/* esl_histogram.c
 * Collecting and displaying histograms.
 * 
 * SRE, Fri Jul  1 13:21:45 2005 [St. Louis]
 * SVN $Id$
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <easel.h>
#include <esl_histogram.h>

#ifdef eslAUGMENT_STATS	 /* stats augmentation gives you goodness-of-fit testing */
#include <esl_stats.h>
#endif



/* Function:  esl_histogram_Create()
 * Incept:    SRE, Fri Jul  1 13:40:26 2005 [St. Louis]
 *
 * Purpose:   Creates and returns a new histogram object, initially
 *            allocated to count scores $>$ <bmin> and $<=$ <xmax> into
 *            bins of width <w>. Thus, a total of <bmax>-<bmin>/<w> bins
 *            are initially created. 
 *            
 *            The bounds <bmin> and <bmax> only need to be initial
 *            guesses.  The histogram object will reallocate itself
 *            dynamically as needed to accomodate scores that exceed
 *            current bounds.
 *
 *            For example, <esl_histogram_Create(-100, 100, 0.5)> 
 *            would init the object to collect scores into 400 bins 
 *            $[-100<x<=-99.5],[-99.5<x<=-99.0]...[99.5<x<=100.0]$.
 *            
 *            <esl_histogram_Create()> creates a simplified histogram
 *            object that collates only the "display" histogram. For
 *            a more complex object that also keeps the raw data samples,
 *            better suited for fitting distributions and goodness-of-fit
 *            testing, use <esl_histogram_CreateFull()>.
 *  
 * Args:      bmin - caller guesses that minimum score will be > bmin
 *            bmax - caller guesses that max score will be <= bmax
 *            w    - size of bins (1.0, for example)
 *            
 * Returns:   ptr to new <ESL_HISTOGRAM> object, which caller is responsible
 *            for free'ing with <esl_histogram_Destroy()>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_HISTOGRAM *
esl_histogram_Create(double bmin, double bmax, double w)
{
  ESL_HISTOGRAM *h = NULL;
  int i;

  if ((h = malloc(sizeof(ESL_HISTOGRAM))) == NULL)    
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  h->xmin      =  DBL_MAX;
  h->xmax      = -DBL_MAX;
  h->n         = 0;
  h->obs       = NULL;		/* briefly... */
  h->expect    = NULL;		/* 'til a Set*() call */
  h->bmin      = bmin;
  h->bmax      = bmax;
  h->nb        = (int)((bmax-bmin)/w);
  h->imin      = h->nb;
  h->imax      = -1;
  h->w         = w;

  h->x         = NULL;
  h->nalloc    = 0;

  h->is_full   = FALSE;
  h->is_sorted = FALSE;

  if ((h->obs = malloc(sizeof(int) * h->nb)) == NULL)
    { esl_histogram_Destroy(h); ESL_ERROR_NULL(eslEMEM, "malloc failed");}
  for (i = 0;  i < h->nb; i++)
    h->obs[i] = 0;
  return h;
}

/* Function:  esl_histogram_CreateFull()
 * Incept:    SRE, Tue Jul 26 13:19:27 2005 [St. Louis]
 *
 * Purpose:   Alternative form of <esl_histogram_Create()> that 
 *            creates a more complex histogram that will contain not just the 
 *            display histogram, but also keeps track of all
 *            the raw sample values. Having a complete vector of raw
 *            samples improves distribution-fitting and goodness-of-fit 
 *            tests. 
 */
ESL_HISTOGRAM *
esl_histogram_CreateFull(double bmin, double bmax, double w)
{
  ESL_HISTOGRAM *h = esl_histogram_Create(bmin, bmax, w);
  if (h == NULL) return NULL;

  h->n      = 0;		/* make sure */
  h->nalloc = 128;		/* arbitrary initial allocation size */
  h->x      = malloc(sizeof(double) * h->nalloc);
  if (h->x == NULL)
    { esl_histogram_Destroy(h); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }

  h->is_full = TRUE;
  return h;
}


/* Function:  esl_histogram_Destroy()
 * Incept:    SRE, Sat Jul  2 19:41:17 2005 [St. Louis]
 *
 * Purpose:   Frees an <ESL_HISTOGRAM> object <h>.
 */
void
esl_histogram_Destroy(ESL_HISTOGRAM *h)
{
  if (h ==  NULL) return;
  if (h->x      != NULL) free(h->x);
  if (h->obs    != NULL) free(h->obs); 
  if (h->expect != NULL) free(h->expect);
  free(h);
  return;
}


/* Function:  esl_histogram_Add()
 * Incept:    SRE, Sat Jul  2 19:41:45 2005 [St. Louis]
 *
 * Purpose:   Adds score <x> to a histogram <h>.
 *           
 *            The histogram will be automatically reallocated as
 *            needed if the score is smaller or larger than the
 *            current allocated bounds.  
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 *            <eslEINVAL> if the histogram was already finished.
 */
int
esl_histogram_Add(ESL_HISTOGRAM *h, double x)
{
  int b,i;			/* what bin we're in                       */
  int nnew;			/* # of new bins created by a reallocation */

  h->is_sorted = FALSE;		/* not any more! */

  /* If we're a full histogram, then we keep the raw x value,
   * reallocating as needed.
   */
  if (h->is_full) 
    {
      if (h->nalloc == h->n) 
	{
	  h->nalloc *= 2;
	  h->x = realloc(h->x, sizeof(double) * h->nalloc);
	  if (h->x == NULL) ESL_ERROR(eslEMEM, "reallocation failed");
	}
      h->x[h->n] = x;
    }

  /* Collate x in the histogram, reallocating the bins
   * if needed.
   */
  b = h->nb - (int) floor((h->bmax - x)/h->w);	

  /* Reallocate below?
   */
  if (b < 0) 
    {				
      nnew = -b*2;	/* overallocate by 2x */
      h->obs = realloc(h->obs, sizeof(int) * (nnew+ h->nb));
      if (h->obs == NULL) { ESL_ERROR(eslEMEM, "reallocation failed"); }
      
      memmove(h->obs+nnew, h->obs, sizeof(int) * h->nb);
      h->nb    += nnew;
      b        += nnew;
      h->bmin  -= nnew*h->w;
      h->imin  += nnew;
      if (h->imax > -1) h->imax += nnew;
      for (i = 0; i < nnew; i++) h->obs[i] = 0;
    }
  /* Reallocate above?
   */
  else if (b >= h->nb)
    {
      nnew = (b-h->nb+1) * 2; /* 2x overalloc */
      
      h->obs = realloc(h->obs, sizeof(int) * (nnew+ h->nb));
      if (h->obs == NULL) { ESL_ERROR(eslEMEM, "reallocation failed"); }

      for (i = h->nb; i < h->nb+nnew; i++)
	h->obs[i] = 0;      
      if (h->imin == h->nb) h->imin+=nnew;
      h->bmax  += nnew*h->w;
      h->nb    += nnew;
    }

  /* Bump the bin counter, and whatever other data we're keeping.
   */
  h->obs[b]++;
  h->n++;
  if (b > h->imax) h->imax = b;
  if (b < h->imin) h->imin = b;
  if (x > h->xmax) h->xmax = x;
  if (x < h->xmin) h->xmin = x;
  return eslOK;
}
  

/* qsort_numerically:
 * this'll be used in the next function.
 */
static int
qsort_numerically(const void *xp1, const void *xp2)
{
  double x1;
  double x2; 
  x1 = * (double *) xp1;
  x2 = * (double *) xp2;
  if (x1 < x2) return -1;
  if (x1 > x2) return 1;
  return 0;
}

/* Function:  esl_histogram_Sort()
 * Incept:    SRE, Thu Aug 18 10:45:46 2005 [St. Louis]
 *
 * Purpose:   Sort the raw scores in a full histogram.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if you try to sort a default, non-full histogram that
 *            doesn't have raw scores stored in it.
 */
int
esl_histogram_Sort(ESL_HISTOGRAM *h)
{
  if (! h->is_full) ESL_ERROR(eslEINVAL, "not a full histogram, nothing to sort");

  if (! h->is_sorted) 
    {
      qsort((void *) h->x, h->n, sizeof(double), qsort_numerically);
      h->is_sorted = TRUE;
    }
  return eslOK;
}


/* Function:  esl_histogram_GetScoreAtRank()
 * Incept:    SRE, Thu Jul 28 08:39:52 2005 [St. Louis]
 *
 * Purpose:   Retrieve the <rank>'th highest score from a 
 *            full, finished histogram <h>. <rank> is 1..n, for
 *            n total samples in the histogram; return it through
 *            <ret_x>.
 *            
 *            If the raw scores aren't sorted, this function sorts
 *            them all (an NlogN operation) before returning the score
 *            you're looking for.
 *
 * Throws:    <eslEINVAL> if the histogram is display-only,
 *            if it isn't finished, or if <rank> isn't in
 *            the range 1..n.
 */
int
esl_histogram_GetScoreAtRank(ESL_HISTOGRAM *h, int rank, double *ret_x)
{
  if (! h->is_full) 
    ESL_ERROR(eslEINVAL, 
	      "esl_histogram_GetRank() needs a full histogram");
  if (rank > h->n)
    ESL_ERROR(eslEINVAL, 
	      "no such rank: not that many scores in the histogram");
  if (rank < 1)
    ESL_ERROR(eslEINVAL, "histogram rank must be a value from 1..n");

  esl_histogram_Sort(h);	/* make sure */
  *ret_x = h->x[h->n - rank + 1];
  return eslOK;
}


/* Function:  esl_histogram_GetBinBounds()
 * Incept:    SRE, Wed Aug 17 08:02:17 2005 [St. Louis]
 *
 * Purpose:   For bin index <whichbin> in the histogram <h>,
 *            retrieve the bounds <ret_low>..<ret_high> for
 *            values in that bin. That is, all values $x$ in bin
 *            <whichbin> satisfy $l < x \leq h$.
 *
 *            Also returns the bin width in <ret_delta>;
 *            $h = l + \delta$.
 *            
 *            All three returned values are optional. Pass NULL
 *            for any of them that you don't want.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_GetBinBounds(ESL_HISTOGRAM *h, int whichbin, 
			   double *ret_low, double *ret_high, double *ret_delta)
{
  if (ret_low  != NULL) *ret_low  = h->w * whichbin     + h->bmin;
  if (ret_high != NULL) *ret_high = h->w * (whichbin+1) + h->bmin;
  if (ret_delta!= NULL) *ret_delta= h->w;
  return eslOK;
}


/* Function:  esl_histogram_Print() 
 * Incept:    SRE, Sat Jul  2 16:03:37 2005 [St. Louis]
 *
 * Purpose:   Print a "prettified" display histogram <h> to a file pointer <fp>.
 *            Deliberately a look-and-feel clone of Bill Pearson's 
 *            excellent FASTA output.
 * 
 * Args:      fp     - open file to print to (stdout works)
 *            h      - histogram to print
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_Print(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;
  int maxbar;
  int imode;
  int units;
  int num;
  char buffer[81];		/* output line buffer */
  int  pos;			/* position in output line buffer */
  int  ilowbound, lowcount;	/* cutoffs on the low side  */
  int  ihighbound, highcount;	/* cutoffs on the high side */
  int  emptybins = 3;

  /* Find out how we'll scale the histogram.  We have 59 characters to
   * play with on a standard 80-column terminal display: leading "%5d
   * %6d %6d|" occupies 20 chars.  Save the peak position, we'll use
   * it later.
   */
  maxbar = 0;
  for (i = 0; i < h->nb; i++)
    if (h->obs[i] > maxbar) 
      {
	maxbar  = h->obs[i];     /* max height    */
	imode   = i;
      }

  /* Truncate histogram display on both sides, ad hoc fashion.
   * Start from the peak; then move out until we see <emptybins> empty bins,
   * and stop.
   */
  for (num = 0, ihighbound = imode; ihighbound < h->imax; ihighbound++)
    {
      if (h->obs[ihighbound] > 0) { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }
  for (num = 0, ilowbound = imode; ilowbound > h->imin; ilowbound--)
    {
      if (h->obs[ilowbound] > 0)  { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }

				/* collect counts outside of bounds */
  for (lowcount = 0, i = h->imin; i < ilowbound; i++)
    lowcount += h->obs[i];
  for (highcount = 0, i = h->imax; i > ihighbound; i--)
    highcount += h->obs[i];

				/* maxbar might need to be raised now; then set our units  */
  if (lowcount  > maxbar) maxbar = lowcount;
  if (highcount > maxbar) maxbar = highcount;
  units = ((maxbar-1)/ 59) + 1;

  /* Print the histogram
   */
  fprintf(fp, "%5s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  fprintf(fp, "%5s %6s %6s\n", "-----", "---", "---");
  buffer[80] = '\0';
  buffer[79] = '\n';
  for (i = h->imin; i <= h->imax; i++)
    {
      memset(buffer, ' ', 79 * sizeof(char));
      x = i*h->w + h->bmin;

      /* Deal with special cases at edges
       */
      if      (i < ilowbound)  continue;
      else if (i > ihighbound) continue;
      else if (i == ilowbound && i != h->imin) 
	{
	  sprintf(buffer, "<%4.0f %6d %6s|", x+h->w, lowcount, "-");
	  if (lowcount > 0) {
	    num = 1+(lowcount-1) / units;
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}
      else if (i == ihighbound && i != h->imax)
	{
	  sprintf(buffer, ">%4.0f %6d %6s|", x, highcount, "-");
	  if (highcount > 0) {
	    num = 1+(highcount-1) / units;
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}

      /* Deal with most cases
       */
      if (h->expect != NULL) 
	sprintf(buffer, "%5.0f %6d %6d|", 
		x, h->obs[i], (int) h->expect[i]);
      else
	sprintf(buffer, "%5.0f %6d %6s|", x, h->obs[i], "-");
      buffer[20] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->obs[i] > 0) {
	num = 1 + (h->obs[i]-1) / units;
	for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       */
      if (h->expect != NULL && h->expect[i] > 0.)
	{
	  pos = 20 + (int)(h->expect[i]-1) / units;
	  if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	  buffer[pos] = '*';
	}

      /* Print the line
       */
      fputs(buffer, fp);
    }

  return eslOK;
}
  

/* Function:  esl_histogram_Plot()
 * Incept:    SRE, Sat Jul  2 19:43:20 2005 [St. Louis]
 *
 * Purpose:   Print a histogram <h> to open file ptr <fp>, 
 *            as an XY file suitable for input to the xmgrace
 *            graphing program.
 *
 */
int
esl_histogram_Plot(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;

  /* First data set is the observed histogram
   */
  for (i = h->imin; i <= h->imax; i++)
    if (h->obs[i] > 0)
      {
	x = i*h->w + h->bmin;
	fprintf(fp, "%f %d\n", x, h->obs[i]);
      }
  fprintf(fp, "&\n");

  /* Second data set is the theoretical histogram
   */
  if (h->expect != NULL)
    {
      for (i = 0; i < h->nb; i++)
	if (h->expect[i] > 0.)
	  {
	    x = i*h->w + h->bmin;
	    fprintf(fp, "%.2f %.2f\n", x, h->expect[i]);
	  }
      fprintf(fp, "&\n");
    }
  return eslOK;
}


/* Function:  esl_histogram_PlotSurvival()
 * Incept:    SRE, Wed Aug 17 08:26:10 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h>, output the empirical survival function
 *            (1-CDF, P(X>x)) to an xmgrace XY file <fp>.
 *
 *            If raw scores are available (in a full histogram) it uses those
 *            for a higher-resolution plot. If not, it uses the binned scores
 *            and produces a lower-resolution plot.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h)
{
  int i;
  int c = 0;
  
  if (h->is_full)    		/* use all (raw) scores? */
    {
      esl_histogram_Sort(h);
      for (i = h->n-1; i >= 0; i--)	/* sorted w/ low score at 0, high at n-1 */
	fprintf(fp, "%f\t%g\n", h->x[i], (double)(h->n-i)/(double)h->n);
    }
  else				/* else, use binned counts */
    {
      for (i = h->imax; i >= h->imin; i--)
	{
	  c   += h->obs[i];
	  fprintf(fp, "%f\t%f\n", i*h->w + h->bmin, (double) c / (double) h->n);
	}
    }
  fprintf(fp, "&\n");
  return eslOK;
}

/* Function:  esl_histogram_PlotQQ()
 * Incept:    SRE, Sat Aug 20 14:15:01 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> containing an empirically observed
 *            distribution, and a pointer to a function <(*invcdf)()>
 *            for an expected inverse cumulative distribution
 *            function conditional on some parameters <params>;
 *            output a Q-Q plot in xmgrace XY format to file <fp>.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotQQ(FILE *fp, ESL_HISTOGRAM *h, 
		     double (*invcdf)(double x, void *params), void *params)
{
  int i;
  int c = 0;
  double cdf;
  double bi;
  
  if (h->is_full)    		/* use all (raw) scores? */
    {
      esl_histogram_Sort(h);

      for (i = 0; i <= h->n-2; i++)   /* avoid last sample where cdf=1.0 */
	{
	  cdf = (double)(i+1)/(double)h->n;
	  fprintf(fp, "%f\t%f\n", h->x[i], (*invcdf)(cdf, params));
	}
    }
  else				/* else, use binned counts */
    {
      for (i = h->imin; i < h->imax; i--) /* again, avoid last bin, cdf=1.0 */
	{
	  c   += h->obs[i];
	  cdf =  (double) c / h->n;
	  esl_histogram_GetBinBounds(h, i, NULL, &bi, NULL);
	  fprintf(fp, "%f\t%f\n", bi, (*invcdf)(cdf, params));
	}
    }
  fprintf(fp, "&\n");

  /* this plots a 45-degree expected QQ line:
   */
  fprintf(fp, "%f\t%f\n", h->xmin, h->xmin);
  fprintf(fp, "%f\t%f\n", h->xmax, h->xmax);
  fprintf(fp, "&\n");

  return eslOK;
}



/* Function:  esl_histogram_SetExpect()
 * Incept:    SRE, Wed Aug 17 17:36:58 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> containing some number of empirically
 *            observed binned counts, and a pointer to a function <(*cdf)()>
 *            that describes the expected cumulative distribution function 
 *            (CDF) conditional on some parameters <params>;
 *            calculate the expected counts in each bin of the histogram,
 *            and hold that information internally in the structure.
 *            
 *            Expected counts (when calculated) are displayed by 
 *            <esl_histogram_Print()> and <esl_histogram_Plot()>.
 *            
 *            The caller provides a function <(*cdf)()> that calculates
 *            the CDF via a generic interface, taking only two
 *            arguments: a quantile <x> and a void pointer to whatever
 *            parameters it needs, which it will cast and interpret.
 *            The <params> void pointer to the given parameters will
 *            just be passed along to the <(*cdf)()> function. The
 *            caller will probably implement this <(*cdf)()> function as
 *            a wrapper around its real CDF function that takes
 *            explicit (non-void-pointer) arguments.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_histogram_SetExpect(ESL_HISTOGRAM *h, 
			double (*cdf)(double x, void *params), void *params)
{
  int    i;
  double x1,x2;

  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      x1 = i * h->w + h->bmin;	   /* scores in bin i are > x1 */
      x2 = (i+1) * h->w + h->bmin; /*     ...and <= x2.        */
      h->expect[i] = h->n * (*cdf)(x2, params) - (*cdf)(x1, params);
    }
  return eslOK;
}


#ifdef eslAUGMENT_STATS
/* Function:  esl_histogram_Goodness()
 * Incept:    SRE, Wed Aug 17 12:46:05 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> with observed counts,
 *            and a pointer to a function <*cdf()> with parameters <params>,
 *            describing the expected cumulative probability distribution 
 *            function (CDF), of which <nfitted> ($\geq 0$) were fitted (and
 *            thus should be subtracted from the degrees of freedom);
 *            Perform a G-test and/or a chi-squared test for goodness of 
 *            fit between observed and expected, and optionally return
 *            the number of bins the data were sorted into
 *            (<ret_bins>), the G statistic and its probability (<ret_G> and
 *            <ret_Gp>), and the X$^2$ statistic and its probability
 *            (<ret_X2> and <ret_X2p>). 
 *            
 *            The function that calculates the CDF has a generic
 *            interface, taking only two arguments: a quantile <x> and
 *            a void pointer to whatever parameters it needs, which it
 *            will cast and interpret.  The <params> void pointer to
 *            the given parameters will just be passed along to the
 *            <(*cdf)()> function. The caller will probably implement
 *            a <(*cdf)()> function as a wrapper around a real CDF
 *            function with explicit (non-void-pointer) arguments.
 *            
 *            If a goodness-of-fit probability is less than some threshold
 *            (usually taken to be 0.01 or 0.05), that is considered to
 *            be evidence that the observed data are unlikely to be consistent
 *            with the tested distribution.
 *            
 *            The two tests should give similar probabilities. In practice,
 *            both goodness-of-fit tests are sensitive to arbitrary choices
 *            in how the data are binned. Neither test is on a sound
 *            theoretical footing.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslERANGE> or <eslECONVERGENCE> may arise on internal
 *            errors in calculating the probability.
 */
int
esl_histogram_Goodness(ESL_HISTOGRAM *h, 
		       double (*cdf)(double x, void *params), void *params,
		       int nfitted, int use_bindata,
		       int *ret_nbins,
		       double *ret_G,  double *ret_Gp,
		       double *ret_X2, double *ret_X2p)
{
  int     *obs;			/* observed in bin i, [0..nb-1]   */
  double  *exp;			/* expected in bin i, [0..nb-1    */
  double  *topx;		/* all values in bin i <= topx[i] */
  int      nb;			/* # of bins                      */
  int      minc;		/* target # of counts/bin         */
  int      i,b;
  int      sum;
  double   G, Gp;
  double   X2, X2p;
  double   tmp;
  int      status;

  /* Figure out how many bins we'd like to have, then allocate.
   * Number of bins for goodness-of-fit tests like G and X^2 
   * is crucial but arbitrary, unfortunately. Some literature suggests
   * using 2*n^{0.4}, which gives:
   *        n    nbins     #/bin
   *    -----    ------   ------
   *     1000      31       32
   *    10000      79      127
   *   100000     200      500
   *  1000000     502     1992
   */
  nb   = 2* (int) pow((double) h->n, 0.4); /* "desired" nb. */
  minc = h->n / (2*nb);		/* arbitrarily set min = 1/2 of the target # */
  ESL_MALLOC(obs,  sizeof(int)    * (nb*2+1)); /* final nb must be <= 2*nb+1 */
  ESL_MALLOC(exp,  sizeof(double) * (nb*2+1));
  ESL_MALLOC(topx, sizeof(double) * (nb*2+1));

  /* Determine the observed counts in each bin;
   * If we have raw counts, sort and use them (unless overridden by use_bindata)
   * If not, use the binned histogram. 
   * In either case, sweep left to right on the histogram bins,
   * collecting sum of counts, dropping the sum into the next bin 
   * whenever we have more than <minc> counts.
   * In the case of the raw counts, be careful that ties all go into
   * the same bin.
   */
  if (! use_bindata && h->is_full)
    {	/* collate raw counts */
      esl_histogram_Sort(h);
      sum = 0;
      i = 0;
      for (b = 0; b < h->n; b++) /* "b" here is a counter in the raw data */
	{
	  sum++;
	  if (sum >= minc) { /* enough? then drop them in bin i */
	    {
	      while (b < h->n-1 && h->x[b+1] == h->x[b]) { sum++; b++; } /* ties */
	      assert(i < nb*2+1);
	      obs[i]  = sum;
	      topx[i] = h->x[b];
	      sum     = 0;
	      i++;
	    }
	  }
	}
      obs[i-1] += sum;		/* add the remaining right tail to the last bin */
      topx[i-1] = h->x[h->n-1]; /* by def'n */
      nb        = i;
    }
  else
    {	/* merge histogram bins */
      sum = 0;
      i   = 0;
      for (b = h->imin; b <= h->imax; b++) 
	{
	  sum += h->obs[b];
	  if (sum >= minc) {	/* if we have enough counts, drop them in i */
	    assert(i < (nb*2+1));
	    obs[i]  = sum;
	    topx[i] = h->w*(b+1) + h->bmin;
	    sum     = 0;
	    i++;
	  }
	}
      obs[i-1]  += sum;		/* add the right tail to our final bin        */
      topx[i-1]  = h->w * (h->imax+1) + h->bmin;
      nb         = i;		/* nb is now the actual # of bins, not target */
    }

  /* Determine the expected counts in each bin.
   *  bin 0    is the left tail, <= topx[0];
   *  bin nb-1 is the right tail, > topx[nb-2],  1-P(<= topx[nb-2]).
   *  others are   P(<= topx[b]) - P(<topx[b-1]).
   */
  exp[0]    = (double) h->n * (*cdf)(topx[0], params);
  exp[nb-1] = (double) h->n * (1 - (*cdf)(topx[nb-2], params));
  for (i = 1; i < nb-1; i++)
    exp[i] = (double) h->n *
      ((*cdf)(topx[i], params) - (*cdf)(topx[i-1], params));
  
  /* Calculate the G statistic: 2 * LLR  */
  G = 0.;
  for (i = 0; i < nb; i++)
    {
      if (exp[i] == 0) G = eslINFINITY;
      else             G += (double) obs[i] * log ((double) obs[i] / exp[i]);
    }
  G *= 2;
  
  /* G is distributed approximately as \chi^2 */
  if (nb-nfitted-1 >= 0 && G != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted-1, G, &Gp);
      if (status != eslOK) return status;
    }
  else Gp = 0.;

  /* Calculate the X^2 statistic: \sum (obs_i - exp_i)^2 / exp_i */
  X2 = 0.;
  for (i = 0; i < nb; i++)
    {
      if (exp[i] == 0) X2 = eslINFINITY;
      else  {
	tmp = obs[i] - exp[i];
	X2 += tmp*tmp / exp[i];
      }
    }

  /* X^2 is distributed approximately chi^2  */
  if (nb-nfitted-1 >= 0 && X2 != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted-1, X2, &X2p);
      if (status != eslOK) return status;
    }
  else X2p = 0.;

  if (ret_nbins != NULL) *ret_nbins = nb;
  if (ret_G     != NULL) *ret_G     = G;
  if (ret_Gp    != NULL) *ret_Gp    = Gp;
  if (ret_X2    != NULL) *ret_X2    = X2;
  if (ret_X2p   != NULL) *ret_X2p   = X2p;
  free(obs);
  free(exp);
  free(topx);
  return eslOK;
}
#endif /* eslAUGMENT_STATS */





/*****************************************************************
 * Example main()
 *****************************************************************/
#ifdef eslHISTOGRAM_EXAMPLE
/*::cexcerpt::histogram_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslHISTOGRAM_EXAMPLE esl_histogram.c esl_random.c easel.c
 * run:     ./example 
 */
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_Create(-100, 100, 1);
  int    nsamples    = 1000;
  double mean        = 20.0;
  double stddev      = 10.0;
  int    i;
  double x;

  for (i = 0; i < nsamples; i++) {
    x = esl_rnd_Gaussian(r, mean, stddev);
    esl_histogram_Add(h, x);
  }
  esl_histogram_Finish(h);

  esl_histogram_Print(stdout, h);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example::end::*/
#endif /*eslHISTOGRAM_EXAMPLE*/
