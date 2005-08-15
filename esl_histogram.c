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

#include <easel.h>
#include <esl_histogram.h>
#ifdef eslAUGMENT_GUMBEL
#include <esl_gumbel.h>
#endif
#ifdef eslAUGMENT_GEV
#include <esl_gev.h>
#endif
#ifdef eslAUGMENT_STATS
#include <esl_stats.h>
#endif

static int qsort_numerically(const void *xp1, const void *xp2);

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
 *            suited for fitting distributions and goodness-of-fit
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

  h->x     = NULL;
  h->xmin  =  DBL_MAX;
  h->xmax  = -DBL_MAX;
  h->n     = 0;
  h->nalloc= 0;

  h->obs   = NULL;		/* briefly... */
  h->expect= NULL;		/* 'til a Set*() call */
  h->bmin  = bmin;
  h->bmax  = bmax;
  h->nb    = (int)((bmax-bmin)/w);
  h->imin  = h->nb;
  h->imax  = -1;
  h->w     = w;

  h->obs2    = NULL;		/* 'til a Finish() call */
  h->expect2 = NULL;            /* 'til a Set*() call   */
  h->topx    = NULL;
  h->nb2     = 0;

  if ((h->obs = malloc(sizeof(int) * h->nb)) == NULL)
    { esl_histogram_Destroy(h); ESL_ERROR_NULL(eslEMEM, "malloc failed");}
  for (i = 0;  i < h->nb; i++)
    h->obs[i] = 0;

  h->is_full     = FALSE;
  h->is_finished = FALSE;

  return h;
}

/* Function:  esl_histogram_CreateFull()
 * Incept:    SRE, Tue Jul 26 13:19:27 2005 [St. Louis]
 *
 * Purpose:   Alternative form of <esl_histogram_Create()> that 
 *            creates a more complex histogram that will contain not just the 
 *            display histogram, but also keeps track of all
 *            the raw sample values. A full histogram can be used
 *            for fitting distributions and goodness-of-fit 
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
  if (h->obs2   != NULL) free(h->obs2);
  if (h->expect2!= NULL) free(h->expect2);
  if (h->topx   != NULL) free(h->topx);
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

  if (h->is_finished)
    ESL_ERROR(eslEINVAL, "Can't add to a finished histogram");

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

  /* Collate x in the display histogram, reallocating the bins
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
  

/* Function:  esl_histogram_Finish()
 * Incept:    SRE, Tue Jul 26 13:35:14 2005 [St. Louis]
 *
 * Purpose:   Called after all data collation is complete, and no
 *            more calls to <esl_histogram_Add()> will be made.
 *            Do any necessary internal bookkeeping before we do
 *            any printing, plotting, distribution fitting or
 *            setting, or goodness-of-fit testing.
 *            
 * Note:      Currently, this involves just building the 
 *            secondary histogram, with roughly equal-sized bins,
 *            by sorting and partitioning the raw sample values <x>.          
 *            
 *            If there are lots of equal values in <x>, this can
 *            result in some bins containing more counts than
 *            expected, and some bins having zero, because a bin must
 *            contain all values $\leq$ the highest value in the bin.
 *            The goodness-of-fit tests must watch out for these empty bins.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation or reallocation failure.
 */
int
esl_histogram_Finish(ESL_HISTOGRAM *h)
{
  int i;
  int j, prv_j;			/* indices into x[], dividing it into bins */

  /* If we're a display-only histogram, not a full histogram
   * with raw sample data, we don't have anything to do; just return success.
   */
  if (! h->is_full) { h->is_finished = TRUE; return eslOK; }

  /* Sort x such that x[0] is smallest, x[n-1] is largest value.
   */
  qsort((void *) h->x, h->n, sizeof(double), qsort_numerically);

  /* Number of bins for goodness-of-fit tests like G and X^2 
   * is crucial but arbitrary, unfortunately. Some literature suggests
   * using 2*n^{0.4}, which gives:
   *        n    nbins     #/bin
   *    -----    ------   ------
   *     1000      31       32
   *    10000      79      127
   *   100000     200      500
   *  1000000     502     1992
   */
  h->nb2 = 2* (int) pow((double) h->n, 0.4);
  
  /* Allocate for the secondary histogram (expect2 will be done later by
   * a Set() function)
   */
  ESL_MALLOC(h->obs2, sizeof(int)    * h->nb2);
  ESL_MALLOC(h->topx, sizeof(double) * h->nb2);

  /* Partition the sorted data, being careful about skipping through ties.
   */
  prv_j = -1;
  for (i = 0; i < h->nb2; i++)
    {
      j = (int) ((i+1)*h->n/h->nb2) - 1;    /* bin i contains (at least) up to x[j] */
      while (j < h->n-1 && h->x[j+1] == h->x[j]) j++;  /* but also include ties */
      
      h->obs2[i] = j - prv_j; 	/* how many observations got put into bin i */
      h->topx[i] = h->x[j];	/* the highest score x[] in bin i */

      prv_j = j;
    }

  h->is_finished = TRUE;
  return eslOK;
}
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



/* Function:  esl_histogram_GetRank()
 * Incept:    SRE, Thu Jul 28 08:39:52 2005 [St. Louis]
 *
 * Purpose:   Retrieve the <rank>'th highest score from a 
 *            full, finished histogram <h>. <rank> is 1..n, for
 *            n total samples in the histogram.
 *
 * Throws:    <eslEINVAL> if the histogram is display-only,
 *            if it isn't finished, or if <rank> isn't in
 *            the range 1..n.
 */
int
esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank)
{
  if (! h->is_full) 
    ESL_ERROR(eslEINVAL, 
	      "esl_histogram_GetRank() needs a full histogram");
  if (! h->is_finished) 
    ESL_ERROR(eslEINVAL, 
	      "esl_histogram_GetRank() needs a finished, sorted histogram");
  if (rank > h->n)
    ESL_ERROR(eslEINVAL, 
	      "no such rank: not that many scores in the histogram");
  if (rank < 1)
    ESL_ERROR(eslEINVAL, "histogram rank must be a value from 1..n");

  return (h->x[h->n - rank + 1]);
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
 * Incept:    SRE, Fri Aug 12 08:40:37 2005 [St. Louis]
 *
 * Purpose:   Given a full histogram <h> - one with sorted list of
 *            raw scores - output the empirical survival function
 *            (1 - CDF) to an xmgrace XY file <fp>.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if histogram hasn't been finished, or
 *            if it isn't full.
 */
int
esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;

  if (! h->is_full)    
    ESL_ERROR(eslEINVAL, "need full histogram to plot empirical distributions");
  if (! h->is_finished)
    ESL_ERROR(eslEINVAL, "histogram not finished and sorted");

  for (i = h->n-1; i >= 0; i--)	/* sorted w/ low score at 0, high at n-1 */
    fprintf(fp, "%f\t%g\n", h->x[i], (double)(h->n-i)/(double)h->n);
  fprintf(fp, "&\n");
  return eslOK;
}
	      


/*****************************************************************
 * Functions for setting expected histogram frequencies
 *****************************************************************/
#ifdef eslAUGMENT_GUMBEL
int
esl_histogram_SetGumbel(ESL_HISTOGRAM *h, double mu, double lambda)
{
  int    i;
  double x1,x2;

  /* Set expectations in the display histogram (fixed-width bins)
   */
  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      x1 = i * h->w + h->bmin;	   /* scores in bin i are > x1 */
      x2 = (i+1) * h->w + h->bmin; /*     ...and <= x2.        */

      h->expect[i] = h->n *
	(esl_gumbel_cdf(x2, mu, lambda) - esl_gumbel_cdf(x1, mu, lambda));
    }

  /* Optionally (in a full histogram), set expectations in
   * the equal-binsize secondary histogram.
   */
  if (h->is_full)
    {
      if (h->expect2 == NULL)
	ESL_MALLOC(h->expect2, sizeof(double) * h->nb2);

      /* bin[0] is a tail, cdf up to topx[0] */
      h->expect2[0] = (double) h->n * esl_gumbel_cdf(h->topx[0], mu, lambda);
      /* bin[n-1] is a tail, survival function after topx[n-2] */
      h->expect2[h->nb2-1] = (double)h->n * esl_gumbel_surv(h->topx[h->nb2-2], mu, lambda);
      /* remaining bins are differences between cdfs */
      for (i = 1; i < h->nb2-1; i++)
	h->expect2[i] = (double) h->n *
	  (esl_gumbel_cdf(h->topx[i],   mu, lambda) - 
	   esl_gumbel_cdf(h->topx[i-1], mu, lambda));
    }

  return eslOK;
}
#endif /*eslAUGMENT_GUMBEL*/

#ifdef eslAUGMENT_GEV
int
esl_histogram_SetGEV(ESL_HISTOGRAM *h, double mu, double lambda, double alpha)
{
  int    i;
  double x1,x2;

  /* Set expectations in the display histogram (fixed-width bins)
   */
  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      x1 = i * h->w + h->bmin;	   /* scores in bin i are > x1 */
      x2 = (i+1) * h->w + h->bmin; /*     ...and <= x2.        */

      h->expect[i] = (double) h->n *
	(esl_gev_cdf(x2, mu, lambda, alpha) - esl_gev_cdf(x1, mu, lambda, alpha));
    }

  /* Optionally (in a full histogram), set expectations in
   * the equal-binsize secondary histogram.
   */
  if (h->is_full)
    {
      if (h->expect2 == NULL)
	ESL_MALLOC(h->expect2, sizeof(double) * h->nb2);

      /* bin[0] is a tail, cdf up to topx[0] */
      h->expect2[0] = (double) h->n * esl_gev_cdf(h->topx[0], mu, lambda, alpha);
      /* bin[n-1] is a tail, survival function after topx[n-2] */
      h->expect2[h->nb2-1] = (double)h->n * esl_gev_surv(h->topx[h->nb2-2], mu, lambda, alpha);
      /* remaining bins are differences between cdfs */
      for (i = 1; i < h->nb2-1; i++)
	h->expect2[i] = (double) h->n *
	  (esl_gev_cdf(h->topx[i],   mu, lambda, alpha) - 
	   esl_gev_cdf(h->topx[i-1], mu, lambda, alpha));
    }

  return eslOK;
}
#endif /*eslAUGMENT_GEV*/


#ifdef eslAUGMENT_STATS
/* Function:  esl_histogram_GTestGoodness()
 * Incept:    SRE, Tue Jul 26 14:06:02 2005 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, with expectations already set by
 *            a prior Set() call, run a G test on the observed vs.
 *            expected numbers in the binned secondary histogram.
 *            
 *            The G statistic is distributed roughly \chi^2, with
 *            at most nbins-1 degrees of freedom. If the expectations
 *            were set using parameters fit to the observed data,
 *            then there are fewer degrees of freedom; <nfitted>
 *            is the number of fitted parameters, where <nfitted> $\geq 0$.
 *            
 *            Returns <ret_nbins>, the number of valid bins counted
 *            toward the G statistic; <ret_G>, the G statistic; and
 *            <ret_p>, the probability of obtaining a statistic of
 *            at least G, assuming a \chi^2 distribution with
 *            <nbins>-1-<ndeg> degrees of freedom. If <p> is small,
 *            the data are a poor fit to the expected distribution.
 *
 * Args:      h         - full histogram, expectations already Set()
 *            nfitted   - >=0; number of fitted parameters 
 *            ret_nbins - optRETURN: # of bins counted toward G statistic
 *            ret_G     - optRETURN: G statistic
 *            ret_p     - optRETURN: P(>=G)
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if it isn't a full histogram, or if it isn't
 *            sorted and finished. 
 *            <eslERANGE> or <eslECONVERGENCE> may arise on internal
 *            errors in calculating the probability.
 */
int
esl_histogram_GTestGoodness(ESL_HISTOGRAM *h, int nfitted,
			    int *ret_nbins, double *ret_G, double *ret_p)
{
  int    i;
  int    nbins;			/* number of bins counted toward G */
  double G = 0.;		/* the G-statistic */
  double Gp;			/* P(test > G) by chi-square distribution */
  int    status;
  
  if (! h->is_full)     ESL_ERROR(eslEINVAL, "goodness tests need full histogram");
  if (! h->is_finished) ESL_ERROR(eslEINVAL, "goodness tests need finished histogram");

  /* Calculate the G statistic = 2 * log likelihood ratio.
   */
  nbins = 0;
  for (i = 0; i < h->nb2; i++)
    if (h->obs2[i] > 0 && h->expect2[i] > 0)
      {
	G += (double) h->obs2[i] * log((double) h->obs2[i]/ h->expect2[i]);
	nbins++;
      }
  G *= 2.;
  
  /* G is distributed approximately as \chi^2
   */
  if (nbins-1-nfitted >= 0) 
    {
      status = esl_stats_ChiSquaredTest(nbins-1-nfitted, G, &Gp);
      if (status != eslOK) return status;
    }
  else Gp = 0.;

  if (ret_nbins != NULL) *ret_nbins = nbins;
  if (ret_G     != NULL) *ret_G     = G;
  if (ret_p     != NULL) *ret_p     = Gp;
  return eslOK;
}

/* Function:  esl_histogram_ChiSquaredGoodness()
 * Incept:    SRE, Tue Jul 26 14:18:02 2005 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, with expectations already set by
 *            a prior Set() call, run a chi-squared test on the observed vs.
 *            expected numbers in the binned secondary histogram.
 *            
 *            The X^2 statistic is distributed roughly \chi^2, with
 *            at most nbins-1 degrees of freedom. If the expectations
 *            were set using parameters fit to the observed data,
 *            then there are fewer degrees of freedom; <nfitted>
 *            is the number of fitted parameters, where <nfitted> $\geq 0$.
 *            
 *            Returns <ret_nbins>, the number of valid bins counted
 *            toward the X^2 statistic; <ret_X2>, the X^2 statistic; and
 *            <ret_p>, the probability of obtaining a statistic of
 *            at least X2, assuming a \chi^2 distribution with
 *            <nbins>-1-<ndeg> degrees of freedom. If <p> is small,
 *            the data are a poor fit to the expected distribution.
 *
 * Args:      h         - full histogram, expectations already Set()
 *            nfitted   - >=0; number of fitted parameters 
 *            ret_nbins - optRETURN: # of bins counted toward X^2 statistic
 *            ret_X2    - optRETURN: X^2 statistic
 *            ret_p     - optRETURN: P(>=X^2)
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if it isn't a full histogram, or if it isn't
 *            sorted and finished. 
 *            <eslERANGE> or <eslECONVERGENCE> may arise on internal
 *            errors in calculating the probability.
 */
int
esl_histogram_ChiSquaredGoodness(ESL_HISTOGRAM *h, int nfitted, 
				 int *ret_nbins, double *ret_X2, double *ret_p)
{
  int    i;
  int    nbins;			/* number of bins counted toward X^2 */
  double chisq = 0.;		/* the X^2 statistic */
  double delta;			/* obs - exp in a bin */
  double chip;		        /* P(test > X^2) by chi-square distribution */
  int    status;
  
  if (! h->is_full)     ESL_ERROR(eslEINVAL, "goodness tests need full histogram");
  if (! h->is_finished) ESL_ERROR(eslEINVAL, "goodness tests need finished histogram");

  /* Calculate the X^2 statistic = \sum_i (obs_i-exp_i)^2/exp, over
   * all bins containing some minimum size (arbitrarily 5)
   */
  nbins = 0;
  for (i = 0; i < h->nb2; i++)
    if (h->obs2[i] > 0 && h->expect2[i] > 0.)
      {
	delta = h->obs2[i] - h->expect2[i];
	chisq += delta * delta / h->expect2[i];
	nbins++;
      }

  
  /* X^2 is distributed approximately as \chi^2
   */
  if (nbins-1-nfitted >= 0)
    {
      status = esl_stats_ChiSquaredTest(nbins-1-nfitted, chisq, &chip);
      if (status != eslOK) return status;
    }
  else chip = 0.;

  if (ret_nbins != NULL) *ret_nbins = nbins;
  if (ret_X2    != NULL) *ret_X2    = chisq;
  if (ret_p     != NULL) *ret_p     = chip;
  return eslOK;
}
#endif /*eslAUGMENT_STATS*/


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
