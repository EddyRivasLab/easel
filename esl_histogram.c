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
  h->bmin      = bmin;
  h->bmax      = bmax;
  h->nb        = (int)((bmax-bmin)/w);
  h->imin      = h->nb;
  h->imax      = -1;
  h->w         = w;

  h->x         = NULL;
  h->nalloc    = 0;

  h->expect    = NULL;		/* 'til a Set*() call */
  h->phi       = 0.;
  h->cmin      = h->imin;
  h->z         = 0;
  h->Nc        = 0;
  h->No        = 0;
  h->Nx        = 0;

  h->is_full       = FALSE;
  h->is_sorted     = FALSE;
  h->dataset_is    = COMPLETE;
  h->fit_describes = COMPLETE_FIT;

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
 *            <eslEINVAL> for cases where something has been done
 *            to the histogram that requires it to be 'finished', and
 *            adding more data is prohibited; for example, 
 *            if censoring information has already been set.
 */
int
esl_histogram_Add(ESL_HISTOGRAM *h, double x)
{
  int b,i;			/* what bin we're in                       */
  int nnew;			/* # of new bins created by a reallocation */

  /* Censoring info must only be set on a finished histogram */
  if (h->dataset_is != COMPLETE)
    ESL_ERROR(eslEINVAL, "can't add more data to this histogram");

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

  /* Which bin will we want to put x into?
   */
  b = esl_histogram_Score2Bin(h,x);

  /* Make sure we have that bin...
   * Reallocate below?
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
      h->cmin  += nnew;
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
      if (h->imin == h->nb) { /* boundary condition of no data yet*/
	h->imin+=nnew; 
	h->cmin+=nnew;
      }
      h->bmax  += nnew*h->w;
      h->nb    += nnew;
    }

  /* Bump the bin counter, and whatever other data we're keeping.
   * Since censoring is "virtual" (we still keep the complete
   * dataset), we just have to keep track of the number of
   * censored samples, z, if it's already a censored histogram.
   */
  h->obs[b]++;
  h->n++;
  h->Nc++;
  h->No++;
  h->Nx++;

  if (b > h->imax) h->imax = b;
  if (b < h->imin) { h->imin = b; h->cmin = b; }
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
 * Purpose:   Sort the raw scores in a full histogram. Has no effect on a
 *            normal histogram, or on a full histogram that is already sorted.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_Sort(ESL_HISTOGRAM *h)
{
  if (! h->is_full) return eslOK;
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
 *            If the raw scores aren't sorted, they are sorted
 *            first (an NlogN operation).
 *
 * Throws:    <eslEINVAL> if the histogram is display-only,
 *            or if <rank> isn't in the range 1..n.
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


/* Function:  esl_histogram_TrueCensoring()
 * Incept:    SRE, Tue Aug 23 10:00:14 2005 [St. Louis]
 *
 * Purpose:   Declare that the dataset collected in <h> is known to be a
 *            censored distribution, where <z> samples were unobserved because
 *            they had values $\leq$ some threshold <phi>.
 *            
 *            No more data can be added to the histogram with <_Add()>
 *            after censoring information has been set.
 *            
 *            This function is for "true" censored datasets, where
 *            the histogram truly contains no observed points
 *            $x \leq phi$. It's the caller's responsibility to
 *            make sure that it didn't <_Add()> any points
 *            $x \leq \phi$ to the histogram.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if you try to set <phi> to a value that is
 *            greater than the minimum <x> stored in the histogram.
 */
int
esl_histogram_TrueCensoring(ESL_HISTOGRAM *h, int z, double phi)
{
  if (phi > h->xmin) ESL_ERROR(eslEINVAL, "no uncensored x can be <= phi");

  h->phi         = phi;
  h->cmin        = h->imin;
  h->z           = z;
  h->Nc          = h->n + z;
  h->No          = h->n;
  h->Nx          = h->n + z;
  h->dataset_is  = TRUE_CENSORED;
  return eslOK;
}

/* Function:  esl_histogram_VirtCensorByValue()
 * Incept:    SRE, Tue Aug 23 09:01:10 2005 [St. Louis]
 *
 * Purpose:   Suggest a censoring threshold <phi> to split a histogram <h>
 *            into "unobserved" data (values $\leq \phi$) and "observed" 
 *            data (values $> \phi$). 
 *
 *            The suggested <phi> is revised downwards to a $\phi$ at the next 
 *            bin lower bound, because operations on binned data in <h>
 *            need to know unambiguously whether all the data in a bin
 *            will be counted as observed or unobserved. 
 *
 *            Any data point $x_i \leq phi$ is then considered to be
 *            in the censored region for purposes of calculating
 *            expected counts and goodness-of-fit tests. 
 *            
 *            No more data can be added to the histogram after
 *            censoring information has been set.
 *            
 *            This function defines a "virtual" censoring: the
 *            histogram actually contains complete data in <h->obs[]> (and
 *            <h->x>, for a full histogram with raw samples), but only
 *            the <h->Nc>-<h->z> samples above the threshold $\phi$
 *            are counted toward fitting distributions, calculating
 *            expected counts, and running goodness-of-fit tests.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_VirtCensorByValue(ESL_HISTOGRAM *h, double phi)
{
  int b;

  /* Usually, put true phi at the next bin lower bound, but
   * watch for a special case where phi is already exactly equal to a 
   * bin upper bound.
   */
  h->cmin = esl_histogram_Score2Bin(h,phi);
  if (phi == esl_histogram_Bin2UBound(h,h->cmin)) h->phi = phi;
  else   h->phi  = esl_histogram_Bin2LBound(h, h->cmin);

  h->z    = 0;
  for (b = h->imin; b < h->cmin; b++)
    h->z += h->obs[b];
  h->Nc         = h->n;		/* (redundant) */
  h->No         = h->n - h->z;
  h->Nx         = h->n;        	/* (redundant) */
  h->dataset_is = VIRTUAL_CENSORED;

  esl_histogram_Sort(h); /* uncensored raw tail now starts at h->x+h->z */
  return eslOK;
}

/* Function:  esl_histogram_VirtCensorByMass()
 * Incept:    SRE, Tue Aug 23 08:10:39 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> (with or without raw data samples),
 *            find a cutoff score that at least a fraction <tfrac> of the samples
 *            exceed. This threshold is stored internally in the histogram
 *            as <h->phi>. The number of virtually censored samples (to the left,
 *            with scores $<= \phi$) is stored internally in <h->z>.
 *            
 *            The identified cutoff score must be a lower bound for some bin
 *            (bins can't be partially censored). The censored mass
 *            will thus usually be a bit greater than <tfrac>, as the
 *            routine will find the highest satisfactory <h->phi>. The
 *            narrower the bin widths, the more accurately the routine
 *            will be able to satisfy the requested <frac>. The caller
 *            can figure out how much tail mass was really left 
 *            by calculating <(h->n-h->z)/h->n>.
 *            
 *            This function defines a virtually censored dataset.
 *            The complete data are still in <h->obs[]> (and possibly
 *            <h->x[]>. <h->n> (the number of stored samples) will
 *            equal <h->Nc> (the size of the complete data), but only
 *            the <Nc-z> "observed" data points are counted towards
 *            fitting distributions, calculating expected counts, and
 *            running goodness-of-fit tests.
 *            
 *            After calling <_CensorByMass()> on a full histogram
 *            (with raw scores), the caller can retrieve the sorted
 *            censored data as <h->x+h->z>, which is a (partial)
 *            vector containing <h->Nc-h->z> numbers, all satisfying
 *            $x_i > \phi$. The caller can then call a censored
 *            distribution fitting method on this dataset.
 *            
 *            Additionally, after calling <_CensorByMass()> on a
 *            histogram, the index of the first uncensored bin is in
 *            <h->cmin>. That is, the censored bins are <0..h->cmin-1>
 *            and the uncensored bins are <h->cmin..nb-1>; or
 *            alternatively, for the range of bins that contain
 *            counts, <h->imin..h->cmin-1> are censored and
 *            <h->cmin..h->imax>.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_VirtCensorByMass(ESL_HISTOGRAM *h, double tfrac)
{
  int b;
  int sum = 0;
	    
  for (b = h->imax; b >= h->imin; b--)
    {
      sum += h->obs[b];
      if (sum >= (tfrac * (double)h->n)) break;
    }

  h->phi         = esl_histogram_Bin2LBound(h,b);
  h->z           = h->n - sum;
  h->cmin        = b;
  h->Nc          = h->n;	/* (redundant) */
  h->No          = h->n - h->z;
  h->Nx         = h->n;        	/* (redundant) */
  h->dataset_is  = VIRTUAL_CENSORED;

  esl_histogram_Sort(h); /* uncensored raw tail now starts at h->x+h->z */
  return eslOK;
}

/* Function:  esl_histogram_SetTailfitting()
 * Incept:    SRE, Tue Aug 23 10:59:11 2005 [St. Louis]
 *
 * Purpose:   Inform a histogram that the expected fit (and subsequent
 *            goodness-of-fit testing and plotting) will only be to the
 *            <h->Nc-h->z> samples in the uncensored tail: that is,
 *            the expected distribution is only appropriate for 
 *            describing the tail, like perhaps an exponential tail. 
 *            
 *            This affects how expected numbers are calculated. If
 *            a tail fit is declared, expected numbers in the tail
 *            are calculated as <h->Nc-h->z> times the expected density.
 *            Otherwise, expected numbers in the tail are calculated
 *            as  <h->Nc> times the expected density.
 */
int
esl_histogram_SetTailfitting(ESL_HISTOGRAM *h)
{
  h->Nx            = h->Nc - h->z; 
  h->fit_describes = TAIL_FIT;
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
 *            Respects any censoring information that has been
 *            set, and whether tail fitting is been declared.
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
  double ai,bi;			/* ai < x <= bi : lower,upper bounds in bin */

  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      ai = esl_histogram_Bin2LBound(h, i);
      bi = esl_histogram_Bin2UBound(h, i);
      if (h->dataset_is == COMPLETE)
	h->expect[i] = h->Nx * ( (*cdf)(bi, params) - (*cdf)(ai, params) );
      else /* either virtual or true censoring: beware the phi limit */
	{
	  if (ai < h->phi) ai = h->phi;
	  if (i >= h->cmin)
	    h->expect[i] = h->Nx * ((*cdf)(bi, params) - (*cdf)(ai, params));
	  else
	    h->expect[i] = 0;
	}
    }
  return eslOK;
}

/* Function:  esl_histogram_Print() 
 * Incept:    SRE, Sat Jul  2 16:03:37 2005 [St. Louis]
 *
 * Purpose:   Print a "prettified" display histogram <h> to a file pointer <fp>.
 *            Deliberately a look-and-feel clone of Bill Pearson's 
 *            excellent FASTA output.
 *            
 *            This will only work well if the bin width (w) is 0.1 or more,
 *            because the score labels are only shown to one decimal point.
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

  /* Find out how we'll scale the histogram.  We have 58 characters to
   * play with on a standard 80-column terminal display: leading "%6.1f
   * %6d %6d|" occupies 21 chars.  Save the peak position, we'll use
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
  units = ((maxbar-1)/ 58) + 1;

  /* Print the histogram
   */
  fprintf(fp, "%6s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  fprintf(fp, "%6s %6s %6s\n", "-----", "---", "---");
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
	  sprintf(buffer, "<%5.1f %6d %6s|", x+h->w, lowcount, "-");
	  if (lowcount > 0) {
	    num = 1+(lowcount-1) / units;
	    for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}
      else if (i == ihighbound && i != h->imax)
	{
	  sprintf(buffer, ">%5.1f %6d %6s|", x, highcount, "-");
	  if (highcount > 0) {
	    num = 1+(highcount-1) / units;
	    for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}

      /* Deal with most cases
       */
      if (h->expect != NULL) 
	sprintf(buffer, "%6.1f %6d %6d|", 
		x, h->obs[i], (int) h->expect[i]);
      else
	sprintf(buffer, "%6.1f %6d %6s|", x, h->obs[i], "-");
      buffer[21] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->obs[i] > 0) {
	num = 1 + (h->obs[i]-1) / units;
	for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       * (The test > 0. also suffices to remove any censored region.)
       */
      if (h->expect != NULL && h->expect[i] > 0.)
	{
	  pos = 21 + (int)(h->expect[i]-1) / units;
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
	x = esl_histogram_Bin2LBound(h,i);
	fprintf(fp, "%f %d\n", x, h->obs[i]);
      }
  fprintf(fp, "&\n");

  /* Second data set is the theoretical (expected) histogram
   */
  if (h->expect != NULL)
    {
      for (i = 0; i < h->nb; i++)
	if (h->expect[i] > 0.)	/* >0 test suffices to remove censored region */
	  {
	    x = esl_histogram_Bin2LBound(h,i);
	    fprintf(fp, "%.2f %g\n", x, h->expect[i]);
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
  double lastx, delta; /* used to reduce the density of points in full plots */
  
  if (h->is_full)    		/* use all (raw) scores? */
    {
      esl_histogram_Sort(h);
      delta = h->w / 20.;

      lastx = DBL_MAX;
      for (i = h->n-1; i >= h->z; i--) /* sorted w/ low score at 0, high at n-1 */
	if (h->x[i] < lastx - delta)
	  {
	    fprintf(fp, "%f\t%g\n", h->x[i], (double)(h->n-i)/(double)h->Nc);
	    lastx = h->x[i];
	  }
    }
  else				/* else, use binned counts */
    {
      for (i = h->imax; i >= h->cmin; i--)
	{
	  c   += h->obs[i];
	  fprintf(fp, "%f\t%f\n", i*h->w + h->bmin, (double) c / (double) h->Nc);
	}
    }
  fprintf(fp, "&\n");
  return eslOK;
}

/* Function:  esl_histogram_PlotTheory()
 * Incept:    SRE, Thu Aug 25 08:09:24 2005 [St. Louis]
 *
 * Purpose:   Plot some theoretical distribution function <(*fx)()> 
 *            (a PDF, CDF, or survival function) that describes the
 *            data in histogram <h>, writing the plot in xmgrace XY
 *            format to open file <fp>. 
 *            
 *            The function is one of the generic API functions (such
 *            as <esl_gumbel_generic_pdf()>) that take a single <void
 *            *> argument to the distribution parameters, <params>,
 *            in addition to a quantile <x>. 
 *            
 *            The x axis (the quantile <x>) is varied from the minimum
 *            to the maximum of the observed data in the histogram.
 *            
 *            If the caller wants a wider range to be plotted (perhaps
 *            an extrapolated tail for larger <x>), it can use the
 *            appropriate plotting function for the specific distribution;
 *            these more specific plotting functions allow you to
 *            control the range of the x-axis.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotTheory(FILE *fp, ESL_HISTOGRAM *h, 
			 double (*fx)(double, void *), void *params)
{
  double x;
  double xmin, xmax, xstep;

  xmin  = (h->dataset_is == COMPLETE) ? h->xmin : h->phi;
  xmax  = h->xmax;
  xstep = h->w / 20.; 	/* plot points at 20x resolution of bin width */

  for (x = xmin; x <= xmax; x += xstep)
    printf("%f\t%f\n", x, (*fx)(x, params));
  printf("&\n");
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
 *            Respects any censoring information that's been set,
 *            or tail fitting that's been declared.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotQQ(FILE *fp, ESL_HISTOGRAM *h, 
		     double (*invcdf)(double x, void *params), void *params)
{
  int    i;
  double cdf;
  double bi;
  int    ibase;
  double sum = 0.;
  double lastx, delta; /* used to reduce the density of points in full plots */

  /* on censored data, fitted to a complete dist, start counting cdf at z,
   * not 0. 
   */
  if ((h->dataset_is == TRUE_CENSORED || h->dataset_is == VIRTUAL_CENSORED)
      && h->fit_describes == COMPLETE_FIT)	  
    sum = h->z; 

  if (h->is_full)    		/* use all (raw) scores? */
    {
      esl_histogram_Sort(h);
      delta = h->w / 20.;

      /* count empirical cdf only on 'observed' & fitted data:
       * so in virtual censored data, skip the first z samples.
       */
      if (h->dataset_is == VIRTUAL_CENSORED) ibase = h->z;
      else                                   ibase = 0;

      /* For each 'observed'/fitted data sample... bump the cdf & print a pt.
       */
      lastx = -DBL_MAX;		/* guarantee first delta test succeeds */
      for (i = ibase; i <= h->n-2; i++)   /* avoid last sample where cdf=1.0 */
	{
	  sum += 1.;
	  if (h->x[i] >= lastx + delta) { /* enforce some minimum spacing */
	    cdf = sum / (double) h->Nx;   /* to reduce the # of points    */
	    fprintf(fp, "%f\t%f\n", h->x[i], (*invcdf)(cdf, params));
	    lastx = h->x[i];
	  }
	}
    }
  else				/* else, use binned counts */
    {
      for (i = h->cmin; i < h->imax; i--) /* again, avoid last bin, cdf=1.0 */
	{
	  sum += (double) h->obs[i];
	  cdf = sum / (double) h->Nx;

	  bi = esl_histogram_Bin2UBound(h, i);
	  fprintf(fp, "%f\t%f\n", bi, (*invcdf)(cdf, params));
	}
    }
  fprintf(fp, "&\n");

  /* this plots a 45-degree expected QQ line:
   */
  if (h->dataset_is != COMPLETE) fprintf(fp, "%f\t%f\n", h->phi,  h->phi);
  else                           fprintf(fp, "%f\t%f\n", h->xmin, h->xmin);
  fprintf(fp, "%f\t%f\n", h->xmax, h->xmax);
  fprintf(fp, "&\n");

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
 *            The two tests should give similar
 *            probabilities. However, both tests are sensitive to
 *            arbitrary choices in how the data are binned, and
 *            neither seems to be on an entirely sound theoretical footing.
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
  int      ibase;

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
   *  
   * remember, h->No is the number of samples 'observed' & fitted.
   */
  nb   = 2* (int) pow((double) h->No, 0.4); /* "desired" nb. */
  minc = 1 + h->No / (2*nb);	/* arbitrarily set min = 1/2 of the target # */
  ESL_MALLOC(obs,  sizeof(int)    * (nb*2+1)); /* final nb must be <= 2*nb+1 */
  ESL_MALLOC(exp,  sizeof(double) * (nb*2+1));
  ESL_MALLOC(topx, sizeof(double) * (nb*2+1));

  /* Determine the observed counts in each bin: that is, partition h->No.
   * If we have raw counts, sort and use them (unless overridden by use_bindata)
   * If not, use the binned histogram. 
   * In either case, sweep left to right on the histogram bins,
   * collecting sum of counts, dropping the sum into the next bin 
   * whenever we have more than <minc> counts.
   * In the case of the raw counts, be careful that ties all go into
   * the same bin.
   * Also be careful to respect virtual censoring.
   */
  if (! use_bindata && h->is_full)
    {	/* collate raw counts */
      esl_histogram_Sort(h);

      /* Iterate over all observed counts:
       */
      if (h->dataset_is == VIRTUAL_CENSORED) ibase = h->z; 
      else                                   ibase = 0;    
      for (sum = 0, b = 0, i = ibase; i < h->n; i++) 
	{
	  sum++;
	  if (sum >= minc) { /* enough? then drop them, and all ties, in bin b */
	    {
	      ESL_DASSERT1( (b < nb*2+1) );
	      while (i < h->n-1 && h->x[i+1] == h->x[i]) { sum++; i++; } /* ties */
	      obs[b]  = sum;
	      topx[b] = h->x[i];
	      sum     = 0;
	      b++;
	    }
	  }
	}
      obs[b-1] += sum;		/* add the remaining right tail to the last bin */
      topx[b-1] = h->x[h->n-1]; /* by def'n */
      nb        = b;
    }
  else
    {	/* merge histogram bins */
      for (sum = 0, i = 0, b = h->cmin; b <= h->imax; b++) 
	{
	  sum += h->obs[b];
	  if (sum >= minc) {	/* if we have enough counts, drop them in bin i */
	    ESL_DASSERT1( (i < (nb*2+1)) );
	    obs[i]  = sum;
	    topx[i] = esl_histogram_Bin2UBound(h,b);
	    sum     = 0;
	    i++;
	  }
	}
      obs[i-1]  += sum;		/* add the right tail to our final bin        */
      topx[i-1]  = esl_histogram_Bin2UBound(h, h->imax);
      nb         = i;		/* nb is now the actual # of bins, not target */
    }


  /* Determine the expected counts in each bin.
   *  bin 0    is the left tail, <= topx[0];
   *  bin nb-1 is the right tail, > topx[nb-2],  1-P(<= topx[nb-2]).
   *  others are   P(<= topx[b]) - P(<topx[b-1]).
   */
  if (h->dataset_is == VIRTUAL_CENSORED || h->dataset_is == TRUE_CENSORED)
    exp[0] = (double) h->Nx * ( (*cdf)(topx[0], params) - (*cdf)(h->phi, params));
  else
    exp[0] = (double) h->Nx * (*cdf)(topx[0], params);

  for (i = 1; i < nb-1; i++)
    exp[i] = (double) h->Nx * ((*cdf)(topx[i], params) - (*cdf)(topx[i-1], params));

  exp[nb-1] = (double) h->Nx * (1 - (*cdf)(topx[nb-2], params));


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
  /* X^2 is distributed approximately chi^2.
   * If # obs = # expected, subtract an extra degree of freedom  */
  if (h->No == h->Nx) i = 1; else i = 0;
  if (nb-nfitted-i >= 0 && X2 != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted-i, X2, &X2p);
      if (status != eslOK) return status;
    }
  else X2p = 0.;


  /* Now, the G test assumes that #exp=#obs (the X^2 test didn't).
   * If that's not true, renormalize to make it so. Note that the
   * sum of exp[i] is not necessarily h->Nx, if we've fit a
   * complete distribution to a censored dataset; we actually have no
   * guarantees on what exp[i] might be in the fitted region.  We
   * do know that the total of obs[i] is h->No; we've always binned
   * all the 'observed' data.
   */
  if (h->No != h->Nx)
    {
      for (tmp = 0., i = 0; i < nb; i++) 
	tmp += exp[i];
      for (i = 0; i < nb; i++)
	exp[i] = exp[i] * (double) h->No / tmp;
    }

  for (tmp = 0., i = 0; i < nb; i++) 
    tmp += exp[i];
  
  /* Calculate the G statistic: 2 * LLR  */
  G = 0.;
  for (i = 0; i < nb; i++)
    {
      if (exp[i] == 0) G = eslINFINITY;
      else             G += (double) obs[i] * log ((double) obs[i] / exp[i]);
    }
  G *= 2;
  
  /* G is distributed approximately as \chi^2.
   * -1 is because total #obs=#exp (which is must be)
   */
  ESL_DASSERT1( (G >= 0.));
  if (nb-nfitted-1 >= 0 && G != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted-1, G, &Gp);
      if (status != eslOK) return status;
    }
  else Gp = 0.;

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
/* compile: gcc -g -Wall -I. -o example -DeslHISTOGRAM_EXAMPLE esl_histogram.c esl_random.c esl_stats.c easel.c -lm
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

  esl_histogram_Print(stdout, h);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example::end::*/
#endif /*eslHISTOGRAM_EXAMPLE*/




/*****************************************************************
 * Test driver code
 *****************************************************************/
#ifdef eslHISTOGRAM_TEST
/* compile: gcc -g -Wall -I. -L. -o test -DeslHISTOGRAM_TEST esl_histogram.c -leasel -lm
 * run:     
 *   ./test -1; ./test -2; ./test -3; ./test -4; ./test -5
 *   
 * Options for showing data for more detailed manual testing:
 *   -b      :  run goodness tests on binned data, not raw
 *   -n <n>  :  run <n> simulation trials, not just 1
 *   -p      :  dump empirical, expected histograms in xmgrace format to test.xy
 *   -P      :  print fancy ASCII histogram to stdout
 *   -q      :  dump QQ plot to test.xy
 *   -s      :  dump empirical, fitted survival plots to test.xy
 *   -v      :  verbose: print params, G tests, X^2 test to stdout
 *   -V      :  for tests -2 or -4: censor data by value, not by tail fraction
 *   
 * Some suggestions for manual testing:
 *   ./test -1 -n 5 -v -s; xmgrace test.xy          
 *        examine 5 survivor plot fits, for -1 
 *        do -2 thru -5 too
 *     
 *   ./test -1 -n 5 -v -q; xmgrace test.xy          
 *        examine 5 QQ plot fits, for -1 
 *        do -2 thru -5 too
 *        
 *   ./test -1 -n 1000 -v > foo
 *   grep "^Estimated" foo | awk '{print $9}' | sort -g > test.xy
 *        Look for straight line fit to G-test p values.
 *        sub $9->$13 for chi-squared
 *        sub Estimated -> Parametric for the parametric fits
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <easel.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_gumbel.h>
#include <esl_exponential.h>

static ESL_HISTOGRAM *
sim_complete_complete(ESL_RANDOMNESS *r, int nsamples, double *p, double *ep)
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-100, 100, 0.1);
  int    i;
  double x;
  
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, p[0], p[1]);
    esl_histogram_Add(h, x);
  }
  esl_gumbel_FitComplete(h->x, h->n, &(ep[0]), &(ep[1]));
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, ep);
  return h;
}
static ESL_HISTOGRAM *
sim_virtcensor_complete(ESL_RANDOMNESS *r, int nsamples, double *p, double *ep, 
			double vm, int by_value) /*v is either a mass or a value*/
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-100, 100, 0.1);
  int    i;
  double x;
  
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, p[0], p[1]);
    esl_histogram_Add(h, x);
  }
  if (by_value)
    esl_histogram_VirtCensorByValue(h, vm);
  else
    esl_histogram_VirtCensorByMass(h, vm);

  esl_gumbel_FitCensored(h->x+h->z, h->No, h->z, h->phi, &(ep[0]), &(ep[1]));
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, ep);
  return h;
}
static ESL_HISTOGRAM *
sim_truecensor_complete(ESL_RANDOMNESS *r, int nsamples, double *p, double *ep, double phi)
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-100, 100, 0.1);
  int    i;
  double x;
  int    z = 0;

  
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, p[0], p[1]);
    if (x > phi)
      esl_histogram_Add(h, x);
    else z++;
  }
  esl_histogram_TrueCensoring(h, z, phi);
  esl_gumbel_FitCensored(h->x, h->No, h->z, h->phi, &(ep[0]), &(ep[1]));
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, ep);
  return h;
}
static ESL_HISTOGRAM *
sim_virtcensor_tail(ESL_RANDOMNESS *r, int nsamples, double *p, double *ep, 
		    double vm, int by_value) /*v is either a mass or a value*/
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-100, 100, 0.1);
  int    i;
  double x;
  
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, p[0], p[1]);
    esl_histogram_Add(h, x);
  }
  if (by_value)
    esl_histogram_VirtCensorByValue(h, vm);
  else
    esl_histogram_VirtCensorByMass(h, vm);

  esl_histogram_SetTailfitting(h);

  ep[0] = h->phi;
  esl_exp_FitComplete(h->x+h->z, h->No, ep[0], &(ep[1]));
  esl_histogram_SetExpect(h, &esl_exp_generic_cdf, ep);
  return h;
}
static ESL_HISTOGRAM *
sim_truecensor_tail(ESL_RANDOMNESS *r, int nsamples, double *p, double *ep, double phi)
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-100, 100, 0.1);
  int    i;
  double x;
  int    z = 0;
  
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, p[0], p[1]);
    if (x > phi) 
      esl_histogram_Add(h, x);
    else z++;
  }
  esl_histogram_TrueCensoring(h, z, phi);
  esl_histogram_SetTailfitting(h);
  ep[0] = phi;
  esl_exp_FitComplete(h->x, h->No, ep[0], &(ep[1]));
  esl_histogram_SetExpect(h, &esl_exp_generic_cdf, ep);
  return h;
}
static int
binmacro_test(void)
{
  ESL_HISTOGRAM *h = esl_histogram_Create(-100, 100, 1.0);
  double trialx[3]  = { -42.42, 0, 42.42 };
  double x, ai, bi;  
  int    i,b;

  /* test bin<->score conversion macros.
   */
  for (i = 0; i < 3; i++)
    {
      x  = trialx[i];
      b  = esl_histogram_Score2Bin(h, x);
      ai = esl_histogram_Bin2LBound(h, b);
      bi = esl_histogram_Bin2UBound(h, b);
      if (x <= ai || x > bi) {
	fprintf(stderr,
		"failed: (ai=%.1f) <= (x=%.2f) < (bi=%.1f) in bin %d, bin macro test\n",
		ai, x, bi, b);
	esl_histogram_Destroy(h);
	return 0;
      }
    }
  esl_histogram_Destroy(h);
  return 1;
}
int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_Create(42);
  ESL_HISTOGRAM *h;
  double p[2] = { 10.0, 1.0 }; /* mu, lambda */
  double ep[2];                /* reestimated params */
  int    j;
  int    nb;
  double G, Gp, X2, X2p;
  double  (*cdf)(double, void *);
  double  (*invcdf)(double, void *);
  double avg_ep[2];
  int    nsamples = 10000;
  double minGp    = 1.;
  double minX2p   = 1.;
  double phi;
  double tfrac;

  int    optind       = 1;
  char  *outfile      = "test.xy";
  FILE  *fp;
  int    show_print   = FALSE;
  int    show_plot    = FALSE;
  int    show_surv    = FALSE;
  int    show_qq      = FALSE;
  int    bin_goodness = FALSE;
  int    by_value     = FALSE;
  int    verbose      = FALSE;
  int    ntrials      = 1;
  int    test_type    = 0;

  while (optind < argc && *argv[optind] == '-')
    {
      if      (strcmp(argv[optind], "-1")  == 0)  test_type    = 1;
      else if (strcmp(argv[optind], "-2")  == 0)  test_type    = 2;
      else if (strcmp(argv[optind], "-3")  == 0)  test_type    = 3;
      else if (strcmp(argv[optind], "-4")  == 0)  test_type    = 4;
      else if (strcmp(argv[optind], "-5")  == 0)  test_type    = 5;
      else if (strcmp(argv[optind], "-b")  == 0)  bin_goodness = TRUE;
      else if (strcmp(argv[optind], "-n")  == 0)  ntrials      = atoi(argv[++optind]);
      else if (strcmp(argv[optind], "-p")  == 0)  show_plot    = TRUE;
      else if (strcmp(argv[optind], "-P")  == 0)  show_print   = TRUE;
      else if (strcmp(argv[optind], "-q")  == 0)  show_qq      = TRUE;
      else if (strcmp(argv[optind], "-s")  == 0)  show_surv    = TRUE;
      else if (strcmp(argv[optind], "-v")  == 0)  verbose      = TRUE;
      else if (strcmp(argv[optind], "-V")  == 0)  by_value     = TRUE;
      optind++;
    }
  if (test_type == 0) {
    fprintf(stderr, "5 available test modes: choose option -1,-2,..,-5\n");
    exit(1);
  }
  avg_ep[0] = 0.;
  avg_ep[1] = 0.;

  fp = fopen(outfile, "w");

  for (j = 0; j < ntrials; j++)
    {
      switch (test_type) {
      case 1: 	/* complete dataset fitted to Gumbel */
	h   = sim_complete_complete(r, nsamples, p, ep);  
	cdf    = &esl_gumbel_generic_cdf;
	invcdf = &esl_gumbel_generic_invcdf;
	break;

      case 2:	/* virtually censored dataset, censored fit to Gumbel */
	phi      = 10.0;
        tfrac    = 0.5;
	if (by_value) h = sim_virtcensor_complete(r, nsamples, p, ep, phi,   TRUE);  
	else          h = sim_virtcensor_complete(r, nsamples, p, ep, tfrac, FALSE);  
	cdf    = &esl_gumbel_generic_cdf;
	invcdf = &esl_gumbel_generic_invcdf;
	break;

      case 3:	/* true censored dataset, censored fit to gumbel */
	phi    = 10.0;
	h = sim_truecensor_complete(r, nsamples, p, ep, phi);
	cdf    = &esl_gumbel_generic_cdf;
	invcdf = &esl_gumbel_generic_invcdf;
	break;

      case 4:	/* virtual censored dataset, tail fit to exponential */
	phi      = 12.5;
        tfrac    = 0.1;
	if (by_value) h = sim_virtcensor_tail(r, nsamples, p, ep, phi,   TRUE);
	else          h = sim_virtcensor_tail(r, nsamples, p, ep, tfrac, FALSE);
	cdf    = &esl_exp_generic_cdf;
	invcdf = &esl_exp_generic_invcdf;
	break;

      case 5:	/* true censored dataset, tail fit to exponential */
	phi      = 12.5;
	h = sim_truecensor_tail(r, nsamples, p, ep, phi);
	cdf    = &esl_exp_generic_cdf;
	invcdf = &esl_exp_generic_invcdf;
	break;
      }
	
      avg_ep[0] += ep[0];
      avg_ep[1] += ep[1];

      if (show_print) esl_histogram_Print(stdout, h);

      /* parametric is always Gumbel
       */
      if (test_type <= 3)
	{
	  esl_histogram_Goodness(h, esl_gumbel_generic_cdf, p, 
			     0, bin_goodness, &nb, &G, &Gp, &X2, &X2p);
	  if (Gp  < minGp)  minGp  = Gp;
	  if (X2p < minX2p) minX2p = X2p;
	  if (verbose) 
	    printf("Parametric: %6.2f %6.4f nb %4d G %g\tGp %g\tX2 %g\tX2p %g\n",
		   p[0], p[1], nb, G, Gp, X2, X2p);
	}

      /* fitted may be Gumbel or exponential; use "cdf" ptr
       */
      esl_histogram_Goodness(h, cdf, ep, 
			     2, bin_goodness, &nb, &G, &Gp, &X2, &X2p);
      if (Gp  < minGp)  minGp  = Gp;
      if (X2p < minX2p) minX2p = X2p;
      if (verbose)
	printf("Estimated:  %6.2f %6.4f nb %4d G %g\tGp %g\tX2 %g\tX2p %g\n",
	       ep[0], ep[1], nb, G, Gp, X2, X2p);

      if (show_plot) esl_histogram_Plot(fp, h);
      if (show_qq)   esl_histogram_PlotQQ(fp, h, invcdf, ep);
      
      if (show_surv) {
	esl_histogram_PlotSurvival(fp, h);
	if (test_type <= 3)
	  esl_gumbel_Plot(fp, ep[0], ep[1], &esl_gumbel_surv, h->xmin-5, h->xmax+5, 0.1);
	else
	  esl_exp_Plot(fp, ep[0], ep[1], &esl_exp_surv, h->phi, h->xmax+5, 0.1);
      }

      esl_histogram_Destroy(h);
    }
  fclose(fp);

  avg_ep[0] /= (double) ntrials;
  avg_ep[1] /= (double) ntrials;

  /* Trap bad fits.
   */
  if (test_type <= 3 && fabs(avg_ep[0] - p[0]) > 0.1)
    ESL_ERROR(eslFAIL, "Something awry with Gumbel mu fit");
  if (fabs(avg_ep[1] - p[1]) > 0.1)
    ESL_ERROR(eslFAIL, "Something awry with lambda fit");

  if (minGp < 1. / (1000. * ntrials))
    ESL_ERROR(eslFAIL, "Something awry with G-test");
  if (minX2p < 1. / (1000. * ntrials))
    ESL_ERROR(eslFAIL, "Something awry with chi squared test");

  /* Smaller final tests
   */
  if (! binmacro_test()) exit(0);
  esl_randomness_Destroy(r);
  return 0;
}
#endif /*eslHISTOGRAM_TEST*/

