/* esl_histogram.c
 * Collecting and displaying histograms.
 *
 * SRE, Fri Jul  1 13:21:45 2005 [St. Louis]
 * SVN $Id$
 */

/*
 * Creating/destroying histograms and collecting data:
 *  ESL_HISTOGRAM *esl_histogram_Create(double xmin, double xmax, double w)
 *  ESL_HISTOGRAM *esl_histogram_CreateFull(double xmin, double xmax, double w)
 *  void           esl_histogram_Destroy(ESL_HISTOGRAM *h)
 *  int            esl_histogram_Add(ESL_HISTOGRAM *h, double x)
 *  int            esl_histogram_Sort(ESL_HISTOGRAM *h)
 *
 * Accessing data samples in a full histogram:
 *  int esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank, double *ret_x)
 *  int esl_histogram_GetData(ESL_HISTOGRAM *h, double **ret_x, int *ret_n)
 *  int esl_histogram_GetTail(ESL_HISTOGRAM *h, double phi, 
 *                            double **ret_x, int *ret_n, int *ret_z)
 *  int esl_histogram_GetTailByMass(ESL_HISTOGRAM *h, double pmass, 
 *                            double **ret_x, int *ret_n, int *ret_z)
 *    
 * Declarations about the binned data before parameter fitting:
 *  int esl_histogram_SetTrueCensoring(ESL_HISTOGRAM *h, int z, double phi)
 *  int esl_histogram_VirtCensor(ESL_HISTOGRAM *h, double phi)
 *  int esl_histogram_VirtCensorByMass(ESL_HISTOGRAM *h, double pmass)
 *  
 * Setting expected binned counts:
 *  int esl_histogram_SetExpect(ESL_HISTOGRAM *h, 
 *	  	  double (*cdf)(double x, void *params), void *params)
 *  int esl_histogram_SetExpectedTail(ESL_HISTOGRAM *h, double base_val, 
 *                double pmass, double (*cdf)(double x, void *params), 
 *  	          void *params)
 *  
 * Output/display of binned data:
 *   int esl_histogram_Print(FILE *fp, ESL_HISTOGRAM *h) 
 *   int esl_histogram_Plot(FILE *fp, ESL_HISTOGRAM *h)
 *   int esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h)
 *   int esl_histogram_PlotQQ(FILE *fp, ESL_HISTOGRAM *h, 
 *		     double (*invcdf)(double x, void *params), void *params)
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <easel.h>
#include <esl_histogram.h>
#include <esl_vectorops.h>

#ifdef eslAUGMENT_STATS	 /* stats augmentation gives goodness-of-fit testing */
#include <esl_stats.h>
#endif


/*****************************************************************
 * Creating/destroying histograms and collecting data.
 *****************************************************************/

/* Function:  esl_histogram_Create()
 * Incept:    SRE, Fri Jul  1 13:40:26 2005 [St. Louis]
 *
 * Purpose:   Creates and returns a new histogram object, initially
 *            allocated to count scores $>$ <xmin> and $<=$ <xmax> into
 *            bins of width <w>. Thus, a total of <xmax>-<xmin>/<w> bins
 *            are initially created. 
 *            
 *            The bounds <xmin> and <xmax> only need to be initial
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
 * Args:      xmin - caller guesses that minimum score will be > xmin
 *            xmax - caller guesses that max score will be <= xmax
 *            w    - size of bins (1.0, for example)
 *            
 * Returns:   ptr to new <ESL_HISTOGRAM> object, which caller is responsible
 *            for free'ing with <esl_histogram_Destroy()>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_HISTOGRAM *
esl_histogram_Create(double xmin, double xmax, double w)
{
  ESL_HISTOGRAM *h = NULL;
  int i;

  if ((h = malloc(sizeof(ESL_HISTOGRAM))) == NULL)    
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  h->xmin      =  DBL_MAX;	/* xmin/xmax are the observed min/max */
  h->xmax      = -DBL_MAX;
  h->n         = 0;
  h->obs       = NULL;		/* briefly... */
  h->bmin      = xmin;		/* bmin/bmax are the allocated bounds */
  h->bmax      = xmax;
  h->nb        = (int)((xmax-xmin)/w);
  h->imin      = h->nb;
  h->imax      = -1;
  h->w         = w;

  h->x         = NULL;
  h->nalloc    = 0;

  h->phi       = 0.;
  h->cmin      = h->imin;	/* sentinel: no observed data yet */
  h->z         = 0;
  h->Nc        = 0;
  h->No        = 0;

  h->expect    = NULL;		/* 'til a Set*() call */
  h->emin      = -1;            /* sentinel: no expected counts yet */
  h->tailbase  = 0.;		/* unused unless is_tailfit TRUE */
  h->tailmass  = 1.0;		/* <= 1.0 if is_tailfit TRUE */

  h->is_full       = FALSE;
  h->is_done       = FALSE;
  h->is_sorted     = FALSE;
  h->is_tailfit    = FALSE;
  h->dataset_is    = COMPLETE;

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
 *            tests, but will consume more memory. 
 */
ESL_HISTOGRAM *
esl_histogram_CreateFull(double xmin, double xmax, double w)
{
  ESL_HISTOGRAM *h = esl_histogram_Create(xmin, xmax, w);
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

  /* Censoring info must only be set on a finished histogram;
   * don't allow caller to add data after configuration has been declared
   */
  if (h->is_done)
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

  /* Bump the bin counter, and all the data sample counters.
   */
  h->obs[b]++;
  h->n++;
  h->Nc++;
  h->No++;

  if (b > h->imax) h->imax = b;
  if (b < h->imin) { h->imin = b; h->cmin = b; }
  if (x > h->xmax) h->xmax = x;
  if (x < h->xmin) h->xmin = x;
  return eslOK;
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
  if (h->is_sorted) return eslOK; /* already sorted, don't do anything */
  if (! h->is_full) return eslOK; /* nothing to sort */
  if (h->is_done) ESL_ERROR(eslEINCONCEIVABLE, "unsorted data can't be done");
  
  esl_vec_DSortIncreasing(h->x, h->n);
  h->is_sorted = TRUE;
  return eslOK;
}


/*****************************************************************
 * Routines for accessing data samples in a full histogram.
 *****************************************************************/

/* Function:  esl_histogram_GetRank()
 * Incept:    SRE, Thu Jul 28 08:39:52 2005 [St. Louis]
 *
 * Purpose:   Retrieve the <rank>'th highest score from a 
 *            full histogram <h>. <rank> is <1..n>, for
 *            <n> total samples in the histogram; return it through
 *            <ret_x>.
 *            
 *            If the raw scores aren't sorted, they are sorted
 *            first (an NlogN operation).
 *            
 *            This can be called at any time, even during data
 *            collection, to see the current <rank>'th highest score.
 *
 * Throws:    <eslEINVAL> if the histogram is display-only,
 *            or if <rank> isn't in the range 1..n.
 */
int
esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank, double *ret_x)
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

/* Function:  esl_histogram_GetData()
 * Incept:    SRE, Fri Jan 27 07:57:21 2006 [St. Louis]
 *
 * Purpose:   Retrieve the raw data values from the histogram <h>.
 *            Return them in the vector <ret_x>, and the number
 *            of values in <ret_n>. The values are indexed <[0..n-1]>,
 *            from smallest to largest (<x[n-1]> is the high score).
 *            
 *            <ret_x> is a pointer to internal memory in the histogram <h>.
 *            The histogram <h> is still responsible for that storage;
 *            its memory will be free'd when you call
 *            <esl_histogram_Destroy()>.
 *            
 *            You can only call this after you have finished collecting
 *            all the data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.
 *            
 * Internal note:
 *            The prohibition against adding more data (by raising
 *            the h->is_done flag) is because we're passing a pointer
 *            to internal data storage back to the caller. Subsequent
 *            calls to Add() will modify that memory -- in the worst case,
 *            if Add() has to reallocate that storage, completely invalidating
 *            the pointer that the caller has a copy of. We want to make
 *            sure that the <ret_x> pointer stays valid.
 *            
 * Args:      h     - histogram to retrieve data values from
 *            ret_x - RETURN: pointer to the data samples, [0..n-1] 
 *            ret_n - RETURN: number of data samples
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram <h> is not a full histogram.
 */
int
esl_histogram_GetData(ESL_HISTOGRAM *h, double **ret_x, int *ret_n)
{
  if (! h->is_full) ESL_ERROR(eslEINVAL, "not a full histogram");
  esl_histogram_Sort(h);

  *ret_x = h->x;
  *ret_n = h->n;

  h->is_done = TRUE;
  return eslOK;
}


/* Function:  esl_histogram_GetTail()
 * Incept:    SRE, Fri Jan 27 07:56:38 2006 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, retrieve all data values 
 *            above the threshold <phi> in the right (high scoring) 
 *            tail, as a ptr <ret_x> to an array of <ret_n> values 
 *            indexed <[0..n-1]> from lowest to highest score. 
 *            Optionally, it also returns the number of values in 
 *            rest of the histogram in <ret_z);
 *            this number is useful if you are going to fit
 *            the tail as a left-censored distribution.
 *            
 *            The test is strictly greater than <phi>, not greater
 *            than or equal to.
 *            
 *            <ret_x> is a pointer to internal memory in the histogram <h>.
 *            The histogram <h> is still responsible for that storage;
 *            its memory will be free'd when you call 
 *            <esl_histogram_Destroy()>.
 *            
 *            You can only call this after you have finished collecting
 *            all the data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.             
 *            
 * Args:      h     - histogram to retrieve the tail from
 *            phi   - threshold: tail is all scores > phi
 *            ret_x - optRETURN: ptr to vector of data values [0..n-1]
 *            ret_n - optRETURN: number of data values in tail
 *            ret_z - optRETURN: number of data values not in tail.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram is not a full histogram.
 */
int
esl_histogram_GetTail(ESL_HISTOGRAM *h, double phi, 
		      double **ret_x, int *ret_n, int *ret_z)
{
  int hi, lo, mid;

  if (! h->is_full) ESL_ERROR(eslEINVAL, "not a full histogram");
  esl_histogram_Sort(h);

  if      (h->n         == 0)   mid = h->n;  /* we'll return NULL, 0, n */  
  else if (h->x[0]       > phi) mid = 0;     /* we'll return x, n, 0    */
  else if (h->x[h->n-1] <= phi) mid = h->n;  /* we'll return NULL, 0, n */
  else /* binary search, faster than a brute force scan */
    {
      lo = 0;
      hi = h->n-1; /* know hi>0, because above took care of n=0 and n=1 cases */
      while (1) {
	mid = (lo + hi + 1) / 2;  /* +1 makes mid round up, mid=0 impossible */
	if      (h->x[mid]  <= phi) lo = mid; /* we're too far left  */
	else if (h->x[mid-1] > phi) hi = mid; /* we're too far right */
	else break;		              /* ta-da! */
      }
    }

  if (ret_x != NULL) *ret_x = h->x + mid;
  if (ret_n != NULL) *ret_n = h->n - mid;
  if (ret_z != NULL) *ret_z = mid;
  h->is_done = TRUE;
  return eslOK;
}


/* Function:  esl_histogram_GetTailByMass()
 * Incept:    SRE, Sun Jan 29 17:56:37 2006 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, retrieve the data values in
 *            the right (high scoring) tail, as a pointer <ret_x>
 *            to an array of <ret_n> values indexed <[0..n-1]> from
 *            lowest to highest score. The tail is defined by a
 *            given mass fraction threshold <pmass>; the mass in the returned
 *            tail is $\leq$ this threshold. <pmass> is a probability,
 *            so it must be $\geq 0$ and $\leq 1$.
 *            
 *            Optionally, the number of values in the rest of the
 *            histogram can be returned in <ret_z>. This is useful
 *            if you are going to fit the tail as a left-censored
 *            distribution.
 *            
 *            <ret_x> is a pointer to internal memory in <h>. 
 *            The histogram <h> remains responsible for its storage,
 *            which will be free'd when you call <esl_histogram_Destroy()>.
 *            As a consequence, you can only call 
 *            <esl_histogram_GetTailByMass()> after you have finished
 *            collecting data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.
 *
 * Args:      h     - histogram to retrieve the tail from
 *            pmass - fractional mass threshold; tail contains <= pmass
 *            ret_x - optRETURN: ptr to vector of data values [0..n-1]
 *            ret_n - optRETURN: number of data values in tail x
 *            ret_z - optRETURN: number of data values not in tail
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram is not a full histogram, 
 *            or <pmass> is not a probability.
 */
int
esl_histogram_GetTailByMass(ESL_HISTOGRAM *h, double pmass,
			    double **ret_x, int *ret_n, int *ret_z)
{
  int n;

  if (! h->is_full) 
    ESL_ERROR(eslEINVAL, "not a full histogram");
  if (pmass < 0. || pmass > 1.) 
    ESL_ERROR(eslEINVAL, "pmass not a probability");

  esl_histogram_Sort(h);

  n = (int) ((float) h->n * pmass); /* rounds down, guaranteeing <= pmass */

  if (ret_x != NULL) *ret_x = h->x + (h->n - n);
  if (ret_n != NULL) *ret_n = n;
  if (ret_z != NULL) *ret_z = h->n - z;
  h->is_done = TRUE;
  return eslOK;
}




/*****************************************************************
 * Declarations about the binned data before parameter fitting
 *****************************************************************/ 

/* Function:  esl_histogram_SetTrueCensoring()
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
 *            $x \leq phi$. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if you try to set <phi> to a value that is
 *            greater than the minimum <x> stored in the histogram.
 */
int
esl_histogram_SetTrueCensoring(ESL_HISTOGRAM *h, int z, double phi)
{
  if (phi > h->xmin) ESL_ERROR(eslEINVAL, "no uncensored x can be <= phi");

  h->phi         = phi;
  h->cmin        = h->imin;
  h->z           = z;
  h->Nc          = h->n + z;
  h->No          = h->n;
  h->dataset_is  = TRUE_CENSORED;
  h->is_done     = TRUE;
  return eslOK;
}

/* Function:  esl_histogram_VirtCensor()
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
 *            in the censored region for purposes of parameter
 *            fitting, calculating expected binned counts,
 *            and binned goodness-of-fit tests. 
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
esl_histogram_VirtCensor(ESL_HISTOGRAM *h, double phi)
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
  h->dataset_is = VIRTUAL_CENSORED;
  h->is_done    = TRUE;

  esl_histogram_Sort(h); /* uncensored raw tail now starts at h->x+h->z */
  return eslOK;
}

/* Function:  esl_histogram_VirtCensorByMass()
 * Incept:    SRE, Tue Aug 23 08:10:39 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> (with or without raw data samples),
 *            find a cutoff score that at least fraction <pmass> of the samples
 *            exceed. This threshold is stored internally in the histogram
 *            as <h->phi>. The number of virtually censored samples (to the 
 *            left, with scores $\leq \phi$) is stored internally in <h->z>.
 *            
 *            The identified cutoff score must be a lower bound for some bin
 *            (bins can't be partially censored). The censored mass
 *            will thus usually be a bit greater than <pmass>, as the
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
 *            After calling <esl_histogram_VirtCensorByMass()> on a
 *            full histogram, the caller can retrieve the sorted
 *            censored data as <h->x+h->z>, which is a (partial)
 *            vector containing <h->Nc-h->z> numbers, all satisfying
 *            $x_i > \phi$.
 *            
 *            Additionally, after calling <esl_histogram_VirtCensorByMass()>
 *            on a histogram, the index of the first uncensored bin is in
 *            <h->cmin>. That is, the censored bins are <0..h->cmin-1>
 *            and the uncensored bins are <h->cmin..nb-1>; or
 *            alternatively, for the range of bins that contain
 *            counts, <h->imin..h->cmin-1> are censored and
 *            <h->cmin..h->imax>.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_VirtCensorByMass(ESL_HISTOGRAM *h, double pmass)
{
  int b;
  int sum = 0;
	    
  for (b = h->imax; b >= h->imin; b--)
    {
      sum += h->obs[b];
      if (sum >= (pmass * (double)h->n)) break;
    }

  h->phi         = esl_histogram_Bin2LBound(h,b);
  h->z           = h->n - sum;
  h->cmin        = b;
  h->Nc          = h->n;	/* (redundant) */
  h->No          = h->n - h->z;
  h->dataset_is  = VIRTUAL_CENSORED;
  h->is_done     = TRUE;

  esl_histogram_Sort(h); /* uncensored raw tail now starts at h->x+h->z */
  return eslOK;
}



/*****************************************************************
 * Setting expected counts
 *****************************************************************/ 

/* Function:  esl_histogram_SetExpect()
 * Incept:    SRE, Wed Aug 17 17:36:58 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> containing some number of empirically
 *            observed binned counts, and a pointer to a function <(*cdf)()>
 *            that describes the expected cumulative distribution function 
 *            (CDF) for the complete data, conditional on some parameters 
 *            <params>; calculate the expected counts in each bin of the 
 *            histogram, and hold that information internally in the structure.
 *            
 *            The caller provides a function <(*cdf)()> that calculates
 *            the CDF via a generic interface, taking only two
 *            arguments: a quantile <x> and a void pointer to whatever
 *            parameters it needs, which it will cast and interpret.
 *            The <params> void pointer to the given parameters is
 *            just passed along to the generic <(*cdf)()> function. The
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
  double ai,bi;			/* ai < x <= bi : lower,upper bounds in bin */

  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      ai = esl_histogram_Bin2LBound(h, i);
      bi = esl_histogram_Bin2UBound(h, i);
      h->expect[i] = h->Nc * ( (*cdf)(bi, params) - (*cdf)(ai, params) );

      if (h->emin != -1 && h->expect[i] > 0.) h->emin = i;
    }

  h->is_done = TRUE;
  return eslOK;
}

/* Function:  esl_histogram_SetExpectedTail()
 * Incept:    SRE, Mon Jan 30 08:57:57 2006 [St. Louis]
 *
 * Purpose:   Given a histogram <h>, and a pointer to a generic function
 *            <(*cdf)()> that describes the expected cumulative
 *            distribution function for the right (high-scoring) tail
 *            starting at <base_val> (all expected <x> $>$ <base_val>) and
 *            containing a fraction <pmass> of the complete data
 *            distribution (<pmass> $\geq 0$ and $\leq 1);
 *            set the expected binned counts for all complete bins
 *            $\geq$ <base_val>. 
 *            
 *            If <base_val> falls within a bin, that bin is considered
 *            to be incomplete, and the next higher bin is the starting
 *            point. 
 *           
 * Args:      h          - finished histogram
 *            base_val   - threshold for the tail: all expected x > base_val
 *            pmass      - fractional mass in the tail: 0 <= pmass <= 1
 *            cdf        - generic-interface CDF function describing the tail
 *            params     - void pointer to parameters for (*cdf)()
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
esl_histogram_SetExpectedTail(ESL_HISTOGRAM *h, double base_val, double pmass,
			      double (*cdf)(double x, void *params), 
			      void *params)
{
  int b;
  double ai, bi;

  h->emin = 1 + esl_histogram_Score2Bin(h, base_val);

  if (h->expect == NULL) 
    ESL_MALLOC(h->expect, sizeof(double) * h->nb);
  for (b = 0; b < h->emin; b++) 
    h->expect[b] = 0.;

  for (b = h->emin; b < h->nb; b++)
    {
      ai = esl_histogram_Bin2LBound(h, b);
      bi = esl_histogram_Bin2UBound(h, b);
      h->expect[b] = pmass * (double) h->Nc * 
	             ( (*cdf)(bi, params) - (*cdf)(ai, params) );
    }
  
  h->tailbase   = base_val;
  h->tailmass   = pmass;
  h->is_tailfit = TRUE;
  h->is_done    = TRUE;
  return eslOK;
}




/*****************************************************************
 * Output and display of binned data.
 *****************************************************************/ 

/* Function:  esl_histogram_Print() 
 * Incept:    SRE, Sat Jul  2 16:03:37 2005 [St. Louis]
 *
 * Purpose:   Print a "prettified" display histogram <h> to a file 
 *            pointer <fp>.
 *            Deliberately a look-and-feel clone of Bill Pearson's 
 *            excellent FASTA output.
 *            
 *            Also displays expected binned counts, if they've been
 *            set.
 *            
 *            Display will only work well if the bin width (w) is 0.1 or more,
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
 * Incept:    SRE, Mon Jan 30 11:09:01 2006 [St. Louis]
 *
 * Purpose:   Print observed (and expected, if set) counts
 *            in a histogram <h> to open file pointer <fp>
 *            in xmgrace XY input file format.
 *
 * Returns:   <eslOK> on success.
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
	if (h->expect[i] > 0.)	/* >0 suffices to remove censored region */
	  {
	    x = esl_histogram_Bin2LBound(h,i);
	    fprintf(fp, "%.2f %g\n", x, h->expect[i]);
	  }
      fprintf(fp, "&\n");
    }
  return eslOK;
}

/* Function:  esl_histogram_PlotSurvival()
 * Incept:    SRE, Mon Jan 30 11:11:05 2006 [St. Louis]
 *
 * Purpose:   Given a histogram <h>, output the observed (and
 *            expected, if available) survival function $P(X>x)$
 *            to file pointer <fp> in xmgrace XY input file format.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h)
{
  int i;
  double c;
  double ai;
  
  /* The observed binned counts:
   */
  c = 0.;
  for (i = h->imax; i >= h->imin; i--)
    {
      if (h->obs[i] > 0) {
	c   += h->obs[i];
	ai = esl_histogram_Bin2LBound(h, i);
	fprintf(fp, "%f\t%f\n", ai, c / (double) h->Nc);
      }
    }
  fprintf(fp, "&\n");

  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      c = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	  if (h->expect[i] > 0.) {
	    c += h->expect[i];
	    ai = esl_histogram_Bin2LBound(h, i);
	    fprintf(fp, "%f\t%f\n", ai, c / (double) h->Nc);
	  }
	}
      fprintf(fp, "&\n");
    }
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
 *            Same domain limits as goodness-of-fit testing: output
 *            is restricted to overlap between observed data (excluding
 *            any censored data) and expected data (which may be limited
 *            if only a tail was fit).
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
  double sum;

  /* on censored data, start counting observed cdf at z, not 0
   */
  if (h->dataset_is == TRUE_CENSORED || h->dataset_is == VIRTUAL_CENSORED)
    sum = h->z; 
  else
    sum = 0.;

  /* Determine smallest bin included in goodness of fit eval
   */
  bbase = h->cmin;
  if (h->is_tailfit && h->emin > bbase) bbase = h->emin;
  for (i = h->cmin; i < bbase; i++) sum += (double) h->obs[i];
  
  /* The q-q plot:
   */
  for (i = bbase; i < h->imax; i++) /* avoid last bin where upper cdf=1.0 */
    {
      sum += (double) h->obs[i];
      cdf = sum / (double) h->Nc;

      if (h->is_tailfit) cdf = (cdf + h->tailmass - 1.) / (h->tailmass);

      bi = esl_histogram_Bin2UBound(h, i);
      fprintf(fp, "%f\t%f\n", bi, (*invcdf)(cdf, params));
    }
  fprintf(fp, "&\n");

  /* Plot a 45-degree expected QQ line:
   */
  bi = esl_histogram_Bin2LBound(h, bbase);
  fprintf(fp, "%f\t%f\n", bi,  bi);
  fprintf(fp, "%f\t%f\n", h->xmax, h->xmax);
  fprintf(fp, "&\n");

  return eslOK;
}




#ifdef eslAUGMENT_STATS
/* Function:  esl_histogram_Goodness()
 * Incept:    SRE, Wed Aug 17 12:46:05 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> with observed and expected counts,
 *            where, for the expected counts, <nfitted> ($\geq 0$)
 *            parameters were fitted (and thus should be subtracted
 *            from the degrees of freedom);
 *            Perform a G-test and/or a chi-squared test for goodness of 
 *            fit between observed and expected, and optionally return
 *            the number of bins the data were sorted into
 *            (<ret_bins>), the G statistic and its probability (<ret_G> and
 *            <ret_Gp>), and the X$^2$ statistic and its probability
 *            (<ret_X2> and <ret_X2p>). 
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
 *            errors in calculating the probability;
 *            <eslEMEM> on internal allocation failure.
 */
int
esl_histogram_Goodness(ESL_HISTOGRAM *h, 
		       int nfitted, int *ret_nbins,
		       double *ret_G,  double *ret_Gp,
		       double *ret_X2, double *ret_X2p)
{
  int     *obs;			/* observed in bin i, [0..nb-1]   */
  double  *exp;			/* expected in bin i, [0..nb-1]   */
  double  *topx;		/* all values in bin i <= topx[i] */
  int      nb;			/* # of re-bins                   */
  int      minc;		/* minimum target # of counts/bin */
  int      i,b;
  double   G, Gp;
  double   X2, X2p;
  double   tmp;
  int      status;
  int      bbase;
  int      hmax;
  int      nobs;
  double   nexp;

  /* Determine the smallest histogram bin included in 
   * the goodness of fit evaluation.
   */
  bbase = h->cmin;		
  if (h->is_tailfit && h->emin > bbase) bbase = h->emin;
  
  /* How many observed total counts are in the evaluated range,
   * and what is the maximum in any given histogram bin?
   */
  nobs = 0;
  hmax = 0;
  for (i = bbase; i <= h->imax; i++)
    {
      nobs += h->obs[i];
      if (h->obs[i] > hmax) hmax = h->obs[i];
    }

  /* Figure out how many eval bins we'd like to have, then allocate
   * for re-binning.
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
   * The most important thing seems to be to get the # of counts
   * in each bin to be roughly equal.
   */
  nb   = 2* (int) pow((double) nobs, 0.4); /* "desired" nb. */
  minc = 1 + nobs / (2*nb);	/* arbitrarily set min = 1/2 of the target # */
  ESL_MALLOC(obs,  sizeof(int)    * (nb*2+1)); /* final nb must be <= 2*nb+1 */
  ESL_MALLOC(exp,  sizeof(double) * (nb*2+1));
  ESL_MALLOC(topx, sizeof(double) * (nb*2+1));

  /* Determine the observed counts in each bin: that is, partition 
   * the <sum> in the evaluated region.
   * Sweep left to right on the histogram bins,
   * collecting sum of counts, dropping the sum into the next re-bin 
   * whenever we have more than <minc> counts.
   */
  nobs = nexp = 0;
  for (i = 0, b = bbase; b <= h->imax; b++) 
    {
      nobs += h->obs[b];
      nexp += h->expect[b];

      /* if we have enough counts, drop into bin i: */
      if (nobs >= minc && nexp >= minc) {
	ESL_DASSERT1( (i < (nb*2+1)) );
	obs[i]  = nobs;
	exp[i]  = nexp;
	topx[i] = esl_histogram_Bin2UBound(h,b);
	nobs = nexp = 0;
	i++;
      }
    }
  obs[i-1]  += nobs;		/* add the right tail to final bin */
  exp[i-1]  += nexp;
  topx[i-1]  = esl_histogram_Bin2UBound(h, h->imax);
  nb         = i;		/* nb is now actual # of bins, not target */

  /* Calculate the X^2 statistic: \sum (obs_i - exp_i)^2 / exp_i */
  X2 = 0.;
  for (i = 0; i < nb; i++)
    {
      tmp = obs[i] - exp[i];
      X2 += tmp*tmp / exp[i];
    }
  /* X^2 is distributed approximately chi^2. */
  if (nb-nfitted >= 0 && X2 != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted, X2, &X2p);
      if (status != eslOK) return status;
    }
  else X2p = 0.;

  /* The G test assumes that #exp=#obs (the X^2 test didn't).
   * If that's not true, renormalize to make it so. 
   */
  nobs = nexp = 0;
  for (i = 0; i < nb; i++) 
    {
      nobs += obs[i];
      nexp += exp[i];
    }
  for (i = 0; i < nb; i++)
    exp[i] = exp[i] * (double) nobs / nexp;
  
  /* Calculate the G statistic: 2 * LLR  */
  G = 0.;
  for (i = 0; i < nb; i++)
    G += (double) obs[i] * log ((double) obs[i] / exp[i]);
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

  esl_histogram_DeclareTailfitting(h);

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
  esl_histogram_DeclareTailfitting(h);
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
			     0, &nb, &G, &Gp, &X2, &X2p);
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

