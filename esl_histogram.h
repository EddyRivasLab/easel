/* esl_histogram.h
 * Collection and display of score histograms.
 * 
 * SRE, Fri Jul  1 13:22:45 2005 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_HISTOGRAM_INCLUDED
#define ESL_HISTOGRAM_INCLUDED


/* Structure: ESL_HISTOGRAM
 * 
 * Keeps a score histogram, in which scores are counted into bins of
 * size (width) w. 
 *   histogram starts at bmin <  floor(xmin/w) * w
 *   histogram ends at   bmax >= ceil(xmax/w)*w
 *   nb = (bmax-bmin)/w
 *   each score x is counted into bin b = nb - (int) (bmax-x)/w
 *   each bin b contains scores bw+bmin < x <= (b+1)w + bmin
 *
 */  
typedef struct {
  /* The histogram is kept as counts in fixed-width bins.
   */
  double  xmin, xmax;	/* smallest, largest sample value observed          */
  int     n;            /* total number of raw data samples                 */
  int    *obs;		/* observed counts in bin b, 0..nb-1 (dynamic)      */
  double *expect;	/* expected counts in bin b, 0..nb-1 (not resized)  */
  double  bmin, bmax;	/* histogram bounds: all x satisfy bmin < x <= bmax */
  int     imin, imax;	/* smallest, largest bin that contain obs[i] > 0    */
  int     nb;           /* number of bins                                   */
  double  w;		/* fixed width of each bin                          */

  /* Optionally, in a "full" h, we can also keep all the raw samples in x.
   */
  double *x;		/* optional: raw sample values x[0..n-1]            */
  int     nalloc;	/* current allocated size of x                      */

  /* Some status flags
   */
  int is_full;		/* TRUE when we're keeping raw data in x   */
  int is_sorted;	/* TRUE if x is sorted smallest-to-largest */
} ESL_HISTOGRAM;


extern ESL_HISTOGRAM *esl_histogram_Create    (double bmin, double bmax, double w);
extern ESL_HISTOGRAM *esl_histogram_CreateFull(double bmin, double bmax, double w);
extern void           esl_histogram_Destroy(ESL_HISTOGRAM *h);

extern int esl_histogram_Add(ESL_HISTOGRAM *h, double x);

extern int esl_histogram_Sort(ESL_HISTOGRAM *h);
extern int esl_histogram_GetScoreAtRank(ESL_HISTOGRAM *h, int rank, double *ret_x);
extern int esl_histogram_GetBinBounds(ESL_HISTOGRAM *h, int whichbin,
				      double *ret_low, double *ret_high,
				      double *ret_delta);

extern int esl_histogram_Print       (FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_Plot        (FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_PlotQQ      (FILE *fp, ESL_HISTOGRAM *h, 
				      double (*invcdf)(double x, void *params),
				      void *params);

extern int esl_histogram_SetExpect(ESL_HISTOGRAM *h, 
				   double (*cdf)(double x, void *params),
				   void *params);

#ifdef eslAUGMENT_STATS
extern int esl_histogram_Goodness(ESL_HISTOGRAM *h, 
				  double (*cdf)(double x, void *params),
				  void *params, int nfitted, int use_bindata,
				  int *ret_nbins,
				  double *ret_G,  double *ret_Gp,
				  double *ret_X2, double *ret_X2p);
#endif



#endif /*!ESL_HISTOGRAM_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
