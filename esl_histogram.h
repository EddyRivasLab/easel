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
 *   histogram starts at bmin <= floor(xmin/w) * w
 *   histogram ends at bmax >= ceil(xmax/w)*w
 *   nb = (bmax-bmin)/w
 *   each score x is counted into bin b = nb - (int) (bmax-x)/w
 *   each bin b contains scores bw+bmin < x <= (b+1)w + bmin
 *
 */  
typedef struct {
  /* The raw data. We always track xmin, xmax, n.
   * Optionally, we keep all the samples, for fitting, for instance.
   * Created dynamically as we use esl_histogram_Add().
   */
  double *x;		/* optional: raw sample values x[0..n-1]            */
  double  xmin, xmax;	/* smallest, largest sample value observed          */
  int     n;            /* total number of raw data samples                 */
  int     nalloc;	/* current allocated size of x                      */

  /* The main (display) histogram is in fixed-width bins, and
   * it's created dynamically as we use esl_histogram_Add().
   */
  int    *obs;		/* observed counts in bin b, 0..nb-1 (dynamic)      */
  double *expect;	/* expected counts in bin b, 0..nb-1 (not resized)  */
  double  bmin, bmax;	/* histogram bounds: all x satisfy bmin <= x < bmax */
  int     imin, imax;	/* smallest, largest bin that contain obs[i] > 0    */
  int     nb;           /* number of bins                                   */
  double  w;		/* fixed width of each bin                          */

  /* A secondary histogram, used for goodness-of-fit testing, is
   * in roughly equal-sized bins of variable width. Requires that 
   * we've kept all the data samples x[]. Created when we call
   * esl_histogram_Finish().
   */
  int    *obs2;		/* obs counts in bin b, 0..nb2-1      */
  double *expect2;	/* expected counts in bin b, 0..nb2-1 */
  double *topx;		/* all values in bin b are <= topx[b] */
  int     nb2;		/* # of bins in secondary histogram   */

  /* Some status flags
   */
  int is_full;		/* TRUE when we're keeping raw data in x */
  int is_finished;	/* TRUE when we've sorted x and binned the 2nd histogram */
} ESL_HISTOGRAM;


extern ESL_HISTOGRAM *esl_histogram_Create    (double bmin, double bmax, double w);
extern ESL_HISTOGRAM *esl_histogram_CreateFull(double bmin, double bmax, double w);
extern void           esl_histogram_Destroy(ESL_HISTOGRAM *h);

extern int    esl_histogram_Add(ESL_HISTOGRAM *h, double x);
extern int    esl_histogram_Finish(ESL_HISTOGRAM *h);

extern int    esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank);

extern int    esl_histogram_Print       (FILE *fp, ESL_HISTOGRAM *h);
extern int    esl_histogram_Plot        (FILE *fp, ESL_HISTOGRAM *h);
extern int    esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h);

#ifdef eslAUGMENT_GUMBEL
extern int    esl_histogram_SetGumbel(ESL_HISTOGRAM *h, double mu, double lambda);
#endif /*eslAUGMENT_GUMBEL*/

#ifdef eslAUGMENT_GEV
extern int    esl_histogram_SetGEV(ESL_HISTOGRAM *h, double mu, double lambda, double alpha);
#endif /*eslAUGMENT_GEV*/

#ifdef eslAUGMENT_STATS
extern int    esl_histogram_GTestGoodness(ESL_HISTOGRAM *h, int ndeg, 
					  int *ret_nbins, double *ret_G, 
					  double *ret_p);
extern int    esl_histogram_ChiSquaredGoodness(ESL_HISTOGRAM *h, int ndeg, 
					       int *ret_nbins, double *ret_X2, 
					       double *ret_p);
#endif /*eslAUGMENT_STATS*/

#endif /*!ESL_HISTOGRAM_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
