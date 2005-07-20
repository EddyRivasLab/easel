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
 * size (width) w. A score of x is counted into bin i = (x-xmin)/w,
 * and bin i contains scores iw + xmin <= x < (i+1)w + xmin.
 */  
typedef struct {
  int    *bin;		/* count bins, 0..nbins-1                                        */
  int     nbins;        /* nbins  = 1 + ((xmax-xmin) / w)                                */
  double  w;		/* width of each bin                                             */
  double  xmin, xmax;	/* allocation bounds: all scores x must satisfy xmin <= x < xmax */
  int     imin, imax;	/* indices of smallest, largest bins with any counts             */

  int     total;	/* total # of observed counts in bins    */
  double *expect;	/* expected counts in bins, 0..nbins-1   */
} ESL_HISTOGRAM;


extern ESL_HISTOGRAM *esl_histogram_Create(double xmin, double xmax, double w);
extern void           esl_histogram_Destroy(ESL_HISTOGRAM *h);

extern int    esl_histogram_Add(ESL_HISTOGRAM *h, double x);
extern int    esl_histogram_Print(FILE *fp, ESL_HISTOGRAM *h);
extern void   esl_histogram_Plot(FILE *fp, ESL_HISTOGRAM *h);

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
