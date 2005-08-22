/* esl_hyperexp.h
 * Hyperexponential (mixture exponential) distributions.
 * 
 * SRE, Mon Aug 15 08:32:33 2005 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_HYPEREXP_INCLUDED
#define ESL_HYPEREXP_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

#ifdef eslAUGMENT_HISTOGRAM
#include <esl_histogram.h>
#endif

typedef struct {
  double *q;			/* mixture coefficients   [0..K-1]*/
  double *lambda;		/* scale params           [0..K-1]*/
  double *wrk;			/* tmp K-vector, for logpdf calc  */
  double  mu;			/* location (x offset) parameter  */
  int     K;			/* # of components                */
  char   *fixlambda;		/* TRUE to constrain a lambda val */
} ESL_HYPEREXP;


extern ESL_HYPEREXP *esl_hyperexp_Create(int K);
extern void          esl_hyperexp_Destroy(ESL_HYPEREXP *h);
extern int           esl_hyperexp_Copy(ESL_HYPEREXP *src, ESL_HYPEREXP *dest);

extern double  esl_hxp_pdf    (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logpdf (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_cdf    (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logcdf (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_surv   (double x, ESL_HYPEREXP *h);
extern double  esl_hxp_logsurv(double x, ESL_HYPEREXP *h);
extern double  esl_hxp_invcdf (double p, ESL_HYPEREXP *h);

extern double  esl_hxp_generic_cdf   (double x, void *params);
extern double  esl_hxp_generic_invcdf(double x, void *params);

extern int esl_hxp_Plot(FILE *fp, ESL_HYPEREXP *h,
			double (*func)(double x, ESL_HYPEREXP *h), 
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
extern double esl_hxp_Sample(ESL_RANDOMNESS *r, ESL_HYPEREXP *h);
#endif

#ifdef eslAUGMENT_MINIMIZER
extern int esl_hxp_FitGuess   (double *x, int n, ESL_HYPEREXP *h);
extern int esl_hxp_FitComplete(double *x, int n, ESL_HYPEREXP *h);
#ifdef eslAUGMENT_HISTOGRAM
extern int esl_hxp_FitGuessBinned   (ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
extern int esl_hxp_FitCompleteBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
#endif
#endif

#endif /*ESL_HYPEREXP_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
