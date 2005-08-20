/* esl_exponential.h
 * Exponential distributions.
 * 
 * SRE, Wed Aug 10 08:32:45 2005 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_EXP_INCLUDED
#define ESL_EXP_INCLUDED

extern double esl_exp_pdf    (double x, double mu, double lambda);
extern double esl_exp_logpdf (double x, double mu, double lambda);
extern double esl_exp_cdf    (double x, double mu, double lambda);
extern double esl_exp_logcdf (double x, double mu, double lambda);
extern double esl_exp_surv   (double x, double mu, double lambda);
extern double esl_exp_logsurv(double x, double mu, double lambda);

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
extern double esl_exp_Sample(ESL_RANDOMNESS *r, double mu, double lambda);
#endif

extern int esl_exp_FitComplete(double *x, int n, double mu, double *ret_lambda);


#endif /*ESL_SXP_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
