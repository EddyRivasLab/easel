/* esl_gev.h
 * Generalized extreme value (GEV) distributions.
 * 
 * SRE, Tue Jul 12 09:15:56 2005
 * SVN $Id$
 */
#ifndef ESL_GEV_INCLUDED
#define ESL_GEV_INCLUDED



extern double esl_gev_pdf    (double x, double mu, double lambda, double alpha);
extern double esl_gev_logpdf (double x, double mu, double lambda, double alpha);
extern double esl_gev_cdf    (double x, double mu, double lambda, double alpha);
extern double esl_gev_logcdf (double x, double mu, double lambda, double alpha);
extern double esl_gev_surv   (double x, double mu, double lambda, double alpha);
extern double esl_gev_logsurv(double x, double mu, double lambda, double alpha);


#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
extern double esl_gev_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double alpha);
#endif

#ifdef eslAUGMENT_MINIMIZER
extern int esl_gev_FitComplete(double *x, int n, 
			       double *ret_mu, double *ret_lambda, 
			       double *ret_alpha);
#endif /*eslAUGMENT_MINIMIZER*/


#endif /*ESL_GEV_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
