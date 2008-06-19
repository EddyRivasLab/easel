/* Statistical routines for normal distributions
 * 
 * SRE, Tue Nov 21 14:29:02 2006 [Janelia]
 * SVN $Id$
 */

#ifndef ESL_NORMAL_INCLUDED
#define ESL_NORMAL_INCLUDED

extern double esl_normal_pdf   (double x, double mu, double sigma);
extern double esl_normal_logpdf(double x, double mu, double sigma);
extern double esl_normal_cdf   (double x, double mu, double sigma);
extern double esl_normal_surv  (double x, double mu, double sigma);

#endif /*ESL_NORMAL_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
