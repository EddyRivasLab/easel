/* esl_stats.h
 * Foundation for the statistics modules.
 * 
 * SRE, Tue Jul 19 11:35:28 2005
 * SVN $Id$
 */
#ifndef ESL_STATS_INCLUDED
#define ESL_STATS_INCLUDED

extern int esl_stats_Mean(const double *x, int n, double *opt_mean, double *opt_var);
extern int esl_stats_LogGamma(double x, double *ret_answer);
extern int esl_stats_Psi(double x, double *ret_answer);
extern int esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax);
extern int esl_stats_ChiSquaredTest(int v, double x, double *ret_answer);
extern int esl_stats_LinearRegression(const double *x, const double *y, const double *sigma, int n,
				      double *opt_a,       double *opt_b,
				      double *opt_sigma_a, double *opt_sigma_b, double *opt_cov_ab,
				      double *opt_cc,      double *opt_Q);
#endif /*ESL_STATS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
