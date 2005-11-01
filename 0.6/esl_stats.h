/* esl_stats.h
 * Foundation for the statistics modules.
 * 
 * SRE, Tue Jul 19 11:35:28 2005
 * SVN $Id$
 */
#ifndef ESL_STATS_INCLUDED
#define ESL_STATS_INCLUDED

extern int esl_stats_Mean(double *x, int n, double *ret_mean, double *ret_var);
extern int esl_stats_LogGamma(double x, double *ret_answer);
extern int esl_stats_IncompleteGamma(double a, double x, double *ret_answer);
extern int esl_stats_ChiSquaredTest(int v, double x, double *ret_answer);

#endif /*ESL_STATS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
