/* dirichlet.h
 * Functions for the Dirichlet density.
 * 
 * SRE, Tue Nov  2 14:35:06 2004 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_DIRICHLET_INCLUDED
#define ESL_DIRICHLET_INCLUDED

#include <easel/random.h>

extern int esl_dirichlet_LogProbData(double *c, double *alpha, int K, double *ret_answer);
extern int esl_dirichlet_sample(ESL_RANDOMNESS *r, double *alpha, int K, double *p);

#endif /*ESL_DIRICHLET_INCLUDED*/
