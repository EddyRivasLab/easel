/* gamma.h
 * SRE, Tue Nov  2 14:17:09 2004 [St. Louis]
 * SVN $Id$
 * 
 * Functions for gamma functions and the gamma density.
 */
#ifndef ESL_GAMMA_INCLUDED
#define ESL_GAMMA_INCLUDED

#include <easel/random.h>

extern int esl_gamma_log(double x, double *ret_answer);
extern int esl_gamma_sample(ESL_RANDOMNESS *r, double a, double *ret_answer);

#endif /*ESL_GAMMA_INCLUDED*/
