/* ratematrix.h
 * Header for evolutionary rate matrix routines in ratemx.c
 * 
 * SRE, Tue Jul 13 16:09:05 2004 [St. Louos]
 * SVN $Id$
 */
#ifndef ESL_RATEMATRIX_INCLUDED
#define ESL_RATEMATRIX_INCLUDED

extern int esl_ratemx_Symm2Q(ESL_DMATRIX *s, double *pi, ESL_DMATRIX **ret_Q);
extern int esl_ratemx_Normalize(ESL_DMATRIX *Q, double *pi, double x);

#endif /*ESL_RATEMATRIX_INCLUDED*/
