/* dmatrix.h
 * 
 * SRE, Tue Jul 13 14:41:07 2004 [St. Louis]
 * SVN $Id$
 */

#ifndef ESL_DMATRIX_INCLUDED
#define ESL_DMATRIX_INCLUDED

#include <stdio.h>

struct esl_dmx_s {
  double **mx;			/* matrix: mx, mx[0] are allocated. */
  int      n;			/* rows    */
  int      m;			/* columns */
};
typedef struct esl_dmx_s ESL_DMATRIX;

struct esl_permutation_s {
  int     *pi;
  int      n;
};
typedef struct esl_permutation_s ESL_PERMUTATION;

extern ESL_DMATRIX *esl_dmx_Alloc(int n, int m);
extern int          esl_dmx_Free(ESL_DMATRIX *A);
extern int          esl_dmx_Copy(ESL_DMATRIX *src, ESL_DMATRIX *dest);
extern int          esl_dmx_MatricesEqual(ESL_DMATRIX *A, ESL_DMATRIX *B, double tol);
extern int          esl_dmx_SetAll(ESL_DMATRIX *A, double x);
extern int          esl_dmx_SetZero(ESL_DMATRIX *A);
extern int          esl_dmx_SetIdentity(ESL_DMATRIX *A);
extern int          esl_dmx_Multiply(ESL_DMATRIX *A, ESL_DMATRIX *B, ESL_DMATRIX *C);
extern int          esl_dmx_Transpose(ESL_DMATRIX *A);
extern int          esl_dmx_Add(ESL_DMATRIX *A, ESL_DMATRIX *B);
extern int          esl_dmx_Scale(ESL_DMATRIX *A, double k);
extern int          esl_dmx_fprintf_alphalabeled(FILE *ofp, ESL_DMATRIX *A, char *alphabet);
extern int          esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P);
extern int          esl_dmx_LU_separate(ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U);
extern int          esl_dmx_Invert(ESL_DMATRIX *A, ESL_DMATRIX *Ai);


extern ESL_PERMUTATION *esl_permutation_Alloc(int n);
extern int esl_permutation_Init(ESL_PERMUTATION *P);
extern int esl_permutation_fprintf_numlabeled(FILE *ofp, ESL_PERMUTATION *P);
extern int esl_permutation_Free(ESL_PERMUTATION *P);
extern int esl_permute_PA(ESL_PERMUTATION *P, ESL_DMATRIX *A, ESL_DMATRIX *B);

#endif /*ESL_DMATRIX_INCLUDED*/
