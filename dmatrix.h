/* dmatrix.h
 * 
 * SRE, Tue Jul 13 14:41:07 2004 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_DMATRIX_INCLUDED
#define ESL_DMATRIX_INCLUDED

#include <stdio.h>


struct esl_dmatrix_s {
  double **mx;			/* matrix: mx, mx[0] are allocated. */
  int      n;			/* rows    */
  int      m;			/* columns */
};
typedef struct esl_dmatrix_s ESL_DMATRIX;

extern ESL_DMATRIX *esl_dmatrix_Create(int n, int m);
extern int          esl_dmatrix_Destroy(ESL_DMATRIX *A);
extern int          esl_dmatrix_Dump(FILE *ofp, ESL_DMATRIX *A, char *rowlabel, char *collabel);
extern int          esl_dmatrix_Copy(ESL_DMATRIX *src, ESL_DMATRIX *dest);
extern int          esl_dmatrix_Compare(ESL_DMATRIX *A, ESL_DMATRIX *B, double tol);
extern int          esl_dmatrix_Set(ESL_DMATRIX *A, double x);
extern int          esl_dmatrix_SetZero(ESL_DMATRIX *A);
extern int          esl_dmatrix_SetIdentity(ESL_DMATRIX *A);


struct esl_permutation_s {
  int     *pi;
  int      n;
};
typedef struct esl_permutation_s ESL_PERMUTATION;

extern ESL_PERMUTATION *esl_permutation_Create(int n);
extern int              esl_permutation_Destroy(ESL_PERMUTATION *P);
extern int              esl_permutation_Reuse(ESL_PERMUTATION *P);
extern int              esl_permutation_Dump(FILE *ofp, ESL_PERMUTATION *P, char *rowlabel, char *collabel);




extern int          esl_dmx_Multiply(ESL_DMATRIX *A, ESL_DMATRIX *B, ESL_DMATRIX *C);
extern int          esl_dmx_Transpose(ESL_DMATRIX *A);
extern int          esl_dmx_Add(ESL_DMATRIX *A, ESL_DMATRIX *B);
extern int          esl_dmx_Scale(ESL_DMATRIX *A, double k);
extern int          esl_dmx_AddScale(ESL_DMATRIX *A, double k, ESL_DMATRIX *B);
extern int          esl_dmx_Permute_PA(ESL_PERMUTATION *P, ESL_DMATRIX *A, ESL_DMATRIX *B);
extern int          esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P);
extern int          esl_dmx_LU_separate(ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U);
extern int          esl_dmx_Invert(ESL_DMATRIX *A, ESL_DMATRIX *Ai);

#endif /*ESL_DMATRIX_INCLUDED*/
