/* dmatrix.c
 * 
 * SRE, Tue Jul 13 14:42:14 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <easel/easel.h>
#include <easel/dmatrix.h>

ESL_DMATRIX *
esl_dmx_Alloc(int n, int m)
{
  ESL_DMATRIX *A = NULL;
  int r;

  if ((A     = malloc(sizeof(ESL_DMATRIX))) == NULL) goto FAILURE;
  A->mx = NULL;
  A->n  = n;
  A->m  = m;

  if ((A->mx = malloc(sizeof(double *) * n)) == NULL) goto FAILURE;
  A->mx[0] = NULL;
  if ((A->mx[0] = malloc(sizeof(double) * n * m)) == NULL) goto FAILURE;

  for (r = 1; r < n; r++)
    A->mx[r] = A->mx[0] + r*n;
  return A;
  
 FAILURE:
  esl_dmx_Free(A);
  ESL_ERROR_VAL(NULL, ESL_EMEM, "malloc failed");
}

int
esl_dmx_Free(ESL_DMATRIX *A)
{
  if (A != NULL && A->mx != NULL && A->mx[0] != NULL) free(A->mx[0]);
  if (A != NULL && A->mx != NULL)                     free(A->mx);
  if (A != NULL)                                      free(A);
  return ESL_OK;
}

int
esl_dmx_Copy(ESL_DMATRIX *src, ESL_DMATRIX *dest)
{
  int i;
  if (dest->n != src->n || dest->m != src->n)
    ESL_ERROR(ESL_EINCOMPAT, "matrices of different size");
  for (i = 0; i < src->n*src->m; i++)
    dest->mx[0][i] = src->mx[0][i];
  return ESL_OK;
}

int
esl_dmx_MatricesEqual(ESL_DMATRIX *A, ESL_DMATRIX *B, double tol)
{
  int i,j;
  if (A->n != B->n) return FALSE;
  if (A->m != B->m) return FALSE;
  for (i = 0; i < A->n; i++)
    for (j = 0; j < A->m; j++)
      if (fabs(A->mx[i][j] - B->mx[i][j]) > tol) return FALSE;
  return TRUE;
}


int
esl_dmx_SetAll(ESL_DMATRIX *A, double x)
{
  int i;
  for (i = 0; i < A->n*A->m; i++) A->mx[0][i] = x;
  return ESL_OK;
}

/* zero an MxN matrix
 */
int
esl_dmx_SetZero(ESL_DMATRIX *A)
{
  int i;
  for (i = 0; i < A->n*A->m; i++) A->mx[0][i] = 0.;
  return ESL_OK;
}
  
/* set a matrix to the identity matrix, a_ii = 1, a_ij = 0 \forall j \neq i
 */
int
esl_dmx_SetIdentity(ESL_DMATRIX *A)
{
  int i;
  
  if (A->n != A->m) ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  esl_dmx_SetZero(A);
  for (i = 0; i < A->n; i++) A->mx[i][i] = 1.;
  return ESL_OK;
}
  
/* Function: esl_dmx_Multiply()
 * 
 * Purpose:  Matrix multiplication.
 *           Multiply AB, giving C.
 *           A is nxm; B is mxp; C is nxp.
 *           Matrix C must be preallocated.
 */
int
esl_dmx_Multiply(ESL_DMATRIX *A, ESL_DMATRIX *B, ESL_DMATRIX *C)
{
  int i, j, k;

  if (A->m != B->n) ESL_ERROR(ESL_EINVAL, "can't multiply those");

  for (i = 0; i < A->n; i++)
    for (j = 0; j < B->m; j++)
      {
	C->mx[i][j] = 0.;
	for (k = 0; k < A->m; k++)
	  C->mx[i][j] += A->mx[i][k] * B->mx[k][j];
      }
  return ESL_OK;
}

/* Transpose A, in place.
 * Matrix must be square.
 */
int
esl_dmx_Transpose(ESL_DMATRIX *A)
{
  int    i,j;
  double swap;

  if (A->n != A->m) ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  for (i = 0; i < A->n; i++)
    for (j = i+1; j < A->m; j++)
      { swap = A->mx[i][j]; A->mx[i][j] = A->mx[j][i]; A->mx[j][i] = swap; }
  return ESL_OK;
}


/* Calculates A + B, leave answer in A.
 */
int
esl_dmx_Add(ESL_DMATRIX *A, ESL_DMATRIX *B)
{
  int i,j;
  
  if (A->n != B->n || A->m != B->n)
    ESL_ERROR(ESL_EINCOMPAT, "matrices of different size");
  for (i = 0; i < A->n; i++)
    for (j = 0; j < A->m; j++)
      A->mx[i][j] +=  B->mx[i][j];
  return ESL_OK;
}

/* Calculates kA, leaves answer in A
 */
int 
esl_dmx_Scale(ESL_DMATRIX *A, double k)
{
  int i,j;
  for (i = 0; i < A->n; i++)
    for (j = 0; j < A->m; j++)
      A->mx[i][j] *=  k;
  return ESL_OK;
}
  


/* Output an alphabet-labeled square matrix of floating point #'s.
 */
int
esl_dmx_fprintf_alphalabeled(FILE *ofp, ESL_DMATRIX *A, char *alphabet)
{
  int a,b;

  if (A->n != A->m) ESL_ERROR(ESL_EINVAL, "matrix isn't square");

  fprintf(ofp, "  ");
  for (b = 0; b < A->n; b++)
    fprintf(ofp, "       %c ", alphabet[b]);
  fprintf(ofp, "\n");
  for (a = 0; a < A->n; a++) {
    fprintf(ofp, "%c ", alphabet[a]);
    for (b = 0; b < A->m; b++)
      fprintf(ofp, "%8.4f ", A->mx[a][b]);
    fprintf(ofp, "\n");
  }
  return ESL_OK;
}


ESL_PERMUTATION *
esl_permutation_Alloc(int n)
{
  ESL_PERMUTATION *P = NULL;
  int i;

  if ((P = malloc(sizeof(ESL_PERMUTATION))) == NULL) goto FAILURE;
  P->pi = NULL;
  P->n  = n;
  if ((P->pi = malloc(sizeof(int) * n)) == NULL) goto FAILURE;
  for (i = 0; i < n; i++)
    P->pi[i] = i;
  return P;

 FAILURE:
  esl_permutation_Free(P);
  ESL_ERROR_VAL(NULL, ESL_EMEM, "malloc failed");
}
  
int
esl_permutation_Init(ESL_PERMUTATION *P)
{
  int i;
  for (i = 0; i < P->n; i++)
    P->pi[i] = i;
}

int 
esl_permutation_fprintf_numlabeled(FILE *ofp, ESL_PERMUTATION *P)
{
  int i,j;

  fprintf(ofp, "    ");
  for (j = 0; j < P->n; j++)
    fprintf(ofp, " %3d ", j);
  fprintf(ofp, "\n");
  for (i = 0; i < P->n; i++) {
    fprintf(ofp, "%3d ", i);
    for (j = 0; j < P->n; j++)
      fprintf(ofp, " %3d ", (j == P->pi[i]) ? 1 : 0);
    fprintf(ofp, "\n");
  }
  return ESL_OK;
}

int
esl_permutation_Free(ESL_PERMUTATION *P)
{
  if (P != NULL && P->pi != NULL) free(P->pi);
  if (P != NULL)                  free(P);
  return ESL_OK;
}

/* Compute B = PA: a row-wise permutation of A
 */
int
esl_permute_PA(ESL_PERMUTATION *P, ESL_DMATRIX *A, ESL_DMATRIX *B)
{
  int i,ip,j;

  for (i = 0; i < A->n; i++)
    {
      ip = P->pi[i];
      for (j = 0; j < A->m; j++)
	B->mx[i][j] = A->mx[ip][j];
    }
  return ESL_OK;
}

/* 
 * Upon return, A is replaced by LU:
 *    U is in upper triangle (inclusive of diagonal)
 *    L is lower triangle (exclusive of diagonal, which is 1's by definition)
 *    
 * Algorithm: Gaussian elimination, with pivoting;
 *            [Cormen, Leiserson, Rivest; _Algorithms_, MIT Press 1999; p.759]
 */
int
esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P)
{
  int    i,j,k,kpiv;
  double max;
  double swap;

  if (A->n != A->m) ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  if (P->n != A->n) ESL_ERROR(ESL_EINVAL, "permutation isn't the right size");
  esl_permutation_Init(P);

  for (k = 0; k < A->n-1; k++)
    {
      /* Identify our pivot;
       * find row with maximum value in col[k].
       */
      max = 0.; 
      for (i = k; i < A->n; i++)
	if (fabs(A->mx[i][k]) > max) {
	  max = fabs(A->mx[i][k]);
	  kpiv = i;
	}
      if (max == 0.) ESL_ERROR(ESL_EDIVZERO, "matrix is singular");
      
      /* Swap those rows (k and kpiv);
       * and keep track of that permutation in P. (misuse j for swapping integers)
       */
      j = P->pi[k]; P->pi[k] = P->pi[kpiv]; P->pi[kpiv] = j;
      for (j = 0; j < A->m; j++)
	{ swap = A->mx[k][j]; A->mx[k][j] = A->mx[kpiv][j]; A->mx[kpiv][j] = swap; }

      /* Gaussian elimination for all rows k+1..n.
       */
      for (i = k+1; i < A->n; i++)
	{
	  A->mx[i][k] /= A->mx[k][k];
	  for (j = k+1; j < A->m; j++)
	    A->mx[i][j] -= A->mx[i][k] * A->mx[k][j];
	}
    }
  return ESL_OK;
}


/* Separate a LU decomposition matrix into its two 
 * triangular matrices L and U.
 */
int
esl_dmx_LU_separate(ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U)
{
  int i,j;

  esl_dmx_SetZero(L);
  esl_dmx_SetZero(U);

  for (i = 0; i < LU->n; i++)
    for (j = i; j < LU->m; j++)
      U->mx[i][j] = LU->mx[i][j];

  for (i = 0; i < LU->n; i++) 
    {
      L->mx[i][i] = 1.;
      for (j = 0; j < i; j++)
	L->mx[i][j] = LU->mx[i][j];
    }
  return ESL_OK;
}


/* Invert an NxN square matrix; leave the result in Ai,
 * which the caller allocated.
 * 
 * Algorithm is LUP decomposition, followed by solving
 * for the inverse by forward/back-substitution.
 * 
 * Reference : [Cormen, Leiserson, Rivest; _Algorithms_, MIT Press 1999; p.753]
 */
int
esl_dmx_Invert(ESL_DMATRIX *A, ESL_DMATRIX *Ai)
{
  ESL_DMATRIX      *LU;
  ESL_PERMUTATION  *P;
  double           *y;		/* column vector, intermediate calculation   */
  double           *b;		/* column vector of permuted identity matrix */
  int               i,j,k;

  if (A->n != A->m)  ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  if (A->n != Ai->n || A->m != Ai->m) ESL_ERROR(ESL_EINVAL, "matrices are different size");
  /* Copy A to LU
   */
  LU = esl_dmx_Alloc(A->n, A->m);
  P  = esl_permutation_Alloc(A->n);
  esl_dmx_Copy(A, LU);
  esl_dmx_LUP_decompose(LU, P);

  /* Now we have:
   *   PA = LU
   *   
   * to invert a matrix A, we want A A^-1 = I;
   * that's PAx = Pb, for columns x of A^-1 and b of the identity matrix;
   * and that's n equations LUx = Pb;
   * 
   * so, solve Ly = Pb for y by forward substitution;
   * then Ux = y by back substitution;
   * x is then a column of A^-1.
   * 
   * Do that for all columns.
   */
  b  = malloc(sizeof(double) * A->n);
  y  = malloc(sizeof(double) * A->n);
  for (k = 0; k < A->m; k++)	/* for each column... */
    {
      /* build Pb for column j of the identity matrix */
      for (i = 0; i < A->n; i++)
	if (P->pi[i] == k) b[i] = 1.; else b[i] = 0.;

      /* forward substitution
       */
      for (i = 0; i < A->n; i++)
	{
	  y[i] = b[i];
	  for (j = 0; j < i; j++) y[i] -= LU->mx[i][j] * y[j];
	}

      /* back substitution
       */
      for (i = A->n-1; i >= 0; i--)
	{
	  Ai->mx[i][k] = y[i];
	  for (j = i+1; j < A->n; j++) Ai->mx[i][k] -= LU->mx[i][j] * Ai->mx[j][k];
	  Ai->mx[i][k] /= LU->mx[i][i];
	}
    }

  free(b);
  free(y);
  esl_dmx_Free(LU);
  esl_permutation_Free(P);
  return ESL_OK;
}

