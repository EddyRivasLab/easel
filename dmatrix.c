/* dmatrix.c
 * Matrix algebra operations in double-precision matrices.
 * 
 * Implements ESL_DMATRIX (double-precision matrix) and 
 * ESL_PERMUTATION (permutation matrix) objects.
 * 
 * Other easel modules required: easel
 * 
 * SRE, Tue Jul 13 14:42:14 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <easel/easel.h>
#include <easel/dmatrix.h>


/* Function:  esl_dmatrix_Create()
 *
 * Purpose:   Creates an <n> x <m> matrix (<n> rows, <m> 
 *            columns); returns a pointer to the new matrix.
 *
 * Args:      <n> - number of rows;    $>= 1$
 *            <m> - number of columns; $>= 1$
 * 
 * Returns:   A new <ESL_DMATRIX> object. Free with <esl_dmatrix_Destroy()>.
 *
 * Throws:    <NULL> if malloc failed.
 */
ESL_DMATRIX *
esl_dmatrix_Create(int n, int m)
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
  esl_dmatrix_Destroy(A);
  ESL_ERROR_VAL(NULL, ESL_EMEM, "esl_dmx_Create(): malloc failed");
}


/* Function:  esl_dmatrix_Destroy()
 *            
 * Purpose:   Frees an <ESL_DMATRIX> object.
 */
int
esl_dmatrix_Destroy(ESL_DMATRIX *A)
{
  if (A != NULL && A->mx != NULL && A->mx[0] != NULL) free(A->mx[0]);
  if (A != NULL && A->mx != NULL)                     free(A->mx);
  if (A != NULL)                                      free(A);
  return ESL_OK;
}

/* Function:  esl_dmatrix_Dump()
 * Incept:    SRE, Mon Nov 29 19:21:20 2004 [St. Louis]
 *
 * Purpose:   Given a matrix <A>, dump it to stream <ofp> in human-readable
 *            format.
 * 
 *            If <rowlabel> or <collabel> are non-NULL, they represent
 *            single-character labels to put on the rows and columns,
 *            respectively. (For example, these might be a sequence
 *            alphabet for a 4x4 or 20x20 rate matrix or substitution
 *            matrix.)  Numbers 1..ncols or 1..nrows are used if
 *            <collabel> or <rowlabel> are NULL.
 *
 * Args:      ofp      -  output file pointer; stdout, for example.
 *            A        -  matrix to dump.
 *            rowlabel -  optional: NULL, or character labels for rows
 *            collabel -  optional: NULL, or character labels for cols
 *
 * Returns:   <ESL_OK> on success.
 */
int
esl_dmatrix_Dump(FILE *ofp, ESL_DMATRIX *A, char *rowlabel, char *collabel)
{
  int a,b;

  fprintf(ofp, "     ");
  if (collabel != NULL) 
    for (b = 0; b < A->n; b++) fprintf(ofp, "       %c ", collabel[b]);
  else
    for (b = 0; b < A->n; b++) fprintf(ofp, "%8d ", b+1);
  fprintf(ofp, "\n");

  for (a = 0; a < A->n; a++) {
    if (rowlabel != NULL)      fprintf(ofp, "    %c ", rowlabel[a]);
    else                       fprintf(ofp, "%5d ",    a+1);

    for (b = 0; b < A->m; b++) fprintf(ofp, "%8.4f ", A->mx[a][b]);
    fprintf(ofp, "\n");
  }
  return ESL_OK;
}


/* Function:  esl_dmatrix_Copy()
 *
 * Purpose:   Copies <src> matrix into <dest> matrix.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EINCOMPAT> if <src>, <dest> are different sizes
 */
int
esl_dmatrix_Copy(ESL_DMATRIX *src, ESL_DMATRIX *dest)
{
  int i;
  if (dest->n != src->n || dest->m != src->n)
    ESL_ERROR(ESL_EINCOMPAT, "matrices of different size");
  for (i = 0; i < src->n*src->m; i++)
    dest->mx[0][i] = src->mx[0][i];
  return ESL_OK;
}

/* Function:  esl_dmatrix_Compare()
 *
 * Purpose:   Compares matrix <A> to matrix <B>. If all elements
 *            $a_{ij}$ and $b_{ij}$ differ by less than <fabs(tol)>,
 *            return <TRUE>; else return <FALSE>. 
 */
int
esl_dmatrix_Compare(ESL_DMATRIX *A, ESL_DMATRIX *B, double tol)
{
  int i,j;
  if (A->n != B->n) return FALSE;
  if (A->m != B->m) return FALSE;
  for (i = 0; i < A->n; i++)
    for (j = 0; j < A->m; j++)
      if (fabs(A->mx[i][j] - B->mx[i][j]) > tol) return FALSE;
  return TRUE;
}


/* Function:  esl_dmatrix_Set()
 *
 * Purpose:   Set all elements $a_{ij}$ in matrix <A> to <x>,
 *            and returns <ESL_OK>.
 */
int
esl_dmatrix_Set(ESL_DMATRIX *A, double x)
{
  int i;
  for (i = 0; i < A->n*A->m; i++) A->mx[0][i] = x;
  return ESL_OK;
}


/* Function:  esl_dmatrix_SetZero()
 *
 * Purpose:   Sets all elements $a_{ij}$ in matrix <A> to 0.0,
 *            and returns <ESL_OK>.
 */
int
esl_dmatrix_SetZero(ESL_DMATRIX *A)
{
  int i;
  for (i = 0; i < A->n*A->m; i++) A->mx[0][i] = 0.;
  return ESL_OK;
}
  

/* Function:  esl_dmatrix_SetIdentity()
 *
 * Purpose:   Given a square matrix <A>, sets all diagonal elements 
 *            $a_{ii}$ to 1, and all off-diagonal elements to 0.
 *            Returns <ESL_OK>.
 *
 * Throws:    <ESL_EINVAL> if the matrix isn't square.
 */
int
esl_dmatrix_SetIdentity(ESL_DMATRIX *A)
{
  int i;
  
  if (A->n != A->m) ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  esl_dmatrix_SetZero(A);
  for (i = 0; i < A->n; i++) A->mx[i][i] = 1.;
  return ESL_OK;
}
  

/* Function:  esl_permutation_Create()
 *
 * Purpose:   Creates a new permutation "matrix" of size <n> for
 *            permuting <n> x <n> square matrices; returns a 
 *            pointer to it.
 *
 *            A permutation matrix consists of 1's and 0's such that
 *            any given row or column contains only one 1. We store it
 *            more efficiently as a vector; each value $p_i$
 *            represents the column $j$ that has the 1. Thus, on
 *            initialization, $p_i = i$ for all $i = 0..n-1$.
 *
 * Returns:   A new <ESL_PERMUTATION> object. Free with 
 *            <esl_permutation_Destroy()>.
 *
 * Throws:    NULL if allocation fails.
 */
ESL_PERMUTATION *
esl_permutation_Create(int n)
{
  ESL_PERMUTATION *P = NULL;
  int i;

  if ((P = malloc(sizeof(ESL_PERMUTATION))) == NULL) goto FAILURE;
  P->pi = NULL;
  P->n  = n;
  if ((P->pi = malloc(sizeof(int) * n)) == NULL) goto FAILURE;

  esl_permutation_Reuse(P);	/* initialize it */
  return P;

 FAILURE:
  esl_permutation_Free(P);
  ESL_ERROR_VAL(NULL, ESL_EMEM, "esl_permutation_Create: malloc failed");
}
  
/* Function:  esl_permutation_Destroy()
 *
 * Purpose:   Frees an <ESL_PERMUTATION> object.
 */
int
esl_permutation_Destroy(ESL_PERMUTATION *P)
{
  if (P != NULL && P->pi != NULL) free(P->pi);
  if (P != NULL)                  free(P);
  return ESL_OK;
}

/* Function:  esl_permutation_Reuse()
 *
 * Purpose:   Resets a permutation matrix to
 *            $p_i = i$ for all $i = 0..n-1$.
 */
int
esl_permutation_Reuse(ESL_PERMUTATION *P)
{
  int i;
  for (i = 0; i < P->n; i++)
    P->pi[i] = i;
}


/* Function:  esl_permutation_Dump()
 *
 * Purpose:   Given a permutation matrix <P>, dump it to stream <ofp>
 *            in human-readable format.
 *            
 *            If <rowlabel> or <collabel> are non-NULL, they represent
 *            single-character labels to put on the rows and columns,
 *            respectively. (For example, these might be a sequence
 *            alphabet for a 4x4 or 20x20 rate matrix or substitution
 *            matrix.)  Numbers 1..ncols or 1..nrows are used if
 *            <collabel> or <rowlabel> are NULL.
 *
 * Args:      ofp      - output file pointer; stdout, for example
 *            P        - permutation matrix to dump
 *            rowlabel - optional: NULL, or character labels for rows
 *            collabel - optional: NULL, or character labels for cols
 *
 * Returns:   <ESL_OK> on success.
 */
int 
esl_permutation_Dump(FILE *ofp, ESL_PERMUTATION *P, char *rowlabel, char *collabel)
{
  int i,j;

  fprintf(ofp, "    ");
  if (collabel != NULL)
    for (j = 0; j < P->n; j++) fprintf(ofp, "  %c ", collabel[j]);
  else
    for (j = 0; j < P->n; j++) fprintf(ofp, "%3d ", j+1);
  fprintf(ofp, "\n");

  for (i = 0; i < P->n; i++) {
    if (rowlabel != NULL) fprintf(ofp, "  %c ", rowlabel[i]);
    else                  fprintf(ofp, "%3d ", i+1);

    for (j = 0; j < P->n; j++)
      fprintf(ofp, "%3d ", (j == P->pi[i]) ? 1 : 0);
    fprintf(ofp, "\n");
  }
  return ESL_OK;
}



/* Function: esl_dmx_Multiply()
 * 
 * Purpose:  Matrix multiplication: calculate <AB>, store result in <C>.
 *           <A> is nxm; <B> is mxp; <C> is nxp.
 *           Matrix <C> must be allocated appropriately by the caller.
 *           
 * Throws:   <ESL_EINVAL> if matrices don't have compatible dimensions. 
 */
int
esl_dmx_Multiply(ESL_DMATRIX *A, ESL_DMATRIX *B, ESL_DMATRIX *C)
{
  int i, j, k;

  if (A->m != B->n) ESL_ERROR(ESL_EINVAL, "can't multiply A,B");
  if (A->n != C->n) ESL_ERROR(ESL_EINVAL, "A,C # of rows not equal");
  if (B->m != C->m) ESL_ERROR(ESL_EINVAL, "B,C # of cols not equal");

  for (i = 0; i < A->n; i++)
    for (j = 0; j < B->m; j++)
      {
	C->mx[i][j] = 0.;
	for (k = 0; k < A->m; k++)
	  C->mx[i][j] += A->mx[i][k] * B->mx[k][j];
      }
  return ESL_OK;
}

/* Function:  esl_dmx_Transpose()
 *
 * Purpose:   Transpose a square matrix <A> in place.
 *
 * Throws:    <ESL_EINVAL> if <A> isn't square.
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


/* Function:  esl_dmx_Add()
 *
 * Purpose:   <A = A+B>; adds matrix <B> to matrix <A> and leaves result
 *            in matrix <A>.
 *
 * Throws:    <ESL_EINVAL> if matrices aren't the same dimensions.
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

/* Function:  esl_dmx_Scale()
 *
 * Purpose:   Calculates <A = kA>: multiply matrix <A> by scalar
 *            <k> and leave answer in <A>.
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


/* Function:  esl_dmx_AddScale()
 * 
 * Purpose:   Calculates <A + kB>, leaves answer in <A>.
 * 
 * Throws:    <ESL_EINVAL> if matrices aren't the same dimensions.
 */
int
esl_dmx_AddScale(ESL_DMATRIX *A, double k, ESL_DMATRIX *B)
{
  int i,j;

  if (A->n != B->n || A->m != B->n)
    ESL_ERROR(ESL_EINCOMPAT, "matrices of different size");
  for (i = 0; i < A->n; i++)
    for (j = 0; j < A->m; j++)
      A->mx[i][j] +=  k * B->mx[i][j];
  return ESL_OK;
}


/* Function:  esl_dmx_Permute_PA()
 *
 * Purpose:   Computes <B = PA>: do a row-wise permutation of a square
 *            matrix <A>, using the permutation matrix <P>, and put
 *            the result in a square matrix <B> that the caller has
 *            allocated.
 *
 * Throws:    <ESL_EINVAL> if A, B, P do not have compatible dimensions.
 */
int
esl_dmx_Permute_PA(ESL_PERMUTATION *P, ESL_DMATRIX *A, ESL_DMATRIX *B)
{
  int i,ip,j;

  if (A->n != P->n || A->n != B->n || A->n != A->m || B->n != B->m)
    ESL_ERROR(ESL_EINVAL, "matrix dimensions not compatible");

  for (i = 0; i < A->n; i++)
    {
      ip = P->pi[i];
      for (j = 0; j < A->m; j++)
	B->mx[i][j] = A->mx[ip][j];
    }
  return ESL_OK;
}

/* Function:  esl_dmx_LUP_decompose()
 *
 * Purpose:   Calculates a permuted LU decomposition of square matrix
 *            <A>; upon return, <A> is replaced by this decomposition,
 *            where <U> is in the lower triangle (inclusive of the 
 *            diagonal) and <L> is the upper triangle (exclusive of
 *            diagonal, which is 1's by definition), and <P> is the
 *            permutation matrix. Caller provides an allocated 
 *            permutation matrix <P> compatible with the square matrix
 *            <A>.
 *            
 *            Implements Gaussian elimination with pivoting; xref
 *            [Cormen, Leiserson, Rivest; "Algorithms", MIT Press,
 *            1999; p.759].
 *
 * Throws:    <ESL_EINVAL> if <A> isn't square, or if <P> isn't the right
 *            size for <A>.
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


/* Function:  esl_dmx_LU_separate()
 *
 * Purpose:   Separate a square <LU> decomposition matrix into its two
 *            triangular matrices <L> and <U>. Caller provides two
 *            allocated <L> and <U> matrices of same size as <LU> for
 *            storing the results.
 *
 * Throws:    <ESL_EINVAL> if <LU>, <L>, <U> are not of compatible dimensions.
 */
int
esl_dmx_LU_separate(ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U)
{
  int i,j;

  if (LU->n != LU->m) ESL_ERROR(ESL_EINVAL, "LU isn't square");
  if (L->n  != L->m)  ESL_ERROR(ESL_EINVAL, "L isn't square");
  if (U->n  != U->m)  ESL_ERROR(ESL_EINVAL, "U isn't square");
  if (LU->n != L->n)  ESL_ERROR(ESL_EINVAL, "LU, L have incompatible dimensions");
  if (LU->n != U->n)  ESL_ERROR(ESL_EINVAL, "LU, U have incompatible dimensions");

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

/* Function:  esl_dmx_Invert()
 *
 * Purpose:   Calculates the inverse of square matrix <A>, and stores the
 *            result in matrix <Ai>. Caller provides an allocated
 *            matrix <Ai> of same dimensions as <A>.
 *            
 *            Peforms the inversion by LUP decomposition followed by 
 *            forward/back-substitution; xref [Cormen, Leiserson, 
 *            Rivest; "Algorithms", MIT Press 1999; p.753].
 *
 * Throws:    <ESL_EINVAL> if <A>, <Ai> do not have same dimensions, or
 *                         if <A> isn't square.
 *            <ESL_EMEM>   if internal allocations (for LU, and some other
 *                         bookkeeping) fail.
 */
int
esl_dmx_Invert(ESL_DMATRIX *A, ESL_DMATRIX *Ai)
{
  ESL_DMATRIX      *LU = NULL;
  ESL_PERMUTATION  *P  = NULL;
  double           *y  = NULL;	/* column vector, intermediate calculation   */
  double           *b  = NULL;	/* column vector of permuted identity matrix */
  int               i,j,k;

  if (A->n != A->m)                   ESL_ERROR(ESL_EINVAL, "matrix isn't square");
  if (A->n != Ai->n || A->m != Ai->m) ESL_ERROR(ESL_EINVAL, "matrices are different size");

  /* Copy A to LU, and do an LU decomposition.
   */
  if ((LU = esl_dmatrix_Create(A->n, A->m)) == NULL) goto FAILURE;
  if ((P  = esl_permutation_Create(A->n))   == NULL) goto FAILURE;
  esl_dmatrix_Copy(A, LU);
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
  if ((b = malloc(sizeof(double) * A->n)) == NULL) goto FAILURE;
  if ((y = malloc(sizeof(double) * A->n)) == NULL) goto FAILURE;

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
  esl_dmatrix_Destroy(LU);
  esl_permutation_Destroy(P);
  return ESL_OK;

 FAILURE:
  if (y  != NULL) free(y);
  if (b  != NULL) free(b);
  if (LU != NULL) esl_dmatrix_Destroy(LU);
  if (P  != NULL) esl_permutation_Destroy(P);
  ESL_ERROR(ESL_EMEM, "allocation failed in esl_dmx_Invert()");
}

