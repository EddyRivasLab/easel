/* dmatrix.c
 * Matrix algebra operations in double-precision matrices.
 * 
 * Implements ESL_DMATRIX (double-precision matrix) and 
 * ESL_PERMUTATION (permutation matrix) objects.
 * 
 * To do:
 *   - table of contents here for .c,.h file, section splits
 *   - eventually probably want additional matrix types
 *   - unit tests poor 
 *   
 *
 * SRE, Tue Jul 13 14:42:14 2004 [St. Louis]
 * SVN $Id$
 */
#include <esl_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <easel.h>
#include <esl_dmatrix.h>


/* Function:  esl_dmatrix_Create()
 *
 * Purpose:   Creates a general <n> x <m> matrix (<n> rows, <m> 
 *            columns).
 *
 * Args:      <n> - number of rows;    $>= 1$
 *            <m> - number of columns; $>= 1$
 * 
 * Returns:   a pointer to a new <ESL_DMATRIX> object. Caller frees
 *            with <esl_dmatrix_Destroy()>.
 *
 * Throws:    <NULL> if an allocation failed.
 */
ESL_DMATRIX *
esl_dmatrix_Create(int n, int m)
{
  ESL_DMATRIX *A = NULL;
  int r;
  int status;

  ESL_ALLOC(A, sizeof(ESL_DMATRIX));
  A->mx = NULL;
  A->n  = n;
  A->m  = m;

  ESL_ALLOC(A->mx, sizeof(double *) * n);
  A->mx[0] = NULL;

  ESL_ALLOC(A->mx[0], sizeof(double) * n * m);
  for (r = 1; r < n; r++)
    A->mx[r] = A->mx[0] + r*n;

  A->type   = eslGENERAL;
  A->ncells = n * m; 
  return A;
  
 ERROR:
  esl_dmatrix_Destroy(A);
  return NULL;
}


/* Function:  esl_dmatrix_CreateUpper()
 * Incept:    SRE, Wed Feb 28 08:45:45 2007 [Janelia]
 *
 * Purpose:   Creates a packed upper triangular matrix of <n> rows and
 *            <n> columns. Caller may only access cells $i \leq j$.
 *            Cells $i > j$ are not stored and are implicitly 0.
 *            
 *            Not all matrix operations in Easel can work on packed
 *            upper triangular matrices.
 *
 * Returns:   a pointer to a new <ESL_DMATRIX> object of type
 *            <eslUPPER>. Caller frees with <esl_dmatrix_Destroy()>.
 *
 * Throws:    <NULL> if allocation fails.
 *
 * Xref:      J1/10
 */
ESL_DMATRIX *
esl_dmatrix_CreateUpper(int n)
{
  int status;
  ESL_DMATRIX *A = NULL;
  int r;			/* counter over rows */
  int nc;			/* cell counter */

  /* matrix structure allocation */
  ESL_ALLOC(A, sizeof(ESL_DMATRIX)); 
  A->mx = NULL;
  A->n  = n;
  A->m  = n;

  /* n row ptrs */
  ESL_ALLOC(A->mx, sizeof(double *) * n); 
  A->mx[0] = NULL;

  /* cell storage */
  ESL_ALLOC(A->mx[0], sizeof(double) * n * (n+1) / 2);
  
  /* row pointers set in a tricksy overlapping way, so
   * mx[i][j] access works normally but only i<=j are valid.
   * xref J1/10.
   */
  nc = n;  /* nc is the number of valid cells assigned to rows so far */
  for (r = 1; r < n; r++) {
    A->mx[r] = A->mx[0] + nc - r; /* -r overlaps this row w/ previous row */
    nc += n-r;
  }
  A->type   = eslUPPER;
  A->ncells = n * (n+1) / 2; 
  return A;

 ERROR:
  esl_dmatrix_Destroy(A);
  return NULL;
}



/* Function:  esl_dmatrix_Destroy()
 *            
 * Purpose:   Frees an <ESL_DMATRIX> object <A>.
 */
int
esl_dmatrix_Destroy(ESL_DMATRIX *A)
{
  if (A != NULL && A->mx != NULL && A->mx[0] != NULL) free(A->mx[0]);
  if (A != NULL && A->mx != NULL)                     free(A->mx);
  if (A != NULL)                                      free(A);
  return eslOK;
}

/* Function:  esl_dmatrix_Dump()
 * Incept:    SRE, Mon Nov 29 19:21:20 2004 [St. Louis]
 *
 * Purpose:   Given a matrix <A>, dump it to output stream <ofp> in human-readable
 *            format.
 * 
 *            If <rowlabel> or <collabel> are non-NULL, they specify a
 *            string of single-character labels to put on the rows and
 *            columns, respectively. (For example, these might be a
 *            sequence alphabet for a 4x4 or 20x20 rate matrix or
 *            substitution matrix.)  Numbers <1..ncols> or <1..nrows> are
 *            used if <collabel> or <rowlabel> are passed as <NULL>.
 *
 * Args:      ofp      -  output file pointer; stdout, for example.
 *            A        -  matrix to dump.
 *            rowlabel -  optional: NULL, or character labels for rows
 *            collabel -  optional: NULL, or character labels for cols
 *
 * Returns:   <eslOK> on success.
 */
int
esl_dmatrix_Dump(FILE *ofp, ESL_DMATRIX *A, char *rowlabel, char *collabel)
{
  int a,b;

  fprintf(ofp, "     ");
  if (collabel != NULL) 
    for (b = 0; b < A->m; b++) fprintf(ofp, "       %c ", collabel[b]);
  else
    for (b = 0; b < A->m; b++) fprintf(ofp, "%8d ", b+1);
  fprintf(ofp, "\n");

  for (a = 0; a < A->n; a++) {
    if (rowlabel != NULL)      fprintf(ofp, "    %c ", rowlabel[a]);
    else                       fprintf(ofp, "%5d ",    a+1);

    for (b = 0; b < A->m; b++) {
      switch (A->type) {
      case eslUPPER:
	if (a > b) 	fprintf(ofp, "%8s ", "");
	else            fprintf(ofp, "%8.4f ", A->mx[a][b]); 
	break;

       default: case eslGENERAL:
	fprintf(ofp, "%8.4f ", A->mx[a][b]); 
	break;
      }
      fprintf(ofp, "\n");
    }
  }
  return eslOK;
}


/* Function:  esl_dmatrix_Copy()
 *
 * Purpose:   Copies <src> matrix into <dest> matrix. <dest> must
 *            be allocated already by the caller.
 * 
 *            You may copy to a matrix of a different type, so long as
 *            the copy makes sense. If <dest> matrix is a packed type
 *            and <src> is not, the values that should be zeros must
 *            be zero in <src>, else the routine throws
 *            <eslEINCOMPAT>. If the <src> matrix is a packed type and
 *            <dest> is not, the values that are implicitly zeros are
 *            set to zeros in the <dest> matrix.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINCOMPAT> if <src>, <dest> are different sizes,
 *            or if their types differ and <dest> cannot represent
 *            <src>.
 */
int
esl_dmatrix_Copy(ESL_DMATRIX *src, ESL_DMATRIX *dest)
{
  int i,j;

  if (dest->n != src->n || dest->m != src->m)
    ESL_EXCEPTION(eslEINCOMPAT, "matrices of different size");

  if (src->type == dest->type)   /* simple case. */
    memcpy(dest->mx[0], src->mx[0], src->ncells * sizeof(double));

  else if (src->type == eslGENERAL && dest->type == eslUPPER)		
    {
      for (i = 1; i < src->n; i++)
	for (j = 0; j < i; j++)
	  if (src->mx[i][j] != 0.) 
	    ESL_EXCEPTION(eslEINCOMPAT, "general matrix isn't upper triangular, can't be copied/packed");
      for (i = 0; i < src->n; i++)
	for (j = i; j < src->m; j++)
	  dest->mx[i][j] = src->mx[i][j];
    }
  
  else if (src->type == eslUPPER && dest->type == eslGENERAL)		
    {
      for (i = 1; i < src->n; i++)
	for (j = 0; j < i; j++)
	  dest->mx[i][j] = 0.;
      for (i = 0; i < src->n; i++)
	for (j = i; j < src->m; j++)
	  dest->mx[i][j] = src->mx[i][j];      
    }

  return eslOK;
}


/* Function:  esl_dmatrix_Duplicate()
 * Incept:    SRE, Tue May  2 14:38:45 2006 [St. Louis]
 *
 * Purpose:   Duplicates matrix <A>, making a copy in newly
 *            allocated space.
 *
 * Returns:   a pointer to the copy. Caller frees with 
 *            <esl_dmatrix_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_DMATRIX *
esl_dmatrix_Duplicate(ESL_DMATRIX *A)
{
  ESL_DMATRIX *new;

  switch (A->type) {
  case eslUPPER:             if ( (new = esl_dmatrix_CreateUpper(A->n))  == NULL) return NULL; break;
  default: case eslGENERAL:  if ( (new = esl_dmatrix_Create(A->n, A->m)) == NULL) return NULL; break;
  }
  esl_dmatrix_Copy(A, new);
  return new;
}


/* Function:  esl_dmatrix_Compare()
 *
 * Purpose:   Compares matrix <A> to matrix <B> element by element,
 *            using <esl_DCompare()> on each cognate element pair, 
 *            with equality defined by a fractional tolerance <tol>.
 *            If all elements are equal, return <eslOK>; if any
 *            elements differ, return <eslFAIL>. 
 *            
 *            <A> and <B> may be of different types; for example,
 *            a packed upper triangular matrix A is compared to
 *            a general matrix B by assuming <A->mx[i][j] = 0.> for
 *            all $i>j$.
 */
int
esl_dmatrix_Compare(ESL_DMATRIX *A, ESL_DMATRIX *B, double tol)
{
  int i,j,c;
  double x1,x2;

  if (A->n != B->n) return eslFAIL;
  if (A->m != B->m) return eslFAIL;

  if (A->type == B->type) 
    {  /* simple case. */
      for (c = 0; c < A->ncells; c++) /* can deal w/ packed or unpacked storage */
	if (esl_DCompare(A->mx[0][c], B->mx[0][c], tol) == eslFAIL) return eslFAIL;
    }
  else 
    { /* comparing matrices of different types */
      for (i = 0; i < A->n; i++)
	for (j = 0; j < A->m; j++)
	  {
	    if (A->type == eslUPPER && i > j) x1 = 0.;
	    else                                         x1 = A->mx[i][j];

	    if (B->type == eslUPPER && i > j) x2 = 0.;
	    else                                         x2 = B->mx[i][j];

	    if (esl_DCompare(x1, x2, tol) == eslFAIL) return eslFAIL;
	  }
    }
  return eslOK;
}


/* Function:  esl_dmatrix_Set()
 *
 * Purpose:   Set all elements $a_{ij}$ in matrix <A> to <x>,
 *            and returns <eslOK>.
 */
int
esl_dmatrix_Set(ESL_DMATRIX *A, double x)
{
  int i;
  for (i = 0; i < A->ncells; i++) A->mx[0][i] = x;
  return eslOK;
}


/* Function:  esl_dmatrix_SetZero()
 *
 * Purpose:   Sets all elements $a_{ij}$ in matrix <A> to 0,
 *            and returns <eslOK>.
 */
int
esl_dmatrix_SetZero(ESL_DMATRIX *A)
{
  int i;
  for (i = 0; i < A->ncells; i++) A->mx[0][i] = 0.;
  return eslOK;
}
  

/* Function:  esl_dmatrix_SetIdentity()
 *
 * Purpose:   Given a square matrix <A>, sets all diagonal elements 
 *            $a_{ii}$ to 1, and all off-diagonal elements $a_{ij},
 *            j \ne i$ to 0. Returns <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the matrix isn't square.
 */
int
esl_dmatrix_SetIdentity(ESL_DMATRIX *A)
{
  int i;
  
  if (A->n != A->m) ESL_EXCEPTION(eslEINVAL, "matrix isn't square");
  esl_dmatrix_SetZero(A);
  for (i = 0; i < A->n; i++) A->mx[i][i] = 1.;
  return eslOK;
}
  


/* Function:  esl_dmx_Max()
 * Incept:    SRE, Thu Mar  1 14:46:48 2007 [Janelia]
 *
 * Purpose:   Returns the maximum value of all the elements $a_{ij}$ in matrix <A>.
 */
double
esl_dmx_Max(ESL_DMATRIX *A)
{
  int    i;
  double best;

  best = A->mx[0][0];
  for (i = 0; i < A->ncells; i++)
    if (A->mx[0][i] > best) best = A->mx[0][i];
  return best;
}

/* Function:  esl_dmx_Min()
 * Incept:    SRE, Thu Mar  1 14:49:29 2007 [Janelia]
 *
 * Purpose:   Returns the minimum value of all the elements $a_{ij}$ in matrix <A>.
 */
double
esl_dmx_Min(ESL_DMATRIX *A)
{
  int    i;
  double best;

  best = A->mx[0][0];
  for (i = 0; i < A->ncells; i++)
    if (A->mx[0][i] < best) best = A->mx[0][i];
  return best;
}


/* Function:  esl_dmx_Sum()
 * Incept:    SRE, Thu Mar  1 16:45:16 2007
 *
 * Purpose:   Returns the scalar sum of all the elements $a_{ij}$ in matrix <A>,
 *            $\sum_{ij} a_{ij}$.
 */
double
esl_dmx_Sum(ESL_DMATRIX *A)
{
  int    i;
  double sum = 0.;

  for (i = 0; i < A->ncells; i++)
    sum += A->mx[0][i];
  return sum;
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
 * Returns:   a pointer to a new <ESL_PERMUTATION> object. Free with 
 *            <esl_permutation_Destroy()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
ESL_PERMUTATION *
esl_permutation_Create(int n)
{
  int status;
  ESL_PERMUTATION *P = NULL;

  ESL_ALLOC(P, sizeof(ESL_PERMUTATION));
  P->pi = NULL;
  P->n  = n;
  ESL_ALLOC(P->pi, sizeof(int) * n);

  esl_permutation_Reuse(P);	/* initialize it */
  return P;

 ERROR:
  esl_permutation_Destroy(P);
  return NULL;
}
  
/* Function:  esl_permutation_Destroy()
 *
 * Purpose:   Frees an <ESL_PERMUTATION> object <P>.
 */
int
esl_permutation_Destroy(ESL_PERMUTATION *P)
{
  if (P != NULL && P->pi != NULL) free(P->pi);
  if (P != NULL)                  free(P);
  return eslOK;
}

/* Function:  esl_permutation_Reuse()
 *
 * Purpose:   Resets a permutation matrix <P> to
 *            $p_i = i$ for all $i = 0..n-1$.
 *            
 * Returns:   <eslOK> on success.           
 */
int
esl_permutation_Reuse(ESL_PERMUTATION *P)
{
  int i;
  for (i = 0; i < P->n; i++)
    P->pi[i] = i;
  return eslOK;
}


/* Function:  esl_permutation_Dump()
 *
 * Purpose:   Given a permutation matrix <P>, dump it to output stream <ofp>
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
 * Returns:   <eslOK> on success.
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
  return eslOK;
}



/* Function: esl_dmx_Multiply()
 * 
 * Purpose:  Matrix multiplication: calculate <AB>, store result in <C>.
 *           <A> is $n times m$; <B> is $m \times p$; <C> is $n \times p$.
 *           Matrix <C> must be allocated appropriately by the caller.
 *
 *           Not supported for anything but general (<eslGENERAL>)
 *           matrix type, at present.
 *           
 * Throws:   <eslEINVAL> if matrices don't have compatible dimensions,
 *           or if any of them isn't a general (<eslGENERAL>) matrix.
 */
int
esl_dmx_Multiply(ESL_DMATRIX *A, ESL_DMATRIX *B, ESL_DMATRIX *C)
{
  int i, j, k;

  if (A->m    != B->n)       ESL_EXCEPTION(eslEINVAL, "can't multiply A,B");
  if (A->n    != C->n)       ESL_EXCEPTION(eslEINVAL, "A,C # of rows not equal");
  if (B->m    != C->m)       ESL_EXCEPTION(eslEINVAL, "B,C # of cols not equal");
  if (A->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "A isn't of type eslGENERAL");
  if (B->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "B isn't of type eslGENERAL");
  if (C->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "B isn't of type eslGENERAL");

  for (i = 0; i < A->n; i++)
    for (j = 0; j < B->m; j++)
      {
	C->mx[i][j] = 0.;
	for (k = 0; k < A->m; k++)
	  C->mx[i][j] += A->mx[i][k] * B->mx[k][j];
      }
  return eslOK;
}


/* Function:  esl_dmx_Transpose()
 *
 * Purpose:   Transpose a square matrix <A> in place.
 *
 *            <A> must be a general (<eslGENERAL>) matrix type.
 *
 * Throws:    <eslEINVAL> if <A> isn't square, or if it isn't
 *            of type <eslGENERAL>.
 */
int
esl_dmx_Transpose(ESL_DMATRIX *A)
{
  int    i,j;
  double swap;

  if (A->n    != A->m)       ESL_EXCEPTION(eslEINVAL, "matrix isn't square");
  if (A->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "A isn't of type eslGENERAL");

  for (i = 0; i < A->n; i++)
    for (j = i+1; j < A->m; j++)
      { swap = A->mx[i][j]; A->mx[i][j] = A->mx[j][i]; A->mx[j][i] = swap; }
  return eslOK;
}


/* Function:  esl_dmx_Add()
 *
 * Purpose:   <A = A+B>; adds matrix <B> to matrix <A> and leaves result
 *            in matrix <A>.
 *
 *            <A> and <B> may be of any type. However, if <A> is a
 *            packed upper triangular matrix (type
 *            <eslUPPER>), all values $i>j$ in <B> must be
 *            zero (i.e. <B> must also be upper triangular, though
 *            not necessarily packed upper triangular).
 *
 * Throws:    <eslEINVAL> if matrices aren't the same dimensions, or
 *            if <A> is <eslUPPER> and any cell $i>j$ in
 *            <B> is nonzero.
 */
int
esl_dmx_Add(ESL_DMATRIX *A, ESL_DMATRIX *B)
{
  int    i,j;
  
  if (A->n    != B->n)              ESL_EXCEPTION(eslEINVAL, "matrices of different size");
  if (A->m    != B->m)              ESL_EXCEPTION(eslEINVAL, "matrices of different size");

  if (A->type == B->type)	/* in this case, can just add cell by cell */
    {
      for (i = 0; i < A->ncells; i++)
	A->mx[0][i] += B->mx[0][i];
    }
  else if (A->type == eslUPPER || B->type == eslUPPER)
    {
      /* Logic is: if either matrix is upper triangular, then the operation is
       * to add upper triangles only. If we try to add a general matrix <B>
       * to packed UT <A>, make sure all lower triangle entries in <B> are zero.
       */
      if (B->type != eslUPPER) {
	for (i = 1; i < A->n; i++)
	  for (j = 0; j < i; j++)
	    if (B->mx[i][j] != 0.) ESL_EXCEPTION(eslEINVAL, "<B> has nonzero cells in lower triangle");
      }
      for (i = 0; i < A->n; i++)
	for (j = i; j < A->m; j++)
	  A->mx[i][j] += B->mx[i][j];
    }
  return eslOK;
}

/* Function:  esl_dmx_Scale()
 *
 * Purpose:   Calculates <A = kA>: multiply matrix <A> by scalar
 *            <k> and leave answer in <A>.
 */
int 
esl_dmx_Scale(ESL_DMATRIX *A, double k)
{
  int i;

  for (i = 0; i < A->ncells; i++)  A->mx[0][i] *=  k;
  return eslOK;
}


/* Function:  esl_dmx_AddScale()
 * 
 * Purpose:   Calculates <A + kB>, leaves answer in <A>.
 * 
 *            Only defined for matrices of the same type (<eslGENERAL>
 *            or <eslUPPER>).
 * 
 * Throws:    <eslEINVAL> if matrices aren't the same dimensions, or
 *            of different types.
 */
int
esl_dmx_AddScale(ESL_DMATRIX *A, double k, ESL_DMATRIX *B)
{
  int i;

  if (A->n    != B->n)    ESL_EXCEPTION(eslEINVAL, "matrices of different size");
  if (A->m    != B->m)    ESL_EXCEPTION(eslEINVAL, "matrices of different size");
  if (A->type != A->type) ESL_EXCEPTION(eslEINVAL, "matrices of different type");

  for (i = 0; i < A->ncells; i++) A->mx[0][i] +=  k * B->mx[0][i];
  return eslOK;
}


/* Function:  esl_dmx_Permute_PA()
 *
 * Purpose:   Computes <B = PA>: do a row-wise permutation of a square
 *            matrix <A>, using the permutation matrix <P>, and put
 *            the result in a square matrix <B> that the caller has
 *            allocated.
 *
 * Throws:    <eslEINVAL> if <A>, <B>, <P> do not have compatible dimensions,
 *            or if <A> or <B> is not of type <eslGENERAL>.
 */
int
esl_dmx_Permute_PA(ESL_PERMUTATION *P, ESL_DMATRIX *A, ESL_DMATRIX *B)
{
  int i,ip,j;

  if (A->n    != P->n)       ESL_EXCEPTION(eslEINVAL, "matrix dimensions not compatible");
  if (A->n    != B->n)       ESL_EXCEPTION(eslEINVAL, "matrix dimensions not compatible");
  if (A->n    != A->m)       ESL_EXCEPTION(eslEINVAL, "matrix dimensions not compatible");
  if (B->n    != B->m)       ESL_EXCEPTION(eslEINVAL, "matrix dimensions not compatible");
  if (A->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "matrix A not of type eslGENERAL");
  if (B->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "matrix B not of type eslGENERAL");

  for (i = 0; i < A->n; i++)
    {
      ip = P->pi[i];
      for (j = 0; j < A->m; j++)
	B->mx[i][j] = A->mx[ip][j];
    }
  return eslOK;
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
 *            Implements Gaussian elimination with pivoting 
 *            \citep[p.~759]{Cormen99}.
 *
 * Throws:    <eslEINVAL> if <A> isn't square, or if <P> isn't the right
 *            size for <A>, or if <A> isn't of general type.
 */
int
esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P)
{
  int    i,j,k,kpiv;
  double max;
  double swap;

  if (A->n    != A->m)       ESL_EXCEPTION(eslEINVAL, "matrix isn't square");
  if (P->n    != A->n)       ESL_EXCEPTION(eslEINVAL, "permutation isn't the right size");
  if (A->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "matrix isn't of general type");
  esl_permutation_Reuse(P);

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
      if (max == 0.) ESL_EXCEPTION(eslEDIVZERO, "matrix is singular");
      
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
  return eslOK;
}


/* Function:  esl_dmx_LU_separate()
 *
 * Purpose:   Separate a square <LU> decomposition matrix into its two
 *            triangular matrices <L> and <U>. Caller provides two
 *            allocated <L> and <U> matrices of same size as <LU> for
 *            storing the results.
 *            
 *            <U> may be an upper triangular matrix in either unpacked
 *            (<eslGENERAL>) or packed (<eslUPPER>) form.
 *            <LU> and <L> must be of <eslGENERAL> type.
 *
 * Throws:    <eslEINVAL> if <LU>, <L>, <U> are not of compatible dimensions,
 *            or if <LU> or <L> aren't of general type. 
 */
int
esl_dmx_LU_separate(ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U)
{
  int i,j;

  if (LU->n    != LU->m)      ESL_EXCEPTION(eslEINVAL, "LU isn't square");
  if (L->n     != L->m)       ESL_EXCEPTION(eslEINVAL, "L isn't square");
  if (U->n     != U->m)       ESL_EXCEPTION(eslEINVAL, "U isn't square");
  if (LU->n    != L->n)       ESL_EXCEPTION(eslEINVAL, "LU, L have incompatible dimensions");
  if (LU->n    != U->n)       ESL_EXCEPTION(eslEINVAL, "LU, U have incompatible dimensions");
  if (LU->type != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "matrix isn't of general type");
  if (L->type  != eslGENERAL) ESL_EXCEPTION(eslEINVAL, "matrix isn't of general type");

  esl_dmatrix_SetZero(L);
  esl_dmatrix_SetZero(U);

  for (i = 0; i < LU->n; i++)
    for (j = i; j < LU->m; j++)
      U->mx[i][j] = LU->mx[i][j];

  for (i = 0; i < LU->n; i++) 
    {
      L->mx[i][i] = 1.;
      for (j = 0; j < i; j++)
	L->mx[i][j] = LU->mx[i][j];
    }
  return eslOK;
}

/* Function:  esl_dmx_Invert()
 *
 * Purpose:   Calculates the inverse of square matrix <A>, and stores the
 *            result in matrix <Ai>. Caller provides an allocated
 *            matrix <Ai> of same dimensions as <A>. Both must be
 *            of type <eslGENERAL>.
 *            
 *            Peforms the inversion by LUP decomposition followed by 
 *            forward/back-substitution \citep[p.~753]{Cormen99}.
 *
 * Throws:    <eslEINVAL> if <A>, <Ai> do not have same dimensions, 
 *                        if <A> isn't square, or if either isn't of
 *                        type <eslGENERAL>.
 *            <eslEMEM>   if internal allocations (for LU, and some other
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
  int               status;

  if (A->n     != A->m)                   ESL_EXCEPTION(eslEINVAL, "matrix isn't square");
  if (A->n     != Ai->n || A->m != Ai->m) ESL_EXCEPTION(eslEINVAL, "matrices are different size");
  if (A->type  != eslGENERAL)             ESL_EXCEPTION(eslEINVAL, "matrix A not of general type");
  if (Ai->type != eslGENERAL)             ESL_EXCEPTION(eslEINVAL, "matrix B not of general type");

  /* Copy A to LU, and do an LU decomposition.
   */
  if ((LU = esl_dmatrix_Create(A->n, A->m)) == NULL) goto ERROR;
  if ((P  = esl_permutation_Create(A->n))   == NULL) goto ERROR;
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
  ESL_ALLOC(b, sizeof(double) * A->n);
  ESL_ALLOC(y, sizeof(double) * A->n);

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
  return eslOK;

 ERROR:
  if (y  != NULL) free(y);
  if (b  != NULL) free(b);
  if (LU != NULL) esl_dmatrix_Destroy(LU);
  if (P  != NULL) esl_permutation_Destroy(P);
  return status;
}


/*****************************************************************
 * Example code
 *****************************************************************/ 

/*   gcc -g -Wall -o example -I. -DeslDMATRIX_EXAMPLE esl_dmatrix.c easel.c -lm
 */
#ifdef eslDMATRIX_EXAMPLE
/*::cexcerpt::dmatrix_example::begin::*/
#include <easel.h>
#include <esl_dmatrix.h>

int main(void)
{
  ESL_DMATRIX *A, *B, *C;

  A = esl_dmatrix_Create(4,4);
  B = esl_dmatrix_Create(4,4);
  C = esl_dmatrix_Create(4,4);
  
  esl_dmatrix_SetIdentity(A);
  esl_dmatrix_Copy(A, B);

  esl_dmx_Multiply(A,B,C);

  esl_dmatrix_Dump(stdout, C, NULL, NULL);

  esl_dmatrix_Destroy(A);
  esl_dmatrix_Destroy(B);
  esl_dmatrix_Destroy(C);
  return 0;
}
/*::cexcerpt::dmatrix_example::end::*/
#endif /*eslDMATRIX_EXAMPLE*/


/*****************************************************************
 * Test driver
 *****************************************************************/ 

/*   gcc -g -Wall -o testdriver -I. -L. -DeslDMATRIX_TESTDRIVE esl_dmatrix.c -leasel -lm
 */
#ifdef eslDMATRIX_TESTDRIVE
#include <easel.h>
#include <esl_dmatrix.h>

int main(void)
{
  ESL_DMATRIX *A, *B, *C;

  A = esl_dmatrix_Create(4,4);
  B = esl_dmatrix_Create(4,4);
  C = esl_dmatrix_Create(4,4);
  
  esl_dmatrix_SetIdentity(A);   /* A=I */
  esl_dmx_Invert(A, B);		/* B=I-1=I */
  esl_dmx_Multiply(A,B,C);	/* C=I */
  esl_dmx_Transpose(A);         /* A=I still */

  esl_dmx_Scale(A, 0.5);	/* A= 0.5I */
  esl_dmx_AddScale(B, -0.5, C);	/* B= 0.5I */
  
  esl_dmx_Add(A,B);		/* A=I */
  esl_dmx_Scale(B, 2.0);	/* B=I */

  if (esl_dmatrix_Compare(A, B, 1e-6) != eslOK) esl_fatal("A != B");
  if (esl_dmatrix_Compare(A, C, 1e-6) != eslOK) esl_fatal("A != C");
  esl_dmatrix_Copy(B, C);
  if (esl_dmatrix_Compare(A, C, 1e-6) != eslOK) esl_fatal("A != copied B");    

  esl_dmatrix_Destroy(A);
  esl_dmatrix_Destroy(B);
  esl_dmatrix_Destroy(C);
  return 0;
}
#endif /*eslDMATRIX_TESTDRIVE*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
