# esl_dmatrix - 2D matrices and linear algebra

The `dmatrix` module implements 2D matrices and linear algebra
operations.

There are two objects. The main one is a `ESL_DMATRIX`, a 2D
real-valued matrix of n rows and m columns. There is also
`ESL_PERMUTATION`, a special matrix used in LU decompositions.

It is straightforward to call standard BLAS and LAPACK linear algebra
routines on the data in an `ESL_DMATRIX`.

## accessing matrix values

The accessible internals of the `ESL_DMATRIX` structure are:

```c
  double **mx;                  /* mx[i][j] is i'th row, j'th col */
  int      n;                   /* rows    */
  int      m;                   /* columns */
  enum { eslGENERAL, eslUPPER } type;
```

The matrix is stored in row-major orientation: the value in cell
$(i,j)$ in row $i$ and column $j$ is in `mx->mx[i][j]`.

Elements are stored in a single array `mx->mx[0]`. This is important
for interoperability with BLAS and LAPACK; see below. The row pointers
`mx->mx[i]` are initialized so that elements may be accessed simply as
`mx->mx[i][j]`, rather than by pointer arithmetic 
`mx->mx[0] + i*mx->m + j`.


## specialized matrix types

Normally matrices are created with `esl_dmatrix_Create()`, which
allocates storage for all $n \times m$ cells. Easel calls this a
matrix of type `eslGENERAL`.

Matrices may have more restricted forms, which may constrain certain
values and may allow packed storage. For example, an upper triangular
matrix is one in which all elements $i>j$ have a value of zero. When we
calculate the minimum in such a matrix with `esl_dmx_Min()`, we
probably don't want to consider the $i>j$ elements. We also can save
almost two-fold in storage by not storing the $i>j$ elements at all.
Other types include square, lower triagonal, and symmetric matrices.

We expect to need to expand Easel's implementation of different matrix
types in the future, but right now, Easel has just one other matrix
type, `eslUPPER`, for packed upper triangular matrices.


### `eslUPPER`: packed upper triangular matrices

An `eslUPPER` matrix is created with `esl_dmatrix_CreateUpper(int n)`.
It is necessarily square $n \times n$, so only one dimension argument
is passed. Most but not all functions in `dmatrix` can operate on
`eslUPPER` matrix types in addition to the usual `eslGENERAL` type.

The caller must not access any cell $i>j$ in an `eslUPPER` matrix.
Setting a cell $i>j$ will corrupt the matrix. Accessing cell $i>j$ will
return an incorrect value, not zero.

The $n (n+1) / 2$ elements of the upper triagonal matrix are packed
into an array `mx->mx[0]`. You can access element $i,j$ by pointer
arithmetic at `mx->mx[j + i(2*mx->m-i-1)/2]` if you like, but it is
easier to access element $i,j$ by the usual `mx->mx[i][j]`. This is
made possible because the row pointers `mx->mx[i]` in an `eslUPPER`
matrix are tricksily initialized in an overlapping fashion so that
`mx->mx[i][j]` does the right thing for $i \leq j$. This overlapping is
also the reason why `mx->mx[i][j]` accesses the wrong element when
$i>j$.


### notes on the current implementation of matrix types

Easel matrix types conflate packing and element validity together. For
example, an upper triangular matrix may be stored either in an
`eslGENERAL` matrix type (in which case elements $i>j$ are set to zero)
or the packed `eslUPPER` matrix type (in which case elements $i>j$
aren't even stored). Using the `eslUPPER` matrix type is 2x more space
efficient, and also, operations like `esl_dmx_Min()` and
`esl_dmx_Max()` will examine all elements in an `eslGENERAL` matrix
(including the zeros), but only the elements $i \leq j$ in a `eslUPPER`
matrix.

This design is provisional. We may adopt a system more closely akin to
BLAS/LAPACK in the future, which distinguish between matrix type and
matrix storage. For example, BLAS has matrices of form `TR` and `TP`
for triangular and packed triangular. Easel's implementation seems
sufficient for the moment, and should also extend to lower diagonal and
symmetric matrices without difficulty when and if they become needed.
In any future development, look to BLAS and LAPACK for guidance.


## interoperability with BLAS and LAPACK

The BLAS and LAPACK libraries provide optimized, standardized linear
algebra routines. The storage in `ESL_DMATRIX` is designed so you can
call routines in these libraries. The `mx->mx[0]` array is a valid
matrix for BLAS and LAPACK so long as you know the right incantations.
These are summarized here:

| Easel type   | `CBLAS_ORDER`    | stride   | `CBLAS_UPLO`  | type   | code                    |
|--------------|------------------|----------|---------------|--------|-------------------------|
| `eslGENERAL` | `CblasRowMajor`  | `mx->m`  | n/a           | double | `GE` (GEneral)          |
| `eslUPPER`   | `CBlasRowMajor`  | `mx->m`  | `CblasUpper`  | double | `TP` (Triangular Packed) |

For example, to call the CBLAS (C implementation of BLAS) for an
operation on an Easel matrix of type `eslGENERAL`, you look for a
routine that starts with prefix `cblas_dge*` (`d` for double, `ge` for
general). An example is `cblas_dgemm()`, the matrix multiplication
(`mm`) routine, which computes $C = \alpha \mathit{op}(A)
\mathit{op}(B) + \beta C$ for matrices $A,B,C$ and scalars
$\alpha,\beta$, where $\mathit{op}(A)$ means $A$, $A^T$ (the
transpose), or $A^H$ (the conjugate transpose). $\mathit{op}(A)$ is an
$M \times K$ matrix, $\mathit{op}(B)$ is $K \times N$ matrix, and the
result $C$ is $M \times N$. The prototype for `cblas_dgemm` is:

```c
    void
    cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta, double *C,
                 const int ldc)
```

The `Order` argument is always `CblasRowMajor` for Easel matrices. The
`TransA` and `TransB` arguments specify $\mathit{op}()$: `CblasNoTrans`
means just the matrix itself. The `ld*` arguments are the major strides
for each matrix: the number of elements in each row, for our row-major
matrices. So, we could call:

```c
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                 A->n, B->m, A->m,
		 1.0, A->mx[0], A->m,
		 B->mx[0], B->m,
		 1.0, C->mx[0], C->m);
```
