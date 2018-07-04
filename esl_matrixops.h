/* Simple matrix operations.
 * 
 * Compare:
 *    esl_vectorops: simple vector operations.
 *    esl_dmatrix:   matrix algebra, double precision
 */
#ifndef eslMATRIXOPS_INCLUDED
#define eslMATRIXOPS_INCLUDED
#include "esl_config.h"

extern double **esl_mat_DCreate(int M, int N);
extern float  **esl_mat_FCreate(int M, int N);
extern int    **esl_mat_ICreate(int M, int N);

extern void     esl_mat_DSet(double **A, int M, int N, double value);
extern void     esl_mat_FSet(float  **A, int M, int N, float  value);
extern void     esl_mat_ISet(int    **A, int M, int N, int    value);

extern void     esl_mat_DCopy(double **src, int M, int N, double **dest);
extern void     esl_mat_FCopy(float  **src, int M, int N, float  **dest);
extern void     esl_mat_ICopy(int    **src, int M, int N, int    **dest);

extern double   esl_mat_DMax(double **A, int M, int N);
extern float    esl_mat_FMax(float  **A, int M, int N);
extern int      esl_mat_IMax(int    **A, int M, int N);

extern int      esl_mat_DCompare(double **A, double **B, int M, int N, double tol);
extern int      esl_mat_FCompare(float  **A, float  **B, int M, int N, float  tol);
extern int      esl_mat_ICompare(int    **A, int    **B, int M, int N);

extern void     esl_mat_DDestroy(double **A, int M);
extern void     esl_mat_FDestroy(float  **A, int M);
extern void     esl_mat_IDestroy(int    **A, int M);

extern int      esl_mat_IDump(int **A, int M, int N);

#endif // eslMATRIXOPS_INCLUDED
