/* esl_vectorops.h
 * Header file for vectorops.c
 * 
 * SRE, Tue Oct  1 15:23:37 2002 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_VECTOROPS_INCLUDED
#define ESL_VECTOROPS_INCLUDED


extern void   esl_vec_DSet(double *vec, int n, double value);
extern void   esl_vec_FSet(float *vec, int n, float value);
extern void   esl_vec_ISet(int *vec, int n, int value);

extern void   esl_vec_DScale(double *vec, int n, double scale);
extern void   esl_vec_FScale(float *vec, int n, float scale);
extern void   esl_vec_IScale(int *vec, int n, int scale);

extern double esl_vec_DSum(double *vec, int n);
extern float  esl_vec_FSum(float *vec, int n);
extern int    esl_vec_ISum(int *vec, int n);

extern void   esl_vec_DAdd(double *vec1, double *vec2, int n);
extern void   esl_vec_FAdd(float *vec1, float *vec2, int n);
extern void   esl_vec_IAdd(int *vec1, int *vec2, int n);

extern void   esl_vec_DCopy(double *vec1, double *vec2, int n);
extern void   esl_vec_FCopy(float *vec1, float *vec2, int n);
extern void   esl_vec_ICopy(int *vec1, int *vec2, int n);

extern double esl_vec_DDot(double *vec1, double *vec2, int n);
extern float  esl_vec_FDot(float *vec1, float *vec2, int n);
extern int    esl_vec_IDot(int *vec1, int *vec2, int n);

extern double esl_vec_DMax(double *vec, int n);
extern float  esl_vec_FMax(float *vec, int n);
extern int    esl_vec_IMax(int *vec, int n);

extern double esl_vec_DMin(double *vec, int n);
extern float  esl_vec_FMin(float *vec, int n);
extern int    esl_vec_IMin(int *vec, int n);

extern int    esl_vec_DArgMax(double *vec, int n);
extern int    esl_vec_FArgMax(float *vec, int n);
extern int    esl_vec_IArgMax(int *vec, int n);

extern int    esl_vec_DArgMin(double *vec, int n);
extern int    esl_vec_FArgMin(float *vec, int n);
extern int    esl_vec_IArgMin(int *vec, int n);

extern void   esl_vec_DNorm(double *vec, int n);
extern void   esl_vec_FNorm(float *vec, int n);

extern void   esl_vec_DLog(double *vec, int n);
extern void   esl_vec_FLog(float *vec, int n);

extern double esl_vec_DEntropy(double *p, int n);
extern float  esl_vec_FEntropy(float *p, int n);

extern void   esl_vec_DExp(double *vec, int n);
extern void   esl_vec_FExp(float *vec, int n);

extern double esl_vec_DLogSum(double *vec, int n);
extern float  esl_vec_FLogSum(float *vec, int n);


#endif /* ESL_VECTOROPS_INCLUDED */
