/* vectorops.c
 * Operations on vectors of floats or doubles.
 * 
 * Can operate on vectors of doubles, floats, or integers - appropriate
 * routine is prefixed with D, F, or I. For example, esl_vec_DSet() is
 * the Set routine for a vector of doubles; esl_vec_ISet() is for integers.
 * 
 * esl_vec_{D,F,I}Set()         - set all items in vector to value.
 * esl_vec_{D,F,I}Scale()       - multiply all items in vector by scale
 * esl_vec_{D,F,I}Sum()         - return sum of values in vector
 * esl_vec_{D,F,I}Add()         - add vec2 to vec1.
 * esl_vec_{D,F,I}Copy()        - set vec1 to be same as vec2. 
 * esl_vec_{D,F,I}Dot()         - return dot product of two vectors.
 * esl_vec_{D,F,I}Max()         - return value of maximum element in vector
 * esl_vec_{D,F,I}Min()         - return value of minimum element in vector 
 * esl_vec_{D,F,I}ArgMax()      - return index of maximum element in vector
 * esl_vec_{D,F,I}ArgMin()      - return index of minimum element in vector
 * 
 * esl_vec_{D,F}Norm()          - normalize a probability vector of length n.
 * esl_vec_{D,F}Log()           - convert to log probabilities 
 * esl_vec_{D,F}Exp()           - convert log p's back to probabilities
 * esl_vec_{D,F}LogSum()        - given vector of log p's; return log of summed p's.
 * esl_vec_{D,F}Entropy()       - return Shannon entropy of probability vector, in bits
 *                        
 * SRE, Tue Oct  1 15:23:25 2002 [St. Louis]
 * SVN $Id$
 */                      
  

#include <math.h>
#include <float.h>
#include <easel/easel.h>
#include <easel/vectorops.h>

void
esl_vec_DSet(double *vec, int n, double value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}
void
esl_vec_FSet(float *vec, int n, float value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}
void
esl_vec_ISet(int *vec, int n, int value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}


void
esl_vec_DScale(double *vec, int n, double scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}
void
esl_vec_FScale(float *vec, int n, float scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}
void
esl_vec_IScale(int *vec, int n, int scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}


double 
esl_vec_DSum(double *vec, int n)
{
  double sum = 0.;
  int    x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}
float 
esl_vec_FSum(float *vec, int n)
{
  float sum = 0.;
  int   x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}
int
esl_vec_ISum(int *vec, int n)
{
  int sum = 0;
  int   x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}


void
esl_vec_DAdd(double *vec1, double *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}
void
esl_vec_FAdd(float *vec1, float *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}
void
esl_vec_IAdd(int *vec1, int *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}


void
esl_vec_DCopy(double *vec1, double *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] = vec2[x];
}
void
esl_vec_FCopy(float *vec1, float *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] = vec2[x];
}
void
esl_vec_ICopy(int *vec1, int *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] = vec2[x];
}


double
esl_vec_DDot(double *vec1, double *vec2, int n)
{
  double result = 0.;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}
float
esl_vec_FDot(float *vec1, float *vec2, int n)
{
  float result = 0.;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}
int
esl_vec_IDot(int *vec1, int *vec2, int n)
{
  int result = 0;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}

double
esl_vec_DMax(double *vec, int n)
{
  int i;
  double best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}
float
esl_vec_FMax(float *vec, int n)
{
  int   i;
  float best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}
int
esl_vec_IMax(int *vec, int n)
{
  int   i;
  int   best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}


double
esl_vec_DMin(double *vec, int n)
{
  int i;
  double best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}
float
esl_vec_FMin(float *vec, int n)
{
  int   i;
  float best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}
int
esl_vec_IMin(int *vec, int n)
{
  int   i;
  int   best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}


int
esl_vec_DArgMax(double *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}
int
esl_vec_FArgMax(float *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}
int
esl_vec_IArgMax(int *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}


int
esl_vec_DArgMin(double *vec, int n)
{
  int i;
  int best = 0;
  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}
int
esl_vec_FArgMin(float *vec, int n)
{
  int   i;
  int   best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}
int
esl_vec_IArgMin(int *vec, int n)
{
  int   i;
  int   best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}


void
esl_vec_DNorm(double *vec, int n)
{
  int    x;
  double sum;

  sum = esl_vec_DSum(vec, n);
  if (sum != 0.0) for (x = 0; x < n; x++) vec[x] /= sum;
  else            for (x = 0; x < n; x++) vec[x] = 1. / (double) n;
}
void
esl_vec_FNorm(float *vec, int n)
{
  int    x;
  float  sum;

  sum = esl_vec_FSum(vec, n);
  if (sum != 0.0) for (x = 0; x < n; x++) vec[x] /= sum;
  else            for (x = 0; x < n; x++) vec[x] = 1. / (float) n;
}


void
esl_vec_DLog(double *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) 
    if (vec[x] > 0.) vec[x] = log(vec[x]);
    else vec[x] = -DBL_MAX;
}
void
esl_vec_FLog(float *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) 
    if (vec[x] > 0.) vec[x] = log(vec[x]);
    else vec[x] = -FLT_MAX;
}


void
esl_vec_DExp(double *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) vec[x] = exp(vec[x]);
}
void
esl_vec_FExp(float *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) vec[x] = exp(vec[x]);
}

double
esl_vec_DLogSum(double *vec, int n)
{
  int x;
  double max, sum;
  
  max = esl_vec_DMax(vec, n);
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      sum += exp(vec[x] - max);
  sum = log(sum) + max;
  return sum;
}
float
esl_vec_FLogSum(float *vec, int n)
{
  int x;
  float max, sum;
  
  max = esl_vec_FMax(vec, n);
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      sum += exp(vec[x] - max);
  sum = log(sum) + max;
  return sum;
}


double
esl_vec_DEntropy(double *p, int n)
{
  int    i;
  double entropy;

  entropy = 0.;
  for(i = 0; i < n; i++)
    if (p[i] > 0.) entropy += p[i] * log(p[i]);
  return(-1.44269504 * entropy); /* converts to bits */
}
float
esl_vec_FEntropy(float *p, int n)
{
  int    i;
  float  entropy;

  entropy = 0.;
  for(i = 0; i < n; i++)
    if (p[i] > 0.) entropy += p[i] * log(p[i]);
  return(-1.44269504 * entropy); /* converts to bits */
}
