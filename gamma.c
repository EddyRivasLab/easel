/* gamma.c
 * Functions relative to the Gamma function and Gamma densities.
 * 
 * SRE, Tue Nov  2 13:46:12 2004 [St. Louis]
 * SVN $Id$
 */

#include <math.h>

#include <easel/easel.h>
#include <easel/random.h>
#include <easel/gamma.h>

/* These three subfunctions are used by esl_gamma_Sample().
 */
static double gamma_ahrens(ESL_RANDOMNESS *r, double a);
static double gamma_integer(ESL_RANDOMNESS *r, unsigned int a);
static double gamma_fraction(ESL_RANDOMNESS *r, double a);

/* Function:  esl_gamma_log()
 * Incept:    SRE, Tue Nov  2 13:47:01 2004 [St. Louis]
 *
 * Purpose:   Returns natural log of Gamma(x), for x > 0.0.
 * 
 * Credit:    Adapted from a public domain implementation in the
 *            NCBI core math library. Thanks to John Spouge and
 *            the NCBI. (According to NCBI, that's Dr. John
 *            "Gammas Galore" Spouge to you, pal.)
 *
 * Args:      x          : argument, x > 0.0
 *            ret_answer : RETURN: the answer
 *
 * Returns:   Put the answer in *ret_answer;
 *            returns eslOK.
 *            
 *            returns eslEINVAL if x <= 0.
 */
int
esl_gamma_log(double x, double *ret_answer)
{
  int i;
  double xx, tx;
  double tmp, value;
  static double cof[11] = {
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10
  };
  
  /* Protect against invalid x<=0
   */
  if (x <= 0.0)  ESL_ERROR(eslEINVAL, "invalid x <= 0 in esl_gamma_LogGamma()");

  xx       = x - 1.0;
  tx = tmp = xx + 11.0;
  value    = 1.0;
  for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
      value += cof[i] / tmp;
      tmp   -= 1.0;
    }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  *ret_answer = value;
  return eslOK;
}


/* Function: esl_gamma_sample()
 * Date:     SRE, Wed Apr 17 13:10:03 2002 [St. Louis]
 *
 * Purpose:  Return a random deviate distributed as Gamma(a, 1.).
 *           
 *           Follows Knuth, vol. 2 Seminumerical Algorithms, pp.133-134.
 *           Also relies on examination of the implementation in
 *           the GNU Scientific Library (libgsl). The implementation
 *           relies on three separate gamma function algorithms:
 *           gamma_ahrens(), gamma_integer(), and gamma_fraction().
 *
 * Args:     r          - random number generation seed
 *           alpha      - order of the gamma function
 *           ret_answer - RETURN: a sample from Gamma(a, 1).
 *
 * Returns:  a gamma-distributed deviate is put in *ret_answer;
 *           returns eslOK.
 *
 *           Returns eslEINVAL for a <= 0.
 */
int
esl_gamma_sample(ESL_RANDOMNESS *r, double a, double *ret_answer)
{
  double aint;

  if (a <= 0.) ESL_ERROR(eslEINVAL, "a <= 0 in esl_gamma_Sample()");

  aint = floor(a);
  if (a == aint && a < 12.) 
    *ret_answer = gamma_integer(r, (unsigned int) a);
  else if (a > 3.) 
    *ret_answer = gamma_ahrens(r, a);
  else if (a < 1.) 
    *ret_answer = gamma_fraction(r, a);
  else 
    *ret_answer = gamma_ahrens(r, aint) + gamma_fraction(r, a-aint);
  return eslOK;
}

static double
gamma_ahrens(ESL_RANDOMNESS *r, double a)	/* for a >= 3 */
{
  double V;			/* uniform deviates */
  double X,Y;
  double test;
  
  do {
    do {				/* generate candidate X */
      Y = tan(eslCONSTANT_PI * esl_random(r)); 
      X = Y * sqrt(2.*a -1.) + a - 1.;
    } while (X <= 0.);
				/* accept/reject X */
    V    = esl_random(r);
    test = (1+Y*Y) * exp( (a-1.)* log(X/(a-1.)) - Y*sqrt(2.*a-1.));
  } while (V > test);
  return X;
}
static double
gamma_integer(ESL_RANDOMNESS *r, unsigned int a)	/* for small integer a, a < 12 */
{
  int    i;
  double U,X;

  U = 1.;
  for (i = 0; i < a; i++) 
    U *= esl_rnd_UniformPositive(r);
  X = -log(U);

  return X;
}
static double
gamma_fraction(ESL_RANDOMNESS *r, double a)	/* for fractional a, 0 < a < 1 */
{				/* Knuth 3.4.1, exercise 16, pp. 586-587 */
  double p, U, V, X, q;
  
  p = eslCONSTANT_E / (a + eslCONSTANT_E);
  do {
    U = esl_random(r);
    V = esl_rnd_UniformPositive(r);
    if (U < p) {
      X = pow(V, 1./a);
      q = exp(-X);
    } else {
      X = 1. - log(V);
      q = pow(X, a-1.);
    }
    U = esl_random(r);
  } while (U >= q);
  return X;
}

