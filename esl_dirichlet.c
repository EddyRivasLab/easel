/* esl_dirichlet.c
 * Functions relevant to Beta, Gamma, and Dirichlet densities,
 * and simple and mixture Dirichlet priors.
 * 
 * SRE, Tue Nov  2 13:42:59 2004 [St. Louis]
 * SVN $Id$
 */

#include <easel.h>
#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_FILEPARSER
#include <esl_fileparser.h>
#endif
#include <esl_vectorops.h>
#include <esl_dirichlet.h>


/* Function:  esl_mixdchlet_Create()
 * Incept:    SRE, Fri Apr  8 10:44:34 2005 [St. Louis]
 *
 * Purpose:   Create a new mixture Dirichlet prior with <N> components,
 *            each with <K> parameters.
 *
 * Returns:   initialized <ESL_MIXDCHLET *> on success.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_MIXDCHLET *
esl_mixdchlet_Create(int N, int K)
{
  ESL_MIXDCHLET *pri = NULL;
  int q;

  if ((pri = malloc(sizeof(ESL_MIXDCHLET))) == NULL) goto FAILURE;
  pri->pq = NULL; 
  pri->alpha = NULL;

  if ((pri->pq = malloc(sizeof(double) * N)) == NULL) goto FAILURE;
  if ((pri->alpha = malloc(sizeof(double *) * N)) == NULL) goto FAILURE;
  pri->alpha[0] = NULL;

  if ((pri->alpha[0] = malloc(sizeof(double) * N * K)) == NULL) goto FAILURE;
  for (q = 1; q < N; q++)
    pri->alpha[q] = pri->alpha[0] + q*K;

  pri->N = N;
  pri->K = K;
  return pri;

 FAILURE:
  esl_mixdchlet_Destroy(pri);
  ESL_ERROR_NULL(eslEMEM, "malloc failed in esl_mixdchlet_Create()");
}

/* Function:  esl_mixdchlet_Destroy()
 * Incept:    SRE, Fri Apr  8 11:00:19 2005 [St. Louis]
 *
 * Purpose:   Free's the mixture Dirichlet <pri>.
 */
void
esl_mixdchlet_Destroy(ESL_MIXDCHLET *pri)
{
  if (pri     == NULL)  return;
  if (pri->pq != NULL)  free(pri->pq);
  if (pri->alpha != NULL) {
    if (pri->alpha[0] != NULL) free(pri->alpha[0]); 
    free(pri->alpha);
  }
  free(pri);
}

/* Function:  
 * Incept:    SRE, Fri Apr  8 12:47:03 2005 [St. Louis]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
ESL_MIXDCHLET *
esl_mixdchlet_Read(ESL_FILEPARSER *efp, char *errbuf)
{
  ESL_MIXCHLET *pri;
  int   K;			/* Dirichlet param vector size */
  int   N;			/* number of mixture components */
  char *tok;			/* ptr to a whitespace-delim, noncomment token */
  int   toklen;			/* length of a parsed token */
  int   status;			/* return status of an Easel call */
  int   q;			/* counter over mixture components (0..N-1) */
  int   i;			/* counter over params (0..K-1) */
  
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
  K = atoi(tok);
  if (K < 1) { sprintf(errbuf, "Bad vector size %.32s\n", tok); goto FAILURE; }
  
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
  N = atoi(tok);
  if (N < 1) { sprintf(errbuf, "Bad mixture number %.32s\n", tok); goto FAILURE; }

  pri = esl_mixdchlet_Create(N, K);
  if (pri == NULL) { sprintf(errbuf, "mxdchlet alloc failed\n", tok); goto FAILURE; }
 
  for (q = 0; q < N; q++)
    {
      if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
      pri->pq[q] = atof(tok);
      if (pri->pq[q] < 0.0 || pri->pq[q] > 1.0) 
	{ sprintf(errbuf, "bad mixture coefficient %.32s\n", tok); goto FAILURE; }      

      for (i = 0; i < K; i++)
	{
	  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
	  pri->alpha[q][i] = atof(tok);
	  if (pri->alpha[q][i] <= 0.0)
	    { sprintf(errbuf, "Dirichlet params must be positive, got %.32s\n", tok); goto FAILURE; }      	  
	}
    }
  esl_vec_DNorm(pri->pq, N);
  return pri;
}


/* Function:  esl_dirichlet_LogProbData()
 * Incept:    SRE, Tue Nov  2 14:22:37 2004 [St. Louis]
 *
 * Purpose:   Given an observed count vector $c[0..K-1]$, 
 *            and a simple Dirichlet density parameterized by
 *            $\alpha[0..K-1]$;
 *            calculate and return $\log P(c \mid \alpha)$.
 *            
 *            This is $\int P(c \mid p) P(p \mid \alpha) dp$,
 *            an integral that can be solved analytically.
 *
 * Args:      c          - count vector, [0..K-1]
 *            alpha      - Dirichlet parameters, [0..K-1]
 *            K          - size of c, alpha vectors
 *            ret_answer - RETURN: log P(c | \alpha)
 *
 * Returns:   $\log P(c \mid \alpha)$.
 */
double
esl_dirichlet_LogProbData(double *c, double *alpha, int K)
{
  double lnp;      
  double sum1, sum2, sum3;
  double a1, a2, a3;
  int   x;

  sum1 = sum2 = sum3 = lnp = 0.0;
  for (x = 0; x < K; x++)
    {
      sum1 += c[x] + alpha[x];
      sum2 += alpha[x];
      sum3 += c[x];
      esl_dirichlet_LogGamma(alpha[x] + c[x], &a1); 
      esl_dirichlet_LogGamma(c[x] + 1.,       &a2);
      esl_dirichlet_LogGamma(alpha[x],        &a3);
      lnp  += a1 - a2 - a3;
    }
  esl_dirichlet_LogGamma(sum1,      &a1);
  esl_dirichlet_LogGamma(sum2,      &a2);
  esl_dirichlet_LogGamma(sum3 + 1., &a3);
  lnp += a2 + a3 - a1;

  return lnp;
}

int
esl_mixdchlet_MPParameters(double *c, int K, MIXDCHLET *pri, double *mix, double *p)
{
  int q;			/* counter over mixture components */
  double val;
  double totc;
  double tota;
  
  if (K != pri-K) ESL_ERROR(eslEINCOMPAT, "cvec's K != mixture Dirichlet's K");

  /* Calculate mix[], the posterior probability
   * P(q | c) of mixture component q given the count vector c.
   */
  for (q = 0; q < pri->nq; q++)
    if (pri->pq[q] > 0.0)  
      {
	mix[q] =  log(pri->pq[q]);
	mix[q] += esl_dirichlet_LogProbData(c, pri->alpha[q], K);
      }
    else
      mix[q] = -HUGE_VAL;
  esl_vec_DLogNorm(mix, pri->nq);
  esl_vec_DExp(mix, pri->nq);	/* mix[q] is now P(q|c) */

  totc = esl_vec_DSum(c, K);
  esl_vec_DSet(p, K, 0.);
  for (x = 0; x < K; x++)
    for (q = 0; q < pri->nq; q++)
      {
	tota = esl_vec_DSum(pri->alpha[q], K);
	p[x] += mix[q] * (c[x] + pri->alpha[q][x]) / (totc + tota);
      }
  /* should be normalized already, but for good measure: */
  esl_vec_DNorm(p, K);
}


/* Function:  esl_dirichlet_LogGamma()
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
 * Returns:   Put the answer in <ret_answer>; returns <eslOK>.
 *            
 * Throws:    <eslEINVAL> if $x <= 0$.
 */
int
esl_dirichlet_LogGamma(double x, double *ret_answer)
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
  if (x <= 0.0)  ESL_ERROR(eslEINVAL, "invalid x <= 0 in esl_dirichlet_LogGamma()");

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





/*****************************************************************
 * Sampling code: 
 * only included when augmented by esl_random module.
 ***************************************************************** 
 */
#ifdef eslAUGMENT_RANDOM

static double gamma_ahrens(ESL_RANDOMNESS *r, double a);
static double gamma_integer(ESL_RANDOMNESS *r, unsigned int a);
static double gamma_fraction(ESL_RANDOMNESS *r, double a);

/* Function:  esl_dirichlet_Sample()
 * Incept:    SRE, Tue Nov  2 14:30:31 2004 [St. Louis]
 *
 * Purpose:   Given a Dirichlet density parameterized by $\alpha[0..K-1]$,
 *            sample a probability vector $p[0..K-1]$ from
 *            $P(p \mid \alpha)$.
 *
 * Args:      r      - random number generation object
 *            alpha  - parameters of Dirichlet density [0..K-1]
 *            K      - vector size
 *            p      - RETURN: sampled probability vector
 *                     (caller allocates 0..K-1).         
 *
 * Returns:   <eslOK>, and <p> will contain the sampled vector.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_dirichlet_Sample(ESL_RANDOMNESS *r, double *alpha, int K, double *p)
{
  int    x;
  for (x = 0; x < K; x++) esl_dirichlet_SampleGamma(r, alpha[x], &(p[x]));
  esl_vec_DNorm(p, K);
}


/* Function: esl_dirichlet_SampleGamma()
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
 *           returns <eslOK>.
 *
 * Throws:   <eslEINVAL> for $a <= 0$.
 */
int
esl_dirichlet_SampleGamma(ESL_RANDOMNESS *r, double a, double *ret_answer)
{
  double aint;

  if (a <= 0.) ESL_ERROR(eslEINVAL, "a <= 0 in esl_dirichlet_SampleGamma()");

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

#endif /*eslAUGMENT_RANDOM*/
