/* esl_dirichlet.c
 * Functions relevant to Beta, Gamma, and Dirichlet densities,
 * and simple and mixture Dirichlet priors.
 * 
 * SRE, Tue Nov  2 13:42:59 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <easel.h>
#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_FILEPARSER
#include <esl_fileparser.h>
#endif
#include <esl_vectorops.h>
#include <esl_stats.h>
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


/* Function:  esl_mixdchlet_MPParameters()
 * Incept:    SRE, Sat Apr  9 14:28:26 2005 [St. Louis]
 *
 * Purpose:   Parameter estimation for a count vector <c> of cardinality
 *            <K>, and a mixture Dirichlet prior <pri>. Calculates
 *            mean posterior estimates for probability parameters, and
 *            returns them in <p>. Also returns the posterior probabilities
 *            of each Dirichlet mixture component, $P(q \mid c)$, in <mix>.
 *            Caller must provide allocated space for <mix> and <p>, both
 *            of length <K>.
 *
 * Returns:   <eslOK> on success; <mix> contains posterior probabilities of
 *            the Dirichlet components, and <p> contains mean posterior 
 *            probability parameter estimates.
 *
 * Throws:    <esl_EINCOMPAT> if <pri> has different cardinality than <c>.
 */
int
esl_mixdchlet_MPParameters(double *c, int K, ESL_MIXDCHLET *pri, double *mix, double *p)
{
  int q;			/* counter over mixture components */
  int x;
  double val;
  double totc;
  double tota;
  
  if (K != pri->K) ESL_ERROR(eslEINCOMPAT, "cvec's K != mixture Dirichlet's K");

  /* Calculate mix[], the posterior probability
   * P(q | c) of mixture component q given the count vector c.
   */
  for (q = 0; q < pri->N; q++)
    if (pri->pq[q] > 0.0)  
      {
	esl_dirichlet_LogProbData(c, pri->alpha[q], K, &val);
	mix[q] =  val + log(pri->pq[q]);
      }
    else
      mix[q] = -HUGE_VAL;
  esl_vec_DLogNorm(mix, pri->N); /* mix[q] is now P(q|c) */

  totc = esl_vec_DSum(c, K);
  esl_vec_DSet(p, K, 0.);
  for (x = 0; x < K; x++)
    for (q = 0; q < pri->N; q++)
      {
	tota = esl_vec_DSum(pri->alpha[q], K);
	p[x] += mix[q] * (c[x] + pri->alpha[q][x]) / (totc + tota);
      }
  /* should be normalized already, but for good measure: */
  esl_vec_DNorm(p, K);
  return eslOK;
}


/* Function:  esl_dirichlet_LogProbData()
 * Incept:    SRE, Tue Nov  2 14:22:37 2004 [St. Louis]
 *
 * Purpose:   Given an observed count vector $c[0..K-1]$, 
 *            and a simple Dirichlet density parameterized by
 *            $\alpha[0..K-1]$;
 *            calculate $\log P(c \mid \alpha)$.
 *            
 *            This is $\int P(c \mid p) P(p \mid \alpha) dp$,
 *            an integral that can be solved analytically.
 *
 * Args:      c          - count vector, [0..K-1]
 *            alpha      - Dirichlet parameters, [0..K-1]
 *            K          - size of c, alpha vectors
 *            ret_answer - RETURN: log P(c | \alpha)
 *
 * Returns:   <eslOK> on success, and puts result $\log P(c \mid \alpha)$
 *            in <ret_answer>.
 */
int
esl_dirichlet_LogProbData(double *c, double *alpha, int K, double *ret_answer)
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
      esl_stats_LogGamma(alpha[x] + c[x], &a1); 
      esl_stats_LogGamma(c[x] + 1.,       &a2);
      esl_stats_LogGamma(alpha[x],        &a3);
      lnp  += a1 - a2 - a3;
    }
  esl_stats_LogGamma(sum1,      &a1);
  esl_stats_LogGamma(sum2,      &a2);
  esl_stats_LogGamma(sum3 + 1., &a3);
  lnp += a2 + a3 - a1;

  *ret_answer = lnp;
  return eslOK;
}


/* Function:  esl_dirichlet_LogProbProbs()
 * Incept:    SRE, Sat Apr  9 14:35:17 2005 [St. Louis]
 *
 * Purpose:   Given Dirichlet parameter vector <alpha> and a probability
 *            vector <p>, both of cardinality <K>; return
 *            $\log P(p \mid alpha)$.
 *            
 * Returns:   <eslOK> on success, and the result is in <ret_answer>.           
 *            
 * Xref:      Sjolander (1996) appendix, lemma 2.
 */
int
esl_dirichlet_LogProbProbs(double *p, double *alpha, int K, double *ret_answer)
{
  double sum;		        /* for Gammln(|alpha|) in Z     */
  double logp;			/* RETURN: log P(p|alpha)       */
  double val;
  int x;

  sum = logp = 0.0;
  for (x = 0; x < K; x++)
    if (p[x] > 0.0)		/* any param that is == 0.0 doesn't exist */
      {
	esl_stats_LogGamma(alpha[x], &val);
	logp -= val;
	logp += (alpha[x]-1.0) * log(p[x]);
	sum  += alpha[x];
      }
  esl_stats_LogGamma(sum, &val);
  logp += val;
  *ret_answer = logp;
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
 */
int
esl_dirichlet_Sample(ESL_RANDOMNESS *r, double *alpha, int K, double *p)
{
  int x;
  int status;  

  for (x = 0; x < K; x++) 
    if ((status = esl_dirichlet_SampleGamma(r, alpha[x], &(p[x]))) != eslOK) return status;
  esl_vec_DNorm(p, K);
  return eslOK;
}

/* Function:  esl_dirichlet_SampleBeta()
 * Incept:    SRE, Sat Oct 25 12:20:31 2003 [Stanford]
 *
 * Purpose:   Samples from a Beta(theta1, theta2) density, leaves answer
 *            in <ret_answer>. (Special case of sampling Dirichlet.)
 *            
 * Returns:   <eslOK>.           
 */
int
esl_dirichlet_SampleBeta(ESL_RANDOMNESS *r, double theta1, double theta2, double *ret_answer)
{
  int status;
  double p, q;

  if ((status = esl_dirichlet_SampleGamma(r, theta1, &p)) != eslOK) return status;
  if ((status = esl_dirichlet_SampleGamma(r, theta2, &q)) != eslOK) return status;
  *ret_answer = p / (p+q);
  return eslOK;
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
 *           <gamma_ahrens()>, <gamma_integer()>, and <gamma_fraction()>.
 *
 * Args:     r          - random number generation seed
 *           alpha      - order of the gamma function
 *           ret_answer - RETURN: a sample from Gamma(a, 1).
 *
 * Returns:  a gamma-distributed deviate is put in <ret_answer>;
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
      Y = tan(eslCONST_PI * esl_random(r)); 
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
  
  p = eslCONST_E / (a + eslCONST_E);
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



/*****************************************************************
 * File input code:
 * only included when augmented by esl_fileparser module.
 ***************************************************************** 
 */
#ifdef eslAUGMENT_FILEPARSER 
/* Function:  esl_mixdchlet_Read()
 * Incept:    SRE, Fri Apr  8 12:47:03 2005 [St. Louis]
 *
 * Purpose:   Reads a mixture Dirichlet from an open stream <efp>, using the 
 *            <ESL_FILEPARSER> token-based parser. 
 *            
 *            The first two tokens are <K>, the length of the Dirichlet parameter
 *            vector(s), and <N>, the number of mixture components. Then for
 *            each of the <N> mixture components <i>, it reads a mixture coefficient
 *            <pq[i]> followed by <K> Dirichlet parameters <alpha[i][0..K-1]>.
 *            
 *            This function may be called more than once on the same open file,
 *            to read multiple different mixture Dirichlets from it (transitions,
 *            match emissions, insert emissions, for example).
 *
 * Returns:   <eslOK> on success, and <ret_pri> contains a new <ESL_MIXDCHLET> object 
 *            that the caller is responsible for free'ing.
 *
 * Throws:    <eslEFORMAT> on parse failure, in which case <efp->errbuf>
 *            contains an informative diagnostic message, and <efp->linenumber>
 *            contains the linenumber at which the parse failed.
 */
int
esl_mixdchlet_Read(ESL_FILEPARSER *efp,  ESL_MIXDCHLET **ret_pri)
{
  ESL_MIXDCHLET *pri;
  int   K;			/* Dirichlet param vector size */
  int   N;			/* number of mixture components */
  char *tok;			/* ptr to a whitespace-delim, noncomment token */
  int   toklen;			/* length of a parsed token */
  int   status;			/* return status of an Easel call */
  int   q;			/* counter over mixture components (0..N-1) */
  int   i;			/* counter over params (0..K-1) */
  
  *ret_pri = pri = NULL;

  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
  K = atoi(tok);
  if (K < 1) { sprintf(efp->errbuf, "Bad vector size %.32s", tok); goto FAILURE; }
  
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
  N = atoi(tok);
  if (N < 1) { sprintf(efp->errbuf, "Bad mixture number %.32s", tok); goto FAILURE; }

  pri = esl_mixdchlet_Create(N, K);
  if (pri == NULL) { sprintf(efp->errbuf, "mxdchlet alloc failed"); goto FAILURE; }
 
  for (q = 0; q < N; q++)
    {
      if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
      pri->pq[q] = atof(tok);
      if (pri->pq[q] < 0.0 || pri->pq[q] > 1.0) 
	{ sprintf(efp->errbuf, "bad mixture coefficient %.32s", tok); goto FAILURE; }      

      for (i = 0; i < K; i++)
	{
	  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto FAILURE;
	  pri->alpha[q][i] = atof(tok);
	  if (pri->alpha[q][i] <= 0.0)
	    { sprintf(efp->errbuf, "Dirichlet params must be positive, got %.32s", tok); goto FAILURE; } 
	}
    }
  esl_vec_DNorm(pri->pq, N);
  *ret_pri = pri;
  return eslOK;

 FAILURE:
  esl_mixdchlet_Destroy(pri);
  return eslEFORMAT;
}
#endif /* eslAUGMENT_FILEPARSER */




/*****************************************************************
 * Example main():
 *****************************************************************/
#ifdef eslDIRICHLET_EXAMPLE
/*::cexcerpt::dirichlet_example::begin::*/
/* compile: 
    gcc -g -Wall -I. -o example -DeslDIRICHLET_EXAMPLE\
      -DeslAUGMENT_RANDOM -DeslAUGMENT_FILEPARSER esl_random.c esl_fileparser.c\
      esl_vectorops.c esl_dirichlet.c easel.c -lm
 * run:     ./example <mixture Dirichlet file>
 */
#include <stdlib.h>
#include <stdio.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_fileparser.h>
#include <esl_vectorops.h>
#include <esl_dirichlet.h>

int
main(int argc, char **argv)
{
  FILE           *fp;
  ESL_FILEPARSER *efp;
  ESL_RANDOMNESS *r;
  ESL_MIXDCHLET  *pri;
  int             c,i,q,qused;
  double         *counts, *probs, *iq, *ip;

  /* Read in a mixture Dirichlet from a file. */
  fp  = fopen(argv[1], "r");
  efp = esl_fileparser_Create(fp);
  if (esl_mixdchlet_Read(efp, &pri) != eslOK) {
    fprintf(stderr, "%s;\ndirichlet file %s parse failed at line %d\n",
	    efp->errbuf, argv[1], efp->linenumber);
    exit(1);
  }
  esl_fileparser_Destroy(efp);
  fclose(fp);  

  /* Allocate some working spaces */
  probs  = malloc(sizeof(double) * pri->K);
  counts = malloc(sizeof(double) * pri->K);
  iq     = malloc(sizeof(double) * pri->N);
  ip     = malloc(sizeof(double) * pri->K);

  /* Sample a probability vector from it. */
  r = esl_randomness_CreateTimeseeded(); /* init the random generator */
  qused = esl_rnd_DChoose(r, pri->pq, pri->N); /* sample a component */
  esl_dirichlet_Sample(r, pri->alpha[qused], pri->K, probs);

  printf("Component %2d: p[] = ", qused);
  for (i = 0; i < pri->K; i++) printf("%.3f ", probs[i]);
  printf("\n");

  /* Sample a count vector from that prob vector. */
  esl_vec_DSet(counts, pri->K, 0.);
  for (c = 0; c < 20; c++)
    counts[esl_rnd_DChoose(r, probs, pri->K)] += 1.;

  printf("              c[] = ");
  for (i = 0; i < pri->K; i++) printf("%5.0f ", counts[i]);
  printf("\n");

  /* Estimate a probability vector (ip) from those counts, and
   * also get back the posterior prob P(q|c) of each component (iq). */
  esl_mixdchlet_MPParameters(counts, pri->K, pri, iq, ip);

  printf("  reestimated p[] = ");
  for (i = 0; i < pri->K; i++) printf("%.3f ", ip[i]);
  printf("\n");

  q = esl_vec_DArgMax(iq, pri->N);
  printf("probably generated by component %d; P(q%d | c) = %.3f\n",
	 q, q, iq[q]);

  esl_mixdchlet_Destroy(pri);
  free(probs); free(counts); free(iq); free(ip);
  return 0;
}
/*::cexcerpt::dirichlet_example::end::*/
#endif /*eslDIRICHLET_EXAMPLE*/

/*****************************************************************
 * Test driver:
 * gcc -g -Wall -I. -o test -DeslDIRICHLET_TESTDRIVE -DeslAUGMENT_FILEPARSER\
 *    -DeslAUGMENT_RANDOM esl_fileparser.c esl_random.c esl_vectorops.c\
 *    esl_dirichlet.c easel.c -lm
 * ./test
 *****************************************************************/
#ifdef eslDIRICHLET_TESTDRIVE
#define NCOMPONENTS 2
#define NALPHA      6		/* dice example, 6 faces */
#define NCOUNTS     1000
#define NTRIALS     100

int
main(void)
{
  ESL_FILEPARSER *efp;
  ESL_RANDOMNESS *r;
  ESL_MIXDCHLET  *d1, *d2;
  FILE *fp;
  char  filename[] = "tmpxxx.pri";
  int   q, i, c, t;

  double pq[NCOMPONENTS] = {0.5, 0.5};
  double alpha[NCOMPONENTS][NALPHA] = { {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
					{0.1, 0.1, 0.1, 0.1, 0.1, 0.1} };
  double counts[NALPHA];
  double probs[NALPHA];
  double iq[NCOMPONENTS];	/* inferred posterior probs over components */
  double ip[NALPHA];		/* inferred probability parameters */
  int    qused;
  int    qguess;		/* inferred guess at which component  */
  double maxdeviation;

  /* Get hold of some reproducible randomness.
   * (It has to be reproducible, because we need to be able
   *  to guarantee that the tests will succeed, even though
   *  we're doing a stochastic sampling procedure.)
   */
  if ((r = esl_randomness_Create(42)) == NULL) abort();

  /* Create a mixture Dirichlet file.
   */
  if ((fp = fopen(filename, "w")) == NULL) abort();
  fprintf(fp, "%d %d\n", NALPHA, NCOMPONENTS);
  for (q = 0; q < NCOMPONENTS; q++)
    {
      fprintf(fp, "%.3f ", pq[q]);
      for (i = 0; i < NALPHA; i++)
	fprintf(fp, "%.3f ", alpha[q][i]);
      fprintf(fp, "\n");
    }
  fclose(fp);
  
  /* Read it back in.
   */
  if ((fp = fopen(filename, "r")) == NULL) abort();
  if ((efp = esl_fileparser_Create(fp)) == NULL) abort();
  if (esl_mixdchlet_Read(efp, &d1) != eslOK) abort();
  esl_fileparser_Destroy(efp);
  fclose(fp);

  /* Make a copy of it - artificially testing the _Create() call.
   */
  if ((d2 = esl_mixdchlet_Create(d1->N, d1->K)) == NULL) abort();
  esl_vec_DCopy(d2->pq, d1->pq, d1->N);
  for (q = 0; q < d1->N; q++)
    esl_vec_DCopy(d2->alpha[q], d1->alpha[q], d1->K);

  /* Sample from it.
   */
  for (t = 0; t < NTRIALS; t++)
    {
      qused = esl_rnd_DChoose(r, d2->pq, d2->N); /* sample a component */
      esl_dirichlet_Sample(r, d2->alpha[qused], d2->K, probs);
      esl_vec_DSet(counts, NALPHA, 0.);
      for (c = 0; c < NCOUNTS; c++)
	{
	  i = esl_rnd_DChoose(r, probs, NALPHA);
	  counts[i] += 1.;
	}
      /* printf("%1d ", qused); */
  
      /* Classify by posterior inference on the sampled probability vector.
       */
      for (q = 0; q < d2->N; q++)
	{
	  esl_dirichlet_LogProbProbs(probs, d2->alpha[q], d2->K, &(iq[q]));
	  iq[q] += log(d2->pq[q]);
	}
      esl_vec_DLogNorm(iq, d2->N);
      qguess = esl_vec_DArgMax(iq, d2->N); /* the MP guess from the probs */
      /* printf("%1d ", qguess); */
      if (qused != qguess) abort();
  
      /* Classify by posterior inference on the sampled count vector;
       * and attempt to estimate the probability vector.
       */
      esl_mixdchlet_MPParameters(counts, d2->K, d2, iq, ip);
      qguess = esl_vec_DArgMax(iq, d2->N); /* the MP guess from the counts */
      /* printf("%1d ", qguess); */
      if (qused != qguess) abort();

      for (i = 0; i < d2->K; i++)
	ip[i] = fabs(ip[i] - probs[i]); /* ip[] is now the differences rel to probs */
      maxdeviation = esl_vec_DMax(ip, d2->K);
      /* printf("%.3f\n", maxdeviation); */
      if (maxdeviation > 0.05) abort();

    }
  esl_mixdchlet_Destroy(d1);
  esl_mixdchlet_Destroy(d2);
  return 0;
}
#endif /*eslDIRICHLET_TESTDRIVE*/
