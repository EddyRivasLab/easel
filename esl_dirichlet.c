/* Functions relevant to Beta, Gamma, and Dirichlet densities,
 * including simple and mixture Dirichlet priors.
 * 
 * Contents:
 *   1. The <ESL_MIXDCHLET> object for mixture Dirichlet priors
 *   2. Dirichlet likelihood functions
 *   3. Sampling from Dirichlets              [with <random>]
 *   4. Reading mixture Dirichlets from files [with <fileparser>]
 *   5. Unit tests
 *   6. Test driver
 *   7. Example
 *   8. Copyright and license information
 * 
 * SRE, Tue Nov  2 13:42:59 2004 [St. Louis]
 * SVN $Id$
 */
#include <esl_config.h>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "easel.h"
#ifdef eslAUGMENT_RANDOM
#include "esl_random.h"
#endif
#ifdef eslAUGMENT_MINIMIZER
#include "esl_minimizer.h"
#endif
#ifdef eslAUGMENT_FILEPARSER
#include "esl_fileparser.h"
#endif
#include "esl_vectorops.h"
#include "esl_stats.h"
#include "esl_dirichlet.h"


/*****************************************************************
 *# 1. The <ESL_MIXDCHLET> object for mixture Dirichlet priors
 *****************************************************************/

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
  int status;
  ESL_MIXDCHLET *pri = NULL;
  int q;

  ESL_ALLOC(pri, sizeof(ESL_MIXDCHLET));
  pri->pq = NULL; 
  pri->alpha = NULL;

  ESL_ALLOC(pri->pq,    sizeof(double)   * N);
  ESL_ALLOC(pri->alpha, sizeof(double *) * N);
  pri->alpha[0] = NULL;

  ESL_ALLOC(pri->alpha[0], sizeof(double) * N * K);
  for (q = 1; q < N; q++)
    pri->alpha[q] = pri->alpha[0] + q*K;

  pri->N = N;
  pri->K = K;
  return pri;

 ERROR:
  esl_mixdchlet_Destroy(pri);
  return NULL;
}

/* Function:  esl_mixdchlet_Compare()
 * Synopsis:  Compare two mixture Dirichlets for equality.
 * Incept:    SRE, Sat May 30 09:37:40 2009 [Stockholm]
 *
 * Purpose:   Compare mixture Dirichlet objects <d1> and <d2> 
 *            for equality. For real numbered values, equality
 *            is defined by <esl_DCompare()> with a fractional
 *            tolerance <tol>.                    
 *
 * Returns:   <eslOK> on equality; <eslFAIL> otherwise.
 */
int
esl_mixdchlet_Compare(ESL_MIXDCHLET *d1, ESL_MIXDCHLET *d2, double tol)
{
  int q;

  if (d1->N != d2->N) return eslFAIL;
  if (d1->K != d2->K) return eslFAIL;

  if (esl_vec_DCompare(d1->pq, d2->pq, d1->N, tol) != eslOK) return eslFAIL;
  
  for (q = 0; q < d1->N; q++)
    if (esl_vec_DCompare(d1->alpha[q], d2->alpha[q], d1->K, tol) != eslOK) return eslFAIL;    

  return eslOK;
}

/* Function:  esl_mixdchlet_Copy()
 * Synopsis:  Copies mixture dirichlet object <d> to <d_dst>.
 *            Both objects are of size <N> and <K>.  
 *            <d> is unchanged. 
 * Incept:    ER, Thu Jun 18 10:30:06 2009 [Janelia]
 *
 * Purpose:                      
 *
 * Returns:   <eslOK> on equality; <eslFAIL> otherwise.
 */
int
esl_mixdchlet_Copy(ESL_MIXDCHLET *d, ESL_MIXDCHLET *d_dst)
{
  int q;

  if (d->N != d_dst->N) return eslFAIL;
  if (d->K != d_dst->K) return eslFAIL;

  esl_vec_DCopy(d->pq, d->N, d_dst->pq);
  
  for (q = 0; q < d->N; q++)
    esl_vec_DCopy(d->alpha[q], d->K, d_dst->alpha[q]);

  return eslOK;
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
  
  if (K != pri->K) ESL_EXCEPTION(eslEINCOMPAT, "cvec's K != mixture Dirichlet's K");

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
/*---------------- end, ESL_MIXDCHLET ---------------------------*/


/*****************************************************************
 *# 2. Dirichlet likelihood functions
 *****************************************************************/

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
  int    x;

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

/* Function:  esl_dirichlet_LogProbData_Mixture()
 * Incept:    ER, Wed Jun 17 14:41:23 EDT 2009 [janelia]
 *
 * Purpose:   Given an observed count vector $c[0..K-1]$, 
 *            and a mixture Dirichlet density parameterized by
 *            $\alpha_1[0..K-1]$ ... $\alpha_N[0..K-1]$;
 *            calculate $\log \sum_i pq_i * P(c \mid \alpha_i)$.
 *            
 *
 * Args:      c          - count vector, [0..K-1]
 *            d      - Dirichlet parameters, [0..K-1]
 *            ret_answer - RETURN: log P(c | \alpha)
 *
 * Returns:   <eslOK> on success, and puts result $\log P(c \mid \alpha)$
 *            in <ret_answer>.
 */
int
esl_dirichlet_LogProbData_Mixture(double *c, ESL_MIXDCHLET *d, double *ret_answer)
{
  double *mixq = NULL;
  double  lnp;
  double  val;
  int     q;             /* counter over mixture components */
  int     status;

  ESL_ALLOC(mixq, sizeof(double)*d->N);

  for (q = 0; q < d->N; q++) {
    esl_dirichlet_LogProbData(c, d->alpha[q], d->K, &val);
    mixq[q] = val + log(d->pq[q]);
  }
  lnp = esl_vec_DLogSum(mixq, d->N);

  free(mixq);

  *ret_answer = lnp;
  return eslOK;

 ERROR:
  free(mixq);
  return status;
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
/*----------- end, Dirichlet likelihood functions ---------------*/

/*****************************************************************
 * Dirichlet Maximum likelihood fit from counts
 *****************************************************************/

#ifdef eslAUGMENT_MINIMIZER
/* This structure is used to sneak the data into minimizer's generic
 * (void *) API for all aux data
 */
struct mixdchlet_data {
  ESL_MIXDCHLET  *d;      /* the dirichlet mixture parameters */
  double        **c;      /* count vector array [0..nc-1][0..alphabet_size(d->K)] */
  int             nc;     /* number of count samples */
};

/*****************************************************************
 * Parameter vector packing/unpacking
 *
 * The conjugate gradient code takes a single parameter vector <p>,
 * where the values are unconstrained real numbers.
 *
 * We have a mixture Dirichlet with two kinds of parameters.
 * pq_i are mixture coefficients, constrained to be >= 0 and
 * \sum_i pq_i = 1.  alpha^i_x are the Dirichlet parameters
 * for component i, constrained to be > 0.
 *
 * Our p's are therefore not only packed into a single vector;
 * they're reparameterized to implement the constraints:
 *   for a Dirichlet parameter:
 *      alpha = exp(p)   p = log(alpha)
 *      (thus, alpha > 0 for all real p)
 *
 *   for a mixture coefficient:
 *      pq = exp(p) / \sum_a exp(p_a)
 *      (thus, 0 < pq < 1, \sum_a pq_a = 1, for all real p)
 *
 * Conjugate gradients optimizes the <p> parameter vector,
 * but we can convert that back out into a Dirichlet answer.
 *
 * The packing order is: the first N terms of a parameter vector are
 * the mixture coefficients pq_i. N different alpha_i vectors follow.
 *
 * [0 ... N-1] [0 ... K-1] [0 ... K-1]  ...
 *     q's      alpha_0     alpha_1     ...
 *
 * In both functions below, p, pq, and alpha are all allocated
 * and free'd by the caller.
 *      p : length N + N*K = N*(K+1)  [0.. N*(K+1)-1]
 *     pq : length N, [0..N-1]
 *  alpha : NxK, [0..N-1][0..K-1].
 */
static void
mixdchlet_pack_paramvector(double *p, int np, ESL_MIXDCHLET *d)
{
  int q;			/* counter over mixture components */
  int x;                        /* counter in alphabet size */

  /* the mixture coeficients */
  for (q = 0; q < d->N; q++)
    p[q] = log(d->pq[q]);

  /* the dirichlet parameters */
  for (q = 0; q < d->N; q++)
    for (x = 0; x < d->K; x++)
      p[d->N + q*d->K + x] = log(d->alpha[q][x]);
 
}

/* Same as above but in reverse: given parameter vector <p>,
 * do appropriate c.o.v. back to desired parameter space, and
 * update the mixdchlet <d>.
 */
static void
mixdchlet_unpack_paramvector(double *p, int np, ESL_MIXDCHLET *d)
{
  int q;			/* counter over mixture components */
  int x;                        /* counter in alphabet size */

  /* the mixture coeficients */
  for (q = 0; q < d->N; q++)
    d->pq[q] = exp(p[q]);
  esl_vec_DNorm(d->pq, d->N);

  /* the dirichlet parameters */
  for (q = 0; q < d->N; q++)
    for (x = 0; x < d->K; x++)
      d->alpha[q][x] = exp(p[d->N + q*d->K + x]);
}

/* The log likelihood function to be optimized by ML fitting:
 *   This needs to be careful of a case where a lambda = inf.
 */
static double
mixdchlet_complete_func(double *p, int np, void *dptr)
{
  struct mixdchlet_data *data = (struct mixdchlet_data *) dptr;
  ESL_MIXDCHLET         *d    = data->d;
  double  logPsample;
  double  logP = 0.;
  int     m;             /* counter over count samples */
 
  mixdchlet_unpack_paramvector(p, np, d);

  for (m = 0; m < data->nc; m++) {
    esl_dirichlet_LogProbData_Mixture(data->c[m], d, &logPsample);
    logP += logPsample;
 }

  if (isnan(logP)) esl_fatal("logP is NaN");

  return -logP;
}

/* The gradient of the NLL w.r.t. each free parameter in p.
 */
static void
mixdchlet_complete_gradient(double *p, int np, void *dptr, double *dp)
{
  struct mixdchlet_data *data = (struct mixdchlet_data *) dptr;
  ESL_MIXDCHLET         *d    = data->d;
  double  sum_alpha;
  double  sum_c;
  double  val1;
  double  val2;
  double  ratio1;
  double  ratio2;
  double  psi1, psi2, psi3, psi4;
  int     m;                     /* counter over count samples */
  int     q;		 	 /* counter over mixture components */
  int     x;                     /* counter in alphabet size */
 
  mixdchlet_unpack_paramvector(p, np, d);

  /* initialize */
  esl_vec_DSet(dp, np, 0.0);
  
  for (q = 0; q < d->N; q++) {
    
    sum_alpha = esl_vec_DSum(d->alpha[q], d->K);
    esl_stats_Psi(sum_alpha, &psi1);
   
    for (m = 0; m < data->nc; m++) {
      sum_c = esl_vec_DSum(data->c[m], d->K);
      esl_stats_Psi(sum_alpha+sum_c, &psi2);
     
      esl_dirichlet_LogProbData(data->c[m], d->alpha[q], d->K, &val1);
      esl_dirichlet_LogProbData_Mixture(data->c[m], d, &val2);
      
      ratio1 = exp(val1 - val2);
      ratio2 = ratio1 * d->pq[q];

      /* derivative respect to the mixture coeficients */
      dp[q] += ratio1;
      
      /* derivative respect to the dirichlet parameters */
      for (x = 0; x < d->K; x++) {
	esl_stats_Psi(d->alpha[q][x]+data->c[m][x], &psi3);
	esl_stats_Psi(d->alpha[q][x],               &psi4);
	
	dp[d->N + q*d->K + x] += ratio2 * (psi1 - psi2 + psi3 - psi4);
     }              
    }
  }
  
  /* check */
  for (q = 0; q < d->N; q++) {
    if (isnan(dp[q])) esl_fatal("dp for pq[%d] is NaN", q);
    for (x = 0; x < d->K; x++) 
      if(isnan(dp[d->N + q*d->K + x])) esl_fatal("dp for alpha[%d][%d] is NaN", q, x);
  }
  
}

/* Function:  esl_mixdchlet_Fit()
 * Incept:    ER, Wed Jun 17 10:58:50 2009 [Janelia]
 *
 * Purpose:   Given a count vector <c>, and an initial guess <d> for
 *            a mixdchlet, find maximum likelihood parameters
 *            by conjugate gradient descent optimization, starting
 *            from <d> and leaving the final optimized solution in
 *            <d>.
 *            
 * Returns:   <eslOK> on success, and <d> contains the fitted 
 *            mixdchlet parameters.
 *            
 * Throws:    <eslEMEM> on allocation error, and <d> is left in
 *            in its initial state.           
 */
int
esl_mixdchlet_Fit(double **c, int nc, ESL_MIXDCHLET *d, int be_verbose)
{
  struct mixdchlet_data data;
  double *p   = NULL;
  double *u   = NULL;
  double *wrk = NULL;
  double  tol;
  int     np;
  double  fx;
  int     i;
  int     status;

  tol = 1e-6;

  /* Allocate parameters
   */
  np = d->N *(d->K+1);
  ESL_ALLOC(p,   sizeof(double) * np);
  ESL_ALLOC(u,   sizeof(double) * np);
  ESL_ALLOC(wrk, sizeof(double) * np * 4);

  /* Copy shared info into the "data" structure
   */
  data.d  = d;
  data.c  = c;
  data.nc = nc;

  /* From d, create the parameter vector.
   */
  mixdchlet_pack_paramvector(p, np, d);

  /* Define the step size vector u.
   */
  for (i = 0; i < np; i++) u[i] = 1.0;

  /* Feed it all to the mighty optimizer.
   */
  status = esl_min_ConjugateGradientDescent(p, u, np, 
					    &mixdchlet_complete_func, 
					    &mixdchlet_complete_gradient,
					    (void *) (&data), tol, wrk, &fx);
  if (status != eslOK) goto ERROR;

  /* Convert the final parameter vector back to a mixdchlet
   */
  mixdchlet_unpack_paramvector(p, np, d);

  free(p);
  free(u);
  free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}
#endif /*eslAUGMENT_MINIMIZER*/
/*----------- end, Dirichlet Maximum likelihood fit from counts ---------------*/


/*****************************************************************
 *# 3. Sampling from Dirichlets: requires <esl_random>
 *****************************************************************/
#ifdef eslAUGMENT_RANDOM
/* Function:  esl_dirichlet_DSample()
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
esl_dirichlet_DSample(ESL_RANDOMNESS *r, double *alpha, int K, double *p)
{
  int x;

  for (x = 0; x < K; x++) 
    p[x] = esl_rnd_Gamma(r, alpha[x]);
  esl_vec_DNorm(p, K);
  return eslOK;
}

/* Function:  esl_dirichlet_FSample()
 * Incept:    SRE, Sat Jan  6 17:09:05 2007 [Casa de Gatos]
 *
 * Purpose:   Same as <esl_dirichlet_DSample()>, except it
 *            works in single-precision floats, not doubles.
 */
int
esl_dirichlet_FSample(ESL_RANDOMNESS *r, float *alpha, int K, float *p)
{
  int x;

  for (x = 0; x < K; x++) 
    p[x] = (float) esl_rnd_Gamma(r, (double) alpha[x]);
  esl_vec_FNorm(p, K);
  return eslOK;
}

/* Function:  esl_dirichlet_DSampleUniform()
 * Incept:    SRE, Thu Aug 11 10:12:49 2005 [St. Louis]
 *
 * Purpose:   Sample a probability vector $p[0..K-1]$ uniformly, by
 *            sampling from a Dirichlet of $\alpha_i = 1.0 \forall i$.
 *
 * Args:      r  - source of random numbers
 *            K  - vector size
 *            p  - RETURN: sampled prob vector, caller alloc'ed 0..K-1
 *
 * Returns:   <eslOK>, and <p> will contain the sampled vector.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_dirichlet_DSampleUniform(ESL_RANDOMNESS *r, int K, double *p)
{
  int x;
  for (x = 0; x < K; x++) 
    p[x] = esl_rnd_Gamma(r, 1.0);
  esl_vec_DNorm(p, K);
  return eslOK;
}

/* Function:  esl_dirichlet_FSampleUniform()
 * Incept:    SRE, Sat Jan  6 17:10:54 2007 [Casa de Gatos]
 *
 * Purpose:   Same as <esl_dirichlet_DSampleUniform()>, except it
 *            works in single-precision floats, not doubles.
 */
int
esl_dirichlet_FSampleUniform(ESL_RANDOMNESS *r, int K, float *p)
{
  int x;
  for (x = 0; x < K; x++) 
    p[x] = (float) esl_rnd_Gamma(r, 1.0);
  esl_vec_FNorm(p, K);
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
  double p, q;

  p = esl_rnd_Gamma(r, theta1);
  q = esl_rnd_Gamma(r, theta2);
  *ret_answer = p / (p+q);
  return eslOK;
}
#endif /*eslAUGMENT_RANDOM*/
/*---------------- end, Dirichlet sampling ----------------------*/


/*****************************************************************
 *# 4. Reading mixture Dirichlets from files [requires esl_fileparser]
 *****************************************************************/
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
 * Note:      One reason this function takes an ESL_FILEPARSER instead of 
 *            a filename or an open FILE pointer is that file format errors
 *            in Easel are non-fatal "normal" errors, and we want to record
 *            an informative error message. The ESL_FILEPARSER has an error
 *            buffer for this purpose. 
 *
 * Returns:   <eslOK> on success, and <ret_pri> contains a new <ESL_MIXDCHLET> object 
 *            that the caller is responsible for free'ing.
 *
 *            <eslEFORMAT> on 'normal' parse failure, in which case <efp->errbuf>
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

  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto ERROR;
  K = atoi(tok);
  if (K < 1) { sprintf(efp->errbuf, "Bad vector size %.32s", tok); goto ERROR; }
  
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto ERROR;
  N = atoi(tok);
  if (N < 1) { sprintf(efp->errbuf, "Bad mixture number %.32s", tok); goto ERROR; }

  pri = esl_mixdchlet_Create(N, K);
  if (pri == NULL) { sprintf(efp->errbuf, "mxdchlet alloc failed"); goto ERROR; }
 
  for (q = 0; q < N; q++)
    {
      if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto ERROR;
      pri->pq[q] = atof(tok);
      if (pri->pq[q] < 0.0 || pri->pq[q] > 1.0) 
	{ sprintf(efp->errbuf, "bad mixture coefficient %.32s", tok); goto ERROR; }      

      for (i = 0; i < K; i++)
	{
	  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) goto ERROR;
	  pri->alpha[q][i] = atof(tok);
	  if (pri->alpha[q][i] <= 0.0)
	    { sprintf(efp->errbuf, "Dirichlet params must be positive, got %.32s", tok); goto ERROR; } 
	}
    }
  esl_vec_DNorm(pri->pq, N);
  *ret_pri = pri;
  return eslOK;

 ERROR:
  esl_mixdchlet_Destroy(pri);
  return eslEFORMAT;
}

/* Function:  esl_mixdchlet_Write()
 * Synopsis:  Write a mixture Dirichlet to an open output stream.
 * Incept:    SRE, Sat May 30 09:30:39 2009 [Stockholm]
 *
 * Purpose:   Write mixture Dirichlet <d> to open output stream <d>.
 *
 * Args:      fp   - open output stream
 *            d    - mixture Dirichlet to write
 * 
 * Returns:   <eslOK> on success.
 */
int
esl_mixdchlet_Write(FILE *fp, ESL_MIXDCHLET *d)
{
  int q,i;

  fprintf(fp, "%d %d\n", d->K, d->N);
  for (q = 0; q < d->N; q++)
    {
      fprintf(fp, "%.3f ", d->pq[q]);
      for (i = 0; i < d->K; i++)
	fprintf(fp, "%.3f ", d->alpha[q][i]);
      fprintf(fp, "\n");
    }
  return eslOK;
}


#endif /* eslAUGMENT_FILEPARSER */
/*-------------- end, reading mixture Dirichlets ----------------*/



/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef eslDIRICHLET_TESTDRIVE

static void
utest_io(ESL_MIXDCHLET *d, double tol)
{
  char           *msg       = "esl_dirichlet: io unit test failed";
  ESL_MIXDCHLET  *d2        = NULL;
  ESL_FILEPARSER *efp       = NULL;
  FILE           *fp        = NULL;
  char            tmpfile[] = "esltmpXXXXXX";

  /* Create a mixture Dirichlet file, as a named tmpfile.  */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal(msg);
  if (esl_mixdchlet_Write(fp, d)      != eslOK) esl_fatal(msg);
  fclose(fp);

  /* Read it back in */
  if ((fp = fopen(tmpfile, "r")) == NULL)        esl_fatal(msg);
  if ((efp = esl_fileparser_Create(fp)) == NULL) esl_fatal(msg);
  if (esl_mixdchlet_Read(efp, &d2) != eslOK)     esl_fatal(msg);
  esl_fileparser_Destroy(efp);
  fclose(fp);

  if (esl_mixdchlet_Compare(d, d2, tol) != eslOK) esl_fatal(msg);

  esl_mixdchlet_Destroy(d2);
  return;
}

static void
utest_inference(ESL_RANDOMNESS *r, ESL_MIXDCHLET *d, int ncounts, int be_verbose)
{
  char   *msg    = "esl_dirichlet: inference unit test failed";
  double *counts = malloc(sizeof(double) * d->K);
  double *probs  = malloc(sizeof(double) * d->K);
  double *iq     = malloc(sizeof(double) * d->N);
  double *ip     = malloc(sizeof(double) * d->K);
  int     qused, qguess;
  int     c, i, q;
  double  maxdeviation;

  /* Sample component, p vector, c vector from mixture Dirichlet */
  qused = esl_rnd_DChoose(r, d->pq, d->N); 
  printf("qused=%1d\n", qused); 
  esl_dirichlet_DSample(r, d->alpha[qused], d->K, probs);
  esl_vec_DSet(counts, d->K, 0.);
  for (c = 0; c < ncounts; c++)
    {
      i = esl_rnd_DChoose(r, probs, d->K);
      counts[i] += 1.;
    }

  /* First inference test: 
   * classify by posterior inference on the sampled probability vector. 
   */
  for (q = 0; q < d->N; q++)
    {
      esl_dirichlet_LogProbProbs(probs, d->alpha[q], d->K, &(iq[q]));
      iq[q] += log(d->pq[q]);
    }
  esl_vec_DLogNorm(iq, d->N);
  qguess = esl_vec_DArgMax(iq, d->N); /* the MP guess from the probs */
  printf("qguess: %1d\n", qguess); 
  if (qused != qguess) esl_fatal(msg);

  /* Second inference test: 
   * classify by posterior inference on the sampled count vector;
   * then attempt to estimate the probability vector.
   */
  esl_mixdchlet_MPParameters(counts, d->K, d, iq, ip);
  qguess = esl_vec_DArgMax(iq, d->N); /* the MP guess from the counts */
  printf("%1d\n", qguess); 
  if (qused != qguess) esl_fatal(msg);

  for (i = 0; i < d->K; i++)
    ip[i] = fabs(ip[i] - probs[i]); /* ip[] is now the differences rel to probs */
  maxdeviation = esl_vec_DMax(ip, d->K);
  printf("maxdev=%.3f\n", maxdeviation);
  if (maxdeviation > 0.05) esl_fatal(msg);

  free(counts);
  free(probs);
  free(iq);
  free(ip);
  return;
}

static void
utest_fit(ESL_RANDOMNESS *r, ESL_MIXDCHLET *d, int ntrials, int ncounts, double tol, int be_verbose)
{
  char           *msg    = "esl_dirichlet: fit unit test failed";
  ESL_MIXDCHLET  *id = NULL;
  double        **counts;
  double         *probs = malloc(sizeof(double) * d->K);
  int             qused;
  int             m;
  int             c;
  int             q;			/* counter over mixture components (0..N-1) */
  int             i;			/* counter over params (0..K-1) */

  counts = malloc(sizeof(double *) * ntrials);
  for (m = 0; m < ntrials; m ++)
    counts[m] = malloc(sizeof(double) * d->K);

  for (m = 0; m < ntrials; m ++) {
    /* Sample component, p vector, c vector from mixture Dirichlet */
    qused = esl_rnd_DChoose(r, d->pq, d->N); 
    esl_dirichlet_DSample(r, d->alpha[qused], d->K, probs);
    esl_vec_DSet(counts[m], d->K, 0.);

    for (c = 0; c < ncounts; c++)
      {
	i = esl_rnd_DChoose(r, probs, d->K);
	counts[m][i] += 1.;
      }
  }
  
  /* Start with a random id, use the counts to infer d by 
   * maximum likelihood gradient descent.
   * Generate a random starting point, alphas range from 0..10. 
   */
  id = esl_mixdchlet_Create(d->N, d->K);
  for (q = 0; q < id->N; q++) {
    id->pq[q] = esl_rnd_UniformPositive(r);

     for (i = 0; i < id->K; i++)
       id->alpha[q][i] = 10.0*esl_rnd_UniformPositive(r);
  }
  esl_vec_DNorm(id->pq, id->N);

  /* optimize id */
  esl_mixdchlet_Fit(counts, ntrials, id, be_verbose);
  printf("\nGiven dirichtlet\n");
  for (q = 0; q < d->N; q++) {
    printf("q[%d] %f\n", q, d->pq[q]);
    for (i = 0; i < d->K; i++)
      printf("alpha[%d][%d] %f\n", q, i, d->alpha[q][i]);
  }
  printf("\nInfered dirichtlet\n");
  for (q = 0; q < id->N; q++) {
    printf("q[%d] %f\n", q, id->pq[q]);
    for (i = 0; i < id->K; i++)
      printf("alpha[%d][%d] %f\n", q, i, id->alpha[q][i]);
  }
  if (esl_mixdchlet_Compare(d, id, tol) != eslOK) esl_fatal(msg);
  
  for (m = 0; m < ntrials; m ++)
    free(counts[m]);
  free(counts);
  free(probs);
  esl_mixdchlet_Destroy(id);

  return;
}

#endif /*eslDIRICHLET_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/



/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef eslDIRICHLET_TESTDRIVE
/*
 * gcc -g -Wall -I. -L. -o esl_dirichlet_utest -DeslDIRICHLET_TESTDRIVE esl_dirichlet.c -leasel -lm
 * ./esl_dirichlet_utest
 */
#include "easel.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_dirichlet.h"
/* Note that the RNG seed of 10 is carefully chosen to make the stochastic 
 * tests work reproducibly. Other choices will tend to fail.
 */
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "10",  NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-t",        eslARG_REAL,  "1e-4",  NULL, NULL,  NULL,  NULL, NULL, "tolerance for real-value equality comparisons",    0 },
  { "-C",        eslARG_INT,      "2",  NULL, NULL,  NULL,  NULL, NULL, "number of components in test mixture D'chlets",    0 },
  { "-K",        eslARG_INT,      "6",  NULL, NULL,  NULL,  NULL, NULL, "alphabet size in test mixture D'chlets",           0 },
  { "-N",        eslARG_INT,   "1000",  NULL, NULL,  NULL,  NULL, NULL, "number of sample counts in mixture D'chlet tests", 0 },
  { "-T",        eslARG_INT,    "100",  NULL, NULL,  NULL,  NULL, NULL, "number of trials of mixture D'chlet tests",        0 },
  { "-v",        eslARG_NONE,    NULL,  NULL, NULL,  NULL,  NULL, NULL, "show verbose output",                              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for dirichlet module";

int
main(int argc, char **argv)
{
  char           *msg          = "esl_dirichlet unit test failed";
  ESL_GETOPTS    *go           = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r            = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_MIXDCHLET  *d            = NULL;
  int             K            = esl_opt_GetInteger(go, "-K");
  int             ncounts      = esl_opt_GetInteger(go, "-N");
  int             ntrials      = esl_opt_GetInteger(go, "-T");
  double          tol          = esl_opt_GetReal   (go, "-t");
  int             be_verbose   = esl_opt_GetBoolean(go, "-v");
  int             t;

  if (be_verbose) printf("rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* Create a two-component mixture Dirichlet for testing */
  if ((d = esl_mixdchlet_Create(2, K)) == NULL) esl_fatal(msg);
  esl_vec_DSet(d->pq,       2, 0.5);
  esl_vec_DSet(d->alpha[0], K, 1.0);
  esl_vec_DSet(d->alpha[1], K, 0.1);

  utest_io(d, tol);
  utest_fit(r, d, ntrials, ncounts, tol, be_verbose);
  for (t = 0; t < ntrials; t++) 
    utest_inference(r, d, ncounts, be_verbose);
 
  esl_randomness_Destroy(r);
  esl_mixdchlet_Destroy(d);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslDIRICHLET_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/



/*****************************************************************
 * 7. Example 
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
#include "easel.h"
#include "esl_random.h"
#include "esl_fileparser.h"
#include "esl_vectorops.h"
#include "esl_dirichlet.h"

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
  r = esl_randomness_Create(0);            /* init the random generator */
  qused = esl_rnd_DChoose(r, pri->pq, pri->N); /* sample a component */
  esl_dirichlet_DSample(r, pri->alpha[qused], pri->K, probs);

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
  esl_randomness_Destroy(r);
  free(probs); free(counts); free(iq); free(ip);
  return 0;
}
/*::cexcerpt::dirichlet_example::end::*/
#endif /*eslDIRICHLET_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
