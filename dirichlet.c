/* dirichlet.c
 * Functions relevant to Dirichlet densities.
 * 
 * SRE, Tue Nov  2 13:42:59 2004 [St. Louis]
 * SVN $Id$
 */

#include <easel/easel.h>
#include <easel/random.h>
#include <easel/gamma.h>

/* Function:  esl_dirichlet_LogProbData()
 * Incept:    SRE, Tue Nov  2 14:22:37 2004 [St. Louis]
 *
 * Purpose:   Given an observed count vector c[0..K-1], 
 *            and a Dirichlet density parameterized by
 *            \alpha[0..K-1];
 *            calculate log P(c | \alpha).
 *            
 *            This is \int P(c | p) P(p | \alpha) dp,
 *            an integral that can be solved analytically.
 *
 * Args:      c          - count vector, [0..K-1]
 *            alpha      - Dirichlet parameters, [0..K-1]
 *            K          - size of c, alpha vectors
 *            ret_answer - RETURN: log P(c | \alpha)
 *
 * Returns:   On success, log P(c | \alpha) is put in *ret_answer;
 *            returns eslOK.
 *
 * Throws:    (no abnormal error conditions)
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
      esl_gamma_log(alpha[x] + c[x], &a1); 
      esl_gamma_log(c[x] + 1.,       &a2);
      esl_gamma_log(alpha[x],        &a3);
      lnp  += a1 - a2 - a3;
    }
  esl_gamma_log(sum1,      &a1);
  esl_gamma_log(sum2,      &a2);
  esl_gamma_log(sum3 + 1., &a3);
  lnp += a2 + a3 - a1;

  *ret_answer = lnp;
  return eslOK;
}


/* Function:  esl_dirichlet_sample()
 * Incept:    SRE, Tue Nov  2 14:30:31 2004 [St. Louis]
 *
 * Purpose:   Given a Dirichlet density parameterized by \alpha[0..K-1],
 *            sample a probability vector p[0..K-1] from
 *            P(p | \alpha).
 *
 * Args:      r      - random number generation object
 *            alpha  - parameters of Dirichlet density [0..K-1]
 *            K      - vector size
 *            p      - RETURN: sampled probability vector
 *                     (caller allocates 0..K-1).         
 *
 * Returns:   eslOK; p contains the sampled vector.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_dirichlet_sample(ESL_RANDOMNESS *r, double *alpha, int K, double *p)
{
  int    x;
  for (x = 0; x < K; x++) esl_gamma_sample(r, alpha[x], &(p[x]));
  esl_vec_DNorm(p, K);
}
