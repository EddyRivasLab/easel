/* ratematrix.c
 * Routines for manipulating evolutionary rate matrices.
 * 
 * SRE, Tue Jul 13 15:51:23 2004 [St. Louis]
 * SVN $Id$
 */

#include <easel/easel.h>
#include <easel/dmatrix.h>
#include <easel/vectorops.h>
#include <easel/ratematrix.h>

/* Function:  esl_ratemx_Symm2Q()
 * Incept:    SRE, Tue Jul 13 15:52:41 2004 [St. Louis]
 *
 * Purpose:   Given a lower triangular matrix (j<i) of 
 *            residue exchangeabilities s, and a statiionary residue
 *            frequency vector pi;
 *            calculates a rate matrix Q.
 *            
 *            Q(j | i) = Q_ij = s_ij * \pi_j
 *            
 * Args:      s     - symmetric residue "exchangeabilities", lower
 *                    triangular, in Easel alphabet order.
 *            pi    - residue frequencies at stationarity, 
 *                    in Easel alphabetic order.
 *            ret_Q - RETURN: rate matrix, square (NxN), 
 *                    in Easel alphabetic order; 
 *                    allocated here, caller must free.          
 *                    
 * Returns:   ESL_OK on success; ret_Q is allocated here.
 * 
 * Xref:      STL8/p56.
 */
int
esl_ratemx_Symm2Q(ESL_DMATRIX *s, double *pi, ESL_DMATRIX **ret_Q)
{
  ESL_DMATRIX *Q;
  int          i,j;
  double       sum;

  if ((Q = esl_dmx_Alloc(s->n, s->n)) == NULL) ESL_ERROR(ESL_EMEM, "malloc failed"); 

  /* Scale all off-diagonals to pi[j] * sij[i][j].
   * Set diagonal to 1 - sum of all j != i.
   */
  for (i = 0; i < s->n; i++)
    {
      for (j = 0; j < i; j++)	/* only look at lower triangle of s. */
	{
	  Q->mx[i][j] = pi[j] * s->mx[i][j]; 
	  Q->mx[j][i] = pi[i] * s->mx[i][j];
	}
      Q->mx[i][i] = 0.;		/* makes the vector sum work for j != i */
      Q->mx[i][i] = 1 - esl_vec_sum(Q->mx[i], Q->n);
    }

  *ret_Q = Q;
  return ESL_OK;
}


/* Function:  esl_ratemx_Normalize()
 * Incept:    SRE, Tue Jul 13 16:05:16 2004 [St. Louis]
 *
 * Purpose:   Normalize a rate matrix Q so that expected substitution
 *            rate per dt is x.
 *
 *            Expected substitution rate is:
 *               \sum_i \sum_j pi_i Q_ij  \forall i \neq j
 *
 *            x typically taken to be 1.0, so time units are substitutions/site.
 *            An exception is PAM, where x = 0.01 for 1 PAM unit.
 *
 * Args:      Q   - rate matrix to normalize
 *            pi  - stationary residue frequencies
 *            x   - expected subsitution rate per dt 
 *                  (1.0 = substitutions/site; 0.01 = PAMs)
 *
 * Returns:   ESL_OK on success;
 *            rate matrix Q is altered.
 *
 * Xref:      STL8/p56.
 */
int
esl_ratemx_Normalize(ESL_DMATRIX *Q, double *pi, double x)
{
  int     i,j;
  double  sum = 0.;

  for (i = 0; i < Q->n; i++)
    for (j = 0; j < Q->n; j++)
      if (i != j) sum += pi[i] * Q->mx[i][j];

  for (i = 0; i < Q->n; i++)
    for (j = 0; j < Q->n; j++)
      Q->mx[i][j] *= (x / sum);

  return ESL_OK;
}
