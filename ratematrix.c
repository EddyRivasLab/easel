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
 *            residue exchangeabilities s, and a stationary residue
 *            frequency vector pi; assuming s_ij = s_ji;
 *            calculates a rate matrix Q.
 *            
 *            Q(j | i) = Q_ij = s_ij * \pi_j
 *            
 *            The resulting Q is not normalized to any particular
 *            number of substitutions/site/time unit; see 
 *            esl_ratemx_Normalize() for that.
 *            
 * Args:      s     - symmetric residue "exchangeabilities";
 *                    only lower triangular entries are used.
 *            pi    - residue frequencies at stationarity. 
 *            Q     - RETURN: rate matrix, square (NxN). 
 *                    Caller allocates the memory for this.
 *                    
 * Returns:   ESL_OK on success; Q is calculated.
 * 
 * Xref:      STL8/p56.
 */
int
esl_ratemx_Symm2Q(ESL_DMATRIX *s, double *pi, ESL_DMATRIX *Q)
{
  int          i,j;
  double       sum;

  /* Scale all off-diagonals to pi[j] * sij[i][j].
   */
  for (i = 0; i < s->n; i++)
    for (j = 0; j < i; j++)	/* only look at lower triangle of s. */
      {
	Q->mx[i][j] = pi[j] * s->mx[i][j]; 
	Q->mx[j][i] = pi[i] * s->mx[i][j];
      }

  /* Set diagonal to  -\sum of all j != i.
   */
  for (i = 0; i < s->n; i++)
    {
      Q->mx[i][i] = 0.;		/* makes the vector sum work for j != i */
      Q->mx[i][i] = -1. * esl_vec_DSum(Q->mx[i], Q->n);
    }

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



/* Function:  esl_ratemx_TaylorExp()
 * Incept:    SRE, Thu Aug 12 08:45:26 2004 [St. Louis]
 *
 * Purpose:   Given rate matrix Q and time t; 
 *            also given pre-allocated space for result matrix P;
 *            calculates conditional substitution matrix
 *            P=e^{tQ}, w/ values P(y | x, t).
 *            
 *            Uses Taylor series approximation:
 *            
 *            e^{tQ} = \sum_{k=0}^{\infty} \frac{t^k Q^k} {k!} =
 *            
 *                               t^2 Q^2   t^3 Q^3
 *                   = I + tQ +  ------- + ------- ...
 *                                   2        3!
 *
 * Args:      Q    - rate matrix to exponentiate
 *            t    - time units
 *            P    - RETURN: substitution matrix       
 *
 * Returns:   ESL_OK on success; conditional probabilities are in P.
 *
 * Xref:      
 */
int
esl_ratemx_TaylorExp(ESL_DMATRIX *Q, double t, ESL_DMATRIX *P)
{
  ESL_DMATRIX *tmp;             /* keeps running product Q^k */
  ESL_DMATRIX *C;
  double       factor;
  int          k;

  if ((tmp = esl_dmx_Alloc(Q->n, Q->n)) == NULL)  ESL_ERROR(ESL_EMEM, "malloc failed");
  if ((C   = esl_dmx_Alloc(Q->n, Q->n)) == NULL)  ESL_ERROR(ESL_EMEM, "malloc failed");
  
  esl_dmx_SetIdentity(P);
  factor = 1;
  esl_dmx_Copy(Q, tmp);		/* tmp is now = Q */

  /* WARNING: no convergence test here. arbitrarily taking Taylor out
   * through 100 terms. Don't leave this this way forever.
   */
  for (k = 1; k < 100; k++)
    {
      factor *= t/k;
      esl_dmx_AddScale(P, factor, tmp);    /* P += factor*tmp */
      esl_dmx_Multiply(tmp, Q, C);         /* C = tmp*Q */
      esl_dmx_Copy(C, tmp);	           /* tmp = C = Q^{k+1} */
    }

  esl_dmx_Free(tmp);
  esl_dmx_Free(C);
  return ESL_OK;
}


/* Function:  esl_ratemx_CreateHKY()
 * Incept:    SRE, Thu Aug 12 08:26:39 2004 [St. Louis]
 *
 * Purpose:   Given base composition f[{ACGT}] and transition/
 *            transversion relative rates \alpha and \beta;
 *            allocate and return an HKY (Hasegawa/Kishino/Yano)
 *            DNA rate matrix, normalized to a unit of
 *            1t= 1.0 substitutions/site.
 *            
 *            Ref: [Hasegawa85]
 *
 * Args:      f      - stationary base composition A..T
 *            alpha  - relative transition rate
 *            beta   - relative transversion rate
 *
 * Returns:   Q      - allocated HKY rate matrix.
 *
 * Xref:      
 */
ESL_DMATRIX *
esl_ratemx_CreateHKY(double *f, double alpha, double beta)
{
  ESL_DMATRIX *Q;
  int i,j;

  if ((Q = esl_dmx_Alloc(4, 4)) == NULL)
    ESL_ERROR_VAL(NULL, ESL_EMEM, "malloc failed");   
  
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
	{
	  if (i != j)  Q->mx[i][j] = ((i+j)%2)? f[j]*beta : f[j]*alpha; /* even=0=transition;odd=1=transversion */
	  else         Q->mx[i][j] = 0.;
	}
      Q->mx[i][i] =  -1. * esl_vec_DSum(Q->mx[i], 4);
    }
  esl_ratemx_Normalize(Q, f, 1.0);
  return Q;
}
