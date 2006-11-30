/* esl_stats.c
 * Foundation for the statistics modules.
 * 
 * SRE, Tue Jul 19 10:57:44 2005
 * SVN $Id$
 */
#include <esl_config.h>

#include <math.h>

#include <easel.h>
#include <esl_stats.h>


/* Function:  esl_stats_Mean()
 * Incept:    SRE, Tue Jul 19 11:04:00 2005 [St. Louis]
 *
 * Purpose:   Calculates the sample mean and s^2, the unbiased
 *            estimator of the population variance, for a
 *            sample of <n> numbers <x[0]..x[n-1]>, and optionally
 *            returns either or both through <ret_mean> and
 *            <ret_var>.
 *
 * Args:      x        - samples x[0]..x[n-1]
 *            n        - number of samples
 *            ret_mean - optRETURN: mean
 *            ret_var  - optRETURN: estimate of population variance       
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stats_Mean(double *x, int n, double *ret_mean, double *ret_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (ret_mean != NULL)  *ret_mean = sum / (double) n;
  if (ret_var  != NULL)  *ret_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}

/* Function:  esl_stats_LogGamma()
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
 * Throws:    <eslERANGE> if $x <= 0$.
 */
int
esl_stats_LogGamma(double x, double *ret_answer)
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
  if (x <= 0.0)  ESL_EXCEPTION(eslERANGE, "invalid x <= 0 in esl_stats_LogGamma()");

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


/* Function:  esl_stats_Psi()
 * Incept:    SRE, Tue Nov 15 13:57:59 2005 [St. Louis]
 *
 * Purpose:   Computes $\Psi(x)$ (the "digamma" function), which is
 *            the derivative of log of the Gamma function:
 *            $d/dx \log \Gamma(x) = \frac{\Gamma'(x)}{\Gamma(x)} = \Psi(x)$.
 *            Argument $x$ is $> 0$. 
 * 
 *            This is J.M. Bernardo's "Algorithm AS103",
 *            Appl. Stat. 25:315-317 (1976).  
 */
int
esl_stats_Psi(double x, double *ret_answer)
{
  double answer = 0.;
  double x2;

  if (x <= 0.0) ESL_EXCEPTION(eslERANGE, "invalid x <= 0 in esl_stats_Psi()");
  
  /* For small x, Psi(x) ~= -0.5772 - 1/x + O(x), we're done.
   */
  if (x <= 1e-5) {
    *ret_answer = -eslCONST_EULER - 1./x;
    return eslOK;
  }

  /* For medium x, use Psi(1+x) = \Psi(x) + 1/x to c.o.v. x,
   * big enough for Stirling approximation to work...
   */
  while (x < 8.5) {
    answer = answer - 1./x;
    x += 1.;
  }
  
  /* For large X, use Stirling approximation
   */
  x2 = 1./x;
  answer += log(x) - 0.5 * x2;
  x2 = x2*x2;
  answer -= (1./12.)*x2;
  answer += (1./120.)*x2*x2;
  answer -= (1./252.)*x2*x2*x2;

  *ret_answer = answer;
  return eslOK;
}



/* Function: esl_stats_IncompleteGamma()
 * 
 * Purpose:  Returns $P(a,x)$ and $Q(a,x) where:
 *           $P(a,x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt$
 *                  $= \frac{\gamma(a,x)}{\Gamma(a)}$
 *           $Q(a,x) = \frac{1}{\Gamma(a)} \int_{x}^{\infty} t^{a-1} e^{-t} dt$
 *                  $= 1 - P(a,x)$
 *
 *           $P(a,x)$ is the CDF of a gamma density with $\lambda = 1$,
 *           and $Q(a,x)$ is the survival function.
 *           
 *           For $x \simeq 0$, $P(a,x) \simeq 0$ and $Q(a,x) \simeq 1$; and
 *           $P(a,x)$ is less prone to roundoff error. 
 *           
 *           The opposite is the case for large $x >> a$, where
 *           $P(a,x) \simeq 1$ and $Q(a,x) \simeq 0$; there, $Q(a,x)$ is
 *           less prone to roundoff error.
 *
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988. 
 *           
 * Args:     a          - for instance, degrees of freedom / 2     [a > 0]
 *           x          - for instance, chi-squared statistic / 2  [x >= 0] 
 *           ret_pax    - RETURN: P(a,x)
 *           ret_qax    - RETURN: Q(a,x)
 *
 * Return:   <eslOK> on success.
 *
 * Throws:   <eslERANGE> if a or x is out of accepted range.
 *           <eslECONVERGENCE> if approximation fails to converge.
 */          
int
esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax)
{
  int    iter;			/* iteration counter */
  double pax;			/* P(a,x) */
  double qax;			/* Q(a,x) */

  if (a <= 0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): a must be > 0");
  if (x <  0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): x must be >= 0");

  /* For x > a + 1 the following gives rapid convergence;
   * calculate Q(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)},
   * using a continued fraction development for \Gamma(a,x).
   */
  if (x > a+1) 
    {
      double oldp;		/* previous value of p    */
      double nu0, nu1;		/* numerators for continued fraction calc   */
      double de0, de1;		/* denominators for continued fraction calc */

      nu0 = 0.;			/* A_0 = 0       */
      de0 = 1.;			/* B_0 = 1       */
      nu1 = 1.;			/* A_1 = 1       */
      de1 = x;			/* B_1 = x       */

      oldp = nu1;
      for (iter = 1; iter < 100; iter++)
	{
	  /* Continued fraction development:
	   * set A_j = b_j A_j-1 + a_j A_j-2
	   *     B_j = b_j B_j-1 + a_j B_j-2
           * We start with A_2, B_2.
	   */
				/* j = even: a_j = iter-a, b_j = 1 */
				/* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
	  nu0 = nu1 + ((double)iter - a) * nu0;
	  de0 = de1 + ((double)iter - a) * de0;
				/* j = odd: a_j = iter, b_j = x */
				/* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
	  nu1 = x * nu0 + (double) iter * nu1;
	  de1 = x * de0 + (double) iter * de1;
				/* rescale */
	  if (de1 != 0.) 
	    { 
	      nu0 /= de1; 
	      de0 /= de1;
	      nu1 /= de1;
	      de1 =  1.;
	    }
				/* check for convergence */
	  if (fabs((nu1-oldp)/nu1) < 1.e-7)
	    {
	      esl_stats_LogGamma(a, &qax);	      
	      qax = nu1 * exp(a * log(x) - x - qax);

	      if (ret_pax != NULL) *ret_pax = 1 - qax;
	      if (ret_qax != NULL) *ret_qax = qax;
	      return eslOK;
	    }

	  oldp = nu1;
	}
      ESL_EXCEPTION(eslECONVERGENCE,
		"esl_stats_IncompleteGamma(): fraction failed to converge");
    }
  else /* x <= a+1 */
    {
      double p;			/* current sum               */
      double val;		/* current value used in sum */

      /* For x <= a+1 we use a convergent series instead:
       *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
       * where
       *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
       * which looks appalling but the sum is in fact rearrangeable to
       * a simple series without the \Gamma functions:
       *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
       * and it's obvious that this should converge nicely for x <= a+1.
       */
      p = val = 1. / a;
      for (iter = 1; iter < 10000; iter++)
	{
	  val *= x / (a+(double)iter);
	  p   += val;
	  
	  if (fabs(val/p) < 1.e-7)
	    {
	      esl_stats_LogGamma(a, &pax);
	      pax = p * exp(a * log(x) - x - pax);

	      if (ret_pax != NULL) *ret_pax = pax;
	      if (ret_qax != NULL) *ret_qax = 1. - pax;
	      return eslOK;
	    }
	}
      ESL_EXCEPTION(eslECONVERGENCE,
		"esl_stats_IncompleteGamma(): series failed to converge");
    }
  /*NOTREACHED*/
  return eslOK;
}


/* Function:  esl_stats_ChiSquaredTest()
 * Incept:    SRE, Tue Jul 19 11:39:32 2005 [St. Louis]
 *
 * Purpose:   Calculate the probability that a chi-squared statistic
 *            with <v> degrees of freedom would exceed the observed
 *            chi-squared value <x>; return it in <ret_answer>. If
 *            this probability is less than some small threshold (say,
 *            0.05 or 0.01), then we may reject the hypothesis we're
 *            testing.
 *
 * Args:      v          - degrees of freedom
 *            x          - observed chi-squared value
 *            ret_answer - RETURN: P(\chi^2 > x)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslERANGE> if <v> or <x> are out of valid range.
 *            <eslECONVERGENCE> if iterative calculation fails.
 */
int
esl_stats_ChiSquaredTest(int v, double x, double *ret_answer)
{
  return esl_stats_IncompleteGamma((double)v/2, x/2, NULL, ret_answer);
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
