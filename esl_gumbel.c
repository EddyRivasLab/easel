/* esl_gumbel.c
 * Statistical routines for Gumvel (type I extreme value) distributions.
 * 
 * SRE, Thu Jun 23 11:48:39 2005
 * SVN $Id$
 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <easel.h>
#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_MINIMIZER
#include <esl_minimizer.h>
#endif
#include <esl_vectorops.h>
#include <esl_gumbel.h>

static void mean_and_variance(double *x, int n, double *ret_mean, double *ret_var);


/*****************************************************************
 * Routines for evaluating densities and distributions
 *****************************************************************/ 

/* Function:  esl_gumbel_pdf()
 * Incept:    SRE, Sun Jun 26 14:08:19 2005 [St. Louis]
 *
 * Purpose:   Calculates the probability density function for the Gumbel,
 *            $P(X=x)$, given quantile <x> and Gumbel location and
 *            scale parameters <mu> and <lambda>.
 *            
 *            Let $y = \lambda(x-\mu)$; for 64-bit doubles,
 *            useful dynamic range is about $-6.5 <= y <= 710$.
 *            Returns 0.0 for smaller $y$, 0.0 for larger $y$.
 */
double
esl_gumbel_pdf(double x, double mu, double lambda)
{
  double y;
  y = lambda * (x - mu);
  return (lambda * exp(-y - exp(-y)));
}


/* Function:  esl_gumbel_logpdf()
 * Incept:    SRE, Sun Jun 26 14:08:19 2005 [St. Louis]
 *
 * Purpose:   Calculates the log probability density function for the Gumbel,
 *            $\log P(X=x)$.
 *                                                     
 *            Let $y = \lambda(x-\mu)$; for 64-bit doubles,
 *            useful dynamic range is about $-708 <= y <= \infty$.
 *            Returns $-\infty$ for smaller or larger $y$.
 */
double
esl_gumbel_logpdf(double x, double mu, double lambda)
{
  double y;
  y = lambda * (x - mu);
  return (log(lambda) -y - exp(-y));
}


/* Function:  esl_gumbel_cdf()
 * Incept:    SRE, Sun Jun 26 10:18:51 2005 [St. Louis]
 *
 * Purpose:   Calculates the cumulative distribution function
 *            for the Gumbel, $P(X \leq x)$.
 *            
 *            Let $y = \lambda(x-\mu)$; for 64-bit doubles,
 *            useful dynamic range for $y$ is about $-6.5 <= y <=36$.
 *            Returns 0.0 for smaller $y$, 1.0 for larger $y$.
 */
double 
esl_gumbel_cdf(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);
  return exp(-exp(-y));
}

/* Function:  esl_gumbel_logcdf()
 * Incept:    SRE, Sun Jun 26 10:18:51 2005 [St. Louis]
 *
 * Purpose:   Calculates the log of the cumulative distribution function
 *            for the Gumbel, $\log P(X \leq x)$.
 *            
 *            Let $y = \lambda(x-\mu)$; for 64-bit doubles,
 *            useful dynamic range for $y$ is about $-708 <= y <= 708$.
 *            Returns $-\infty$ for smaller $y$, 0.0 for larger $y$.
 */
double 
esl_gumbel_logcdf(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);
  return (-exp(-y));
}

/* Function:  esl_gumbel_surv()
 * Incept:    SRE, Sun Jun 26 09:54:31 2005 [St. Louis]
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ for a Gumbel 
 *            (that is, 1-cdf), the right tail's probability mass.
 * 
 *            Let $y=\lambda(x-\mu)$; for 64-bit doubles, 
 *            useful dynamic range for $y$ is $-3.6 <= y <= 708$.
 *            Returns 1.0 for $y$ below lower limit, and 0.0
 *            for $y$ above upper limit.
 */
double
esl_gumbel_surv(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);

  /* Near P=0.0 we'll use the approx P(S>=x) ~ e^-y, for "large" y,
   * because 1-e^-x ~ x for small x.  What's a sufficiently large y to
   * use P(S>=x) ~ e^-y?  For sure, we need to use the approx if
   * exp(-e^-y) > 1-epsilon, that is, y > -log(epsilon), which is
   * about 36. But the accuracy of the full calculation starts
   * breaking down even lower than that; 0.5 is an arbitrary factor,
   * tested empirically, where the two approaches give the most
   * similar answers (to within about 10^-8). This crossover is at y >
   * 18.
   */
  if (y > -0.5*log(DBL_EPSILON)) return exp(-y);
  return (1 - exp(-exp(-y)));
}

/* Function:  esl_gumbel_logsurv()
 * Incept:    SRE, Sun Jun 26 13:45:52 2005 [St. Louis]
 *
 * Purpose:   Calculates $\log P(X>x)$ for a Gumbel (that is, $\log$(1-cdf)):
 *            the log of the right tail's probability mass.
 * 
 *            Let $y=\lambda(x-\mu)$; for 64-bit doubles, 
 *            useful dynamic range for $y$ is $-6.5 <= y <= \infty$.
 *            Returns 0.0 for smaller $y$.
 */
double
esl_gumbel_logsurv(double x, double mu, double lambda)
{
  double y;
  y = lambda*(x-mu);

   /* For "large" y, we can use 1-e^-a = a to get log P(S>=y) ~ -y;
   * .5 log(DBL_EPSILON) is an arbitrary crossover for 64-bit IEEE doubles.
   */
  if (y > -0.5*log(DBL_EPSILON)) return -y;

  /* For "small y", we can use log(1-x) ~ -x. 
   * -2.9 is an arbitrary crossover, tested for 64-bit IEEE doubles
   */
  if (y < -2.9) return (-exp(-exp(-y)));

   return log(1-exp(-exp(-y)));
}
/*------------------ end of densities and distributions --------------------*/




/*****************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 *****************************************************************/ 

#ifdef eslAUGMENT_RANDOM
/* Function:  esl_gumbel_Sample()
 * Incept:    SRE, Thu Jun 23 11:38:39 2005 [St. Louis]
 *
 * Purpose:   Sample a Gumbel-distributed random variate
 *            by the transformation method.
 */
double
esl_gumbel_Sample(ESL_RANDOMNESS *r, double mu, double lambda)
{
  double p;
  p = esl_rnd_UniformPositive(r); 
  return mu - log(-1. * log(p)) / lambda;
} 
#endif /*eslAUGMENT_RANDOM*/

/*------------------------ end of sampling --------------------------------*/



/*****************************************************************
 * Routines for maximum likelihood fitting Gumbels to data
 * (fitting truncated distributions requires augmentation w/ minimizer module)
 *****************************************************************/ 

/*****************************************************************
 * Complete data, maximum a posteriori parameters
 *****************************************************************/ 

/* lawless416()
 * SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 * 
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to Gumbel lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *           
 * Args:     x      - array of sample values 
 *           n      - number of samples 
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *           
 * Return:   (void)
 */ 
static void
lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df)
{
  double esum;			/* \sum e^(-lambda xi)      */
  double xesum;			/* \sum xi e^(-lambda xi)   */
  double xxesum;		/* \sum xi^2 e^(-lambda xi) */
  double xsum;			/* \sum xi                  */
  int i;

  esum = xesum = xsum  = xxesum = 0.;
  for (i = 0; i < n; i++)
    {
      xsum   += x[i];
      xesum  += x[i] * exp(-1. * lambda * x[i]);
      xxesum += x[i] * x[i] * exp(-1. * lambda * x[i]);
      esum   += exp(-1. * lambda * x[i]);
    }
  *ret_f  = (1./lambda) - (xsum / n)  + (xesum / esum);
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));
}

/* Function: esl_gumbel_FitComplete()
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis] - HMMER's EVDMaxLikelyFit()
 * 
 * Purpose:  Given an array of Gumbel-distributed samples <x[0]..x[n-1]>,
 *           find maximum likelihood parameters <mu> and <lambda>.
 *           
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *            for lambda using Newton/Raphson iterations,
 *            then substitutes lambda into Lawless' equation 4.1.5
 *            to get mu. 
 *           
 * Args:     x          - list of Gumbel distributed samples
 *           n          - number of samples
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslECONVERGENCE> if the fit doesn't converge.
 */
int
esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda)
{
  double  mean, variance;
  double  lambda, mu;
  double  fx;			/* f(x)  */
  double  dfx;			/* f'(x) */
  double  esum;                 /* \sum e^(-lambda xi) */ 
  double  tol = 1e-5;
  int     i;

  /* 1. Find an initial guess at lambda
   *    (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
   */
  mean_and_variance(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);

  /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      lawless416(x, n, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

  /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      double left, right, mid;
      ESL_DPRINTF1(("esl_gumbel_FitComplete(): Newton/Raphson failed; switchover to bisection"));

      /* First bracket the root */
      left  = 0.;	                 	/* for sure */
      right = eslCONST_PI / sqrt(6.*variance);  /* an initial guess */
      lawless416(x, n, lambda, &fx, &dfx);
      while (fx > 0.) 
	{		
	  right *= 2.;		/* arbitrary leap to the right */
	  if (right > 100.) /* no reasonable lambda should be > 100, we assert */
	    ESL_ERROR(eslECONVERGENCE, "Failed to bracket root in esl_gumbel_FitComplete().");
	  lawless416(x, n, right, &fx, &dfx);
	}

      /* Now, bisection search in left/right interval */
      for (i = 0; i < 100; i++)
	{
	  mid = (left + right) / 2.; 
	  lawless416(x, n, mid, &fx, &dfx);
	  if (fabs(fx) < tol) break;             /* success */
	  if (fx > 0.)	left = mid;
	  else          right = mid;
	}
      if (i == 100) 
	ESL_ERROR(eslECONVERGENCE, "Even bisection search failed in esl_gumbel_FitComplete().");

      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.1.5 to find mu
   */
  esum = 0.;
  for (i = 0; i < n; i++)
    esum  += exp(-1 * lambda * x[i]);
  mu = -1. * log(esum / n) / lambda;

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return 1;
}

/* direct_mv_fit()
 * SRE, Wed Jun 29 08:23:47 2005
 * 
 * Purely for curiousity: a complete data fit using the
 * simple direct method, calculating mu and lambda from mean
 * and variance.
 */
static int
direct_mv_fit(double *x, int n, double *ret_mu, double *ret_lambda)
{
  double mean, variance;

  mean_and_variance(x, n, &mean, &variance);
  *ret_lambda = eslCONST_PI / sqrt(6.*variance);
  *ret_mu     = mean - 0.57722/(*ret_lambda);
  return eslOK;
}


/*------------------- end of complete data fit ---------------------------------*/


/*****************************************************************
 * Censored data, MAP/ML parameters
 *****************************************************************/ 
/* lawless422()
 * SRE, Mon Nov 17 09:42:48 1997 [St. Louis]
 * 
 * Purpose:  Equation 4.2.2 from [Lawless82], pg. 169, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to Gumbel lambda parameter
 *           for Type I censored data. 
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *           
 * Args:     x      - array of observed sample values 
 *           n      - number of observed samples 
 *           z      - number of censored samples = N-n 
 *           phi    - censoring value; all observed x_i >= phi         
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.2.2 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.2.2 evaluated at lambda
 *           
 * Return:   (void)
 */ 
static void
lawless422(double *x, int n, int z, double phi,
	   double lambda, double *ret_f, double *ret_df)
{
  double esum;			/* \sum e^(-lambda xi)      + z term    */
  double xesum;			/* \sum xi e^(-lambda xi)   + z term    */
  double xxesum;		/* \sum xi^2 e^(-lambda xi) + z term    */
  double xsum;			/* \sum xi                  (no z term) */
  int i;

  esum = xesum = xsum  = xxesum = 0.;
  for (i = 0; i < n; i++)
    {
      xsum   += x[i];
      esum   +=               exp(-1. * lambda * x[i]);
      xesum  +=        x[i] * exp(-1. * lambda * x[i]);
      xxesum += x[i] * x[i] * exp(-1. * lambda * x[i]);
    }

  /* Add z terms for censored data
   */
  esum   += (double) z *             exp(-1. * lambda * phi);
  xesum  += (double) z * phi *       exp(-1. * lambda * phi);
  xxesum += (double) z * phi * phi * exp(-1. * lambda * phi);

  *ret_f  = 1./lambda - xsum / n + xesum / esum;
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));

  return;
}

/* Function: esl_gumbel_FitCensored()
 * Date:     SRE, Mon Nov 17 10:01:05 1997 [St. Louis]
 * 
 * Purpose: Given a left-censored array of Gumbel-distributed samples
 *          <x[0]..x[n-1]>, the number of censored samples <z>, and the
 *          censoring value <phi> (all <x[i]> $\geq$ <phi>).
 *          Find maximum likelihood parameters <mu> and <lambda>.
 *           
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *            for lambda using Newton/Raphson iterations;
 *            then substitutes lambda into Lawless' equation 4.2.3
 *            to get mu. 
 *           
 * Args:     x          - array of Gumbel-distributed samples, 0..n-1
 *           n          - number of observed samples
 *           z          - number of censored samples
 *           phi        - censoring value (all x_i >= phi)
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *           
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslECONVERGENCE> if the fit doesn't converge.
 */
int
esl_gumbel_FitCensored(double *x, int n, int z, double phi, 
		       double *ret_mu, double *ret_lambda)
{
  double mean, variance;
  double lambda, mu;
  double fx;			/* f(x)  */
  double dfx;			/* f'(x) */
  double esum;                  /* \sum e^(-lambda xi) */ 
  double tol = 1e-5;
  int    i;

  /* 1. Find an initial guess at lambda
   *    (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
   */
  mean_and_variance(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);

  /* 2. Use Newton/Raphson to solve Lawless 4.2.2 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      lawless422(x, n, z, phi, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

 /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume (!?) that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      double left, right, mid;
      ESL_DPRINTF1(("esl_gumbel_FitCensored(): Newton/Raphson failed; switched to bisection"));

      /* First bracket the root */
      left  = 0.;		/* we know that's the left bound */
      right = eslCONST_PI / sqrt(6.*variance); /* start from here, move "right"... */
      lawless422(x, n, z, phi, right, &fx, &dfx);
      while (fx > 0.)
	{
	  right *= 2.;
	  if (right > 100.) /* no reasonable lambda should be > 100, we assert */
	    ESL_ERROR(eslECONVERGENCE, "Failed to bracket root in esl_gumbel_FitCensored().");
	  lawless422(x, n, z, phi, right, &fx, &dfx);
	}

      /* Now we bisection search in left/right interval */
      for (i = 0; i < 100; i++)
	{
	  mid = (left + right) / 2.; 
	  lawless422(x, n, z, phi, mid, &fx, &dfx);
	  if (fabs(fx) < tol) break;             /* success */
	  if (fx > 0.)	left = mid;
	  else          right = mid;
	}
      if (i == 100) 
	ESL_ERROR(eslECONVERGENCE, "Even bisection search failed in esl_gumbel_FitCensored().");
      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.2.3 to find mu
   */
  esum = 0.;
  for (i = 0; i < n; i++)
    esum  += exp(-lambda * x[i]);
  esum += z * exp(-1. * lambda * phi);    /* term from censored data */
  mu = -log(esum / n) / lambda;        

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return 1;
}


/*****************************************************************
 * Truncated data, MAP parameters (requires minimizer augmentation)
 *****************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER
/* Easel's conjugate gradient descent code allows a single void ptr to
 * point to any necessary fixed data, so we'll put everything into one
 * structure:
 */
struct tevd_data {
  double *x;	/* data: n observed samples from a Gumbel */
  int     n;	/* number of observed samples */
  double  phi;	/* truncation threshold: all observed x_i >= phi */
};

/* tevd_func()
 * 
 * Called by the optimizer: evaluate the objective function
 * for the negative posterior log probability of a particular choice 
 * of parameters mu and lambda, given truncated Gumbel samples.
 */
static double 
tevd_func(double *p, int nparam, void *dptr)
{
  double mu, w, lambda;
  struct tevd_data *data;
  double *x;
  int     n;
  double  phi;
  double  logL;
  int     i;
  
  /* unpack what the optimizer gave us; nparam==2 always
   */
  mu     = p[0];
  w      = p[1];
  lambda = exp(w);
  data   = (struct tevd_data *) dptr;
  x      = data->x;
  n      = data->n;
  phi    = data->phi;

  /* The log likelihood equation
   */
  logL   = n * log(lambda);
  for (i = 0; i < n; i++) 
    logL -= lambda * (x[i] - mu);
  for (i = 0; i < n; i++)
    logL -= exp(-1. * lambda * (x[i] - mu));
  logL -= n * esl_gumbel_logsurv(phi, mu, lambda);    

  return -1.0 * logL;		/* objective: minimize the NLP */
}

/* tevd_grad()
 * 
 * Called by the optimizer: evaluate the gradient of the objective
 * function (the negative posterior log probability of the parameters
 * mu and w, where w = log(lambda), at a particular choice of mu and
 * lambda.
 */
static void
tevd_grad(double *p, int nparam, void *dptr, double *dp)
{
  double mu, lambda, w;
  struct tevd_data *data;
  double *x;
  int     n;
  double  phi;
  double  dmu, dw;
  double  coeff;
  int     i;
  
  /* unpack what the optimizer gave us; nparam==2 always
   */
  mu     = p[0];
  w      = p[1];
  lambda = exp(w);
  data   = (struct tevd_data *) dptr;
  x      = data->x;
  n      = data->n;
  phi    = data->phi;

  /* Both partials include a coefficient that
   * basically looks like P(S=phi) / P(S>=phi); pre-calculate it.
   * Watch out when phi >> mu, which'll give us 0/0; instead,
   * recognize that for phi >> mu, coeff converges to \lambda.
   */
  if (lambda*(phi-mu) > 50.)	/* arbitrary crossover. */
    coeff = lambda;
  else
    coeff = esl_gumbel_pdf(phi, mu, lambda) / esl_gumbel_surv(phi, mu, lambda); 

  /* Partial derivative w.r.t. mu.
   */
  dmu = n * lambda;
  for (i = 0; i < n; i++) 
    dmu -= lambda * exp(-1. * lambda * (x[i] - mu));
  dmu -= n * coeff;    /*!!*/

  /* Partial derivative w.r.t. w=log(lambda).
   */
  dw = n;
  for (i = 0; i < n; i++) dw -= (x[i] - mu) * lambda;
  for (i = 0; i < n; i++) dw += (x[i] - mu) * lambda * exp(-1. * lambda * (x[i] - mu));
  dw += n * (phi - mu) * coeff;    /*!!*/

  /* Return the negative, because we're minimizing NLP, not maximizing.
   */
  dp[0] = -1. * dmu;	/* negative because we're minimizing NLP, not maximizing */
  dp[1] = -1. * dw;
  return;
}
  
/* Function:  esl_gumbel_FitTruncated()
 * Incept:    SRE, Wed Jun 29 14:14:17 2005 [St. Louis]
 *
 * Purpose:   Given a left-truncated array of Gumbel-distributed
 *            samples <x[0]..x[n-1]> and the truncation threshold
 *            <phi> (such that all <x[i]> $\geq$ <phi>).
 *            Find maximum likelihood parameters <mu> and <lambda>.
 *            
 *            <phi> should not be much greater than <mu>, the
 *            mode of the Gumbel, or the fit will become unstable
 *            or may even fail to converge. The problem is
 *            that for <phi> $>$ <mu>, the tail of the Gumbel
 *            becomes a scale-free exponential, and <mu> becomes
 *            undetermined.
 *            
 * Algorithm: Uses conjugate gradient descent to optimize the
 *            log likelihood of the data. Follows a general
 *            approach to fitting missing data problems outlined
 *            in [Gelman95].
 *
 * Args:      x          - observed data samples [0..n-1]
 *            n          - number of samples
 *            phi        - truncation threshold; all x[i] >= phi
 *            ret_mu     - RETURN: ML estimate of mu       
 *            ret_lambda - RETURN: ML estimate of lambda
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECONVERGENCE> if the fit doesn't converge.
 */
int
esl_gumbel_FitTruncated(double *x, int n, double phi, 
			double *ret_mu, double *ret_lambda)
{
  struct tevd_data data;
  double wrk[8];		/* workspace for CG: 4 tmp vectors of size 2 */
  double p[2];			/* mu, w;  lambda = e^w */
  double u[2];			/* max initial step size for mu, lambda */
  int    status;
  double mean, variance;
  double mu, lambda;
  double fx;
  
  data.x   = x;
  data.n   = n;
  data.phi = phi;

  /* The source of the following magic is Evans/Hastings/Peacock, 
   * Statistical Distributions, 3rd edition (2000), p.86, which gives
   * eq's for the mean and variance of a Gumbel in terms of mu and lambda;
   * we turn them around to get mu and lambda in terms of the mean and variance.
   * These would be reasonable estimators if we had a full set of Gumbel
   * distributed variates. They'll be off for a truncated sample, but
   * close enough to be a useful starting point.
   */
  mean_and_variance(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);
  mu     = mean - 0.57722/lambda;

  p[0] = mu;
  p[1] = log(lambda);		/* c.o.v. because lambda is constrained to >0 */

  u[0] = 2.0;
  u[1] = 0.1;

  /* Pass the problem to the optimizer. The work is done by the
   * equations in tevd_func() and tevd_grad().
   */
  status = esl_min_ConjugateGradientDescent(p, u, 2, 
					    &tevd_func, &tevd_grad,(void *)(&data),
					    1e-4, wrk, &fx);
  
  *ret_mu     = p[0];
  *ret_lambda = exp(p[1]);	/* reverse the c.o.v. */
  return status;
}
#endif /*eslAUGMENT_MINIMIZER*/

/* mean_and_variance()
 * 
 * Return the mean and s^2, the unbiased estimator
 * of the population variance, for a sample of numbers <x>.
 */
static void
mean_and_variance(double *x, int n, double *ret_mean, double *ret_var)
{
  double sum = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) { 
    sum   += x[i];
    sqsum += x[i]*x[i];
  }
  *ret_mean = sum / (double) n;
  *ret_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
}
/*------------------------ end of fitting --------------------------------*/



/*****************************************************************
 * Example, test, and stats drivers
 *****************************************************************/ 
/* Example main()
 */
#ifdef eslGUMBEL_EXAMPLE
/*::cexcerpt::gumbel_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslGUMBEL_EXAMPLE -DeslAUGMENT_RANDOM esl_gumbel.c esl_random.c esl_vectorops.c easel.c -lm
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r = esl_randomness_CreateTimeseeded();;
  int     n         = 10000; 	/* simulate 10,000 samples */
  double  mu        = -20.0;       /* with mu = -20 */ 
  double  lambda    = 0.4;         /* and lambda = 0.4 */
  double  min       =  9999.;
  double  max       = -9999.;
  double *x         = malloc(sizeof(double) * n);
  double  z, est_mu, est_lambda;
  int     i;

  for (i = 0; i < n; i++)	/* generate the 10,000 samples */
    { 
      x[i] = esl_gumbel_Sample(r, mu, lambda);
      if (x[i] < min) min = x[i];
      if (x[i] > max) max = x[i];
    }

  z = esl_gumbel_surv(max, mu, lambda);           /* right tail p~1e-4 >= max */
  printf("max = %6.1f  P(>=max) = %g\n", max, z);
  z = esl_gumbel_cdf(min, mu, lambda);             /* left tail p~1e-4 < min */
  printf("min = %6.1f  P(<min)  = %g\n", min, z);

  esl_gumbel_FitComplete(x, n, &est_mu, &est_lambda); /* fit params to the data */
  
  z = 100. * fabs((est_mu - mu) / mu);
  printf("Parametric mu     = %6.1f.  Estimated mu     = %6.2f.  Difference = %.1f%%.\n",
	 mu, est_mu, z);
  z = 100. * fabs((est_lambda - lambda) /lambda);
  printf("Parametric lambda = %6.1f.  Estimated lambda = %6.2f.  Difference = %.1f%%.\n",
	 lambda, est_lambda, z);

  free(x);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::gumbel_example::end::*/
#endif /*eslGUMBEL_EXAMPLE*/


/* Test driver's main()
 */
#ifdef eslGUMBEL_TESTDRIVE
/* compile: gcc -g -Wall -I. -o test -DeslGUMBEL_TESTDRIVE -DeslAUGMENT_RANDOM -DeslAUGMENT_MINIMIZER esl_gumbel.c esl_random.c esl_minimizer.c esl_vectorops.c easel.c -lm
 * run:     ./test
 */
#include <stdio.h>

#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r;
  int    totalN;
  int    n;
  double phi;		/* truncation threshold. */
  int    i;
  double *x;
  double  mu, lambda;
  double  est_mu, est_lambda;

  totalN = 10000;
  mu     = -20.;
  lambda = 0.4;
  phi    = -20.;

  r = esl_randomness_Create(42); /* make the sims reproducible */
  x = malloc(sizeof(double) * totalN);
  
  /* Test complete data fitting on simulated data.
   * Don't tolerate more than 1% error in mu, 3% in lambda.
   */
  for (i = 0; i < totalN; i++)
    x[i] = esl_gumbel_Sample(r, mu, lambda);
  esl_gumbel_FitComplete(x, totalN, &est_mu, &est_lambda);
  if (fabs((est_mu    -mu)    /mu)     > 0.01) abort();
  if (fabs((est_lambda-lambda)/lambda) > 0.03) abort();

  /* Test censored fitting on simulated data, for 
   * the right tail mass above the mode.
   * Don't tolerate more than 1% error in mu, 4% in lambda.
   */
  for (n=0, i = 0; i < totalN; i++)
    if (x[i] >= phi) x[n++] = x[i];
  esl_gumbel_FitCensored(x, n, totalN-n, phi, &est_mu, &est_lambda);
  if (fabs((est_mu    -mu)    /mu)     > 0.01) abort();
  if (fabs((est_lambda-lambda)/lambda) > 0.04) abort();

  /* Test truncated fitting on simulated data.
   * Don't tolerate more than 5% error in mu, 8% in lambda.
   */
#ifdef eslAUGMENT_MINIMIZER
  esl_gumbel_FitTruncated(x, n, phi, &est_mu, &est_lambda);
  if (fabs((est_mu    -mu)    /mu)     > 0.05) abort();
  if (fabs((est_lambda-lambda)/lambda) > 0.08) abort();
#endif /*eslAUGMENT_MINIMIZER*/


  free(x);
  esl_randomness_Destroy(r);
  return 0;
}
#endif /*eslGUMBEL_TESTDRIVE*/


/* Stats driver's main()
 */
#ifdef eslGUMBEL_STATS
/* compile: gcc -g -Wall -I. -o stats -DeslGUMBEL_STATS -DeslAUGMENT_RANDOM -DeslAUGMENT_MINIMIZER esl_gumbel.c esl_random.c esl_minimizer.c esl_vectorops.c easel.c -lm
 * run:     ./stats > stats.out
 * process w/ lines like:
 *    grep "complete    100" stats.out | awk '{$i = 100*($5-$4)/$4; if ($i < 0) $i = -$i; print $i}' | avg
 *    grep "complete    100" stats.out | awk '{$i = 100*($7-$6)/$6; if ($i < 0) $i = -$i; print $i}' | avg
 * to get accuracy summary (in %) for mu, lambda; first part of the grep pattern may be "complete", "censored", or
 * "truncated", second part may be "    100", "   1000", "  10000", or " 100000".
 *
 * This is the routine that collects the accuracy statistics that appear
 * in tables in the Gumbel chapter of the guide, esl_gumbel.tex.
 */
#include <stdio.h>

#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r;
  int    totalN[4] = {100, 1000, 10000, 100000}; /*biggest last; one malloc*/
  int    nexps = 4;
  int    exp;
  int    trial, ntrials;
  double phi;		/* truncation threshold. */
  int    i;
  int    n;
  double *x;
  double  mu, lambda;
  double  est_mu, est_lambda;
  double  val;
  int     do_complete, do_censored, do_truncated;

  ntrials = 500;
  mu      = -20.0;
  lambda  = 0.4;
  phi     = -15.;

  do_complete  = FALSE;		/* Flip these on/off as desired */
  do_censored  = FALSE;
  do_truncated = TRUE;

  r = esl_randomness_CreateTimeseeded();
  x = malloc(sizeof(double) * totalN[nexps-1]);
  
  /* Fitting to simulated complete datasets
   */
  if (do_complete) {
    for (exp = 0; exp < nexps; exp++)
      {
	for (trial = 0; trial < ntrials; trial++)
	  {
	    for (i = 0; i < totalN[exp]; i++)
	      x[i] = esl_gumbel_Sample(r, mu, lambda);

	    /*direct_mv_fit(x, totalN[exp], &est_mu, &est_lambda);*/
	    esl_gumbel_FitComplete(x, totalN[exp], &est_mu, &est_lambda);
	  
	    printf("complete %6d %6d %9.5f %9.5f %8.6f %8.6f\n", 
		   totalN[exp], totalN[exp], mu, est_mu, lambda, est_lambda);
	  }
	printf("\n");
      }
  }

  /* Fitting to simulated censored datasets
   */
  if (do_censored) {
    for (exp = 0; exp < nexps; exp++)
      {
	for (trial = 0; trial < ntrials; trial++)
	  {
	    for (n = 0, i = 0; i < totalN[exp]; i++)
	      {
		val = esl_gumbel_Sample(r, mu, lambda);
		if (val >= phi) x[n++] = val;
	      }
	    esl_gumbel_FitCensored(x, n, totalN[exp]-n, phi, &est_mu, &est_lambda);
	  
	    printf("censored %6d %6d %9.5f %9.5f %8.6f %8.6f\n", 
		   totalN[exp], n, mu, est_mu, lambda, est_lambda);
	  }
	printf("\n");
      }
  }

  /* Fitting to simulated truncated datasets
   */
#ifdef eslAUGMENT_MINIMIZER
  if (do_truncated) {
    for (exp = 0; exp < nexps; exp++)
      {
	for (trial = 0; trial < ntrials; trial++)
	  {
	    for (n = 0, i = 0; i < totalN[exp]; i++)
	      {
		val = esl_gumbel_Sample(r, mu, lambda);
		if (val >= phi) x[n++] = val;
	      }
	    esl_gumbel_FitTruncated(x, n, phi, &est_mu, &est_lambda);
	  
	    printf("truncated %6d %6d %9.5f %9.5f %8.6f %8.6f\n", 
		   totalN[exp], n, mu, est_mu, lambda, est_lambda);
	  }
	printf("\n");
      }
  }
#endif /*eslAUGMENT_MINIMIZER*/

  esl_randomness_Destroy(r);
  free(x);
  return 0;
}
#endif /*eslGUMBEL_STATS*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
