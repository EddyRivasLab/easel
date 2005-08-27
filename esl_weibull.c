/* esl_weibull.c
 * Statistical routines for Weibull distributions.
 * 
 * SRE, Tue Aug  9 10:35:06 2005 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <math.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_weibull.h>

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_MINIMIZER
#include <esl_minimizer.h>
#endif

/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 

/* Function:  esl_wei_pdf()
 * Incept:    SRE, Tue Aug  9 13:42:17 2005 [St. Louis]
 *
 * Purpose:   Calculates the probability density function for the 
 *            Weibull pdf $P(X=x)$, given quantile <x>, offset <mu>, 
 *            and parameters <lambda> and <tau>.
 */
double
esl_wei_pdf(double x, double mu, double lambda, double tau)
{
  double y    = lambda * (x-mu);
  double val;

  if (x < mu) return 0.;

  val = lambda * tau * 
    exp((tau-1)*log(y)) *
    exp(- exp(tau * log(y)));
  return val;
}

/* Function:  esl_wei_logpdf()
 * Incept:    SRE, Tue Aug  9 13:42:07 2005 [St. Louis]
 *
 * Purpose:   Calculates the log probability density function for the
 *            Weibull, $\log P(X=x)$, given quantile <x>,
 *            offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_wei_logpdf(double x, double mu, double lambda, double tau)
{
  double y = lambda * (x-mu);
  double val;

  if (x < mu) return -eslINFINITY;
  val = log(tau) + tau*log(lambda) + (tau-1)*log(x-mu) - exp(tau * log(y));
  return val;
}

/* Function:  esl_wei_cdf()
 * Incept:    SRE, Tue Aug  9 15:19:06 2005 [St. Louis]
 *
 * Purpose:   Calculates the cumulative distribution function for the
 *            Weibull, $P(X \leq x)$, given quantile <x>,
 *            offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_wei_cdf(double x, double mu, double lambda, double tau)
{
  double y   = lambda*(x-mu);
  double tly = tau * log(y);

  if      (x <= mu)                return 0.0;
  else if (fabs(tly) < eslSMALLX1) return exp(tly);
  else                             return 1 - exp(-exp(tly));
}

/* Function:  esl_wei_logcdf()
 * Incept:    SRE, Tue Aug  9 15:21:52 2005 [St. Louis]
 *
 * Purpose:   Calculates the log of the cumulative distribution function for a
 *            Weibull, $P(X \leq x)$, given quantile <x>,
 *            offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_wei_logcdf(double x, double mu, double lambda, double tau)
{
  double y   = lambda*(x-mu);
  double tly = tau * log(y);

  if      (fabs(tly) < eslSMALLX1)              return tly;
  else if (fabs(exp(-exp(tly))) < eslSMALLX1)   return -exp(-exp(tly)); 
  else                                          return log(1 - exp(-exp(tly)));
}


/* Function:  esl_wei_surv()
 * Incept:    SRE, Tue Aug  9 15:23:06 2005 [St. Louis]
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ (that is, 1-CDF,
 *            the right tail probability mass) for a Weibull
 *            distribution, given quantile <x>, offset <mu>, and parameters
 *            <lambda> and <tau>.
 */
double
esl_wei_surv(double x, double mu, double lambda, double tau)
{
  double y   = lambda*(x-mu);
  double tly = tau * log(y);

  return exp(-exp(tly));
}

/* Function:  esl_wei_logsurv()
 * Incept:    SRE, Tue Aug  9 15:33:53 2005 [St. Louis]
 *
 * Purpose:   Calculates the log survivor function, $\log P(X>x)$ (that is, 
 *            log(1-CDF), the right tail log probability mass) for a 
 *            Weibull distribution, given quantile <x>, offset <mu>,
 *            and parameters <lambda> and <tau>.
 */
double
esl_wei_logsurv(double x, double mu, double lambda, double tau)
{
  double y   = lambda*(x-mu);
  double tly = tau * log(y);

  return -exp(tly);
}

/* Function:  esl_wei_invcdf()
 * Incept:    SRE, Sun Aug 21 14:50:00 2005 [St. Louis]
 *
 * Purpose:   Calculates the inverse CDF for a Weibull distribution
 *            with parameters <mu>, <lambda>, and <tau>, returning
 *            the quantile <x> at which the CDF is <p>, for $0<p<1$.
 */
double
esl_wei_invcdf(double p, double mu, double lambda, double tau)
{
  return mu + 1/lambda * exp(1/tau * log(-log((1.-p))));
}
/*-------------------- end densities & distributions ------------------------*/




/****************************************************************************
 * Generic API routines: for general interface w/ histogram module
 ****************************************************************************/ 

/* Function:  esl_wei_generic_pdf()
 * Incept:    SRE, Thu Aug 25 08:04:48 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_wei_pdf()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_wei_generic_pdf(double x, void *params)
{
  double *p = (double *) params;
  return esl_wei_pdf(x, p[0], p[1], p[2]);
}

/* Function:  esl_wei_generic_cdf()
 * Incept:    SRE, Fri Aug 19 09:34:26 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_wei_cdf()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_wei_generic_cdf(double x, void *params)
{
  double *p = (double *) params;
  return esl_wei_cdf(x, p[0], p[1], p[2]);
}

/* Function:  esl_wei_generic_surv()
 * Incept:    SRE, Fri Aug 19 09:34:26 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_wei_surv()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_wei_generic_surv(double x, void *params)
{
  double *p = (double *) params;
  return esl_wei_surv(x, p[0], p[1], p[2]);
}

/* Function:  esl_wei_generic_invcdf()
 * Incept:    SRE, Sun Aug 21 14:51:33 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_wei_invcdf()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_wei_generic_invcdf(double p, void *params)
{
  double *v = (double *) params;
  return esl_wei_invcdf(p, v[0], v[1], v[2]);
}
/*------------------------ end generic API ---------------------------------*/



/****************************************************************************
 * Routines for dumping plots for files
 ****************************************************************************/ 

/* Function:  esl_wei_Plot()
 * Incept:    SRE, Fri Aug 19 09:38:02 2005 [St. Louis]
 *
 * Purpose:   Plot some Weibull function <func> (for instance, <esl_wei_pdf()>)
 *            for Weibull parameters <mu>, <lambda>, and <tau>, for a range of
 *            quantiles x from <xmin> to <xmax> in steps of <xstep>;
 *            output to an open stream <fp> in xmgrace XY input format.
 *
 * Returns:   <eslOK>.
 */
int
esl_wei_Plot(FILE *fp, double mu, double lambda, double tau,
	     double (*func)(double x, double mu, double lambda, double tau), 
	     double xmin, double xmax, double xstep)
{
  double x;
  for (x = xmin; x <= xmax; x += xstep)
    fprintf(fp, "%f\t%g\n", x, (*func)(x, mu, lambda, tau));
  fprintf(fp, "&\n");
  return eslOK;
}
/*-------------------- end plot dumping routines ---------------------------*/





/****************************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_wei_Sample()
 * Incept:    SRE, Tue Aug  9 13:42:28 2005 [St. Louis]
 *
 * Purpose:   Sample a Weibull random variate,
 *            by the transformation method.
 */
double
esl_wei_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau)
{
  double p;
  p = esl_rnd_UniformPositive(r); 
  return esl_wei_invcdf(p, mu, lambda, tau);
} 
#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/


/****************************************************************************
 * Maximum likelihood fitting
 ****************************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER
/* Easel's conjugate gradient descent code allows a single void ptr to
 * point to any necessary fixed data, so we put everything into one
 * structure:
 */
struct wei_data {
  double *x;	        /* data: n observed samples    */
  int     n;		/* number of observed samples  */
  double  mu;		/* mu is considered to be known, not fitted */
};

/* wei_func():
 * Returns the negative log likelihood of a complete data sample,
 * in the API of the conjugate gradient descent optimizer in esl_minimizer.
 */
static double
wei_func(double *p, int nparam, void *dptr)
{
  double lambda, tau;
  struct wei_data *data;
  double logL;
  int    i; 
    
  /* Unpack what the optimizer gave us.
   */
  lambda = exp(p[0]); /* see below for c.o.v. notes */
  tau    = exp(p[1]);
  data   = (struct wei_data *) dptr;

  logL = 0.;
  for (i = 0; i < data->n; i++)
    logL += esl_wei_logpdf(data->x[i], data->mu, lambda, tau);
  return -logL;			/* goal: minimize NLL */
}

/* Function:  esl_wei_FitComplete()
 * Incept:    SRE, Tue Aug  9 13:55:37 2005 [St. Louis]
 *
 * Purpose:   Given an array of <n> samples <x[0]..x[n-1>, fit
 *            them to a stretched exponential distribution starting
 *            at lower bound <mu> (all $x_i > \mu$), and 
 *            return maximum likelihood parameters <ret_lambda>
 *            and <ret_tau>.
 *            
 * Args:      x          - complete GEV-distributed data [0..n-1]
 *            n          - number of samples in <x>
 *            mu         - lower bound of the distribution (all x_i > mu)
 *            ret_lambda - RETURN: maximum likelihood estimate of lambda
 *            ret_tau    - RETURN: maximum likelihood estimate of tau
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECONVERGENCE> if the fit doesn't converge.
 *
 * Xref:      STL9/136-137
 */
int
esl_wei_FitComplete(double *x, int n, double mu,
		    double *ret_lambda, double *ret_tau)
{
  struct wei_data data;
  double p[2];			/* parameter vector                  */
  double u[2];			/* max initial step size vector      */
  double wrk[8];		/* 4 tmp vectors of length 2         */
  double mean, variance;
  double lambda, tau;      	/* initial param guesses             */
  double tol = 1e-6;		/* convergence criterion for CG      */
  double fx;			/* f(x) at minimum; currently unused */
  int    status;

  data.x   = x;
  data.n   = n;
  data.mu  = mu;

  /* Make a good initial guess, based on exponential fit.
   */
  esl_stats_Mean(x, n, &mean, &variance);
  lambda = 1 / (mean - mu);
  tau    = 0.9;

  /* Change of variables;
   *   lambda > 0, so c.o.v.  lambda = exp^w,  w = log(lambda);
   *   tau > 0, same c.o.v.
   */
  p[0] = log(lambda);		
  p[1] = log(tau);

  u[0] = 1.0;
  u[1] = 1.0;

  /* pass problem to the optimizer
   */
  status = esl_min_ConjugateGradientDescent(p, u, 2, 
					    &wei_func, NULL,
					    (void *)(&data),
					    tol, wrk, &fx);
  *ret_lambda = exp(p[0]);
  *ret_tau    = exp(p[1]);
  return status;
}

#ifdef eslAUGMENT_HISTOGRAM
struct wei_binned_data {
  ESL_HISTOGRAM *h;	/* contains the binned observed data        */
  double  mu;		/* mu is considered to be known, not fitted */
};

/* wei_binned_func():
 * Returns the negative log likelihood of a binned data sample,
 * in the API of the conjugate gradient descent optimizer in esl_minimizer.
 */
static double
wei_binned_func(double *p, int nparam, void *dptr)
{
  struct wei_binned_data *data = (struct wei_binned_data *) dptr;
  ESL_HISTOGRAM          *h    = data->h;
  double lambda, tau;
  double logL;
  double ai,bi;
  int    i; 
  double tmp;
    
  /* Unpack what the optimizer gave us.
   */
  lambda = exp(p[0]); /* see below for c.o.v. notes */
  tau    = exp(p[1]);

  logL = 0.;
  for (i = h->cmin; i <= h->imax; i++)
    {
      if (h->obs[i] == 0) continue;

      ai = esl_histogram_Bin2LBound(h,i);
      bi = esl_histogram_Bin2UBound(h,i);
      if (ai < data->mu) ai = data->mu;

      tmp = esl_wei_cdf(bi, data->mu, lambda, tau) -
            esl_wei_cdf(ai, data->mu, lambda, tau);

      if (tmp == 0.) return eslINFINITY;
      ESL_DASSERT1( (tmp > 0.)); 

      logL += h->obs[i] * log(tmp);
    }
  return -logL;			/* goal: minimize NLL */
}

/* Function:  esl_wei_FitCompleteBinned()
 * Incept:    SRE, Sun Aug 21 15:17:45 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <g> with binned observations, where each
 *            bin i holds some number of observed samples x with values from 
 *            lower bound l to upper bound u (that is, $l < x \leq u$), and
 *            <mu>, the known offset (minimum value) of the distribution;
 *            return maximum likelihood parameters <ret_lambda>
 *            and <ret_tau>.
 *            
 * Args:      x          - complete GEV-distributed data [0..n-1]
 *            n          - number of samples in <x>
 *            ret_mu     - lower bound of the distribution (all x_i > mu)
 *            ret_lambda - RETURN: maximum likelihood estimate of lambda
 *            ret_tau    - RETURN: maximum likelihood estimate of tau
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECONVERGENCE> if the fit doesn't converge.
 *
 * Xref:      STL9/136-137
 */
int
esl_wei_FitCompleteBinned(ESL_HISTOGRAM *h, double *ret_mu,
			  double *ret_lambda, double *ret_tau)
{
  struct wei_binned_data data;
  double p[2];			/* parameter vector                  */
  double u[2];			/* max initial step size vector      */
  double wrk[8];		/* 4 tmp vectors of length 2         */
  double mean;
  double mu, lambda, tau;      	/* initial param guesses             */
  double tol = 1e-6;		/* convergence criterion for CG      */
  double fx;			/* f(x) at minimum; currently unused */
  int    status;
  int    i;
  double ai;

  /* Make a good initial guess, based on exponential fit.
   */
  mean = 0.;
  for (i = h->cmin; i <= h->imax; i++) 
    { 
      ai = esl_histogram_Bin2LBound(h, i);
      ai += 0.5*h->w;		/* midpoint in bin */
      mean += (double)h->obs[i] * ai;
    }
  mean  /= h->No;
  mu     = esl_histogram_Bin2LBound(h,h->cmin);
  lambda = 1 / (mean - mu);
  tau    = 0.9;

  /* load the data structure */
  data.h   = h;
  data.mu  = mu;

  /* Change of variables;
   *   lambda > 0, so c.o.v.  lambda = exp^w,  w = log(lambda);
   *   tau > 0, same c.o.v.
   */
  p[0] = log(lambda);		
  p[1] = log(tau);

  u[0] = 1.0;
  u[1] = 1.0;

  /* pass problem to the optimizer
   */
  status = esl_min_ConjugateGradientDescent(p, u, 2, 
					    &wei_binned_func, NULL,
					    (void *)(&data),
					    tol, wrk, &fx);
  *ret_mu     = mu;
  *ret_lambda = exp(p[0]);
  *ret_tau    = exp(p[1]);
  return status;
}
#endif /*eslAUGMENT_HISTOGRAM*/
#endif /*eslAUGMENT_MINIMIZER*/
/*--------------------------- end fitting ----------------------------------*/


/****************************************************************************
 * Example, test, and stats drivers
 ****************************************************************************/ 
/* Example main()
 */
#ifdef eslWEI_EXAMPLE
/*::cexcerpt::wei_example::begin::*/
/* compile: 
     gcc -g -Wall -I. -o example -DeslWEI_EXAMPLE -DeslAUGMENT_RANDOM\
       -DeslAUGMENT_MINIMIZER esl_weibull.c esl_random.c esl_minimizer.c\
       esl_vectorops.c easel.c -lm
   or:
     gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o example -DeslWEI_EXAMPLE\
       esl_weibull.c -leasel -lm
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_weibull.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r = esl_randomness_Create(42);
  double  mu        = -2.1;         
  double  lambda    =  1.0;         
  double  tau       =  0.8;	   
  double  min       = -2.0;
  double  max       = 100.0;
  double  step      = 0.1;
  int     n         = 100000; 
  double *x;
  double  pdf;
  int     i;
  double  est_lambda, est_tau;

  /*
  for (x = min; x <= max; x += step)
    { 
      pdf = esl_wei_pdf(x, mu, lambda, tau);
      printf("%.4f  %g\n", x, pdf);
    }
  printf("&\n");
  */
  
  x = malloc(sizeof(double) * n);

  pdf  = 0.;
  for (i = 0; i < n; i++)
    {
      x[i]   = esl_wei_Sample(r, mu, lambda, tau);
      pdf   += log(esl_wei_pdf(x[i], mu, lambda, tau));
    }
  printf("At parametric: -logL=%f\n", -pdf);
  
  esl_wei_FitComplete(x, n, mu, &est_lambda, &est_tau);

  printf("estimated lambda = %g\n", est_lambda);
  printf("estimated tau    = %g\n", est_tau);

  esl_randomness_Destroy(r);
  free(x);
  return 0;
}
/*::cexcerpt::wei_example::end::*/
#endif /*eslWEI_EXAMPLE*/


#ifdef eslWEI_STATS
/* compile: 
     gcc -g -Wall -I. -o stats -DeslWEI_STATS -DeslAUGMENT_RANDOM\
       -DeslAUGMENT_MINIMIZER esl_stretchexp.c esl_random.c esl_minimizer.c\
       esl_vectorops.c easel.c -lm
  or:
     gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o stats -DeslWEI_STATS\
       esl_stretchexp.c -leasel -lm
 * run:     ./stats <test#>...
 * e.g. 
 *          ./stats 1 2 3
 * would run tests 1, 2, 3.
 */
#include <stdio.h>
#include <math.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_stretchexp.h>

#define MAX_STATS_TESTS 10
static void test_range(FILE *fp, double mu, double lambda, double tau);

int
main(int argc, char **argv)
{
  FILE *fp;
  double  mu        = 0.0;
  double  lambda    = 1.0;  
  double  tau       = 0.5;      
  double  xmin      = 0.1;
  double  xmax      = 40.;
  double  xstep     = .1; 
  double  x;
  int     do_test[MAX_STATS_TESTS+1];
  int     i;

  if (argc == 1) {
    printf("Diagnostic test output driver for stretchexp module.\n");
    printf("Usage: ./stats <#> [<#>...]\n");
    printf("Available test numbers:\n");
    printf("#     Description        Output format   Output file\n");
    printf("--  ------------------   -------------   -----------\n");
    printf("1    pdf plot            xmgrace xy       stats.1   \n");
    printf("2    log pdf plot        xmgrace xy       stats.2   \n");
    printf("3    cdf plot            xmgrace xy       stats.3   \n");
    printf("4    log cdf plot        xmgrace xy       stats.4   \n");
    printf("5    survivor plot       xmgrace xy       stats.5   \n");
    printf("6    log surv plot       xmgrace xy       stats.6   \n");
    printf("7    range tests         R table          stats.7   \n");
    printf("----------------------------------------------------\n");
    printf("Using mu = %f, lambda = %f, tau = %f\n", mu, lambda, tau);
    return 0;
  }

  for (i = 0; i <= MAX_STATS_TESTS; i++) do_test[i] = 0;
  for (i = 1; i < argc; i++)
    do_test[atoi(argv[i])] = 1;

  /* stats.1: density plot, xmgrace xy file */
  if (do_test[1]) {
    if ((fp = fopen("stats.1", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_wei_pdf(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.2: log density plot, xmgrace xy file */
  if (do_test[2]) {
    if ((fp = fopen("stats.2", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %gf\n", x, esl_wei_logpdf(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.3: CDF plot, xmgrace xy file */
  if (do_test[3]) {
    if ((fp = fopen("stats.3", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_wei_cdf(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.4: log CDF plot, xmgrace xy file */
  if (do_test[4]) {
    if ((fp = fopen("stats.4", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep) 
      fprintf(fp, "%.1f  %g\n", x, esl_wei_logcdf(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }
  
  /* stats.5: survivor plot (right tail), xmgrace xy file */
  if (do_test[5]) {
    if ((fp = fopen("stats.5", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_wei_surv(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }
    
  /* stats.6: log survivor plot, xmgrace xy file */
  if (do_test[6]) {
    if ((fp = fopen("stats.6", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep) 
      fprintf(fp, "%.1f  %g\n", x, esl_wei_logsurv(x, mu, lambda, tau));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.7: test extreme range of x
   */
  if (do_test[7]) {
    if ((fp = fopen("stats.7", "w")) == NULL) abort();
    test_range(fp, mu, lambda, tau);
    fclose(fp);
  }
  return 0;
}

static void
test_range(FILE *fp, double mu, double lambda, double tau)
{
  double xpoints[] = { 0.,     1e-100, 1e-10,  1.0,   
                       10.,    100,     200,   300,
                       400,    500,     1000,   1e4,     
                       1e100,  1e300};
  double n = sizeof(xpoints)/sizeof(double);
  int    i;
  double x;

  fprintf(fp, "%14s %14s %14s %14s %14s %14s %14s\n",
	  "", "pdf", "logpdf", "cdf", "logcdf", "surv", "logsurv");
  for (i = 0; i < n; i++)
    {
      x = xpoints[i];
      fprintf(fp, "%14g ", x);
      fprintf(fp, "%14g ", esl_wei_pdf    (x, mu, lambda, tau));
      fprintf(fp, "%14g ", esl_wei_logpdf (x, mu, lambda, tau));
      fprintf(fp, "%14g ", esl_wei_cdf    (x, mu, lambda, tau));
      fprintf(fp, "%14g ", esl_wei_logcdf (x, mu, lambda, tau));
      fprintf(fp, "%14g ", esl_wei_surv   (x, mu, lambda, tau));
      fprintf(fp, "%14g ", esl_wei_logsurv(x, mu, lambda, tau));
      fprintf(fp, "\n");
    }
}

#endif /*eslWEI_STATS*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
