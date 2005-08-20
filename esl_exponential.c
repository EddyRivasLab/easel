/* esl_exponential.c
 * Statistical routines for exponential distributions.
 * 
 * SRE, Wed Aug 10 08:15:57 2005   xref:STL9/138  [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_vectorops.h>
#include <esl_exponential.h>
#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 
/* gotchas: 
 *   - any lambda > 0 is valid... including infinity! (fitting code
 *     may try to test such lambdas, and it must get back valid numbers,
 *     never an NaN, or it will fail.) Note that IEEE754 allows us
 *     to calculate log(inf) = inf, exp(-inf) = 0, and exp(inf) = inf.
 *     But remember that inf-inf = NaN, so Don't Do That.
 */

/* Function:  esl_exp_pdf()
 * Incept:    SRE, Wed Aug 10 08:30:46 2005 [St. Louis]
 *
 * Purpose:   Calculates the probability density function for the
 *            exponential, $P(X=x)$, given quantile <x>, offset <mu>,
 *            and decay parameter <lambda>.
 */
double
esl_exp_pdf(double x, double mu, double lambda)
{
  if (x < mu) return 0.;
  return (lambda * exp(-lambda*(x-mu)));
}

/* Function:  esl_exp_logpdf()
 * Incept:    SRE, Wed Aug 10 08:35:06 2005 [St. Louis]
 *
 * Purpose:   Calculates the log probability density function for the
 *            exponential, $P(X=x)$, given quantile <x>, offset <mu>,
 *            and decay parameter <lambda>.
 */
double
esl_exp_logpdf(double x, double mu, double lambda)
{
  if (x < mu) return -eslINFINITY;

  if (lambda == eslINFINITY) 
    {	/* limit as lambda->inf: avoid inf-inf! */
      if (x == mu) return  eslINFINITY;
      else         return -eslINFINITY;
    }
  return (log(lambda) - lambda*(x-mu));
}

/* Function:  esl_exp_cdf()
 * Incept:    SRE, Wed Aug 10 08:36:04 2005 [St. Louis]
 *
 * Purpose:   Calculates the cumulative distribution function for the
 *            exponential, $P(X \leq x)$, given quantile <x>, offset <mu>,
 *            and decay parameter <lambda>.
 */
double
esl_exp_cdf(double x, double mu, double lambda)
{
  double y = lambda*(x-mu);	/* y>=0 because lambda>0 and x>=mu */

  if (x < mu) return 0.;

  /* 1-e^-y ~ y for small |y| */
  if (y < eslSMALLX1) return y;
  else                return 1 - exp(-y);
}

/* Function:  esl_exp_logcdf()
 * Incept:    SRE, Wed Aug 10 10:03:28 2005 [St. Louis]
 *
 * Purpose:   Calculates the log of the cumulative distribution function
 *            for the exponential, $log P(X \leq x)$, given quantile <x>,
 *            offset <mu>, and decay parameter <lambda>.
 */
double
esl_exp_logcdf(double x, double mu, double lambda)
{
  double y  = lambda * (x-mu);
  double ey = exp(-y);

  if (x < mu) return -eslINFINITY;

  /* When y is small, 1-e^-y = y, so answer is log(y);
   * when y is large, exp(-y) is small, log(1-exp(-y)) = -exp(-y).
   */
  if      (y == 0)           return -eslINFINITY; /* don't allow NaN */
  else if (y  < eslSMALLX1)  return log(y);
  else if (ey < eslSMALLX1)  return -ey;
  else                       return log(1-ey);
}

/* Function:  esl_exp_surv()
 * Incept:    SRE, Wed Aug 10 10:14:49 2005 [St. Louis]
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ (that is, 1-CDF,
 *            the right tail probability mass) for an exponential distribution,
 *            given quantile <x>, offset <mu>, and decay parameter <lambda>.
 */
double
esl_exp_surv(double x, double mu, double lambda)
{
  if (x < mu) return 1.0;
  return exp(-lambda * (x-mu));
}

/* Function:  esl_exp_logsurv()
 * Incept:    SRE, Wed Aug 10 10:14:49 2005 [St. Louis]
 *
 * Purpose:   Calculates the log survivor function, $\log P(X>x)$ (that is,
 *            log(1-CDF), the log of the right tail probability mass) for an 
 *            exponential distribution, given quantile <x>, offset <mu>, and 
 *            decay parameter <lambda>.
 */
double
esl_exp_logsurv(double x, double mu, double lambda)
{
  if (x < mu) return 0.0;
  return -lambda * (x-mu);
}


/****************************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_exp_Sample()
 * Incept:    SRE, Wed Aug 10 10:46:51 2005 [St. Louis]
 *
 * Purpose:   Sample an exponential random variate
 *            by the transformation method, given offset <mu>
 *            and decay parameter <lambda>.
 */
double
esl_exp_Sample(ESL_RANDOMNESS *r, double mu, double lambda)
{
  double p, x;
  p = esl_rnd_UniformPositive(r); 

  x = mu - 1/lambda * log(p);	/* really log(1-p), but if p is uniform on 0..1 
				 * then so is 1-p. 
                                 */
  return x;
} 
#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/




/****************************************************************************
 * Maximum likelihood fitting
 ****************************************************************************/ 

/* Function:  esl_exp_FitComplete()
 * Incept:    SRE, Wed Aug 10 10:53:47 2005 [St. Louis]
 *
 * Purpose:   Given an array of <n> samples <x[0]..x[n-1]>, fit
 *            them to an exponential distribution starting at a
 *            known lower bound <mu> (all $x_i \geq \mu$). 
 *            Return maximum likelihood decay parameter <ret_lambda>.
 *
 * Args:      x          - complete exponentially-distributed data [0..n-1]
 *            n          - number of samples in <x>
 *            mu         - lower bound of the distribution (all x_i >= mu)
 *            ret_lambda - RETURN: maximum likelihood estimate of lambda
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      STL9/138.
 */
int
esl_exp_FitComplete(double *x, int n, double mu, double *ret_lambda)
{
  double mean;
  int    i;

  mean = 0.;
  for (i = 0; i < n; i++) mean += x[i] - mu;
  mean /= (double) n;

  *ret_lambda = 1./mean;	/* ML estimation is trivial in this case */
  return eslOK;
}


/****************************************************************************
 * "stats" code driver for dumping plots and tables for verification
 ****************************************************************************/ 
#ifdef eslEXP_STATS
/* compile: 
     gcc -g -Wall -I. -o stats -DeslEXP_STATS -DeslAUGMENT_RANDOM\
       -DeslAUGMENT_MINIMIZER esl_exponential.c esl_random.c esl_minimizer.c\
       esl_vectorops.c easel.c -lm
  or:
     gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o stats -DeslEXP_STATS\
       esl_exponential.c -leasel -lm
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
#include <esl_exponential.h>

#define MAX_STATS_TESTS 10
static void test_range(FILE *fp, double mu, double lambda);

int
main(int argc, char **argv)
{
  FILE *fp;
  double  mu        = 0.0;
  double  lambda    = 1.0;  
  double  xmin      = 0.;
  double  xmax      = 40.;
  double  xstep     = .1; 
  double  x;
  int     do_test[MAX_STATS_TESTS+1];
  int     i;

  if (argc == 1) {
    printf("Diagnostic test output driver for exponential module.\n");
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
    printf("Using mu = %f, lambda = %f\n", mu, lambda);
    return 0;
  }

  
 for (i = 0; i <= MAX_STATS_TESTS; i++) do_test[i] = 0;
  for (i = 1; i < argc; i++)
    do_test[atoi(argv[i])] = 1;

  /* stats.1: density plot, xmgrace xy file */
  if (do_test[1]) {
    if ((fp = fopen("stats.1", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_exp_pdf(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.2: log density plot, xmgrace xy file */
  if (do_test[2]) {
    if ((fp = fopen("stats.2", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %gf\n", x, esl_exp_logpdf(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.3: CDF plot, xmgrace xy file */
  if (do_test[3]) {
    if ((fp = fopen("stats.3", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_exp_cdf(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.4: log CDF plot, xmgrace xy file */
  if (do_test[4]) {
    if ((fp = fopen("stats.4", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep) 
      fprintf(fp, "%.1f  %g\n", x, esl_exp_logcdf(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }
  
  /* stats.5: survivor plot (right tail), xmgrace xy file */
  if (do_test[5]) {
    if ((fp = fopen("stats.5", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep)
      fprintf(fp, "%.1f  %g\n", x, esl_exp_surv(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }
    
  /* stats.6: log survivor plot, xmgrace xy file */
  if (do_test[6]) {
    if ((fp = fopen("stats.6", "w")) == NULL) abort();
    for (x = xmin; x <= xmax; x+= xstep) 
      fprintf(fp, "%.1f  %g\n", x, esl_exp_logsurv(x, mu, lambda));
    fprintf(fp, "&\n");
    fclose(fp);
  }

  /* stats.7: test extreme range of x
   */
  if (do_test[7]) {
    if ((fp = fopen("stats.7", "w")) == NULL) abort();
    test_range(fp, mu, lambda);
    fclose(fp);
  }
  return 0;
}

static void
test_range(FILE *fp, double mu, double lambda)
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
      fprintf(fp, "%14g ", esl_exp_pdf    (x, mu, lambda));
      fprintf(fp, "%14g ", esl_exp_logpdf (x, mu, lambda));
      fprintf(fp, "%14g ", esl_exp_cdf    (x, mu, lambda));
      fprintf(fp, "%14g ", esl_exp_logcdf (x, mu, lambda));
      fprintf(fp, "%14g ", esl_exp_surv   (x, mu, lambda));
      fprintf(fp, "%14g ", esl_exp_logsurv(x, mu, lambda));
      fprintf(fp, "\n");
    }
}
#endif /*eslEXP_STATS*/






/*****************************************************************
 * @LICENSE@
 *****************************************************************/
