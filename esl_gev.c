/* esl_gev.c
 * Statistical routines for generalized extreme value (GEV) distributions.
 *
 * GEV distribution 
 *     G(x) = exp{ -[1 + \alpha \lambda(x - \mu)]^{-1/\alpha} }
 * where:
 *    \mu     = location parameter
 *    \lambda = scale parameter   (\lambda = 1/\sigma, in [Coles01] notation)
 *    \alpha  = shape parameter   (\alpha  = \xi, in [Coles01] notation) 
 * 
 * lim_{\alpha -> 0} is a type I EVD (Gumbel)
 * \alpha > 0  is a Type II  EVD (Frechet)
 * \alpha < 0  is a Type III EVD (Weibull)
 * 
 * Reference: 
 *   [Coles01] 
 *   S. Coles, An Introduction to Statistical Modeling of Extreme Values, 
 *   Springer, 2001.
 *            
 * Xref: 
 *   STL9/118, 2005/0712-easel-gev-impl. Verified against evd package in R.
 *
 * SRE, Tue Jul 12 09:02:08 2005
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
#include <esl_gev.h>

extern void
ConjugateGradientDescent(double *p, int n, double ftol, int *ret_iter,
			 double *ret_min, double (*func)(double *),
			 void (*dfunc)(double *, double *),
			 void (*pfunc)(double *, double, double *));


/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 

/* Function:  esl_gev_pdf()
 * Incept:    SRE, Tue Jul 12 10:53:19 2005 [St. Louis]
 *
 * Purpose:   Calculates the probability density function for the 
 *            generalized extreme value distribution, $P(X=x)$, given
 *            quantile <x> and GEV location, scale, shape parameters 
 *            <mu>, <lambda>, <alpha>.
 */
double
esl_gev_pdf(double x, double mu, double lambda, double alpha)
{
  double y     = lambda * (x-mu);
  double ya1   = 1. + alpha * y;
  double lya1;  

  /* Special case: if alpha is tiny, approximate by a Gumbel */
  if (fabs(y*alpha) < 1e-12) return (lambda * exp(-y - exp(-y)));

  /* Else, use GEV; but use log/exp to avoid a pow() call,
   * as that's almost 2x faster (on my machine anyway).
   */
  if (ya1 <= 0) return 0.;
  lya1 = log(ya1);
  return (lambda * exp(-(1.+ 1./alpha)*lya1 - exp(-lya1/alpha)));
}

/* Function:  esl_gev_logpdf()
 * Incept:    SRE, Tue Jul 12 16:43:13 2005 [St. Louis]
 *
 * Purpose:   Calculates the log probability density function for the
 *            generalized extreme value distribution, $\log P(X=x)$, 
 *            given quantile <x> and GEV location, scale, shape
 *            parameters <mu>, <lambda>, <alpha>.
 */
double
esl_gev_logpdf(double x, double mu, double lambda, double alpha)
{
  double y          = lambda *(x-mu);
  double ya1        = 1. + alpha*y;
  double lya1;

  /* Special case: if alpha is tiny, approx by a Gumbel */
  if (fabs(y*alpha) < 1e-12) return (log(lambda) -y - exp(-y));

  /* It's important not to return NaN for this domain error;
   * minimizer relies on being able to compare logL's for any parameter,
   * and you can't compare NaN to anything.
   */
  if (ya1 <= 0) return -eslINFINITY;

  lya1 = log(ya1);
  return (log(lambda) - (1.+1./alpha)*lya1 - exp(-lya1/alpha));
}


/* Function:  esl_gev_cdf()
 * Incept:    SRE, Tue Jul 12 16:45:55 2005 [St. Louis]
 *
 * Purpose:   Calculates the cumulative distribution function for the
 *            generalized extreme value distribution, $P(X \leq x)$, 
 *            given quantile <x> and GEV location, scale, shape
 *            parameters <mu>, <lambda>, <alpha>.
 */
double
esl_gev_cdf(double x, double mu, double lambda, double alpha)
{
  double y          = lambda *(x-mu);
  double ya1        = 1. + alpha*y;
  double lya1;

  /* Special case: if alpha is tiny, approx by a Gumbel */
  if (fabs(y*alpha) < 1e-12) return (exp(-exp(-y)));

  if (ya1 <= 0) {
    if (x < mu) return 0.0; /* the frechet case */
    else        return 1.0; /* the weibull case */
  }
  lya1 = log(ya1);
  return (exp(-exp(-lya1/alpha)));
}



/* Function:  esl_gev_logcdf()
 * Incept:    SRE, Tue Jul 12 17:15:49 2005 [St. Louis]
 *
 * Purpose:   Calculates the log of the cumulative distribution function for the
 *            generalized extreme value distribution, $\log P(X \leq x)$, 
 *            given quantile <x> and GEV location, scale, shape
 *            parameters <mu>, <lambda>, <alpha>.
 */
double
esl_gev_logcdf(double x, double mu, double lambda, double alpha)
{
  double y          = lambda *(x-mu);
  double ya1        = 1. + alpha*y;
  double lya1;

  /* Special case: if alpha is tiny, approx by a Gumbel */
  if (fabs(y*alpha) < 1e-12) return (-exp(-y));

  if (ya1 <= 0) {
    if (x < mu) return -eslINFINITY;    /* Frechet  */
    else        return 0.0;     	/* Weibull  */
  }

  lya1 = log(ya1);
  return (-exp(-lya1/alpha));
}


/* Function:  esl_gev_surv()
 * Incept:    SRE, Wed Jul 13 07:41:12 2005 [St. Louis]
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ (that is, 1-cdf),
 *            the right tail's probability mass,  given quantile <x> and
 *            GEV location, scale, shape parameters <mu>, <lambda>, <alpha>.
 */
double
esl_gev_surv(double x, double mu, double lambda, double alpha)
{
   double y          = lambda *(x-mu);
   double ya1        = 1. + alpha*y;
   double lya1;

   /* Special case: for tiny alpha, use Gumbel (xref esl_gumbel.c) */
   if (fabs(y*alpha) < 1e-12) 
     return ((y > -0.5*log(DBL_EPSILON)) ? exp(-y) : (1 - exp(-exp(-y))));
   
   if (ya1 <= 0) {
     if (x < mu) return 1.0;	/* the frechet case */
     else        return 0.0;	/* the weibull case */
   }
   lya1 = log(ya1)/alpha;
   return ((lya1 > -0.5*log(DBL_EPSILON)) ? exp(-lya1) : (1 - exp(-exp(-lya1))));
}


/* Function:  esl_gev_logsurv()
 * Incept:    SRE, Wed Jul 13 08:14:48 2005 [St. Louis]
 *
 * Purpose:   Calculates the log survivor function $\log P(X>x)$ for a 
 *            generalized extreme value distribution (that is, 
 *            $\log (1 - \mbox(cdf))$); log of the right tail's probability
 *            mass; given quantile <x> and GEV location, scale, shape
 *            parameters <mu>, <lambda>, <alpha>.
 */
double
esl_gev_logsurv(double x, double mu, double lambda, double alpha)
{
   double y          = lambda *(x-mu);
   double ya1        = 1. + alpha*y;
   double lya1;

   /* Special case: for tiny alpha, use Gumbel (xref esl_gumbel.c) */
   if (fabs(y*alpha) < 1e-12) 
     {
       if      (y > -0.5 * log(DBL_EPSILON)) return (-y);
       else if (y < -2.9)                    return (-exp(-exp(-y)));
       else                                  return (log(1-exp(-exp(-y))));
     }
   
   /* See esl_gumbel.c for analysis of the crossovers in
    * the three cases (small, large, and ok lya1)
    */
   if (ya1 <= 0) {
     if (x < mu) return 1.0;        	/* Frechet case */
     else        return -eslINFINITY;   /* Weibull case */
   }

   lya1 = log(ya1)/alpha;
   if      (lya1 > -0.5 * log(DBL_EPSILON)) return (-lya1);
   else if (lya1 < -2.9)                    return (-exp(-exp(-lya1)));
   else                                     return (log(1-exp(-exp(-lya1))));
}
/*-------------------- end densities & distributions ---------------------------*/


/****************************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM
/* Function:  esl_gev_Sample()
 * Incept:    SRE, Wed Jul 13 08:30:49 2005 [St. Louis]
 *
 * Purpose:   Sample a GEV-distributed random variate
 *            by the transformation method.
 */
double
esl_gev_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double alpha)
{
  double p, x;
  p = esl_rnd_UniformPositive(r); 

  /* failover to Gumbel sample, for tiny alpha */
  if (fabs(alpha) < 1e-12) return (mu - log(-1. * log(p)) / lambda);

  x = mu + (exp(-alpha*log(-log(p))) - 1.) / (alpha * lambda) ;
  return x;
} 
#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/




/****************************************************************************
 * Maximum likelihood fitting to GEV distributions
 ****************************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER

/* Easel's conjugate gradient descent code allows a single void ptr to
 * point to any necessary fixed data, so we'll put everything into one
 * structure:
 */
struct gev_data {
  double *x;	/* data: n observed samples from a Gumbel */
  int     n;	/* number of observed samples */
  double  phi;	/* censoring or truncation threshold: all observed x_i >= phi */
  int     z;	/* # of censored samples */
};

/* gev_complete_func():
 * Returns the negative log likelihood of a complete GEV data sample;
 * in the API of the conjugate gradient descent optimizer in esl_minimizer.
 */
static double
gev_complete_func(double *p, int nparam, void *dptr)
{
  double mu, w, lambda, alpha;
  struct gev_data *data;
  double logL;
  int    i; 
    
  /* Unpack what the optimizer gave us.
   */
  mu     = p[0];
  w      = p[1];   /* w is a c.o.v. to allow unconstrained opt of lambda>0 */
  lambda = exp(w);
  alpha  = p[2];
  data   = (struct gev_data *) dptr;

  logL = 0.;
  for (i = 0; i < data->n; i++)
    logL += esl_gev_logpdf(data->x[i], mu, lambda, alpha);
  return -logL;			/* goal: minimize NLL */
}


static void
gev_numeric_grad(double *p, int nparam, void *dptr, double *dp)
{
  double mu, w, lambda, alpha;
  struct gev_data *data;
  double *x;
  double dmu, dw, dalpha;
  double fx1,fx2;
  double delta;
    
  /* Unpack what the optimizer gave us */
  mu     = p[0];
  w      = p[1];   /* w is a c.o.v. to allow unconstrained opt of lambda>0 */
  lambda = exp(w);
  alpha  = p[2];
  data   = (struct gev_data *) dptr;
  x      = data->x;

  delta = 0.001;
  fx1 =  gev_complete_func(p, 3, dptr);

  p[0]  = mu + delta*mu;
  fx2   = gev_complete_func(p, 3, dptr);
  dmu   = (fx2-fx1)/(delta*mu);
  p[0]  = mu;

  p[1]  = w + delta*w;
  fx2   = gev_complete_func(p, 3, dptr);
  dw    = (fx2-fx1)/(delta*w);
  p[1]  = w;

  p[2]  = alpha + delta*alpha;
  fx2   = gev_complete_func(p, 3, dptr);
  dalpha= (fx2-fx1)/(delta*alpha);
  p[2]  = alpha;

  dp[0] = dmu;
  dp[1] = dw;
  dp[2] = dalpha;
}


/* gev_complete_grad()
 * Computes the gradient of the negative log likelihood of a complete
 * GEV sample; in the API of the CG optimizer.
 */
static void
gev_complete_grad(double *p, int nparam, void *dptr, double *dp)
{
  double mu, w, lambda, alpha;
  struct gev_data *data;
  double *x;
  int    i; 
  double dmu, dw, dalpha;
  double y, ay, ay1, lay1;
    
  /* Unpack what the optimizer gave us */
  mu     = p[0];
  w      = p[1];   /* w is a c.o.v. to allow unconstrained opt of lambda>0 */
  lambda = exp(w);
  alpha  = p[2];
  data   = (struct gev_data *) dptr;
  x      = data->x;

  dmu    = 0.;
  dw     = data->n; /* d/dw, term1 */
  dalpha = 0.;

  for (i = 0; i < data->n; i++)
    {
      y    = lambda * (x[i]-mu);
      ay   = alpha*y;
      ay1  = 1+ay;		/* 1+ay=1, for ay < DBL_EPSILON */
      lay1 = log(ay1);
      
      /* d/dmu, term1. (will become 1, for small alpha.) */
      dmu += (alpha+1) / ay1;
      
      /* d/dmu, term2. For tiny ay, use log(1+x) ~ x to simplify. */
      if (fabs(ay) < 1e-12) dmu -= exp(-y);
      else                  dmu -= exp(-(1+1/alpha) * lay1);

      /* d/dw, term2. converges to -y, for small alpha. */
      dw -= y*(1+alpha) / ay1;

      /* d/dw, term2. For tiny ay, use log(1+x) ~ x to simplify. */
      if (fabs(ay) < 1e-12) dw += y*exp(-y);
      else                  dw += y*exp(-(1+1/alpha) * lay1);

      /* d/dalpha, term1
       */
      dalpha -= (1 + 1/alpha) * y/ay1;

      /* d/dalpha, terms 2,3,4: for tiny ay, simplify.
       * d/dalpha will go to +/-inf for alpha ~ 0, so watch out.
       */
      if (fabs(ay) < 1e-12) {
	dalpha += y/alpha;
	dalpha += y*exp(-y) / (alpha*ay1);
	dalpha -= y*exp(-y) / alpha;
      } else {
	dalpha += lay1 / (alpha*alpha);
	dalpha += y*exp( (-1/alpha)*lay1)/ (alpha*ay1);
	dalpha -= lay1 * exp( (-1/alpha)*lay1) / (alpha*alpha);
      }
    }
  dmu *= lambda;

  /* Return the negative gradient, because we're minimizing NLL,
   * not maximizing LL.
   */
  dp[0] = -dmu;
  dp[1] = -dw;
  dp[2] = -dalpha;
  return;
}

/*****************************************************************
 * temporary, while using NR to debug
 */
#if 0
static struct gev_data nrdata;
static double 
nr_func(double *p)
{
  return (gev_complete_func(p, 3, (void*) &nrdata));
}
static void
nr_dfunc(double *p, double *dp)
{
  gev_complete_grad(p, 3, (void *) &nrdata, dp);
}
#endif

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

int
esl_gev_FitComplete(double *x, int n, 
		    double *ret_mu, double *ret_lambda, double *ret_alpha)
{
  struct gev_data data;
  double fx;
  double p[3];
  double u[3];
  double wrk[12];		/* 4 vectors of length 3 */
  int    status;
  double mean, variance;
  double mu, lambda, alpha;

  data.x   = x;
  data.n   = n;
  data.phi = -DBL_MAX;
  data.z   = 0;

  mean_and_variance(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);
  mu     = mean - 0.57722/lambda;
  alpha  = 0.0001;

 /* make 'em up, for now */
  p[0] = mu;
  p[1] = log(lambda);
  p[2] = alpha;

  /* pass to the optimizer
   */
#if 0
  int    niter;
  nrdata.x   = x;
  nrdata.n   = n;
  nrdata.phi = -DBL_MAX;
  nrdata.z   = 0;


  ConjugateGradientDescent(p, 3, 1e-5, &niter, &fx, &nr_func, &nr_dfunc, NULL);

  /* make 'em up, for now */
  p[0] = mu;
  p[1] = log(lambda);
  p[2] = alpha;
#endif



  /* initial step sizes */
  u[0] = 1.0;
  u[1] = fabs(log(0.02));
  u[2] = 0.02;

  status = esl_min_ConjugateGradientDescent(p, u, 3, 
					    &gev_complete_func, 
					    &gev_complete_grad,
					    (void *)(&data),
					    1e-7, wrk, &fx);
  *ret_mu     = p[0];
  *ret_lambda = exp(p[1]);
  *ret_alpha  = p[2];
  return status;
}
#endif /*eslAUGMENT_MINIMIZER*/
/*--------------------------- end fitting ----------------------------------*/


/****************************************************************************
 * Example, test, and stats drivers
 ****************************************************************************/ 
/* Example main()
 */
#ifdef eslGEV_EXAMPLE
/*::cexcerpt::gev_example::begin::*/
/* compile: 
     gcc -g -Wall -I. -o example -DeslGEV_EXAMPLE -DeslAUGMENT_RANDOM\
       -DeslAUGMENT_MINIMIZER esl_gev.c esl_random.c esl_minimizer.c\
       esl_vectorops.c easel.c -lm
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_gev.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r = esl_randomness_CreateTimeseeded();;
  int     n         = 10000; 	   /* simulate 10,000 samples */
  double  mu        = -20.0;       /* with mu = -20    */ 
  double  lambda    = 0.4;         /* and lambda = 0.4 */
  double  alpha     = 0.1;	   /* and alpha = 0.1  */
  double  min       =  9999.;
  double  max       = -9999.;
  double *x         = malloc(sizeof(double) * n);
  double est_mu, est_lambda, est_alpha;
  double  z;
  double  nll;
  int     i;

  nll = 0.;
  for (i = 0; i < n; i++)	/* generate the 10,000 samples */
    { 
      x[i] = esl_gev_Sample(r, mu, lambda, alpha);
      nll -= esl_gev_logpdf(x[i], mu, lambda, alpha);
      if (x[i] < min) min = x[i];
      if (x[i] > max) max = x[i];
    }

  z = esl_gev_surv(max, mu, lambda, alpha);       /* right tail p~1e-4 >= max */
  printf("max = %6.1f  P(>max)  = %g   E=%6.3f\n", max, z, z*(double)n);
  z = esl_gev_cdf(min, mu, lambda, alpha);        /* left tail p~1e-4 < min */
  printf("min = %6.1f  P(<=min) = %g   E=%6.3f\n", min, z, z*(double)n);

  esl_gev_FitComplete(x, n, &est_mu, &est_lambda, &est_alpha);
 
  z = 100. * fabs((est_mu - mu) / mu);
  printf("Parametric mu     = %6.1f.  Estimated mu     = %6.2f.  Difference = %.1f%%.\n",
	 mu, est_mu, z);
  z = 100. * fabs((est_lambda - lambda) /lambda);
  printf("Parametric lambda = %6.2f.  Estimated lambda = %6.2f.  Difference = %.1f%%.\n",
	 lambda, est_lambda, z);
  z = 100. * fabs((est_alpha - alpha) /alpha);
  printf("Parametric alpha  = %6.4f.  Estimated alpha  = %6.4f.  Difference = %.1f%%.\n",
	 alpha, est_alpha, z);

  z = mu + (exp(-alpha*log(1/(double)n)) - 1 ) / (alpha*lambda);/* x at E=1*/
  z = (double) n * esl_gev_surv(z, est_mu, est_lambda, est_alpha); /* E at x */
  printf("Estimated E of x at true E=1: %6.4f\n", z);

  printf("NLL at true parameters: %6.4f\n", nll);

  free(x);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::gev_example::end::*/
#endif /*eslGEV_EXAMPLE*/


#ifdef eslGEV_STATS
/* compile: 
     gcc -g -Wall -I. -o stats -DeslGEV_STATS -DeslAUGMENT_RANDOM\
       -DeslAUGMENT_MINIMIZER esl_gev.c esl_random.c esl_minimizer.c\
       esl_vectorops.c easel.c -lm
 * run:     ./example
 */
#include <stdio.h>
#include <math.h>
#include <easel.h>
#include <esl_random.h>
#include <esl_minimizer.h>
#include <esl_gev.h>

static void stats_sample(FILE *fp);

int
main(int argc, char **argv)
{
  FILE *fp;
  double  mu        = -20.0;       /* with mu = -20    */ 
  double  lambda    = 0.4;         /* and lambda = 0.4 */
  double  xmin      = -40.;
  double  xmax      = 40.;
  double  xstep     = 0.1; 
  double  x,z;

  /* stats.1: xmgrace xy file w/ densities for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.1", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_pdf(x, mu, lambda, 0.0));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_pdf(x, mu, lambda, 0.6));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_pdf(x, mu, lambda, -0.6));
  fprintf(fp, "&\n");
  fclose(fp);

  /* stats.2: xmgrace xy file w/ log densities for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.2", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logpdf(x, mu, lambda, 0.0);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logpdf(x, mu, lambda, 0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logpdf(x, mu, lambda, -0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  fclose(fp);

  /* stats.3: xmgrace xy file w/ CDF for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.3", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_cdf(x, mu, lambda, 0.0));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_cdf(x, mu, lambda, 0.6));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_cdf(x, mu, lambda, -0.6));
  fprintf(fp, "&\n");
  fclose(fp);
  
  /* stats.4: xmgrace xy file w/ logCDF for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.4", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logcdf(x, mu, lambda, 0.0);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logcdf(x, mu, lambda, 0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logcdf(x, mu, lambda, -0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  fclose(fp);

 /* stats.5: xmgrace xy file w/ surv for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.5", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_surv(x, mu, lambda, 0.0));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_surv(x, mu, lambda, 0.6));
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep)
    fprintf(fp, "%.1f  %9.7f\n", x, esl_gev_surv(x, mu, lambda, -0.6));
  fprintf(fp, "&\n");
  fclose(fp);

 /* stats.6: xmgrace xy file w/ logsurv for Gumbel, Frechet, Weibull */
  if ((fp = fopen("stats.6", "w")) == NULL) abort();
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logsurv(x, mu, lambda, 0.0);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logsurv(x, mu, lambda, 0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  for (x = xmin; x <= xmax; x+= xstep) {
    z = esl_gev_logsurv(x, mu, lambda, -0.2);
    if (finite(z)) fprintf(fp, "%.1f  %9.7f\n", x, z);
  }
  fprintf(fp, "&\n");
  fclose(fp);

  /* stats.7. R input file of 10,000 random GEV samples.
   */
  if ((fp = fopen("stats.7", "w")) == NULL) abort();  
  stats_sample(fp);
  fclose(fp);
  return 0;
}

/* stats_sample()
 * Creates an R input table containing 10,000 random samples
 * each in columns labeled "gumbel", "frechet", "weibull".
 * To process in R (remember that R uses 1/lambda for scale):
     library(ismev)
     library(evd)
     z=read.table("stats.7")
     x1 <- sort(z$gumbel,  decreasing=T)
     x2 <- sort(z$frechet, decreasing=T)
     x3 <- sort(z$weibull, decreasing=T)
     q1 <- qgumbel(ppoints(10000), -20., 1./0.4)
     q2 <- qgev(ppoints(10000), -20., 1./0.4, 0.2)
     q3 <- qgev(ppoints(10000), -20., 1./0.4, -0.2)
     xax<- seq(-40,40,by=0.1)
     a1 <- dgumbel(xax, -20, 1/0.4)
     a2 <- dgev(xax, -20, 1/0.4, 0.2)
     a3 <- dgev(xax, -20, 1/0.4, -0.2)
     qqplot(x1,q1); abline(0,1)
     qqplot(x2,q2); abline(0,1)
     qqplot(x3,q3); abline(0,1)
     plot(density(x1,bw=0.2)); lines(xax,a1)
     plot(density(x2,bw=0.2)); lines(xax,a2)
     plot(density(x3,bw=0.2)); lines(xax,a3)
 */
static void
stats_sample(FILE *fp)
{
  ESL_RANDOMNESS *r;
  double mu     = -20.;
  double lambda = 0.4;
  int    n      = 10000;
  double a,b,c;
  int    i;

  r = esl_randomness_Create(42);
  fprintf(fp, "         gumbel  \t  frechet\t  weibull\n");
  for (i = 1; i <= n; i++)
    {
      a  = esl_gev_Sample(r, mu, lambda, 0.0);
      b  = esl_gev_Sample(r, mu, lambda, 0.2);
      c  = esl_gev_Sample(r, mu, lambda, -0.2);
      fprintf(fp, "%d\t%8.4f\t%8.4f\t%8.4f\n", i, a,b,c);
    }
  esl_randomness_Destroy(r);
}
#endif /*eslGEV_STATS*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
