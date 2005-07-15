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

  if (ya1 <= 0) return -DBL_MAX;
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

  if (ya1 <= 0) return 0.;
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

  if (ya1 <= 0) return -DBL_MAX;
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
   
   if (ya1 <= 0) return 1.0;
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
   if (ya1 <= 0) return 0.;
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
  if (alpha < 1e-12) return (mu - log(-1. * log(p)) / lambda);

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
      if (ay < 1e-12) dmu -= exp(-y);
      else            dmu -= exp(-(1+1/alpha) * lay1);

      /* d/dw, term2. converges to -y, for small alpha. */
      dw -= y*(1+alpha) / ay1;

      /* d/dw, term2. For tiny ay, use log(1+x) ~ x to simplify. */
      if (ay < 1e-12) dw += y*exp(-y);
      else            dw += y*exp(-(1+1/alpha) * lay1);

      /* d/dalpha, term1
       */
      dalpha -= (1 + 1/alpha) * y/ay1;

      /* d/dalpha, terms 2,3,4: for tiny ay, simplify.
       * d/dalpha will go to +/-inf for alpha ~ 0, so watch out.
       */
      if (ay < 1e-12) {
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
  int    niter;
  double mean, variance;
  double mu, lambda, alpha;

  nrdata.x   = x;
  nrdata.n   = n;
  nrdata.phi = -DBL_MAX;
  nrdata.z   = 0;

  data.x   = x;
  data.n   = n;
  data.phi = -DBL_MAX;
  data.z   = 0;

  mean_and_variance(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);
  mu     = mean - 0.57722/lambda;
  alpha  = 0.01;

 /* make 'em up, for now */
  p[0] = mu;
  p[1] = log(lambda);
  p[2] = alpha;

  /* pass to the optimizer
   */
  ConjugateGradientDescent(p, 3, 1e-5, &niter, &fx, &nr_func, &nr_dfunc, NULL);
  printf("fx = %f\n", fx);

  /* make 'em up, for now */
  p[0] = mu;
  p[1] = log(lambda);
  p[2] = alpha;

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
/* compile: gcc -g -Wall -I. -o example -DeslGEV_EXAMPLE esl_gev.c esl_vectorops.c easel.c -lm
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <esl_gev.h>

int
main(int argc, char **argv)
{
  double x;
  double a,b,c,d;
  double mu     = 0.;
  double lambda = 1.0;
  double alpha  = .1;
  double min    = -5.;
  double max    = 100.;
  double step   = 0.2;

  for (x = min; x <= max; x += step)
    {
      /*pdf   = log(esl_gev_pdf(x, mu, lambda, alpha)); */
      a  = esl_gev_cdf(x, mu, lambda, alpha); 
      b  = esl_evd_cdf(x, mu, lambda); 
      printf("%g    %g    %g   \n", x, a, b);
    }
  printf ("&\n");
  return 0;
}
/*::cexcerpt::gev_example::end::*/
#endif /*eslGEV_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/