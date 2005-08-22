/* esl_stretchexp.c
 * Statistical routines for stretched exponential distributions.
 * 
 * SRE, Fri Aug 19 11:15:21 2005 [St. Louis] 
 * xref STL9/146
 * SVN $Id
 */

#include <stdio.h>
#include <math.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_vectorops.h>  /* one DMin() call that could be removed*/
#include <esl_stretchexp.h>

#include <esl_dirichlet.h>   /* just need SampleGamma() */

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_MINIMIZER
#include <esl_minimizer.h>
#endif
#ifdef eslAUGMENT_HISTOGRAM
#include <esl_histogram.h>
#endif
/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 

/* Function:  esl_sxp_pdf()
 * Incept:    SRE, Fri Aug 19 11:17:47 2005 [St. Louis]
 *
 * Purpose:   Calculates the probability density function for the 
 *            stretched exponential pdf, $P(X=x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_sxp_pdf(double x, double mu, double lambda, double tau)
{
  double y    = lambda * (x-mu);
  double val;
  double gt;
  
  if (x < mu) return 0.;
  esl_stats_LogGamma(1/tau, &gt);
  val = (lambda * tau / exp(gt)) * exp(- exp(tau * log(y)));
  return val;
}

/* Function:  esl_sxp_logpdf()
 * Incept:    SRE, Fri Aug 19 11:27:32 2005 [St. Louis]
 *
 * Purpose:   Calculates the log probability density function for the 
 *            stretched exponential pdf, $\log P(X=x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double 
esl_sxp_logpdf(double x, double mu, double lambda, double tau)
{
  double y    = lambda * (x-mu);
  double gt;
  double val;

  if (x < mu) return -eslINFINITY;
  esl_stats_LogGamma(1/tau, &gt);
  val = log(lambda) + log(tau) - gt - exp(tau*log(y));
  return val;
}

/* Function:  esl_sxp_cdf()
 * Incept:    SRE, Fri Aug 19 11:30:55 2005 [St. Louis]
 *
 * Purpose:   Calculates the cumulative distribution function for the 
 *            stretched exponential pdf, $P(X \leq x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_sxp_cdf(double x, double mu, double lambda, double tau)
{
  double y = lambda * (x-mu);
  double val;

  if (x <= mu) return 0.;
  esl_stats_IncompleteGamma(1/tau, exp(tau * log(y)), &val);
  
  ESL_DASSERT1 (( !isnan(val)));
  return (1-val);
}

/* Function:  esl_sxp_logcdf()
 * Incept:    SRE, Fri Aug 19 11:37:20 2005 [St. Louis]
 *
 * Purpose:   Calculates the log of the cumulative distribution function for the 
 *            stretched exponential pdf, $\log P(X \leq x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_sxp_logcdf(double x, double mu, double lambda, double tau)
{
  double y = lambda * (x-mu);
  double val;

  if (x <= mu) return -eslINFINITY;
  esl_stats_IncompleteGamma(1./tau, exp(tau * log(y)), &val);
  return log(1 - val);
}



/* Function:  esl_sxp_surv()
 * Incept:    SRE, Fri Aug 19 11:38:24 2005 [St. Louis]
 *
 * Purpose:   Calculates the survival function for the 
 *            stretched exponential pdf, $P(X > x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_sxp_surv(double x, double mu, double lambda, double tau)
{
  double y = lambda * (x-mu);
  double val;

  if (x <= mu) return 1.0;

  esl_stats_IncompleteGamma(1./tau, exp(tau * log(y)), &val);
  return val;
}

/* Function:  esl_sxp_logsurv()
 * Incept:    SRE, Fri Aug 19 11:38:24 2005 [St. Louis]
 *
 * Purpose:   Calculates the log survival function for the 
 *            stretched exponential pdf, $\log P(X > x)$, given
 *            quantile <x>, offset <mu>, and parameters <lambda> and <tau>.
 */
double
esl_sxp_logsurv(double x, double mu, double lambda, double tau)
{
  double y = lambda * (x-mu);
  double val;

  if (x <= mu) return 0.0;

  esl_stats_IncompleteGamma(1./tau, exp(tau * log(y)), &val);
  return log(val);
}

/* Function:  esl_sxp_invcdf()
 * Incept:    SRE, Sat Aug 20 14:42:06 2005 [St. Louis]
 *
 * Purpose:   Calculates the inverse CDF for a stretched exponential
 *            with parameters <mu>, <lambda>, and <tau>, returning
 *            the quantile <x> at which the CDF is <p>.
 *            
 *            The inverse CDF of the stretched exponential has no
 *            analytical expression as far as I'm aware. The calculation
 *            here is a computationally expensive, brute force bisection
 *            search in <x> using the CDF function. It will suffice for
 *            a small number of calls (for plotting applications, for example),
 *            but it is not sufficient for a large number of calls.
 */
double
esl_sxp_invcdf(double p, double mu, double lambda, double tau)
{
  double x1, x2, xm;		/* low, high guesses at x */
  double f1, f2, fm;
  double tol = 1e-6;

  x1 = mu;
  f1 = 0.;
  x2 = mu + 1.;
  do {				/* bracket */
    x2 = x2 + 2.*(x2-x1);
    f2 = esl_sxp_cdf(x2, mu, lambda, tau);
  } while (f2 < p);

  do {				/* bisection */
    xm = (x1+x2) / 2.;
    fm = esl_sxp_cdf(xm, mu, lambda, tau);
    
    if      (fm > p) x2 = xm;
    else if (fm < p) x1 = xm;
    else return xm;		/* unlikely case of fm==cdf */
  } while ( (x2-x1)/(x1+x2-2*mu) > tol);

  xm = (x1+x2) / 2.;
  return xm;
}
/*-------------------- end densities & distributions ------------------------*/
	  



/****************************************************************************
 * Generic API routines: for general interface w/ histogram module
 ****************************************************************************/ 

/* Function:  esl_sxp_generic_cdf()
 * Incept:    SRE, Fri Aug 19 13:54:26 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_sxp_cdf()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_sxp_generic_cdf(double x, void *params)
{
  double *p = (double *) params;
  return esl_sxp_cdf(x, p[0], p[1], p[2]);
}

/* Function:  esl_sxp_generic_invcdf()
 * Incept:    SRE, Sat Aug 20 14:46:55 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_sxp_invcdf()>, taking
 *            a void ptr to a double array containing $\mu$, $\lambda$,
 *            $\tau$ parameters.
 */
double
esl_sxp_generic_invcdf(double p, void *params)
{
  double *v = (double *) params;
  return esl_sxp_invcdf(p, v[0], v[1], v[2]);
}
/*------------------------ end generic API ---------------------------------*/



/****************************************************************************
 * Routines for dumping plots for files
 ****************************************************************************/ 

/* Function:  esl_sxp_Plot()
 * Incept:    SRE, Fri Aug 19 11:48:27 2005 [St. Louis]
 *
 * Purpose:   Plot some stretched exponential function <func> (for instance,
 *            <esl_sxp_pdf()>) for parameters <mu>, <lambda>, and <tau>, for
 *            a range of quantiles x from <xmin> to <xmax> in steps of <xstep>;
 *            output to an open stream <fp> in xmgrace XY input format.
 *
 * Returns:   <eslOK>.
 */
int
esl_sxp_Plot(FILE *fp, double mu, double lambda, double tau,
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
 * Routines for sampling (requires augmentation w/ random, dirichlet modules)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_sxp_Sample()
 * Incept:    SRE, Fri Aug 19 13:39:36 2005 [St. Louis]
 *
 * Purpose:   Sample a stretched exponential random variate,
 *            by a change of variable from a Gamma sample.
 */
double
esl_sxp_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau)
{
  double t,x;

  esl_dirichlet_SampleGamma(r, 1./tau, &t); 
  x = mu + 1./lambda * exp(1./tau * log(t));
  return x;
} 
#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/



/****************************************************************************
 * Maximum likelihood fitting
 ****************************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER
/* This structure is used to sneak the data into minimizer's generic
 * (void *) API for all aux data
 */
struct sxp_data {
  double *x;
  int     n;
  double  mu;
};

static double
sxp_complete_func(double *p, int np, void *dptr)
{
  struct sxp_data *data = (struct sxp_data *) dptr;
  double lambda, tau;
  double logL = 0.;
  int    i;

  lambda = exp(p[0]);
  tau    = exp(p[1]);

  for (i = 0; i < data->n; i++)
    logL += esl_sxp_logpdf(data->x[i], data->mu, lambda, tau);
  return -logL;
}

/* Function:  esl_sxp_FitComplete()
 * Incept:    SRE, Fri Aug 19 15:25:42 2005 [St. Louis]
 *
 * Purpose:   Given a vector of <n> observed data samples <x[]>,
 *            find maximum likelihood parameters by conjugate gradient 
 *            descent optimization.
 */
int
esl_sxp_FitComplete(double *x, int n,
		    double *ret_mu, double *ret_lambda, double *ret_tau)

{
  struct sxp_data data;
  double p[2], u[2], wrk[8];
  double mu, tau, lambda;
  double tol = 1e-6;
  double fx;
  int    status;

  /* initial guesses */
  mu     = esl_vec_DMin(x, n); 
  lambda = 1.0;
  tau    = 0.42;

  /* load data structure, param vector, and step vector */
  data.x  = x;
  data.n  = n;
  data.mu = mu;
  p[0]    = log(lambda);
  p[1]    = log(tau);
  u[0]    = 1.0;
  u[1]    = 1.0;

  /* hand it off */
  status =  esl_min_ConjugateGradientDescent(p, u, 2, 
					     &sxp_complete_func, 
					     NULL,
					     (void *) (&data), tol, wrk, &fx);
  
  *ret_mu     = mu;
  *ret_lambda = exp(p[0]);
  *ret_tau    = exp(p[1]);
  return eslOK;
}

#ifdef eslAUGMENT_HISTOGRAM
struct sxp_binned_data {
  ESL_HISTOGRAM *g;	/* contains the binned data    */
  double mu;		/* mu is not a learnable param */
};

static double 
sxp_complete_binned_func(double *p, int np, void *dptr)
{
  struct sxp_binned_data *data = (struct sxp_binned_data *) dptr;
  ESL_HISTOGRAM          *g    = data->g;
  double logL = 0.;
  double ai, bi;		/* lower, upper bounds on bin */
  double lambda, tau;
  int    i;
  double tmp;

  lambda = exp(p[0]);
  tau    = exp(p[1]);  

  ESL_DASSERT1(( ! isnan(lambda) ));
  ESL_DASSERT1(( ! isnan(tau) ));
  
  for (i = g->imin; i <= g->imax; i++) /* for each occupied bin */
    {
      if (g->obs[i] == 0) continue;
      
      esl_histogram_GetBinBounds(g, i, &ai, &bi, NULL);
      if (ai < data->mu) ai = data->mu; /* careful at leftmost bound */

      tmp = esl_sxp_cdf(bi, data->mu, lambda, tau) -
            esl_sxp_cdf(ai, data->mu, lambda, tau);
      if      (tmp == 0.) return eslINFINITY;
      logL += g->obs[i] * log(tmp);
    }
  return -logL;			/* minimizing NLL */
}

/* Function:  esl_sxp_FitCompleteBinned()
 * Incept:    SRE, Sat Aug 20 13:28:00 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <g> with binned observations, where each
 *            bin i holds some number of observed samples x with values from 
 *            lower bound l to upper bound u (that is, $l < x \leq u$);
 *            find maximum likelihood parameters mu, lambda, tau by conjugate
 *            gradient descent optimization.
 */
int
esl_sxp_FitCompleteBinned(ESL_HISTOGRAM *g,
			  double *ret_mu, double *ret_lambda, double *ret_tau)

{
  struct sxp_binned_data data;
  double p[2], u[2], wrk[8];
  double mu, tau, lambda;
  double tol = 1e-6;
  double fx;
  int    status;

  /* initial guesses are arbitrary */
  mu     = g->xmin; 	/* fix mu here; no point in optimizing it */
  lambda = 1.0;
  tau    = 0.42;

  /* load data structure, param vector, and step vector */
  data.g  = g;
  data.mu = mu;
  p[0]    = log(lambda);
  p[1]    = log(tau);
  u[0]    = 1.0;
  u[1]    = 1.0;

  /* hand it off */
  status =  esl_min_ConjugateGradientDescent(p, u, 2, 
					     &sxp_complete_binned_func, 
					     NULL,
					     (void *) (&data), tol, wrk, &fx);
  *ret_mu     = mu;
  *ret_lambda = exp(p[0]);
  *ret_tau    = exp(p[1]);
  return eslOK;
}
#endif /*eslAUGMENT_HISTOGRAM*/
#endif /*eslAUGMENT_MINIMIZER*/

/****************************************************************************
 * Example, test, and stats drivers
 ****************************************************************************/ 
/* Example main()
 */
#ifdef eslSXP_EXAMPLE
/*::cexcerpt::sxp_example::begin::*/
/* compile:
   gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o example -DeslSXP_EXAMPLE\
     esl_stretchexp.c -leasel -lm
 */
#include <stdio.h>
#include <easel.h>
#include <esl_stats.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_stretchexp.h>

int
main(int argc, char **argv)
{
  ESL_HISTOGRAM  *h;
  ESL_RANDOMNESS *r;
  double mu     = -50.0;
  double lambda = 2.5;
  double tau    = 0.7;
  double emu, elam, etau;
  int    n      = 10000;
  int    i;
  double x;

  r = esl_randomness_CreateTimeseeded();
  h = esl_histogram_CreateFull(mu, 100., 0.1);
  for (i = 0; i < n; i++)
    {
      x = esl_sxp_Sample(r, mu, lambda, tau);
      esl_histogram_Add(h, x);
    }
  esl_histogram_Sort(h);

  /* Plot the empirical (sampled) and expected survivals */
  esl_histogram_PlotSurvival(stdout, h);
  esl_sxp_Plot(stdout, mu, lambda, tau,
	       &esl_sxp_surv,  h->xmin, h->xmax, 0.1);

  /* ML fit to complete data, and plot fitted survival curve */
  esl_sxp_FitComplete(h->x, h->n, &emu, &elam, &etau);
  esl_sxp_Plot(stdout, emu, elam, etau,
	       &esl_sxp_surv,  h->xmin, h->xmax, 0.1);

  /* ML fit to binned data, plot fitted survival curve  */
  esl_sxp_FitCompleteBinned(h, &emu, &elam, &etau);
  esl_sxp_Plot(stdout, emu, elam, etau,
	       &esl_sxp_surv,  h->xmin, h->xmax, 0.1);

  return 0;
}
/*::cexcerpt::sxp_example::end::*/
#endif /*eslSXP_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
