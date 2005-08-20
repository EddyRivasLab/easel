/* esl_hyperexp.c 
 * Statistical routines for hyperexponential distributions.
 * 
 * SRE, Mon Aug 15 08:29:45 2005  xref:STL9/140  [St. Louis]
 * SVN $Id$
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <easel.h>
#include <esl_vectorops.h>
#include <esl_exponential.h>
#include <esl_hyperexp.h>

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

#ifdef eslAUGMENT_MINIMIZER
#include <esl_histogram.h>
#include <esl_minimizer.h>
#endif

/****************************************************************************
 * Routines for the ESL_HYPEREXP object
 ****************************************************************************/ 

/* Function:  esl_hyperexp_Create()
 * Incept:    SRE, Mon Aug 15 08:40:44 2005 [St. Louis]
 *
 * Purpose:   Creates an object to hold parameters for a <K>-component
 *            hyperexponential. 
 *
 *            Parameters in the object are initialized 
 *            ($q_k = \frac{1}{K}$, $\lambda_k = 1$, $\mu = 0$), but
 *            the caller will want to set these according to its own
 *            purposes.
 *
 * Args:      K  - number of components in the mixture
 *
 * Returns:   ptr to newly allocated/initialized <ESL_HYPEREXP> object.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_HYPEREXP *
esl_hyperexp_Create(int K)
{
  ESL_HYPEREXP *h;
  int           k;

  if ((h = malloc(sizeof(ESL_HYPEREXP))) == NULL)         goto FAILURE;
  h->q = h->lambda = h->wrk = NULL;
  h->fixlambda = NULL;
  h->K = K;

  if ((h->q        = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((h->lambda   = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((h->wrk      = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((h->fixlambda= malloc(sizeof(char)   * K)) == NULL) goto FAILURE;

  for (k = 0; k < K; k++)
    {
      h->q[k]        = 1. / (double) K;
      h->lambda[k]   = 1.;
      h->fixlambda[k]= 0;
    }
  h->mu = 0.;
  return h;
  
 FAILURE:
  esl_hyperexp_Destroy(h);
  ESL_ERROR_NULL(eslEMEM, "malloc failed");
}

/* Function:  esl_hyperexp_Destroy()
 * Incept:    SRE, Mon Aug 15 08:53:50 2005 [St. Louis]
 *
 * Purpose:   Deallocates the hyperexponential parameter object <h>.
 *
 * Args:      h  - ptr to the object to be deallocated.
 *
 * Returns:   (void).
 */
void
esl_hyperexp_Destroy(ESL_HYPEREXP *h)
{
  if (h == NULL) return;

  if (h->q        != NULL) free(h->q);
  if (h->lambda   != NULL) free(h->lambda);
  if (h->wrk      != NULL) free(h->wrk);
  if (h->fixlambda!= NULL) free(h->fixlambda);
  free(h);
}
  

/* Function:  esl_hyperexp_Copy()
 * Incept:    SRE, Mon Aug 15 08:58:17 2005 [St. Louis]
 *
 * Purpose:   Makes a copy of the hyperexponential parameter object <src>
 *            in <dest>. Caller must have already allocated <dest> to have
 *            (at least) the same number of components as <src>.
 *
 * Args:      src   - object to be copied
 *            dest  - allocated object to copy <src> into
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINCOMPAT> if <dest> isn't allocated with enough
 *            components to hold a copy of <src>.
 */
int
esl_hyperexp_Copy(ESL_HYPEREXP *src, ESL_HYPEREXP *dest)
{
  int k;

  if (dest->K < src->K) 
    ESL_ERROR(eslEINCOMPAT, "hyperexponential too small to copy into");

  for (k = 0; k < src->K; k++)
    {
      dest->q[k]        = src->q[k];
      dest->lambda[k]   = src->lambda[k];
      dest->fixlambda[k]= src->fixlambda[k];
    }
  dest->mu = src->mu;
  dest->K  = src->K;
  return eslOK;
}
/*----------------- end ESL_HYPEREXP object maintenance --------------------*/






/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 

/* Function:  esl_hxp_pdf()
 * Incept:    SRE, Mon Aug 15 09:17:34 2005 [St. Louis]
 *
 * Purpose:   Returns the probability density function $P(X=x)$ for
 *            quantile <x>, given hyperexponential parameters <h>.
 */
double
esl_hxp_pdf(double x, ESL_HYPEREXP *h)
{
  double pdf = 0.;
  int    k;

  if (x < h->mu) return 0.;

  for (k = 0; k < h->K; k++)
    pdf += h->q[k] * esl_exp_pdf(x, h->mu, h->lambda[k]);
  return pdf;
}


/* Function:  esl_hxp_logpdf()
 * Incept:    SRE, Mon Aug 15 09:25:45 2005 [St. Louis]
 *
 * Purpose:   Returns the log of the PDF ($\log P(X=x)$) for quantile <x>,
 *            given hyperexponential parameters <h>.
 */
double
esl_hxp_logpdf(double x, ESL_HYPEREXP *h)
{
  int    k;
  double z;

  if (x < h->mu) return -eslINFINITY;

  for (k = 0; k < h->K; k++)
    if (h->q[k] == 0.0) 
      h->wrk[k] = -eslINFINITY;	
    else
      h->wrk[k] = log(h->q[k]) + esl_exp_logpdf(x, h->mu, h->lambda[k]);

  z = esl_vec_DLogSum(h->wrk, h->K);
  return z;
}

/* Function:  esl_hxp_cdf()
 * Incept:    SRE, Mon Aug 15 09:48:44 2005 [St. Louis]
 *
 * Purpose:   Returns the cumulative distribution function $P(X \leq x)$
 *            for quantile <x>, given hyperexponential parameters <h>.
 */
double
esl_hxp_cdf(double x, ESL_HYPEREXP *h)
{
  double cdf = 0.;
  int    k;
  
  if (x < h->mu) return 0.;

  for (k = 0; k < h->K; k++)
    cdf += h->q[k] * esl_exp_cdf(x, h->mu, h->lambda[k]);
  return cdf;
}

/* Function:  esl_hxp_logcdf()
 * Incept:    SRE, Mon Aug 15 09:52:31 2005 [St. Louis]
 *
 * Purpose:   Returns the log of the CDF $\log P(X \leq x)$
 *            for quantile <x>, given hyperexponential parameters <h>.
 */
double
esl_hxp_logcdf(double x, ESL_HYPEREXP *h)
{
  int k;

  if (x < h->mu) return -eslINFINITY;

  for (k = 0; k < h->K; k++)
    if (h->q[k] == 0.0) 
      h->wrk[k] = -eslINFINITY;
    else
      h->wrk[k] = log(h->q[k]) + esl_exp_logcdf(x, h->mu, h->lambda[k]);

  return esl_vec_DLogSum(h->wrk, h->K);
}


/* Function:  esl_hxp_surv()
 * Incept:    SRE, Mon Aug 15 09:57:39 2005 [St. Louis]
 *
 * Purpose:   Returns the survivor function $P(X > x)$ (1-CDF)
 *            for quantile <x>, given hyperexponential parameters <h>.
 */
double
esl_hxp_surv(double x, ESL_HYPEREXP *h)
{
  double srv = 0.;
  int    k;
  
  if (x < h->mu) return 1.0;

  for (k = 0; k < h->K; k++)
    srv += h->q[k] * esl_exp_surv(x, h->mu, h->lambda[k]);
  return srv;
}

  
/* Function:  esl_hxp_logsurv()
 * Incept:    SRE, Mon Aug 15 10:00:46 2005 [St. Louis]
 *
 * Purpose:   Returns the log survivor function $\log P(X > x)$ (log(1-CDF))
 *            for quantile <x>, given hyperexponential parameters <h>.
 */
double
esl_hxp_logsurv(double x, ESL_HYPEREXP *h)
{
  int k;
  
  if (x < h->mu) return 0.0;

  for (k = 0; k < h->K; k++)
    if (h->q[k] == 0.0) 
      h->wrk[k] = -eslINFINITY;
    else
      h->wrk[k] = log(h->q[k]) + esl_exp_logsurv(x, h->mu, h->lambda[k]);
  
  return esl_vec_DLogSum(h->wrk, h->K);
}
/*-------------------- end densities & distributions ---------------------------*/




/****************************************************************************
 * Generic API routines: for general interface w/ histogram module
 ****************************************************************************/ 

/* Function:  esl_hxp_generic_cdf()
 * Incept:    SRE, Wed Aug 17 13:17:08 2005 [St. Louis]
 *
 * Purpose:   Generic-API version of CDF call, for passing to histogram module's
 *            <SetExpected()> and <Goodness()>.
 */
double
esl_hxp_generic_cdf(double x, void *params)
{
  ESL_HYPEREXP *h = (ESL_HYPEREXP *) params;
  return esl_hxp_cdf(x, h);
}

/*------------------------ end generic API ---------------------------------*/






/****************************************************************************
 * Routines for dumping plots for files
 ****************************************************************************/ 

/* Function:  esl_hxp_Plot()
 * Incept:    SRE, Mon Aug 15 10:06:35 2005 [St. Louis]
 *
 * Purpose:   Plot some function <func> (for instance, <esl_hxp_pdf()>)
 *            for hyperexponential parameters <h>, for a range of
 *            quantiles x from <xmin> to <xmax> in steps of <xstep>;
 *            output to an open stream <fp> in xmgrace XY input format.
 *
 * Returns:   <eslOK>.
 */
int
esl_hxp_Plot(FILE *fp, ESL_HYPEREXP *h,
	     double (*func)(double x, ESL_HYPEREXP *h), 
	     double xmin, double xmax, double xstep)
{
  double x;
  for (x = xmin; x <= xmax; x += xstep)
    fprintf(fp, "%f\t%g\n", x, (*func)(x, h));
  fprintf(fp, "&\n");
  return eslOK;
}
/*-------------------- end plot dumping routines ---------------------------*/





/****************************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_hxp_Sample()
 * Incept:    SRE, Mon Aug 15 10:10:26 2005 [St. Louis]
 *
 * Purpose:   Sample a random variate x from a hyperexponential <h>, 
 *            given random number source <r>.
 */
double
esl_hxp_Sample(ESL_RANDOMNESS *r, ESL_HYPEREXP *h)
{
  int k;	
  k = esl_rnd_DChoose(r, h->q, h->K);
  return esl_exp_Sample(r, h->mu, h->lambda[k]);
}
#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/




/****************************************************************************
 * Maximum likelihood fitting, complete unbinned data
 ****************************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER
/* This structure is used to sneak the data into minimizer's generic
 * (void *) API for all aux data
 */
struct hyperexp_data {
  double *x;
  int     n;
  ESL_HYPEREXP *h;
};

/* Given hyperexponential parameters in <h>;
 * do appropriate c.o.v.'s to unconstrained real parameters
 * and fill in the packed parameter vector <p>.
 * 
 * <p> must be allocated for at least (2K-1) doubles: K-1 mixture 
 * coefficients and K lambda parameters. (mu is not a free param).
 *
 * First K-1 are Q_1..Q_K-1 mixture coefficient parameters; Q_0 implicitly 0;
 *  cov is q_k = e^{Q_k} / \sum_j e^{Q_j};    Q_k = log(q_k) - log(q_0).
 * Then K lambda params;
 * lambda cov is lambda = e^w, w = log(lambda).
 */
static void
hyperexp_pack_paramvector(double *p, int np, ESL_HYPEREXP *h)
{
  int    i;			/* counter in parameter vector p */
  int    k;			/* counter in mixture components */
  double z;			/* tmp variable */

  /* mixture coefficients */
  z = log(h->q[0]);
  i = 0;
  for (k = 1; k < h->K; k++) 
    p[i++] = log(h->q[k]) - z;
  
  /* exponential parameters */
  for (k = 0; k < h->K; k++)
    if (! h->fixlambda[k])
      p[i++] = log(h->lambda[k]);
  /* you can assert(i==np) in debugging if you want */
}

/* Same as above but in reverse: given parameter vector <p>,
 * <np> = 2K-1, do appropriate c.o.v. back to desired parameter space, and
 * fill in the hyperexponential <h>.
 */
static void
hyperexp_unpack_paramvector(double *p, int np, ESL_HYPEREXP *h)
{
  int    i;			/* counter in parameter vector p */
  int    k;			/* counter in mixture components */
  double z;			/* tmp variable  */

  /* Fetch the params in their c.o.v. space first
   */
  i = 0;
  h->q[0] = 0;	/* implicitly */
  for (k = 1; k < h->K; k++) 
    h->q[k] = p[i++]; 
  for (k = 0; k < h->K; k++)
    if (! h->fixlambda[k]) 
      h->lambda[k] = p[i++];
  /* assert(i==np) */
  
  /* Convert mix coefficients back to probabilities;
   * their  c.o.v. is q_k = e^{Q_k} / \sum_k e^{Q_k}
   * which rearranges to exp(Q_k - log[\sum_k e^Q_k]),
   * and we have the DLogSum() function to compute the log sum.
   */
  z = esl_vec_DLogSum(h->q, h->K);
  for (k = 0; k < h->K; k++)
    h->q[k] = exp(h->q[k] - z);
  
  /* lambda c.o.v. is \lambda = e^w */
  for (k = 0; k < h->K; k++)
    if (! h->fixlambda[k]) 
      h->lambda[k] = exp(h->lambda[k]);
}

/* The log likelihood function to be optimized by ML fitting:
 *   This needs to be careful of a case where a lambda = inf.
 */
static double
hyperexp_complete_func(double *p, int np, void *dptr)
{
  struct hyperexp_data *data = (struct hyperexp_data *) dptr;
  ESL_HYPEREXP         *h    = data->h;
  double logL = 0.;
  int    i;

  hyperexp_unpack_paramvector(p, np, h);
  for (i = 0; i < data->n; i++)
    logL += esl_hxp_logpdf(data->x[i], h);
  return -logL;
}

/* The gradient of the NLL w.r.t. each free parameter in p.
 */
static void
hyperexp_complete_gradient(double *p, int np, void *dptr, double *dp)
{
  struct hyperexp_data *data = (struct hyperexp_data *) dptr;
  ESL_HYPEREXP         *h    = data->h;
  double pdf;
  int i,k;
  int pidx;			
  
  hyperexp_unpack_paramvector(p, np, h);
  esl_vec_DSet(dp, np, 0.);
  for (i = 0; i < data->n; i++)
    {
      /* FIXME: I think the calculation below needs to be done
       * in log space, to avoid underflow errors; see complete_binned_gradient()
       */
      /* Precalculate q_k PDF_k(x) terms, and their sum */
      for (k = 0; k < h->K; k++)
	h->wrk[k] = h->q[k] * esl_exp_pdf(data->x[i], h->mu, h->lambda[k]);
      pdf = esl_vec_DSum(h->wrk, h->K);

      pidx = 0;
      for (k = 1; k < h->K; k++) /* generic d/dQ solution for mixture models */
	dp[pidx++] -= h->wrk[k]/pdf - h->q[k];
      
      for (k = 0; k < h->K; k++)
	if (! h->fixlambda[k])
	  dp[pidx++] -= (1.-h->lambda[k]*(data->x[i]-h->mu))*h->wrk[k]/pdf; /* d/dw */
    }
}


/* Function:  esl_hyperexp_FitGuess()
 * Incept:    SRE, Mon Aug 15 14:02:02 2005 [St. Louis]
 *
 * Purpose:   Given a sorted vector of <n> observed data samples <x[]>,
 *            from smallest <x[0]> to largest <x[n-1]>, calculate a
 *            very crude guesstimate of a fit -- suitable only as a starting
 *            point for further optimization -- and return those parameters
 *            in <h>.
 *
 *            Assigns $q_k \propto \frac{1/k}$ and  $\mu = \min_i x_i$;
 *            splits $x$ into $K$ roughly equal-sized bins, and
 *            and assigns $\lambda_k$ as the ML estimate from bin $k$.
 */
int
esl_hxp_FitGuess(double *x, int n, ESL_HYPEREXP *h)
{
  double tmu;			/* current mu */
  double mean;			/* mean (x-tmu) in a bin */
  int    i,k;
  int    imin, imax;

  h->mu = x[0];
  for (k = 0; k < h->K; k++)
    {
      h->q[k] = 1 / (double)(k+1); /* priors ~ 1, 1/2, 1/3... */

      imin = (int) ((double)(k*n)/(double)h->K);
      imax = (int) ((double)((k+1)*n)/(double)h->K);
      tmu = x[imin];
      mean = 0.;
      for (i = imin; i < imax; i++)
	mean += x[i] - tmu;
      mean /= (double)(imax-imin);
      h->lambda[k] = 1 / mean;
    }
  esl_vec_DNorm(h->q, h->K);
  return eslOK;
}

/* Function:  esl_hxp_FitComplete()
 * Incept:    SRE, Mon Aug 15 14:11:19 2005 [St. Louis]
 *
 * Purpose:   Given a vector of <n> observed data samples <x[]> 
 *            (sorted or unsorted), and an initial guess <h> for
 *            a hyperexponential, find maximum likelihood parameters
 *            by conjugate gradient descent optimization, starting
 *            from <h> and leaving the final optimized solution in
 *            <h>.
 */
int
esl_hxp_FitComplete(double *x, int n, ESL_HYPEREXP *h)
{
  struct hyperexp_data data;
  int     status;
  double *p;
  double *u;
  double *wrk;
  double  tol;
  int     np;
  double  fx;
  int     i;

  tol = 1e-6;

  /* Determine number of free parameters and allocate 
   */
  np = h->K-1;	              /* K-1 mix coefficients...     */
  for (i = 0; i < h->K; i++)  /* ...and up to K lambdas free */
    if (! h->fixlambda[i]) np++;	
  p   = malloc(sizeof(double) * np);
  u   = malloc(sizeof(double) * np);
  wrk = malloc(sizeof(double) * np * 4);

  /* Copy shared info into the "data" structure
   */
  data.x   = x;
  data.n   = n;
  data.h   = h;

  /* From h, create the parameter vector.
   */
  hyperexp_pack_paramvector(p, np, h);

  /* Define the step size vector u.
   */
  for (i = 0; i < np; i++) u[i] = 1.0;

  /* Feed it all to the mighty optimizer.
   */
  status = esl_min_ConjugateGradientDescent(p, u, np, 
					    &hyperexp_complete_func, 
					    &hyperexp_complete_gradient,
					    (void *) (&data), tol, wrk, &fx);

  /* Convert the final parameter vector back to a hyperexponential
   */
  hyperexp_unpack_paramvector(p, np, h);
  
  free(p);
  free(u);
  free(wrk);
  return status;
}


/****************************************************************************
 * Maximum likelihood fitting, complete binned data         xref STL9/143-144
 ****************************************************************************/ 
/* minimizer API only allows us one generic void ptr to pass
 * our data through:
 */
struct hyperexp_binned_data {
  ESL_HISTOGRAM *g;	
  ESL_HYPEREXP  *h;
};
  
static double 
hyperexp_complete_binned_func(double *p, int np, void *dptr)
{
  struct hyperexp_binned_data *data = (struct hyperexp_binned_data *) dptr;
  ESL_HISTOGRAM               *g    = data->g;
  ESL_HYPEREXP                *h    = data->h;
  double logL = 0.;
  double ai, delta;
  int    i,k;

  hyperexp_unpack_paramvector(p, np, h);
  for (i = g->imin; i <= g->imax; i++) /* counting over occupied histogram bins */
    {
      if (g->obs[i] == 0) continue; /* skip unoccupied ones */

      esl_histogram_GetBinBounds(g, i, &ai, NULL, &delta);
      if (ai < h->mu) ai = h->mu; /* careful about the left boundary: no x < h->mu */

      for (k = 0; k < h->K; k++)
	{
	  h->wrk[k] = log(h->q[k]) - h->lambda[k]*(ai-h->mu);
	  if (delta * h->lambda[k] < eslSMALLX1) 
	    h->wrk[k] += log(delta * h->lambda[k]);
	  else
	    h->wrk[k] += log(1 - exp(-delta * h->lambda[k]));
	}
      logL += g->obs[i] * esl_vec_DLogSum(h->wrk, h->K);
    }
  return -logL;
}

static void
hyperexp_complete_binned_gradient(double *p, int np, void *dptr, double *dp)
{
  struct hyperexp_binned_data *data = (struct hyperexp_binned_data *) dptr;
  ESL_HISTOGRAM               *g    = data->g;
  ESL_HYPEREXP                *h    = data->h;
  int i,k;
  int pidx;			
  double z;
  double tmp;
  double ai, delta;
  
  hyperexp_unpack_paramvector(p, np, h);
  esl_vec_DSet(dp, np, 0.);
  for (i = g->imin; i <= g->imax; i++)
    {
      if (g->obs[i] == 0) continue;
      esl_histogram_GetBinBounds(g, i, &ai, NULL, &delta);
      if (ai < h->mu) ai = h->mu; /* careful about the left boundary: no x < h->mu */

      /* Calculate log (q_m alpha_m(a_i) terms
       */
      for (k = 0; k < h->K; k++)
	{
	  h->wrk[k] = log(h->q[k]) - h->lambda[k]*(ai-h->mu);
	  if (delta * h->lambda[k] < eslSMALLX1) 
	    h->wrk[k] += log(delta * h->lambda[k]);
	  else
	    h->wrk[k] += log(1 - exp(-delta * h->lambda[k]));
	}
      z = esl_vec_DLogSum(h->wrk, h->K); /* z= log \sum_k q_k alpha_k(a_i) */

      /* Bump the gradients for Q_1..Q_{K-1} */
      pidx = 0;
      for (k = 1; k < h->K; k++)
	dp[pidx++] -= g->obs[i] * (exp(h->wrk[k] - z) - h->q[k]);
	
      /* Bump the gradients for w_0..w_{K-1}
       */
      for (k = 0; k < h->K; k++)
	if (! h->fixlambda[k])
	  {
	    tmp  = log(h->q[k]) + log(h->lambda[k])- h->lambda[k]*(ai-h->mu);
	    tmp  = exp(tmp - z);
	    tmp *= (ai + delta - h->mu) * exp(-delta * h->lambda[k]) - (ai - h->mu);
	    dp[pidx++] -= g->obs[i] * tmp;
	  }
    }  
}

/* Function:  esl_hyperexp_FitGuessBinned()
 * Incept:    SRE, Mon Aug 15 14:02:02 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <g> with binned observations;
 *            obtain a very crude guesstimate of a fit -- suitable only 
 *            as a starting point for further optimization -- and return 
 *            those parameters in <h>.
 *
 *            Assigns $q_k \propto \frac{1/k}$ and  $\mu = \min_i x_i$;
 *            splits $x$ into $K$ roughly equal-sized bins, and
 *            and assigns $\lambda_k$ as the ML estimate from bin $k$.
 */
int
esl_hxp_FitGuessBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h)
{
  double sum;
  int    n;
  int    i,k;
  int    nb;
  double lowval;

  h->mu = g->xmin;
  nb    = g->imax - g->imin + 1;
  k     = h->K-1;
  sum   = 0;
  n     = 0;
  for (i = g->imax; i >= g->imin; i--)
    {
      lowval  = g->w*(double)i+g->bmin; /* low bound in this bin */
      if (lowval < g->xmin) lowval = g->xmin;
      n      += g->obs[i];
      sum    += g->obs[i] * lowval;
      
      if (i == g->imin + (k*nb)/h->K)
	h->lambda[k--] = 1 / ((sum/(double) n) - lowval);
    }

  for (k = 0; k < h->K; k++)
    h->q[k] = 1 / (double) h->K;

  return eslOK;
}


/* Function:  esl_hxp_FitCompleteBinned()
 * Incept:    SRE, Tue Aug 16 13:30:43 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <g> with binned observations, where each
 *            bin i holds some number of observed samples x with values from 
 *            lower bound l to upper bound u (that is, $l < x \leq u$),
 *            and given a starting guess <h> for hyperexponential parameters;
 *
 *            Find maximum likelihood parameters <h> by conjugate gradient
 *            descent, starting from the initial <h> and leaving the
 *            optimized solution in <h>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int 
esl_hxp_FitCompleteBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h)
{
  struct hyperexp_binned_data data;
  int     status;
  double *p;
  double *u;
  double *wrk;
  double  fx;
  int     i;
  double  tol = 1e-6;
  int     np;

  np = h->K-1;                /* K-1 mix coefficients...      */
  for (i = 0; i < h->K; i++)  /* ...and up to K lambdas free. */
    if (! h->fixlambda[i]) np++;

  p   = malloc(sizeof(double) * np);
  u   = malloc(sizeof(double) * np);
  wrk = malloc(sizeof(double) * np * 4);

  /* Copy shared info into the "data" structure  */
  data.g     = g;
  data.h     = h;

  /* From h, create the parameter vector. */
  hyperexp_pack_paramvector(p, np, h);

  /* Define the step size vector u.
   */
  for (i = 0; i < np; i++) u[i] = 1.0;

  /* Feed it all to the mighty optimizer.
   */
  status = esl_min_ConjugateGradientDescent(p, u, np, 
					    &hyperexp_complete_binned_func, 
					    &hyperexp_complete_binned_gradient,
					    (void *) (&data), tol, wrk, &fx);

  /* Convert the final parameter vector back to a hyperexponential
   */
  hyperexp_unpack_paramvector(p, np, h);
  
  free(p);
  free(u);
  free(wrk);
  return status;
}
#endif /*eslAUGMENT_MINIMIZER*/


/****************************************************************************
 * Example main()
 ****************************************************************************/ 

#ifdef eslHYPEREXP_EXAMPLE
/*::cexcerpt::hyperexp_example::begin::*/
/* compile: 
   gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o example -DeslHYPEREXP_EXAMPLE\
     esl_hyperexp.c esl_exponential.c eslesl_vectorops.c esl_histogram.c \
     esl_random.c esl_minimizer.c easel.c -lm

   gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o example -DeslHYPEREXP_EXAMPLE\
     esl_hyperexp.c -leasel -lm
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <esl_hyperexp.h>
#include <esl_histogram.h>
#include <esl_random.h>

int
main(int argc, char **argv)
{
  FILE *fp;
  ESL_RANDOMNESS *r;		/* source of random numbers        */
  ESL_HISTOGRAM  *h;		/* histogram to store the data     */
  ESL_HYPEREXP *hxp;		/* hyperexponential to sample from */
  ESL_HYPEREXP *ehxp;		/* estimated hyperexponential      */
  double      x;		/* sampled data point              */
  int         n = 100000;	/* number of samples               */
  int         i;
  int         k;
  double      nll;
  double      min, max;

  fp = fopen("data.xy", "w");

  r   = esl_randomness_CreateTimeseeded();
  h   = esl_histogram_CreateFull(-3, 100, 1.0);
  hxp = esl_hyperexp_Create(3);
  hxp->mu = -2.0;
  hxp->q[0]      = 0.6;    hxp->q[1]      = 0.3;   hxp->q[2]      = 0.1; 
  hxp->lambda[0] = 1.0;    hxp->lambda[1] = 0.3;   hxp->lambda[2] = 0.1;

  nll = 0.;
  for (i = 0; i < n; i++)
    {
      x    = esl_hxp_Sample(r, hxp);
      nll -= esl_hxp_logpdf(x, hxp);
      esl_histogram_Add(h, x);
    }
  esl_histogram_Finish(h);
  min = h->x[0];
  max = h->x[n-1];
  printf("NLL of known hyperexp: %g\n", nll);

  esl_histogram_PlotSurvival(fp, h);
  esl_hxp_Plot(fp, hxp, &esl_hxp_surv, hxp->mu, max+5, 0.1);

  ehxp = esl_hyperexp_Create(3);
  esl_hxp_FitGuess(h->x, n, ehxp);
  /* esl_hyperexp_Copy(hxp, ehxp); */
  esl_hxp_Plot(fp, ehxp, &esl_hxp_surv, hxp->mu, max+5, 0.1);
  printf("Guessed:\n");
  printf("Component   q      lambda\n");
  for (k=0; k < 3; k++)
    printf("%d\t%7.4f\t%7.4f\n",
	   k, ehxp->q[k], ehxp->lambda[k]);
  printf("and mu = %f\n", ehxp->mu);
  nll = 0.;
  for (i = 0; i < n; i++)
    nll -= esl_hxp_logpdf(h->x[i], ehxp);
  printf("NLL of guessed fit: %g\n", nll);

  esl_hxp_FitComplete(h->x, n, ehxp);
  esl_hxp_Plot(fp, ehxp, &esl_hxp_surv, hxp->mu, max+5, 0.1);
  printf("Optimized:\n");
  printf("Component   q      lambda\n");
  for (k=0; k < 3; k++)
    printf("%d\t%7.4f\t%7.4f\n",
	   k, ehxp->q[k], ehxp->lambda[k]);
  printf("and mu = %f\n", ehxp->mu);
  nll = 0.;
  for (i = 0; i < n; i++)
    nll -= esl_hxp_logpdf(h->x[i], ehxp);
  printf("NLL of optimized fit: %g\n", nll);

  fclose(fp);
  esl_hyperexp_Destroy(hxp);
  esl_hyperexp_Destroy(ehxp);
  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}

/*::cexcerpt::hyperexp_example::end::*/
#endif /*eslHYPEREXP_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
