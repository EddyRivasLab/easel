/* Finding roots.
 * 
 * Contents:
 *   1. The ESL_ROOTFINDER object.
 *   2. One-dimensional root finding.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Examples.
 *   6. Copyright and license information.
 * 
 * SRE, Fri Apr  6 09:14:13 2007 [Janelia]
 * SVN $Id$
 */

#include <math.h>

#include <easel.h>
#include <esl_rootfinder.h>

/*****************************************************************
 * 1. The ESL_ROOTFINDER object.
 *****************************************************************/

/* Function:  esl_rootfinder_CreateBracketer()
 * Synopsis:  Creates ESL_ROOTFINDER for bisection method.
 * Incept:    SRE, Tue Apr 10 19:54:09 2007 [Janelia]
 *
 * Purpose:   Create a bisection rootfinder that will find a value x in
 *            the open interval <xl..xr> such that $f(x) = 0$, where a
 *            pointer to the function $f(x)$ is passed as <func()>.
 *            
 *            Caller provides a <func()> that takes two arguments:
 *            <x>, and a void pointer to any additional parameters that
 *            $f(x)$ depends on; it returns $f(x)$. 
 *            
 *            Caller is responsible for being sure that a root
 *            actually exists in the open interval <xl..xr>.
 *            
 *            The bisection method may be slow, but it is guaranteed
 *            to converge to the root, provided that the brackets
 *            enclose one.
 *
 * Args:      (*func)() - ptr to function that evaluates f(x)
 *            params    - ptr to parameters to be passed to (*func)()
 *            xl,xr     - interval to seek root x in
 *
 * Returns:   pointer to a new ESL_ROOTFINDER structure.
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINVAL> if <xl..xr> cannot bracket a root because $f(x_l)$
 *            $f(x_r)$ have the same sign.
 */
ESL_ROOTFINDER *
esl_rootfinder_CreateBracketer(double (*func)(double, void*), void *params, double xl, double xr)
{
  int status;
  ESL_ROOTFINDER *R = NULL;

  ESL_ALLOC(R, sizeof(ESL_ROOTFINDER));
  R->func          = func;
  R->fdf           = NULL;
  R->params        = params;
  R->xl            = xl;
  R->fl            = (*func)(xl, params);
  R->xr            = xr;
  R->fr            = (*func)(xr, params);
  R->x0            = 0.;	
  R->f0            = 0.;	
  R->x             = 0.;	/* not set yet */
  R->fx            = 0.;	/* not set yet */
  R->dfx           = 0.;	/* unused */
  R->iter          = 0;
  R->abs_tolerance = 1e-15;
  R->rel_tolerance = 1e-15;
  R->residual_tol  = 0.;
  R->max_iter      = 100; 
  
  if (R->fl * R->fr > 0) ESL_XEXCEPTION(eslEINVAL, "xl, xr don't bracket the root");

  return R;

 ERROR:
  esl_rootfinder_Destroy(R);
  return NULL;
}


/* Function:  esl_rootfinder_CreatePolisher()
 * Synopsis:  Creates ESL_ROOTFINDER for Newton/Raphson method.
 * Incept:    SRE, Tue Apr 10 20:47:42 2007 [Janelia]
 *
 * Purpose:   Create a Newton/Raphson rootfinder that will find 
 *            a root of a function $f(x) = 0$ using first derivative
 *            information and an initial guess at the root $x$. 
 *            
 *            Caller provides a pointer <*fdf()> to a function that
 *            takes four arguments. The first two are the current <x>
 *            value, and a void pointer to any additional parameters
 *            that $f(x)$ depends on. <*fdf()> calculates the function
 *            $f(x)$ and the derivative $f'(x)$ and returns them
 *            through the remaining two arguments.
 *            
 *            The Newton/Raphson algorithm is not guaranteed to
 *            succeed. Caller needs to worry about this, and should
 *            provide an initial guess for $x$ that is "suitably
 *            close" to the root: hence the name "polisher" for this
 *            sort of root finding.
 *
 * Args:      (*fdf)() - ptr to function that returns f(x) and f'(x)
 *            params   - ptr to parameters to be passed to (*fdf)()
 *            guess    - initial guess 
 *
 * Returns:   
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      
 */
ESL_ROOTFINDER *
esl_rootfinder_CreatePolisher(void (*fdf)(double, void*, double*, double*), void *params, double guess)
{
  int status;
  ESL_ROOTFINDER *R = NULL;

  ESL_ALLOC(R, sizeof(ESL_ROOTFINDER));
  R->func          = NULL;
  R->fdf           = fdf;
  R->params        = params;
  R->xl            = -eslINFINITY;
  R->fl            = 0.;	/* unused */
  R->xr            = eslINFINITY;
  R->fr            = 0.;	/* unused */
  R->x0            = 0.;	
  R->f0            = 0.;	
  R->x             = guess;
  R->iter          = 0;
  R->abs_tolerance = 1e-15;
  R->rel_tolerance = 1e-15;
  R->residual_tol  = 0.;
  R->max_iter      = 100; 
  (*fdf)(R->x, params, &(R->fx), &(R->dfx));

  return R;

 ERROR:
  esl_rootfinder_Destroy(R);
  return NULL;
}

int
esl_rootfinder_SetBrackets(ESL_ROOTFINDER *R, double xl, double xr)
{
  R->xl = xl;
  R->xr = xr;
  return eslOK;
}

int
esl_rootfinder_SetAbsoluteTolerance(ESL_ROOTFINDER *R, double tol)
{
  R->abs_tolerance = tol;
  return eslOK;
}

int
esl_rootfinder_SetRelativeTolerance(ESL_ROOTFINDER *R, double tol)
{
  R->rel_tolerance = tol;
  return eslOK;
}

int
esl_rootfinder_SetResidualTolerance(ESL_ROOTFINDER *R, double tol)
{
  R->residual_tol = tol;
  return eslOK;
}

int
esl_rootfinder_SetMaxIterations(ESL_ROOTFINDER *R, int maxiter)
{
  R->max_iter = maxiter;
  return eslOK;
}


void
esl_rootfinder_Destroy(ESL_ROOTFINDER *R)
{
  if (R == NULL) return;
  free(R);
}


/*****************************************************************
 * 2. One-dimensional root finding.
 *****************************************************************/

int
esl_root_Bisection(ESL_ROOTFINDER *R, double *ret_x)
{
  double xmag;

  while (1) {
    R->iter++;
    if (R->iter > R->max_iter) ESL_EXCEPTION(eslENOHALT, "failed to converge in Bisection");

    /* Bisect and evaluate the function */
    R->x  = (R->xl+R->xr)/2.; 	          
    R->fx = (*R->func)(R->x, R->params);

    /* Test for convergence */
    xmag = (R->xl < 0. && R->xr > 0.) ?  0. : R->x;
    if (((R->xr-R->xl)  <  R->abs_tolerance + R->rel_tolerance*xmag) || fabs(R->fx) < R->residual_tol) break;

    /* Narrow the bracket; pay attention to directionality */
    if (R->fl > 0.) {
      if   (R->fx > 0.) { R->xl = R->x; R->fl = R->fx; }
      else              { R->xr = R->x; R->fr = R->fx; }
    } else {
      if   (R->fx < 0.) { R->xl = R->x; R->fl = R->fx; }
      else              { R->xr = R->x; R->fr = R->fx; }      
    }
  }
  
  *ret_x = R->x;
  return eslOK;
}

int
esl_root_Newton(ESL_ROOTFINDER *R, double *ret_x)
{
  while (1) {
    R->iter++;
    if (R->iter > R->max_iter) ESL_EXCEPTION(eslENOHALT, "failed to converge in Newton");

    printf("current: x=%20g   f(x) = %20g   f'(x) = %20g\n", R->x, R->fx, R->dfx);

    /* Take a Newton/Raphson step. */
    R->x0  = R->x;
    R->f0  = R->fx;
    R->x   = R->x - R->fx / R->dfx;
    (*R->fdf)(R->x, R->params, &(R->fx), &(R->dfx));  

    /* Test for convergence. */
    if ( (fabs(R->x - R->x0) < R->abs_tolerance + R->rel_tolerance*R->x) || fabs(R->fx < R->residual_tol)) break;
  }

  *ret_x = R->x;
  return eslOK;
}


/*****************************************************************
 * Examples.
 *****************************************************************/

/* An example of bisection.
 *   gcc -g -Wall -o example -I. -DeslROOTFINDER_EXAMPLE esl_rootfinder.c easel.c -lm
 */
#ifdef eslROOTFINDER_EXAMPLE
/*::cexcerpt::rootfinder_example::begin::*/
#include <easel.h>
#include <esl_rootfinder.h>

struct polyparams { double a,b,c; };

double quadratic_f(double x, void *params)
{
  struct polyparams *p = (struct polyparams *) params;
  return (p->a * x * x + p->b * x + p->c);
}

int main(void)
{
  ESL_ROOTFINDER *R = NULL;
  struct polyparams p;
  double x;

  p.a = 5.;
  p.b = 2.;
  p.c = -1.;

  R = esl_rootfinder_CreateBracketer(quadratic_f, &p, 0., 100.);
  esl_root_Bisection(R, &x);

  printf("Find an x such that f(x) = %.0fx^2 + %.0fx + %.0f = 0 ...\n", p.a, p.b, p.c);
  printf("x = %f (f(x) = %f)\n", x, quadratic_f(x, &p));

  esl_rootfinder_Destroy(R);
  return 0;
}
/*::cexcerpt::rootfinder_example::end::*/
#endif /*eslROOTFINDER_EXAMPLE*/


/* An example of Newton/Raphson.
 *   gcc -g -Wall -o example -I. -DeslROOTFINDER_EXAMPLE2 esl_rootfinder.c easel.c -lm
 */
#ifdef eslROOTFINDER_EXAMPLE2
/*::cexcerpt::rootfinder_example2::begin::*/
#include <easel.h>
#include <esl_rootfinder.h>

struct polyparams { double a,b,c; };

void quadratic_fdf(double x, void *params, double *ret_fx, double *ret_dfx)
{
  struct polyparams *p = (struct polyparams *) params;
  
  *ret_fx  = (p->a * x * x + p->b * x + p->c);
  *ret_dfx =  (2 * p->a) * x + p->b;
  return;
}

int main(void)
{
  ESL_ROOTFINDER *R = NULL;
  struct polyparams p;
  double x;

  p.a = 5.;
  p.b = 2.;
  p.c = -1.;

  R = esl_rootfinder_CreatePolisher(quadratic_fdf, &p, -1);
  esl_root_Newton(R, &x);

  printf("Find an x such that f(x) = %.0fx^2 + %.0fx + %.0f = 0 ...\n", p.a, p.b, p.c);
  printf("x = %f\n", x);

  esl_rootfinder_Destroy(R);
  return 0;
}
/*::cexcerpt::rootfinder_example2::end::*/
#endif /*eslROOTFINDER_EXAMPLE2*/




/*****************************************************************
 * @LICENSE@
 *****************************************************************/
