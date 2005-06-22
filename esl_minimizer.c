/* esl_minimizer.c
 * Multidimensional optimization using conjugate gradient descent.
 * 
 * SRE, Wed Jun 22 11:37:29 2005
 * SVN $Id$
 */

#include <math.h>
#include <float.h>

#include <easel.h>
#include <esl_vectorops.h>
#include <esl_minimizer.h>

/* Function:  esl_min_Bracket()
 * Incept:    SRE, Mon Jun 20 16:47:51 2005 [St. Louis]
 *
 * Purpose:   Brackets a minimum: starting at initial vector <a>
 *            and minimizing along the gradient <d>, for vectors
 *            of length <n>, determines two additional vectors <b>
 *            and <c> on that line such that $b = a + bx * d$,
 *            $c = a + cx * d$, $a < b < c$, and $f(b> < f(a),f(c)$.
 *            
 *            Caller passes an objective function <*func()>,
 *            and a void pointer to any necessary conditional 
 *            parameters <prm>. The objective function will
 *            be evaluated at a point <x> by calling
 *            <(*func)(x, n, prm)>. The caller's function
 *            is responsible to casting <prm> to whatever it's
 *            supposed to be, which might be a ptr to a structure,
 *            for example.
 *         
 *            The magnitude of the initial step along the gradient
 *            is determined by <xinit>: <xinit>*<d>. In principle, this
 *            can be any positive number, and the bracketing will work;
 *            but if the caller has some rough knowledge about where
 *            the minimum might lie, it can set <xinit>*<d> to try
 *            to jump there, perhaps accelerating bracketing and
 *            subsequent line minimization.
 *
 * Args:      a         - initial point; an n-vector [0..n-1]
 *            d         - gradient evaluated at a; an n-vector [0..n-1]
 *            n         - size of all vectors
 *            (*func)() - objective function
 *            prm       - void ptr to any conditional info that *func needs
 *            xinit     - initial step along gradient is xinit*d; [xinit>0]
 *            ret_fa    - optRETURN: function evaluation at a
 *            b         - RETURN: midpoint; [0..n-1], caller-provided space 
 *            ret_bx    - optRETURN: scalar that determines b = a + bx*d
 *            ret_fb    - optRETURN: function evaluation at b
 *            c         - RETURN: farpoint; [0..n-1], caller-provided space
 *            ret_cx    - optRETURN: scalar that determines c = a + cx*d
 *            ret_fc    - optRETURN: function evaluation at c
 *
 * Returns:   <eslOK> on success; and
 *            returns scalars <ret_bx> and <ret_cx>, and sets the vectors
 *            <b> and <c> they determine (caller provides allocated
 *            space for these vectors). Also returns the function
 *            evaluations at all three points in <ret_fa>, <ret_fb>,
 *            and <ret_fc>.
 *            
 * Throws:    <eslECONVERGENCE> if bracketing is not achieved in a certain
 *            number of iterations (set by a compile-time constant,
 *            <MAXITERATIONS>).
 *            
 * Xref:      STL9/109.
 */
int
esl_min_Bracket(double *a, double *d, int n, 
		double (*func)(double *, int, void *), void *prm, 
		double xinit, double *ret_fa,
		double *b, double *ret_bx, double *ret_fb,
		double *c, double *ret_cx, double *ret_fc)
{
  double fa, fb, fc;		/* f(a), f(b), f(c); f(b) < f(a), f(c) */
  double bx, cx;		/* scalar multipliers that give the b,c vectors */
  int    i;			/* counter over iterations */

  fa = (*func)(a, n, prm);

  /* Reach out by xinit. For now, set both b, c vectors the same here.
   * If f(c) >= f(a), then we'll look for minimum b in between them;
   *                  since d was supposed to be a gradient, we think
   *                  we're guaranteed that such b must exist.
   * If f(b) < f(a),  then it's a candidate for b, and we'll look for c 
   *                  further out.
   */ 
  bx = cx = xinit;
  esl_vec_DCopy(b, a, n); 
  esl_vec_DAddScaled(b, d, bx, n);
  esl_vec_DCopy(c, b, n);
  fb = fc = (*func)(b, n, prm);

  /* Now, from there, iteratively search for bracketing triplet a,b,c.
   */
  if (fc >= fa) /* case 1: where we have c and we looking in between for b... */
    {
      for (i = 0; i < MAXITERATIONS; i++)
	{
	  if (i > 0) 	/* make prev b new c. unnecessary on first iteration */
	    { esl_vec_DCopy(c, b, n); fc = fb; cx = bx; }
	  bx /= eslCONST_GOLD;	        /* maintain the golden ratio */
	  esl_vec_DCopy(b, a, n); 
	  esl_vec_DAddScaled(b, d, bx, n); /* find new b */
	  fb = (*func)(b, n, prm); 
	  if (fb < fa) break;
	} 
      if (i == MAXITERATIONS) 
	ESL_ERROR(eslECONVERGENCE, "Failed to bracket a minimum.");
    }
  else 		/* case 2: where we have b, and we look further out for c... */
    {
      for (i = 0; i < MAXITERATIONS; i++)
	{
	  if (i > 0)   /* make prev c new b. wasted on iter 1. */
	    { esl_vec_DCopy(b, c, n);  fb = fc; bx = cx; }
	  cx *= eslCONST_GOLD;     /* maintain the golden ratio */
	  esl_vec_DCopy(c, a, n);
	  esl_vec_DAddScaled(c, d, cx, n); /* find new c */
	  fc = (*func)(c, n, prm); 
	  if (fc > fb) break;
	}
      if (i == MAXITERATIONS) 
	ESL_ERROR(eslECONVERGENCE, "Failed to bracket a minimum.");
    }

  /* a,b,c vectors now bracket a minimum along gradient d,
   * where b = a + bx*d, c = a + cx*d, and f(b) < f(a), f(c);
   * set up our optional returned data, and return.
   */
  if (ret_bx != NULL) *ret_bx = bx;
  if (ret_cx != NULL) *ret_cx = cx;
  if (ret_fa != NULL) *ret_fa = fa;
  if (ret_fb != NULL) *ret_fb = fb;
  if (ret_fc != NULL) *ret_fc = fc;
  return eslOK;
}


/* Function:  esl_min_LineSearch()
 * Incept:    SRE, Tue Jun 21 15:31:30 2005 [St. Louis]
 *
 * Purpose:   Minimization along a gradient in n-dimensional space by 
 *            golden section search.
 *            
 *            We know there is a minimum on a line starting
 *            at <ori> in the direction of the gradient <d>,
 *            where these are vectors of dimension <n>.
 *            
 *            We provide the objective function <*func()>, and a void
 *            pointer to any necessary data and conditional parameters
 *            <prm>. The objective function will be evaluated at
 *            vectors <x> by calling <(*func)(x, n, prm)>. The
 *            provided function is responsible to casting <prm> to
 *            whatever it's supposed to be, which might be a ptr to a
 *            structure or a data vector, for example.
 *            
 *            We also provide an allocated <n>-vector <b> as
 *            temporary work space for the routine.
 *            
 *            The routine executes an iterative golden section
 *            search along this line, and narrows in on a
 *            suitably infinestimal interval in which the
 *            minimum must lie.
 *            
 *            It then returns <x>, the vector at the minimum (caller
 *            provides this allocated memory); returns <ret_xx),
 *            the scalar multiplier used to find <x> from <ori>
 *            and the gradient; and <ret_fx>, the value of the
 *            objective function at the minimum <x>.
 *            
 *            <ori> and <d> are unchanged by this procedure.
 *            <b> is workspace; it will contain numbers upon return,
 *            but don't use them for anything.
 *            
 * Args:      ori	- n-vector at origin
 *            d         - gradient from ori
 *            n         - dimensionality of all vectors
 *            (*func)   - pointer to caller's objective function
 *            prm       - void pointer to any data that (*func) needs
 *            b         - n-vector to be used for workspace
 *            x         - RETURN: vector at the minimum
 *            ret_xx    - optRETURN: scalar multiplier that gave vector x
 *            ret_fx    - optRETURN: f(x)
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      STL9/110.
 */
int
esl_min_LineSearch(double *ori, double *d, int n,
		   double (*func)(double *, int, void *), void *prm,
		   double *b, 
		   double *x, double *ret_xx, double *ret_fx)
{
  double xx, swap;
  double fa, fb, fx;
  double tolerance;
  double ax, bx, cx;

  /* rationale for sqrt() below is a Taylor expansion around b; 
   * xref Numerical Recipes */
  tolerance = sqrt(DBL_EPSILON);
  
  /* Bracket the minimum along line d first; this gives us scalar
   * points $0 <= ax < bx < cx$ relative to vector d. That is, vectors
   * at these points are determined as, for example, $a = ori + ax*d$.
   */
  ax = 0.;
  esl_min_Bracket(ori, d, n, func, prm, 1.0, &fa, b, &bx, &fb, x, &cx, &fx);

  /* the main do loop is guaranteed to terminate; 
   * no need to count iterations 
   */
  do {	
    if ((bx - ax) > (cx - bx))	/* find the larger interval */
      {				/* a..b larger: set up a..x.b..c */
	xx = (bx + ax*eslCONST_GOLD) / (1.+eslCONST_GOLD);
	swap = xx; xx = bx; bx = swap; /* now a..b.x..c */
      }
    else	/* b..c larger: set up a..b.x..c */
      xx = (cx + bx*eslCONST_GOLD) / (1.+eslCONST_GOLD);
      
    /* Calculate new vectors at b and x, along line d from ori. */
    esl_vec_DCopy(b, ori, n);  esl_vec_DAddScaled(b, d, bx, n);
    esl_vec_DCopy(x, ori, n);  esl_vec_DAddScaled(x, d, xx, n);
  
    /* Calculate new function values at those points */
    fb = (*func)(b, n, prm);
    fx = (*func)(x, n, prm);

    if (fb < fx) 
      { /* then a..b.x is the new bracketing; discard c */ 				
	cx = xx;
	fx = fb;
      }
    else /* else fx <= fb; discard a; b.x..c is the new bracketing */
      {
	ax = bx;
	bx = xx;
      }
    /* now we have a new a.b.c bracket; in principle we also
     * know which subinterval is larger, but it's not expensive
     * to figure it out again when we loop now. */
  } while ((cx-ax)/bx > tolerance); 

  /* Minimum is at b. 
   * Make sure we have x vector there: wasteful 50% of the time */
  esl_vec_DCopy(x, ori, n); esl_vec_DAddScaled(x, d, bx, n); 
  if (ret_xx != NULL) *ret_xx = bx;
  if (ret_fx != NULL) *ret_fx = fx;
  return eslOK;
}

		   

/* Function:  esl_min_ConjugateGradientDescent()
 * Incept:    SRE, Wed Jun 22 08:49:42 2005 [St. Louis]
 *
 * Purpose:   n-dimensional minimization by conjugate gradient descent.
 *           
 *            An initial point is provided by <x>, a vector of <n>
 *            components. The caller also provides a function <*func()> that 
 *            compute the objective function f(x) when called as 
 *            <(*func)(x, n, prm)>, and a function <*dfunc()> that can
 *            compute the gradient <dx> at <x> when called as 
 *            <(*dfunc)(x, n, prm, dx)>, given an allocated vector <dx>
 *            to put the derivative in. Any additional data or fixed
 *            parameters that these functions require are passed by
 *            the void pointer <prm>.
 *            
 *            The caller also provides four allocated n-vectors as
 *            workspace.  It is not expected that the contents of these
 *            vectors will be used upon return. However, if you
 *            care: upon return, <dx> contains the gradient at the
 *            minimum (which ought to be close to zero), and <cg>
 *            contains the conjugate direction that would've been
 *            followed next if we hadn't terminated the
 *            minimization. <w1> and <w2> have no particular meaning.
 *            
 *            Upon return, <x> is the minimum, and <ret_fx> is f(x),
 *            the function value at <x>.
 *            
 * Args:      x        - an initial guess n-vector; RETURN: x at the minimum
 *            n        - dimensionality of all vectors
 *            *func()  - function for computing objective function f(x)
 *            *dfunc() - function for computing a gradient at x
 *            prm      - void ptr to any data/params func,dfunc need 
 *            dx       - allocated n-vector for workspace; current gradient
 *            cg       - allocated n-vector for workspace; current direction
 *            w1       - allocated n-vector for workspace
 *            w2       - allocated n-vector for workspace
 *            ret_fx   - optRETURN: f(x) at the minimum
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECONVERGENCE> if it fails to converge in MAXITERATIONS.
 *
 * Xref:      STL9/101.
 */
int
esl_min_ConjugateGradientDescent(double *x, int n, 
       				 double (*func)(double *, int, void *),
				 void (*dfunc)(double *, int, void *, double *),
				 void *prm,
				 double *dx, double *cg, double *w1, double *w2,
				 double *ret_fx)
{
  double oldfx;
  double coeff;
  int    i, i1;
  double tolerance;

  oldfx = (*func)(x, n, prm);	/* init the objective function */
  (*dfunc)(x, n, prm, dx);	/* find the current negative gradient, - df(x)/dxi  */
  esl_vec_DScale(dx, n, -1.0);
  esl_vec_DCopy(cg, dx, n);	/* and make that the first conjugate direction, cg  */

  tolerance = DBL_EPSILON;

  for (i = 0; i < MAXITERATIONS; i++)
    {
      /* Minimize along the line given by the conjugate gradient <cg> */
      esl_min_LineSearch(x, cg, n, func, prm, w1, w2, NULL, ret_fx);
      esl_vec_DCopy(x, w2, n);

      /* Find the negative gradient at that point (temporarily in w1) */
      (*dfunc)(x, n, prm, w1);
      esl_vec_DScale(w1, n, -1.0);

      /* Calculate the Polak-Ribiere coefficient */
      for (coeff = 0., i1 = 0; i1 < n; i1++)
	coeff += (w1[i1] - dx[i1]) * w1[i1];
      coeff /= esl_vec_DDot(dx, dx, n);
      
      /* Calculate the next conjugate gradient direction in w2 */
      esl_vec_DCopy(w2, w1, n);
      esl_vec_DAddScaled(w2, cg, coeff, n);

      /* Finishing set up for next iteration: */
      esl_vec_DCopy(dx, w1, n);
      esl_vec_DCopy(cg, w2, n);

      /* Now: x is the current point; 
       *      *ret_fx is the function value at that point;
       *      dx is the current gradient at x;
       *      cg is the current conjugate gradient direction. 
       */

      /* Convergence test.
       */
      if (2.0 * fabs((oldfx-*ret_fx)) / (1e-9 + fabs(oldfx) + fabs(*ret_fx)) <= tolerance) break;
      oldfx = *ret_fx;
    }
  if (i == MAXITERATIONS) 
    ESL_ERROR(eslECONVERGENCE, "Failed to converge in ConjugateGradientDescent()");

  return eslOK;
}





/*****************************************************************
 * Example main()
 *****************************************************************/
#ifdef eslMINIMIZER_EXAMPLE
/*::cexcerpt::minimizer_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslMINIMIZER_EXAMPLE esl_minimizer.c esl_vectorops.c easel.c -lm
 * run:     ./example 
 */
#include <stdio.h>

#include <easel.h>
#include <esl_vectorops.h>
#include <esl_minimizer.h>

/* a simple multidimensional quadratic w/ a minimum at 0:
 *    $f(x) = a_1 x_1^2 + ... a_n x_n^2$
 */ 
static double
example_func(double *x, int n, void *prm)
{
  double *a;
  double  fx;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (fx = 0., i = 0; i < n; i++)
    fx += a[i] * x[i] * x[i];
  return fx;
}
/* gradient of the f(x): d/dx_i = 2 a_i x_i
 */
static void
example_dfunc(double *x, int n, void *prm, double *dx)
{
  double *a;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (i = 0; i < n; i++)
    dx[i] = 2.0 * a[i] * x[i];
}
int
main(int argc, char **argv)
{
  int    n = 6;
  double a[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double x[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double v1[6], v2[6], v3[6], v4[6]; 
  double fx;
  int    i;

  esl_min_ConjugateGradientDescent(x, n, 
				   &example_func, &example_dfunc, (void *) a, 
				   v1, v2, v3, v4,
				   &fx);

  printf("At minimum: f(x) = %g\n", fx);
  printf("vector x = ");
  for (i = 0; i < 6; i++) printf("%g  ", x[i]);
  printf("\n");

  return 0;
}
/*::cexcerpt::minimizer_example::end::*/
#endif /*eslMINIMIZER_EXAMPLE*/






/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
