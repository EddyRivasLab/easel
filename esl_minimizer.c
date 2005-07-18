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
 *            The initial step is key. It shouldn't be too small (or
 *            we'll spend a long time finding the other side where it
 *            increases), and definitely shouldn't be too large (else
 *            we may even blow up the objective function). The caller
 *            provides us with <u>, which indicates the "natural step
 *            size" in each dimension, as if we were marking the axes
 *            of a graph of f(). The initial step cannot exceed u[i]
 *            for any dimension [i].
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
 * Args:      a         - initial point; an n-vector [0..n-1]
 *            d         - gradient evaluated at a; an n-vector [0..n-1]
 *            u         - "units": max initial step size for any dimension i
 *            n         - size of all vectors
 *            (*func)() - objective function
 *            prm       - void ptr to any conditional info that *func needs
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
esl_min_Bracket(double *a, double *d, double *u, int n, 
		double (*func)(double *, int, void *), void *prm, 
		double *ret_fa,
		double *b, double *ret_bx, double *ret_fb,
		double *c, double *ret_cx, double *ret_fc)
{
  double fa, fb, fc;		/* f(a), f(b), f(c); f(b) < f(a), f(c) */
  double bx, cx;		/* scalar multipliers that give the b,c vectors */
  int    i;			/* counter over iterations */

  fa = (*func)(a, n, prm);

  /* Figure out the initial step size.
   */
  bx = fabs(u[0] / d[0]);
  for (i = 1; i < n; i++)
    {
      cx = fabs(u[i] / d[i]);
      if (cx < bx) bx = cx;
    }
  bx = 0.01; /*SRE!*/

  /* Reach out by that step. For now, set both b, c vectors the same here.
   * If f(c) >= f(a), then we'll look for minimum b in between them;
   *                  since d was supposed to be a gradient, we think
   *                  we're guaranteed that such b must exist.
   * If f(b) < f(a),  then it's a candidate for b, and we'll look for c 
   *                  further out.
   */ 
  cx = bx;
  esl_vec_DCopy(b, a, n); 
  esl_vec_DAddScaled(b, d, bx, n);
  esl_vec_DCopy(c, b, n);
  fb = fc = (*func)(b, n, prm);

  /* Now, from there, iteratively search for bracketing triplet a,b,c.
   */
  if (fc >= fa) /* case 1: where we have c and we look in between for b... */
    {
      for (i = 0; i < MAXITERATIONS; i++)
	{
	  if (i > 0) 	/* make prev b new c. unnecessary on first iteration */
	    { esl_vec_DCopy(c, b, n); fc = fb; cx = bx; }
	  bx /= eslCONST_GOLD;
	  esl_vec_DCopy(b, a, n); 
	  esl_vec_DAddScaled(b, d, bx, n); /* find new b */
	  fb = (*func)(b, n, prm); 
	  if (fb <= fa) break;
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
	  cx *= eslCONST_GOLD; 
	  esl_vec_DCopy(c, a, n);
	  esl_vec_DAddScaled(c, d, cx, n); /* find new c */
	  fc = (*func)(c, n, prm); 
	  if (fc >= fb) break;
	}
      if (i == MAXITERATIONS) 
	ESL_ERROR(eslECONVERGENCE, "Failed to bracket a minimum.");
    }

  ESL_DPRINTF2(("\nesl_min_Bracket(): %d iterations\n", i));
  ESL_DPRINTF2(("esl_min_Bracket(): triplet is %g  %g  %g along current direction\n", 0.0, bx, cx));
  ESL_DPRINTF2(("esl_min_Bracket(): f()'s there are: %g  %g  %g\n\n", fa, fb, fc));

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

#if eslDEBUGLEVEL >= 2
static void
show_line(double *ori, double *d, double *u, double n,
	  double ax, double cx,
	  double (*func)(double *, int, void *), void *prm, double *wrk)
{
  double bx, xx, fx;
  int i;

  bx = fabs(u[0] / d[0]);
  for (i = 1; i < n; i++) {
    xx = fabs(u[i] / d[i]);
    if (xx < bx) bx = xx;
  }
  for (xx = ax; xx <= cx; xx += bx)
    {
      esl_vec_DCopy(wrk, ori, n); 
      esl_vec_DAddScaled(wrk, d, xx, n); 
      fx = (*func)(wrk, n, prm);
      printf("%g  %g\n", xx, fx);
    }
}
#endif



/* brent():
 * SRE, Sun Jul 10 19:07:05 2005 [St. Louis]
 *
 * Purpose:   Quasi-one-dimensional minimization of a function <*func()>
 *            in <n>-dimensions, along vector <dir> starting from a
 *            point <ori>. Identifies a scalar $x$ that approximates
 *            the position of the minimum along this direction, in a
 *            given bracketing interval (<a,b>).  The minimum must
 *            have been bracketed by the caller in the <(a,b)>
 *            interval.  <a> is often 0, because we often start at the
 *            <ori>.
 *
 *            A quasi-1D scalar coordinate $x$ (such as <a> or <b>) is
 *            transformed to a point $\mathbf{p}$ in n-space as:
 *            $\mathbf{p} = \mathbf{\mbox{ori}} + x
 *            \mathbf{\mbox{dir}}$.
 *
 *            Any extra (fixed) data needed to calculate <func> can be
 *            passed through the void <prm> pointer.
 *
 *            <eps> and <t> define the relative convergence tolerance,
 *            $\mbox{tol} = \mbox{eps} |x| + t$. <eps> should not be
 *            less than the square root of the machine precision.  The
 *            <DBL_EPSILON> is 2.2e-16 on many machines with 64-bit
 *            doubles, so <eps> is on the order of 1e-8 or more. <t>
 *            is a yet smaller number, used to avoid nonconvergence in
 *            the pathological case $x=0$.
 *
 *            Upon convergence (which is guaranteed), returns <xvec>,
 *            the n-dimensional minimum. Optionally, will also return
 *            <ret_x>, the scalar <x> that resulted in that
 *            n-dimensional minimum, and <ret_fx>, the objective
 *            function <*func(x)> at the minimum.
 *
 *            This is an implementation of the R.P. Brent (1973)
 *            algorithm for one-dimensional minimization without
 *            derivatives (modified from Brent's ALGOL60 code). Uses a
 *            combination of bisection search and parabolic
 *            interpolation; should exhibit superlinear convergence in
 *            most functions.
 *
 *
 * Args:      ori     - n-vector at origin
 *            dir     - direction vector (gradient) we're following from ori
 *            n       - dimensionality of ori, dir, and xvec
 *            (*func) - ptr to caller's objective function
 *            prm     - ptr to any additional data (*func)() needs
 *            a,b     - minimum is bracketed on interval [a,b]
 *            eps     - tol = eps |x| + t; eps >= 2 * relative machine precision
 *            t       - additional factor for tol to avoid x=0 case.
 *            xvec    - RETURN: minimum, as an n-vector (caller allocated)
 *            ret_x   - optRETURN: scalar multiplier that gave xvec
 *            ret_fx  - optRETURN: f(x)
 *
 * Returns:   (void)
 *
 * Reference: See [Brent73], Chapter 5. My version is derived directly
 *            from Brent's description and his ALGOL60 code. I've
 *            preserved his variable names as much as possible, to
 *            make the routine follow his published description
 *            closely. The Brent algorithm is also discussed in
 *            Numerical Recipes [Press88].
 */
static void
brent(double *ori, double *dir, int n,
      double (*func)(double *, int, void *), void *prm,
      double a, double b, double eps, double t,
      double *xvec, double *ret_x, double *ret_fx)
{
  double w,x,v,u;               /* with [a,b]: Brent's six points     */
  double m;                     /* midpoint of current [a,b] interval */
  double tol;                   /* tolerance = eps|x| + t */
  double fu,fv,fw,fx;           /* function evaluations */
  double p,q;                   /* numerator, denominator of parabolic interpolation */
  double r;
  double d,e;                   /* last, next-to-last values of p/q  */
  double c = 1. - (1./eslCONST_GOLD); /* Brent's c; 0.381966; golden ratio */

  x=v=w= a + c*(b-a);           /* initial guess of x by golden section */
  esl_vec_DCopy(xvec, ori, n);  /* build xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  fx = (*func)(xvec, n, prm);   /* initial function evaluation */

  e = 0.;
  while (1) /* algorithm is guaranteed to converge. */
    {
      m   = 0.5 * (a+b);
      tol = eps*fabs(x) + t;
      if (fabs(x-m) <= 2*tol - 0.5*(b-a)) break; /* convergence test. */

      p = q = r = 0.;
      if (fabs(e) > tol)
        { /* Compute parabolic interpolation, u = x + p/q */
          r = (x-w)*(fx-fv);
          q = (x-v)*(fx-fw);
          p = (x-v)*q - (x-w)*r;
          q = 2*(q-r);
          if (q > 0) { p = -p; } else {q = -q;}
          r = e;
          e=d;                  /* e is now the next-to-last p/q  */
        }

      if (fabs(p) < fabs(0.5*q*r) || p < q*(a-x) || p < q*(b-x))
        { /* Seems well-behaved? Use parabolic interpolation to compute new point u */
          d = p/q;              /* d remembers last p/q */
          u = x+d;              /* trial point, for now... */

          if (2.0*(u-a) < tol || 2.0*(b-u) < tol) /* don't evaluate func too close to a,b */
            d = (x < m)? tol : -tol;
        }
      else /* Badly behaved? Use golden section search to compute u. */
        {
          e = (x<m)? b-x : a-x;  /* e = largest interval */
          d = c*e;
        }

      /* Evaluate f(), but not too close to x.  */
      if      (fabs(d) >= tol) u = x+d;
      else if (d > 0)          u = x+tol;
      else                     u = x-tol;
      esl_vec_DCopy(xvec, ori, n);  /* build xvec from ori, dir, u */
      esl_vec_DAddScaled(xvec, dir, u, n);
      fu = (*func)(xvec, n, prm);   /* f(u) */

      /* Bookkeeping.  */
     if (fu <= fx)
        {
          if (u < x) b = x; else a = x;
          v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        }
      else
        {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x)
            { v = w; fv = fw; w = u; fw = fu; }
          else if (fu <= fv || v==x || v ==w)
            { v = u; fv = fu; }
        }
    }

  /* Return.
   */
  esl_vec_DCopy(xvec, ori, n);  /* build final xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  if (ret_x  != NULL) *ret_x  = x;
  if (ret_fx != NULL) *ret_fx = fx;
  ESL_DPRINTF1(("xx=%10.8f fx=%10.1f\n", x, fx));
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
 *            It stops when the objective function has not changed
 *            by more than a fraction <tol>. This should not be smaller
 *            than sqrt(DBL_EPSILON). 
 *            
 *            Upon convergence, it sets <x>, the vector at the minimum (caller
 *            provides this allocated memory); sets <ret_xx),
 *            the scalar multiplier used to find <x> from <ori>
 *            and the gradient; and sets <ret_fx>, the value of the
 *            objective function at the minimum <x>.
 *            
 *            <ori> and <d> are unchanged by this procedure.
 *            <b> is a temporary workspace.
 *            
 * Args:      ori	- n-vector at origin
 *            d         - gradient from ori
 *            u         - natural units for each dimension i
 *            n         - dimensionality of all vectors
 *            (*func)   - pointer to caller's objective function
 *            prm       - void pointer to any data that (*func) needs
 *            tol       - tolerance: test for convergence
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
esl_min_LineSearch(double *ori, double *d, double *u, int n,
		   double (*func)(double *, int, void *), void *prm,
		   double tol, double *b, 
		   double *x, double *ret_xx, double *ret_fx)
{
  double xx, swap;
  double fa, fb, fx;
  double ax, bx, cx;
  int    niter;

  /* Bracket the minimum along line d first; this gives us scalar
   * points $0 <= ax < bx < cx$ relative to vector d. That is, vectors
   * at these points are determined as, for example, $a = ori + ax*d$.
   */
  ax = 0.;
  esl_min_Bracket(ori, d, u, n, func, prm, &fa, b, &bx, &fb, x, &cx, &fx);

#if eslDEBUGLEVEL >= 2
   show_line(ori, d, u, n, ax, cx, func, prm, b); */
#endif

  /* the main do loop is guaranteed to terminate; 
   * iterations are counted only for debugging output
   */
  niter = 0;
  do {	
    niter++;
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

    /* printf("golden section: %g %g %g\n", ax,bx,cx); */

  } while ((cx-ax)/bx > tol); 

  ESL_DPRINTF1(("xx=%10.8f fx=%10.1f\n", bx, fx));
  ESL_DPRINTF2(("\nesl_min_LineSearch(): %d iterations\n", niter));
  ESL_DPRINTF2(("esl_min_LineSearch(): moved to %g along curr direction\n", bx));
  ESL_DPRINTF2(("esl_min_LineSearch(): new f() is %g\n\n", fx));
  
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
 *            The first step of each iteration is to try to bracket
 *            the minimum along the current direction. The initial step
 *            size is controlled by <u[]>; the first step will not exceed 
 *            <u[i]> for any dimension <i>. (You can think of <u> as
 *            being the natural "units" to use along a graph axis, if
 *            you were plotting the objective function.)
 *
 *            The caller also provides an allocated workspace sufficient to
 *            hold four allocated n-vectors. (4 * sizeof(double) * n).
 *
 *            Iterations continue until the objective function has changed
 *            by less than a fraction <tol>. This should not be set to less than
 *            sqrt(DBL_EPSILON). 
 *
 *            Upon return, <x> is the minimum, and <ret_fx> is f(x),
 *            the function value at <x>.
 *            
 * Args:      x        - an initial guess n-vector; RETURN: x at the minimum
 *            u        - "units": maximum initial step size along gradient when bracketing.
 *            n        - dimensionality of all vectors
 *            *func()  - function for computing objective function f(x)
 *            *dfunc() - function for computing a gradient at x
 *            prm      - void ptr to any data/params func,dfunc need 
 *            tol      - convergence criterion applied to f(x)
 *            wrk      - allocated 4xn-vector for workspace
 *            ret_fx   - optRETURN: f(x) at the minimum
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECONVERGENCE> if it fails to converge in MAXITERATIONS.
 *
 * Xref:      STL9/101.
 */
int
esl_min_ConjugateGradientDescent(double *x, double *u, int n, 
       				 double (*func)(double *, int, void *),
				 void (*dfunc)(double *, int, void *, double *),
				 void *prm, double tol, double *wrk, double *ret_fx)
{
  double oldfx;
  double coeff;
  int    i, i1;
  double *dx, *cg, *w1, *w2;
  double cvg;
  double fa,fb,fc;
  double ax,bx,cx;

  dx = wrk;
  cg = wrk + n;
  w1 = wrk + 2*n;
  w2 = wrk + 3*n;

  oldfx = (*func)(x, n, prm);	/* init the objective function */
  (*dfunc)(x, n, prm, dx);	/* find the current negative gradient, - df(x)/dxi  */

  esl_vec_DScale(dx, n, -1.0);
  esl_vec_DCopy(cg, dx, n);	/* and make that the first conjugate direction, cg  */

  for (i = 0; i < MAXITERATIONS; i++)
    {
      ESL_DPRINTF1(("p: %5.2f %7.4f %7.4f  cg: %9.2f %9.2f %9.2f ", 
		    x[0], x[1], x[2], cg[0], cg[1], cg[2]));

      ax = 0.;
      esl_min_Bracket(x, cg, u, n, func, prm, &fa, w1, &bx, &fb, w2, &cx, &fc);

      /* Minimize along the line given by the conjugate gradient <cg> */
      /*esl_min_LineSearch(x, cg, u, n, func, prm, 1e-4, w1, w2, NULL, ret_fx);*/
      brent(x, cg, n, func, prm, ax,cx, 1e-7, 1e-8, w2, NULL, ret_fx);
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

      /* Convergence test. 1e-9 factor is fudging the case where our
       * minimum is at exactly f()=0.
       */
      cvg = 2.0 * fabs((oldfx-*ret_fx)) / (1e-9 + fabs(oldfx) + fabs(*ret_fx));

#if eslDEBUGLEVEL >= 2
      printf("\nesl_min_ConjugateGradientDescent():\n");
      printf("new point:     ");
      for (i = 0; i < n; i++)
	printf("%g ", x[i]);

      printf("\nnew gradient:  ");
      for (i = 0; i < n; i++)
	printf("%g ", dx[i]);

      printf("\nnew direction: ");
      for (i = 0; i < n; i++)
	printf("%g ", cg[i]);

      printf("\nOld f() = %g    New f() = %g    Convergence = %g\n\n", oldfx, *ret_fx, cvg);
#endif

      if (cvg <= tol) break;
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
  double wrk[24];
  double fx;
  int    i;

  esl_min_ConjugateGradientDescent(x, n, 
				   &example_func, &example_dfunc, (void *) a, 
				   0.0001, wrk, &fx);

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
