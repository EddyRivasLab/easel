/* esl_minimizer.h
 * Multidimensional optimization by conjugate gradient descent.
 * 
 * SRE, Wed Jun 22 09:53:05 2005
 * SVN $Id$
 */


#ifndef ESL_MINIMIZER_INCLUDED

#define MAXITERATIONS 100

extern int esl_min_Bracket(double *a, double *d, int n, 
			   double (*func)(double *, int, void *), void *prm, 
			   double xinit, double *ret_fa,
			   double *b, double *ret_bx, double *ret_fb,
			   double *c, double *ret_cx, double *ret_fc);
extern int esl_min_LineSearch(double *ori, double *d, int n,
			      double (*func)(double *, int, void *), void *prm,
			      double *b, 
			      double *x, double *ret_xx, double *ret_fx);
extern int esl_min_ConjugateGradientDescent(double *x, int n, 
					    double (*func)(double *, int, void *),
					    void (*dfunc)(double *, int, void *, double *),
					    void *prm,
					    double *dx, double *cg, double *w1, double *w2,
					    double *ret_fx);

#endif /*ESL_MINIMIZER_INCLUDED*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
