/* Finding roots of functions.
 * 
 * SRE, Fri Apr  6 10:01:43 2007 [Janelia]
 * SVN $Id$
 */
#ifndef ESL_ROOTFINDER_INCLUDED
#define ESL_ROOTFINDER_INCLUDED

typedef struct {
  double (*func)(double, void*);
  void   (*fdf) (double, void*, double*, double*);
  void   *params;

  double xl;
  double fl;
  double xr;
  double fr;

  double x0;
  double f0;

  double x;
  double fx;
  double dfx;
  int    iter;

  double abs_tolerance;
  double rel_tolerance;
  double residual_tol;
  int    max_iter;
} ESL_ROOTFINDER;


extern ESL_ROOTFINDER *esl_rootfinder_CreateBracketer(double (*func)(double, void*), void *params, double xl, double xr);
extern ESL_ROOTFINDER *esl_rootfinder_CreatePolisher(void (*fdf)(double, void*, double*, double*), void *params, double guess);
extern int esl_rootfinder_SetBrackets(ESL_ROOTFINDER *R, double xl, double xr);
extern int esl_rootfinder_SetAbsoluteTolerance(ESL_ROOTFINDER *R, double tol);
extern int esl_rootfinder_SetRelativeTolerance(ESL_ROOTFINDER *R, double tol);
extern int esl_rootfinder_SetResidualTolerance(ESL_ROOTFINDER *R, double tol);
extern int esl_rootfinder_SetMaxIterations(ESL_ROOTFINDER *R, int maxiter);
extern void esl_rootfinder_Destroy(ESL_ROOTFINDER *R);

extern int esl_root_Bisection(ESL_ROOTFINDER *R, double *ret_x);
extern int esl_root_Newton(ESL_ROOTFINDER *R, double *ret_x);


#endif /*ESL_ROOTFINDER_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
