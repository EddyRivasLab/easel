/* Statistical routines for normal distributions
 * (This module is currently incomplete and not fully tested,
 *  though it compiles cleanly.)
 * 
 * SRE, Tue Nov 21 14:12:59 2006 [Janelia]
 * SVN $Id$
 */

#include "esl_config.h"

#include <math.h>

#include "easel.h"
#include "esl_normal.h"

/*****************************************************************
 * 1. Densities and distributions.
 *****************************************************************/

/* Function:  esl_normal_pdf()
 * Incept:    SRE, Tue Nov 21 14:15:43 2006 [Janelia]
 *
 * Purpose:   Calculates the normal (Gaussian) probability density
 *            function $P(X=x)$ for a normal distribution, given
 *            value <x>, mean <mu>, and standard deviation <sigma>.
 * 
 * Xref:      STL11/94.
 */
double 
esl_normal_pdf(double x, double mu, double sigma)
{
  double z;
  
  z = (x - mu) / sigma;
  return  exp(-z*z) / (sigma * sqrt(2. * eslCONST_PI));
}

/* Function:  esl_normal_logpdf()
 * Incept:    SRE, Tue Jan  9 20:43:52 2007 [Casa de Gatos]
 *
 * Purpose:   Calculates the log of the probabiility density function
 *            for the normal (Gaussian), $\log P(X=x)$, given value
 *            <x>, mean <mu>, and standard deviation <sigma>.
 *
 * Xref:      STL11/94.
 */
double
esl_normal_logpdf(double x, double mu, double sigma)
{
  double z;

  z = (x - mu) / sigma;
  return  (-z*z) - log(sigma) - log(sqrt(2.*eslCONST_PI));
}

/* Function:  esl_normal_cdf()
 * Incept:    SRE, Tue Jan  9 20:59:04 2007 [Casa de Gatos]
 *
 * Purpose:   Calculates the cumulative distribution function for the
 *            normal, $P(X \leq x)$, given value <x>, mean <mu>,
 *            and standard deviation <sigma>.
 *
 * Xref:      STL11/94.
 */
double
esl_normal_cdf(double x, double mu, double sigma)
{
  double z;

  z = (x - mu) / sigma;
  return 0.5 + 0.5 * erf(z / sqrt(2.));
}

/* Function:  esl_normal_surv()
 * Incept:    SRE, Thu Jan 11 20:16:23 2007 [Casa de Gatos]
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ (that is,
 *            1-CDF, the right tail probability mass) for a normal
 *            distribution, given value <x>, mean <mu>, and standard
 *            deviation <sigma>.
 *
 * Xref:      STL11/94
 */
double
esl_normal_surv(double x, double mu, double sigma)
{
  double z;

  z = (x - mu) / sigma;
  return erfc(z / sqrt(2.));
}


/*****************************************************************
 * Unit tests.
 *****************************************************************/
#ifdef eslNORMAL_TESTDRIVE
static int
utest_pdf(void)
{
  double mu    = 0.;
  double sigma = 1.;
  double x;
  double newpdf, lastpdf;

  x = 0.;
  newpdf = esl_normal_pdf(x, mu, sigma);
  do {
    x += 1.;
    lastpdf = newpdf;
    newpdf = esl_normal_pdf(x, mu, sigma);
  } while (newpdf > 0.);

  if (lastpdf > 1e-300) esl_fatal("Dynamic range of esl_normal_pdf insufficient");
  return eslOK;
}

static int
utest_logpdf(void)
{
  double mu    = 0.;
  double sigma = 1.;
  double x;
  double old, new;

  x = 0.;
  new = esl_normal_logpdf(x, mu, sigma);
  do {
    x += 1.;
    old = new;
    new = esl_normal_logpdf(x, mu, sigma);
  } while (new > 0.);

  if (old > 1e-300) esl_fatal("Dynamic range of esl_normal_logpdf insufficient");

  old = esl_normal_pdf(42, -5., 2.1);
  new = exp(esl_normal_logpdf(42, -5., 2.1));
  if (esl_DCompare(old, new, 1e-9) != eslOK)
    esl_fatal("logpdf and pdf aren't giving the same answer");
  return eslOK;
}

static int
utest_cdf(void)
{
  double mu    = 0.;
  double sigma = 1.;
  double x;
  double new;

  for (x = 0.; x > -100.; x -= 1.0)
    {
      new = esl_normal_cdf(x, mu, sigma);
      printf("%.0f %g\n", x, new);
    }
  return eslOK;
}


static int
utest_surv(void)
{
  double mu    = 0.;
  double sigma = 1.;
  double x;
  double new;

  for (x = 0.; x < 100.; x += 1.0)
    {
      new = esl_normal_surv(x, mu, sigma);
      printf("%.0f %g\n", x, new);
    }
  return eslOK;
}
  


#endif /*eslNORMAL_TESTDRIVE*/




/*****************************************************************
 * Test driver.
 *****************************************************************/
#ifdef eslNORMAL_TESTDRIVE
/* Compile:
   gcc -g -Wall -I. -L. -o test -DeslNORMAL_TESTDRIVE esl_normal.c -leasel -lm
*/
#include <stdio.h>
#include <math.h>
#include "easel.h"
#include "esl_normal.h"

int
main(int argc, char **argv)
{
  utest_pdf();
  utest_logpdf();
  /* utest_cdf(); */
  utest_surv();

  return eslOK;
}
#endif /*eslNORMAL_TESTDRIVE*/

/*****************************************************************
 * Example.
 *****************************************************************/

#ifdef eslNORMAL_EXAMPLE
/*::cexcerpt::normal_example::begin::*/
/* compile:
   gcc -g -Wall -I. -o example -DeslNORMAL_EXAMPLE esl_normal.c easel.c -lm
 */
#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_normal.h"

int
main(int argc, char **argv)
{
  double z;

  z = sqrt(2 * eslCONST_PI);
  printf("%.60f\n", z);
  printf("%.60f\n", eslCONST_PI);
  printf("%.60f\n", (1. + sqrt(5.)) / 2.);
  return 0;
}



#endif /*eslNORMAL_EXAMPLE*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
