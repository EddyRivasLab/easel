/* Statistical routines for normal distributions
 * 
 * SRE, Tue Nov 21 14:12:59 2006 [Janelia]
 * SVN $Id$
 */

#include <esl_config.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_normal.h>

/*****************************************************************
 * 1. Densities and distributions.
 *****************************************************************/

/* Function:  esl_normal_pdf()
 * Incept:    SRE, Tue Nov 21 14:15:43 2006 [Janelia]
 *
 * Purpose:   Calculates the probability density function for
 *            a normal distribution, given 
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */


double 
esl_normal_pdf(double x, double mu, double sigma)
{
  double z;

  z = (x - mu) / sigma;
  

  return z;
}



/*****************************************************************
 * Example.
 *****************************************************************/

#ifdef eslNORMAL_EXAMPLE
/*::cexcerpt::normal_example::begin::*/
/* compile:
   gcc -g -Wall -I. -o example -DeslNORMAL_EXAMPLE\
     -DeslAUGMENT_HISTOGRAM -DeslAUGMENT_RANDOM -DeslAUGMENT_STATS\
     esl_normal.c esl_histogram.c esl_random.c esl_stats.c esl_vectorops.c easel.c -lm
 */
#include <stdio.h>
#include <math.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_normal.h>

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
