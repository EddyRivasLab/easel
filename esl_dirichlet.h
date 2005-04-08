/* esl_dirichlet.h
 * Functions relevant to Beta, gamma, and Dirichlet densities,
 * and simple and mixture Dirichlet priors.
 * 
 * SRE, Tue Nov  2 14:35:06 2004 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_DIRICHLET_INCLUDED
#define ESL_DIRICHLET_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <easel/random.h>
#endif

/* Structure: MIXDCHLET
 * 
 * A mixture Dirichlet density, usually used as a prior 
 * for a multinomial model (turning count vectors into probability
 * parameters).
 */
typedef struct {
  double  *pq;			/* mixture coefficients pq[0..N-1]          */
  double **alpha;               /* Dirichlet params alpha[0..N-1][0..K-1]   */
  int      N;			/* number of mixtures, e.g. 9 for Sjolander */
  int      K;			/* alphabet size, e.g. 20                   */
} ESL_MIXDCHLET;

extern ESL_MIXDCHLET *esl_mixdchlet_Create(int N, int K);
extern void           esl_mixdchlet_Destroy(ESL_MIXDCHLET *pri);

extern int esl_dirichlet_LogProbData(double *c, double *alpha, int K, double *ret_answer);
extern int esl_dirichlet_LogGamma(double x, double *ret_answer);


/* Optional sampling code, when augmented by random module.
 */
#ifdef eslAUGMENT_RANDOM
extern int esl_dirichlet_SampleGamma(ESL_RANDOMNESS *r, double a, double *ret_answer)
extern int esl_dirichlet_Sample(ESL_RANDOMNESS *r, double *alpha, int K, double *p);
#endif /*eslAUGMENT_RANDOM*/

#endif /*ESL_DIRICHLET_INCLUDED*/
