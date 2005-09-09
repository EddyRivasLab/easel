/* random.h
 * 
 * Header file for random.c: 
 * Easel's portable, threadsafe random number generator.
 * 
 * SRE, Wed Jul 14 11:23:57 2004 [St. Louis]
 * SVN $Id: esl_random.h 35 2005-02-24 23:18:34Z eddy $
 */
#ifndef ESL_RANDOM_INCLUDED
#define ESL_RANDOM_INCLUDED

struct esl_randomness_s {
  long  seed; 		 /* if >0; reseed with this value */
  long  rnd1;		 /* random number from LCG1       */
  long  rnd2;            /* random number from LCG2       */
  long  rnd;             /* random number we return       */
  long  tbl[64];	 /* table for Bays/Durham shuffle */
};
typedef struct esl_randomness_s ESL_RANDOMNESS;

/* esl_rnd_Choose(a) chooses a uniformly distributed integer
 * in the range 0..a-1, given an initialized ESL_RANDOMNESS r.
 */
#define esl_rnd_Choose(r, a)    ((int) (esl_random(r) * (a)))


extern ESL_RANDOMNESS *esl_randomness_Create(long seed);
extern ESL_RANDOMNESS *esl_randomness_CreateTimeseeded(void);
extern void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
extern int             esl_randomness_Init(ESL_RANDOMNESS *r, long seed);

extern double esl_random(ESL_RANDOMNESS *r);

extern double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
extern double esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev);
extern double esl_rnd_Gamma(ESL_RANDOMNESS *r, double a);
extern int    esl_rnd_DChoose(ESL_RANDOMNESS *r, double *p, int N);
extern int    esl_rnd_FChoose(ESL_RANDOMNESS *r, float *p, int N);
extern char  *esl_rnd_IID(ESL_RANDOMNESS *r, char *alphabet, double *p, int n, 
			  int len);

#endif /*ESL_RANDOM_INCLUDED*/
