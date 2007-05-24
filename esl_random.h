/* random.h
 * 
 * Header file for random.c: 
 * Easel's portable, threadsafe random number generator.
 * 
 * SRE, Wed Jul 14 11:23:57 2004 [St. Louis]
 * SVN $Id$
 */
#ifndef ESL_RANDOM_INCLUDED
#define ESL_RANDOM_INCLUDED

typedef struct {
  long  seed;           /* reseed with this value, >0    */
  long  rnd1;           /* random number from LCG1       */
  long  rnd2;           /* random number from LCG2       */
  long  rnd;            /* random number we return       */
  long  tbl[64];        /* table for Bays/Durham shuffle */
  int   reseeding;	/* TRUE if seed is new           */
} ESL_RANDOMNESS;

/* esl_rnd_Choose(a) chooses a uniformly distributed integer
 * in the range 0..a-1, given an initialized ESL_RANDOMNESS r.
 */
#define esl_rnd_Choose(r, a)    ((int) (esl_random(r) * (a)))


/* 1. The ESL_RANDOMNESS object.
 */
extern ESL_RANDOMNESS *esl_randomness_Create(long seed);
extern ESL_RANDOMNESS *esl_randomness_CreateTimeseeded(void);
extern void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
extern int             esl_randomness_Init(ESL_RANDOMNESS *r, long seed);
extern long            esl_randomness_GetSeed(const ESL_RANDOMNESS *r);

/* 2. The generator, esl_random().
 */
extern double esl_random(ESL_RANDOMNESS *r);

/* 3. Other fundamental sampling (including Gaussian, gamma).
 */
extern double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
extern double esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev);
extern double esl_rnd_Gamma(ESL_RANDOMNESS *r, double a);

/* 4. Multinomial sampling from discrete probability n-vectors.
 */
extern int    esl_rnd_DChoose(ESL_RANDOMNESS *r, const double *p, int N);
extern int    esl_rnd_FChoose(ESL_RANDOMNESS *r, const float *p, int N);

/* 5. Generating iid sequences, either text or digital mode.
 */
extern int esl_rnd_IID  (ESL_RANDOMNESS *r, const char *alphabet, const double *p, int K, int L, char *s);
extern int esl_rnd_fIID (ESL_RANDOMNESS *r, const char *alphabet, const float  *p, int K, int L, char *s);
extern int esl_rnd_xIID (ESL_RANDOMNESS *r, const double *p, int K, int L, ESL_DSQ *dsq);
extern int esl_rnd_xfIID(ESL_RANDOMNESS *r, const float  *p, int K, int L, ESL_DSQ *dsq);

/* 6. Randomizing sequences.
 */
extern int esl_rnd_CShuffle  (ESL_RANDOMNESS *r, const char *s, char *shuffled);
extern int esl_rnd_CShuffleDP(ESL_RANDOMNESS *r, const char *s, char *shuffled);
extern int esl_rnd_CMarkov0  (ESL_RANDOMNESS *r, const char *s, char *markoved);
extern int esl_rnd_CMarkov1  (ESL_RANDOMNESS *r, const char *s, char *markoved);
extern int esl_rnd_CReverse  (const char *s, char *rev);
extern int esl_rnd_CShuffleWindows(ESL_RANDOMNESS *r, const char *s, int w, char *shuffled);

extern int esl_rnd_XShuffle  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, ESL_DSQ *shuffled);
extern int esl_rnd_XShuffleDP(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *shuffled);
extern int esl_rnd_XMarkov0  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved);
extern int esl_rnd_XMarkov1  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved);
extern int esl_rnd_XReverse(const ESL_DSQ *dsq, int L, ESL_DSQ *rev);
extern int esl_rnd_XShuffleWindows(ESL_RANDOMNESS *r, const ESL_DSQ *s, int L, int w, ESL_DSQ *shuffled);

/* 7. Randomizing alignments.
 */

#endif /*ESL_RANDOM_INCLUDED*/

