/* random.c
 * 
 * Easel's portable, threadsafe random number generator.
 * 
 * SRE, Wed Jul 14 10:54:46 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <easel/easel.h>
#include <easel/random.h>

/* Function:  esl_random_seed()
 * Incept:    SRE, Wed Jul 14 13:02:18 2004 [St. Louis]
 *
 * Purpose:   Initialize the random number generator using
 *            a given random seed. Seed must be >0.
 *            Returns an ESL_RANDOMNESS object that 
 *            the esl_random() generator will use.
 *            
 * Args:      seed >= 0.
 *
 * Returns:   initialized ESL_RANDOMNESS * on success;
 *            caller free's with esl_random_free().
 *              
 *            NULL on failure.
 * 
 * Xref:      STL8/p57.
 */
ESL_RANDOMNESS *
esl_random_seed(long seed)
{
  ESL_RANDOMNESS *r;
  int             burnin = 7;

  if (seed <= 0)                                    ESL_ERROR_VAL(NULL, ESL_EINVAL, "bad seed");
  if ((r = malloc(sizeof(ESL_RANDOMNESS))) == NULL) ESL_ERROR_VAL(NULL, ESL_EMEM,   "malloc failed");
  r->seed = seed;
  /* we observe that the first random number
   * isn't very random, if closely spaced seeds are used, like what
   * we get with using time().
   * So, "burn in" the random chain just a little.
   */
  while (burnin--) esl_random(r);
  return r;
}

/* Function:  esl_random_seed_time()
 * Incept:    SRE, Wed Jul 14 11:22:54 2004 [St. Louis]
 *
 * Purpose:   Like esl_random_seed(), but initializes the
 *            the random number generator using a time() call 
 *            (# of msec since the POSIX epoch).
 *
 * Args:      (void)
 *
 * Returns:   initialized ESL_RANDOMNESS * on success;
 *            caller free's with esl_random_free().
 *              
 *            NULL on failure.
 * 
 * Xref:      STL8/p57.
 */
ESL_RANDOMNESS *
esl_random_seed_time(void)
{
  ESL_RANDOMNESS *r;
  int             burnin = 7;

  if ((r = malloc(sizeof(ESL_RANDOMNESS))) == NULL) ESL_ERROR_VAL(NULL, ESL_EMEM,   "malloc failed");
  r->seed = time ((time_t *) NULL);
  while (burnin--) esl_random(r);
  return r;
}

/* Function:  esl_random_free()
 * Incept:    SRE, Wed Jul 14 13:19:08 2004 [St. Louis]
 *
 * Purpose:   Frees an ESL_RANDOMNESS object.
 */
int
esl_random_free(ESL_RANDOMNESS *r)
{
  free(r);
  return ESL_OK;
}


/* Function:  esl_random_reseed()
 * Incept:    SRE, Wed Jul 14 13:13:05 2004 [St. Louis]
 *
 * Purpose:   Reset and reinitialize an existing ESL_RANDOMNESS
 *            object. (Not generally recommended.)
 *
 * Args:      r     - randomness object
 *            seed  - new seed to use; >0.
 *
 * Returns:   ESL_OK on success.
 * Xref:      STL8/p57.
 */
int
esl_random_reseed(ESL_RANDOMNESS *r, long seed)
{
  int burnin = 7;
  if (seed <= 0) ESL_ERROR(ESL_EINVAL, "bad seed");
  r->seed = seed;
  while (burnin--) esl_random(r);
  return ESL_OK;
}



/* Function: esl_random()
 * 
 * Purpose:  Returns a uniform deviate x, 0.0 <= x < 1.0.
 * 
 *           The "randomness object" r contains the information
 *           we need for the generator; keeping it in an object
 *           (as opposed to static variables) makes us threadsafe.
 *           This object is created by esl_random_seed() or
 *           esl_random_seed_time(); it must be free'd eventually
 *           using esl_random_free().
 *
 *           If the seed in r is >0, that's a flag to reset
 *           and reinitialize the generator.
 *           
 * Method:   Implements L'Ecuyer's algorithm for combining output
 *           of two linear congruential generators, plus a Bays-Durham
 *           shuffle. This is essentially ran2() from Numerical Recipes,
 *           sans their nonhelpful Rand/McNally-esque code obfuscation.
 *           
 *           Overflow errors are avoided by Schrage's algorithm:
 *               az % m = a(z%q) - r(z/q) (+m if <0)
 *           where q=m/a, r=m%a
 *
 *           Requires that long int's have at least 32 bits.
 *           
 *           Reliable and portable, but slow. Benchmarks on wrasse,
 *           using Linux gcc and Linux glibc rand() (see randspeed, in Testsuite):
 *           sre_random():    0.5 usec/call
 *           rand():          0.2 usec/call
 *           
 * Reference: Press et al. Numerical Recipes in C, 1992. 
 */
double
esl_random(ESL_RANDOMNESS *r)
{
  long x,y;
  int i;
  /* Magic numbers a1,m1, a2,m2 from L'Ecuyer, for 2 LCGs.
   * q,r derive from them (q=m/a, r=m%a) and are needed for Schrage's algorithm.
   */
  long a1 = 40014;		
  long m1 = 2147483563;		
  long q1 = 53668;
  long r1 = 12211;

  long a2 = 40692;
  long m2 = 2147483399;
  long q2 = 52774;
  long r2 = 3791;

  if (r->seed > 0) 
    {
      r->rnd1 = r->seed;
      r->rnd2 = r->seed;
				/* Fill the table for Bays/Durham */
      for (i = 0; i < 64; i++) {
	x    = a1*(r->rnd1%q1);   /* LCG1 in action... */
	y    = r1*(r->rnd1/q1);
	r->rnd1 = x-y;
	if (r->rnd1 < 0) r->rnd1 += m1;

	x    = a2*(r->rnd2%q2);   /* LCG2 in action... */
	y    = r2*(r->rnd2/q2);
	r->rnd2 = x-y;
	if (r->rnd2 < 0) r->rnd2 += m2;

	r->tbl[i] = r->rnd1 - r->rnd2;
	if (r->tbl[i] < 0) r->tbl[i] += m1;
      }
      r->seed = 0;		/* drop the flag. */
    }/* end of initialization*/


  x    = a1*(r->rnd1%q1);   /* LCG1 in action... */
  y    = r1*(r->rnd1/q1);
  r->rnd1 = x-y;
  if (r->rnd1 < 0) r->rnd1 += m1;

  x    = a2*(r->rnd2%q2);   /* LCG2 in action... */
  y    = r2*(r->rnd2/q2);
  r->rnd2 = x-y;
  if (r->rnd2 < 0) r->rnd2 += m2;

   			/* Choose our random number from the table... */
  i   = (int) (((double) r->rnd / (double) m1) * 64.);
  r->rnd = r->tbl[i];
			/* and replace with a new number by L'Ecuyer. */
  r->tbl[i] = r->rnd1 - r->rnd2;
  if (r->tbl[i] < 0) r->tbl[i] += m1;

  return ((double) r->rnd / (double) m1);  
}


/* Function: esl_random_positive()
 * Incept:   SRE, Wed Jul 14 13:31:23 2004 [St. Louis]
 *
 * Purpose:  Same as esl_random(), but assure 0 < x < 1;
 *           (positive uniform deviate).
 */
double
esl_random_positive(ESL_RANDOMNESS *r)
{
  double x;
  do { x = esl_random(r); } while (x == 0.0);
  return x;
}


/* Function: esl_random_exponential()
 * Date:     SRE, Mon Sep  6 21:24:29 1999 [St. Louis]
 *
 * Purpose:  Pick an exponentially distributed random variable
 *           0 > x >= infinity
 *           
 * Args:     r   - an ESL_RANDOMNESS object
 *
 * Returns:  x
 */
double
esl_random_exponential(ESL_RANDOMNESS *r)
{
  return -log(esl_random_positive(r));
}    


/* Function:  esl_random_gaussian()
 * Incept:    SRE, Wed Jul 14 13:50:36 2004 [St. Louis]
 *
 * Purpose:   Pick a Gaussian-distributed random variable
 *            with some mean and standard deviation, and
 *            return it.
 * 
 * Method:    Based on RANLIB.c gennor() public domain implementation.
 *            Thanks to the authors, Barry W. Brown and James Lovato,
 *            University of Texas, M.D. Anderson Cancer Center, Houston TX.
 *            Their implementation is from Ahrens and Dieter, "Extensions 
 *            of Forsythe's method for random sampling from the normal
 *            distribution", Math. Comput. 27:927-937 (1973).
 *
 *            Impenetrability of the code is to be blamed on 
 *            FORTRAN/f2c lineage.
 *
 * Args:      r      - ESL_RANDOMNESS object
 *            mean   - mean of the Gaussian we're sampling from
 *            stddev - standard deviation of the Gaussian     
 *
 * Returns:   x
 */
double
esl_random_gaussian(ESL_RANDOMNESS *r, double mean, double stddev)
{
  long   i;
  double snorm,u,s,ustar,aa,w,y,tt;

  /* These static's are threadsafe: they are magic constants
   * we will not touch.
   */
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,    
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };

  u = esl_random(r);
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
   * START CENTER
   */
  ustar = u-(double)i;
  aa = a[i-1];
S40:
  if (ustar <= t[i-1]) goto S60;
  w = (ustar - t[i-1]) * h[i-1];
S50:
  /*
   * EXIT   (BOTH CASES)
   */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return (stddev*snorm + mean);
S60:
  /*
   * CENTER CONTINUED
   */
  u = esl_random(r);
  w = u*(a[i]-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = esl_random(r);
S80:
  if(ustar > tt) goto S50;
  u = esl_random(r);
  if(ustar >= u) goto S70;
  ustar = esl_random(r);
  goto S40;
S100:
  /*
   * START TAIL
   */
  i = 6;
  aa = a[31];
  goto S120;
S110:
  aa += d[i-1];
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u*d[i-1];
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = esl_random(r);
  if(ustar > tt) goto S50;
  u = esl_random(r);
  if(ustar >= u) goto S150;
  u = esl_random(r);
  goto S140;
}


/* Functions: esl_random_DChoose(), esl_random_FChoose()
 *
 * Purpose:   Make a random choice from a normalized distribution.
 *            DChoose() is for double-precision vectors;
 *            FChoose() is for single-precision float vectors.
 *            Returns the number of the choice.
 *
 * Limitation:
 *            All p's must be >> FLT_EPSILON or DBL_EPSILON.
 */
int
esl_random_DChoose(ESL_RANDOMNESS *r, double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }
  /* See comment on this next line in FChoose() */
  do { i = (int) (esl_random(r) * N); } while (p[i] == 0.);
  return i;
}
int
esl_random_FChoose(ESL_RANDOMNESS *r, float *p, int N)
{
  float roll;                   /* random fraction */
  float sum;                    /* integrated prob */
  int   i;                      /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }

  /* Very rarely, because of machine floating point representation,
   * our roll is "impossibly" >= total sum, even though any roll of
   * esl_random() is < 1.0 and the total sum is supposed to be 1.0 by
   * definition. This can happen when the total_sum is not really 1.0,
   * but something just less than that in the machine representation,
   * and the roll happens to also be very very close to 1. I have not
   * examined this analytically. Empirically, it occurs at a frequency
   * of about 1/10^8, as measured for bug #sq5. To work around, choose
   * one of the *nonzero* p[i]'s at random.  (If you chooose *any*
   * p[i] you get bug #sq5; routines like StrMarkov0() fail because
   * they choose impossible residues.)
   */
  do { i = (int) (esl_random(r) * N); } while (p[i] == 0.);
  return i;
}


/* Function: esl_random_sequence()
 * Incept:   SRE, Thu Aug  5 09:03:03 2004 [St. Louis]
 *
 * Purpose:  Generate an iid symbol sequence according
 *           to some alphabet, alphabet_size, probability
 *           distribution, and length. Return the
 *           sequence.
 *
 * Args:     r         - ESL_RANDOMNESS object
 *           alphabet  - e.g. "ACGT"
 *           p         - probability distribution [0..n-1]
 *           n         - number of symbols in alphabet
 *           len       - length of generated sequence
 *
 * Return:   ptr to random sequence, or NULL on failure.
 *           sequence is allocated here, must be free'd by caller.
 */
char *
esl_random_sequence(ESL_RANDOMNESS *r, char *alphabet, double *p, int n, int len)
{
  char *s;
  int   x;

  if ((s = (char *) malloc(sizeof(char) * (len+1))) == NULL) 
    ESL_ERROR_VAL(NULL, ESL_EMEM, "malloc failed");

  for (x = 0; x < len; x++)
    s[x] = alphabet[esl_random_DChoose(r,p,n)];
  s[x] = '\0';
  return s;
}

