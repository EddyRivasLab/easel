/* A portable, threadsafe random number generator.
 *
 *  1. The ESL_RANDOMNESS object.
 *  2. The generator, esl_random().
 *  3. Other fundamental sampling (including Gaussian, gamma).
 *  4. Multinomial sampling from discrete probability n-vectors.
 *  5. Unit tests.
 *  6. Test driver.
 *  7. An example of using the random module.
 *  
 * See http://csrc.nist.gov/rng/ for the NIST random number
 * generation test suite.
 * 
 * SRE, Wed Jul 14 10:54:46 2004 [St. Louis]
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_random.h"


/*****************************************************************
 *# 1. The <ESL_RANDOMNESS> object.
 *****************************************************************/

/* Function:  esl_randomness_Create()
 * Synopsis:  Create an RNG with a given seed.
 * Incept:    SRE, Wed Jul 14 13:02:18 2004 [St. Louis]
 *
 * Purpose:   Create a random number generator using
 *            a given random seed. Seed must be $>0$.
 *            
 * Args:      seed $>= 0$.
 *
 * Returns:   an initialized <ESL_RANDOMNESS *> on success.
 *            Caller free's with <esl_randomness_Destroy()>.
 *              
 * Throws:    <NULL> on failure.
 * 
 * Xref:      STL8/p57.
 */
ESL_RANDOMNESS *
esl_randomness_Create(long seed)
{
  ESL_RANDOMNESS *r      = NULL;
  int             burnin = 7;
  int             status;

  if (seed <= 0) ESL_XEXCEPTION(eslEINVAL, "bad seed");
  ESL_ALLOC(r, sizeof(ESL_RANDOMNESS));
  r->seed      = seed;
  r->reseeding = TRUE;

  /* we observe that the first random number isn't very random, if
   * closely spaced seeds are used, like what we get with using
   * time().  So, "burn in" the random chain just a little.
   */
  while (burnin--) esl_random(r);
  return r;

 ERROR:
  return NULL;
}

/* Function:  esl_randomness_CreateTimeseeded()
 * Synopsis:  Create an RNG with a quasirandom seed.
 * Incept:    SRE, Wed Jul 14 11:22:54 2004 [St. Louis]
 *
 * Purpose:   Like <esl_randomness_Create()>, but it initializes the
 *            the random number generator using a POSIX <time()> call 
 *            (number of seconds since the POSIX epoch).
 *
 * Returns:   an initialized <ESL_RANDOMNESS *> on success.
 *            Caller free's with <esl_randomness_Destroy()>.
 *              
 * Throws:    <NULL> on failure.
 * 
 * Xref:      STL8/p57.
 */
ESL_RANDOMNESS *
esl_randomness_CreateTimeseeded(void)
{
  ESL_RANDOMNESS *r      = NULL;
  int             burnin = 7;
  int             status;

  ESL_ALLOC(r, sizeof(ESL_RANDOMNESS));
  r->seed      = time ((time_t *) NULL);
  r->reseeding = TRUE;
  while (burnin--) esl_random(r);
  return r;

 ERROR:
  return NULL;
}

/* Function:  esl_randomness_Destroy()
 * Synopsis:  Free an RNG.            
 * Incept:    SRE, Wed Jul 14 13:19:08 2004 [St. Louis]
 *
 * Purpose:   Frees an <ESL_RANDOMNESS> object.
 */
void
esl_randomness_Destroy(ESL_RANDOMNESS *r)
{
  free(r);
  return;
}


/* Function:  esl_randomness_Init()
 * Synopsis:  Reinitialize an RNG.           
 * Incept:    SRE, Wed Jul 14 13:13:05 2004 [St. Louis]
 *
 * Purpose:   Reset and reinitialize an existing <ESL_RANDOMNESS>
 *            object. 
 *            
 *            (Not generally recommended. This does not make a
 *            sequence of numbers more random, and may make it less
 *            so.)
 *
 * Args:      r     - randomness object
 *            seed  - new seed to use; >0.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if seed is $<= 0$.
 *
 * Xref:      STL8/p57.
 */
int
esl_randomness_Init(ESL_RANDOMNESS *r, long seed)
{
  int burnin = 7;
  if (seed <= 0) ESL_EXCEPTION(eslEINVAL, "bad seed");
  r->seed      = seed;
  r->reseeding = TRUE;
  while (burnin--) esl_random(r);
  return eslOK;
}

/* Function:  esl_randomness_GetSeed()
 * Synopsis:  Returns the value of RNG's seed.
 * Incept:    SRE, Wed May 23 17:02:59 2007 [Janelia]
 *
 * Purpose:   Return the value of the seed. 
 * 
 *            (You already know what the seed was if you used
 *            <esl_randomness_Create()>, but not if you used
 *            <esl_randomness_CreateTimeseeded()>. It is often useful
 *            to record what the seed was, in order to be able to
 *            exactly reproduce results.)
 */
long
esl_randomness_GetSeed(const ESL_RANDOMNESS *r)
{
  return r->seed;
}
/*----------- end of ESL_RANDOMNESS object functions --------------*/



/*****************************************************************
 *# 2. The generator, <esl_random()>
 *****************************************************************/  

/* Function: esl_random()
 * Synopsis: Generate a uniform random deviate $0.0 <= x < 1.0$.
 *            
 * Purpose:  Returns a uniform deviate x, $0.0 <= x < 1.0$, given
 *           RNG <r>, using L'Ecuyer's algorithm for combining output
 *           of two linear congruential generators, plus a Bays-Durham
 *           shuffle \citep{Press93}.
 *           
 * Method:   This is essentially ran2() from Numerical Recipes,
 *           rewritten, sans their nonhelpful Rand/McNally-esque code
 *           obfuscation.
 *           
 *           Overflow errors are avoided by Schrage's algorithm:
 *               az % m = a(z%q) - r(z/q) (+m if <0)
 *           where q=m/a, r=m%a
 *
 *           Requires that long int's have at least 32 bits.
 *           
 *           Reliable and portable, but slow. Benchmarks on wrasse,
 *           using Linux gcc and Linux glibc rand() (see randspeed, Testsuite):
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

  if (r->reseeding) 
    {
      r->rnd1 = r->seed;
      r->rnd2 = r->seed;

      /* Fill the table for Bays/Durham; first 64 (0..63)
       * random #'s are for our table, 65th is to init r->rnd.
       */
      for (i = 0; i <= 64; i++) {
	x    = a1*(r->rnd1%q1);   /* LCG1 in action... */
	y    = r1*(r->rnd1/q1);
	r->rnd1 = x-y;
	if (r->rnd1 < 0) r->rnd1 += m1;

	x    = a2*(r->rnd2%q2);   /* LCG2 in action... */
	y    = r2*(r->rnd2/q2);
	r->rnd2 = x-y;
	if (r->rnd2 < 0) r->rnd2 += m2;

	if (i < 64) {
	  r->tbl[i] = r->rnd1 - r->rnd2;
	  if (r->tbl[i] < 0) r->tbl[i] += m1;
	} else {
	  r->rnd = r->rnd1 - r->rnd2;
	  if (r->rnd < 0) r->rnd += m1;
	}
      }
      r->reseeding = FALSE;	/* drop the flag. */
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
/*----------- end of esl_random() --------------*/



/*****************************************************************
 *# 3. Other fundamental sampling (including Gaussian, gamma)
 *****************************************************************/ 

/* Function: esl_rnd_UniformPositive()
 * Synopsis: Generate a uniform positive random deviate $0 < x < 1$.
 * Incept:   SRE, Wed Jul 14 13:31:23 2004 [St. Louis]
 *
 * Purpose:  Same as <esl_random()>, but assure $0 < x < 1$;
 *           (positive uniform deviate).
 */
double
esl_rnd_UniformPositive(ESL_RANDOMNESS *r)
{
  double x;
  do { x = esl_random(r); } while (x == 0.0);
  return x;
}


/* Function:  esl_rnd_Gaussian()
 * Synopsis:  Generate a Gaussian-distributed sample.
 * Incept:    SRE, Wed Jul 14 13:50:36 2004 [St. Louis]
 *
 * Purpose:   Pick a Gaussian-distributed random variable
 *            with a given <mean> and standard deviation <stddev>, and
 *            return it. 
 *            
 *            Implementation is derived from the public domain
 *            RANLIB.c <gennor()> function, written by Barry W. Brown
 *            and James Lovato (M.D. Anderson Cancer Center, Texas
 *            USA) using the method described in
 *            \citep{AhrensDieter73}.
 * 
 * Method:    Impenetrability of the code is to be blamed on 
 *            FORTRAN/f2c lineage.
 *
 * Args:      r      - ESL_RANDOMNESS object
 *            mean   - mean of the Gaussian we're sampling from
 *            stddev - standard deviation of the Gaussian     
 */
double
esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev)
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



/* subfunctions that esl_rnd_Gamma() is going to call:
 */
static double
gamma_ahrens(ESL_RANDOMNESS *r, double a)	/* for a >= 3 */
{
  double V;			/* uniform deviates */
  double X,Y;
  double test;
  
  do {
    do {				/* generate candidate X */
      Y = tan(eslCONST_PI * esl_random(r)); 
      X = Y * sqrt(2.*a -1.) + a - 1.;
    } while (X <= 0.);
				/* accept/reject X */
    V    = esl_random(r);
    test = (1+Y*Y) * exp( (a-1.)* log(X/(a-1.)) - Y*sqrt(2.*a-1.));
  } while (V > test);
  return X;
}
static double
gamma_integer(ESL_RANDOMNESS *r, unsigned int a)	/* for small integer a, a < 12 */
{
  int    i;
  double U,X;

  U = 1.;
  for (i = 0; i < a; i++) 
    U *= esl_rnd_UniformPositive(r);
  X = -log(U);

  return X;
}
static double
gamma_fraction(ESL_RANDOMNESS *r, double a)	/* for fractional a, 0 < a < 1 */
{				/* Knuth 3.4.1, exercise 16, pp. 586-587 */
  double p, U, V, X, q;
  
  p = eslCONST_E / (a + eslCONST_E);
  do {
    U = esl_random(r);
    V = esl_rnd_UniformPositive(r);
    if (U < p) {
      X = pow(V, 1./a);
      q = exp(-X);
    } else {
      X = 1. - log(V);
      q = pow(X, a-1.);
    }
    U = esl_random(r);
  } while (U >= q);
  return X;
}


/* Function: esl_rnd_Gamma()
 * Synopsis: Returns a random deviate from a Gamma(a, 1) distribution.
 * Incept:   SRE, Wed Apr 17 13:10:03 2002 [St. Louis]
 *
 * Purpose:  Return a random deviate distributed as Gamma(a, 1.)
 *           \citep[pp. 133--134]{Knu-81a}.
 *           
 *           The implementation follows not only Knuth \citep{Knu-81a},
 *           but also relied on examination of the implementation in
 *           the GNU Scientific Library (libgsl) \citep{Galassi06}.
 *
 * Args:     r      - random number generation seed
 *           a      - order of the gamma function; a > 0
 *
 * Throws:   <eslEINVAL> for $a <= 0$.
 */
double
esl_rnd_Gamma(ESL_RANDOMNESS *r, double a)
{
  double aint;

  aint = floor(a);
  if (a == aint && a < 12.) 
    return gamma_integer(r, (unsigned int) a);
  else if (a > 3.) 
    return gamma_ahrens(r, a);
  else if (a < 1.) 
    return gamma_fraction(r, a);
  else 
    return gamma_integer(r, aint) + gamma_fraction(r, a-aint);
  return eslOK;
}


/*****************************************************************
 *# 4. Multinomial sampling from discrete probability n-vectors
 *****************************************************************/ 

/* Function:  esl_rnd_DChoose()
 * Synopsis:  Return random choice from discrete multinomial distribution.          
 *
 * Purpose:   Make a random choice from a normalized discrete
 *            distribution <p> of <N> elements, where <p>
 *            is double-precision. Returns the index of the
 *            selected element, $0..N-1$.
 *            
 *            <p> must be a normalized probability distribution
 *            (i.e. must sum to one). Sampling distribution is
 *            undefined otherwise: that is, a choice will always
 *            be returned, but it might be an arbitrary one.
 *
 *            All $p_i$ must be $>>$ <DBL_EPSILON> in order to 
 *            have a non-zero probability of being sampled.
 *
 *            <esl_rnd_FChoose()> is the same, but for floats in <p>.
 *
 * Note:      Why the while (1) loop? Very rarely, because of machine
 *            floating point representation, our roll is "impossibly" 
 *            >= total sum, even though any roll of esl_random() is 
 *            < 1.0 and the total sum is supposed to be 1.0 by
 *            definition. This can happen when the total_sum is not
 *            really 1.0, but something just less than that in the 
 *            machine representation, and the roll happens to also be 
 *            very very close to 1. I have not examined this analytically, 
 *            but empirically, it occurs at a frequency of about 1/10^8
 *            as measured for bug #sq5... which suggests it is on the
 *            order of machine epsilon (not surprisingly). The while 
 *            loop makes you go around and try again; it must eventually
 *            succeed.
 *            
 *            The while() loop then makes the function vulnerable to
 *            an infinite loop if <p> sums to <=0 -- which shouldn't
 *            happen, but we shouldn't infinite loop if it does,
 *            either.  That's why there's a check on the sum of
 *            <p>. We return -1 in this case, a non-standard error code
 *            for Easel.
 * 
 * Throws:    -1 on failure. (This is a non-standard error code for Easel,
 *            but the only way an error can happen is if <p> isn't a 
 *            normalized probability distribution.)
 */
int
esl_rnd_DChoose(ESL_RANDOMNESS *r, const double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;

  while (1) {	/* see note in header about this while() */
    for (i = 0; i < N; i++)
      {
	sum += p[i];
	if (roll < sum) return i;  /* success! */
      }
    if (sum < 0.99) ESL_EXCEPTION(-1, "unnormalized distribution");    /* avoid inf loop */
  }
  /*UNREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}
int
esl_rnd_FChoose(ESL_RANDOMNESS *r, const float *p, int N)
{
  float  roll;                  /* random fraction */
  float  sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;

  while (1) {	/* see note in header about this while() */
    for (i = 0; i < N; i++)
      {
	sum += p[i];
	if (roll < sum) return i; /* success */
      }
    if (sum < 0.99) ESL_EXCEPTION(-1, "unnormalized distribution");    /* avoid inf loop */
  }
  /*UNREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}





/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/

#ifdef eslRANDOM_TESTDRIVE
#include "esl_vectorops.h"
#include "esl_stats.h"
#include "esl_dirichlet.h"
    
  
/* The esl_random() unit test:
 * a binned frequency test.
 */
static void
utest_random(long seed, int n, int nbins, int be_verbose)
{
  ESL_RANDOMNESS *r      = NULL;
  int            *counts = NULL;
  double          X2p    = 0.;
  int             i;
  double          X2, exp, diff;

  if ((counts = malloc(sizeof(int) * nbins)) == NULL) esl_fatal("malloc failed");
  esl_vec_ISet(counts, nbins, 0);

  /* This contrived call sequence exercises CreateTimeseeded() and
   * Init(), while leaving us a reproducible chain. Because it's
   * reproducible, we know this test succeeds, despite being
   * statistical in nature.
   */
  if ((r = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("randomness create failed");
  if (esl_randomness_Init(r, seed)            != eslOK) esl_fatal("randomness init failed");

  for (i = 0; i < n; i++)
    counts[esl_rnd_Roll(r, nbins)]++;

  /* X^2 value: \sum (o_i - e_i)^2 / e_i */
  for (X2 = 0., i = 0; i < nbins; i++) {
    exp  = (double) n / (double) nbins;
    diff = (double) counts[i] - exp;
    X2 +=  diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal("chi squared eval failed");
  if (be_verbose) printf("random():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal("chi squared test failed");

  esl_randomness_Destroy(r);
  free(counts);
  return;
}

/* The DChoose() and FChoose() unit tests.
 */
static void
utest_choose(ESL_RANDOMNESS *r, int n, int nbins, int be_verbose)
{
  double *pd = NULL;
  float  *pf = NULL;
  int    *ct = NULL;
  int     i;
  double  X2, diff, exp, X2p;

  if ((pd = malloc(sizeof(double) * nbins)) == NULL) esl_fatal("malloc failed"); 
  if ((pf = malloc(sizeof(float)  * nbins)) == NULL) esl_fatal("malloc failed");
  if ((ct = malloc(sizeof(int)    * nbins)) == NULL) esl_fatal("malloc failed");

  /* Sample a random multinomial probability vector.  */
  if (esl_dirichlet_DSampleUniform(r, nbins, pd) != eslOK) esl_fatal("dirichlet sample failed");
  esl_vec_D2F(pd, nbins, pf);

  /* Sample observed counts using DChoose(). */
  esl_vec_ISet(ct, nbins, 0);
  for (i = 0; i < n; i++)
    ct[esl_rnd_DChoose(r, pd, nbins)]++;

  /* X^2 test on those observed counts. */
  for (X2 = 0., i=0; i < nbins; i++) {
    exp = (double) n * pd[i];
    diff = (double) ct[i] - exp;
    X2 += diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal("chi square eval failed");
  if (be_verbose) printf("DChoose():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal("chi squared test failed");

  /* Repeat above for FChoose(). */
  esl_vec_ISet(ct, nbins, 0);
  for (i = 0; i < n; i++)
    ct[esl_rnd_FChoose(r, pf, nbins)]++;
  for (X2 = 0., i=0; i < nbins; i++) {
    exp = (double) n * pd[i];
    diff = (double) ct[i] - exp;
    X2 += diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal("chi square eval failed");
  if (be_verbose) printf("FChoose():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal("chi squared test failed");
  
  free(pd);
  free(pf);
  free(ct);
  return;
}
 

#endif /*eslRANDOM_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/


/*****************************************************************
 * 6. Test driver.
 *****************************************************************/
#ifdef eslRANDOM_TESTDRIVE
/* gcc -g -Wall -o random_utest -L. -I. -DeslRANDOM_TESTDRIVE esl_random.c -leasel -lm
 */
#include "esl_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  {"-b",  eslARG_INT,      "20", NULL, "n>0",NULL, NULL, NULL, "number of test bins",               0},
  {"-n",  eslARG_INT, "1000000", NULL, "n>0",NULL, NULL, NULL, "number of samples",                 0},
  {"-r",  eslARG_NONE,     NULL, NULL, NULL, NULL, NULL, NULL, "use arbitrary random number seed",  0},
  {"-s",  eslARG_INT,      "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",     0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose output",               0},
  {"--bitfile",eslARG_STRING,NULL,NULL,NULL, NULL, NULL, NULL, "save bit file for NIST benchmark",  0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for random module";

static int save_bitfile(char *bitfile, ESL_RANDOMNESS *r, int n);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r          = NULL;
  char           *bitfile    = esl_opt_GetString (go, "--bitfile");
  int             nbins      = esl_opt_GetInteger(go, "-b");
  int             n          = esl_opt_GetInteger(go, "-n");
  int             be_verbose = esl_opt_GetBoolean(go, "-v");
  int             seed       = esl_opt_GetInteger(go, "-s");

  if (esl_opt_GetBoolean(go, "-r")) r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(seed);

  utest_random(seed, n, nbins, be_verbose);
  utest_choose(r,    n, nbins, be_verbose);

  if (bitfile != NULL) save_bitfile(bitfile, r, n);

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

static int
save_bitfile(char *bitfile, ESL_RANDOMNESS *r, int n)
{
  FILE *fp = NULL;
  int b,i;
  long x;

  /* Open the file. 
   */
  if ((fp = fopen(bitfile, "w")) == NULL) 
    esl_fatal("failed to open %s for writing", bitfile);

  /* Sample <n> random numbers, output 31n random bits to the file.
   */
  for (i = 0; i < n; i++)
    {
      esl_random(r);
      x = r->rnd;		/* peek inside, get the 31 bit random long */

      for (b = 0; b < 31; b++)  /* don't print the sign bit. */
	{
	  if (x & 01) fprintf(fp, "1");
	  else        fprintf(fp, "0");
	  x >>= 1;
	}
      fprintf(fp, "\n");
    }
  fclose(fp);
  return eslOK;
}
#endif /*eslRANDOM_TESTDRIVE*/



/*****************************************************************
 * 7. An example of using the random module.
 *****************************************************************/
#ifdef eslRANDOM_EXAMPLE
/*::cexcerpt::random_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslRANDOM_EXAMPLE esl_random.c easel.c -lm
 * run:     ./example
 */
#include <stdio.h>
#include "easel.h"
#include "esl_random.h"

int 
main(void)
{
  ESL_RANDOMNESS *r = esl_randomness_Create(42); 
  int             n = 10;

  printf("A sequence of %d pseudorandom numbers:\n", n);
  while (n--)  printf("%f\n", esl_random(r));

  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::random_example::end::*/
#endif /*eslRANDOM_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/


