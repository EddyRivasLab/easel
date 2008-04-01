/* A portable, threadsafe random number generator.
 *
 *  1. The ESL_RANDOMNESS object.
 *  2. The generator, esl_random().
 *  3. Other fundamental sampling (including Gaussian, gamma).
 *  4. Multinomial sampling from discrete probability n-vectors.
 *  5. Generating iid sequences, either text or digital mode.
 *  6. Randomizing sequences.
 *  7. Unit tests.
 *  8. The test driver.
 *  9. An example of using the random module.
 *  
 * See http://csrc.nist.gov/rng/ for the NIST random number
 * generation test suite.
 * 
 * SRE, Wed Jul 14 10:54:46 2004 [St. Louis]
 * SVN $Id$
 */
#include <esl_config.h>

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
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
  r->seed = seed;
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
 *# 5. Generating iid sequences, either text or digital mode.
 *****************************************************************/ 

/* Function: esl_rnd_IID()
 * Synopsis: Generate an iid random text sequence.
 * Incept:   SRE, Thu Aug  5 09:03:03 2004 [St. Louis]
 *
 * Purpose:  Generate a <NUL>-terminated i.i.d. symbol string of length <L>,
 *           $0..L-1$, and leave it in <s>. The symbol alphabet is given
 *           as a string <alphabet> of <K> total symbols, and the iid
 *           probability of each residue is given in <p>. The caller
 *           must provide an <s> that is allocated for at least
 *           <(L+1)*sizeof(char)>, room for <L> residues and the <NUL> terminator.
 *           
 *           <esl_rnd_fIID()> does the same, but for a floating point
 *           probability vector <p>, rather than a double precision
 *           vector.
 *
 * Args:     r         - ESL_RANDOMNESS object
 *           alphabet  - e.g. "ACGT"
 *           p         - probability distribution [0..n-1]
 *           K         - number of symbols in alphabet
 *           L         - length of generated sequence
 *           s         - the generated sequence.
 *                       Caller allocated, >= (L+1) * sizeof(char).
 *            
 * Return:   <eslOK> on success.
 */
int
esl_rnd_IID(ESL_RANDOMNESS *r, const char *alphabet, const double *p, int K, int L, char *s)
{
  int   x;

  for (x = 0; x < L; x++)
    s[x] = alphabet[esl_rnd_DChoose(r,p,K)];
  s[x] = '\0';
  return eslOK;
}
int
esl_rnd_fIID(ESL_RANDOMNESS *r, const char *alphabet, const float *p, int K, int L, char *s)
{
  int   x;

  for (x = 0; x < L; x++)
    s[x] = alphabet[esl_rnd_FChoose(r,p,K)];
  s[x] = '\0';
  return eslOK;
}


/* Function: esl_rnd_xIID()
 * Synopsis: Generate an iid random digital sequence.
 * Incept:   SRE, Sat Feb 17 16:39:01 2007 [Casa de Gatos]
 *
 * Purpose:  Generate an i.i.d. digital sequence of length <L> (1..L) and
 *           leave it in <dsq>. The i.i.d. probability of each residue is
 *           given in the probability vector <p>, and the number of
 *           possible residues (the alphabet size) is given by <K>.
 *           (Only the alphabet size <K> is needed here, as opposed to
 *           a digital <ESL_ALPHABET>, but the caller presumably
 *           has a digital alphabet.) The caller must provide a <dsq>
 *           allocated for at least <L+2> residues of type <ESL_DSQ>,
 *           room for <L> residues and leading/trailing digital sentinel bytes.
 *           
 *           <esl_rnd_xfIID()> does the same, but for a
 *           single-precision float vector <p> rather than a
 *           double-precision vector <p>.
 *
 * Args:     r         - ESL_RANDOMNESS object
 *           p         - probability distribution [0..n-1]
 *           K         - number of symbols in alphabet
 *           L         - length of generated sequence
 *           ret_s     - RETURN: the generated sequence. 
 *                       (Caller-allocated, >= (L+2)*ESL_DSQ)
 *
 * Return:   <eslOK> on success.
 */
int
esl_rnd_xIID(ESL_RANDOMNESS *r, const double *p, int K, int L, ESL_DSQ *dsq)
{
  int   x;

  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
  for (x = 1; x <= L; x++) 
    dsq[x] = esl_rnd_DChoose(r,p,K);
  return eslOK;
}
int
esl_rnd_xfIID(ESL_RANDOMNESS *r, const float *p, int K, int L, ESL_DSQ *dsq)
{
  int   x;

  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
  for (x = 1; x <= L; x++) 
    dsq[x] = esl_rnd_FChoose(r,p,K);
  return eslOK;
}

/*****************************************************************
 *# 6. Randomizing sequences.
 *****************************************************************/

/* Function:  esl_rnd_CShuffle()
 * Synopsis:  Shuffle a text sequence.
 * Incept:    SRE, Fri Feb 23 08:17:50 2007 [Casa de Gatos]
 *
 * Purpose:   Returns a shuffled version of <s> in <shuffled>, given
 *            a source of randomness <r>.
 *            
 *            Caller provides allocated storage for <shuffled>, for at
 *            least the same length as <s>.
 *
 *            <shuffled> may also point to the same storage as <s>,
 *            in which case <s> is shuffled in place.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_rnd_CShuffle(ESL_RANDOMNESS *r, const char  *s, char *shuffled)
{
  int  L, i;
  char c;

  L = strlen(s);
  if (shuffled != s) strcpy(shuffled, s);
  while (L > 1) {
    i             = esl_rnd_Choose(r, L);
    c             = shuffled[i];
    shuffled[i]   = shuffled[L-1];
    shuffled[L-1] = c;
    L--;
  }
  return eslOK;
}

/* Function:  esl_rnd_CShuffleDP()
 * Synopsis:  Shuffle a text sequence, preserving diresidue composition.
 * Incept:    SRE, Fri Feb 23 08:56:03 2007 [Casa de Gatos]
 *
 * Purpose:   Given string <s>, and a source of randomness <r>,
 *            returns shuffled version in <shuffled>. The shuffle
 *            is a "doublet-preserving" (DP) shuffle which
 *            shuffles a sequence while exactly preserving both mono-
 *            and di-symbol composition. 
 *            
 *            <s> may only consist of alphabetic characters [a-zA-Z].
 *            The shuffle is done case-insensitively. The shuffled
 *            string result is all upper case.
 *
 *            Caller provides storage in <shuffled> of at least the
 *            same length as <s>.
 *            
 *            <shuffled> may also point to the same storage as <s>,
 *            in which case <s> is shuffled in place.
 *            
 *            The algorithm does an internal allocation of a
 *            substantial amount of temporary storage, on the order of
 *            <26 * strlen(s)>, so an allocation failure is possible
 *            if <s> is long enough.
 *
 *            The algorithm is a search for a random Eulerian walk on
 *            a directed multigraph \citep{AltschulErickson85}.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains nonalphabetic characters.
 *            <eslEMEM> on allocation failure.
 */
int
esl_rnd_CShuffleDP(ESL_RANDOMNESS *r, const char *s, char *shuffled)
{
  int    status;          /* Easel return status code */
  int    len;	          /* length of s */
  int    pos;	          /* a position in s or shuffled */
  int    x,y;             /* indices of two characters */
  char **E  = NULL;       /* edge lists: E[0] is the edge list from vertex A */
  int   *nE = NULL;       /* lengths of edge lists */
  int   *iE = NULL;       /* positions in edge lists */
  int    n;	          /* tmp: remaining length of an edge list to be shuffled */
  char   sf;              /* last character in shuffled */
  char   Z[26];           /* connectivity in last edge graph Z */ 
  int    keep_connecting; /* flag used in Z connectivity algorithm */
  int    is_eulerian;	  /* flag used for when we've got a good Z */
  
  /* First, verify that the string is entirely alphabetic. */
  len = strlen(s);
  for (pos = 0; pos < len; pos++)
    if (! isalpha((int) s[pos]))
      ESL_EXCEPTION(eslEINVAL, "String contains nonalphabetic characters");

  /* Allocations. */
  ESL_ALLOC(E,  sizeof(char *) * 26);   for (x = 0; x < 26; x++) E[x] = NULL;
  ESL_ALLOC(nE, sizeof(int)    * 26);   for (x = 0; x < 26; x++) nE[x] = 0;
  ESL_ALLOC(iE, sizeof(int)    * 26);   for (x = 0; x < 26; x++) iE[x] = 0; 
  for (x = 0; x < 26; x++) 
    ESL_ALLOC(E[x], sizeof(char) * (len-1));

  /* "(1) Construct the doublet graph G and edge ordering E
   *      corresponding to S."
   * 
   * Note that these also imply the graph G; and note,
   * for any list x with nE[x] = 0, vertex x is not part
   * of G.
   */
  x = toupper((int) s[0]) - 'A';
  for (pos = 1; pos < len; pos++)
    {
      y = toupper((int) s[pos]) - 'A';
      E[x][nE[x]] = y;
      nE[x]++;
      x = y;
    }
  
  /* Now we have to find a random Eulerian edge ordering. */
  sf = toupper((int) s[len-1]) - 'A'; 
  is_eulerian = 0;
  while (! is_eulerian)
    {
      /* "(2) For each vertex s in G except s_f, randomly select
       *      one edge from the s edge list of E(S) to be the
       *      last edge of the s list in a new edge ordering."
       *
       * select random edges and move them to the end of each 
       * edge list.
       */
      for (x = 0; x < 26; x++)
	{
	  if (nE[x] == 0 || x == sf) continue;
	  pos           = esl_rnd_Choose(r, nE[x]);
	  ESL_SWAP(E[x][pos], E[x][nE[x]-1], char);
	}

      /* "(3) From this last set of edges, construct the last-edge
       *      graph Z and determine whether or not all of its
       *      vertices are connected to s_f."
       * 
       * a probably stupid algorithm for looking at the
       * connectivity in Z: iteratively sweep through the
       * edges in Z, and build up an array (confusing called Z[x])
       * whose elements are 1 if x is connected to sf, else 0.
       */
      for (x = 0; x < 26; x++) Z[x] = 0;
      Z[(int) sf] = keep_connecting = 1;

      while (keep_connecting) {
	keep_connecting = 0;
	for (x = 0; x < 26; x++) {
	  if (nE[x] == 0) continue;
	  y = E[x][nE[x]-1];            /* xy is an edge in Z */
	  if (Z[x] == 0 && Z[y] == 1) {  /* x is connected to sf in Z */
	    Z[x] = 1;
	    keep_connecting = 1;
	  }
	}
      }

      /* if any vertex in Z is tagged with a 0, it's
       * not connected to sf, and we won't have a Eulerian
       * walk.
       */
      is_eulerian = 1;
      for (x = 0; x < 26; x++) {
	if (nE[x] == 0 || x == sf) continue;
	if (Z[x] == 0) {
	  is_eulerian = 0;
	  break;
	}
      }

      /* "(4) If any vertex is not connected in Z to s_f, the
       *      new edge ordering will not be Eulerian, so return to
       *      (2). If all vertices are connected in Z to s_f, 
       *      the new edge ordering will be Eulerian, so
       *      continue to (5)."
       *      
       * e.g. note infinite loop while is_eulerian is FALSE.
       */
    }

  /* "(5) For each vertex s in G, randomly permute the remaining
   *      edges of the s edge list of E(S) to generate the s
   *      edge list of the new edge ordering E(S')."
   *      
   * Essentially a StrShuffle() on the remaining nE[x]-1 elements
   * of each edge list; unfortunately our edge lists are arrays,
   * not strings, so we can't just call out to StrShuffle().
   */
  for (x = 0; x < 26; x++)
    for (n = nE[x] - 1; n > 1; n--)
      {
	pos       = esl_rnd_Choose(r, n);
	ESL_SWAP(E[x][pos], E[x][n-1], char);
      }

  /* "(6) Construct sequence S', a random DP permutation of
   *      S, from E(S') as follows. Start at the s_1 edge list.
   *      At each s_i edge list, add s_i to S', delete the
   *      first edge s_i,s_j of the edge list, and move to
   *      the s_j edge list. Continue this process until
   *      all edge lists are exhausted."
   */ 
  pos = 0; 
  x = toupper((int) s[0]) - 'A';
  while (1) 
    {
      shuffled[pos++] = 'A'+ x; /* add s_i to S' */
      
      y = E[x][iE[x]];
      iE[x]++;			/* "delete" s_i,s_j from edge list */
  
      x = y;			/* move to s_j edge list. */

      if (iE[x] == nE[x])
	break;			/* the edge list is exhausted. */
    }
  shuffled[pos++] = 'A' + sf;
  shuffled[pos]   = '\0';  

  /* Reality checks.
   */
  if (x   != sf)  ESL_XEXCEPTION(eslEINCONCEIVABLE, "hey, you didn't end on s_f.");
  if (pos != len) ESL_XEXCEPTION(eslEINCONCEIVABLE, "hey, pos (%d) != len (%d).", pos, len);
  
  /* Free and return.
   */
  esl_Free2D((void **) E, 26);
  free(nE);
  free(iE);
  return eslOK;

 ERROR:
  esl_Free2D((void **) E, 26);
  if (nE != NULL) free(nE);
  if (iE != NULL) free(nE);
  return status;
}


/* Function:  esl_rnd_CMarkov0()
 * Synopsis:  Generate new text string of same 0th order Markov properties.
 * Incept:    SRE, Sat Feb 24 08:47:43 2007 [Casa de Gatos]
 *
 * Purpose:   Makes a random string <markoved> with the same length and
 *            0-th order Markov properties as <s>, given randomness
 *            source <r>.
 *            
 *            <s> and <markoved> can be point to the same storage, in which
 *            case <s> is randomized in place, destroying the original
 *            string.
 *            
 *            <s> must consist only of alphabetic characters [a-zA-Z].
 *            Statistics are collected case-insensitively over 26 possible
 *            residues. The random string is generated all upper case.
 *
 * Args:      s         - input string
 *            markoved  - randomly generated string 
 *                        (storage allocated by caller, at least strlen(s)+1)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains nonalphabetic characters.
 */
int 
esl_rnd_CMarkov0(ESL_RANDOMNESS *r, const char *s, char *markoved)
{
  int    L;
  int    i; 
  double p[26];		/* initially counts, then probabilities */
  int    x;

  /* First, verify that the string is entirely alphabetic. */
  L = strlen(s);
  for (i = 0; i < L; i++)
    if (! isalpha((int) s[i])) 
      ESL_EXCEPTION(eslEINVAL, "String contains nonalphabetic characters");

  /* Collect zeroth order counts and convert to frequencies. 
   */
  for (x = 0; x < 26; x++) p[x] = 0.;
  for (i = 0; i < L; i++)
    p[(int)(toupper((int) s[i]) - 'A')] += 1.0;
  if (L > 0)
    for (x = 0; x < 26; x++) p[x] /= (double) L;

  /* Generate a random string using those p's. */
  for (i = 0; i < L; i++)
    markoved[i] = esl_rnd_DChoose(r, p, 26) + 'A';
  markoved[i] = '\0';

  return eslOK;
}

/* Function:  esl_rnd_CMarkov1()
 * Synopsis:  Generate new text string of same 1st order Markov properties.
 * Incept:    SRE, Sat Feb 24 09:21:46 2007 [Casa de Gatos]
 *
 * Purpose:   Makes a random string <markoved> with the same length and
 *            1st order (di-residue) Markov properties as <s>, given
 *            randomness source <r>.
 *            
 *            <s> and <markoved> can be point to the same storage, in which
 *            case <s> is randomized in place, destroying the original
 *            string.
 *            
 *            <s> must consist only of alphabetic characters [a-zA-Z].
 *            Statistics are collected case-insensitively over 26 possible
 *            residues. The random string is generated all upper case.
 *
 * Args:      s         - input string
 *            markoved  - new randomly generated string 
 *                        (storage allocated by caller, at least strlen(s)+1)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains nonalphabetic characters.
 */
int 
esl_rnd_CMarkov1(ESL_RANDOMNESS *r, const char *s, char *markoved) 
{
  int    L;
  int    i; 
  int    x,y;
  int    i0;			/* initial symbol */
  double p[26][26];		/* conditional probabilities p[x][y] = P(y | x) */
  double p0[26];		/* marginal probabilities P(x), just for initial residue. */

  /* First, verify that the string is entirely alphabetic. */
  L = strlen(s);
  for (i = 0; i < L; i++)
    if (! isalpha((int) s[i])) 
     ESL_EXCEPTION(eslEINVAL, "String contains nonalphabetic characters");

  /* Collect first order counts and convert to frequencies. */
  for (x = 0; x < 26; x++) 
    for (y = 0; y < 26; y++) 
      p[x][y] = 0.;

  i0 = x = toupper((int) s[0]) - 'A';
  for (i = 1; i < L; i++) 
    {
      y = toupper((int) s[i]) - 'A';
      p[x][y] += 1.0;
      x = y;
    }
  p[x][i0] += 1.0; 		/* "circularized": avoids a bug; see markov1_bug utest */

  for (x = 0; x < 26; x++) 
    {
      p0[x] = 0.;
      for (y = 0; y < 26; y++)
	p0[x] += p[x][y];	/* now p0[x] = marginal counts of x, inclusive of 1st residue */

      for (y = 0; y < 26; y++) 
	p[x][y] = (p0[x] > 0. ? p[x][y] / p0[x] : 0.); /* now p[x][y] = P(y | x) */
      
      p0[x] /= (double) L;	/* now p0[x] = marginal P(x) */
    }

  /* Generate a random string using those p's. */
  x = esl_rnd_DChoose(r, p0, 26);
  markoved[0] = x + 'A';
  for (i = 1; i < L; i++)
    {
      y           = esl_rnd_DChoose(r, p[x], 26);
      markoved[i] = y + 'A';
      x           = y;
    } 
  markoved[L] = '\0';

  return eslOK;
}

/* Function:  esl_rnd_CReverse()
 * Synopsis:  Reverse a string.
 * Incept:    SRE, Sat Feb 24 10:06:34 2007 [Casa de Gatos]
 *
 * Purpose:   Returns a reversed version of <s> in <rev>. 
 * 
 *            There are no restrictions on the symbols that <s>
 *            might contain.
 * 
 *            Caller provides storage in <rev> for at least
 *            <(strlen(s)+1)*sizeof(char)>.
 *            
 *            <s> and <rev> can point to the same storage, in which
 *            case <s> is reversed in place.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_rnd_CReverse(const char *s, char *rev)
{
  int  L, i;
  char c;
  
  L = strlen(s);
  for (i = 0; i < L/2; i++)
    {				/* swap ends */
      c          = s[L-i-1];
      rev[L-i-1] = s[i];
      rev[i]     = c;
    }
  if (L%2) { rev[i] = s[i]; } /* don't forget middle residue in odd-length s */
  rev[L] = '\0';
  return eslOK;
}

/* Function: esl_rnd_CShuffleWindows()
 * Synopsis: Shuffle local windows of a text string.
 * Incept:   SRE, Sat Feb 24 10:17:59 2007 [Casa de Gatos]
 * 
 * Purpose:  Given string <s>, shuffle residues in nonoverlapping
 *           windows of width <w>, and put the result in <shuffled>.
 *           See [Pearson88].
 *
 *           <s> and <shuffled> can be identical to shuffle in place.
 * 
 *           Caller provides storage in <shuffled> for at least
 *           <(strlen(s)+1)*sizeof(char)>.
 *
 * Args:     s        - string to shuffle in windows
 *           w        - window size (typically 10 or 20)      
 *           shuffled - allocated space for window-shuffled result.
 *           
 * Return:   <eslOK> on success.
 */
int
esl_rnd_CShuffleWindows(ESL_RANDOMNESS *r, const char *s, int w, char *shuffled)
{
  int  L;
  char c;
  int  i, j, k;

  L = strlen(s);
  if (shuffled != s) strcpy(shuffled, s);
  for (i = 0; i < L; i += w)
    for (j = ESL_MIN(L-1, i+w-1); j > i; j--)
      {
	k             = i + esl_rnd_Choose(r, j-i);
	c             = shuffled[k];  /* semantics of a j,k swap, because we might be shuffling in-place */
	shuffled[k]   = shuffled[j];
	shuffled[j]   = c;
      }
  return eslOK;
}




/* Function:  esl_rnd_XShuffle()
 * Synopsis:  Shuffle a digital sequence.
 * Incept:    SRE, Fri Feb 23 08:24:20 2007 [Casa de Gatos]
 *
 * Purpose:   Given a digital sequence <dsq> of length <L> residues,
 *            shuffle it, and leave the shuffled version in <shuffled>.
 *            
 *            Caller provides allocated storage for <shuffled> for at
 *            least the same length as <dsq>. 
 * 
 *            <shuffled> may also point to the same storage as <dsq>,
 *            in which case <dsq> is shuffled in place.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_rnd_XShuffle(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, ESL_DSQ *shuffled)
{
  int     i;
  ESL_DSQ x;

  if (dsq != shuffled) esl_abc_dsqcpy(dsq, L, shuffled);
  while (L > 1) {
    i           = 1 + esl_rnd_Choose(r, L);
    x           = shuffled[i];
    shuffled[i] = shuffled[L];
    shuffled[L] = x;
    L--;
  }
  return eslOK;
}

/* Function:  esl_rnd_XShuffleDP()
 * Synopsis:  Shuffle a digital sequence, preserving diresidue composition.
 * Incept:    SRE, Fri Feb 23 09:23:47 2007 [Casa de Gatos]
 *
 * Purpose:   Same as <esl_rnd_CShuffleDP()>, except for a digital
 *            sequence <dsq> of length <L>, encoded in a digital alphabet
 *            of <K> residues. 
 *            
 *            <dsq> may only consist of residue codes <0..K-1>; if it
 *            contains gaps, degeneracies, or missing data, pass the alphabet's
 *            <Kp> size, not its canonical <K>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains digital residue codes
 *            outside the range <0..K-1>.
 *            <eslEMEM> on allocation failure.
 */
int
esl_rnd_XShuffleDP(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *shuffled)
{
  int     status;           /* Easel return status code */
  int     i;	            /* a position in dsq or shuffled */
  ESL_DSQ x,y;              /* indices of two characters */
  ESL_DSQ **E  = NULL;      /* edge lists: E[0] is the edge list from vertex A */
  int     *nE  = NULL;      /* lengths of edge lists */
  int     *iE  = NULL;      /* positions in edge lists */
  int      n;	            /* tmp: remaining length of an edge list to be shuffled */
  ESL_DSQ  sf;              /* last character in shuffled */
  ESL_DSQ *Z;               /* connectivity in last edge graph Z */ 
  int      keep_connecting; /* flag used in Z connectivity algorithm */
  int      is_eulerian;	    /* flag used for when we've got a good Z */
  
  /* First, verify that we can deal with all the residues in dsq. */
  for (i = 1; i <= L; i++)
    if (dsq[i] >= K)
      ESL_EXCEPTION(eslEINVAL, "dsq contains unexpected residue codes");

  /* Allocations. */
  ESL_ALLOC(nE, sizeof(int)       * K);  for (x = 0; x < K; x++) nE[x] = 0;
  ESL_ALLOC(E,  sizeof(ESL_DSQ *) * K);  for (x = 0; x < K; x++) E[x]  = NULL;
  ESL_ALLOC(iE, sizeof(int)       * K);  for (x = 0; x < K; x++) iE[x] = 0; 
  ESL_ALLOC(Z,  sizeof(ESL_DSQ)   * K);
  for (x = 0; x < K; x++) 
    ESL_ALLOC(E[x], sizeof(ESL_DSQ) * (L-1));

  /* "(1) Construct the doublet graph G and edge ordering E... */
  x = dsq[1];
  for (i = 2; i <= L; i++) {
    E[x][nE[x]] = dsq[i];
    nE[x]++;
    x = dsq[i];
  }
  
  /* Now we have to find a random Eulerian edge ordering. */
  sf = dsq[L];
  is_eulerian = 0;
  while (! is_eulerian)
    {
      for (x = 0; x < K; x++) {
	if (nE[x] == 0 || x == sf) continue;
	i           = esl_rnd_Choose(r, nE[x]);
	ESL_SWAP(E[x][i], E[x][nE[x]-1], ESL_DSQ);
      }

      for (x = 0; x < K; x++) Z[x] = 0;
      Z[(int) sf] = keep_connecting = 1;
      while (keep_connecting) {
	keep_connecting = 0;
	for (x = 0; x < K; x++) {
	  if (nE[x] == 0) continue;
	  y = E[x][nE[x]-1];            /* xy is an edge in Z */
	  if (Z[x] == 0 && Z[y] == 1) {  /* x is connected to sf in Z */
	    Z[x] = 1;
	    keep_connecting = 1;
	  }
	}
      }

      is_eulerian = 1;
      for (x = 0; x < K; x++) {
	if (nE[x] == 0 || x == sf) continue;
	if (Z[x] == 0) {
	  is_eulerian = 0;
	  break;
	}
      }
    }

  /* "(5) For each vertex s in G, randomly permute... */
  for (x = 0; x < K; x++)
    for (n = nE[x] - 1; n > 1; n--)
      {
	i       = esl_rnd_Choose(r, n);
	ESL_SWAP(E[x][i], E[x][n-1], ESL_DSQ);
      }

  /* "(6) Construct sequence S'... */
  i = 1; 
  x = dsq[1];
  while (1) {
    shuffled[i++] = x; 
    y = E[x][iE[x]++];
    x = y;			
    if (iE[x] == nE[x]) break;
  }
  shuffled[i++] = sf;
  shuffled[i]   = eslDSQ_SENTINEL;
  shuffled[0]   = eslDSQ_SENTINEL;

  /* Reality checks. */
  if (x != sf)   ESL_XEXCEPTION(eslEINCONCEIVABLE, "hey, you didn't end on s_f.");
  if (i != L+1)  ESL_XEXCEPTION(eslEINCONCEIVABLE, "hey, i (%d) overran L+1 (%d).", i, L+1);
  
  esl_Free2D((void **) E, K);
  free(nE);
  free(iE);
  free(Z);
  return eslOK;

 ERROR:
  esl_Free2D((void **) E, K);
  if (nE != NULL) free(nE);
  if (iE != NULL) free(nE);
  if (Z  != NULL) free(Z);
  return status;
}


/* Function:  esl_rnd_XMarkov0()
 * Synopsis:  Generate new digital sequence of same 0th order Markov properties.
 * Incept:    SRE, Sat Feb 24 09:12:32 2007 [Casa de Gatos]
 *
 * Purpose:   Same as <esl_rnd_CMarkov0()>, except for a digital
 *            sequence <dsq> of length <L>, encoded in a digital 
 *            alphabet of <K> residues; caller provides storage
 *            for the randomized sequence <markoved> for at least 
 *            <L+2> <ESL_DSQ> residues, including the two flanking
 *            sentinel bytes.
 *            
 *            <dsq> therefore may only consist of residue codes
 *            in the range <0..K-1>. If it contains gaps,
 *            degeneracies, or missing data, pass the alphabet's
 *            <Kp> size, not its canonical <K>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains digital residue codes outside
 *            the range <0..K-1>.
 *            <eslEMEM> on allocation failure.
 */
int 
esl_rnd_XMarkov0(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved)
{
  int     status;
  int     i; 
  double *p = NULL;	/* initially counts, then probabilities */
  int     x;

  /* First, verify that the string is entirely alphabetic. */
  for (i = 1; i <= L; i++)
    if (dsq[i] >= K)
      ESL_XEXCEPTION(eslEINVAL, "String contains unexpected residue codes");

  ESL_ALLOC(p, sizeof(double) * K);
  for (x = 0; x < K; x++) p[x] = 0.;

  for (i = 1; i <= L; i++)
    p[(int) dsq[i]] += 1.0;
  if (L > 0)
    for (x = 0; x < K; x++) p[x] /= (double) L;

  for (i = 1; i <= L; i++)
    markoved[i] = esl_rnd_DChoose(r, p, K);
  markoved[0]   = eslDSQ_SENTINEL;
  markoved[L+1] = eslDSQ_SENTINEL;

  free(p);
  return eslOK;

 ERROR:
  if (p != NULL) free(p);
  return status;
}



/* Function:  esl_rnd_XMarkov1()
 * Synopsis:  Generate new digital sequence of same 1st order Markov properties.
 * Incept:    SRE, Sat Feb 24 09:46:09 2007 [Casa de Gatos]
 *
 * Purpose:   Same as <esl_rnd_CMarkov1()>, except for a digital
 *            sequence <dsq> of length <L>, encoded in a digital 
 *            alphabet of <K> residues. Caller provides storage
 *            for the randomized sequence <markoved> for at least 
 *            <L+2> <ESL_DSQ> residues, including the two flanking
 *            sentinel bytes.
 *            
 *            <dsq> and <markoved> can be point to the same storage, in which
 *            case <dsq> is randomized in place, destroying the original
 *            string.
 *            
 *            <dsq> therefore may only consist of residue codes
 *            in the range <0..K-1>. If it contains gaps,
 *            degeneracies, or missing data, pass the alphabet's
 *            <Kp> size, not its canonical <K>.
 *
 * Args:      dsq       - input digital sequence 1..L
 *            L         - length of dsq
 *            K         - residue codes in dsq are in range 0..K-1
 *            markoved  - new randomly generated digital sequence;
 *                        storage allocated by caller, at least (L+2)*ESL_DSQ;
 *                        may be same as dsq to randomize in place.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <s> contains digital residue codes outside
 *            the range <0..K-1>.
 *            <eslEMEM> on allocation failure.
 */
int 
esl_rnd_XMarkov1(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved) 
{
  int      status;
  int      i; 
  ESL_DSQ  x,y;
  ESL_DSQ  i0;		/* initial symbol */
  double **p;		/* conditional probabilities p[x][y] = P(y | x) */
  double  *p0;		/* marginal probabilities P(x), just for initial residue. */

  /* validate the input string */
  for (i = 1; i <= L; i++)
    if (dsq[i] >= K)
      ESL_XEXCEPTION(eslEINVAL, "String contains unexpected residue codes");

  /* allocations */
  ESL_ALLOC(p0, sizeof(double)   * K);  for (x = 0; x < K; x++) p0[x] = 0.;
  ESL_ALLOC(p,  sizeof(double *) * K);  for (x = 0; x < K; x++) p[x]  = NULL;
  for (x = 0; x < K; x++)
    { ESL_ALLOC(p[x], sizeof(double) * K); for (y = 0; y < K; y++) p[x][y] = 0.; }
  
  /* Collect first order counts and convert to frequencies. */
  i0 = x = dsq[1];
  for (i = 2; i <= L; i++) 
    {
      y = dsq[i];
      p[x][y] += 1.0;
      x = y;
    }
  p[x][i0] += 1.0;	/* "circularized": avoids a bug; see markov1_bug utest */

  for (x = 0; x < K; x++) 
    {
      p0[x] = 0.;
      for (y = 0; y < K; y++)
	p0[x] += p[x][y];	/* now p0[x] = marginal counts of x, inclusive of 1st residue */

      for (y = 0; y < K; y++) 
	p[x][y] = (p0[x] > 0. ? p[x][y] / p0[x] : 0.);	/* now p[x][y] = P(y | x) */
      
      p0[x] /= (double) L;	/* now p0[x] = marginal P(x) inclusive of 1st residue */
    }

  /* Generate a random string using those p's. */
  markoved[1] = esl_rnd_DChoose(r, p0, K);
  for (i = 2; i <= L; i++)
    markoved[i] = esl_rnd_DChoose(r, p[markoved[i-1]], K);

  markoved[0]   = eslDSQ_SENTINEL;
  markoved[L+1] = eslDSQ_SENTINEL;

  esl_Free2D((void**)p, K);
  free(p0);
  return eslOK;

 ERROR:
  esl_Free2D((void**)p, K);
  if (p0 != NULL) free(p0);
  return status;
}


/* Function:  esl_rnd_XReverse()
 * Synopsis:  Reverse a digital sequence.
 * Incept:    SRE, Sat Feb 24 10:13:30 2007 [Casa de Gatos]
 *
 * Purpose:   Given a digital sequence <dsq> of length <L>, return
 *            reversed version of it in <rev>. 
 * 
 *            Caller provides storage in <rev> for at least
 *            <(L+2)*sizeof(ESL_DSQ)>.
 *            
 *            <s> and <rev> can point to the same storage, in which
 *            case <s> is reversed in place.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_rnd_XReverse(const ESL_DSQ *dsq, int L, ESL_DSQ *rev)
{
  int     i;
  ESL_DSQ x;
  
  for (i = 1; i <= L/2; i++)
    {				/* swap ends */
      x          = dsq[L-i+1];
      rev[L-i+1] = dsq[i];
      rev[i]     = x;
    }
  if (L%2) { rev[i] = dsq[i]; } /* don't forget middle residue in odd-length dsq */
  rev[0]   = eslDSQ_SENTINEL;
  rev[L+1] = eslDSQ_SENTINEL;
  return eslOK;
}


/* Function: esl_rnd_XShuffleWindows()
 * Synopsis: Shuffle local windows of a digital sequence.
 * Incept:   SRE, Sat Feb 24 10:51:31 2007 [Casa de Gatos]
 * 
 * Purpose:  Given a digital sequence <dsq> of length <L>, shuffle
 *           residues in nonoverlapping windows of width <w>, and put
 *           the result in <shuffled>.  See [Pearson88].
 *
 *           Caller provides storage in <shuffled> for at least
 *           <(L+2)*sizeof(ESL_DSQ)>.
 *           
 *           <dsq> and <shuffled> can be identical to shuffle in place.
 *
 * Args:     dsq      - digital sequence to shuffle in windows
 *           L        - length of <dsq>
 *           w        - window size (typically 10 or 20)      
 *           shuffled - allocated space for window-shuffled result.
 *           
 * Return:   <eslOK> on success.
 */
int
esl_rnd_XShuffleWindows(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int w, ESL_DSQ *shuffled)
{
  ESL_DSQ x;
  int  i, j, k;

  if (dsq != shuffled) esl_abc_dsqcpy(dsq, L, shuffled);
  for (i = 1; i <= L; i += w)
    for (j = ESL_MIN(L, i+w-1); j > i; j--)
      {
	k           = i + esl_rnd_Choose(r, j-i+1);
	x           = shuffled[k];  /* semantics of a j,k swap, because we might be shuffling in-place */
	shuffled[k] = shuffled[j];
	shuffled[j] = x;
      }
  return eslOK;
}




/*****************************************************************
 * 7. Unit tests.
 *****************************************************************/

#ifdef eslRANDOM_TESTDRIVE
#include <esl_alphabet.h>
#include <esl_vectorops.h>
#include <esl_stats.h>
#include <esl_dirichlet.h>
    
  
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
    counts[esl_rnd_Choose(r, nbins)]++;

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
 
/* count c(x) monoresidue and c(xy) diresidue composition
 * used for sequence shuffling unit tests
 * mono, di allocated by caller for 26 and 26x26, respectively.
 */
static int
composition(char *s, int L, int *mono, int **di)
{
  int i, x, y;

  for (x = 0; x < 26; x++) {
    mono[x] = 0;
    for (y = 0; y < 26; y++)
      di[x][y] = 0;
  }

  for (i = 0; s[i] != '\0'; i++) { 
    if (!isalpha(s[i])) esl_fatal("bad residue %d", i);
    y = toupper(s[i]) - 'A';
    mono[y]++;
    if (i > 0) {
      x = toupper(s[i-1] - 'A');
      di[x][y]++;
    }
  }
  if (i != L) esl_fatal("sequence length didn't match expected %d", L);
  return eslOK;
}

/* same, but for digital seq., with alphabet size K */
static int
xcomposition(ESL_DSQ *dsq, int L, int K, int *mono, int **di)
{
  int i, x, y;

  for (x = 0; x < K; x++) {
    mono[x] = 0;
    for (y = 0; y < K; y++)
      di[x][y] = 0;
  }

  for (i = 1; dsq[i] != eslDSQ_SENTINEL; i++) { 
    if (dsq[i] > K) esl_fatal("bad residue %d", i);
    if (i > 1) di[(int) dsq[i-1]][(int) dsq[i]]++;
    mono[(int) dsq[i]]++;
  }
  if (i != L+1) esl_fatal("sequence length didn't match expected %d", L);
  return eslOK;
}

static int
composition_allocate(int K, int **ret_mono, int ***ret_di)
{
  int  status;
  int *mono = NULL;
  int **di  = NULL;
  int  x;

  ESL_ALLOC(mono, sizeof(int)   * K);
  ESL_ALLOC(di,   sizeof(int *) * K); for (x = 0; x < K; x++) di[x] = NULL;
  for (x = 0; x < K; x++)
    ESL_ALLOC(di[x], sizeof(int) * K);
  *ret_mono = mono;
  *ret_di   = di;
  return eslOK;

 ERROR:
  esl_Free2D((void **) di, K);
  if (mono != NULL) free(mono);
  *ret_mono = NULL;
  *ret_di   = NULL;
  return status;
}

/* compare compositions before/after.
 * either mono (m1,m2) or di (d1,d2) may be NULL, to compare only the other one */
static int
composition_compare(int *m1, int **di1, int *m2, int **di2, int K)
{
  int x,y;

  for (x = 0; x < K; x++) {
    if (m1 != NULL && m1[x] != m2[x]) return eslFAIL;
    if (di1 != NULL) 
      for (y = 0; y < K; y++) 
	if (di1[x][y] != di2[x][y])   return eslFAIL;
  }
  return eslOK;
}

/* Unit tests for:
 *     esl_rnd_CShuffle()
 *     esl_rnd_CShuffleDP()
 *     esl_rnd_CShuffleWindows()
 *     esl_rnd_CReverse()
 * 
 * All of these exactly preserve residue composition, which is
 * the basis of the unit tests.
 */
static void
utest_CShufflers(ESL_RANDOMNESS *r, int L, char *alphabet, int K)
{
  char   *logmsg  = "Failure in one of the CShuffle* unit tests";
  int     status;
  char   *s   = NULL;
  char   *s2  = NULL;
  int    *m1  = NULL,
         *m2  = NULL;	    /* mono, before and after */
  int   **di1 = NULL,
        **di2 = NULL;       /* di, before and after */
  double  *p;		    
  int      w = 12;   	    /* window width for CShuffleWindows() */

  /* allocations */
  ESL_ALLOC(s,   sizeof(char)   * (L+1));
  ESL_ALLOC(s2,  sizeof(char)   * (L+1));
  ESL_ALLOC(p,   sizeof(double) * K);
  if (composition_allocate(26, &m1, &di1) != eslOK) esl_fatal(logmsg);
  if (composition_allocate(26, &m2, &di2) != eslOK) esl_fatal(logmsg);

  /* generate the string we'll start shuffling */
  if (esl_dirichlet_DSampleUniform(r, K, p) != eslOK) esl_fatal(logmsg);
  if (esl_rnd_IID(r, alphabet, p, K, L, s)  != eslOK) esl_fatal(logmsg);

  /* esl_rnd_CShuffle: mono composition should stay exactly the same, di may change */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s,   L, m1, di1)                != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffle(r, s, s2)                  != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, 26) != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 

  /* esl_rnd_CShuffle, in place */
  strcpy(s, s2);
  if (composition(s2, L, m1, di1)                 != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffle(r, s2, s2)                 != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, 26) != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 

  /* esl_rnd_CShuffleDP: mono and di compositions stay exactly the same */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s, L, m1,  di1)                 != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffleDP(r, s, s2)                != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, di1, m2, di2, 26)   != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 

  /* esl_rnd_CShuffleDP, in place */
  strcpy(s, s2);
  if (composition(s2, L, m1, di1)                 != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffleDP(r, s2, s2)               != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, di1, m2, di2, 26)   != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 
  
  /* esl_rnd_CShuffleWindows(): mono composition stays the same */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s,   L, m1, di1)                != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffleWindows(r, s, w, s2)        != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, 26) != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 
  
  /* esl_rnd_CShuffleWindows(), in place */
  strcpy(s, s2);
  if (composition(s2, L, m1, di1)                 != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CShuffleWindows(r, s2, w, s2)       != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, 26) != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 
  
  /* esl_rnd_CReverse(): two reverses (one in place) give the same seq back */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s,   L, m1, di1)                != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CReverse(s, s2)                     != eslOK) esl_fatal(logmsg);      
  if (composition(s2, L, m2, di2)                 != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, 26) != eslOK) esl_fatal(logmsg);
  if (strcmp(s2, s) == 0)                                   esl_fatal(logmsg); 
  if (esl_rnd_CReverse(s2, s2)                    != eslOK) esl_fatal(logmsg);      
  if (strcmp(s2, s) != 0)                                   esl_fatal(logmsg); 

  free(s);
  free(s2);
  free(p);
  free(m1);
  free(m2);
  esl_Free2D((void **) di1, 26);
  esl_Free2D((void **) di2, 26);
  return;
  
 ERROR:
  esl_fatal(logmsg);
}

/* Unit tests for:
 *    esl_rnd_CMarkov0()
 *    esl_rnd_CMarkov1()
 * 
 * Testing these is less robust than the shufflers, because it's hard
 * to concoct deterministic tests. Instead the test is a weak one,
 * that zero probability events get zero counts.
 */
static void
utest_CMarkovs(ESL_RANDOMNESS *r, int L, char *alphabet)
{
  char   *logmsg = "Failure in a CMarkov*() unit test";
  int     status;
  char   *s   = NULL;
  char   *s2  = NULL;
  float  *p   = NULL;
  int     K;
  int     pzero;
  int    *m1  = NULL,
         *m2  = NULL;	    /* mono, before and after */
  int   **di1 = NULL,
        **di2 = NULL;       /* di, before and after */
  int     i,x;

  K = strlen(alphabet);
  ESL_ALLOC(p,   sizeof(float)  * K);
  ESL_ALLOC(s,   sizeof(char)   * (L+1));
  ESL_ALLOC(s2,  sizeof(char)   * (L+1));
  if (composition_allocate(26, &m1, &di1) != eslOK) esl_fatal(logmsg);
  if (composition_allocate(26, &m2, &di2) != eslOK) esl_fatal(logmsg);

  /* generate string with a random letter prob set to 0  */
  pzero = esl_rnd_Choose(r, K);
  if (esl_dirichlet_FSampleUniform(r, K, p)  != eslOK) esl_fatal(logmsg);
  p[pzero] = 0;
  esl_vec_FNorm(p, K);
  if (esl_rnd_fIID(r, alphabet, p, K, L, s)  != eslOK) esl_fatal(logmsg);

  /* esl_rnd_CMarkov0()  */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s,   L, m1, di1)  != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CMarkov0(r, s, s2)    != eslOK) esl_fatal(logmsg);
  if (composition(s2, L, m2, di2)   != eslOK) esl_fatal(logmsg);  
  if (m1[pzero]                     != 0)     esl_fatal(logmsg);  
  if (m2[pzero]                     != 0)     esl_fatal(logmsg);  
  if (strcmp(s2, s)                 == 0)     esl_fatal(logmsg);  
  
  /* esl_rnd_CMarkov0(), in place */
  strcpy(s, s2);
  if (esl_rnd_CMarkov0(r, s2, s2)   != eslOK) esl_fatal(logmsg);
  if (composition(s2, L, m2, di2)   != eslOK) esl_fatal(logmsg);  
  if (m2[pzero]                     != 0)     esl_fatal(logmsg);  
  if (strcmp(s2, s)                 == 0)     esl_fatal(logmsg);  
  
  /* generate string with all homodiresidues set to 0 */
  if (esl_dirichlet_FSampleUniform(r, K, p)  != eslOK) esl_fatal(logmsg);
  do {
    if (esl_rnd_fIID(r, alphabet, p, K, L, s)  != eslOK) esl_fatal(logmsg);  
    for (i = 1; i < L; i++)
      if (s[i] == s[i-1]) /* this incantation will rotate letter forward in alphabet: */
	s[i] = alphabet[(1+strchr(alphabet,s[i])-alphabet)%K];
  } while (s[0] == s[L-1]);	/* lazy: reject strings where circularization would count a homodimer */
  
  /* esl_rnd_CMarkov1()  */
  memset(s2, 0, (L+1)*sizeof(char));
  if (composition(s,   L, m1, di1)  != eslOK) esl_fatal(logmsg);
  if (esl_rnd_CMarkov1(r, s, s2)    != eslOK) esl_fatal(logmsg);
  if (composition(s2, L, m2, di2)   != eslOK) esl_fatal(logmsg);  
  for (x = 0; x < K; x++) {
    if (di1[x][x]                   != 0)     esl_fatal(logmsg);  
    if (di2[x][x]                   != 0)     esl_fatal(logmsg);  
  }
  if (strcmp(s2, s)                 == 0)     esl_fatal(logmsg);  

  /* esl_rnd_CMarkov1(), in place  */
  strcpy(s, s2);
  if (esl_rnd_CMarkov1(r, s2, s2)  != eslOK)   esl_fatal(logmsg);
  if (composition(s2, L, m2, di2)  != eslOK) esl_fatal(logmsg);  
  for (x = 0; x < K; x++) {
    if (di1[x][x]                   != 0)     esl_fatal(logmsg);  
    if (di2[x][x]                   != 0)     esl_fatal(logmsg);  
  }
  if (strcmp(s2, s)                 == 0)     esl_fatal(logmsg);  
  
  free(s);
  free(s2);
  free(p);
  free(m1);
  free(m2);
  esl_Free2D((void **) di1, 26);
  esl_Free2D((void **) di2, 26);
  return;
  
 ERROR:
  esl_fatal(logmsg);
}


/* Unit tests for:
 *     esl_rnd_XShuffle()
 *     esl_rnd_XShuffleDP()
 *     esl_rnd_XShuffleWindows()
 *     esl_rnd_XReverse()
 * Same ideas as testing the C* versions, adapted for digital sequences. 
 */
static void
utest_XShufflers(ESL_RANDOMNESS *r, int L, int K)
{
  char    *logmsg  = "Failure in one of the XShuffle* unit tests";
  int      status;
  ESL_DSQ *dsq   = NULL;
  ESL_DSQ *ds2   = NULL;
  int     *m1    = NULL,
          *m2    = NULL;    /* mono, before and after */
  int    **di1   = NULL,
         **di2   = NULL;    /* di, before and after */
  float   *p     = NULL;
  int      w = 12;   	    /* window width for XShuffleWindows() */

  /* allocations */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  ESL_ALLOC(ds2, sizeof(ESL_DSQ) * (L+2));
  ESL_ALLOC(p,   sizeof(double)  * K);
  if (composition_allocate(K, &m1, &di1) != eslOK) esl_fatal(logmsg);
  if (composition_allocate(K, &m2, &di2) != eslOK) esl_fatal(logmsg);

  /* generate the string we'll test shuffling on, keep its composition stats */
  if (esl_dirichlet_FSampleUniform(r, K, p) != eslOK) esl_fatal(logmsg);
  if (esl_rnd_xfIID(r, p, K, L, dsq)        != eslOK) esl_fatal(logmsg);

  /* esl_rnd_XShuffle: mono composition should stay exactly the same, di may change */
  memset(ds2, eslDSQ_SENTINEL, (L+2));
  if (xcomposition(dsq, L, K, m1, di1)           != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffle(r, dsq, L, ds2)           != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, K) != eslOK) esl_fatal(logmsg);

  /* esl_rnd_XShuffle, in place */
  if (esl_abc_dsqcpy(ds2, L, dsq)                != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m1,  di1)          != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffle(r, ds2, L, ds2)           != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, K) != eslOK) esl_fatal(logmsg);

  /* esl_rnd_XShuffleDP: mono and di compositions stay exactly the same */
  memset(ds2, eslDSQ_SENTINEL, (L+2));
  if (xcomposition(dsq, L, K, m1,  di1)          != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffleDP(r, dsq, L, K, ds2)      != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, di1, m2, di2, K)   != eslOK) esl_fatal(logmsg);

  /* esl_rnd_XShuffleDP, in place */
  if (esl_abc_dsqcpy(ds2, L, dsq)                != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m1, di1)           != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffleDP(r, ds2, L, K, ds2)      != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, di1, m2, di2, K)   != eslOK) esl_fatal(logmsg);
  
  /* esl_rnd_XShuffleWindows(): mono composition stays the same */
  memset(ds2, eslDSQ_SENTINEL, (L+2));
  if (xcomposition(dsq, L, K, m1, di1)           != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffleWindows(r, dsq, L, w, ds2) != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, K) != eslOK) esl_fatal(logmsg);
  
  /* esl_rnd_XShuffleWindows(), in place */
  if (esl_abc_dsqcpy(ds2, L, dsq)                != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m1,  di1)          != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XShuffleWindows(r, ds2, L, w, ds2) != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)           != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, K) != eslOK) esl_fatal(logmsg);
  
  /* esl_rnd_XReverse(): two reverses (one in place) give the same seq back */
  memset(ds2, eslDSQ_SENTINEL, (L+2));
  if (xcomposition(dsq, L, K, m1, di1)            != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XReverse(dsq, L, ds2)               != eslOK) esl_fatal(logmsg);      
  if (xcomposition(ds2, L, K, m2, di2)            != eslOK) esl_fatal(logmsg);
  if (composition_compare(m1, NULL, m2, NULL, K)  != eslOK) esl_fatal(logmsg);
  if (memcmp((void *) ds2, (void *) dsq, sizeof(ESL_DSQ)*(L+2)) == 0) esl_fatal(logmsg); 
  if (esl_rnd_XReverse(ds2, L, ds2)               != eslOK) esl_fatal(logmsg);      
  if (memcmp((void *) ds2, (void *) dsq, sizeof(ESL_DSQ)*(L+2)) != 0) esl_fatal(logmsg); 

  free(dsq);
  free(ds2);
  free(p);
  free(m1);
  free(m2);
  esl_Free2D((void **) di1, K);
  esl_Free2D((void **) di2, K);
  return;
  
 ERROR:
  esl_fatal(logmsg);
}

/* Unit tests for:
 *    esl_rnd_XMarkov0()
 *    esl_rnd_XMarkov1()
 * Same ideas as in the C* versions, but for digital sequences.
 */
static void
utest_XMarkovs(ESL_RANDOMNESS *r, int L, int K)
{
  char    *logmsg = "Failure in an XMarkov*() unit test";
  int      status;
  ESL_DSQ *dsq = NULL;
  ESL_DSQ *ds2 = NULL;
  int     *m1  = NULL, 
          *m2  = NULL;    /* mono, before and after */
  int    **di1 = NULL,
         **di2 = NULL;    /* di, before and after */
  float   *p   = NULL;
  int      pzero;
  int      i,x;

  /* allocations */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  ESL_ALLOC(ds2, sizeof(ESL_DSQ) * (L+2));
  ESL_ALLOC(p,   sizeof(double)  * K);
  if (composition_allocate(K, &m1, &di1) != eslOK) esl_fatal(logmsg);
  if (composition_allocate(K, &m2, &di2) != eslOK) esl_fatal(logmsg);

  /* generate sequence with a random letter prob set to 0  */
  pzero = esl_rnd_Choose(r, K);
  if (esl_dirichlet_FSampleUniform(r, K, p)  != eslOK) esl_fatal(logmsg);
  p[pzero] = 0.;
  esl_vec_FNorm(p, K);
  if (esl_rnd_xfIID(r, p, K, L, dsq)         != eslOK) esl_fatal(logmsg);

  /* esl_rnd_XMarkov0()  */
  memset(ds2, eslDSQ_SENTINEL, (L+2)*sizeof(ESL_DSQ));
  if (xcomposition(dsq, L, K, m1, di1)        != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XMarkov0(r, dsq, L, K, ds2)     != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m2, di2)        != eslOK) esl_fatal(logmsg);  
  if (m1[pzero]                               != 0)     esl_fatal(logmsg);  
  if (m2[pzero]                               != 0)     esl_fatal(logmsg);  
  if (memcmp(ds2, dsq, sizeof(ESL_DSQ)*(L+2)) == 0)     esl_fatal(logmsg);  
  
  /* esl_rnd_CMarkov0(), in place */
  if (esl_abc_dsqcpy(ds2, L, dsq)             != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XMarkov0(r, ds2, L, K, ds2)     != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m2, di2)        != eslOK) esl_fatal(logmsg);  
  if (m2[pzero]                               != 0)     esl_fatal(logmsg);  
  if (memcmp(ds2, dsq, sizeof(ESL_DSQ)*(L+2)) == 0)     esl_fatal(logmsg);  
  
  /* generate string with all homodiresidues set to 0 */
  if (esl_dirichlet_FSampleUniform(r, K, p)   != eslOK) esl_fatal(logmsg);
  do {
    if (esl_rnd_xfIID(r, p, K, L, dsq)          != eslOK) esl_fatal(logmsg);  
    for (i = 2; i <= L; i++)
      if (dsq[i] == dsq[i-1]) /* this incantation will rotate letter forward in alphabet: */
	dsq[i] = (dsq[i]+1)%K;
  } while (dsq[1] == dsq[L]);	/* lazy. reject strings where circularization would count a homodimer */
    
  /* esl_rnd_XMarkov1()  */
  memset(ds2, eslDSQ_SENTINEL, (L+2)*sizeof(ESL_DSQ));
  if (xcomposition(dsq, L, K, m1, di1)        != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XMarkov1(r, dsq, L, K, ds2)     != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m2, di2)        != eslOK) esl_fatal(logmsg);  
  for (x = 0; x < K; x++) {
    if (di1[x][x]                             != 0)     esl_fatal(logmsg);  
    if (di2[x][x]                             != 0)     esl_fatal(logmsg);  
  }
  if (memcmp(ds2, dsq, sizeof(ESL_DSQ)*(L+2)) == 0)     esl_fatal(logmsg);  

  /* esl_rnd_XMarkov1(), in place  */
  if (esl_abc_dsqcpy(ds2, L, dsq)             != eslOK) esl_fatal(logmsg);
  if (esl_rnd_XMarkov1(r, ds2, L, K, ds2)     != eslOK) esl_fatal(logmsg);
  if (xcomposition(ds2, L, K, m2, di2)        != eslOK) esl_fatal(logmsg);  
  for (x = 0; x < K; x++) {
    if (di1[x][x]                             != 0)     esl_fatal(logmsg);  
    if (di2[x][x]                             != 0)     esl_fatal(logmsg);  
  }
  if (memcmp(ds2, dsq, sizeof(ESL_DSQ)*(L+2)) == 0)     esl_fatal(logmsg);  
  
  free(dsq);
  free(ds2);
  free(p);
  free(m1);
  free(m2);
  esl_Free2D((void **) di1, K);
  esl_Free2D((void **) di2, K);
  return;
  
 ERROR:
  esl_fatal(logmsg);
}

/* utest_markov1_bug()
 * 
 * Given a sequence like AAAAAAAAAT, where a residue only occurs once
 * and at the end of the sequence, a bug can appear: a Markov chain
 * can transit to T, but can't leave. Easel handles this by 
 * counting Markov statistics as if the input sequence were circular.
 */
static void
utest_markov1_bug(ESL_RANDOMNESS *r)
{
  char    logmsg[]  = "Failure in markov1_bug test (zero/absorbing transition)";
  char    testseq[] = "AAAAAAAAAT";
  char   *seq       = NULL;
  ESL_DSQ testdsq[] = { eslDSQ_SENTINEL,0,0,0,0,0,0,0,0,0,3,eslDSQ_SENTINEL};
  ESL_DSQ *dsq      = NULL;
  int     L         = strlen(testseq);
  int    *mono      = NULL;
  int   **di        = NULL;
  int     N         = 100;         
  int     i;

  if ((seq = malloc(sizeof(char)    * (L+1))) == NULL)    esl_fatal(logmsg);
  if ((dsq = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL)    esl_fatal(logmsg);

  if (composition_allocate(4, &mono, &di)       != eslOK) esl_fatal(logmsg);
  for (i = 0; i < N; i++) {
    if (esl_rnd_XMarkov1(r, testdsq, L, 4, dsq) != eslOK) esl_fatal(logmsg);
    if (xcomposition(testdsq, L, 4, mono, di)   != eslOK) esl_fatal(logmsg);
    if (mono[0] + mono[3] != L)                           esl_fatal(logmsg);
  }
  esl_Free2D((void **) di, 4);
  free(mono);

  if (composition_allocate(26, &mono, &di) != eslOK) esl_fatal(logmsg);
  for (i = 0; i < N; i++) {
    if (esl_rnd_CMarkov1(r, testseq, seq)  != eslOK) esl_fatal(logmsg);
    if (composition(seq, L, mono, di)      != eslOK) esl_fatal(logmsg);
    if (mono[0] + mono['T'-'A'] != L)                esl_fatal(logmsg);
  }
  esl_Free2D((void **) di, 26);
  free(mono);
  free(seq);
  free(dsq);
}

#endif /*eslRANDOM_TESTDRIVE*/


/*****************************************************************
 * 8. The test driver.
 *****************************************************************/

/* gcc -g -Wall -o testdrive -L. -I. -DeslRANDOM_TESTDRIVE esl_random.c -leasel -lm
 */
#ifdef eslRANDOM_TESTDRIVE

#include <stdio.h>
#include "easel.h"
#include "esl_getopts.h"
#include "esl_dirichlet.h"
#include "esl_vectorops.h"
#include "esl_random.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  {"-b",  eslARG_INT,      "20", NULL, "n>0",NULL, NULL, NULL, "number of test bins",               0},
  {"-n",  eslARG_INT, "1000000", NULL, "n>0",NULL, NULL, NULL, "number of samples",                 0},
  {"-s",  eslARG_INT,      "42", NULL, "n>0",NULL, NULL, NULL, "random number seed",                0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose output",               0},
  {"-L",  eslARG_INT,    "1000", NULL, NULL, NULL, NULL, NULL, "length of random shuffled seqs",    0},
  {"--bitfile",eslARG_STRING,NULL,NULL,NULL, NULL, NULL, NULL, "save bit file for NIST benchmark",  0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[] = "Usage: ./testdrive [-options]";

static int save_bitfile(char *bitfile, ESL_RANDOMNESS *r, int n);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go = NULL;
  ESL_RANDOMNESS *r  = NULL;
  int             n, nbins, seed, be_verbose;
  char           *bitfile;
  int             L;
  char           *alphabet = "ACGT";
  int             K;

  K = strlen(alphabet);

  /* Command line parsing
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("%s", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("%s", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage); 
    puts("\n  where options are:");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2=indentation; 80=width */
    return 0;
  }
  nbins      = esl_opt_GetInteger(go, "-b");
  n          = esl_opt_GetInteger(go, "-n");
  seed       = esl_opt_GetInteger(go, "-s");
  be_verbose = esl_opt_GetBoolean(go, "-v");
  L          = esl_opt_GetInteger(go, "-L");
  bitfile    = esl_opt_GetString (go, "--bitfile");
  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return 1;
  }

  /* Initialization
   */
  if ((r = esl_randomness_Create(seed)) == NULL) esl_fatal("randomness creation failed");

  /* Unit tests
   */
  utest_random(seed, n, nbins, be_verbose);
  utest_choose(r,    n, nbins, be_verbose);

  utest_CShufflers(r, L, alphabet, K);
  utest_CMarkovs  (r, L, alphabet);
  utest_XShufflers(r, L, K);
  utest_XMarkovs  (r, L, K);

  utest_markov1_bug(r);

  /* Optional datafiles.
   */
  if (bitfile != NULL) save_bitfile(bitfile, r, n);

  /* Exit.
   */
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
 * 9. An example of using the random module.
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


