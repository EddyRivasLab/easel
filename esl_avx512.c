/* Vectorized routines for x86 AVX-512 instructions
 * 
 * Most speed-critical code is in the .h file, to facilitate inlining.
 * 
 * Contents:
 *    1. Vector expf() and logf() functions
 *    2. Debugging/development routines
 *    3. Unit tests
 *    4. Test driver
 *    
 * This code is conditionally compiled, only when <eslENABLE_AVX512> was
 * set in <esl_config.h> by the configure script, and that will only
 * happen on x86 platforms. When <eslENABLE_AVX512> is not set, we
 * include some dummy code to silence compiler and ranlib warnings
 * about empty translation units and no symbols, and dummy drivers
 * that do nothing but declare success.
 */
#include <esl_config.h>
#ifdef eslENABLE_AVX512

#include <stdio.h>
#include <x86intrin.h>		

#include "easel.h"
#include "esl_avx512.h"

/*****************************************************************
 * 1. AVX512 SIMD logf(), expf()
 *****************************************************************/ 

/* Function:  esl_avx512_logf()
 * Synopsis:  <r[z] = log x[z]>
 * Incept:    NPC, April 29 2024
 *
 * Purpose:   Given a vector <x> containing sixteen floats, returns a
 *            vector <r> in which each element <r[z] = logf(x[z])>.
 *            
 *            Valid in the domain $x_z > 0$ for normalized IEEE754
 *            $x_z$.
 *
 *            For <x> $< 0$, including -0, returns <NaN>. For <x> $==
 *            0$ or subnormal <x>, returns <-inf>. For <x = inf>,
 *            returns <inf>. For <x = NaN>, returns <NaN>. For 
 *            subnormal <x>, returns <-inf>.
 *
 * Xref:      J2/71.
 * 
 * Note:      Derived from an SSE1 implementation by Julian
 *            Pommier. Converted to SSE2 and added handling
 *            of IEEE754 specials.  Extended to AVX512 by Nick Carter
 */
__m512 
esl_avx512_logf(__m512 x) 
{
  static float cephes_p[9] = {  7.0376836292E-2f, -1.1514610310E-1f,  1.1676998740E-1f,
				-1.2420140846E-1f, 1.4249322787E-1f, -1.6668057665E-1f,
				2.0000714765E-1f, -2.4999993993E-1f,  3.3333331174E-1f };
  __m512  onev = _mm512_set1_ps(1.0f);          /* all elem = 1.0 */
  __m512  v0p5 = _mm512_set1_ps(0.5f);          /* all elem = 0.5 */
  __m512i vneg = _mm512_set1_epi32(0x80000000); /* all elem have IEEE sign bit up */
  __m512i vexp = _mm512_set1_epi32(0x7f800000); /* all elem have IEEE exponent bits up */
  __m512i ei;
  __m512  e;
  __mmask16  invalid_mask, zero_mask, inf_mask;            /* masks used to handle special IEEE754 inputs */
  __mmask16  mask;
  __m512  origx;
  __m512  tmp;
  __m512  y;
  __m512  z;

  /* first, split x apart: x = frexpf(x, &e); */
  ei           = _mm512_srli_epi32( _mm512_castps_si512(x), 23);	                                        /* shift right 23: IEEE754 floats: ei = biased exponents     */
  invalid_mask = _mm512_cmp_epi32_mask( _mm512_and_si512(_mm512_castps_si512(x), vneg), vneg, 0);  /* mask any elem that's negative; these become NaN           */
  // 0 = compare for equality
  zero_mask    = _mm512_cmp_epi32_mask(ei, _mm512_setzero_si512(), 0);                          /* mask any elem zero or subnormal; these become -inf        */
  inf_mask     = _mm512_cmp_epi32_mask( _mm512_and_si512(_mm512_castps_si512(x), vexp), vexp,0);  /* mask any elem inf or NaN; log(inf)=inf, log(NaN)=NaN      */
  origx        = x;			                                                                /* store original x, used for log(inf) = inf, log(NaN) = NaN */

  x  = _mm512_and_ps(x, _mm512_castsi512_ps(_mm512_set1_epi32(~0x7f800000))); /* x now the stored 23 bits of the 24-bit significand        */
  x  = _mm512_or_ps (x, v0p5);                                          /* sets hidden bit b[0]                                      */

  ei = _mm512_sub_epi32(ei, _mm512_set1_epi32(126));                       /* -127 (ei now signed base-2 exponent); then +1             */
  e  = _mm512_cvtepi32_ps(ei);

  /* now, calculate the log */
  mask = _mm512_cmp_ps_mask(x, _mm512_set1_ps(0.707106781186547524f), 1);
  // 1 = code for less-than comparison
  /* avoid conditional branches.           */
  tmp  = _mm512_mask_blend_ps(mask, _mm512_setzero_ps(), x);
  /* tmp contains x values < 0.707, else 0 */
  x    = _mm512_sub_ps(x, onev);
  e    = _mm512_sub_ps(e, _mm512_mask_blend_ps(mask, _mm512_setzero_ps(), onev));
  x    = _mm512_add_ps(x, tmp);
  z    = _mm512_mul_ps(x,x);

  y =               _mm512_set1_ps(cephes_p[0]);    y = _mm512_mul_ps(y, x); 
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[1]));   y = _mm512_mul_ps(y, x);    
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[2]));   y = _mm512_mul_ps(y, x);   
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[3]));   y = _mm512_mul_ps(y, x);   
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[4]));   y = _mm512_mul_ps(y, x);    
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[5]));   y = _mm512_mul_ps(y, x);   
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[6]));   y = _mm512_mul_ps(y, x); 
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[7]));   y = _mm512_mul_ps(y, x);  
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[8]));   y = _mm512_mul_ps(y, x);
  y = _mm512_mul_ps(y, z);

  tmp = _mm512_mul_ps(e, _mm512_set1_ps(-2.12194440e-4f));
  y   = _mm512_add_ps(y, tmp);

  tmp = _mm512_mul_ps(z, v0p5);
  y   = _mm512_sub_ps(y, tmp);

  tmp = _mm512_mul_ps(e, _mm512_set1_ps(0.693359375f));
  x = _mm512_add_ps(x, y);
  x = _mm512_add_ps(x, tmp);

  /* IEEE754 cleanup: */
  x = esl_avx512_select_ps(x, origx,                     inf_mask);  /* log(inf)=inf; log(NaN)      = NaN  */
  x = esl_avx512_select_ps(x, _mm512_set1_ps(eslNaN), invalid_mask);                                 /* log(x<0, including -0,-inf) = NaN  */
  x = esl_avx512_select_ps(x, _mm512_set1_ps(-eslINFINITY), zero_mask); /* x zero or subnormal         = -inf */
  return x;
}


/* Function:  esl_avx512_expf()
 * Synopsis:  <r[z] = exp x[z]>
 * Incept:    NPC, April 30 2024
 *
 * Purpose:   Given a vector <x> containing sixteen floats, returns a
 *            vector <r> in which each element <r[z] = expf(x[z])>.
 *            
 *            Valid for all IEEE754 floats $x_z$.
 *            
 * Xref:      J2/71
 *            J10/62: bugfix, minlogf/maxlogf range was too wide; 
 *                    (k+127) must be >=0 and <=255, so (k+127)<<23
 *                    is a valid IEEE754 float, without touching 
 *                    the sign bit. Pommier had this right in the
 *                    first place, and I didn't understand.
 * 
 * Note:      Derived from an SSE1 implementation by Julian
 *            Pommier. Converted to SSE2. Converted to AVX512 by
 *            Nick Carter.
 *            
 *            Note on maxlogf/minlogf, which are close to but not
 *            exactly 127.5/log2 [J10/63]. We need -127<=k<=128, so
 *            k+127 is 0..255, a valid IEEE754 8-bit exponent
 *            (0..255), so the bit pattern (k+127)<<23 is IEEE754
 *            single-precision for 2^k.  If k=-127, we get IEEE754 0.
 *            If k=128, we get IEEE754 +inf.  If k<-127, k+127 is
 *            negative and we get screwed up.  If k>128, k+127
 *            overflows the 8-bit exponent and sets the sign bit.  So
 *            for x' (base 2) < -127.5 we must definitely return e^x ~
 *            0; for x' < 126.5 we're going to calculate 0 anyway
 *            (because k=floor(-126.5-epsilon+0.5) = -127).  So any
 *            minlogf between -126.5 log2 ... -127.5 log2 will suffice
 *            as the cutoff. Ditto for 126.5 log2 .. 127.5log2.
 *            That's 87.68312 .. 88.3762655.  I think Pommier's
 *            thinking is, you don't want to get to close to the
 *            edges, lest fp roundoff error screw you (he may have
 *            consider 1 ulp carefully, I can't tell), but otherwise
 *            you may as well put your bounds close to the outer edge;
 *            so 
 *              maxlogf =  127.5 log(2) - epsilon 
 *              minlogf = -127.5 log(2) + epsilon 
 *            for an epsilon that happen to be ~ 3e-6.
 */
__m512 
esl_avx512_expf(__m512 x) 
{
  static float cephes_p[6] = { 1.9875691500E-4f, 1.3981999507E-3f, 8.3334519073E-3f, 
			       4.1665795894E-2f, 1.6666665459E-1f, 5.0000001201E-1f };
  static float cephes_c[2] = { 0.693359375f,    -2.12194440e-4f };
  static float maxlogf     =  88.3762626647949f;  /* 127.5 log(2) - epsilon. above this, 0.5+x/log2 gives k>128 and breaks 2^k "float" construction, because (k+127)<<23 must be a valid IEEE754 exponent 0..255 */
  static float minlogf     = -88.3762626647949f;  /*-127.5 log(2) + epsilon. below this, 0.5+x/log2 gives k<-127 and breaks 2^k, see above */
  __m512i k;
  __m512  tmp, fx, z, y;
  __mmask16 minmask, maxmask;
  
  /* handle out-of-range and special conditions */
  maxmask = _mm512_cmp_ps_mask(x, _mm512_set1_ps(maxlogf), 14);
  // 14 = greater than
  minmask = _mm512_cmp_ps_mask(x, _mm512_set1_ps(minlogf), 2);
  // 2 = less than or equal
  
  /* range reduction: exp(x) = 2^k e^f = exp(f + k log 2); k = floorf(0.5 + x / log2): */
  fx = _mm512_mul_ps(x,  _mm512_set1_ps(eslCONST_LOG2R));
  fx = _mm512_add_ps(fx, _mm512_set1_ps(0.5f));
  fx = _mm512_floor_ps(fx);
  k = _mm512_cvtps_epi32(fx);
  
  /* polynomial approx for e^f for f in range [-0.5, 0.5] */
  tmp = _mm512_mul_ps(fx, _mm512_set1_ps(cephes_c[0]));
  z   = _mm512_mul_ps(fx, _mm512_set1_ps(cephes_c[1]));
  x   = _mm512_sub_ps(x, tmp);
  x   = _mm512_sub_ps(x, z);
  z   = _mm512_mul_ps(x, x);
  
  y =               _mm512_set1_ps(cephes_p[0]);    y = _mm512_mul_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[1]));   y = _mm512_mul_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[2]));   y = _mm512_mul_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[3]));   y = _mm512_mul_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[4]));   y = _mm512_mul_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(cephes_p[5]));   y = _mm512_mul_ps(y, z);
  y = _mm512_add_ps(y, x);
  y = _mm512_add_ps(y, _mm512_set1_ps(1.0f));

  /* build 2^k by hand, by creating a IEEE754 float */
  k  = _mm512_add_epi32(k, _mm512_set1_epi32(127));
  k  = _mm512_slli_epi32(k, 23);
  fx = _mm512_castsi512_ps(k);
  
  /* put 2^k e^f together (fx = 2^k,  y = e^f) and we're done */
  y = _mm512_mul_ps(y, fx);	

  /* special/range cleanup */
  y = esl_avx512_select_ps(y, _mm512_set1_ps(eslINFINITY), maxmask); /* exp(x) = inf for x > log(2^128)  */
  y = esl_avx512_select_ps(y, _mm512_set1_ps(0.0f),        minmask); /* exp(x) = 0   for x < log(2^-149) */
  return y;
}

/*****************************************************************
 * 2. Debugging/development routines
 *****************************************************************/

void 
esl_avx512_dump_512i_hex8(__m512i v)
{
  uint64_t *val = (uint64_t*) &v;
  printf("%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", 
	 val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void
esl_avx512_dump_ps(FILE *fp, __m512 v)
{
  float *p = (float *)&v;
  fprintf(fp, "[%13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g, %13.8g]", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef eslAVX512_TESTDRIVE

#include "esl_random.h"
#include "esl_getopts.h"
#include <float.h>
/* utest_logf():  Test range/domain of logf */
static void utest_logf(ESL_GETOPTS *go)
{
  __m512 x;                            /* test input  */
  union { __m512 v; float x[16]; } r;   /* test output */

  /* Test IEEE754 specials:
   *    log(-inf) = NaN     log(x<0)  = NaN  log(-0)   = NaN
   *    log(0)    = -inf    log(inf)  = inf  log(NaN)  = NaN
   */
  x   = _mm512_set_ps(0.0, -0.0, -1.0, -eslINFINITY, 0.0, -0.0, -1.0, -eslINFINITY, 0.0, -0.0, -1.0, -eslINFINITY, 0.0, -0.0, -1.0, -eslINFINITY); /* set_ps() is in order 7 6 5 4 3 2 1 0 */
  r.v =  esl_avx512_logf(x);
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("logf");
    esl_avx512_dump_ps(stdout, x);    printf(" ==> ");
    esl_avx512_dump_ps(stdout, r.v);  printf("\n");
  }
  if (! isnan(r.x[0]))                 esl_fatal("logf(-inf) should be NaN");
  if (! isnan(r.x[1]))                 esl_fatal("logf(-1)   should be NaN");
  if (! isnan(r.x[2]))                 esl_fatal("logf(-0)   should be NaN");
  if (! (r.x[3] < 0 && isinf(r.x[3]))) esl_fatal("logf(0)    should be -inf");
  if (! isnan(r.x[4]))                 esl_fatal("logf(-inf) should be NaN");
  if (! isnan(r.x[5]))                 esl_fatal("logf(-1)   should be NaN");
  if (! isnan(r.x[6]))                 esl_fatal("logf(-0)   should be NaN");
  if (! (r.x[7] < 0 && isinf(r.x[7]))) esl_fatal("logf(0)    should be -inf");
  if (! isnan(r.x[8]))                 esl_fatal("logf(-inf) should be NaN");
  if (! isnan(r.x[9]))                 esl_fatal("logf(-1)   should be NaN");
  if (! isnan(r.x[10]))                 esl_fatal("logf(-0)   should be NaN");
  if (! (r.x[11] < 0 && isinf(r.x[11]))) esl_fatal("logf(0)    should be -inf");
  if (! isnan(r.x[12]))                 esl_fatal("logf(-inf) should be NaN");
  if (! isnan(r.x[13]))                 esl_fatal("logf(-1)   should be NaN");
  if (! isnan(r.x[14]))                 esl_fatal("logf(-0)   should be NaN");
  if (! (r.x[15] < 0 && isinf(r.x[15]))) esl_fatal("logf(0)    should be -inf");
  x   = _mm512_set_ps(FLT_MAX, FLT_MIN, eslNaN, eslINFINITY, FLT_MAX, FLT_MIN, eslNaN, eslINFINITY, FLT_MAX, FLT_MIN, eslNaN, eslINFINITY, FLT_MAX, FLT_MIN, eslNaN, eslINFINITY);
  r.v = esl_avx512_logf(x);
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("logf");
    esl_avx512_dump_ps(stdout, x);    printf(" ==> ");
    esl_avx512_dump_ps(stdout, r.v);  printf("\n");
  }
  if (! isinf(r.x[0]))  esl_fatal("logf(inf)  should be inf");
  if (! isnan(r.x[1]))  esl_fatal("logf(NaN)  should be NaN");
  if (! isinf(r.x[4]))  esl_fatal("logf(inf)  should be inf");
  if (! isnan(r.x[5]))  esl_fatal("logf(NaN)  should be NaN");
  if (! isinf(r.x[8]))  esl_fatal("logf(inf)  should be inf");
  if (! isnan(r.x[9]))  esl_fatal("logf(NaN)  should be NaN");
  if (! isinf(r.x[12]))  esl_fatal("logf(inf)  should be inf");
  if (! isnan(r.x[13]))  esl_fatal("logf(NaN)  should be NaN");

}


/* utest_expf():  Test range/domain of expf */
static void
utest_expf(ESL_GETOPTS *go)
{
  __m512 x;			       /* test input  */
  union { __m512 v; float x[16]; } r;   /* test output */
  
  /* exp(-inf) = 0    exp(-0)  = 1   exp(0) = 1  exp(inf) = inf   exp(NaN)  = NaN */
  x = _mm512_set_ps(eslINFINITY, 0.0, -0.0, -eslINFINITY, eslINFINITY, 0.0, -0.0, -eslINFINITY, eslINFINITY, 0.0, -0.0, -eslINFINITY, eslINFINITY, 0.0, -0\
.0, -eslINFINITY); /* set_ps() is in order 7 6 5 4 3 2 1 0 */
  r.v =  esl_avx512_expf(x); 
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("expf");
    esl_avx512_dump_ps(stdout, x);    printf(" ==> ");
    esl_avx512_dump_ps(stdout, r.v);  printf("\n");
  }
  if (r.x[0] != 0.0f)   esl_fatal("expf(-inf) should be 0");
  if (r.x[1] != 1.0f)   esl_fatal("expf(-0)   should be 1");
  if (r.x[2] != 1.0f)   esl_fatal("expf(0)    should be 1");
  if (! isinf(r.x[3]))  esl_fatal("expf(inf)  should be inf");
  if (r.x[4] != 0.0f)   esl_fatal("expf(-inf) should be 0");
  if (r.x[5] != 1.0f)   esl_fatal("expf(-0)   should be 1");
  if (r.x[6] != 1.0f)   esl_fatal("expf(0)    should be 1");
  if (! isinf(r.x[7]))  esl_fatal("expf(inf)  should be inf");
  if (r.x[8] != 0.0f)   esl_fatal("expf(-inf) should be 0");
  if (r.x[9] != 1.0f)   esl_fatal("expf(-0)   should be 1");
  if (r.x[10] != 1.0f)   esl_fatal("expf(0)    should be 1");
  if (! isinf(r.x[11]))  esl_fatal("expf(inf)  should be inf");
  if (r.x[12] != 0.0f)   esl_fatal("expf(-inf) should be 0");
  if (r.x[13] != 1.0f)   esl_fatal("expf(-0)   should be 1");
  if (r.x[14] != 1.0f)   esl_fatal("expf(0)    should be 1");
  if (! isinf(r.x[15]))  esl_fatal("expf(inf)  should be inf");
  /* exp(NaN) = NaN    exp(large)  = inf   exp(-large) = 0  exp(1) = exp(1) */
  x = _mm512_set_ps(1.0f, -666.0f, 666.0f, eslNaN, 1.0f, -666.0f, 666.0f, eslNaN, 1.0f, -666.0f, 666.0f, eslNaN, 1.0f, -666.0f, 666.0f, eslNaN); /* set_ps() is in order 7 6 5 4 3 2 1 0 */
  r.v =  esl_avx512_expf(x); 
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("expf");
    esl_avx512_dump_ps(stdout, x);    printf(" ==> ");
    esl_avx512_dump_ps(stdout, r.v);  printf("\n");
  }
  if (! isnan(r.x[0]))  esl_fatal("expf(NaN)      should be NaN");
  if (! isinf(r.x[1]))  esl_fatal("expf(large x)  should be inf");
  if (r.x[2] != 0.0f)   esl_fatal("expf(-large x) should be 0");
    if (! isnan(r.x[4]))  esl_fatal("expf(NaN)      should be NaN");
  if (! isinf(r.x[5]))  esl_fatal("expf(large x)  should be inf");
  if (r.x[6] != 0.0f)   esl_fatal("expf(-large x) should be 0");
  if (! isnan(r.x[8]))  esl_fatal("expf(NaN)      should be NaN");
  if (! isinf(r.x[9]))  esl_fatal("expf(large x)  should be inf");
  if (r.x[10] != 0.0f)   esl_fatal("expf(-large x) should be 0");
    if (! isnan(r.x[12]))  esl_fatal("expf(NaN)      should be NaN");
  if (! isinf(r.x[13]))  esl_fatal("expf(large x)  should be inf");
  if (r.x[14] != 0.0f)   esl_fatal("expf(-large x) should be 0");
  /* Make sure we are correct around the problematic ~minlogf boundary:
   *  (1) e^x for x < -127.5 log2 + epsilon is 0, because that's our minlogf barrier.
   *  (2) e^x for  -127.5 log2 < x < -126.5 log2 is 0 too, but is actually calculated
   *  (3) e^x for  -126.5 log2 < x should be finite (and close to FLT_MIN)
   *
   *  minlogf = -127.5 log(2) + epsilon = -88.3762626647949;
   *        and -126.5 log(2)           = -87.68311834
   *  so for
   *     (1): expf(-88.3763)  => 0
   *     (2): expf(-88.3762)  => 0
   *     (3): expf(-87.6832)   => 0
   *     (4): expf(-87.6831)   => <FLT_MIN (subnormal) : ~8.31e-39 (may become 0 in flush-to-zero mode for subnormals)
   */
  x   = _mm512_set_ps(-88.3763, -88.3762, -87.6832, -87.6831, -88.3763, -88.3762, -87.6832, -87.6831, -88.3763, -88.3762, -87.6832, -87.6831, -88.3763, -88.3762, -87.6832, -87.6831);
  r.v = esl_avx512_expf(x); 
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("expf");
    esl_avx512_dump_ps(stdout, x);    printf(" ==> ");
    esl_avx512_dump_ps(stdout, r.v);  printf("\n");
  }
  if ( r.x[0] >= FLT_MIN) esl_fatal("expf( -126.5 log2 + eps) should be around FLT_MIN");
  if ( r.x[1] != 0.0f)    esl_fatal("expf( -126.5 log2 - eps) should be 0.0 (by calculation)");
  if ( r.x[2] != 0.0f)    esl_fatal("expf( -127.5 log2 + eps) should be 0.0 (by calculation)");
  if ( r.x[3] != 0.0f)    esl_fatal("expf( -127.5 log2 - eps) should be 0.0 (by min bound): %g", r.x[0]);
   if ( r.x[4] >= FLT_MIN) esl_fatal("expf( -126.5 log2 + eps) should be around FLT_MIN");
  if ( r.x[5] != 0.0f)    esl_fatal("expf( -126.5 log2 - eps) should be 0.0 (by calculation)");
  if ( r.x[6] != 0.0f)    esl_fatal("expf( -127.5 log2 + eps) should be 0.0 (by calculation)");
  if ( r.x[7] != 0.0f)    esl_fatal("expf( -127.5 log2 - eps) should be 0.0 (by min bound): %g", r.x[7]);
  if ( r.x[8] >= FLT_MIN) esl_fatal("expf( -126.5 log2 + eps) should be around\
 FLT_MIN");
  if ( r.x[9] != 0.0f)    esl_fatal("expf( -126.5 log2 - eps) should be 0.0 (b\
y calculation)");
  if ( r.x[10] != 0.0f)    esl_fatal("expf( -127.5 log2 + eps) should be 0.0 (b\
y calculation)");
  if ( r.x[11] != 0.0f)    esl_fatal("expf( -127.5 log2 - eps) should be 0.0 (b\
y min bound): %g", r.x[0]);
   if ( r.x[12] >= FLT_MIN) esl_fatal("expf( -126.5 log2 + eps) should be aroun\
d FLT_MIN");
  if ( r.x[13] != 0.0f)    esl_fatal("expf( -126.5 log2 - eps) should be 0.0 (b\
y calculation)");
  if ( r.x[14] != 0.0f)    esl_fatal("expf( -127.5 log2 + eps) should be 0.0 (b\
y calculation)");
  if ( r.x[15] != 0.0f)    esl_fatal("expf( -127.5 log2 - eps) should be 0.0 (b\
y min bound): %g", r.x[15]);
}


/* utest_odds():  test accuracy of logf, expf on odds ratios,
 * our main intended use.
 */
static void
utest_odds(ESL_GETOPTS *go, ESL_RANDOMNESS *r)
{
  int    N            = esl_opt_GetInteger(go, "-N");
  int    verbose      = esl_opt_GetBoolean(go, "-v");
  int    very_verbose = esl_opt_GetBoolean(go, "--vv");
  int    i;
  float  p1, p2, odds;
  union { __m512 v; float x[16]; } r1;   
  union { __m512 v; float x[16]; } r2;   
  float  scalar_r1, scalar_r2;
  double  err1, maxerr1 = 0.0, avgerr1 = 0.0; /* errors on logf() */
  double  err2, maxerr2 = 0.0, avgerr2 = 0.0; /* errors on expf() */

  for (i = 0; i < N; i++)
    {
      p1    = esl_rnd_UniformPositive(r);
      p2    = esl_rnd_UniformPositive(r);
      odds  = p1 / p2;

      if (odds == 0.0) esl_fatal("whoa, odds ratio can't be 0!\n");

      r1.v      = esl_avx512_logf(_mm512_set1_ps(odds));  /* r1.x[z] = log(p1/p2) */
      scalar_r1 = log(odds);

      err1       = (r1.x[0] == 0. && scalar_r1 == 0.) ? 0.0 : 2 * fabs(r1.x[0] - scalar_r1) / fabs(r1.x[0] + scalar_r1);
      if (err1 > maxerr1) maxerr1 = err1;
      avgerr1   += err1 / (float) N;
      if (isnan(avgerr1)) esl_fatal("whoa, what?\n");

      r2.v      = esl_avx512_expf(r1.v);        /* and back to odds */
      scalar_r2 = exp(r1.x[0]);

      err2       = (r2.x[0] == 0. && scalar_r2 == 0.) ? 0.0 : 2 * fabs(r2.x[0] - scalar_r2) / fabs(r2.x[0] + scalar_r2);
      if (err2 > maxerr2) maxerr2 = err2;
      avgerr2   += err2 / (float) N;

      if (very_verbose) 
	printf("%13.7g  %13.7g  %13.7g  %13.7g  %13.7g  %13.7g  %13.7g\n", odds, scalar_r1, r1.x[0], scalar_r2, r2.x[0], err1, err2);
    }

  if (verbose) {
    printf("Average [max] logf() relative error in %d odds trials:  %13.8g  [%13.8g]\n", N, avgerr1, maxerr1);
    printf("Average [max] expf() relative error in %d odds trials:  %13.8g  [%13.8g]\n", N, avgerr2, maxerr2);
    printf("(random seed : %" PRIu32 ")\n", esl_randomness_GetSeed(r));
  }

  if (avgerr1 > 1e-8) esl_fatal("average error on logf() is intolerable\n");
  if (maxerr1 > 1e-6) esl_fatal("maximum error on logf() is intolerable\n");
  if (avgerr2 > 1e-8) esl_fatal("average error on expf() is intolerable\n");
  if (maxerr2 > 1e-6) esl_fatal("maximum error on expf() is intolerable\n");
}

static void
utest_hmax_epu8(ESL_RANDOMNESS *rng)
{
  union { __m512i v; uint8_t x[64]; } u;
  uint8_t r1, r2;
  int     i,z;

  for (i = 0; i < 100; i++)
    {
      r1 = 0;
      for (z = 0; z < 64; z++) 
        {
          u.x[z] = (uint8_t) (esl_rnd_Roll(rng, 256));  // 0..255
          if (u.x[z] > r1) r1 = u.x[z];
        }
      r2 = esl_avx512_hmax_epu8(u.v);
      if (r1 != r2) esl_fatal("hmax_epu8 utest failed");
    }
}

static void
utest_hmax_epi8(ESL_RANDOMNESS *rng)
{
  union { __m512i v; int8_t x[64]; } u;
  int8_t r1, r2;
  int    i,z;

  for (i = 0; i < 100; i++) 
    {
      r1 = 0;
      for (z = 0; z < 64; z++) 
        {
          u.x[z] = (int8_t) (esl_rnd_Roll(rng, 256) - 128);  // -128..127
          if (u.x[z] > r1) r1 = u.x[z];
        }
      r2 = esl_avx512_hmax_epi8(u.v);
      if (r1 != r2) esl_fatal("hmax_epi8 utest failed");
    }
}


static void
utest_hmax_epi16(ESL_RANDOMNESS *rng)
{
  union { __m512i v; int16_t x[32]; } u;
  int16_t r1, r2;
  int     i,z;

  for (i = 0; i < 100; i++) 
    {
      r1 = -32768;
      for (z = 0; z < 32; z++) 
        {
          u.x[z] = (int16_t) (esl_rnd_Roll(rng, 65536) - 32768);  // -32768..32767
          if (u.x[z] > r1) r1 = u.x[z];
        }
      r2 = esl_avx512_hmax_epi16(u.v);
      if (r1 != r2) esl_fatal("hmax_epi16 utest failed %d %d", r1, r2);
    }
}

#endif /*eslAVX512_TESTDRIVE*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/

#ifdef eslAVX512_TESTDRIVE
#include <esl_config.h>

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_cpu.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_avx512.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0",  NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
    { "-N",        eslARG_INT,  "10000",  NULL, NULL,  NULL,  NULL, NULL, "numbe\
r of random test points",                     0 },
  { "-v",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "be ve\
rbose: show test report",                     0 },
  { "--vv",      eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "be ve\
ry verbose: show individual test samples",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for avx512 module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  if (esl_cpu_has_avx512())
    {
      utest_hmax_epu8(rng);
      utest_hmax_epi8(rng);
      utest_hmax_epi16(rng);
      utest_expf(go);
      utest_logf(go);
      utest_odds(go, rng);
    }
  else
    {
      fprintf(stderr, "processor does not support our AVX-512 code; skipping tests.\n");
      fprintf(stderr, "  (we need KNL's F,CD,ER,PF subsets, plus DQ,BW)\n");
    }

  fprintf(stderr, "#  status = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*eslAVX512_TESTDRIVE*/



#else  // ! eslENABLE_AVX512
#include <stdio.h>
void esl_avx512_silence_hack(void) { return; }
#if defined eslAVX512_TESTDRIVE || eslAVX512_EXAMPLE || eslAVX512_BENCHMARK
int main(void) { fprintf(stderr, "# AVX512 support not compiled.\n"); return 0; }
#endif
#endif // eslENABLE_AVX512
