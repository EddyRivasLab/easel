/* Vectorized routines for ARM, using NEON technology.
 *
 * This header file, unusually, provides many complete function
 * implementations; this is so that they can be inlined by the
 * compiler, for maximum efficiency.
 * 
 * Contents:
 *    1. Function declarations (from esl_neon.c)
 *    2. Inlined utilities for ps vectors (4 floats in __arm128f)
 *    3. Inlined utilities for epu8 vectors (16 uchars in __arm128i)
 */
#include "easel.h"

#include <stdio.h>
#include "arm_vector.h"

/*****************************************************************
 * 1. Function declarations (from esl_sse.c)
 *****************************************************************/
/*extern __armm128f  esl_neon_logf(__arm128f x);
extern __arm128f  esl_neon_expf(__arm128f x);
extern void    esl_neon_dump_ps(FILE *fp, __arm128i v);
*/

/*****************************************************************
 * 2. Inline utilities for ps vectors (4 floats in __arm128f)
 *****************************************************************/

/* Function:  esl_neon_select_ps()
 * Synopsis:  NEON equivalent of <vec_sel()>
 *
 * Purpose:   Vector select. Returns a vector <r[z] = a[z]> where <mask[z]>
 *            is all 0's; <r[z] = b[z]> where <mask[z]> is all 1's.
 *            
 *            Useful for avoiding conditional branches. For example,
 *            to implement \ccode{if (a > 0) a += a;}:
 *            
 *            \begin{cchunk}
 *              mask = _mm_cmpgt_ps(a, _mm_setzero_ps());
 *              twoa = _mm_add_ps(a, a);
 *              a    = esl_sse_select_ps(a, twoa, mask);
 *            \end{cchunk}
 *
 */

/* SO it looks like AArch32 does not support floating-point vector 
 * bitwise and in its Advanced SIMD instruction set, so we can 
 * either work around this or move to AArch64. Back burner for now. */
/*
static inline __arm128f
esl_neon_select_ps(__arm128f a, __arm128f b, __arm128f mask)
{
  b.s64x2 = vandq_s64(b.s64x2, mask.s64x2);
  a.s64x2 = vbicq_s64(mask.s64x2, a.s64x2);
  printf("After bicq: %x is lane 0\n", vgetq_lane_s32(a.s32x4, 0));
  __arm128i ret; 
  ret.s64x2 = vorrq_s64(a.s64x2,b.s64x2);  
  return ret; 
}
*/

/* Function:  esl_neon_any_gt_float()
 * Synopsis:  Returns TRUE if any a[z] > b[z]
 *
 * Purpose:   Returns TRUE if any a[z] > b[z] in two
 *            <ps> vectors of floats.
 *
 * Xref:      From Apple Altivec/SSE migration guide.
 */


static inline int 
esl_neon_any_gt_float(__arm128f a, __arm128f b)
{
  __arm128i mask;
  int l0, l1;
  mask.u32x4    = vcgtq_f32(a.f32x4,b.f32x4);
  l0 = vgetq_lane_u64(mask.u64x2, 0);
  l1 = vgetq_lane_u64(mask.u64x2, 1);
  int maskbits = l0 | l1;
  return maskbits != 0;
}


/* Function:  esl_neon_hmax_float()
 * Synopsis:  Find the maximum of elements in a vector.
 *
 * Purpose:   Find the maximum valued element in the four float elements
 *            in <a>, and return that maximum value in <*ret_max>.
 *            
 * Xref:      J3/90 for benchmarking of some alternative implementations.
 */

/* Appears not to be used in SSE implementation of HMMER; skipped for now */

/*
static inline void
esl_neon_hmax_ps(__arm128f a, float *ret_max)
{
  float l0, l1, l2, l3;
  l0 = vgetq_lane_f32(a.f32x4, 0);
  l1 = vgetq_lane_f32(a.f32x4, 0);
  l2 = vgetq_lane_f32(a.f32x4, 0);
  l3 = vgetq_lane_f32(a.f32x4, 0);
  l0 = (l0 > l1) ? l0 : l1;
  l0 = (l0 > l2) ? l0 : l2;
  l0 = (l0 > l3) ? l0 : l3;
  _mm_store_ss(ret_max, a);
}
*/

/* Function:  esl_sse_hmin_ps()
 * Synopsis:  Find the minimum of elements in a vector.
 *
 * Purpose:   Find the minimum valued element in the four float elements
 *            in <a> and return that minimum value in <*ret_min>.
 */
/*
static inline void
esl_sse_hmin_ps(__m128 a, float *ret_min)
{
  a = _mm_min_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
  a = _mm_min_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(ret_min, a);
}
*/
/* Function:  esl_sse_hsum_float()
 * Synopsis:  Takes the horizontal sum of elements in a vector.
 *
 * Purpose:   Add the four float elements in vector <a>; return
 *            that sum in <*ret_sum>.
 */


static inline void
esl_sse_hsum_float(__arm128f a, float *ret_sum)
{
  __arm256f_composite fvec;  
  a.f32x4 = vaddq_f32(a.f32x4, vrev64q_f32(a.f32x4));
  
  fvec.f32x4x2 = vtrnq_f32(a.f32x4, a.f32x4);
  a.f32x4 = vaddq_f32(fvec.f32x4x2.val[0], fvec.f32x4x2.val[1]);
  vst1q_lane_f32(ret_sum, a.f32x4, 0);
}


/* Function:  esl_sse_rightshift_ps()
 * Synopsis:  Shift vector elements to the right.
 *
 * Purpose:   Returns a vector containing
 *            <{ b[0] a[0] a[1] a[2] }>:
 *            i.e. shift the values in <a> to the
 *            right, and load the first value of 
 *            <b> into the first slot.
 */

static inline __arm128f 
esl_neon_rightshift_float(__arm128f a, __arm128f b)
{
  register __arm128f v;
  float floats[5];
  vst1q_f32(&floats[1], a.f32x4); /* Store a[0] a[1] a[2] a[3] */
  vst1q_lane_f32(floats, b.f32x4, 0); /* Store b[0] just below a[0] */
  v.f32x4 = vld1q_f32(floats);
  return v;

}

/* Function:  esl_neon_leftshift_float()
 * Synopsis:  Shift vector elements to the left.
 *
 * Purpose:   Returns a vector containing
 *            <{ a[1] a[2] a[3] b[0]}>:
 *            i.e. shift the values in <a> to the
 *            left and load the first value of 
 *            <b> into the first slot.
 */

static inline __arm128f
esl_sse_leftshift_float(__arm128f a, __arm128f b)
{
  register __arm128f v;
  __arm128i ia, ib, iv;
  ia.u32x4 = vreinterpretq_u32_f32(a.f32x4); /* reinterpret as int for select */
  ib.u32x4 = vreinterpretq_u32_f32(b.f32x4);
  iv.u32x4 = vreinterpretq_u32_f32(v.f32x4); 
  iv.u32x4 = vextq_u32(ia.u32x4, ib.u32x4, 1);/* now a[1] a[2] a[3] b[0] */
  v.f32x4 = vreinterpretq_f32_u32(iv.u32x4);
  return v;
}



/*****************************************************************
 * 3. Inlined utilities for epu8 vectors (16 uchars in __m128i)
 *****************************************************************/ 

/* Function:  esl_sse_any_gt_epu8()
 * Synopsis:  Returns TRUE if any a[z] > b[z].
 *
 * Purpose:   Return TRUE if any <a[z] > b[z]> for <z=0..15>
 *            in two <epu8> vectors of unsigned chars.
 *            
 *            We need this incantation because SSE provides
 *            no <cmpgt_epu8> instruction.
 *            
 *            For equality tests, note that <cmpeq_epi8> works fine
 *            for unsigned ints though there is no <cmpeq_epu8>
 *            instruction either).
 * 
 *            See vec_any_gt
 */

//static inline int 
//esl_sse_any_gt_epu8(__m128i a, __m128i b)
//{
//  __m128i mask    = _mm_cmpeq_epi8(_mm_max_epu8(a,b), b); /* anywhere a>b, mask[z] = 0x0; elsewhere 0xff */
//  int   maskbits  = _mm_movemask_epi8(_mm_xor_si128(mask,  _mm_cmpeq_epi8(mask, mask))); /* the xor incantation is a bitwise inversion */
//  return maskbits != 0;
//}
/*
static inline int 
esl_sse_any_gt_epi16(__m128i a, __m128i b)
{
  return (_mm_movemask_epi8(_mm_cmpgt_epi16(a,b)) != 0); 
}
*/

/* Function:  esl_sse_hmax_u8()
 * Synopsis:  Return the max of the 16 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <epu8> vector.
 */
static inline uint8_t
esl_sse_hmax_u8(__arm128i a)
{
  a.u8x16 = vmaxq_u8(a.u8x16, vrev64q_u8(a.u8x16));
  a.u8x16 = vmaxq_u8(a.u8x16, vreinterpretq_u8_u32(vrev64q_u32(a.u32x4)));
  a.u8x16 = vmaxq_u8(a.u8x16, vrev64q_u8(a.u8x16)); 
  a.u8x16 = vmaxq_u8(a.u8x16, vreinterpretq_u8_u32(vrev64q_u32(a.u32x4)));
  return vgetq_lane_u8(a.u8x16, 15);
}

/* Function:  esl_sse_hmax_s16()
 * Synopsis:  Return the max of the 8 elements in s16 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <s8> vector.
 */
static inline int16_t
esl_sse_hmax_s16(__arm128i a)
{
  a.s16x8 = vmaxq_s16(a.s16x8, vrev64q_s16(a.s16x8));
  a.s16x8 = vmaxq_s16(a.s16x8, vreinterpretq_s16_s32(vrev64q_s32(a.s32x4)));
  a.s16x8 = vmaxq_s16(a.s16x8, vrev64q_s16(a.s16x8));
  return vgetq_lane_s16(a.s16x8, 7);
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
