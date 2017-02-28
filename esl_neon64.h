/* Vectorized routines for ARM, using NEON technology.
 *
 * This header file, unusually, provides many complete function
 * implementations, so they can be inlined by the compiler.
 * 
 * Contents:
 *    1. Data structures for ARM/Intel intrinsics compatibility
 *    2. Function declarations for esl_neon64
 *    3. Inlined utilities for float vectors (4 floats in esl_neon_128f_t)
 *    4. Inlined utilities for epu8 vectors (16 uchars in esl_neon_128i_t)
 */

#ifdef  HAVE_NEON64
#ifndef eslNEON64_INCLUDED
#define eslNEON64_INCLUDED

#include "easel.h"
#include <stdio.h>
#include <arm_neon.h>
#include "esl_neon.h"

/*****************************************************************
 * 1. Function declarations (from esl_neon64.c)
 *****************************************************************/

extern esl_neon_128f_t  esl_neon64_logf(esl_neon_128f_t x);
extern esl_neon_128f_t  esl_neon64_expf(esl_neon_128f_t x);
extern void             esl_neon64_dump_float(FILE *fp, esl_neon_128f_t v);


/*****************************************************************
 * 2. Inline utilities for ps vectors (4 floats in esl_neon_128f_t)
 *****************************************************************/

/* Function:  esl_neon64_select_float()
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
static inline esl_neon_128f_t
esl_neon64_select_float(esl_neon_128f_t a, esl_neon_128f_t b, esl_neon_128f_t mask)
{
  esl_neon_128i_t aview, bview, maskview, masknot;
  esl_neon_128f_t ret;

  maskview.s64x2 = vreinterpretq_s64_f32(mask.f32x4);
  bview.s64x2    = vreinterpretq_s64_f32(b.f32x4);
  aview.s64x2    = vreinterpretq_s64_f32(a.f32x4);
  bview.s64x2    = vandq_s64(bview.s64x2, maskview.s64x2);
  masknot.s32x4  = vmvnq_s32(maskview.s32x4);
  aview.s64x2    = vandq_s64(masknot.s64x2, aview.s64x2);
  ret.f32x4      = vreinterpretq_f32_s64(vorrq_s64(aview.s64x2,bview.s64x2));  
  return ret; 
}


/* Function:  esl_neon64_any_gt_float()
 * Synopsis:  Returns TRUE if any a[z] > b[z]
 *
 * Purpose:   Returns TRUE if any a[z] > b[z] in two
 *            <ps> vectors of floats.
 *
 * Note:      Ported from esl_sse.c::esl_sse_any_gt_float().
 */
static inline int 
esl_neon64_any_gt_float(esl_neon_128f_t a, esl_neon_128f_t b)
{
  esl_neon_128i_t mask;
  int             l0, l1;
  int             maskbits;

  mask.u32x4 = vcgtq_f32(a.f32x4,b.f32x4);
  l0         = vgetq_lane_u64(mask.u64x2, 0);
  l1         = vgetq_lane_u64(mask.u64x2, 1);
  maskbits   = l0 | l1;
  return maskbits != 0;
}


/* Function:  esl_neon64_hsum_float()
 * Synopsis:  Takes the horizontal sum of elements in a vector.
 *
 * Purpose:   Add the four float elements in vector <a>; return
 *            that sum in <*ret_sum>.
 */
static inline void
esl_neon64_hsum_float(esl_neon_128f_t a, float *ret_sum)
{
  *ret_sum = vaddvq_f32(a.f32x4);
}


/* Function:  esl_neon64_rightshift_float()
 * Synopsis:  Shift vector elements to the right.
 *
 * Purpose:   Returns a vector containing
 *            <{ b[0] a[0] a[1] a[2] }>:
 *            i.e. shift the values in <a> to the
 *            right, and load the first value of 
 *            <b> into the first slot.
 */
static inline esl_neon_128f_t 
esl_neon64_rightshift_float(esl_neon_128f_t a, esl_neon_128f_t b)
{
  register esl_neon_128f_t v;

  v.f32x4 = vrev64q_f32(b.f32x4);            /* b[1] b[0] b[3] b[2] */
  v.f32x4 = vextq_f32(v.f32x4, v.f32x4, 2);  /* b[3] b[2] b[1] b[0] */
  v.f32x4 = vextq_f32(v.f32x4, a.f32x4, 3);  /* b[0] a[0] a[1] a[2] */
  return v; 
}

/* Function:  esl_neon64_leftshift_float()
 * Synopsis:  Shift vector elements to the left.
 *
 * Purpose:   Returns a vector containing
 *            <{ a[1] a[2] a[3] b[0]}>:
 *            i.e. shift the values in <a> to the
 *            left and load the first value of 
 *            <b> into the first slot.
 */
static inline esl_neon_128f_t
esl_neon64_leftshift_float(esl_neon_128f_t a, esl_neon_128f_t b)
{
  register esl_neon_128f_t v;
  v.f32x4 = vextq_f32(a.f32x4, b.f32x4, 1);/* now a[1] a[2] a[3] b[0] */
  return v;
}


/*****************************************************************
 * 4. Inlined utilities for epu8 vectors (16 uchars in __m128i)
 *****************************************************************/ 

/* Function:  esl_neon64_any_gt_s16()
 * Synopsis:  Returns TRUE if any a[z] > b[z].
 *
 * Purpose:   Return TRUE if any <a[z] > b[z]> for <z=0..15>
 *            in two <s16> vectors.
 */
static inline int 
esl_neon64_any_gt_s16(esl_neon_128i_t a, esl_neon_128i_t b)
{
  esl_neon_128i_t mask;
  int64_t         l0, l1;
  int64_t         maskbits;

  mask.u16x8 = vcgtq_s16(a.s16x8,b.s16x8);
  l0         = vgetq_lane_u64(mask.u64x2, 0);
  l1         = vgetq_lane_u64(mask.u64x2, 1);
  maskbits   = l0 | l1;
  return maskbits != 0;
}

/* Function:  esl_neon64_hmax_u8()
 * Synopsis:  Return the max of the 16 elements in u8 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            a <u8> vector.
 */
static inline uint8_t
esl_neon64_hmax_u8(esl_neon_128i_t a)
{
  return vmaxvq_u8(a.u8x16);
}

/* Function:  esl_neon64_hmax_s16()
 * Synopsis:  Return the max of the 8 elements in s16 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <s8> vector.
 */
static inline int16_t
esl_neon64_hmax_s16(esl_neon_128i_t a)
{
  return vmaxvq_s16(a.s16x8);
}

#endif /* peslINCLUDED */
#endif /* HAVE_NEON64 */
