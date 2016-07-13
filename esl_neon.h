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

#ifdef HAVE_NEON
#ifndef eslNEON_INCLUDED
#define eslNEON_INCLUDED

#include "easel.h"

#include <stdio.h>
#include "arm_vector.h"




/*****************************************************************
 * 1. Function declarations (from esl_sse.c)
 *****************************************************************/
extern __arm128f  esl_neon_logf(__arm128f x);
extern __arm128f  esl_neon_expf(__arm128f x);
extern void    esl_neon_dump_float(FILE *fp, __arm128f v);


/*****************************************************************
 * 2. Inline utilities for ps vectors (4 floats in __arm128f)
 *****************************************************************/

/* Function:  esl_neon_select_float()
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

static inline __arm128f
esl_neon_select_float(__arm128f a, __arm128f b, __arm128f mask)
{
  
  __arm128i aview, bview, maskview, masknot;
  __arm128f ret;
  maskview.s64x2 = vreinterpretq_s64_f32(mask.f32x4);
  bview.s64x2 = vreinterpretq_s64_f32(b.f32x4);
  aview.s64x2 = vreinterpretq_s64_f32(a.f32x4);
  bview.s64x2 = vandq_s64(bview.s64x2, maskview.s64x2);
  masknot.s32x4 = vmvnq_s32(maskview.s32x4);
  aview.s64x2 = vandq_s64(masknot.s64x2, aview.s64x2);
  ret.f32x4 = vreinterpretq_f32_s64(vorrq_s64(aview.s64x2,bview.s64x2));  
  return ret; 
}


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


/* Function:  esl_neon_hsum_float()
 * Synopsis:  Takes the horizontal sum of elements in a vector.
 *
 * Purpose:   Add the four float elements in vector <a>; return
 *            that sum in <*ret_sum>.
 */


static inline void
esl_neon_hsum_float(__arm128f a, float *ret_sum)
{
  __arm128f fvec;  
  a.f32x4 = vaddq_f32(a.f32x4, vrev64q_f32(a.f32x4));
  fvec.f32x4 = vextq_f32(a.f32x4, a.f32x4, 2);
  a.f32x4 = vaddq_f32(a.f32x4, fvec.f32x4);
  vst1q_lane_f32(ret_sum, a.f32x4, 0);
}


/* Function:  esl_neon_rightshift_ps()
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
  union { __arm128f v; float f[4]; }tmp;
  tmp.v = a;
  tmp.v = b;
  register __arm128f v;
  v.f32x4 = vrev64q_f32(b.f32x4); /* b[1] b[0] b[3] b[2] */
  v.f32x4 = vextq_f32(v.f32x4, v.f32x4, 2);  /* b[3] b[2] b[1] b[0] */
  v.f32x4 = vextq_f32(v.f32x4, a.f32x4, 3); /* b[0] a[0] a[1] a[2] */
  tmp.v = v;
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
esl_neon_leftshift_float(__arm128f a, __arm128f b)
{
  register __arm128f v;
  v.f32x4 = vextq_f32(a.f32x4, b.f32x4, 1);/* now a[1] a[2] a[3] b[0] */
  return v;
}



/*****************************************************************
 * 3. Inlined utilities for epu8 vectors (16 uchars in __m128i)
 *****************************************************************/ 

/* Function:  esl_neon_any_gt_s16()
 * Synopsis:  Returns TRUE if any a[z] > b[z].
 *
 * Purpose:   Return TRUE if any <a[z] > b[z]> for <z=0..15>
 *            in two <s16> vectors.
 *            
 */

static inline int 
esl_neon_any_gt_s16(__arm128i a, __arm128i b)
{
   __arm128i mask;
  int64_t l0, l1;
  mask.u16x8 = vcgtq_s16(a.s16x8,b.s16x8);
  l0 = vgetq_lane_u64(mask.u64x2, 0);
  l1 = vgetq_lane_u64(mask.u64x2, 1);
  int64_t maskbits = l0 | l1;
  return maskbits != 0;
}

/* Function:  esl_neon_hmax_u8()
 * Synopsis:  Return the max of the 16 elements in u8 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            a <u8> vector.
 */
static inline uint8_t
esl_neon_hmax_u8(__arm128i a)
{
  register __arm128i tempv;
  tempv.u8x16 = vreinterpretq_u8_u32(vextq_u32(a.u32x4, a.u32x4, 2));
  a.u8x16 = vmaxq_u8(a.u8x16, tempv.u8x16);
  tempv.u8x16 = vreinterpretq_u8_u32(vextq_u32(a.u32x4, a.u32x4, 1));
  a.u8x16 = vmaxq_u8(a.u8x16, tempv.u8x16);
  tempv.u8x16 = vreinterpretq_u8_u16(vrev64q_u16(a.u16x8));
  a.u8x16 = vmaxq_u8(a.u8x16, tempv.u8x16);
  tempv.u8x16 = vrev64q_u8(a.u8x16);
  a.u8x16 = vmaxq_u8(a.u8x16, tempv.u8x16);

  return vgetq_lane_u8(a.u8x16, 15);
}

/* Function:  esl_neon_hmax_s16()
 * Synopsis:  Return the max of the 8 elements in s16 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <s8> vector.
 */
static inline int16_t
esl_neon_hmax_s16(__arm128i a)
{
  a.s16x8 = vmaxq_s16(a.s16x8, vrev64q_s16(a.s16x8));
  a.s16x8 = vmaxq_s16(a.s16x8, vreinterpretq_s16_s32(vrev64q_s32(a.s32x4)));
  a.s16x8 = vmaxq_s16(a.s16x8, vreinterpretq_s16_s32(vextq_s32(a.s32x4, a.s32x4, 2)));
  return vgetq_lane_s16(a.s16x8, 7);
}

#endif /* eslNEON_INCLUDED */
#endif /* HAVE_NEON */
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
