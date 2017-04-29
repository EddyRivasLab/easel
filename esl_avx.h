/* Vectorized utility routines for Intel AVX instructions and 
 * compatible processors.
 *
 * This header file, unusually, provides many complete function
 * implementations; this is so that they can be inlined by the
 * compiler, for maximum efficiency.
 * 
 * Contents:
 *    1. Inlined horizontal functions for 8 and 16-bit quantities
 *       in 256-bit vectors (__m256i)
 */
#ifndef eslAVX_INCLUDED
#define eslAVX_INCLUDED
#include "esl_config.h"
#ifdef  eslENABLE_AVX

#include "easel.h"

#include <stdio.h>
#include <x86intrin.h>		


/* Function:  esl_avx_hmax_epu8()
 * Synopsis:  Return the unsigned max of the 32 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 32 elements in
 *            an <epu8> vector.
 */
static inline uint8_t
esl_avx_hmax_epu8(__m256i a)
{
      __m256i temp1_AVX = _mm256_permute2x128_si256(a, a, 0x01);
      // Swap the 128-bit halves from a into temp1

      __m256i temp2_AVX = _mm256_max_epu8(temp1_AVX, a); // each 8-bit field in temp2_AVX now has the max of the
      //corresponding fields in the high and low halves of a

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of a
      temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 8-bit fields from the 64-bit quarters of a

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
      temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 8 bit fields from the 32-bit eighths of a

      temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp2_AVX
      temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of a

      uint8_t temp_stash = _mm256_extract_epi8(temp2_AVX, 1);
      temp1_AVX = _mm256_insert_epi8(temp2_AVX, temp_stash, 0);  // low byte of temp1_AVX now has byte 2 of temp2_AVX
      temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of Dmaxv_AVX
      return(_mm256_extract_epi8(temp2_AVX, 0));  // get low byte of temp2_AVX
}

/* Function:  esl_avx_hmax_epi16()
 * Synopsis:  Return the signed max of the 16 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 32 elements in
 *            an <epu8> vector.
 */
static inline uint16_t
esl_avx_hmax_epi16(__m256i a)
{
      __m256i temp1_AVX = _mm256_permute2x128_si256(a, a, 0x01);
      // Swap the 128-bit halves from a into temp1

      __m256i temp2_AVX = _mm256_max_epi16(temp1_AVX, a); // each 8-bit field in temp2_AVX now has the max of the
      //corresponding fields in the high and low halves of a

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of a
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 8-bit fields from the 64-bit quarters of a

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 8 bit fields from the 32-bit eighths of a

      temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of a

      return(_mm256_extract_epi16(temp2_AVX, 0));  // get low 16 bits of temp2_AVX
}


/* naming conventions:  The left/right in the names of these functions refers to the direction of the SSE shift instruction
that they emulate, because that's what the first filters to be ported to AVX used.  The esl_sse_(left/right)shift functions 
are vector-logical, meaning that, on x86, they mirror the function of the shift instruction with the opposite name.  For
self-consistency, I'm sticking with names that match the direction of the instruction, even though this means that the SSE
and AVX filters call different functions.
*/

// shifts vector left by one byte
static inline __m256i esl_avx_leftshift_one(__m256i vector)
{
   register __m256i temp_mask_AVX = _mm256_permute2x128_si256(vector, vector, _MM_SHUFFLE(0,0,3,0) );
   return(_mm256_alignr_epi8(vector, temp_mask_AVX,15));
}

// shifts vector left by two bytes
static inline __m256i esl_avx_leftshift_two(__m256i vector)
{
   register __m256i temp_mask_AVX = _mm256_permute2x128_si256(vector, vector, _MM_SHUFFLE(0,0,3,0) );
   return(_mm256_alignr_epi8(vector, temp_mask_AVX,14));
}
// shifts vector left by four bytes (one float)
static inline __m256 esl_avx_leftshift_ps(__m256 vector)
{
   register __m256i temp_mask_AVX = _mm256_permute2x128_si256((__m256i) vector, (__m256i) vector, _MM_SHUFFLE(0,0,3,0) );
   return((__m256) _mm256_alignr_epi8((__m256i) vector, temp_mask_AVX,12));
}

// shifts vector right by four bytes (one float)
static inline __m256 esl_avx_rightshift_ps(__m256 vector)
{
   register __m256i temp1 = _mm256_permute2x128_si256((__m256i) vector, (__m256i) vector, 0x81);  //result has vector[255:128] in low 128 bits, 0 in high 128
   return((__m256) _mm256_alignr_epi8(temp1, (__m256i) vector,4));
}

/* Function:  esl_avx_hsum_ps()
 * Synopsis:  Takes the horizontal sum of elements in a vector.
 *
 * Purpose:   Add the four float elements in vector <a>; return
 *            that sum in <*ret_sum>.
 */
static inline void
esl_avx_hsum_ps(__m256 a, float *ret_sum)
{
 __m256 temp1_AVX = (__m256) _mm256_permute2x128_si256((__m256i) a, (__m256i) a, 0x01);
      // Swap the 128-bit halves from a into temp1

 __m256 temp2_AVX = _mm256_add_ps(a, temp1_AVX);
 // low 128 bits of temp2_AVX have the sum of the corresponding floats from the high, low
 // 128 bits of a

   temp1_AVX = (__m256) _mm256_shuffle_epi32((__m256i) temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of a
   temp2_AVX = _mm256_add_ps(temp1_AVX, temp2_AVX);  // low 64 bits of temp2_AVX now have the sums of the
   // corresponding floats from the quarters of a

   temp1_AVX = (__m256) _mm256_shuffle_epi32((__m256i) temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
   temp2_AVX = _mm256_add_ps(temp1_AVX, temp2_AVX);  // low 32 bits of temp2_AVX now have the sum of the floats in a

   int *retint_ptr = (int *) ret_sum;  // This is a horrible hack because there isn't an intrinsic to extract a float from
   // an __m256.  Do this to avoid casting an int back to a float and screwing it up
   *retint_ptr = _mm256_extract_epi32((__m256i) temp2_AVX, 0);
}


extern void esl_avx_dump_256i_hex4(__m256i v);

#endif /*eslAVX_INCLUDED*/
#endif // eslENABLE_AVX
