
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

#ifndef eslAVX512_INCLUDED
#define eslAVX512_INCLUDED
#ifdef HAVE_AVX512 // make sure we don't include AVX-512 instructions on architectures that can't handle them
#include "easel.h"

#include <stdio.h>
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#include <immintrin.h>		/* AVX, AVX-512 */

//temporary helper function.  Should not be in release version
 static void print_512(__m512 var){
  float *val = (float*) &var;
    printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
           val[15], val[14], val[13], val[12], val[11], val[10], 
           val[9], val[8],val[7], val[6], val[5], val[4], val[3], val[2], 
           val[1], val[0]);
 }


 /* Function:  esl_avx_512_hmax_epu8()
 * Synopsis:  Return the unsigned max of the 64 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 64 elements in
 *            an <epu8> vector.
 */
static inline uint8_t
esl_avx_512_hmax_epu8(__m512i a)
{
    //Use AVX instructions for this because AVX-512 can't extract 8-bit quantities
    // Intel has stated that there will be no performance penalty for switching between AVX-512 and AVX
      __m256i temp3 = _mm512_extracti32x8_epi32(a, 0);
      __m256i temp4 = _mm512_extracti32x8_epi32(a, 1);
      temp3 = _mm256_max_epu8(temp3, temp4);
      temp4 = _mm256_permute2x128_si256(temp3, temp3, 0x01);
      // Swap the 128-bit halves from xEv_AVX into temp1

      temp3 = _mm256_max_epu8(temp4, temp3); // each 8-bit field in temp3 now has the max of the
      //corresponding fields in the high and low halves of Dmaxv_AVX

      temp4 = _mm256_shuffle_epi32(temp3, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 8-bit fields from the 64-bit quarters of Dmaxv_AVX

      temp4 = _mm256_shuffle_epi32(temp3, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 8 bit fields from the 32-bit eighths of Dmaxv_AVX

      temp4 = _mm256_shufflelo_epi16(temp3, 0xb1); // bottom 32-bits of temp4 now contain the swapped 16-bit halves
      // of the low 32 bits of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  //bottom 16 bits of temp3 now contain the max of the 16-bit fields of Dmaxv_AVX

      uint8_t temp_stash2 = _mm256_extract_epi8(temp3, 1);
      temp4 = _mm256_insert_epi8(temp3, temp_stash2, 0);  // low byte of temp4 now has byte 2 of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  //bottom 16 bits of temp3 now contain the max of the 16-bit fields of Dmaxv_AVX
      return(_mm256_extract_epi8(temp3, 0));  // get low byte of temp3
 
}

/* Function:  esl_avx_512_hmax_epi16()
 * Synopsis:  Return the signed max of the 16 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 32 elements in
 *            an <epu8> vector.
 */
static inline uint16_t
esl_avx_512_hmax_epi16(__m512i a)
{
 	__m512i temp1_AVX_512 = _mm512_shuffle_i64x2(a, a, 0x4e);
    temp1_AVX_512 = _mm512_max_epi16(a, temp1_AVX_512);  // get max of corresponding 16-bit quantities in high and 
      // low halves of a

    __m256i temp3_AVX = _mm512_extracti64x4_epi64(temp1_AVX_512, 0);  //shift to normal AVX for 16-bit operations
    __m256i temp4_AVX = _mm256_permute2x128_si256(temp3_AVX, temp3_AVX, 0x01);
      // Swap the 128-bit halves from temp3 into temp4

    temp3_AVX = _mm256_max_epi16(temp3_AVX, temp4_AVX); // each 16-bit field in temp3_AVX now has the max of the
      //corresponding fields in the high and low halves of a

    temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp3_AVX
    temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  // Each 64-bit quantity in temp4 now has the max of the corresponding
      // 16-bit fields from the 64-bit eighths of a

    temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp3_AVX
    temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 16 bit fields from the 32-bit sixteenths of a

    temp4_AVX = _mm256_shufflelo_epi16(temp3_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp3_AVX
    temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of xEv_AVX

    return(_mm256_extract_epi16(temp3_AVX, 0));  // return those low 16 bits
      
}

/* Function: esl_avx_512_hsum_ps()
 * Synopsis: sums the floating-point values in an __m512 vector
 *           returns the result in ret_sum
 * Purpose:  To compute the sum of the 32-bit float elements of a 512-bit vector 
 */
static inline void
esl_avx_512_hsum_ps(__m512 a, float *ret_sum){
//  printf("In esl_avx_512_hsum_ps.  Input vector is: ");
 // print_512(a);
 __m512 temp1_AVX_512 = _mm512_shuffle_f32x4(a, a, 0x4e);  //swap high and low halves of a
 __m512 temp2_AVX_512 = _mm512_add_ps(a, temp1_AVX_512); // sum corresponding floats in the high, low halves
 
 temp1_AVX_512 = _mm512_shuffle_f32x4(temp2_AVX_512, temp2_AVX_512, 0xb1);  //swap high and low quarters of each half of temp2
 temp2_AVX_512 = _mm512_add_ps(temp2_AVX_512, temp1_AVX_512); // sum corresponding floats in the high, low quarters
 
 temp1_AVX_512 = _mm512_shuffle_ps(temp2_AVX_512, temp2_AVX_512, 0x4e);  //swap high and low eigths of each quarter of a
 temp2_AVX_512 = _mm512_add_ps(temp2_AVX_512, temp1_AVX_512); // sum corresponding floats in the high, low eighths
 
 temp1_AVX_512 = _mm512_shuffle_ps(temp2_AVX_512, temp2_AVX_512, 0xb1);  //swap high and low sixteenths of each eighth 
 temp2_AVX_512 = _mm512_add_ps(temp2_AVX_512, temp1_AVX_512); // each element of temp2_AVX_512 now contains the sum of all the floats in a
 
 __m256 temp3_AVX = _mm512_extractf32x8_ps(temp2_AVX_512, 0); //Grab the low half of temp2_AVX_512 
 // because AVX-512 doesn't provide an operation to extract one float from a 512-bit vector
// printf("output sum vector is: ");
// print_512(temp2_AVX_512);
 int *retint_ptr = (int *) ret_sum;  // This is a horrible hack because there isn't an intrinsic to extract a float from
   // an __m256.  Do this to avoid casting an int back to a float and screwing it up
  *retint_ptr = _mm256_extract_epi32((__m256i) temp3_AVX, 0);
}

// shifts vector left by one byte.  Uses a similar technique to the AVX macro, but is complicated by the 
// lack of permute2x128 instruction in AVX-512

static inline __m512i 
esl_avx_512_leftshift_one(__m512i vector)
{
	__mmask16 zero_low_128 = 0xfff0;   // maskz_shuffle zeroes out 32-bit fields where the mask is 0, so this will
	// zero out the low 128-bits
    __m512i temp_mask_AVX_512 = _mm512_maskz_shuffle_i32x4(zero_low_128, vector, vector, 0x90);
       	
    //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
    return(_mm512_alignr_epi8(vector, temp_mask_AVX_512,15));
	  
}
// shifts vector left by two bytes
static inline __m512i esl_avx_512_leftshift_two(__m512i vector){
   // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
  __mmask16 zero_low_128 = 0xfff0;   // maskz_shuffle zeroes out 32-bit fields where the mask is 0, so this will
  // zero out the low 128-bits
    __m512i temp_mask_AVX_512 = _mm512_maskz_shuffle_i32x4(zero_low_128, vector, vector, 0x90);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by two bytes
     return(_mm512_alignr_epi8(vector, temp_mask_AVX_512,14));
}

// shifts vector left by four bytes (one float)
static inline __m512 esl_avx_512_leftshift_ps(__m512 vector){
   // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
  __mmask16 zero_low_128 = 0xfff0;   // maskz_shuffle zeroes out 32-bit fields where the mask is 0, so this will
  // zero out the low 128-bits
    __m512i temp_mask_AVX_512 = _mm512_maskz_shuffle_i32x4(zero_low_128, (__m512i) vector, (__m512i) vector, 0x90);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by two bytes
    __m512 retval = (__m512) _mm512_alignr_epi8((__m512i) vector, temp_mask_AVX_512,12);
return(retval);
  }

// shifts vector right by four bytes (one float)
static inline __m512 esl_avx_512_rightshift_ps(__m512 vector){
__mmask16 zero_high_128 = 0x0fff;   // maskz_shuffle zeroes out 32-bit fields where the mask is 0, so this will
  // zero out the low 128-bits
    __m512i temp_mask_AVX_512 = _mm512_maskz_shuffle_i32x4(zero_high_128, (__m512i) vector, (__m512i) vector, 0x39);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by two bytes
     return((__m512) _mm512_alignr_epi8(temp_mask_AVX_512, (__m512i) vector,4));
  }
#endif //HAVE_AVX512
#endif //eslAVX512_INCLUDED