/* Data structures for the ARM AArch32/AArch64 architectures NEON technology.
 * 
 * The data structures in this file exist for compatibility between Intel's
 * vector intrinsics (SSE/SSE2/SSE3/AVX) and ARM NEON intrinsics. Intel's 
 * vectorization code utilizes a single type for each view of its vector
 * registers; for example:
 * 
 * __m128 a = _mm_and_ps(...)
 * 
 * would be used for any combination of element sizes and lane numbers 
 * for some Intel vector register mapped to the C variable 'a'.
 *
 * By contrast, ARM requires the programmer to specify both the element 
 * size and the number of lanes when mapping a C variable onto a NEON 
 * register:
 *
 * uint32x4_t a = vdupq_n_s32(...)
 *
 * For compatibility reasons, and to simplify the code porting process and 
 * code maintainability, we define here a union type for each different view 
 * of the 128-bit registers.
 * 
 */ 

#ifndef ARM_VECTOR
#define ARM_VECTOR

#include <arm_neon.h>

/* Union type for vectorized integers. 
 *
 * Fields are named according to the following scheme in keeping with standard 
 * ARM NEON intrinsic naming/type conventions: 
 * 
 * <signed/unsigned><element size>x<number of lanes> 
 * 
 * For example:
 *
 * __arm128i vector.u64x2 
 * 
 * views the 128-bit register as 2 lanes of 64-bit integers. 
 * 
 */
typedef union arm128i
{
	int8x16_t s8x16;
	int16x8_t s16x8;
	int32x4_t s32x4;
	int64x2_t s64x2;
	int8x8x2_t s8x8x2;
	uint8x16_t u8x16;
	uint16x8_t u16x8;
	uint32x4_t u32x4;
	uint64x2_t u64x2;
	uint8x8x2_t u8x8x2;
} __arm128i;

typedef union arm64i
{
	int8x8_t s8x8;
	uint8x8_t u8x8;
	int64x1_t s64x1;
	uint64x1_t u64x1;
} __arm64i;
/* Union type for vectorized floating point values. Note: AArch32 does not 
 * allow double-precision floating-point vector operations; this was newly 
 * introduced in AArch64. */
typedef union arm64f
{
	float16x4_t f16x4;
	float32x2_t f32x2;
} __arm64f;

typedef union arm128f
{
	float32x4_t f32x4;
} __arm128f;
/* Union type for polynomial values. */
typedef union arm128p
{
	poly8x16_t p8x16;
	poly16x8_t p16x8;
} __arm128p;

/* Composite types */
typedef union arm128i_composite
{
	int8x8x2_t s8x8x2;
	int16x4x2_t s16x4x2;
	int32x2x2_t s32x2x2;
	uint8x8x2_t u8x8x2;
	uint16x4x2_t u16x4x2;
	uint32x2x2_t u32x2x2;
	uint64x1_t u64x1; /* useful for loading constants */
} __arm128i_composite;

typedef union arm256i_composite
{
	int8x16x2_t s8x16x2;
	int16x8x2_t s16x8x2;
	int32x4x2_t s32x4x2;
	uint8x16x2_t u8x16x2;
	uint16x8x2_t u16x8x2;
	uint32x4x2_t u32x4x2;
} __arm256i_composite;

typedef union arm128f_composite
{
	float32x2x2_t f32x2x2;
} __arm128f_composite;

typedef union arm256f_composite
{
	float32x4x2_t f32x4x2;
} __arm256f_composite;

#endif
