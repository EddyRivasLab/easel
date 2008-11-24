/* Vectorized routines for Intel/AMD, using Streaming SIMD Extensions (SSE).
 * 
 * SRE, Sun Dec 16 10:01:41 2007 [Janelia]
 * SVN $Id$
 */
#ifdef HAVE_SSE2
#ifndef ESL_SSE_INCLUDED
#define ESL_SSE_INCLUDED

#include "easel.h"

#include <stdio.h>
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */


extern __m128 esl_sse_logf(__m128 x);
extern __m128 esl_sse_expf(__m128 x);
extern __m128 esl_sse_select_ps(__m128 a, __m128 b, __m128 mask);
extern void   esl_sse_dump_ps(FILE *fp, __m128 v);

#endif /*ESL_SSE_INCLUDED*/
#endif /*HAVE_SSE2*/
