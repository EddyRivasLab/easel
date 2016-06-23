#include <stdio.h>
#ifdef HAVE_AVX2
#include <immintrin.h>		/* AVX2 */
#endif
#include "esl_avx.h"


// This file is just a dummy target to make sure that the functions defined in esl_avx.h get included in the library version of hmmer