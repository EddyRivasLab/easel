#include <stdio.h>
#include <immintrin.h>		/* AVX2 */

#include "esl_avx.h"
#include "easel.h"
#ifdef HAVE_AVX2
 void esl_print_256(__m256i var){
  uint64_t *val = (uint64_t*) &var;
    printf("%016lx %016lx %016lx %016lx\n",   val[3], val[2], val[1], val[0]);
}
#endif

