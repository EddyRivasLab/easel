#include "esl_config.h"

#ifdef HAVE_AVX512
#include <stdio.h>
#include <immintrin.h>		/* AVX2 */
#include "esl_avx_512.h"


 void esl_print_512(__m512i var){
  uint64_t *val = (uint64_t*) &var;
    printf("%016lx %016lx %016lx %016lx %016lx %016lx %016lx %016lx \n", 
           val[7], val[6], val[5], val[4], val[3], val[2], 
           val[1], val[0]);
 }

blert
#endif
