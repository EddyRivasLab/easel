#include "esl_config.h"
#ifdef eslENABLE_AVX512

#include <stdio.h>
#include <x86intrin.h>		
#include "esl_avx512.h"

void 
esl_avx512_dump_512i_hex8(__m512i v)
{
  uint64_t *val = (uint64_t*) &v;
  printf("%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", 
	 val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

#else  // ! eslENABLE_AVX512
void esl_avx512_silence_hack(void) { return; }
#endif // eslENABLE_AVX512
