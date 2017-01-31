#include "esl_config.h"
#ifdef eslENABLE_AVX

#include <stdio.h>
#include <x86intrin.h>	

#include "esl_avx.h"
#include "easel.h"

void 
esl_avx_dump_256i_hex4(__m256i v)
{
  uint64_t *val = (uint64_t *) &v;
  printf("%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", 
	 val[3], val[2], val[1], val[0]);
}

#else  // !eslENABLE_AVX 
void esl_avx_silence_hack(void) { return; }
#endif // eslENABLE_AVX

