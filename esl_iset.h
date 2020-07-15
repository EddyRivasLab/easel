/* Generalized single linkage clustering.
 * 
 * SRE, Mon Jan  7 09:40:06 2008 [Janelia]
 */
#ifndef eslISET_INCLUDED
#define eslISET_INCLUDED
#include "esl_config.h"
#include "esl_random.h"

extern int esl_iset_Cobalt(void *base, size_t n, size_t size, 
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *workspace, int *assignments, ESL_RANDOMNESS *r);


extern int shuffle_array(ESL_RANDOMNESS *r, int a[], int n);
#endif /*eslISET_INCLUDED*/