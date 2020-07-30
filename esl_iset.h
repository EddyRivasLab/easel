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


extern int esl_bi_iset_Cobalt(void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *workspace, int *assignments, int *ret_larger, ESL_RANDOMNESS *r);

extern int	esl_iset_Blue(void *base, size_t n, size_t size,
							  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
							  int *workspace, int *assignments, ESL_RANDOMNESS *r);

extern int esl_bi_iset_Blue(void *base, size_t n, size_t size,
											  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
											  int *workspace, int *assignments, int *ret_larger, ESL_RANDOMNESS *r);

static int shuffle_array(ESL_RANDOMNESS *r, int a[], int n);

static void print_array(int array[], int n);

static int check_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments);

static int check_bi_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments);

static int i_select(void *base, size_t n, size_t size, int k,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
         int *dec_o, int *label_o, int *status_d, int *to_add, int *ret_lta);

static void i_update_workspace(int *dec_o, int *label_o, int *status_d, int *to_add, int *assignments, size_t n, int *k, 
	int *lta, ESL_RANDOMNESS *r);

static void bi_update_workspace_blue(int *dec_o, int *label_o, int *status_d, int *to_add, int *elig, int *assignments, int n, 
	int *d, int *l, int *lta1, int *lta2, int *nb1, int *nb2, ESL_RANDOMNESS *r);

static int
update_2_elig(int j, void *base, int n, size_t size,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,int *label_o, int *status_d, int *to_add, int *elig, const int lta1);

static int bi_select_blue(void *base, int n, size_t size, 
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
         int *dec_o, int *label_o, int *status_d, int *to_add, int *elig, const int d, const int l, int *ret_lta1, int *ret_lta2);



#endif /*eslISET_INCLUDED*/
