/* Clustering sequences in an MSA by % identity.
 */
#ifndef eslMSA_ISET_INCLUDED
#define eslMSA_ISET_INCLUDED
#include "esl_config.h"
#include "esl_msa.h"
#include "esl_random.h"

extern int esl_msa_iset_Cobalt(const ESL_MSA *msa, double maxid, 
			     int **opt_c, int **opt_nin, ESL_RANDOMNESS *r);

#endif /*eslMSA_ISET_INCLUDED*/
