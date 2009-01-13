/* esl_msaweight.h
 * Sequence weighting algorithms.
 * 
 * SVN $Id$
 * SRE, Sun Nov  5 09:11:13 2006 [Janelia]
 */
#ifndef ESL_MSAWEIGHT_INCLUDED
#define ESL_MSAWEIGHT_INCLUDED

#include <esl_msa.h>

extern int esl_msaweight_GSC(ESL_MSA *msa);
extern int esl_msaweight_PB(ESL_MSA *msa);
extern int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid);



#endif /*ESL_MSAWEIGHT_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
