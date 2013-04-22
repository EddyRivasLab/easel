/* I/O of multiple sequence alignments in PSI-BLAST format
 */
#ifndef eslMSAFILE_PSIBLAST_INCLUDED
#define eslMSAFILE_PSIBLAST_INCLUDED

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_psiblast_SetInmap     (ESLX_MSAFILE *afp);
extern int esl_msafile_psiblast_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
extern int esl_msafile_psiblast_Read         (ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_psiblast_Write        (FILE *fp, const ESL_MSA *msa);

#endif /* eslMSAFILE_PSIBLAST_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
