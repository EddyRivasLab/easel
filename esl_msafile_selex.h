/* i/o of multiple sequence alignment files in SELEX format
 */
#ifndef eslMSAFILE_SELEX_INCLUDED
#define eslMSAFILE_SELEX_INCLUDED

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_selex_SetInmap(ESLX_MSAFILE *afp);
extern int esl_msafile_selex_CheckFileFormat(ESL_BUFFER *bf);
extern int esl_msafile_selex_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_selex_Write(FILE *fp, const ESL_MSA *msa);

#endif /* eslMSAFILE_SELEX_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
