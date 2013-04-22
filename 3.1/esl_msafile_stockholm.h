/* I/O of multiple alignment files in Stockholm/Pfam format
 */
#ifndef eslMSAFILE_STOCKHOLM_INCLUDED
#define eslMSAFILE_STOCKHOLM_INCLUDED

extern int esl_msafile_stockholm_SetInmap     (ESLX_MSAFILE *afp);
extern int esl_msafile_stockholm_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
extern int esl_msafile_stockholm_Read         (ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_stockholm_Write        (FILE *fp, const ESL_MSA *msa, int fmt);

#endif /*eslMSAFILE_STOCKHOLM_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
