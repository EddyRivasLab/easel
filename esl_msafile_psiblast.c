/* I/O of multiple sequence alignments in PSI-BLAST format
 * 
 * Contents:
 *   1. API for reading/writing PSI-BLAST format
 *   
 */
#include "esl_config.h"

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#include "esl_mem.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_psiblast.h"

/*****************************************************************
 *# 1. API for reading/writing PSI-BLAST format
 *****************************************************************/

/* Function:  esl_msafile_psiblast_Write()
 * Synopsis:  Write an MSA to a stream in PSI-BLAST format
 *
 * Purpose:   Write alignment <msa> in NCBI PSI-BLAST format to 
 *            stream <fp>.
 *            
 *            The <msa> should have a valid reference line <msa->rf>,
 *            with alphanumeric characters marking consensus (match)
 *            columns, and non-alphanumeric characters marking
 *            nonconsensus (insert) columns. 
 *            
 *            As a fallback, if the <msa> does not have a reference
 *            line, a default one is created. In this case,
 *            <esl_msa_MarkFragments()> will be called with
 *            <fragthresh=0.5> to heuristically define sequence
 *            fragments in the alignment, then
 *            <esl_msa_ReasonableRF()> will be called with
 *            <symfrac=0.5> to mark consensus columns.  This will
 *            modify the alignment data in <msa> (converting
 *            leading/trailing gap symbols on "fragments" to '~'
 *            missing data symbols) in addition to adding <#=RF>
 *            annotation to it. If caller doesn't want the alignment
 *            data to be modified, it must provide its own <msa->rf>
 *            line, or provide a copy of <msa>, or it must avoid
 *            writing alignments in PSI-BLAST format.
 *            
 *            PSI-BLAST format allows only one symbol ('-') for gaps,
 *            and cannot represent missing data symbols (Easel's
 *            '~'). Any missing data symbols are converted to gaps.
 *
 * Args:      fp  - open output stream
 *            msa - MSA to write       
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msafile_psiblast_Write(FILE *fp, const ESL_MSA *msa)
{
  int     cpl = 60;
  char    buf[61];
  int     i;
  int     sym;
  int     pos;
  int     bpos;
  int     acpl;
  int     maxnamewidth;
  int     is_consensus;
  int     is_residue;
  int     status;

  maxnamewidth = (int) maxwidth(msa->sqname, msa->nseq);

  /* if <msa> lacks an RF line, make a default one */
  if (msa->rf == NULL) 
    {
      ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
      if ((status = esl_msa_MarkFragments(msa, 0.5))          !=  eslOK) return status;
      if ((status = esl_msa_ReasonableRF (msa, 0.5, msa->rf)) !=  eslOK) return status;
    }

  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      for (i = 0; i < msa->nseq; i++)
	{
	  acpl =  (msa->alen - pos > cpl)? cpl : msa->alen - pos;

#ifdef eslAUGMENT_ALPHABET
	  if (msa->flags & eslMSA_DIGITAL)
	    {
	      for (bpos = 0; bpos < acpl; bpos++)
		{
		  sym          = msa->abc->sym[msa->ax[i][pos + bpos + 1]];
		  is_consensus = (isalnum(msa->rf[pos + bpos + 1]) ? TRUE : FALSE);
		  is_residue   = esl_abc_XIsResidue(msa->abc, msa->ax[i][pos+bpos+1]);
				      
		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos] = toupper(sym);
		      else            buf[bpos] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if (is_residue) buf[bpos] = tolower(sym);
		      else            buf[bpos] = '-';     /* includes conversion of missing data to - */
		    }
		}
	    }
#endif
	  if (! (msa->flags & eslMSA_DIGITAL))
	    {
	      for (bpos = 0; bpos < acpl; bpos++)
		{
		  sym          = msa->aseq[i][pos + bpos];
		  is_consensus = (isalnum(msa->rf[pos + bpos]) ? TRUE : FALSE);
		  is_residue   = isalnum(sym);

		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos] = toupper(sym);
		      else            buf[bpos] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if (is_residue) buf[bpos] = tolower(sym);
		      else            buf[bpos] = '-';     /* includes conversion of missing data to - */
		    }
		}
	    }
	  buf[acpl] = '\0';	      
	  fprintf(fp, "%-*s  %s\n", maxnamewidth, msa->sqname[i], buf);
	}  /* end loop over sequences */

      if (pos + cpl < msa->alen) fputc('\n', fp);
    } /* end loop over alignment blocks */
  return eslOK;
  
 ERROR:
  return status;
}


}
/*----------- end, API for i/o of psi-blast format --------------*/


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
