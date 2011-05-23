/* i/o of multiple sequence alignment files in "aligned FASTA" format
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_mem.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_afa.h"


/* Function:  esl_msafile_afa_SetInmap()
 * Synopsis:  Finishes configuring input map for aligned FASTA format.
 *
 * Purpose:   We tolerate spaces in input lines of aligned FASTA format;
 *            map ' ' to <eslDSQ_IGNORED>.
 */
int
esl_msafile_afa_SetInmap(ESLX_MSAFILE *afp)
{
  afp->inmap[' '] = eslDSQ_IGNORED;
  return eslOK;
}

/* Function:  esl_msafile_afa_Read()
 * Synopsis:  Read in an aligned FASTA format alignment.
 *
 * Purpose:   Read an MSA from an open <ESLX_MSAFILE> <afp>,
 *            parsing for aligned FASTA format. Create
 *            a new MSA, and return a ptr to that alignment
 *            in <*ret_msa>. Caller is responsible for free'ing 
 *            this <ESL_MSA>.
 *
 * Args:      afp     - open <ESL_MSAFILE>
 *            ret_msa - RETURN: newly parsed <ESL_MSA>
 *
 * Returns:   <eslOK> on success.
 * 
 *            In the event of a parse error, returns <eslEFORMAT>, and
 *            set <afp->errmsg> to an appropriately informative error
 *            message that can be shown to the user. 
 *
 *            If no alignment is found at all, returns <eslEOF>,
 *            and <afp->errmsg> is blank.
 *
 * Throws:    <eslEMEM> - an allocation failed.
 *            <eslESYS> - a system call such as fread() failed
 *            <eslEINCONCEIVABLE> - "impossible" corruption 
 */
int
esl_msafile_afa_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa  = NULL;
  int       idx  = 0;
  int64_t   alen = 0;
  int64_t   this_alen = 0;
  char     *p, *tok;
  esl_pos_t n, ntok;
  int       status;

  afp->errmsg[0] = '\0';	                                  /* Blank the error message. */
  if (afp->msa_cache) return eslx_msafile_Decache(afp, ret_msa);  /* Check the <afp>'s cache first */

  if (afp->abc   &&  (msa = esl_msa_CreateDigital(afp->abc, 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
  if (! afp->abc &&  (msa = esl_msa_Create(                 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
  
  do {   /* skip leading blank lines until we see a > name/desc line */ 
    if ( (status = esl_buffer_GetLine(afp->bf, &p, &n)) != eslOK) goto ERROR; /* includes EOF */
    if (afp->linenumber != -1) afp->linenumber++;
  } while (esl_memspn(p, n, " \t\r\n") == n); /* idiomatic for "blank line" */
  while (n && isspace(*p)) { p++; n--; }     /* tolerate sloppy space at line start*/

  do {
    if (n <= 1 || *p != '>' ) ESL_XFAIL(eslEFORMAT, afp->errmsg, "expected aligned FASTA name/desc line starting with >");    
    p++; n--;			/* advance past > */

    if ( (status = esl_memtok(&p, &n, " \t", &tok, &ntok)) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errmsg, "no name found for aligned FASTA record");
    if (idx >= msa->sqalloc && (status = esl_msa_Expand(msa))    != eslOK) goto ERROR;

    if (     (status = esl_msa_SetSeqName       (msa, idx, tok, ntok)) != eslOK) goto ERROR;
    if (n && (status = esl_msa_SetSeqDescription(msa, idx, p,   n))    != eslOK) goto ERROR;

    /* The code below will do a realloc on every line. Possible optimization: once you know
     * alen (from first sequence), allocate subsequent seqs once, use noalloc versions of
     * esl_strmapcat/esl_abc_dsqcat(). Requires implementing protection against overrun, if
     * input is bad and a sequence is too long. Could gain ~25% or so that way (quickie
     * test on PF00005 Full)
     */
    this_alen = 0;
    while ((status = esl_buffer_GetLine(afp->bf, &p, &n)) == eslOK)
      {
	if (afp->linenumber != -1) afp->linenumber++;
	while (n && isspace(*p)) { p++; n--; } /* skip leading whitespace on line */
	if (n  == 0)   continue;
	if (*p == '>') break;

#ifdef eslAUGMENT_ALPHABET
	if (msa->abc)   { status = esl_abc_dsqcat(afp->inmap, &(msa->ax[idx]),   &this_alen, p, n); } 
#endif
	if (! msa->abc) { status = esl_strmapcat (afp->inmap, &(msa->aseq[idx]), &this_alen, p, n); }
	if (status == eslEINVAL)   ESL_XFAIL(eslEFORMAT, afp->errmsg, "one or more invalid sequence characters");
	else if (status != eslOK)  goto ERROR;
      }
    if (alen && alen != this_alen) ESL_XFAIL(eslEFORMAT, afp->errmsg, "sequence %s has alen %" PRId64 "; expected %" PRId64, msa->sqname[idx], this_alen, alen);

    alen = this_alen;
    idx++;
  } while (status == eslOK);	/* normally ends on eslEOF. */

  msa->nseq = idx;
  msa->alen = alen;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;

}

/* Function:  esl_msafile_afa_Write()
 * Synopsis:  Write an aligned FASTA format alignment file to a stream.
 *
 * Purpose:   Write alignment <msa> in aligned FASTA format to a stream
 *            <fp>.
 * 
 *            If <msa> is in text mode, residues and gaps are written
 *            exactly as they appear in the data structure. If <msa>
 *            is digital, residues are in uppercase and all gaps are
 *            -.
 *
 * Args:      fp  - open stream to write to
 *            msa - MSA to write
 *
 * Returns:   <eslOK> on success.
 */
int
esl_msafile_afa_Write(FILE *fp, const ESL_MSA *msa)
{
  int     i;
  int64_t pos;
  char    buf[61];
  int     acpl;       /* actual number of character per line */
  
  for (i = 0; i < msa->nseq; i++)
    {
      fprintf(fp, ">%s", msa->sqname[i]);
      if (msa->sqacc  != NULL && msa->sqacc[i]  != NULL) fprintf(fp, " %s", msa->sqacc[i]);
      if (msa->sqdesc != NULL && msa->sqdesc[i] != NULL) fprintf(fp, " %s", msa->sqdesc[i]);
      fputc('\n', fp);

      pos = 0;
      while (pos < msa->alen)
	{
	  acpl = (msa->alen - pos > 60)? 60 : msa->alen - pos;
#ifdef eslAUGMENT_ALPHABET
	  if (msa->abc)   esl_abc_TextizeN(msa->abc, msa->ax[i] + pos + 1, acpl, buf);
#endif
	  if (! msa->abc) strncpy(buf, msa->aseq[i] + pos, acpl);

	  buf[acpl] = '\0';
	  fprintf(fp, "%s\n", buf);	      
	  pos += 60;
	}
    } 
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
