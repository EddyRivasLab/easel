/* i/o of multiple sequence alignment files in Clustal-like formats
 *
 * This module is responsible for i/o of:
 *    eslMSAFILE_CLUSTAL
 *    eslMSAFILE_MUSCLE
 *
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
#include "esl_msafile_clustal.h"
#include "esl_recorder.h"

static int make_text_consensus_line(const ESL_MSA *msa, char **ret_consline);
#ifdef eslAUGMENT_ALPHABET
static int make_digital_consensus_line(const ESL_MSA *msa, char **ret_consline);
#endif



/* Function:  esl_msafile_clustal_Read()
 * Synopsis:  Read in a CLUSTAL or CLUSTAL-like alignment.
 *
 * Purpose:   Read an open <ESL_MSAFILE> <afp>, starting from the
 *            position of its current point, parsing for CLUSTAL-like
 *            formats. (<afp->format> is expected to be
 *            <eslMSAFILE_CLUSTAL> or <eslMSAFILE_MUSCLE>.) Create a
 *            new multiple alignment, and return a ptr to that
 *            alignment in <*ret_msa>.  Caller is responsible for
 *            free'ing this <ESL_MSA>.
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
 *            and <afp->errmsg> is set to a message as above.
 *
 * Throws:    <eslEMEM> - an allocation failed.
 *            <eslESYS> - a system call such as fread() failed
 *            <eslEINVAL> - anchoring call failed in esl_buffer code
 *            <eslEINCONCEIVABLE> - "impossible" corruption 
 */
int
esl_msafile_clustal_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa     = esl_msa_Create(16, -1); /* a growable MSA. */
  char     *p       = NULL;
  esl_pos_t n       = 0;
  int       nblocks = 0;
  int       idx     = 0;
  esl_pos_t alen    = 0;
  esl_pos_t pos;
  esl_pos_t name_start, name_len;
  esl_pos_t seq_start, seq_len;
  esl_pos_t block_seq_start, block_seq_len;
  int       status;

  afp->errmsg[0] = '\0';
  if (! msa) { status = eslEMEM; goto ERROR; }

  /* skip leading blank lines in file */
  do {
    if ( (status = esl_buffer_GetLine(afp->bf, &p, &n)) != eslOK) goto ERROR; /* includes EOF */
    if (afp->linenumber != -1) afp->linenumber++;
  } while (esl_memspn(p, n, " \t\r\n") == n); /* idiomatic for "blank line" */
  /* now p[0..n-1] is the first non-blank line; point is at the start of the next line. */
    
  /* That first line says something like: "CLUSTAL W (1.83) multiple sequence alignment" */
  switch (afp->format) {
  case eslMSAFILE_CLUSTAL: if (! esl_memstrpfx(p, n, "CLUSTAL")) ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing CLUSTAL header"); break;
  case eslMSAFILE_MUSCLE:  if (! esl_memstrpfx(p, n, "MUSCLE"))  ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing MUSCLE header");  break;
  default:                 ESL_XEXCEPTION(eslEINCONCEIVABLE, "format %s is not clustal-like", esl_msa_DecodeFormat(afp->format));
  }

  /* skip blank lines again */
  do {
    status = esl_buffer_GetLine(afp->bf, &p, &n);
    if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "no alignment data following header");
    else if (status != eslOK) goto ERROR;
    if (afp->linenumber != -1) afp->linenumber++;
  } while (esl_memspn(p, n, " \t\r\n") == n); /* idiom for "blank line" */

  /* Read the file a line at a time. */
  nblocks  = 0;
  do {  /* p, n is now the first line of a block. for each line in a block of lines: */
    idx = 0;
    do {
      for (pos = 0;     pos < n; pos++) if (! isspace(p[pos])) break;  name_start = pos; 
      for (pos = pos+1; pos < n; pos++) if (  isspace(p[pos])) break;  name_len   = pos - name_start;
      for (pos = pos+1; pos < n; pos++) if (! isspace(p[pos])) break;  seq_start  = pos;      
      if (pos >= n) ESL_XFAIL(eslEFORMAT, afp->errmsg, "invalid alignment line");
      for (pos = n-1; pos > 0; pos--)   if (! isspace(p[pos])) break;  seq_len    = pos - seq_start - 1;

      if (idx == 0) {
	block_seq_start = seq_start;
	block_seq_len   = seq_len;
      } else {
	if (seq_start != block_seq_start) ESL_XFAIL(eslEFORMAT, afp->errmsg, "sequence start is misaligned");
	if (seq_len   != block_seq_len)   ESL_XFAIL(eslEFORMAT, afp->errmsg, "sequence end is misaligned");
      }

      /* Store the sequence name. */
      if (nblocks == 0)	{
	/* make sure we have room for another sequence */
	if (idx >= msa->sqalloc &&  (status = esl_msa_Expand(msa))      != eslOK) goto ERROR;
	if ( (status = esl_memstrdup(p+name_start, name_len, &(msa->sqname[idx]))) != eslOK) goto ERROR;
	msa->nseq++;
      } else {
	if (! esl_memstrcmp(p+name_start, name_len, msa->sqname[idx]))
	  ESL_XFAIL(eslEFORMAT, afp->errmsg, "expected sequence %s on this line, but saw %.*s", msa->sqname[idx], (int) name_len, p+name_start);
      }

      /* Append the sequence. */
#ifdef eslAUGMENT_ALPHABET
      if (msa->flags & eslMSA_DIGITAL)
	{ status = esl_abc_dsqcat(msa->abc, &(msa->ax[idx]), &(msa->sqlen[idx]), p+seq_start, seq_len); }
#endif
      if (! (msa->flags & eslMSA_DIGITAL)) 
	{ 
	  status = esl_strcat(&(msa->aseq[idx]), msa->sqlen[idx], p+seq_start, seq_len);
	  msa->sqlen[idx] += seq_len;
	}

      /* get next line. if it's a consensus line, we're done with the block */
      status = esl_buffer_GetLine(afp->bf, &p, &n);
      if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "alignment block did not end with consensus line");
      else if (status != eslOK)  goto ERROR;
      if (afp->linenumber != -1) afp->linenumber++;

      idx++;
    } while (esl_memspn(p, n, " .:*") < n); /* end loop over a block */
    
    if (idx != msa->nseq) ESL_XFAIL(eslEFORMAT, afp->errmsg, "last block didn't contain same # of seqs as earlier blocks");

    /* skip blank lines until we find start of next block, or EOF */
    do {
      status = esl_buffer_GetLine(afp->bf, &p, &n);
      if      (status == eslEOF) break;
      else if (status != eslOK)  goto ERROR;
      if (afp->linenumber != -1) afp->linenumber++;
    } while (esl_memspn(p, n, " \t\r\n") == n); 
    
    alen += block_seq_len;
    nblocks++;
  } while (status == eslOK);	/* normal end has status == EOF after last block. */

  msa->alen = alen;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}  






/* Function:  esl_msafile_clustal_Write()
 * Synopsis:  Write a CLUSTAL format alignment file to a stream.
 *
 * Purpose:   Write alignment <msa> in CLUSTAL W 1.83 format to
 *            output stream <fp>. 
 *            
 *            The alignment is written in blocks of 60 aligned
 *            residues at a time.
 *            
 * Args:      fp  - open output stream
 *            msa - alignment to write      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msafile_clustal_Write(FILE *fp, const ESL_MSA *msa)
{
  int   i;
  char *consline = NULL;
  char  buf[61];
  int   maxnamelen = 0;
  int   namelen;
  esl_pos_t apos;
  int   status;


  /* Find the maximum name length; determines width of name block  */
  for (i = 0; i < msa->nseq; i++)
    {
      namelen = strlen(msa->sqname[i]);
      maxnamelen = ESL_MAX(namelen, maxnamelen);
    }

  /* Make a CLUSTAL-like consensus line */
#ifdef eslAUGMENT_ALPHABET 
  if (msa->flags & eslMSA_DIGITAL) 
    status = make_digital_consensus_line(msa, &consline);
  else 
#endif
  status = make_text_consensus_line(msa, &consline);
  if (status != eslOK) goto ERROR;


  /* The magic header */
  fprintf(fp, "CLUSTAL W (1.83) multiple sequence alignment\n");


  /* The alignment */
  buf[60] = '\0';
  for (apos = 0; apos < msa->alen; apos += 60)
    {
      fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
#ifdef eslAUGMENT_ALPHABET 
	  if (msa->flags & eslMSA_DIGITAL) 
	    esl_abc_TextizeN(msa->abc, msa->ax[i]+apos+1, 60, buf);
	  else
#endif
	  strncpy(buf, msa->aseq[i]+apos, 60);
	  fprintf(fp, "%-*s %s\n", maxnamelen, msa->sqname[i], buf);
	}
      strncpy(buf, consline+apos, 60);
      fprintf(fp, "%-*s %s\n", maxnamelen, "", buf);
    }

  free(consline);
  return eslOK;

 ERROR:
  if (consline != NULL) free(consline);
  return status;
}


/* make_text_consensus_line()
 * 
 * Given a text mode <msa>, allocate and create a CLUSTAL-style
 * consensus line; return it in <*ret_consline>. Caller is responsible
 * for free'ing this strong.
 * 
 * The consensus line is numbered 0..alen-1, and is NUL-terminated.
 * 
 * Currently this only does a subset of what CLUSTAL consensus lines
 * look like; it only uses '*' for completely conserved positions,
 * and elsewise uses ' '.
 * 
 * Returns <eslOK> on success.
 * No normal failure codes.
 * Throws <eslEMEM> on allocation error.
 */
static int
make_text_consensus_line(const ESL_MSA *msa, char **ret_consline)
{
  char  *consline = NULL;
  int   *ct       = NULL;
  int    i, apos, x;
  int    nseen;
  int    status;

  ESL_ALLOC(consline, sizeof(char) * (msa->alen+1));
  ESL_ALLOC(ct,       sizeof(int)  * 27);

  for (x = 0; x <= 26; x++) ct[x] = 0;
				      
  for (apos = 0; apos < msa->alen; apos++)
    {
      for (i = 0; i < msa->nseq; i++)
	{
	  x = toupper(msa->aseq[i][apos]) - 'A';
	  if (x >= 0 && x < 26) ct[x]++;
	  else                  ct[26]++;
	}
	  
      for (nseen = 0, x = 0; x < 26; x++) /* not including gaps */
	if (ct[x] > 0) nseen++;

      if (nseen == 1) consline[apos] = '*';
      else            consline[apos] = ' ';
    }

  consline[msa->alen] = '\0';
  *ret_consline = consline;
  free(ct);
  return eslOK;

 ERROR:
  if (ct != NULL) free(ct);
  *ret_consline = NULL;
  return status;
}


/* make_digital_consensus_line()
 * 
 * Exactly the same as make_text_consensus_line(), except for
 * digital mode <msa>.
 */
#ifdef eslAUGMENT_ALPHABET
static int
make_digital_consensus_line(const ESL_MSA *msa, char **ret_consline)
{
  char  *consline = NULL;
  float *ct       = NULL;
  int    i, apos, x;
  int    nseen;
  int    status;

  ESL_ALLOC(consline, sizeof(char)  * (msa->alen+1));
  ESL_ALLOC(ct,       sizeof(float) * (msa->abc->K+1));  

  for (x = 0; x <= msa->abc->K+1; x++) ct[x] = 0.;
				      
  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (i = 0; i < msa->nseq; i++)
	esl_abc_FCount(msa->abc, ct, msa->ax[i][apos], 1.0);
	
      for (nseen = 0, x = 0; x < msa->abc->K; x++) /* not including gaps */
	if (ct[x] > 0.) nseen++;
      
      if (nseen == 1) consline[apos-1] = '*';
      else            consline[apos-1] = ' ';
    }

  consline[msa->alen] = '\0';
  *ret_consline = consline;
  free(ct);
  return eslOK;

 ERROR:
  if (ct != NULL) free(ct);
  *ret_consline = NULL;
  return status;
}
#endif /*eslAUGMENT_ALPHABET*/


/*****************************************************************
 * Example.
 *****************************************************************/

#ifdef eslMSAFILE_CLUSTAL_EXAMPLE
/* An example of reading an MSA in text mode, and handling any returned errors.
   gcc -g -Wall -o esl_msafile_clustal_example -I. -DeslMSAFILE_CLUSTAL_EXAMPLE esl_msafile_clustal.c esl_msa.c easel.c 
   ./esl_msafile_clustal_example <msafile>
 */

/*::cexcerpt::msafile_clustal_example::begin::*/
#include <stdio.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_clustal.h"

int 
main(int argc, char **argv)
{
  char        *filename = argv[1];
  int          fmt      = eslMSAFILE_UNKNOWN;
  ESLX_MSAFILE *afp      = NULL;
  ESL_MSA     *msa      = NULL;
  int          status;

  if ( (status = eslx_msafile_Open(filename, fmt, NULL, &afp)) != eslOK) 
    eslx_msafile_OpenFailure(afp, status);

  if ( (status = esl_msafile_clustal_Read(afp, &msa))         != eslOK)
    eslx_msafile_ReadFailure(afp, status);

  printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	 1, msa->name, msa->nseq, (int) msa->alen);

  esl_msafile_clustal_Write(stdout, msa);
  esl_msa_Destroy(msa);
  eslx_msafile_Close(afp);
  exit(0);
}
/*::cexcerpt::msafile_clustal_example::end::*/
#endif /*eslMSAFILE_CLUSTAL_EXAMPLE*/
/*--------------------- end of example --------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
