/* i/o of multiple sequence alignment files in Clustal-like formats
 *
 * This module is responsible for i/o of:
 *    eslMSAFILE_CLUSTAL
 *    eslMSAFILE_MUSCLE
 *
 * SRE, Thu Mar 11 12:37:44 2010 [UA 916 Seattle->Dulles]
 * SVN $Id$
 * SVN $URL$
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_clustal.h"
#include "esl_recorder.h"

static int make_text_consensus_line(ESL_MSA *msa, char **ret_consline);
#ifdef eslAUGMENT_ALPHABET
static int make_digital_consensus_line(ESL_MSA *msa, char **ret_consline);
#endif

/* Function:  esl_msafile_clustal_Read()
 * Synopsis:  Read in a CLUSTAL or CLUSTAL-like alignment.
 * Incept:    SRE, Thu Mar 11 10:36:31 2010 [UA 916 from Seattle]
 *
 * Purpose:   Read an open <ESL_MSAFILE> <afp>, parsing for CLUSTAL-like
 *            formats (<afp->format> is <eslMSAFILE_CLUSTAL> or 
 *            <eslMSAFILE_MUSCLE>), and create a new multiple
 *            alignment; return a ptr to that alignment in <*ret_msa>.
 *            Caller is responsible for free'ing this <ESL_MSA>.
 *
 * Args:      afp     - open <ESL_MSAFILE>
 *            ret_msa - RETURN: newly parsed <ESL_MSA>
 *
 * Returns:   <eslOK> on success.
 * 
 *            In the event of a parse error, returns <eslEFORMAT>, and
 *            set <afp->errbuf> to an appropriately informative error
 *            message that can be shown to the user. These messages
 *            look like "parse failed (line %d): blah blah blah" with
 *            no trailing newline. The caller may prefix (with the
 *            filename, perhaps) or suffix this message as it pleases.
 *
 *            If no alignment is found at all, returns <eslEOF>,
 *            and <afp->errbuf> is set to a message as above.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_msafile_clustal_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa   = NULL;
  char     *buf   = NULL;
  char    **lines = NULL;
  int       n;
  int       nblocks;
  char     *tok1 = NULL;	/* first token on line  = name     */
  char     *tok2 = NULL;	/* second token on line = subseq   */
  char      len1, len2;		/* lengths of the two tokens       */
  int       li;			/* line counter in block[0..n-1]   */
  int       startline;		/* line number of start of block   */
  int       leftpos_each_line;	/* position of first aligned chars */
  int       nres_each_line;	/* # of residues on each line      */
  int       status;

  afp->errbuf[0] = '\0';

  /* Set <buf> to the first nonblank line. */
  do { status = esl_recorder_Read(afp->rc, &buf); } while (status == eslOK && esl_str_IsBlank(buf));
  if (status == eslEOF)  ESL_XFAIL(eslEOF, afp->errbuf, "parse failed (line %d): end of file; no alignment?", 1+esl_recorder_GetCurrent(afp->rc));
  if (status == eslEMEM) return status;

  /* That first line says something like: "CLUSTAL W (1.83) multiple sequence alignment" */
  if (afp->format == eslMSAFILE_CLUSTAL) 
    {
      status = strncmp(buf, "CLUSTAL", 7);
      if (status) ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing CLUSTAL header", 1+esl_recorder_GetCurrent(afp->rc));
    }
  else if (afp->format == eslMSAFILE_MUSCLE) 
    {
      status = strncmp(buf, "MUSCLE",  6);
      if (status)  ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing MUSCLE header", 1+esl_recorder_GetCurrent(afp->rc));
    }
  else ESL_XEXCEPTION(eslEINCONCEIVABLE, "format %d is not clustal-like", afp->format);

  /* Loop over all alignment blocks: one block parsed at a time. */
  nblocks = 0;
  status  = eslOK;
  while (status == eslOK)
    {
      /* skip ahead to first line of block */
      do { status = esl_recorder_Read(afp->rc, &buf); } while (status == eslOK && esl_str_IsBlank(buf));
      if (status == eslEOF)  ESL_XFAIL(eslEOF, afp->errbuf, "parse failed (line %d): end of file; no alignment?", 1+esl_recorder_GetCurrent(afp->rc));
      if (status == eslEMEM) return status;
  
      /* read block */
      startline = esl_recorder_GetCurrent(afp->rc);
      esl_recorder_MarkBlock(rc, startline);
      do { status = esl_recorder_Read(afp->rc, &buf); } while (status == eslOK && ! esl_str_IsBlank(buf));
      if (status == eslEMEM) return status; /* both OK and EOF is acceptable here */
      esl_recorder_GetBlock(rc, &block, NULL, NULL, &n);
      nblocks++;
      if (status == eslOK) n--;	/* i.e. status of last _Read(), EOF vs. OK; if OK, last line is blank */

      /* Now we have a block of <n> lines, block[0]..block[n-1]. 
       * In CLUSTAL format, the last line [n-1] is a consensus line,
       * so there should be n-1 aligned sequences.
       */
      if (nblocks == 1) 
	{
	  if ((msa = esl_msa_Create(n-1, -1)) == NULL) { status = eslEMEM; goto ERROR; }
	  msa->nseq = n-1;
	  msa->alen = 0;
	}
      else if (n-1 != msa->nseq)
	ESL_XFAIL(eslEFORMAT, afp->errbuf, 
		  "parse failed (block %d starting at line %d): expected %d seqs in block, saw %d",
		  nblocks, startline+1, msa->nseq, n-1);

      for (li = 0; li < msa->nseq; li++)
	{
	  s = block[li];
	  status = esl_strtok_adv(&s, " \t",     &tok1, &len1, NULL);
	  if (status != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, 
					 "parse failed (line %d): expected to find a seq name but didn't", startline+li+1);

	  status = esl_strtok_adv(&s, " \t\n\r", &tok2, &len2, NULL);
	  if (status != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, 
					 "parse failed (line %d): expected to find aligned sequence but didn't", startline+li+1);
      
	  for (; *s; s++)
	    if (! isspace(*s))
	      ESL_XFAIL(eslEFORMAT, afp->errbuf, 
			"parse failed (line %d) : expected only name/seq on line", 
			startline+li+1);

	  if (li == 0) leftpos_each_line = tok-block[i];
	  else if (leftpos_each_line != tok-block[i]) 
	    ESL_XFAIL(eslEFORMAT, afp->errbuf,
		      "parse failed (line %d; ali block %d) : aligned seqs aren't flush",
		      startline+li+1, nblocks);
      
	  if (li == 0) nres_each_line = len2;
	  else if (nres_each_line != len2) 
	    ESL_XFAIL(eslEFORMAT, afp->errbuf, 
		      "parse failed (line %d, ali block %d) : expected %d residues on line, saw %d", 
		      startline+li, nblocks, nres_each_line, len2);
      
	  if (nblock == 1) if ((status = esl_strdup(tok, ntok, &(msa->sqname[li]))) != eslOK) goto ERROR;
	  else if (esl_strcmp(tok, msa->sqname[li]) != 0) 
	    ESL_XFAIL(eslEFORMAT, afp->errbuf,
		      "parse failed (line %d, ali block %d) : expected seq name %s here, saw %s", 
		      startline+li+1, nblocks, msa->sqname[li], tok1);

	  ESL_RALLOC(msa->aseq[i], tmpp, sizeof(char) * (msa->alen + len2 + 1)); 	  
	  memcpy(msa->aseq[i]+msa->alen, tok2, sizeof(char) * len2);
	  msa->aseq[i][msa->alen+len2] = '\0';
	} /* end loop over lines in a block */
      msa->alen += nres_each_line;

      /* if we were superparanoid, we could validate the consensus line block[n-1] here */

    } /* end of a block; ready to look for next block */

#ifdef eslAUGMENT_ALPHABET 
  if (afp->do_digital && (status = esl_msa_Digitize(afp->abc, msa, afp->errbuf)) != eslOK) goto ERROR;
#endif  

  afp->linenumber = esl_recorder_GetCurrent(afp->rc);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  afp->linenumber = esl_recorder_GetCurrent(afp->rc);
  esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}

/* Function:  esl_msafile_clustal_Write()
 * Synopsis:  Write a CLUSTAL format file to a stream
 * Incept:    SRE, Thu Mar 11 11:59:02 2010 [UA916 from Seattle]
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
  buf[60] = '\0'
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
      fprintf(fp, "%-*s %s\n", "", consline);
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
make_text_consensus_line(ESL_MSA *msa, char **ret_consline)
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


/* make_digital_consensus_line()
 * 
 * Exactly the same as make_text_consensus_line(), except for
 * digital mode <msa>.
 */
#ifdef eslAUGMENT_ALPHABET
static int
make_digital_consensus_line(ESL_MSA *msa, char **ret_consline)
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
#include "esl_msafile_clustal.h"

int 
main(int argc, char **argv)
{
  char        *filename = argv[1];
  int          fmt      = eslMSAFILE_CLUSTAL;
  ESL_MSAFILE *afp      = NULL;
  ESL_MSA     *msa      = NULL;
  int          status;

  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s not found or not readable\n", filename);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of %s\n",  filename);
  else if (status != eslOK)        esl_fatal("Alignment file open failed (error %d)\n", status);

  status = esl_msafile_clustal_Read(afp, &msa);
  if      (status == eslEFORMAT) esl_fatal("alignment file %s: %s\n",                    afp->fname, afp->errbuf);
  else if (status == eslEOF)     esl_fatal("alignment file %s appears empty?\n",         afp->fname);
  else if (status != eslOK)      esl_fatal("alignment file %s: read failed, error %d\n", afp->fname, status);

  printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	 nali, msa->name, msa->nseq, (int) msa->alen);

  esl_msafile_clustal_Write(stdout, msa);
  esl_msa_Destroy(msa);

  esl_msafile_Close(afp);
  exit(0);
}
/*::cexcerpt::msafile_clustal_example::end::*/
#endif /*eslMSAFILE_CLUSTAL_EXAMPLE*/
/*--------------------- end of example --------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
