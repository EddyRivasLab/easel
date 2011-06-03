/* i/o of multiple sequence alignment files in Clustal-like formats
 *
 * Contents:
 *   1. API for reading/writing Clustal and Clustal-like formats
 *   2. Example.
 *   3. Copyright and license information.
 *   
 * This module is responsible for i/o of both eslMSAFILE_CLUSTAL and
 * eslMSAFILE_CLUSTALLIKE alignment formats.
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
#include "esl_mem.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_clustal.h"


/*****************************************************************
 *# 1. API for reading/writing Clustal and Clustal-like formats
 *****************************************************************/

static int make_text_consensus_line(const ESL_MSA *msa, char **ret_consline);
#ifdef eslAUGMENT_ALPHABET
static int make_digital_consensus_line(const ESL_MSA *msa, char **ret_consline);
#endif


/* Function:  esl_msafile_clustal_SetInmap()
 * Synopsis:  Finishes configuring input map for CLUSTAL, CLUSTALLIKE formats.
 *
 * Purpose:   This is a no-op; the default input map is fine for
 *            clustal, clustallike formats. (We don't allow spaces 
 *            interior of input lines, for example.)
 */
int
esl_msafile_clustal_SetInmap(ESLX_MSAFILE *afp)
{
  return eslOK;
}


/* Function:  esl_msafile_clustal_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open Clustal MSA input.
 *
 * Purpose:   Guess the alpbabet of the sequences in open
 *            A2M format MSA file <afp>.
 *            
 *            On a normal return, <*ret_type> is set to <eslDNA>,
 *            <eslRNA>, or <eslAMINO>, and <afp> is reset to its
 *            original position.
 *
 * Args:      afp      - open A2M format MSA file
 *            ret_type - RETURN: <eslDNA>, <eslRNA>, or <eslAMINO>       
 *
 * Returns:   <eslOK> on success.
 *            <eslENOALPHABET> if alphabet type can't be determined.
 *            In either case, <afp> is rewound to the position it
 *            started at.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslESYS> on failures of fread() or other system calls
 */
int 
esl_msafile_clustal_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type)
{
  int       alphatype     = eslUNKNOWN;
  esl_pos_t anchor        = -1;
  int       threshold[3]  = { 500, 5000, 50000 }; /* we check after 500, 5000, 50000 residues; else we go to EOF */
  int       nsteps        = 3;
  int       step          = 0;
  int       nres          = 0;
  int       x;
  int64_t   ct[26];
  char     *p, *tok;
  esl_pos_t n,  toklen, pos;
  int       status;

  for (x = 0; x < 26; x++) ct[x] = 0;

  anchor = esl_buffer_GetOffset(afp->bf);
  if ((status = esl_buffer_SetAnchor(afp->bf, anchor)) != eslOK) { status = eslEINCONCEIVABLE; goto ERROR; } /* [eslINVAL] can't happen here */

  /* Ignore the first nonblank line, which says "CLUSTAL W (1.83) multiple sequence alignment" or some such */
  while ( (status = eslx_msafile_GetLine(afp, &p, &n)) == eslOK  && esl_memspn(afp->line, afp->n, " \t") == afp->n) ;
  if      (status == eslEOF) ESL_XFAIL(eslENOALPHABET, afp->errmsg, "can't determine alphabet: no alignment data found");
  else if (status != eslOK)  goto ERROR;
  
  while ( (status = eslx_msafile_GetLine(afp, &p, &n)) == eslOK)
    {
      if ((status = esl_memtok(&p, &n, " \t", &tok, &toklen)) != eslOK) continue; /* ignore blank lines */
      /* p now points to the rest of the sequence line, after a name */
      
      /* count characters into ct[] array */
      for (pos = 0; pos < n; pos++)
	if (isalpha(p[pos])) {
	  x = toupper(p[pos]) - 'A';
	  ct[x]++;
	  nres++; 	  
	}

      /* try to stop early, checking after 500, 5000, and 50000 residues: */
      if (step < nsteps && nres > threshold[step]) {
	if ((status = esl_abc_GuessAlphabet(ct, &alphatype)) == eslOK) goto DONE; /* (eslENOALPHABET) */
	step++;
      }
    }
  if (status != eslEOF) goto ERROR; /* [eslEMEM,eslESYS,eslEINCONCEIVABLE] */
  status = esl_abc_GuessAlphabet(ct, &alphatype); /* (eslENOALPHABET) */

 DONE:
  esl_buffer_SetOffset(afp->bf, anchor);   /* Rewind to where we were. */
  esl_buffer_RaiseAnchor(afp->bf, anchor);
  *ret_type = alphatype;
  return status;

 ERROR:
  if (anchor != -1) {
    esl_buffer_SetOffset(afp->bf, anchor);
    esl_buffer_RaiseAnchor(afp->bf, anchor);
  }
  *ret_type = eslUNKNOWN;
  return status;
}


/* Function:  esl_msafile_clustal_Read()
 * Synopsis:  Read in a CLUSTAL or CLUSTALLIKE alignment.
 *
 * Purpose:   Read an MSA from an open <ESLX_MSAFILE> <afp>, parsing
 *            for Clustal or Clustal-like format, starting from the 
 *            current point. (<afp->format> is expected to be
 *            <eslMSAFILE_CLUSTAL> or <eslMSAFILE_CLUSTALLIKE>.) Create a
 *            new multiple alignment, and return a ptr to that
 *            alignment in <*ret_msa>.  Caller is responsible for
 *            free'ing this <ESL_MSA>.
 *
 * Args:      afp     - open <ESL_MSAFILE>
 *            ret_msa - RETURN: newly parsed <ESL_MSA>
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEOF> if no (more) alignment data are found in
 *            <afp>, and <afp> is returned at EOF. 
 *
 *            <eslEFORMAT> on a parse error. <*ret_msa> is set to
 *            <NULL>. <afp> contains information sufficient for
 *            constructing useful diagnostic output: 
 *            | <afp->errmsg>       | user-directed error message     |
 *            | <afp->linenumber>   | line # where error was detected |
 *            | <afp->line>         | offending line (not NUL-term)   |
 *            | <afp->n>            | length of offending line        |
 *            | <afp->bf->filename> | name of the file                |
 *            and <afp> is poised at the start of the following line,
 *            so (in principle) the caller could try to resume
 *            parsing.
 *
 * Throws:    <eslEMEM> - an allocation failed.
 *            <eslESYS> - a system call such as fread() failed
 *            <eslEINCONCEIVABLE> - "impossible" corruption 
 */
int
esl_msafile_clustal_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa     = NULL;
  char     *p       = NULL;
  esl_pos_t n       = 0;
  char     *tok     = NULL;
  esl_pos_t ntok    = 0;
  int       nblocks = 0;
  int       idx     = 0;
  int       nseq    = 0;
  int64_t   alen    = 0;
  int64_t   cur_alen;
  esl_pos_t pos;
  esl_pos_t name_start, name_len;
  esl_pos_t seq_start, seq_len;
  esl_pos_t block_seq_start, block_seq_len;
  int       status;

  afp->errmsg[0] = '\0';
  
#ifdef eslAUGMENT_ALPHABET
  if (afp->abc   &&  (msa = esl_msa_CreateDigital(afp->abc, 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
#endif
  if (! afp->abc &&  (msa = esl_msa_Create(                 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }

  /* skip leading blank lines in file */
  while ( (status = eslx_msafile_GetLine(afp, &p, &n)) == eslOK  && esl_memspn(afp->line, afp->n, " \t") == afp->n) ;
  if      (status != eslOK)  goto ERROR; /* includes normal EOF */
    
  /* That first line says something like: "CLUSTAL W (1.83) multiple sequence alignment" */
  if (esl_memtok(&p, &n, " \t", &tok, &ntok) != eslOK)                             ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing CLUSTAL header");
  if (afp->format == eslMSAFILE_CLUSTAL && ! esl_memstrpfx(tok, ntok, "CLUSTAL"))  ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing CLUSTAL header"); 
  if (! esl_memstrcontains(p, n, "multiple sequence alignment"))                   ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing CLUSTAL header");

  /* skip blank lines again */
  do {
    status = eslx_msafile_GetLine(afp, &p, &n);
    if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "no alignment data following header");
    else if (status != eslOK) goto ERROR;
  } while (esl_memspn(afp->line, afp->n, " \t") == afp->n); /* idiom for "blank line" */

  /* Read the file a line at a time. */
  do { 		/* afp->line, afp->n is now the first line of a block... */
    idx = 0;
    do {
      for (pos = 0;     pos < n; pos++) if (! isspace(p[pos])) break;  name_start = pos; 
      for (pos = pos+1; pos < n; pos++) if (  isspace(p[pos])) break;  name_len   = pos - name_start;
      for (pos = pos+1; pos < n; pos++) if (! isspace(p[pos])) break;  seq_start  = pos;      
      if (pos >= n) ESL_XFAIL(eslEFORMAT, afp->errmsg, "invalid alignment line");
      for (pos = n-1; pos > 0; pos--)   if (! isspace(p[pos])) break;  seq_len    = pos - seq_start + 1;

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
	if (idx >= msa->sqalloc &&  (status = esl_msa_Expand(msa))           != eslOK) goto ERROR;
	if ( (status = esl_msa_SetSeqName(msa, idx, p+name_start, name_len)) != eslOK) goto ERROR;
	nseq++;
      } else {
	if (! esl_memstrcmp(p+name_start, name_len, msa->sqname[idx]))
	  ESL_XFAIL(eslEFORMAT, afp->errmsg, "expected sequence %s on this line, but saw %.*s", msa->sqname[idx], (int) name_len, p+name_start);
      }

      /* Append the sequence. */
      cur_alen = alen;
#ifdef eslAUGMENT_ALPHABET
      if (msa->abc)    { status = esl_abc_dsqcat(afp->inmap, &(msa->ax[idx]),   &(cur_alen), p+seq_start, seq_len); }
#endif
      if (! msa->abc)  { status = esl_strmapcat (afp->inmap, &(msa->aseq[idx]), &(cur_alen), p+seq_start, seq_len); }
      if      (status == eslEINVAL)    ESL_XFAIL(eslEFORMAT, afp->errmsg, "one or more invalid sequence characters");
      else if (status != eslOK)        goto ERROR;
      if (cur_alen - alen != seq_len) ESL_XFAIL(eslEFORMAT, afp->errmsg, "unexpected number of seq characters");

      /* get next line. if it's a consensus line, we're done with the block */
      status = eslx_msafile_GetLine(afp, &p, &n);
      if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "alignment block did not end with consensus line");
      else if (status != eslOK)  goto ERROR;

      idx++;
    } while (esl_memspn(afp->line, afp->n, " .:*") < afp->n); /* end loop over a block */
    
    if (idx != nseq) ESL_XFAIL(eslEFORMAT, afp->errmsg, "last block didn't contain same # of seqs as earlier blocks");

    /* skip blank lines until we find start of next block, or EOF */
    do {
      status = eslx_msafile_GetLine(afp, &p, &n);
      if      (status == eslEOF) break;
      else if (status != eslOK)  goto ERROR;
    } while (esl_memspn(p, n, " \t") == n); 
    
    alen += block_seq_len;
    nblocks++;
  } while (status == eslOK);	/* normal end has status == EOF after last block. */

  msa->nseq = nseq;
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
 * Purpose:   Write alignment <msa> to output stream <fp>, in
 *            format <fmt>. If <fmt> is <eslMSAFILE_CLUSTAL>,
 *            write strict CLUSTAL 2.1 format. If <fmt>
 *            is <eslMSAFILE_CLUSTALLIKE>, put "EASEL (VERSION)"
 *            in the header.
 *            
 *            The alignment is written in blocks of 60 aligned
 *            residues at a time.
 *            
 *            Constructing the CLUSTAL consensus line properly
 *            requires knowing the alphabet. If the <msa> is in text
 *            mode, we don't know the alphabet, so then we use a
 *            simplified consensus line, with '*' marking completely
 *            conserved columns, ' ' on everything else. If the <msa>
 *            is in digital mode and of type <eslAMINO>, then we also
 *            use Clustal's "strong" and "weak" residue group
 *            annotations, ':' and '.'.  Strong groups are STA, NEQK,
 *            NHQK, NDEQ, QHRK, MILV, MILF, HY, and FYW. Weak groups
 *            are CSA, ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK,
 *            NEQHRK, FVLIM, and HFY.
 *            
 * Args:      fp  - open output stream, writable
 *            msa - alignment to write      
 *            fmt - eslMSAFILE_CLUSTAL or eslMSAFILE_CLUSTALLIKE      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msafile_clustal_Write(FILE *fp, const ESL_MSA *msa, int fmt)
{
  int       cpl        = 60;
  int       maxnamelen = 0;
  int       namelen;
  char     *consline   = NULL;
  char     *buf        = NULL;
  int64_t   apos;
  int       i;
  int       status;

  ESL_ALLOC(buf, sizeof(char) * (cpl+1));
  buf[cpl] = '\0';
  for (i = 0; i < msa->nseq; i++)
    {
      namelen = strlen(msa->sqname[i]);
      maxnamelen = ESL_MAX(namelen, maxnamelen);
    }

  /* Make a CLUSTAL-like consensus line */
#ifdef eslAUGMENT_ALPHABET 
  if (msa->abc &&  (status = make_digital_consensus_line(msa, &consline)) != eslOK) goto ERROR;
#endif
  if (! msa->abc && (status = make_text_consensus_line(msa, &consline))   != eslOK) goto ERROR;

  /* The magic header */
  if      (fmt == eslMSAFILE_CLUSTAL)     fprintf(fp, "CLUSTAL 2.1 multiple sequence alignment\n");
  else if (fmt == eslMSAFILE_CLUSTALLIKE) fprintf(fp, "EASEL (%s) multiple sequence alignment\n", EASEL_VERSION);

  /* The alignment */
  for (apos = 0; apos < msa->alen; apos += cpl)
    {
      fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
#ifdef eslAUGMENT_ALPHABET 
	  if (msa->abc)   esl_abc_TextizeN(msa->abc, msa->ax[i]+apos+1, cpl, buf);
#endif
	  if (! msa->abc) strncpy(buf, msa->aseq[i]+apos, cpl);
	  fprintf(fp, "%-*s %s\n", maxnamelen, msa->sqname[i], buf);
	}
      strncpy(buf, consline+apos, cpl);
      fprintf(fp, "%-*s %s\n", maxnamelen, "", buf);
    }

  free(buf);
  free(consline);
  return eslOK;

 ERROR:
  if (buf)      free(buf);
  if (consline) free(consline);
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
  int    status;

  ESL_ALLOC(consline, sizeof(char) * (msa->alen+1));
  ESL_ALLOC(ct,       sizeof(int)  * 27);

  for (apos = 0; apos < msa->alen; apos++)
    {
      for (x = 0; x <= 26; x++) ct[x] = 0;

      for (i = 0; i < msa->nseq; i++)
	{
	  x = toupper(msa->aseq[i][apos]) - 'A';
	  if (x >= 0 && x < 26) ct[x]++;
	  else                  ct[26]++;
	}
	  
      consline[apos] = ' ';
      for (x = 0; x < 26; x++) /* not including gaps */
	if (ct[x] == msa->nseq) { consline[apos] = '*'; break; }
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
static char
matches_group_digital(ESL_ALPHABET *abc, double *ct, double nseq, char *residues)
{
  double total = 0.;
  char *c;

  for (c = residues; *c; c++) 
    total += ct[ (int) esl_abc_DigitizeSymbol(abc, *c) ];
  return (total == nseq ? TRUE : FALSE); /* easily changed in the future to be some threshold fraction of nseq */
}

static int
make_digital_consensus_line(const ESL_MSA *msa, char **ret_consline)
{
  char   *consline = NULL;
  double *ct       = NULL;
  int    i, apos, x;
  int    status;

  ESL_ALLOC(consline, sizeof(char)   * (msa->alen+1));
  ESL_ALLOC(ct,       sizeof(double) * (msa->abc->K+1));  

  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (x = 0; x <= msa->abc->K; x++) ct[x] = 0.;

      for (i = 0; i < msa->nseq; i++)
	esl_abc_DCount(msa->abc, ct, msa->ax[i][apos], 1.0);
	
      consline[apos-1] = ' ';
      for (x = 0; x < msa->abc->K; x++)
	if (ct[x] >= (double) msa->nseq) { consline[apos-1] = '*'; break; }

      /* clustalw's "strong groups" */
      if (msa->abc->type == eslAMINO && consline[apos-1] == ' ') {
	if (matches_group_digital(msa->abc, ct, (double) msa->nseq, "STA")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "NEQK") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "NHQK") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "NDEQ") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "QHRK") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "MILV") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "MILF") ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "HY")   ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "FYW"))
	  consline[apos-1] = ':';
      }

      /* clustalw's "weak groups" */
      if (msa->abc->type == eslAMINO && consline[apos-1] == ' ') {
	if (matches_group_digital(msa->abc, ct, (double) msa->nseq, "CSA")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "ATV")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "SAG")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "STNK")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "STPA")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "SGND")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "SNDEQK")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "NDEQHK")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "NEQHRK")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "FVLIM")  ||
	    matches_group_digital(msa->abc, ct, (double) msa->nseq, "HFY"))
	  consline[apos-1] = '.';
      }
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
 * 2. Example.
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

  if ( (status = eslx_msafile_Open(NULL, filename, fmt, NULL, &afp)) != eslOK) 
    eslx_msafile_OpenFailure(afp, status);

  if ( (status = esl_msafile_clustal_Read(afp, &msa))         != eslOK)
    eslx_msafile_ReadFailure(afp, status);

  printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	 1, msa->name, msa->nseq, (int) msa->alen);

  esl_msafile_clustal_Write(stdout, msa, eslMSAFILE_CLUSTAL);
  esl_msa_Destroy(msa);
  eslx_msafile_Close(afp);
  exit(0);
}
/*::cexcerpt::msafile_clustal_example::end::*/
#endif /*eslMSAFILE_CLUSTAL_EXAMPLE*/
/*--------------------- end of example --------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
