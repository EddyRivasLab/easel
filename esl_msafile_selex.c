/* I/O of multiple sequence alignment files in SELEX format
 *
 * Contents:
 *   1. API for reading/writing SELEX input.
 *   2. Internal functions for a block of input lines.
 *   3. Internal functions for parsing SELEX input.
 *   4. License and copyright.
 *   
 * Notes:
 *   In SELEX, a tricky and unusual issue is that spaces are allowed
 *   as gaps, and can even overlap names. Alignments like this are
 *   legitimate:
 *        seq1_longname ACCCGGT
 *        seq2      AAAAACCCGGTT
 *  
 *  You can't determine the aligned length of any sequence in the
 *  block without seeing the whole block.  We define an internal
 *  object (an ESL_SELEX_BLOCK) and some local functions to handle
 *  reading a block of input lines from an input buffer.
 *
 *  Even though spaces are allowed as gaps in input files, Easel
 *  disallows them internally, even in text-mode alignments. Any
 *  spaces are mapped to '.'.
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
#include "esl_msafile_selex.h"

#define eslSELEX_LINE_SQ 1
#define eslSELEX_LINE_RF 2
#define eslSELEX_LINE_CS 3
#define eslSELEX_LINE_SS 4
#define eslSELEX_LINE_SA 5

typedef struct {
  char     **line;		/* line[0..nlines-1][0..llen-1]: memory lines in input buffer */
  esl_pos_t *llen;		/* length of line[] in bytes                                  */
  esl_pos_t *offsets;		/* offset of start of each line in input buffer               */
  int64_t   *linenum;		/* line number of each line[] in input                        */
  int       *ltype;		/* code for line type: eslSELEX_LINE_SQ, etc.                 */
  esl_pos_t *lpos;		/* leftmost position of seq data on line[], 0..llen-1 [or -1] */
  esl_pos_t *rpos;              /* rightmost pos of seq data on line[], 0..llen-1 [or -1]     */
  int        nlines;		/* number of lines in this block                              */
  int        nalloc;		/* number of lines allocated for (>=nlines)                   */
  esl_pos_t  anchor;		/* input buffer anchor set at the start of the block          */
} ESL_SELEX_BLOCK;

static ESL_SELEX_BLOCK *selex_block_Create(int nalloc);
static int              selex_block_Grow(ESL_SELEX_BLOCK *b);
static void             selex_block_Destroy(ESL_SELEX_BLOCK *b);

static int selex_ErrorInBlock(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, int idx);
static int selex_read_block  (ESLX_MSAFILE *afp, ESL_SELEX_BLOCK **block_p);
static int selex_first_block (ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA **ret_msa);
static int selex_other_block (ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA *msa);
static int selex_append_block(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA *msa);

/*****************************************************************
 * 1. API for reading/writing SELEX input
 *****************************************************************/

/* Function:  esl_msafile_selex_SetInmap()
 * Synopsis:  Finishes configuring input map for SELEX format
 *
 * Purpose:   SELEX not only tolerates spaces in input, it
 *            allows a space as a gap character. (Which significantly
 *            complicates parsing.)
 *            
 *            The inmap may not contain <eslDSQ_IGNORED> mappings.
 *            Annotation lines are parsed literally: every character
 *            is copied. If some characters of the aligned sequence
 *            were ignored, we'd be misaligned with the annotation.
 *            In general, because of this, it seems unlikely that any
 *            alignment format would use <eslDSQ_IGNORED> mappings.
 */
int
esl_msafile_selex_SetInmap(ESLX_MSAFILE *afp)
{
#ifdef eslAUGMENT_ALPHABET
  if (  afp->abc) afp->inmap[' '] = esl_abc_XGetGap(afp->abc);
#endif
  if (! afp->abc) afp->inmap[' '] = '.';   /* Easel does not allow spaces as gap characters. */
  return eslOK;
}

/* Function:  esl_msafile_selex_CheckFileFormat()
 * Synopsis:  Checks whether an input source appears to be in SELEX format.
 *
 * Purpose:   Check whether the input source <bf> appears to be a 
 *            SELEX-format alignment file, starting from the current point,
 *            to the end of the input. Return <eslOK> if so, <eslFAIL>
 *            if not.
 *            
 *            This is a SELEX-specific plugin for <esl_msafile_GuessFormat()>.
 */
int
esl_msafile_selex_CheckFileFormat(ESL_BUFFER *bf)
{
  esl_pos_t start_offset = -1;
  int       block_nseq   = 0;	 /* Number of seqs in each block is checked     */
  int       nseq         = 0;
  esl_pos_t block_nres   = 0;	 /* Number of residues in each line is checked  */
  char     *firstname    = NULL; /* First seq name of every block is checked    */
  esl_pos_t namelen      = 0;
  int       blockidx     = 0;
  int       in_block     = FALSE;
  char     *p, *tok;
  esl_pos_t n,  toklen;
  int       status;

  /* Anchor at the start of the input, so we can rewind */
  start_offset = esl_buffer_GetOffset(bf);
  if ( (status = esl_buffer_SetAnchor(bf, start_offset)) != eslOK) goto ERROR;

  while ( (status = esl_buffer_GetLine(bf, &p, &n)) == eslOK)
    {
      /* Some automatic giveaways of SELEX format */
      if (esl_memstrpfx(p, n, "#=RF")) { status = eslOK; goto DONE; }
      if (esl_memstrpfx(p, n, "#=CS")) { status = eslOK; goto DONE; }
      if (esl_memstrpfx(p, n, "#=SS")) { status = eslOK; goto DONE; }
      if (esl_memstrpfx(p, n, "#=SA")) { status = eslOK; goto DONE; }
      
      /* skip comments */
      if (esl_memstrpfx(p, n, "#"))    continue;
      
      /* blank lines: end block, reset block counters */
      if (esl_memspn(p, n, " \t\r\n") == n)
	{
	  if (nseq && block_nseq != nseq) { status = eslFAIL; goto DONE;} /* each block has same # of seqs */
	  if (in_block) blockidx++;
	  if (blockidx >= 3) { status = eslOK; goto DONE; } /* stop after three blocks; we're pretty sure by now */
	  in_block   = FALSE;
	  block_nres = 0;
	  block_nseq = nseq;
	  nseq       = 0;
	  continue;
	}

      /* else we're a "sequence" line. test for two and only two non-whitespace
       * fields; test that second field has same length; test that each block
       * starts with the same seq name..
       */
      in_block = TRUE;
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) != eslOK) goto ERROR; /* there's at least one token - we already checked for blank lines */
      if (nseq == 0)	/* check first seq name in each block */
	{
	  if (blockidx == 0) { firstname = tok; namelen = toklen; } /* First block: set the name we check against. */
	  else if (toklen != namelen || memcmp(tok, firstname, toklen) != 0) { status = eslFAIL; goto DONE; } /* Subsequent blocks */
	}
      if (esl_memtok(&p, &n, " \t", &tok, &toklen) != eslOK) { status = eslFAIL; goto DONE; }
      if (block_nres && toklen != block_nres)                { status = eslFAIL; goto DONE; }
      block_nres = toklen;
      if (esl_memtok(&p, &n, " \t", &tok, &toklen) == eslOK) { status = eslFAIL; goto DONE; }
      nseq++;
    }
  if (status != eslEOF) goto ERROR;  /* EOF is expected and good; anything else is bad */

  if (in_block) blockidx++; 
  if (! blockidx) status = eslFAIL;
  
 DONE:
  if (start_offset != -1) { 
    if (esl_buffer_SetOffset(bf, start_offset)   != eslOK) goto ERROR;
    if (esl_buffer_RaiseAnchor(bf, start_offset) != eslOK) goto ERROR;
    start_offset = -1;
  }
  return status;

 ERROR:
  if (start_offset != -1) { 
    esl_buffer_SetOffset(bf, start_offset);
    esl_buffer_RaiseAnchor(bf, start_offset);
  }
  return status;
}



/* Function:  esl_msafile_selex_Read()
 * Synopsis:  Read in a SELEX format alignment.
 *
 * Purpose:   Read an MSA from an open <ESLX_MSAFILE> <afp>, 
 *            parsing for SELEX format, starting from the
 *            current point. (<afp->format> is expected to
 *            be <eslMSAFILE_SELEX>.) Create a new multiple
 *            alignment and return it via <*ret_msa>. 
 *            Caller is responsible for free'ing this
 *            <ESL_MSA>.
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
esl_msafile_selex_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA         *msa     = NULL;
  ESL_SELEX_BLOCK *b       = NULL;
  int32_t          nblocks = 0;
  int              status;

  afp->errmsg[0] = '\0';

  /* Check the <afp>'s cache first */
  if (afp->msa_cache) return eslx_msafile_Decache(afp, ret_msa);

  while ( (status = selex_read_block(afp, &b)) == eslOK)
    {
      if      (! nblocks &&  (status = selex_first_block(afp, b, &msa)) != eslOK) goto ERROR;
      else if (  nblocks &&  (status = selex_other_block(afp, b, msa))  != eslOK) goto ERROR;

      if ((status = selex_append_block(afp, b, msa)) != eslOK) goto ERROR;

      esl_buffer_RaiseAnchor(afp->bf, b->anchor);
      b->anchor = -1;

      nblocks++;
    }
  /* selex_read_block took care of destroying the block */
  if (status != eslEOF || nblocks == 0) goto ERROR;

  msa->offset = 0;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (b) {
    if (b->anchor != -1) esl_buffer_RaiseAnchor(afp->bf, b->anchor);
    selex_block_Destroy(b);
  }
  *ret_msa = NULL;
  return status;
}

/* Function:  esl_msafile_selex_Write()
 * Synopsis:  Write a SELEX format alignment to a stream
 *
 * Purpose:   Write alignment <msa> to output stream <fp>,
 *            in SELEX format. The alignment is written
 *            in blocks of 60 aligned residues at a time.
 *
 * Args:      fp  - open output stream, writable
 *            msa - alignment to write
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msafile_selex_Write(FILE *fp, const ESL_MSA *msa)
{
  int     cpl        = 60;
  int     maxnamelen = 4;		/* init to 4 because minimum name field is #=CS, etc. */
  int     namelen;
  char   *buf        = NULL;
  int     i;
  int64_t apos;
  int     status;

  ESL_ALLOC(buf, sizeof(char) * (cpl+1));
  buf[cpl] = '\0';
  for (i = 0; i < msa->nseq; i++) {
    namelen    = strlen(msa->sqname[i]);
    maxnamelen = ESL_MAX(namelen, maxnamelen);
  }

  for (apos = 0; apos < msa->alen; apos += cpl)
    {
      if (apos)         fprintf(fp, "\n");
      if (msa->ss_cons) fprintf(fp, "%-*s %.*s\n", maxnamelen, "#=CS", cpl, msa->ss_cons+apos);
      if (msa->rf)      fprintf(fp, "%-*s %.*s\n", maxnamelen, "#=RF", cpl, msa->rf+apos);

      for (i = 0; i < msa->nseq; i++)
	{
#ifdef eslAUGMENT_ALPHABET 
	  if (msa->abc)   esl_abc_TextizeN(msa->abc, msa->ax[i]+apos+1, cpl, buf);
#endif
	  if (! msa->abc) strncpy(buf, msa->aseq[i]+apos, cpl);
	  fprintf(fp, "%-*s %s\n", maxnamelen, msa->sqname[i], buf);	  

	  if (msa->ss && msa->ss[i]) fprintf(fp, "%-*s %.*s\n", maxnamelen, "#=SS", cpl, msa->ss[i]+apos);
	  if (msa->sa && msa->sa[i]) fprintf(fp, "%-*s %.*s\n", maxnamelen, "#=SA", cpl, msa->sa[i]+apos);
	}
    }

  free(buf);
  return eslOK;
  
 ERROR:
  if (buf) free(buf);
  return status;
}
/*--------------------- end, SELEX i/o API ----------------------*/



/*****************************************************************
 * 2. Internal functions handling a block of input lines.
 *****************************************************************/

static ESL_SELEX_BLOCK *
selex_block_Create(int nalloc)
{
  ESL_SELEX_BLOCK *b = NULL;
  int              idx;
  int              status;
  
  ESL_ALLOC(b,       sizeof(ESL_SELEX_BLOCK));
  b->line    = NULL;
  b->llen    = NULL;
  b->offsets = NULL;
  b->linenum = NULL;
  b->ltype   = NULL;
  b->lpos    = NULL;
  b->rpos    = NULL;
  b->nlines  = 0;
  b->anchor  = -1;		/* -1 is a flag for "unused" */

  ESL_ALLOC(b->line,    sizeof(char *)    * nalloc);
  ESL_ALLOC(b->llen,    sizeof(esl_pos_t) * nalloc);
  ESL_ALLOC(b->offsets, sizeof(esl_pos_t) * nalloc);
  ESL_ALLOC(b->linenum, sizeof(int64_t)   * nalloc);
  ESL_ALLOC(b->ltype,   sizeof(int)       * nalloc);
  ESL_ALLOC(b->lpos,    sizeof(esl_pos_t) * nalloc);
  ESL_ALLOC(b->rpos,    sizeof(esl_pos_t) * nalloc);
  for (idx = 0; idx < nalloc; idx++) 
    { 
      b->line[idx]    = NULL; 
      b->llen[idx]    = 0; 
      b->offsets[idx] = 0; 
      b->linenum[idx] = 0; 
      b->ltype[idx]   = 0; 
      b->lpos[idx]    = 0; 
      b->rpos[idx]    = 0; 
    }
  b->nalloc = nalloc;
  return b;

 ERROR:
  if (b) selex_block_Destroy(b);
  return NULL;
}

static int
selex_block_Grow(ESL_SELEX_BLOCK *b)
{
  int idx;
  int status;

  ESL_REALLOC(b->line,    sizeof(char *)    * b->nalloc * 2);
  ESL_REALLOC(b->llen,    sizeof(esl_pos_t) * b->nalloc * 2);
  ESL_REALLOC(b->offsets, sizeof(esl_pos_t) * b->nalloc * 2);
  ESL_REALLOC(b->linenum, sizeof(int64_t)   * b->nalloc * 2);
  ESL_REALLOC(b->ltype,   sizeof(int)       * b->nalloc * 2);
  ESL_REALLOC(b->lpos,    sizeof(esl_pos_t) * b->nalloc * 2);
  ESL_REALLOC(b->rpos,    sizeof(esl_pos_t) * b->nalloc * 2);
  for (idx = b->nalloc; idx < b->nalloc*2; idx++) 
    { 
      b->line[idx]    = NULL; 
      b->llen[idx]    = 0; 
      b->offsets[idx] = 0; 
      b->linenum[idx] = 0; 
      b->ltype[idx]   = 0; 
      b->lpos[idx]    = 0;
      b->rpos[idx]    = 0; 
    }	
  b->nalloc  *= 2;
  return eslOK;

 ERROR:
  return status;
}
  
static void
selex_block_Destroy(ESL_SELEX_BLOCK *b)
{
  if (!b) return;
  if (b->line)    free(b->line);
  if (b->llen)    free(b->llen);
  if (b->offsets) free(b->offsets);
  if (b->linenum) free(b->linenum);
  if (b->ltype)   free(b->ltype);
  if (b->lpos)    free(b->lpos);
  if (b->rpos)    free(b->rpos);
  return;
}

/*------- end, internal functions for input line blocks ---------*/




/*****************************************************************
 * 3. Internal functions for parsing SELEX input.
 *****************************************************************/

/* Before we return a parse error,
 * reset the <afp> so its current line is the one at fault.
 */
static int
selex_ErrorInBlock(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, int which)
{
  afp->line       = b->line[which];
  afp->n      = b->llen[which];
  afp->lineoffset = b->offsets[which];
  afp->linenumber = b->linenum[which];
  return esl_buffer_SetOffset(afp->bf, b->offsets[which] + b->llen[which]);
}

/* selex_read_block:  read one block of alignment data.
 *
 * Note that line numbers aren't necessarily consecutive,
 * because we're stripping out comment lines here. On a parse error
 * on a specific line, we're going to reset the buffer to that line,
 * and we'll need the linenumber to do that reset.
 * 
 * The <afp> detected the end of the block by reading a blank line, or EOF.
 * Thus its point is at the next line after that blank, or at EOF.
 * 
 * The <afp> has a stable anchor set at (*block_p)->anchor.
 * Caller must raise this anchor when it's done parsing the block.
 * 
 * Returns: <eslOK> on success.
 *          
 *          <eslEOF> if no more blocks are found in the input.
 *          <eslEFORMAT> on failure, if a subsequent block has a
 *          different number of data lines than the first block.
 *          On normal errors, all the references are returned set to NULL.
 *       
 * Throws:  <eslEMEM> on allocation failure.      
 */
static int
selex_read_block(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK **block_p) 
{
  ESL_SELEX_BLOCK *b      = *block_p; /* now b==NULL if first block; or on subsequent blocks, reuse prev block storage. */
  int              idx    = 0;
  int              status;

  /* Advance past blank lines until we have the first line of next
   * block.  We may hit a normal EOF here, in which case we return
   * EOF, we're done.
   */
  do { 
    if ( ( status = eslx_msafile_GetLine(afp)) != eslOK) goto ERROR;                               /* EOF here is a normal EOF   */
  } while (esl_memspn(afp->line, afp->n, " \t\r\n") == afp->n ||                                   /* idiomatic for "blank line" */
	   (esl_memstrpfx(afp->line, afp->n, "#") && ! esl_memstrpfx(afp->line, afp->n, "#=")));   /* a SELEX comment line       */

  /* if this is first block, allocate block; subsequent blocks reuse it */
  if (!b && (b = selex_block_Create(16)) == NULL) goto ERROR; 
  
  /* Anchor stably at this point. */
  b->anchor = afp->lineoffset;
  if ((status = esl_buffer_SetStableAnchor(afp->bf, b->anchor)) != eslOK) goto ERROR;

  /* Parse for a block of lines. */
  do {
    if (b->nalloc && idx == b->nalloc && (status = selex_block_Grow(b)) != eslOK) goto ERROR;

    b->line[idx]     = afp->line;
    b->llen[idx]     = afp->n;
    b->offsets[idx]  = afp->lineoffset;
    b->linenum[idx]  = afp->linenumber;   /* ltype, lpos, rpos aren't set yet */
    idx++;
    
    /* Get next non-comment line; this can be next line of block, blank (end of block), or EOF. */
    do { 
      status = eslx_msafile_GetLine(afp); 
    } while ( status == eslOK && (esl_memstrpfx(afp->line, afp->n, "#") && ! esl_memstrpfx(afp->line, afp->n, "#=")));

  } while (status == eslOK && esl_memspn(afp->line, afp->n, " \t\r\n") < afp->n); /* end of block on EOF or blank line */
  
  if (*block_p && b->nlines != idx) 
    ESL_XFAIL(eslEFORMAT, afp->errmsg, "expected %d lines in block, saw %d", b->nlines, idx);

  b->nlines  = idx;
  *block_p   = b;
  return eslOK;	/* EOF status gets turned into OK: we've read a block successfully and hit EOF. Next call will generate the EOF */

 ERROR:
  if (b && b->anchor != -1) esl_buffer_RaiseAnchor(afp->bf, b->anchor);
  if (b) selex_block_Destroy(b);
  *block_p = NULL;
  return status;
}


/* selex_first_block()
 * 
 * 1. Determine and store line types, in b->ltype[0..b->nlines-1].
 * 2. From the number of eslSELEX_LINE_SQ lines, we know nseq.
 * 3. From nseq, we can allocate a new MSA.
 * 4. Parse each line for sequence names, and store them.
 * 5. Determine lpos[] for each line.
 */
static int
selex_first_block(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa = NULL;
  int       nrf, ncs, nss, nsa, nseq;
  int       has_ss, has_sa;
  char     *p, *tok;
  esl_pos_t n,  ntok;
  int       idx, seqi;
  int       status;

  afp->errmsg[0] = '\0';

  nrf = ncs = nss = nsa = nseq = 0;
  has_ss = has_sa = FALSE;
  for (idx = 0; idx < b->nlines; idx++)
    {
      if      (esl_memstrpfx(b->line[idx], b->llen[idx], "#=RF")) { b->ltype[idx] = eslSELEX_LINE_RF; nrf++; }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=CS")) { b->ltype[idx] = eslSELEX_LINE_CS; ncs++; }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=SS")) { b->ltype[idx] = eslSELEX_LINE_SS; nss++; has_ss = TRUE; }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=SA")) { b->ltype[idx] = eslSELEX_LINE_SA; nsa++; has_sa = TRUE; }
      else                                                        { b->ltype[idx] = eslSELEX_LINE_SQ; nseq++; nss = nsa = 0; }

      if (nss && !nseq) { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "#=SS must follow a sequence");   }
      if (nsa && !nseq) { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "#=SA must follow a sequence");   }
      if (nrf > 1)      { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "Too many #=RF lines for block"); }
      if (ncs > 1)      { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "Too many #=CS lines for block"); }
      if (nss > 1)      { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "Too many #=SS lines for seq");   }
      if (nsa > 1)      { selex_ErrorInBlock(afp, b, idx); ESL_XFAIL(eslEFORMAT, afp->errmsg, "Too many #=SA lines for seq");   }
    }

#ifdef eslAUGMENT_ALPHABET
  if ( afp->abc && (msa = esl_msa_CreateDigital(afp->abc, nseq, -1)) == NULL) { status = eslEMEM; goto ERROR; } /* a growable MSA */
#endif
  if (!afp->abc && (msa = esl_msa_Create(                 nseq, -1)) == NULL) { status = eslEMEM; goto ERROR; } 
  if (has_ss) {
    ESL_ALLOC(msa->ss, sizeof(char *) * nseq);
    for (seqi = 0; seqi < nseq; seqi++) msa->ss[seqi] = NULL;
  }
  if (has_sa) {
    ESL_ALLOC(msa->sa, sizeof(char *) * nseq);
    for (seqi = 0; seqi < nseq; seqi++) msa->sa[seqi] = NULL;
  }
  msa->nseq = nseq;
  msa->alen = 0;

  for (seqi = 0, idx = 0; idx < b->nlines; idx++)
    {
      p = b->line[idx];
      n = b->llen[idx];
      if ( esl_memtok(&p, &n, " \t", &tok, &ntok) != eslOK) ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen"); /* because a block by definition consists of non-blank lines */
      if (b->ltype[idx] == eslSELEX_LINE_SQ) /* otherwise, first token is #=XX marking annotation of some sort */
	{
	  if ( esl_msa_SetSeqName(msa, seqi, tok, ntok)   != eslOK) goto ERROR;
	  seqi++;
	}
      b->lpos[idx] = (n ? p-b->line[idx] : -1);  /* set lpos[] to position of first seq or annotation residue */
    }

  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}


/* selex_other_block()
 * We've already parsed the first block.
 * So we know the order of line types, nseq, and sequence names.
 * Validate that a subsequent block has the same.
 */
static int
selex_other_block(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA *msa)
{
  char     *p, *tok;
  esl_pos_t n, ntok;
  int       idx, seqi;

  /* Validate line types */
  for (idx = 0; idx < b->nlines; idx++)
    {
      if      (esl_memstrpfx(b->line[idx], b->llen[idx], "#=RF")) { if (b->ltype[idx] != eslSELEX_LINE_RF) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "#=RF line isn't in expected order in block"); } }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=CS")) { if (b->ltype[idx] != eslSELEX_LINE_CS) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "#=CS line isn't in expected order in block"); } }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=SS")) { if (b->ltype[idx] != eslSELEX_LINE_SS) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "#=SS line isn't in expected order in block"); } }
      else if (esl_memstrpfx(b->line[idx], b->llen[idx], "#=SA")) { if (b->ltype[idx] != eslSELEX_LINE_SA) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "#=SA line isn't in expected order in block"); } }
      else                                                        { if (b->ltype[idx] != eslSELEX_LINE_SQ) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "sequence line isn't in expected order in block"); } }
    }
  
  /* Validate seq names, and set lpos */
  for (seqi = 0, idx = 0; idx < b->nlines; idx++)
    {
      p = b->line[idx];
      n = b->llen[idx];
      if ( esl_memtok(&p, &n, " \t", &tok, &ntok) != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "can't happen"); /* because a block by definition consists of non-blank lines */      
      if (b->ltype[idx] == eslSELEX_LINE_SQ)
	{
	  if (! esl_memstrcmp(tok, ntok, msa->sqname[seqi]))  { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "expected sequence %s at this line of block", msa->sqname[seqi]); }
	  seqi++;
	}
      while (n && isspace(*p)) { p++; n--; }     /* advance past whitespace to first residue */
      b->lpos[idx] = (n ? p-b->line[idx] : -1);  /* set lpos[] to position of first seq or annotation residue */
    }
  return eslOK;
}

static int
selex_append_block(ESLX_MSAFILE *afp, ESL_SELEX_BLOCK *b, ESL_MSA *msa)
{
  char     *p;
  esl_pos_t pos;
  int       idx, seqi;
  esl_pos_t leftmost, rightmost;
  int64_t   nadd;		/* width of this sequence block, in aligned columns added to msa */
  esl_pos_t nleft, ntext, nright;
  int64_t   alen;
  int       status;
  
  /* Determine rpos for each line.  */
  for (idx = 0; idx < b->nlines; idx++)
    {
      p   = b->line[idx];
      pos = b->llen[idx] - 1;
      while (pos>=0 && isspace(p[pos])) pos--;
      b->rpos[idx] = ( (pos < b->lpos[idx]) ? -1 : pos); /* -1: a completely blank seq line is valid */
    }

  /* Determine leftmost and rightmost positions for entire block */
  leftmost  = b->lpos[0];
  rightmost = b->rpos[0];
  for (idx = 1; idx < b->nlines; idx++) {
    leftmost  = (b->lpos[idx] == -1) ? leftmost  : ESL_MIN(leftmost,  b->lpos[idx]);
    rightmost = (b->rpos[idx] == -1) ? rightmost : ESL_MAX(rightmost, b->rpos[idx]);
  }
  if (rightmost == -1) return eslOK; /* super special case: no sequence or annotation data in this block at all! */
  nadd = rightmost - leftmost + 1;

  /* Appends */
  for (seqi = 0, idx = 0; idx < b->nlines; idx++)
    {
      nleft  = ((b->lpos[idx] != -1) ? b->lpos[idx] - leftmost         : nadd); /* watch special case of all whitespace on data line, lpos>rpos */
      ntext  = ((b->lpos[idx] != -1) ? b->rpos[idx] - b->lpos[idx] + 1 : 0);
      nright = ((b->lpos[idx] != -1) ? rightmost    - b->rpos[idx]     : 0);

      if      (b->ltype[idx] == eslSELEX_LINE_SQ)
	{
#if eslAUGMENT_ALPHABET
	  if (msa->abc)
	    {			/* digital sequence append - mapped, preallocated */
	      ESL_REALLOC(msa->ax[seqi],   sizeof(ESL_DSQ) * (msa->alen + nadd + 2)); 
	      for (alen = msa->alen; alen < msa->alen+nleft;  alen++) msa->ax[seqi][alen+1] = esl_abc_XGetGap(msa->abc);

	      status = esl_abc_dsqcat_noalloc(afp->inmap, msa->ax[seqi], &alen, b->line[idx] + b->lpos[idx], ntext);
	      if      (status == eslEINVAL) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "illegal residue(s) in sequence line"); }
	      else if (status != eslOK)     { selex_ErrorInBlock(afp, b, idx); goto ERROR; }
	      if (alen != msa->alen + nleft + ntext) { selex_ErrorInBlock(afp, b, idx); ESL_EXCEPTION(eslEINCONCEIVABLE, afp->errmsg, "unexpected inconsistency appending a sequence"); };

	      for (; alen < msa->alen+nadd;   alen++) msa->ax[seqi][alen+1] = esl_abc_XGetGap(msa->abc);
	      msa->ax[seqi][alen+1] = eslDSQ_SENTINEL;
	    }
#endif
	  if (! msa->abc)
	    {			/* text mode sequence append - mapped, preallocated */
	      ESL_REALLOC(msa->aseq[seqi], sizeof(char)    * (msa->alen + nadd + 1)); 
	      for (alen = msa->alen; alen < msa->alen+nleft; alen++) msa->aseq[seqi][alen] = '.';

	      status = esl_strmapcat_noalloc(afp->inmap, msa->aseq[seqi], &alen, b->line[idx] + b->lpos[idx], ntext);
	      if      (status == eslEINVAL) { selex_ErrorInBlock(afp, b, idx); ESL_FAIL(eslEFORMAT, afp->errmsg, "illegal residue(s) in input line"); }
	      else if (status != eslOK)     { selex_ErrorInBlock(afp, b, idx); goto ERROR; }
	      if (alen != msa->alen + nleft + ntext) { selex_ErrorInBlock(afp, b, idx); ESL_EXCEPTION(eslEINCONCEIVABLE, afp->errmsg, "unexpected inconsistency appending a sequence"); };

	      for  (; alen < msa->alen+nadd;  alen++) msa->aseq[seqi][alen] = '.';
	      msa->aseq[seqi][alen] = '\0';
	    }
	  seqi++;
	}
      else 
	{			/* annotation append: not mapped, characters are copied exactly as they are */
	  if      (b->ltype[idx] == eslSELEX_LINE_RF) { ESL_REALLOC(msa->rf,         sizeof(char) * (msa->alen + nadd + 1)); p = msa->rf;         }
	  else if (b->ltype[idx] == eslSELEX_LINE_CS) { ESL_REALLOC(msa->ss_cons,    sizeof(char) * (msa->alen + nadd + 1)); p = msa->ss_cons;    }
	  else if (b->ltype[idx] == eslSELEX_LINE_SS) { ESL_REALLOC(msa->ss[seqi-1], sizeof(char) * (msa->alen + nadd + 1)); p = msa->ss[seqi-1]; }
	  else if (b->ltype[idx] == eslSELEX_LINE_SA) { ESL_REALLOC(msa->sa[seqi-1], sizeof(char) * (msa->alen + nadd + 1)); p = msa->sa[seqi-1]; }

	  for (alen = msa->alen; alen < msa->alen+nleft; alen++) p[alen] = '.';
	  if (ntext) memcpy(p+msa->alen+nleft, b->line[idx]+b->lpos[idx], sizeof(char)*ntext);
	  for (alen = msa->alen+nleft+ntext; alen < msa->alen+nadd; alen++) p[alen] = '.';
	  p[alen] = '\0';
	}
    }
  msa->alen += nadd;
  return eslOK;

 ERROR:
  return status;
}    
    
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
