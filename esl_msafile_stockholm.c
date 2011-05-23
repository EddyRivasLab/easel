/* I/O of multiple sequence alignment files in Stockholm format
 * 
 * Contents:
 *   1. API for reading/writing Stockholm input.
 *   2. Internal functions for parsing Stockholm.
 *   x. License and copyright.
 */
#include "esl_config.h"

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#include "esl_mem.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile_stockholm.h"

#define eslSTOCKHOLM_LINE_NEWBLOCK 0
#define eslSTOCKHOLM_LINE_SQ 1
#define eslSTOCKHOLM_LINE_GF 2 
#define eslSTOCKHOLM_LINE_GS 3
#define eslSTOCKHOLM_LINE_GC 4
#define eslSTOCKHOLM_LINE_GR 5


typedef struct {
  /* Having to do with the expected order of lines in each Stockholm block: */
  char      *blinetype;		/* blinetype[bi=0..npb-1] = code for linetype on parsed block line [bi]: GC, GR, or seq  */
  int       *bidx;		/* bidx[bi=0.npb-1] = seq index si=0..nseq-1 of seq or GR on parsed block line [bi]; or -1 for GC lines */
  int       npb;		/* number of lines per block. Increments during 1st; constant thereafter */
  int       balloc;		/* number of lines per block currently allocated for. */
  int       bi;			/* index of current line in a block, 0..npb-1  */
  int       nblock;		/* current block number (starting at 0 while in first block) */

  /* Having to do with the growing lengths (and numbers) of sequences and annotations in <msa>: */
  int64_t   *sqlen;		/* current lengths of ax[0..nseq-1] or aseq[0..nseq-1]  */
  int64_t   *sslen;		/* current lengths of ss[0..nseq-1] */
  int64_t   *salen;		/* current lengths of sa[0..nseq-1] */
  int64_t   *pplen;		/* current lengths of pp[0..nseq-1] */
  int64_t   *ogc_len;		/* current lengths of unparsed gc[0..ngc-1]  */
  int64_t  **ogr_len;		/* current lengths of unparsed gr[0..ngr-1][0..nseq-1] */
  int        ngc;		/* number of unparsed GC tag names (synced to msa->ngc) */
  int        ngr;		/* number of unparsed GR tag names (synced to msa->ngr) */
  int        nseq;		/* # of sqnames currently stored, sqname[0..nseq-1]          */
  int        salloc;		/* # of sqnames currently allocated for (synced to msa->sqalloc) */
  int        si;		/* current (or next expected) sequence index   */

  /* Having to do with the expected length of new sequence or annotation on each line: */
  int64_t   alen;		/* alignment length not including current block being parsed */
  int64_t   block_addlen; 	/* residues added by each seq field in curr block            */
  esl_pos_t block_textlen;	/* width of each sequence field in current block             */
} ESL_STOCKHOLM_PARSEDATA;

/*****************************************************************
 *# 1. Reading/writing Stockholm input.
 *****************************************************************/

int
esl_msafile_stockholm_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA *msa      = NULL;
  int      nseq     = 0;
  int64_t  alen     = 0;
  int      nblocks  = 0;
  int      in_block = FALSE;
  int      seqi     = 0;
  int      status;

  afp->errmsg[0] = '\0';

  /* Check the <afp>'s cache first */
  if (afp->msa_cache) return eslx_msafile_Decache(afp, ret_msa);

  /* Allocate a growable MSA. */
#ifdef eslAUGMENT_ALPHABET
  if (afp->abc   &&  (msa = esl_msa_CreateDigital(afp->abc, 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
#endif
  if (! afp->abc &&  (msa = esl_msa_Create(                 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Skip leading blank lines in file. EOF here is a normal EOF return. */
  do { 
    if ( ( status = eslx_msafile_GetLine(afp)) != eslOK) goto ERROR;     /* EOF? return a normal EOF     */
  } while (esl_memspn(afp->line, afp->n, " \t\r\n") == afp->n ||         /* skip blank lines             */
	   (esl_memstrpfx(afp->line, afp->n, "#")                        /* and skip comment lines       */
	    && ! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM")));      /* but stop on Stockholm header */

  /* Check for the magic Stockholm header */
  if (! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM 1."))  ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing Stockholm header");

  nblocks = 0;
  while ( (status = eslx_msafile_GetLine(afp)) == eslOK)
    {
      p = afp->line;
      n = afp->n;
      while (n && ( *p == ' ' || *p == '\t')) { p++; n--; } /* skip leading whitespace */

      if (! n || *p == '\r' || *p == '\n') 
	{ /* a blank line triggers some end-of-block logic */
	  if (in_block) {
	    nblocks++; 
	    if (nseq && nseq != seqi) ESL_XFAIL(eslEFORMAT, afp->errmsg, "number of seqs in this block did not match number in earlier block(s)\n");
	    nseq     = seqi;
	    seqi     = 0; 
	    in_block = FALSE;
	  }
	  continue; 
	}
      if      (esl_memstrpfx(p, n, "//"))   break; /* end-of-record marker */

      if (*p == '#') 
	{
	  if      (esl_memstrpfx(p, n, "#=GF")) status = stockholm_parse_gf(msa, p, n, afp->errmsg);
	  else if (esl_memstrpfx(p, n, "#=GS")) status = stockholm_parse_gs(msa, p, n, afp->errmsg);
	  else if (esl_memstrpfx(p, n, "#=GC")) status = stockholm_parse_gc(msa, p, n, afp->errmsg);
	  else if (esl_memstrpfx(p, n, "#=GR")) status = stockholm_parse_gr(msa, p, n, afp->errmsg);
	  else                                  status = stockholm_parse_comment(msa, p, n, afp->errmsg);
	}
      else status = stockholm_parse_sequence(msa, p, n, afp->errmsg);
    }
  if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing // terminator after MSA");
  else if (status != eslOK)  goto ERROR;
  if (! nblocks)             ESL_XFAIL(eslEFORMAT, afp->errmsg, "no alignment data followed Stockholm header");
  if (seqi != nseq)          ESL_XFAIL(eslEFORMAT, afp->errmsg, "number of seqs in last block did not match number in earlier block(s)\n");

  msa->nseq = nseq;
  msa->alen = alen;
  *ret_msa  = msa;
  return eslOK;

 ERROR:
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}

/*****************************************************************
 * 2. Internal functions for parsing Stockholm.
 *****************************************************************/

/* An important desideratum: the auxiliary parse data is sufficient
 * to validate each line as we see it, so we can immediately report
 * any errors and the line number they occur on. We do not want to 
 * detect errors in some later validation step, after we've lost track
 * of original line numbers of the input.
 */

static ESL_STOCKHOLM_PARSEDATA *
stockholm_parsedata_Create(ESL_MSA *msa)
{
  ESL_STOCKHOLM_PARSEDATA *pd = NULL;
  int z;

  ESL_ALLOC(pd, sizeof(ESL_STOCKHOLM_PARSEDATA));
  pd->blinetype     = NULL;
  pd->bidx          = NULL;
  pd->npb           = 0;
  pd->balloc        = 0;
  pd->bi            = 0;
  pd->nblock        = 0;

  pd->sqlen         = NULL;
  pd->sslen         = NULL;
  pd->salen         = NULL;
  pd->pplen         = NULL;
  pd->ogc_len       = NULL;
  pd->ogr_len       = NULL;
  pd->ngc           = 0;
  pd->ngr           = 0;
  pd->nseq          = 0;
  pd->salloc        = 0;
  pd->si            = 0;

  pd->alen          = 0;
  pd->block_addlen  = 0;
  pd->block_textlen = 0;

  ESL_ALLOC(pd->blinetype, sizeof(char) * 16);
  ESL_ALLOC(pd->bidx,      sizeof(int)  * 16);
  pd->balloc = 16;

  ESL_ALLOC(pd->sqlen,     sizeof(int64_t) * msa->sqalloc);
  for (z = 0; z < msa->sqalloc; z++) pd->sqlen[z] = 0;
  pd->salloc = msa->sqalloc;

  return pd;

 ERROR:
  stockholm_parsedata_Destroy(pd);
  return NULL;
}

static int
stockholm_parsedata_ExpandSeq(ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa)
{
  int tagidx;
  int z;
  int status;

  ESL_REALLOC(pd->sqlen, sizeof(int64_t) * msa->sqalloc);  
  for (z = pd->salloc; z < msa->sqalloc; z++) pd->sqlen[z] = 0; 

  if (pd->sslen) {
    ESL_REALLOC(pd->sslen,   sizeof(int64_t) * msa->sqalloc);
    for (z = pd->salloc; z < msa->sqalloc; z++) pd->sslen[z] = 0;
  }

  if (pd->salen) {
    ESL_REALLOC(pd->salen,   sizeof(int64_t) * msa->sqalloc);
    for (z = pd->salloc; z < msa->sqalloc; z++) pd->salen[z] = 0;
  }

  if (pd->pplen) {
    ESL_REALLOC(pd->pplen,   sizeof(int64_t) * msa->sqalloc);
    for (z = pd->salloc; z < msa->sqalloc; z++) pd->pplen[z] = 0;
  }

  if (pd->ogr_len) {
    for (tagidx = 0; tagidx < pd->ngr; tagidx++) 
      if (pd->ogr_len[tagidx]) {
	ESL_REALLOC(pd->ogr_len[idx], sizeof(int64_t) * msa->sqalloc);
	for (z = pd->alloc; z < msa->sqalloc; z++) pd->ogr_len[tagidx][z] = 0;
      }
  }

  pd->salloc = msa->sqalloc;
  return eslOK;

 ERROR:
  return status;
}
  

static int
stockholm_parsedata_ExpandBlock(ESL_STOCKHOLM_PARSEDATA *pd)
{
  int status;

  ESL_REALLOC(pd->blinetype, sizeof(char) * (pd->balloc * 2));
  ESL_REALLOC(pd->bidx,      sizeof(int)  * (pd->balloc * 2));
  pd->balloc *= 2;
  return eslOK;

 ERROR:
  return status;
}


static void
stockholm_parsedata_Destroy(ESL_STOCKHOLM_PARSEDATA *pd)
{
  int i;
  if (! pd) return;

  if (pd->blinetype) free(pd->blinetype);
  if (pd->bidx)      free(pd->bidx);

  if (pd->sqlen)     free(pd->sqlen);
  if (pd->sslen)     free(pd->sslen);
  if (pd->salen)     free(pd->salen);
  if (pd->pplen)     free(pd->pplen);

  if (pd->ogc_len)   free(pd->ogc_len);
  if (pd->ogr_len) {
    for (i = 0; i < pd->ngr; i++)
      if (pd->ogr_len[i]) free(pd->ogr_len[i]);
    free(pd->ogr_len);
  }
  
  return;
}

/* stockholm_get_seqidx()
 * 
 * Find the index of a given sequence <name>,<n> in a growing <msa>
 * with associated parse data <pdat>. If caller has a good guess (for
 * instance, the sequences are coming in a previously seen order in a
 * block of seqs or annotation), the caller can pass this information
 * in <guess>, or -1 if it has no guess.
 * 
 * If the name does not already exist in the MSA, then it
 * is assumed to be a new sequence name that we need to store.
 * seqidx is set to pdat->nseq, the MSA is Expand()'ed if necessary
 * to make room, the name is stored in msa->sqname[pdat->nseq],
 * (and in the hash table, if we're keyhash augmented)
 * and pdat->nseq is incremented.
 *
 * Returns:  <eslOK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws: <eslEMEM> on allocation failure
 *         <eslEINVAL> if we try to add a name to a non-growable MSA.
 *         <eslEINCONCEIVABLE> on internal coding errors
 */
static int
stockholm_get_seqidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *name, esl_pos_t n, char linetype, int *ret_idx)
{
  int seqidx;
  int status;
  
  if (pd->si > 0 && linetype == eslSTOCKHOLM_LINE_GR &&
      esl_memstrcmp(msa->sqname[pd->si-1], name, n))
    { *ret_idx = pd->si-1; return eslOK; }
  if (pd->si < pd->nseq &&
      esl_memstrcmp(msa->sqname[pd->si], name, n)) 
    { *ret_idx = pd->si;   return eslOK; }

  /* if we're keyhash-augmented, try to find it in the hash.
   * if we're not, try to find it the hard way.
   */
#ifdef eslAUGMENT_KEYHASH
  status = esl_keyhash_Store(msa->index, name, n, &seqidx);
  if (status == eslEDUP) { *ret_idx = seqidx; return eslOK; }
  if (status != eslOK)   goto ERROR;
#else
  for (seqidx = 0; seqidx < pd->nseq; seqidx++)
    if (esl_memstrcmp(name, n, msa->sqname[seqidx])) 
      { *ret_idx = seqidx; return eslOK; }
#endif
  
  /* if we get here, this is a new name we're adding */
  if (seqidx >= msa->sqalloc) {
    if ( (status = esl_msa_Expand(msa))                    != eslOK) goto ERROR;
    if ( (status = stockholm_parsedata_ExpandSeq(pd, msa)) != eslOK) goto ERROR;
  }  

  if ( (status = esl_msa_SetSeqName(msa, seqidx, name, n)) != eslOK) goto ERROR;
  pd->nseq++;

  *ret_idx = seqidx;
  return eslOK;

 ERROR:
  *ret_idx = -1;
  return status;
}

/* stockholm_parse_gf()
 * Line format is:
 *   #=GF <tag> <text>
 */
static int
stockholm_parse_gf(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gf,  *tag,   *tok;
  esl_pos_t gflen, taglen, toklen;
  int       status;

  if ( (status = esl_memtok(&p, &n, " \t", &gr,  &grlen))  != eslOK) ESL_XFAIL();
  if ( (status = esl_memtok(&p, &n, " \t", &tag, &taglen)) != eslOK) ESL_XFAIL();

  if      (esl_memstrcmp(tag, taglen, "ID")) status = esl_msa_SetName     (msa, p, n);
  else if (esl_memstrcmp(tag, taglen, "AC")) status = esl_msa_SetAccession(msa, p, n);
  else if (esl_memstrcmp(tag, taglen, "DE")) status = esl_msa_SetDesc     (msa, p, n);
  else if (esl_memstrcmp(tag, taglen, "AU")) status = esl_msa_SetAuthor   (msa, p, n);
  else if (esl_memstrcmp(tag, taglen, "GA"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_GA1] = atof[tok];
	msa->cutset[eslMSA_GA1] = TRUE;
      } else ESL_XFAIL();
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_GA2] = atof[tok];
	msa->cutset[eslMSA_GA2] = TRUE;
      } 
      status = eslOK;
    }
  else if (esl_memstrcmp(tag, taglen, "NC"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_NC1] = atof[tok];
	msa->cutset[eslMSA_NC1] = TRUE;
      } else ESL_XFAIL();
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_NC2] = atof[tok];
	msa->cutset[eslMSA_NC2] = TRUE;
      } 
      status = eslOK;
    }
  else if (esl_memstrcmp(tag, taglen, "TC"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_TC1] = atof[tok];
	msa->cutset[eslMSA_TC1] = TRUE;
      } else ESL_XFAIL();
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	msa->cutoff[eslMSA_TC2] = atof[tok];
	msa->cutset[eslMSA_TC2] = TRUE;
      } 
      status = eslOK;
    }
  else 
    status = esl_msa_AddGF(msa, tag, taglen, p, n);

  return status;
}

/* stockholm_parse_sequence():
 * Format of line is:
 *   <seqname>  <aligned text>
 */
static int
stockholm_parse_sequence(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char     *seqname, *text;
  esl_pos_t seqnamelen, textlen;
  int       seqidx = pd->si;
  int64_t   seqlen;
  int       status;
  
  if (esl_memtok(&p, &n, " \t", &seqname, &seqnamelen) != eslOK) ESL_XFAIL();
  if (esl_memtok(&p, &n, " \t", &text,    &textlen)    != eslOK) ESL_XFAIL();

  /* Which seqidx is this?
   * In first block:
   *    1. If #=GS lines set sqname[] completely, then it's pd->si.
   *    2. If #=GS lines set sqname[] partially or out of order, then name is in the keyhash.
   *    3. If we haven't seen name before, then we'll add it: seqidx = pd->nseq, add name to keyhash, possibly reallocate.
   * In subsequent blocks, use recorded indices and linetypes:
   *    4. seqidx = saved bidx[]; should be expecting a SQ line; name should match expected name.
   */
  if (! pd->nblock) /* First block: we're setting npb, bidx[], and blinetype[] as we see them */
    {
      if (pd->si < pd->nseq && ! esl_memstrcmp(msa->sqname[seqidx], name, n)) {
	status = stockholm_get_seqidx(msa, pd, seqname, seqnamelen, &seqidx);
      }
      
      if (pd->bi == pd->balloc &&
	  (status = stockholm_parsedata_ExpandBlock(pd)) != eslOK) goto ERROR;
      pd->bidx[pd->bi]      = seqidx;
      pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_SQ;
      pd->npb++;
    }
  else 
    {				/* subsequent block(s) */
      if (pd->bi >= pd->npb)                               ESL_XFAIL();
      seqidx = pd->bidx[pd->bi];
      if (  pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_SQ) ESL_XFAIL();
      if (! esl_memstrcmp(msa->sqname[seqidx], name, n))   ESL_XFAIL();
    }

#ifdef eslAUGMENT_ALPHABET 
  if (  afp->abc) status = esl_abc_dsqcat(afp->inmap, &(msa->ax[seqidx]),   &(pd->sqlen[seqidx]), text, textlen);
#endif
  if (! afp->abc) status = esl_strmapcat (afp->inmap, &(msa->aseq[seqidx]), &(pd->sqlen[seqidx]), text, textlen);

  if (pd->bi)			/* in a block already? */
    {
      if (textlen           != pd->block_textlen)           ESL_XFAIL();
      if (pd->sqlen[seqidx] != pd->block_addlen + pd->alen) ESL_XFAIL();
    }
  else				/* first line of this block */
    {
      pd->block_textlen = textlen;
      pd->block_addlen  = pd->sqlen[seqidx] - pd->alen;
    }
  
  pd->bi++;
  pd->si = seqidx+1;
  return eslOK;
}

/* A GR line is
 *   #=GR <seqname> <featurename> <text>
 * recognized featurenames: { SS | SA | PP }
 * 
 */
static int
parse_gr(ESL_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gr,   *name,    *tag,   *text;
  esl_pos_t  grlen, namelen,  taglen, textlen;
  int        seqidx, tagidx;
  int        z;

  if (esl_memtok(&p, &n, " \t", &gr,   &grlen)    != eslOK) ESL_XFAIL();
  if (esl_memtok(&p, &n, " \t", &name, &namelen)  != eslOK) ESL_XFAIL();  
  if (esl_memtok(&p, &n, " \t", &tag,  &taglen)   != eslOK) ESL_XFAIL();  
  if (esl_memtok(&p, &n, " \t", &text, &textlen)  != eslOK) ESL_XFAIL();  

  /* Which seqidx is this? likely to be either pd->si-1 (#=GR following a seq) or 
   * pd->si (#=GR preceding a seq) 
   */
  if (! pd->nblock) /* First block: we're setting npb, bidx[], and blinetype[] as we see them */
    {
      if      (pd->si >= 1       && esl_memstrcmp(msa->sqname[pd->si-1], name, namelen)) seqidx = pd->si-1;
      else if (pd->si < pd->nseq && esl_memstrcmp(msa->sqname[pd->si],   name, namelen)) seqidx = pd->si;
      else  status = stockholm_get_seqidx(msa, pd, name, namelen, &seqidx);
      
      if (pd->bi == pd->balloc && (status = stockholm_parsedata_ExpandBlock(pd)) != eslOK) goto ERROR;
      pd->bidx[pd->bi]      = seqidx;
      pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GR;
      pd->npb++;
    }
  else 
    {				/* subsequent block(s) */
      if (pd->bi >= pd->npb)                                   ESL_XFAIL();
      seqidx = pd->bidx[pd->bi];
      if (  pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GR)     ESL_XFAIL();
      if (! esl_memstrcmp(msa->sqname[seqidx], name, namelen)) ESL_XFAIL();
    }

  /* Append the annotation where it belongs
   */
  if (esl_memstrcmp("SS", tag, taglen)) 
    {
      if (! msa->ss) {
	ESL_ALLOC(msa->ss,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->sslen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->ss[z] = NULL; pd->sslen[z] = 0; }
      }
      status = esl_strmapcat(afp->inmap, &(msa->ss[seqidx]), &(pd->sslen[seqidx]), text, textlen);
    }
  else if (esl_memstrcmp("PP", tag, taglen)) 
    {
      if (! msa->pp) {
	ESL_ALLOC(msa->pp,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->pplen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->pp[z] = NULL; pd->pplen[z] = 0; }
      }
      status = esl_strmapcat(afp->inmap, &(msa->pp[seqidx]), &(pd->pplen[seqidx]), text, textlen);
    }
  else if (esl_memstrcmp("SA", tag, taglen)) 
    {
      if (! msa->sa) {
	ESL_ALLOC(msa->sa,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->salen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->sa[z] = NULL; pd->salen[z] = 0; }
      }
      status = esl_strmapcat(afp->inmap, &(msa->sa[seqidx]), &(pd->salen[seqidx]), text, textlen);
    }
  else
    {
      status = esl_msa_StoreGRTag(msa, tag, taglen, &tagidx);
      if (status == eslOK) {
	ESL_REALLOC(pd->ogr_len, sizeof(int64_t *) * msa->ngr);
	for (z = pd->ngr; z < msa->ngr; z++) pd->ogr_len = 0;
      } else if (status != eslEDUP) return status;

      status = esl_strmapcat(afp->inmap, &(msa->gr[tagidx]), &(pd->ogr_len[tagidx]), text, textlen);
    }

  if (pd->bi)			/* in a block already? */
    {
      if (textlen           != pd->block_textlen)           ESL_XFAIL();
      if (pd->sqlen[seqidx] != pd->block_addlen + pd->alen) ESL_XFAIL();
    }
  else				/* first line of this block */
    {
      pd->block_textlen = textlen;
      pd->block_addlen  = pd->sqlen[seqidx] - pd->alen;
    }
  
  pd->bi++;
  return eslOK;
}
  
  

}

static int
parse_comment(ESL_MSA *msa, char *p, esl_pos_t n)
{
  if (n && *p == '#')      { p++; n--; }
  while (n && isspace(*p)) { p++; n--; }

  return esl_msa_AddComment(msa, p, n);
}
  


/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      
 */
static int
esl_msa_StoreGRTag(ESL_MSA *msa, char *tag, esl_pos_t taglen, int *ret_tagidx)
{
  int fstatus;
  int tagidx;

  /* Find the tag, if we have it; else, add it, at tagidx = msa->ngr */
#ifdef eslAUGMENT_KEYHASH
  fstatus = esl_keyhash_Store(msa->gr_idx, tag, taglen, &tagidx); 
  if (fstatus != eslOK && fstatus != eslEDUP) return fstatus; /* eslEMEM */
#else
  for (tagidx = 0; tagidx < msa->ngr; tagidx++)
    if (esl_memstrcmp(tag, taglen, msa->gr_tag[tagidx])) break;
  fstatus = (tagidx == msa->ngr ? eslOK : eslEDUP);
#endif

  /* tagidx is now the index of the new tag. Allocate for it, if needed */
  if (tagidx == msa->ngr)
    {
      ESL_REALLOC(msa->gr_tag,       sizeof(char *)  * (msa->ngr+1)); /* +1, we're allocated one new tag at a time, as needed */
      ESL_REALLOC(msa->gr,           sizeof(char **) * (msa->ngr+1));
      ESL_REALLOC(msa->gr[msa->ngr], sizeof(char *)  *  msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->gr[msa->ngr][i] = NULL;
      if ( (status = esl_memstrdup(tag, taglen, &(msa->gr_tag[msa->ngr]))) != eslOK) return status; /* eslEMEM */
    }
  return fstatus; /* eslOK | eslEDUP */
}

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
