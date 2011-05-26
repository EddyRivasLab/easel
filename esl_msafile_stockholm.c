/* I/O of multiple sequence alignment files in Stockholm format
 * 
 * Contents:
 *   1. API for reading/writing Stockholm input.
 *   2. Internal: ESL_STOCKHOLM_PARSEDATA auxiliary structure.
 *   3. Internal: parsing Stockholm line types.
 *   4. Internal: looking up seq, tag indices.
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

/* Valid line types in an alignment block */
#define eslSTOCKHOLM_LINE_SQ        1
#define eslSTOCKHOLM_LINE_GC_SSCONS 2
#define eslSTOCKHOLM_LINE_GC_SACONS 3
#define eslSTOCKHOLM_LINE_GC_PPCONS 4
#define eslSTOCKHOLM_LINE_GC_RF     5
#define eslSTOCKHOLM_LINE_GC_OTHER  6
#define eslSTOCKHOLM_LINE_GR_SS     7
#define eslSTOCKHOLM_LINE_GR_SA     8
#define eslSTOCKHOLM_LINE_GR_PP     9
#define eslSTOCKHOLM_LINE_GR_OTHER  10

typedef struct {
  /* information about the size of the growing alignment parse */
  int       nseq;		/* # of sqnames currently stored, sqname[0..nseq-1]. Becomes msa->nseq when done */
  int64_t   alen;		/* alignment length not including current block being parsed. Becomes msa->alen when done */

  /* Having to do with the expected order of lines in each Stockholm block: */
  int       in_block;		/* TRUE if we're in a block (GC, GR, or sequence lines) */
  char     *blinetype;		/* blinetype[bi=0..npb-1] = code for linetype on parsed block line [bi]: GC, GR, or seq  */
  int      *bidx;		/* bidx[bi=0.npb-1] = seq index si=0..nseq-1 of seq or GR on parsed block line [bi]; or -1 for GC lines */
  int       npb;		/* number of lines per block. Set by bi in 1st block; checked against bi thereafter */
  int       bi;			/* index of current line in a block, 0..npb-1  */
  int       si;		        /* current (next expected) sequence index, 0..nseq */
  int       balloc;		/* number of lines per block currently allocated for. */

  /* Other information kept per block */
  int       nblock;		/* current block number (starting at 0 while in first block) */
  int       nseq_b;             /* number of sequences seen in this block so far */
  int64_t   alen_b;     	/* residues added by each seq field in curr block            */

  /* Having to do with the growing lengths (and numbers) of sequences and annotations in <msa>: */
  /* yes, needed: used to catch dup lines in a block, such as seq1 xxx, seq1 xxx.               */
  int64_t    ssconslen;		/* current length of #=GC SS_cons annotation */
  int64_t    saconslen;		/* current length of #=GC SA_cons annotation */
  int64_t    ppconslen;		/* current length of #=GC PP_cons annotation */
  int64_t    rflen;		/* current length of #=GC RF annotation */
  int64_t   *sqlen;		/* current lengths of ax[0..nseq-1] or aseq[0..nseq-1]  */
  int64_t   *sslen;		/* current lengths of ss[0..nseq-1] */
  int64_t   *salen;		/* current lengths of sa[0..nseq-1] */
  int64_t   *pplen;		/* current lengths of pp[0..nseq-1] */
  int64_t   *ogc_len;		/* current lengths of unparsed gc[0..ngc-1]  */
  int64_t  **ogr_len;		/* current lengths of unparsed gr[0..ngr-1][0..nseq-1] */
  int        salloc;		/* # of sqnames currently allocated for (synced to msa->sqalloc) */
} ESL_STOCKHOLM_PARSEDATA;

static ESL_STOCKHOLM_PARSEDATA *stockholm_parsedata_Create(ESL_MSA *msa);
static int                      stockholm_parsedata_ExpandSeq  (ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa);
static int                      stockholm_parsedata_ExpandBlock(ESL_STOCKHOLM_PARSEDATA *pd);
static void                     stockholm_parsedata_Destroy    (ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa);

static int stockholm_parse_gf(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n);
static int stockholm_parse_gs(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n);
static int stockholm_parse_gc(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n);
static int stockholm_parse_gr(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n);
static int stockholm_parse_sq(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n);
static int stockholm_parse_comment(ESL_MSA *msa, char *p, esl_pos_t n);

static int stockholm_get_seqidx   (ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *name, esl_pos_t n,      int *ret_idx);
static int stockholm_get_gr_tagidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *tag,  esl_pos_t taglen, int *ret_tagidx);
static int stockholm_get_gc_tagidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *tag,  esl_pos_t taglen, int *ret_tagidx);




/*****************************************************************
 *# 1. API for reading/writing Stockholm input.
 *****************************************************************/

/* Function:  esl_msafile_stockholm_Read()
 * Synopsis:  Read an alignment in Stockholm format.
 *
 * Purpose:   Read an MSA from open <ESLX_MSAFILE> <afp>, 
 *            parsing for Stockholm format. Create a new
 *            MSA, and return it by reference through 
 *            <*ret_msa>. Caller is responsible for freeing
 *            this <ESL_MSA>.
 *            
 * Args:      <afp>     - open <ESL_MSAFILE> to read from
 *            <ret_msa> - RETURN: newly parsed, created <ESL_MSA>
 *
 * Returns:   <eslOK> on success. <*ret_msa> contains the newly
 *            allocated MSA. <afp> is poised at start of next
 *            alignment record, or is at EOF.
 *
 *            <eslEFORMAT> on a parsing error. <*ret_msa> is returned
 *            <NULL>. <afp> contains information sufficient for
 *            constructing diagnostic output: <afp->errmsg> contains a
 *            user-directed error message; <afp->linenumber> contains
 *            the line number; <afp->bf->filename> contains the name
 *            of the file (if it was a file that was being parsed);
 *            <afp->line> is the offending line where the parse error
 *            was detected, and <afp->n> is its length (in general,
 *            <afp->line> is not NUL-terminated). The state of the
 *            <afp>'s buffer is poised at the start of the next line
 *            following the offending one, so that (in principle) the
 *            caller could try to resume parsing.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslESYS> if a system call fails, such as fread().
 *            <*ret_msa> is returned <NULL>.
 */
int
esl_msafile_stockholm_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA                 *msa      = NULL;
  ESL_STOCKHOLM_PARSEDATA *pd       = NULL;
  char                    *p;
  esl_pos_t                n;
  int                      status;

  afp->errmsg[0] = '\0';

  /* Allocate a growable MSA, and auxiliary parse data coupled to the MSA allocation */
#ifdef eslAUGMENT_ALPHABET
  if (afp->abc   &&  (msa = esl_msa_CreateDigital(afp->abc, 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
#endif
  if (! afp->abc &&  (msa = esl_msa_Create(                 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
  if ( (pd = stockholm_parsedata_Create(msa))                        == NULL) { status = eslEMEM; goto ERROR; }

  /* Skip leading blank lines in file. EOF here is a normal EOF return. */
  do { 
    if ( ( status = eslx_msafile_GetLine(afp)) != eslOK) goto ERROR;     /* eslEOF is OK here - end of input (eslEOF) [eslEMEM|eslESYS] */
  } while (esl_memspn(afp->line, afp->n, " \t\r\n") == afp->n ||         /* skip blank lines             */
	   (esl_memstrpfx(afp->line, afp->n, "#")                        /* and skip comment lines       */
	    && ! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM")));      /* but stop on Stockholm header */

  /* Check for the magic Stockholm header */
  if (! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM 1."))  ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing Stockholm header");

  while ( (status = eslx_msafile_GetLine(afp)) == eslOK) /* (eslEOF) [eslEMEM|eslESYS] */
    {
      p = afp->line;
      n = afp->n;
      while (n && ( *p == ' ' || *p == '\t')) { p++; n--; } /* skip leading whitespace */

      if (!n || esl_memstrpfx(p, n, "//"))
	{ /* blank lines and the Stockholm end-of-record // trigger end-of-block logic */
	  if (pd->in_block) {
	    if (pd->nblock) { if (pd->nseq_b != pd->nseq) ESL_XFAIL(eslEFORMAT, afp->errmsg, "number of seqs in block did not match number in earlier block(s)");     }
	    else            { if (pd->nseq_b < pd->nseq)  ESL_XFAIL(eslEFORMAT, afp->errmsg, "number of seqs in block did not match number annotated by #=GS lines"); };
	    if (pd->nblock) { if (pd->bi != pd->npb)      ESL_XFAIL(eslEFORMAT, afp->errmsg, "unexpected number of lines in alignment block"); }

	    pd->nseq     = pd->nseq_b;
	    pd->alen    += pd->alen_b;
	    pd->in_block = FALSE;
	    pd->npb      = pd->bi;
	    pd->bi       = 0;
	    pd->si       = 0;
	    pd->nblock  += 1;
	    pd->nseq_b   = 0;
	    pd->alen_b   = 0;
	  }
	  if   (esl_memstrpfx(p, n, "//"))   break; /* Stockholm end-of-record marker */
	  else continue;			    /* else, on to next block */
	}

      if (*p == '#') 
	{
	  if      (esl_memstrpfx(p, n, "#=GF")) { if ((status = stockholm_parse_gf     (afp, pd, msa, p, n)) != eslOK) goto ERROR; }
	  else if (esl_memstrpfx(p, n, "#=GS")) { if ((status = stockholm_parse_gs     (afp, pd, msa, p, n)) != eslOK) goto ERROR; }
	  else if (esl_memstrpfx(p, n, "#=GC")) { if ((status = stockholm_parse_gc     (afp, pd, msa, p, n)) != eslOK) goto ERROR; }
	  else if (esl_memstrpfx(p, n, "#=GR")) { if ((status = stockholm_parse_gr     (afp, pd, msa, p, n)) != eslOK) goto ERROR; }
	  else                                  { if ((status = stockholm_parse_comment(         msa, p, n)) != eslOK) goto ERROR; }
	}
      else if (                                       (status = stockholm_parse_sq     (afp, pd, msa, p, n)) != eslOK) goto ERROR;
    }
  if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing // terminator after MSA");
  else if (status != eslOK)  goto ERROR;
  if (pd->nblock == 0)       ESL_XFAIL(eslEFORMAT, afp->errmsg, "no alignment data followed Stockholm header");

  msa->nseq = pd->nseq;
  msa->alen = pd->alen;
  stockholm_parsedata_Destroy(pd, msa);
  *ret_msa  = msa;
  return eslOK;

 ERROR:
  if (pd)  stockholm_parsedata_Destroy(pd, msa);
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}


/* Function:  esl_msafile_stockholm_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open Stockholm MSA file.
 *
 * Purpose:   Guess the alpbabet of the sequences in open
 *            Stockholm-format MSA file <afp>.
 *            
 *            On a normal return, <*ret_type> is set to <eslDNA>,
 *            <eslRNA>, or <eslAMINO>, and <afp> is reset to its
 *            original position.
 *
 * Args:      afp      - open Stockholm-format MSA file
 *            ret_type - RETURN: <eslDNA>, <eslRNA>, or <eslAMINO>       
 *
 * Returns:   <eslOK> on success.
 *            <eslENOALPHABET> if alphabet type can't be determined.
 *            In either case, <afp> is rewound to the position it
 *            started at.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_msafile_stockholm_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type)
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
  esl_pos_t n,  toklen;
  int       status;

  for (x = 0; x < 26; x++) ct[x] = 0;

  anchor = esl_buffer_GetOffset(afp->bf);
  if ((status = esl_buffer_SetAnchor(afp->bf, anchor)) != eslOK) { status = eslEINCONCEIVABLE; goto ERROR; } /* [eslINVAL] can't happen here */

  while ( (status = eslx_msafile_GetLine(afp)) == eslOK)
    {
      p = afp->line;
      n = afp->n;
      if ((status = esl_memtok(&p, &n, " \t", &tok, &toklen)) != eslOK || *tok == '#') continue; /* blank lines, annotation, comments */
      /* p now points to the rest of the sequence line */
      
      /* count characters into ct[] array */
      while (n)
	{
	  x = toupper(*p) - 'A';
	  if (x >= 0 && x < 26) { ct[x]++; nres++; }
	  p++;
	  n--;
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
/*--------------- end, api for stockholm i/o --------------------*/




/*****************************************************************
 * 2. Internal: ESL_STOCKHOLM_PARSEDATA auxiliary structure 
 *****************************************************************/

/* The auxiliary parse data is sufficient to validate each line as we
 * see it. Our design requires that we immediately report any errors
 * and the line number they occur on. We do not want to detect errors
 * in some later validation step, after we've lost track of original
 * line numbers of the input. 
 */

static ESL_STOCKHOLM_PARSEDATA *
stockholm_parsedata_Create(ESL_MSA *msa)
{
  ESL_STOCKHOLM_PARSEDATA *pd = NULL;
  int z;
  int status;

  ESL_ALLOC(pd, sizeof(ESL_STOCKHOLM_PARSEDATA));
  pd->nseq          = 0;
  pd->alen          = 0;

  pd->in_block      = FALSE;
  pd->blinetype     = NULL;
  pd->bidx          = NULL;
  pd->npb           = 0;
  pd->bi            = 0;
  pd->si            = 0;
  pd->balloc        = 0;

  pd->nblock        = 0;
  pd->nseq_b        = 0;
  pd->alen_b        = 0;

  pd->ssconslen     = 0;
  pd->saconslen     = 0;
  pd->ppconslen     = 0;
  pd->rflen         = 0;
  pd->sqlen         = NULL;
  pd->sslen         = NULL;
  pd->salen         = NULL;
  pd->pplen         = NULL;
  pd->ogc_len       = NULL;
  pd->ogr_len       = NULL;
  pd->salloc        = 0;

  ESL_ALLOC(pd->blinetype, sizeof(char) * 16);
  ESL_ALLOC(pd->bidx,      sizeof(int)  * 16);
  pd->balloc = 16;

  ESL_ALLOC(pd->sqlen,     sizeof(int64_t) * msa->sqalloc);
  for (z = 0; z < msa->sqalloc; z++) 
    pd->sqlen[z] = 0;
  pd->salloc = msa->sqalloc;
  return pd;

 ERROR:
  stockholm_parsedata_Destroy(pd, msa);
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

  /* don't need to reallocate ogc_len here: it's [0..ngc-1], not by seq */

  if (pd->ogr_len) {
    for (tagidx = 0; tagidx < msa->ngr; tagidx++) 
      if (pd->ogr_len[tagidx]) {
	ESL_REALLOC(pd->ogr_len[tagidx], sizeof(int64_t) * msa->sqalloc);
	for (z = pd->salloc; z < msa->sqalloc; z++) pd->ogr_len[tagidx][z] = 0;
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
stockholm_parsedata_Destroy(ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa)
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
    for (i = 0; i < msa->ngr; i++)
      if (pd->ogr_len[i]) free(pd->ogr_len[i]);
    free(pd->ogr_len);
  }
  return;
}
/*------------------ end, ESL_STOCKHOLM_PARSEDATA auxiliary structure -------------*/

/*****************************************************************
 * 3. Internal: parsing Stockholm line types
 *****************************************************************/ 

/* stockholm_parse_gf()
 * Line format is:
 *   #=GF <tag> <text>
 * recognized featurenames: { ID | AC | DE | AU | GA | NC | TC }
 */
static int
stockholm_parse_gf(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gf,  *tag,   *tok;
  esl_pos_t gflen, taglen, toklen;
  int       status;

  if ( (status = esl_memtok(&p, &n, " \t", &gf,  &gflen))  != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "EOL can't happen here.");
  if ( (status = esl_memtok(&p, &n, " \t", &tag, &taglen)) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GF line is missing <tag>, annotation");
  if (! esl_memstrcmp(gf, gflen, "#=GF"))                            ESL_FAIL(eslEFORMAT, afp->errmsg, "faux #=GF line?");

  if      (esl_memstrcmp(tag, taglen, "ID")) 
    {
      if ((status = esl_memtok(&p, &n, " \t", &tok, &toklen)) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "No name found on #=GF ID line");
      if (n)                                                            ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GF ID line should have only one name (no whitespace allowed)");
      if ( (status = esl_msa_SetName (msa, tok, toklen))      != eslOK) return status; /* [eslEMEM] */
    }
  else if (esl_memstrcmp(tag, taglen, "AC")) 
    {
      if ((status = esl_memtok(&p, &n, " \t", &tok, &toklen)) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "No accession found on #=GF AC line");
      if (n)                                                            ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GF AC line should have only one accession (no whitespace allowed)");
      if ((status = esl_msa_SetAccession(msa, p, n))          != eslOK) return status; /* [eslEMEM] */
    }
  else if (esl_memstrcmp(tag, taglen, "DE")) 
    {
      if ((status = esl_msa_SetDesc     (msa, p, n))          != eslOK) return status; /* [eslEMEM] */
    }
  else if (esl_memstrcmp(tag, taglen, "AU")) 
    {
      if ((status = esl_msa_SetAuthor   (msa, p, n))          != eslOK) return status; /* [eslEMEM] */
    }
  else if (esl_memstrcmp(tag, taglen, "GA"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for GA1 value on #=GF GA line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_GA1])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_GA1] = TRUE;
      } else ESL_FAIL(eslEFORMAT, afp->errmsg, "No GA threshold value found on #=GF GA line");
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for GA2 value on #=GF GA line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_GA2])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_GA2] = TRUE;
      } 
    }
  else if (esl_memstrcmp(tag, taglen, "NC"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for NC1 value on #=GF NC line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_NC1])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_NC1] = TRUE;
      } else ESL_FAIL(eslEFORMAT, afp->errmsg, "No NC threshold value found on #=GF NC line");
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for NC2 value on #=GF NC line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_NC2])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_NC2] = TRUE;
      } 
    }
  else if (esl_memstrcmp(tag, taglen, "TC"))
    {
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for TC1 value on #=GF TC line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_TC1])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_TC1] = TRUE;
      } else ESL_FAIL(eslEFORMAT, afp->errmsg, "No TC threshold value found on #=GF TC line");
      if ( (status = esl_memtok(&p, &n, " \t", &tok, &toklen)) == eslOK) {
	if (! esl_mem_IsReal(tok, toklen)) ESL_FAIL(eslEFORMAT, afp->errmsg, "Expected a real number for TC2 value on #=GF TC line");
	if (  esl_memtof(tok, toklen, &(msa->cutoff[eslMSA_TC2])) != eslOK) return status; /* [eslEMEM] */
	msa->cutset[eslMSA_TC2] = TRUE;
      } 
    }
  else 
    {
      if ((status = esl_msa_AddGF(msa, tag, taglen, p, n)) != eslOK) return status;
    }

  return eslOK;
}


/* stockholm_parse_gs()
 * Format:
 *   #=GS <seqname> <tag> <text>
 * recognized featurenames: { WT | AC | DE }
 */
static int
stockholm_parse_gs(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gs,   *seqname,   *tag,   *tok;
  esl_pos_t  gslen, seqnamelen, taglen, toklen;
  int        seqidx;
  int        status;
  
  if (esl_memtok(&p, &n, " \t", &gs,      &gslen)      != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "EOL can't happen here.");
  if (esl_memtok(&p, &n, " \t", &seqname, &seqnamelen) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GS line missing <seqname>, <tag>, annotation");
  if (esl_memtok(&p, &n, " \t", &tag,     &taglen)     != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GS line missing <tag>, annotation");
  if (! esl_memstrcmp(gs, gslen, "#=GS"))                        ESL_FAIL(eslEFORMAT, afp->errmsg, "faux #=GS line?");

  seqidx = pd->si;
  if (seqidx == pd->nseq || ! esl_memstrcmp(seqname, seqnamelen, msa->sqname[seqidx])) {
    status = stockholm_get_seqidx(msa, pd, seqname, seqnamelen, &seqidx);
  }

  if (esl_memstrcmp(tag, taglen, "WT")) 
    {
      if (esl_memtok(&p, &n, " \t", &tok, &toklen) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "no weight value found on #=GS <seqname> WT line");
      if (msa->wgt[seqidx] != -1.0)                          ESL_FAIL(eslEFORMAT, afp->errmsg, "sequence has more than one #=GS <seqname> WT line");
      if (n)                                                 ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GS <seqname> WT line should have only one field, the weight");
      if (! esl_mem_IsReal(tok, toklen))                     ESL_FAIL(eslEFORMAT, afp->errmsg, "value on #=GS <seqname> WT line isn't a real number");
      if ( esl_memtod(tok, toklen, &(msa->wgt[seqidx])))     return status; /* eslEMEM */
      msa->flags |= eslMSA_HASWGTS;
    }
  else if (esl_memstrcmp(tag, taglen, "AC"))
    {
      if (esl_memtok(&p, &n, " \t", &tok, &toklen) != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "no accession found on #=GS <seqname> AC line");
      if (msa->sqacc && msa->sqacc[seqidx])                  ESL_FAIL(eslEFORMAT, afp->errmsg, "sequence has more than one #=GS <seqname> AC accession line");
      if (n)                                                 ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GS <seqname> AC line should have only one field, the accession");
      if ((status = esl_msa_SetSeqAccession(msa, seqidx, tok, toklen)) != eslOK) return status; /* eslEMEM */
    }
  else if (esl_memstrcmp(tag, taglen, "DE"))
    {
      if (msa->sqdesc && msa->sqdesc[seqidx]) ESL_FAIL(eslEFORMAT, afp->errmsg, "sequence has more than one #=GS <seqname> DE accession line");
      if ((status = esl_msa_SetSeqDescription(msa, seqidx, p, n)) != eslOK) return status; /* eslEMEM */
    }
  else
    {
      if ((status = esl_msa_AddGS(msa, tag, taglen, seqidx, p, n)) != eslOK) return status;
    }

  pd->si = seqidx+1;	/* set guess for next sequence index */
  return eslOK;
}  
  
/* stockholm_parse_gc()
 * Format of line is:
 *   #=GC <tag> <aligned text>
 * recognized featurenames: { SS_cons | SA_cons | PP_cons | RF }
 */
static int
stockholm_parse_gc(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gc,    *tag,   *text;
  esl_pos_t  gclen,  taglen, textlen;
  int        tagidx;
  int        status;

  if (esl_memtok(&p, &n, " \t", &gc,   &gclen)    != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "EOL can't happen here.");
  if (esl_memtok(&p, &n, " \t", &tag,  &taglen)   != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GC line missing <tag>, annotation");
  if (esl_memtok(&p, &n, " \t", &text, &textlen)  != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GC line missing annotation");  
  if (! esl_memstrcmp(gc, gclen, "#=GC"))                   ESL_FAIL(eslEFORMAT, afp->errmsg, "faux #=GC line?");
  
  if (pd->nblock) 		/* Subsequent blocks */
    {
      if      (esl_memstrcmp(tag, taglen, "SS_cons") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GC_SSCONS) ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GC SS_cons line; lines in earlier block(s) were in different order?");
      else if (esl_memstrcmp(tag, taglen, "SA_cons") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GC_SACONS) ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GC SA_cons line; lines in earlier block(s) were in different order?");
      else if (esl_memstrcmp(tag, taglen, "PP_cons") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GC_PPCONS) ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GC PP_cons line; lines in earlier block(s) were in different order?");

      else if (esl_memstrcmp(tag, taglen, "RF")      && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GC_RF)     ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GC RF line; lines in earlier block(s) were in different order?");
      else if (                                         pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GC_OTHER)  ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GC line; lines in earlier block(s) were in different order?");
    }
  else				/* First block */
    {
      if (pd->bi == pd->balloc && (status = stockholm_parsedata_ExpandBlock(pd)) != eslOK) return status;

      if      (esl_memstrcmp(tag, taglen, "SS_cons"))  pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GC_SSCONS;
      else if (esl_memstrcmp(tag, taglen, "SA_cons"))  pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GC_SACONS;
      else if (esl_memstrcmp(tag, taglen, "PP_cons"))  pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GC_PPCONS;
      else if (esl_memstrcmp(tag, taglen, "RF"))       pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GC_RF;
      else                                             pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GC_OTHER;
      pd->bidx[pd->bi]      = -1;
    }

  if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GC_SSCONS)
    {
      if (pd->ssconslen != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GC SS_cons line in block");
      if ( (status = esl_strcat(&(msa->ss_cons), pd->ssconslen, text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->ssconslen += textlen;
    }
  else if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GC_SACONS)
    {
      if (pd->saconslen != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GC SA_cons line in block");
      if ((status = esl_strcat(&(msa->sa_cons), pd->saconslen, text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->saconslen += textlen;
    }
  else if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GC_PPCONS)
    {
      if (pd->ppconslen != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GC PP_cons line in block");
      if ((status = esl_strcat(&(msa->pp_cons), pd->ppconslen, text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->ppconslen += textlen;
    }
  else if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GC_RF)
    {
      if (pd->rflen != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GC RF line in block");
      if ((status = esl_strcat(&(msa->rf), pd->rflen, text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->rflen += textlen;
    }
  else
    {
      if ((status = stockholm_get_gc_tagidx(msa, pd, tag, taglen, &tagidx)) != eslOK) return status;
      
      if (pd->ogc_len[tagidx] != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GC %.*s line in block", (int) taglen, tag);
      if ((status = esl_strcat(&(msa->gc[tagidx]), pd->ogc_len[tagidx], text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->ogc_len[tagidx] += textlen;
    }

  if (pd->bi && textlen != pd->alen_b) ESL_FAIL(eslEFORMAT, afp->errmsg, "unexpected number of aligned annotation characters in #=GC %.*s line", (int) taglen, tag); 
  pd->alen_b   = textlen;
  pd->in_block = TRUE;
  pd->bi++;
  return eslOK;
}

/* A GR line is
 *   #=GR <seqname> <featurename> <text>
 * recognized featurenames: { SS | SA | PP }
 * 
 */
static int
stockholm_parse_gr(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char      *gr,   *name,    *tag,   *text;
  esl_pos_t  grlen, namelen,  taglen, textlen;
  int        seqidx, tagidx;
  int        z;
  int        status;

  if (esl_memtok(&p, &n, " \t", &gr,   &grlen)    != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "EOL can't happen here.");
  if (esl_memtok(&p, &n, " \t", &name, &namelen)  != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GR line missing <seqname>, <tag>, annotation");
  if (esl_memtok(&p, &n, " \t", &tag,  &taglen)   != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GR line missing <tag>, annotation");
  if (esl_memtok(&p, &n, " \t", &text, &textlen)  != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "#=GR line missing annotation");
  if (! esl_memstrcmp(gr, grlen, "#=GR"))                   ESL_FAIL(eslEFORMAT, afp->errmsg, "faux #=GR line?");

  /* Which seqidx is this? likely to be either pd->si-1 (#=GR following a seq) or 
   * pd->si (#=GR preceding a seq) 
   */
  if (! pd->nblock) /* First block: we're setting bidx[], blinetype[] as we see them */
    {
      if      (pd->si >= 1       && esl_memstrcmp(name, namelen, msa->sqname[pd->si-1])) seqidx = pd->si-1;
      else if (pd->si < pd->nseq && esl_memstrcmp(name, namelen, msa->sqname[pd->si]))   seqidx = pd->si;
      else  status = stockholm_get_seqidx(msa, pd, name, namelen, &seqidx);
      
      if (pd->bi == pd->balloc && (status = stockholm_parsedata_ExpandBlock(pd)) != eslOK) return status;

      if      (esl_memstrcmp(tag, taglen, "SS")) pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GR_SS;
      else if (esl_memstrcmp(tag, taglen, "SA")) pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GR_SA;
      else if (esl_memstrcmp(tag, taglen, "PP")) pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GR_PP;
      else                                       pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_GR_OTHER;
      pd->bidx[pd->bi]      = seqidx;
    }
  else 
    {				/* subsequent block(s) */
      if (pd->bi >= pd->npb) ESL_FAIL(eslEFORMAT, afp->errmsg, "more lines than expected in this alignment block; earlier blocks had fewer");

      if      (esl_memstrcmp(tag, taglen, "SS") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GR_SS)    ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GR <seqname> SS line; lines in earlier block(s) were in different order?");  
      else if (esl_memstrcmp(tag, taglen, "SA") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GR_SA)    ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GR <seqname> SA line; lines in earlier block(s) were in different order?");  
      else if (esl_memstrcmp(tag, taglen, "PP") && pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GR_PP)    ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GR <seqname> PP line; lines in earlier block(s) were in different order?");  
      else if (                                    pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_GR_OTHER) ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a #=GR line; lines in earlier block(s) were in different order?");  

      seqidx = pd->bidx[pd->bi];
      if (! esl_memstrcmp(name, namelen, msa->sqname[seqidx])) ESL_FAIL(eslEFORMAT, afp->errmsg, "unexpected sequence name %.*s; expected %s from order of earlier blocks", (int) namelen, name, msa->sqname[seqidx]);
    }

  /* Append the annotation where it belongs  */
  if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GR_SS)
    {
      if (! msa->ss) {
	ESL_ALLOC(msa->ss,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->sslen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->ss[z] = NULL; pd->sslen[z] = 0; }
      }
      if (pd->sslen[seqidx] != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GR %.*s SS line in block", (int) namelen, name);
      if (( status = esl_strcat(&(msa->ss[seqidx]), pd->sslen[seqidx], text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->sslen[seqidx] += textlen;
    }
  else if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GR_PP)
    {
      if (! msa->pp) {
	ESL_ALLOC(msa->pp,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->pplen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->pp[z] = NULL; pd->pplen[z] = 0; }
      }
      if (pd->pplen[seqidx] != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GR %.*s PP line in block", (int) namelen, name);
      if ((status = esl_strcat(&(msa->pp[seqidx]), pd->pplen[seqidx], text, textlen)) != eslOK) return status; /* [eslEMEM] */
      pd->pplen[seqidx] += textlen;
    }
  else if (pd->blinetype[pd->bi] == eslSTOCKHOLM_LINE_GR_SA)
    {
      if (! msa->sa) {
	ESL_ALLOC(msa->sa,   sizeof(char *)  * msa->sqalloc);
	ESL_ALLOC(pd->salen, sizeof(int64_t) * msa->sqalloc);
	for (z = 0; z < msa->sqalloc; z++) { msa->sa[z] = NULL; pd->salen[z] = 0; }
      }
      if (pd->salen[seqidx] != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GR %.*s SA line in block", (int) namelen, name);
      if ((status = esl_strcat(&(msa->sa[seqidx]), pd->salen[seqidx], text, textlen)) != eslOK) return status;
      pd->salen[seqidx] += textlen;
    }
  else
    {
      if ((status = stockholm_get_gr_tagidx(msa, pd, tag, taglen, &tagidx)) != eslOK) return status; /* [eslEMEM] */

      if (pd->ogr_len[tagidx][seqidx] != pd->alen) ESL_FAIL(eslEFORMAT, afp->errmsg, "more than one #=GR %.*s %.*s line in block", (int) namelen, name, (int) taglen, tag);
      if ((status = esl_strcat(&(msa->gr[tagidx][seqidx]), pd->ogr_len[tagidx][seqidx], text, textlen)) != eslOK) return status;
      pd->ogr_len[tagidx][seqidx] += textlen;
    }

  if (pd->bi && textlen != pd->alen_b) ESL_FAIL(eslEFORMAT, afp->errmsg, "unexpected number of aligned annotation characters in #=GR %.*s line", (int) taglen, tag); 
  pd->alen_b   = textlen;
  pd->in_block = TRUE;
  pd->bi++;
  return eslOK;

 ERROR:
  return status;
}
  

/* stockholm_parse_sq():
 * Format of line is:
 *   <seqname>  <aligned text>
 */
static int
stockholm_parse_sq(ESLX_MSAFILE *afp, ESL_STOCKHOLM_PARSEDATA *pd, ESL_MSA *msa, char *p, esl_pos_t n)
{
  char     *seqname, *text;
  esl_pos_t seqnamelen, textlen;
  int       seqidx = pd->si;
  int       status;
  
  if (esl_memtok(&p, &n, " \t", &seqname, &seqnamelen) != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "EOL can't happen here.");
  if (esl_memtok(&p, &n, " \t", &text,    &textlen)    != eslOK) ESL_FAIL(eslEFORMAT, afp->errmsg, "line has sequence name, but no sequence");

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
      if (pd->si < pd->nseq && esl_memstrcmp(seqname, seqnamelen, msa->sqname[seqidx])) seqidx = pd->si;
      else if ((status = stockholm_get_seqidx(msa, pd, seqname, seqnamelen, &seqidx)) != eslOK) return status; /* [eslEMEM] */

      if (pd->bi == pd->balloc && (status = stockholm_parsedata_ExpandBlock(pd)) != eslOK) return status;

      pd->blinetype[pd->bi] = eslSTOCKHOLM_LINE_SQ;
      pd->bidx[pd->bi]      = seqidx;
    }
  else 
    {				/* subsequent block(s) */
      if (pd->bi >= pd->npb)                               ESL_FAIL(eslEFORMAT, afp->errmsg, "more lines than expected in this alignment block; earlier blocks had fewer");
      if (  pd->blinetype[pd->bi] != eslSTOCKHOLM_LINE_SQ) ESL_FAIL(eslEFORMAT, afp->errmsg, "didn't expect a sequence line; lines in earlier block(s) were in different order?");
      seqidx = pd->bidx[pd->bi];

      if (! esl_memstrcmp(seqname, seqnamelen, msa->sqname[seqidx])) ESL_FAIL(eslEFORMAT, afp->errmsg, "unexpected sequence name %.*s; expected %s from order of earlier blocks", (int) seqnamelen, seqname, msa->sqname[seqidx]);
    }

#ifdef eslAUGMENT_ALPHABET 
  if (  afp->abc ) {
    status = esl_abc_dsqcat(afp->inmap, &(msa->ax[seqidx]),   &(pd->sqlen[seqidx]), text, textlen);
    if      (status == eslEINVAL) ESL_FAIL(eslEFORMAT, afp->errmsg, "invalid sequence character(s) on line");
    else if (status != eslOK)     return status;
  }
#endif
  if (! afp->abc) {
    status = esl_strmapcat (afp->inmap, &(msa->aseq[seqidx]), &(pd->sqlen[seqidx]), text, textlen);
    if      (status == eslEINVAL) ESL_FAIL(eslEFORMAT, afp->errmsg, "invalid sequence character(s) on line");
    else if (status != eslOK)     return status;
  }

  if (pd->bi && textlen != pd->alen_b)         ESL_FAIL(eslEFORMAT, afp->errmsg, "unexpected number of aligned residues parsed on line");
  if (pd->sqlen[seqidx] - pd->alen != textlen) ESL_EXCEPTION(eslEINCONCEIVABLE, "implementation assumes that no symbols are ignored in inmap; else GR, GC text annotations are messed up");
  pd->alen_b++;
  pd->nseq_b++;
  pd->in_block = TRUE;
  pd->bi++;
  pd->si = seqidx+1;
  return eslOK;
}

  

static int
stockholm_parse_comment(ESL_MSA *msa, char *p, esl_pos_t n)
{
  if (n && *p == '#')      { p++; n--; }
  while (n && isspace(*p)) { p++; n--; }

  return esl_msa_AddComment(msa, p, n);
}
/*------------- end, parsing Stockholm line types ---------------*/  



/*****************************************************************
 * 4. Internal: looking up seq, tag indices
 *****************************************************************/


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
stockholm_get_seqidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *name, esl_pos_t n, int *ret_idx)
{
  int seqidx;
  int status;
  
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

static int
stockholm_get_gr_tagidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *tag, esl_pos_t taglen, int *ret_tagidx)
{
  int tagidx;
  int z;
  int status;


  /* Find the tag, if we have it; else, add it, at tagidx = msa->ngr */
#ifdef eslAUGMENT_KEYHASH
  status = esl_keyhash_Store(msa->gr_idx, tag, taglen, &tagidx); 
  if      (status == eslEDUP) { *ret_tagidx = tagidx; return eslOK; }
  else if (status != eslOK)   goto ERROR;
#else
  for (tagidx = 0; tagidx < msa->ngr; tagidx++)
    if (esl_memstrcmp(tag, taglen, msa->gr_tag[tagidx])) 
      { *ret_tagidx = tagidx; return eslOK; }
#endif

  /* if we get here, this is a new tag we're adding. */
  ESL_REALLOC(msa->gr_tag,       sizeof(char *)    * (msa->ngr+1)); /* +1, we're allocated one new tag at a time, as needed */
  ESL_REALLOC(msa->gr,           sizeof(char **)   * (msa->ngr+1));
  ESL_REALLOC(pd->ogr_len,       sizeof(int64_t *) * (msa->ngr+1));
  msa->gr_tag[tagidx] = NULL;
  msa->gr[tagidx]     = NULL;
  pd->ogr_len[tagidx] = NULL;
  ESL_ALLOC(msa->gr[tagidx],     sizeof(char *)    * msa->sqalloc);
  ESL_ALLOC(pd->ogr_len[tagidx], sizeof(int64_t)   * msa->sqalloc);
  for (z = 0; z < msa->sqalloc; z++) {
    msa->gr[tagidx][z] = NULL;
    pd->ogr_len[tagidx][z] = 0;
  }
   
  if ( (status = esl_memstrdup(tag, taglen, &(msa->gr_tag[tagidx]))) != eslOK) goto ERROR; /* eslEMEM */
  msa->ngr++;	
  *ret_tagidx = tagidx;
  return eslOK;

 ERROR:
  *ret_tagidx = -1;
  return status;
}

static int
stockholm_get_gc_tagidx(ESL_MSA *msa, ESL_STOCKHOLM_PARSEDATA *pd, char *tag, esl_pos_t taglen, int *ret_tagidx)
{
  int tagidx;
  int status;

  /* Find the tag, if we have it; else, add it, at tagidx = msa->ngc */
#ifdef eslAUGMENT_KEYHASH
  status = esl_keyhash_Store(msa->gc_idx, tag, taglen, &tagidx); 
  if      (status == eslEDUP) { *ret_tagidx = tagidx; return eslOK; }
  else if (status != eslOK)   goto ERROR; /* eslEMEM */
#else
  for (tagidx = 0; tagidx < msa->ngc; tagidx++)
    if (esl_memstrcmp(tag, taglen, msa->gc_tag[tagidx])) 
      { *ret_tagidx = tagidx; return eslOK; }
#endif

  /* if we get here, this is a new tag we're adding. */
  ESL_REALLOC(msa->gc_tag, sizeof(char *)  * (msa->ngc+1)); /* +1, we're allocated one new tag at a time, as needed */
  ESL_REALLOC(msa->gc,     sizeof(char *)  * (msa->ngc+1));
  ESL_REALLOC(pd->ogc_len, sizeof(int64_t) * (msa->ngc+1));
  msa->gc_tag[tagidx] = NULL;
  msa->gc[tagidx]     = NULL;
  pd->ogc_len[tagidx] = 0;

  if ( (status = esl_memstrdup(tag, taglen, &(msa->gc_tag[tagidx]))) != eslOK) return status; /* eslEMEM */
  msa->ngc++;
  *ret_tagidx = tagidx;
  return eslOK;

 ERROR:
  *ret_tagidx = -1;
  return status;
}
/*------------ end, looking up seq, tag indices -----------------*/



/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
