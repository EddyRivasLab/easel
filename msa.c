/* msa.c
 * Multiple sequence alignment file i/o.
 * 
 * SRE, Thu Jan 20 08:50:43 2005 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/msa.h>


/******************************************************************************
 * Functions for the ESL_MSA object                                           *
 *     esl_msa_Create()                                                       *
 *     esl_msa_Destroy()                                                      *
 *     esl_msa_Expand()                                                       *
 *****************************************************************************/
/* Forward declarations of private MSA functions
 */
static int get_seqidx(ESL_MSA *msa, char *name, int guess, int *ret_idx);
static int set_seq_accession(ESL_MSA *msa, int seqidx, char *acc);
static int set_seq_description(ESL_MSA *msa, int seqidx, char *desc);
static int add_comment(ESL_MSA *msa, char *s);
static int add_gf(ESL_MSA *msa, char *tag, char *value);
static int add_gs(ESL_MSA *msa, char *tag, int sqidx, char *value);
static int append_gc(ESL_MSA *msa, char *tag, char *value);
static int append_gr(ESL_MSA *msa, char *tag, int sqidx, char *value);
static int verify_parse(ESL_MSA *msa, char *errbuf);

/* Function: esl_msa_Create()
 * Incept:   SRE, Sun Jan 23 08:25:26 2005 [St. Louis]
 *
 * Purpose:  Creates and initializes an <ESL_MSA> object, and returns a
 *           pointer to it. Designed to be used in three ways:
 * 
 *           1) We know exactly the dimensions of the alignment:
 *              both <nseq> and <alen>, e.g.:
 *                \begin{cchunk}
 *                   msa = esl_msa_Create(nseq, alen)
 *                \end{cchunk}
 *                    
 *           2) We know the number of sequences but not alen.
 *              (We add sequences later.) e.g.:
 *                \begin{cchunk}
 *                   msa = esl_msa_Create(nseq, 0)
 *                \end{cchunk}
 *              
 *           3) We don't even know the number of sequences, so
 *              we'll have to dynamically expand allocations.
 *              We provide an initial <nseq> that will be 
 *              expanded (by doubling) when needed. e.g.:
 *                \begin{cchunk}
 *                  msa = esl_msa_Create(16, 0);
 *                  if (msa->nseq == msa->sqalloc) esl_msa_Expand(msa);
 *                \end{cchunk}   
 *                
 *           A created <msa> can only be <Expand()>'ed if <alen> is 0
 *           (i.e. case 2,3). 
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or 0      
 *
 * Returns:  pointer to new MSA object, w/ all values initialized.
 *           Note that msa->nseq is initialized to 0, even though space
 *           is allocated.
 *           
 * Throws:   NULL on allocation failure.          
 *
 * Xref:     squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_Create(int nseq, int alen)
{
  ESL_MSA *msa;
  int  i;

  if ((msa = malloc(sizeof(ESL_MSA))) == NULL)
    ESL_ERROR_NULL(ESL_EMEM, "malloc failed");

  msa->aseq    = NULL;
  msa->sqname  = NULL;
  msa->wgt     = NULL;
  msa->alen    = alen;		/* if 0, then we're growable. */
  msa->nseq    = 0;
  msa->flags   = 0;
  msa->type    = 0;		/* =eslUNKNOWN; no alphabet type yet */
  msa->name    = NULL;
  msa->desc    = NULL;
  msa->acc     = NULL;
  msa->au      = NULL;
  msa->ss_cons = NULL;
  msa->sa_cons = NULL;
  msa->rf      = NULL;
  msa->sqacc   = NULL;
  msa->sqdesc  = NULL;
  msa->ss      = NULL;
  msa->sa      = NULL;
  for (i = 0; i < eslMSA_NCUTS; i++) {
    msa->cutoff[i] = 0.;
    msa->cutset[i] = FALSE;
  }
#ifdef ESL_KEYHASH_INCLUDED	/* SRE: CHECK ME */
  msa->index   = GKIInit();
#endif
  msa->sqalloc = nseq;
  msa->sqlen   = NULL;
  msa->sslen   = NULL;
  msa->salen   = NULL;
  msa->lastidx = 0;

  /* Unparsed markup, including comments and Stockholm tags.
   * GS, GC, and GR Stockholm tags require keyhash augmentation
   */
  msa->comment        = NULL;
  msa->ncomment       = 0;
  msa->alloc_ncomment = 0;

  msa->gf_tag         = NULL;
  msa->gf             = NULL;
  msa->ngf            = 0;
  msa->alloc_ngf      = 0;

#ifdef ESL_KEYHASH_INCLUDED
  msa->gs_tag         = NULL;
  msa->gs             = NULL;
  msa->gs_idx         = NULL;
  msa->ngs            = 0;

  msa->gc_tag         = NULL;
  msa->gc             = NULL;
  msa->gc_idx         = NULL;
  msa->ngc            = 0;

  msa->gr_tag         = NULL;
  msa->gr             = NULL;
  msa->gr_idx         = NULL;
  msa->ngr            = 0;
#endif /*keyhash augmentation*/

  /* Allocation, round 2.
   */
  if ((msa->aseq   = malloc(sizeof(char *) * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(ESL_EMEM, "malloc failed"); }
  if ((msa->sqname = malloc(sizeof(char *) * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(ESL_EMEM, "malloc failed"); }
  if ((msa->wgt    = malloc(sizeof(float)  * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(ESL_EMEM, "malloc failed"); }
  if ((msa->sqlen  = malloc(sizeof(int)    * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(ESL_EMEM, "malloc failed"); }

  /* Initialize at the second level.
   */
  for (i = 0; i < nseq; i++)
    {
      msa->aseq[i]   = NULL;
      msa->sqname[i] = NULL;
      msa->sqlen[i]  = 0;
      msa->wgt[i]    = -1.0;	/* "unset so far" */
    }

  /* Allocation, round 3.
   */
  if (alen != 0)
    {
      for (i = 0; i < nseq; i++)
	if ((msa->aseq[i] = malloc(sizeof(char) * (alen+1))) == NULL) 
	  { 
	    esl_msa_Destroy(msa); 
	    ESL_ERROR_NULL(ESL_EMEM, "malloc failed"); 
	  }
    }      

  return msa;
}


/* Function:  esl_msa_Destroy()
 * Incept:    SRE, Sun Jan 23 08:26:02 2005 [St. Louis]
 *
 * Purpose:   Destroys <msa>.
 *
 * Xref:      squid's MSADestroy().
 */
void
esl_msa_Destroy(ESL_MSA *msa)
{
  if (msa == NULL) return;

  esl_Free2D((void **) msa->aseq,   msa->nseq);
  esl_Free2D((void **) msa->sqname, msa->nseq);
  esl_Free2D((void **) msa->sqacc,  msa->nseq);
  esl_Free2D((void **) msa->sqdesc, msa->nseq);
  esl_Free2D((void **) msa->ss,     msa->nseq);
  esl_Free2D((void **) msa->sa,     msa->nseq);

  if (msa->sqlen   != NULL) free(msa->sqlen);
  if (msa->wgt     != NULL) free(msa->wgt);

  if (msa->name    != NULL) free(msa->name);
  if (msa->desc    != NULL) free(msa->desc);
  if (msa->acc     != NULL) free(msa->acc);
  if (msa->au      != NULL) free(msa->au);
  if (msa->ss_cons != NULL) free(msa->ss_cons);
  if (msa->sa_cons != NULL) free(msa->sa_cons);
  if (msa->rf      != NULL) free(msa->rf);
  if (msa->sslen   != NULL) free(msa->sslen);
  if (msa->salen   != NULL) free(msa->salen);
  
  esl_Free2D((void **) msa->comment, msa->ncomment);
  esl_Free2D((void **) msa->gf_tag,  msa->ngf);
  esl_Free2D((void **) msa->gf,      msa->ngf);
#ifdef ESL_KEYHASH_INCLUDED
  esl_Free2D((void **) msa->gs_tag,  msa->ngs);
  esl_Free3D((void ***)msa->gs,      msa->ngs, msa->nseq);
  esl_Free2D((void **) msa->gc_tag,  msa->ngc);
  esl_Free2D((void **) msa->gc,      msa->ngc);
  esl_Free2D((void **) msa->gr_tag,  msa->ngr);
  esl_Free3D((void ***)msa->gr,      msa->ngr, msa->nseq);
  GKIFree(msa->index);
  GKIFree(msa->gs_idx);
  GKIFree(msa->gc_idx);
  GKIFree(msa->gr_idx);
#endif /* keyhash augmentation */  

  free(msa);
  return;
}



/* Function:  esl_msa_Expand()
 * Incept:    SRE, Sun Jan 23 08:26:30 2005 [St. Louis]
 *
 * Purpose:   Double the current sequence allocation in <msa>.
 *            Typically used when we're reading an alignment sequentially 
 *            from a file, so we don't know nseq 'til we're done.
 *
 * Returns:   <ESL_OK> on success.
 * 
 * Throws:    <ESL_EMEM> on reallocation failure; <msa> is undamaged,
 *            and the caller may attempt to recover from the error.
 *            
 *            <ESL_EINVAL> if <msa> is not growable: its <alen> field
 *            must be 0 to be growable.
 *
 * Xref:      squid's MSAExpand(), 1999.
 */
int
esl_msa_Expand(ESL_MSA *msa)
{
  int   old, new;		/* old & new allocation sizes */
  void *p;			/* tmp ptr to realloc'ed memory */
  int   i;

  if (msa->alen > 0) 
    ESL_ERROR(ESL_EINVAL, "that MSA is not growable");

  old = msa->sqalloc;
  new = 2*old;

  ESL_REALLOC(msa->aseq,   p, sizeof(char *) * new);
  ESL_REALLOC(msa->sqname, p, sizeof(char *) * new);
  ESL_REALLOC(msa->wgt,    p, sizeof(float)  * new);
  ESL_REALLOC(msa->sqlen,  p, sizeof(int)    * new);

  if (msa->ss != NULL) 
    {
      ESL_REALLOC(msa->ss,    p, sizeof(char *) * new);
      ESL_REALLOC(msa->sslen, p, sizeof(int) * new);
    }
  
  if (msa->sa != NULL) 
    {
      ESL_REALLOC(msa->sa,    p, sizeof(char *) * new);
      ESL_REALLOC(msa->salen, p, sizeof(int) * new);
    }

  if (msa->sqacc != NULL)
    ESL_REALLOC(msa->sqacc,  p, sizeof(char *) * new);

  if (msa->sqdesc != NULL)
    ESL_REALLOC(msa->sqdesc, p, sizeof(char *) * new);

  for (i = old; i < new; i++)
    {
      msa->aseq[i]   = NULL;
      msa->sqname[i] = NULL;
      msa->wgt[i]    = -1.0;	/* -1.0 means "unset so far" */
      msa->sqlen[i]  = 0;

      if (msa->ss != NULL) 
	{
	  msa->ss[i] = NULL;
	  msa->sslen[i] = 0;
	}
      if (msa->sa != NULL) 
	{ 
	  msa->sa[i] = NULL;
	  msa->salen[i] = 0;
	}
      if (msa->sqacc  != NULL) msa->sqacc[i]  = NULL;
      if (msa->sqdesc != NULL) msa->sqdesc[i] = NULL;
    }

#ifdef ESL_KEYHASH_INCLUDED
  /* AUGMENTATION-SPECIFIC: optional stockholm tag section
   */
  /* Reallocate and re-init for unparsed #=GS tags, if we have some.
   * gs is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gs != NULL)
    for (i = 0; i < msa->ngs; i++)
      {
	if (msa->gs[i] != NULL)
	  {
	    ESL_REALLOC(msa->gs[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gs[i][j] = NULL;
	  }
      }
  /* Reallocate and re-init for unparsed #=GR tags, if we have some.
   * gr is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gr != NULL)
    for (i = 0; i < msa->ngr; i++)
      {
	if (msa->gr[i] != NULL)
	  {
	    ESL_REALLOC(msa->gr[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gr[i][j] = NULL;
	  }
      }
#endif /* keyhash augmentation */

  msa->sqalloc = new;
  return ESL_OK;
}

/* get_seqidx()
 * 
 * Find the index of a given sequence <name> in an <msa>.
 * If caller has a good guess (for instance, the sequences are
 * coming in a previously seen order in a block of seqs or annotation),
 * the caller can pass this information in <guess>, or -1 if
 * it has no guess.
 * 
 * This function behaves differently depending on whether
 * keyhash augmentation is available or not. Without keyhashing,
 * the name is identified by bruteforce search of the names
 * in the <msa>. With keyhashing, we hash search, which should
 * improve performance for large alignments.
 * 
 * If the name does not already exist in the MSA, then it
 * is assumed to be a new sequence name that we need to store.
 * seqidx is set to msa:nseq, the MSA is Expand()'ed if necessary
 * to make room, the name is stored in msa:sqname[msa:nseq],
 * (and in the hash table, if we're keyhash augmented)
 * and msa:nseq is incremented.
 *
 * Returns:  <ESL_OK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws:   <ESL_EMEM> if we try to add a name and allocation fails.
 *           <ESL_EINVAL> if we try to add a name to a non-growable MSA.
 */
static int
get_seqidx(ESL_MSA *msa, char *name, int guess, int *ret_idx)
{
  int seqidx;
  int status;

  *ret_idx = -1;

  /* can we guess? */
  if (guess >= 0 && 
      guess < msa->nseq && 
      strcmp(name, msa->sqname[guess]) == 0) 
    { *ret_idx = guess; return ESL_OK; }

  /* Else look it up - either brute force
   * or, if we're keyhash-augmented, by hashing.
   */
#ifdef ESL_KEYHASH_INCLUDED                  
  if ((seqidx = GKIKeyIndex(msa->index, name)) >= 0)
    { *ret_idx = seqidx; return ESL_OK; }
				/* else, it's a new name */
  seqidx = GKIStoreKey(msa->index, name);
#else
  for (seqidx = 0; seqidx < msa->nseq; seqidx++)
    if (strcmp(msa->sqname[seqidx], name) == 0) break;
  if (seqidx < msa->nseq) 
    { *ret_idx = seqidx; return ESL_OK; }
#endif

  /* If we reach here, then this is a new name that we're
   * adding.
   */
  if (seqidx >= msa->sqalloc &&  
     (status = esl_msa_Expand(msa)) != ESL_OK)
    return status; 
    
  status = esl_strdup(name, -1, &(msa->sqname[seqidx]));
  msa->nseq++;
  if (ret_idx != NULL) *ret_idx = seqidx;
  return status;
}

/* set_seq_accession()
 *
 * Sets the sequence accession field for sequence
 * number <seqidx> in an alignment <msa>, by
 * duplicating the string <acc>.
 *
 * Returns:   <ESL_OK> on success.
 * Throws:    <ESL_EMEM> on allocation failure.
 */
static int
set_seq_accession(ESL_MSA *msa, int seqidx, char *acc)
{
  int i;

  /* If this is the first acccession, we have to
   * initialize the whole optional array.
   */
  if (msa->sqacc == NULL) 
    {
      ESL_MALLOC(msa->sqacc, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqacc[i] = NULL;
    }
  /* If we already had an accession, that's weird, but free it. 
   */
  if (msa->sqacc[seqidx] != NULL) free(msa->sqacc[seqidx]);
  return (esl_strdup(acc, -1, &(msa->sqacc[seqidx])));
}

/* set_seq_description()
 *
 * Set the sequence description field for sequence number
 * <seqidx> in an alignment <msa> by copying the string <desc>.
 *
 * Returns:  <ESL_OK> on success.
 * 
 * Throws:   <ESL_EMEM> on allocation failure.
 */
int
set_seq_description(ESL_MSA *msa, int seqidx, char *desc)
{
  int i;

  if (msa->sqdesc == NULL) 
    {
      ESL_MALLOC(msa->sqdesc, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqdesc[i] = NULL;
  }
  if (msa->sqdesc[i] != NULL) free(msa->sqdesc[i]);
  return (esl_strdup(desc, -1, &(msa->sqdesc[seqidx])));
}


/* add_comment()
 * SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Add an (unparsed) comment line to the MSA structure, allocating as
 * necessary.
 *
 * Args:     msa - a multiple alignment
 *           s   - comment line to add
 *
 * Returns:  <ESL_OK> on success.
 */
static int
add_comment(ESL_MSA *msa, char *s)
{
  void *p;
  int   status;

  /* If this is our first recorded comment, we need to malloc();
   * and if we've filled available space, we need to realloc().
   */
  if (msa->comment == NULL) {
    ESL_MALLOC(msa->comment, sizeof(char *) * 16);
    msa->alloc_ncomment = 16;
  }
  if (msa->ncomment == msa->alloc_ncomment) {
    ESL_REALLOC(msa->comment, p, sizeof(char *) * msa->alloc_ncomment * 2);
    msa->alloc_ncomment *= 2;
  }
  status = esl_strdup(s, -1, &(msa->comment[msa->ncomment]));
  msa->ncomment++;
  return status;
}


/* add_gf()
 * 
 * Add an unparsed #=GF markup line to the MSA, allocating
 * as necessary. <tag> is the GF markup tag; <value> is
 * the text associated w/ that tag.
 * 
 * Returns ESL_OK on success. 
 * Throws ESL_EMEM on allocation failure.
 */
static int
add_gf(ESL_MSA *msa, char *tag, char *value)
{  
  void *p;
  int   n;
  int   status;

  /* If this is our first recorded unparsed #=GF line, we need to malloc().
   */
  if (msa->gf_tag == NULL) {
    ESL_MALLOC(msa->gf_tag, sizeof(char *) * 16);
    ESL_MALLOC(msa->gf,     sizeof(char *) * 16);
    msa->alloc_ngf = 16;
  }
  /* or if we're out of room for new GF's, realloc() by doubling
   */
  if (msa->ngf == msa->alloc_ngf) {
    n = msa->alloc_ngf * 2;
    ESL_REALLOC(msa->gf_tag, p, sizeof(char *) * n);
    ESL_REALLOC(msa->gf,     p, sizeof(char *) * n);
    msa->alloc_ngf = n;
  }

  status = esl_strdup(tag, -1, &(msa->gf_tag[msa->ngf]));
  if (status != ESL_OK) return status;
  status = esl_strdup(value, -1, &(msa->gf[msa->ngf]));
  if (status != ESL_OK) return status;
  msa->ngf++;

  return ESL_OK;
}


/* add_gs()
 *
 * Adds an unparsed #=GS markup line to the MSA structure, allocating
 * as necessary.
 *           
 * It's possible that we could get more than one of the same type of
 * GS tag per sequence; for example, "DR PDB;" structure links in
 * Pfam.  Hack: handle these by appending to the string, in a \n
 * separated fashion.
 *
 * Args:     msa    - multiple alignment structure
 *           tag    - markup tag (e.g. "AC")
 *           sqidx  - index of sequence to assoc markup with (0..nseq-1)
 *           value  - markup (e.g. "P00666")
 *
 * Returns:  <ESL_OK> on success
 * Throws:   <ESL_EMEM> on allocation failure
 */
int
add_gs(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int   tagidx;
  int   i;
  int   status;

  /* first GS tag? init w/ malloc  */
  if (msa->gs_tag == NULL)	
    {
#ifdef ESL_KEYHASH_INCLUDED
      msa->gs_idx = GKIInit();
      tagidx      = GKIStoreKey(msa->gs_idx, tag);
      SQD_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_MALLOC(msa->gs_tag, sizeof(char *));  /* one at a time. */
      ESL_MALLOC(msa->gs,     sizeof(char **));
      ESL_MALLOC(msa->gs[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->gs[0][i] = NULL;
    }
  else 
    {
      /* Get a tagidx for this GS tag.
       * tagidx < ngs; we already saw this tag;
       * tagidx == ngs; this is a new one.
       */
#ifdef ESL_KEYHASH_INCLUDED
      tagidx  = GKIKeyIndex(msa->gs_idx, tag); 
      if (tagidx < 0) {
	tagidx = GKIStoreKey(msa->gs_idx, tag);
	SQD_DASSERT1((tagidx == msa->ngs));
      }
#else
      for (tagidx = 0; tagidx < msa->ngs; tagidx++)
	if (strcmp(msa->gs_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocation (in blocks of 1) */
      if (tagidx == msa->ngs ) 
	{
	  ESL_REALLOC(msa->gs_tag, p, (msa->ngs+1) * sizeof(char *));
	  ESL_REALLOC(msa->gs,     p, (msa->ngs+1) * sizeof(char **));
	  ESL_MALLOC(msa->gs[msa->ngs], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gs[msa->ngs][i] = NULL;
	}
    }

  /* Store the tag, if it's new.
   */
  if (tagidx == msa->ngs) 
    {
      status = esl_strdup(tag, -1, &(msa->gs_tag[tagidx]));
      if (status != ESL_OK) return status;
      msa->ngs++;
    }
  
  /* Store the annotation on the sequence.
   * If seq is unannotated, dup the value; if
   * seq already has a GS annotation, cat a \n, then cat the value.
   */
  if (msa->gs[tagidx][sqidx] == NULL)
    {
      status = esl_strdup(value, -1, &(msa->gs[tagidx][sqidx]));
      if (status != ESL_OK) return status;
    }
  else 
    {			
      int n1,n2;
      n1 = strlen(msa->gs[tagidx][sqidx]);
      n2 = strlen(value);
      ESL_REALLOC(msa->gs[tagidx][sqidx], p, sizeof(char) * (n1+n2+2));
      msa->gs[tagidx][sqidx][n1] = '\n';
      strcpy(msa->gs[tagidx][sqidx]+n1+1, value);
    }
  return ESL_OK;
} 

/* append_gc()
 *
 * Add an unparsed #=GC markup line to the MSA structure, allocating
 * as necessary.
 *           
 * When called multiple times for the same tag, appends value strings
 * together -- used when parsing multiblock alignment files, for
 * example.
 *
 * Args:     msa   - multiple alignment structure
 *           tag   - markup tag (e.g. "CS")
 *           value - markup, one char per aligned column      
 *
 * Returns:  <ESL_OK> on success
 * 
 * Throws:   <ESL_EMEM> on allocation failure
 */
static int
append_gc(ESL_MSA *msa, char *tag, char *value)
{
  int   tagidx;
  void *p;

  /* Is this an unparsed tag name that we recognize?
   * If not, handle adding it to index, and reallocating
   * as needed.
   */
  if (msa->gc_tag == NULL)	/* first tag? init w/ malloc  */
    {
#ifdef ESL_KEYHASH_INCLUDED
      msa->gc_idx = GKIInit();
      tagidx      = GKIStoreKey(msa->gc_idx, tag);      
      SQD_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_MALLOC(msa->gc_tag, sizeof(char **));
      ESL_MALLOC(msa->gc,     sizeof(char **));
      msa->gc[0]  = NULL;
    }
  else
    {			/* new tag? */
      /* get tagidx for this GC tag. existing tag: <ngc; new: == ngc. */
#ifdef ESL_KEYHASH_INCLUDED
      tagidx  = GKIKeyIndex(msa->gc_idx, tag); 
      if (tagidx < 0) 
	{		
	  tagidx = GKIStoreKey(msa->gc_idx, tag);
	  SQD_DASSERT1((tagidx == msa->ngc));
	}
#else
      for (tagidx = 0; tagidx < msa->ngc; tagidx++)
	if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocate, in block of one tag at a time
       */
      ESL_REALLOC(msa->gc_tag, p, (msa->ngc+1) * sizeof(char **));
      ESL_REALLOC(msa->gc,     p, (msa->ngc+1) * sizeof(char **));
      msa->gc[tagidx] = NULL;
    }
  /* new tag? store it.
   */
  if (tagidx == msa->ngc) 
    {
      if (esl_strdup(tag, -1, &(msa->gc_tag[tagidx])) != ESL_OK) 
	return ESL_EMEM;
      msa->ngc++;
    }
  return (esl_strcat(&(msa->gc[tagidx]), -1, value, -1));
}

/* append_gr()
 * SRE, Thu Jun  3 06:34:38 1999 [Madison]
 *
 * Add an unparsed #=GR markup line to the MSA structure, allocating
 * as necessary.
 *           
 * When called multiple times for the same tag, appends value strings
 * together -- used when parsing multiblock alignment files, for
 * example.
 *
 * Args:     msa    - multiple alignment structure
 *           tag    - markup tag (e.g. "SS")
 *           sqidx  - index of seq to assoc markup with (0..nseq-1)
 *           value  - markup, one char per aligned column      
 *
 * Returns:  <ESL_OK> on success.
 * 
 * Throws:   <ESL_EMEM> on allocation failure.
 */
static int
append_gr(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int tagidx;
  int i;

  if (msa->gr_tag == NULL)	/* first tag? init w/ malloc  */
    {
#ifdef ESL_KEYHASH_INCLUDED
      msa->gr_idx = GKIInit();
      tagidx      = GKIStoreKey(msa->gr_idx, tag);
      SQD_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_MALLOC(msa->gr_tag, sizeof(char *));
      ESL_MALLOC(msa->gr,     sizeof(char **));
      ESL_MALLOC(msa->gr[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) 
	msa->gr[0][i] = NULL;
    }
  else 
    {
      /* get tagidx for this GR tag. existing<ngr; new=ngr.
       */
#ifdef ESL_KEYHASH_INCLUDED
      tagidx  = GKIKeyIndex(msa->gr_idx, tag); 
      if (tagidx < 0) 
	{	
	  tagidx = GKIStoreKey(msa->gr_idx, tag);
	  SQD_DASSERT1((tagidx == msa->ngr));
	}
#else
      for (tagidx = 0; tagidx < msa->ngr; tagidx++)
	if (strcmp(msa->gr_tag[tagidx], tag) == 0) break;
#endif
      /* if a new tag, realloc for it */      
      if (tagidx == msa->ngr)
	{ 
	  ESL_REALLOC(msa->gr_tag, p, (msa->ngr+1) * sizeof(char *));
	  ESL_REALLOC(msa->gr,     p, (msa->ngr+1) * sizeof(char **));
	  ESL_MALLOC(msa->gr[msa->ngr], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gr[msa->ngr][i] = NULL;
	}
    }

  if (tagidx == msa->ngr) 
    {
      if (esl_strdup(tag, -1, &(msa->gr_tag[tagidx])) != ESL_OK)
	return ESL_EMEM;
      msa->ngr++;
    }
  return (esl_strcat(&(msa->gr[tagidx][sqidx]), -1, value, -1));
}

/* verify_parse()
 *
 * Last function called after a multiple alignment parser thinks it's
 * done. Checks that parse was successful; makes sure required
 * information is present; makes sure required information is
 * consistent. Some fields that are only use during parsing may be
 * freed (sqlen, for example), and some fields are finalized now
 * (msa::alen is set, for example). 
 * 
 * <errbuf> is a place to sprintf an informative message about the
 * reason for a parse error. The caller provides an <errbuf>
 * of at least 512 bytes.
 *
 * Returns:  <ESL_OK>, and errbuf is set to an empty string.
 *           
 * Throws:   <ESL_EFORMAT> if a problem is detected, and an
 *           informative message about the failure is in errbuf.
 */
static int
verify_parse(ESL_MSA *msa, char *errbuf)
{
  int idx;

  if (msa->nseq == 0) 
    {
      sprintf(errbuf, 
	      "MSA parse error: no sequences were found for alignment %.128s",
	      msa->name != NULL ? msa->name : "");
      return ESL_EFORMAT;
    }

  /* set alen, until proven otherwise; we'll check that the other seqs
   * have the same length later.
   */
  msa->alen = msa->sqlen[0];

  /* We can rely on msa->sqname[] being valid for any index,
   * because of the way the line parsers always store any name
   * they add to the index.
   */
  for (idx = 0; idx < msa->nseq; idx++)
    {
      /* aseq is required. */
      if (msa->aseq[idx] == NULL) 
	{
	  sprintf(errbuf,
		  "MSA %.128s parse error: no sequence for %.128s",
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]); 
	  return ESL_EFORMAT;
	}

      /* either all weights must be set, or none of them */
      if ((msa->flags & eslMSA_HASWGTS) && msa->wgt[idx] == -1.0)
	{
	  sprintf(errbuf,
		  "MSA %.128s parse error: expected a weight for seq %.128s", 
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]);
	  return ESL_EFORMAT;
	}

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	{
	  sprintf(errbuf,
	  "MSA %.128s parse error: sequence %.128s: length %d, expected %d",
		  msa->name != NULL ? msa->name : "",
		  msa->sqname[idx], msa->sqlen[idx], msa->alen);
	  return ESL_EFORMAT;
	}

      /* if individual SS is present, it must have length right too */
      if (msa->ss != NULL &&
	  msa->ss[idx] != NULL && 
	  msa->sslen[idx] != msa->alen) 
	{
	  sprintf(errbuf,
	  "MSA %.128s parse error: GR SS for %.128s: length %d, expected %d",
		  msa->name != NULL ? msa->name : "",
		  msa->sqname[idx], msa->sslen[idx], msa->alen);
	  return ESL_EFORMAT;
	}
				/* if SA is present, must have length right */
      if (msa->sa != NULL && 
	  msa->sa[idx] != NULL && 
	  msa->salen[idx] != msa->alen) 
	{
	  sprintf(errbuf,
	  "MSA %.128s parse error: GR SA for %.128s: length %d, expected %d",
		  msa->name != NULL ? msa->name : "",
		  msa->sqname[idx], msa->salen[idx], msa->alen);
	  return ESL_EFORMAT;
	}
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC SS_cons markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->ss_cons), msa->alen);
      return ESL_EFORMAT;
    }

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC SA_cons markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->sa_cons), msa->alen);
      return ESL_EFORMAT;
    }

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC RF markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->rf), msa->alen);
      return ESL_EFORMAT;
    }

  /* If no weights were set, set 'em all to 1.0 */
  if (!(msa->flags & eslMSA_HASWGTS))
    for (idx = 0; idx < msa->nseq; idx++)
      msa->wgt[idx] = 1.0;

  /* Clean up a little from the parser */
  if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
  if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
  if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }

  return ESL_OK;
}



/*---------------------- end of ESL_MSA functions ---------------------------*/




/******************************************************************************
 * Functions for an ESL_MSAFILE object                                        *
 *     esl_msafile_Open()                                                     *
 *     esl_msafile_Close()                                                    *
 *****************************************************************************/

/* Function: esl_msafile_Open()
 * Date:     SRE, Sun Jan 23 08:30:33 2005 [St. Louis]
 *
 * Purpose:  Open an alignment database file <filename> and prepare for
 *           reading one alignment, or sequentially in the case of *
 *           multiple MSA * databases (e.g. Stockholm format); returns
 *           the opened file pointer in <ret_msafp>.
 *          
 *           There are one or two special cases for <filename>. If
 *           <filename> is "-", then the alignment is read from
 *           <stdin>. If <filename> * ends in ".gz", then the file is
 *           assumed to be compressed * by gzip, and it is opened as a
 *           pipe from <gunzip -dc>. (Auto-decompression of gzipp'ed files
 *           is only available on POSIX-compliant systems, when 
 *           <ESL_POSIX_AUGMENTATION> is defined at compile-time.)
 *          
 *           If <env> is non-NULL, then we look for <filename> in
 *           one or more directories in a colon-delimited list
 *           that is the value of the environment variable <env>.
 *           (For example, if we had 
 *              <setenv HMMERDB /nfs/db/Pfam:/nfs/db/Rfam> 
 *           in the environment, a profile HMM application
 *           might pass "HMMERDB" as <env>.)
 *          
 *          The file is asserted to be in format <fmt>, which is
 *          either a known format like <eslMSAFILE_STOCKHOLM>, or
 *          <eslMSAFILE_UNKNOWN>; if <fmt> is <eslMSAFILE_UNKNOWN>,
 *          then format autodetection is invoked.
 *
 * Returns:  <ESL_OK> on success, and <ret_msafp> is set to point at
 *           an open <MSAFILE *>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           Returns <ESL_ENOTFOUND> if <filename> cannot be opened,
 *           or <ESL_EFORMAT> if autodetection is attempted and 
 *           format cannot be determined.
 *           
 * Throws:   <ESL_EMEM> on allocation failure.
 *           <ESL_EINVAL> if format autodetection is attempted on 
 *           stdin or a gunzip pipe.
 *
 * Xref:     squid's MSAFileOpen(), 1999.
 */
int
esl_msafile_Open(char *filename, int format, char *env, 
		 ESL_MSAFILE **ret_msafp)
{
  ESL_MSAFILE *afp;
  char *ssifile;
  int  n;
  int  status;
  
  if (ret_msafp != NULL) *ret_msafp = NULL;

  if ((afp = malloc(sizeof(ESL_MSAFILE))) == NULL)
    ESL_ERROR(ESL_EMEM, "malloc failed");
  afp->f          = NULL;
  afp->fname      = NULL;
  afp->linenumber = 0;
  afp->errbuf[0]  = '\0';
  afp->buf        = NULL;
  afp->buflen     = 0;
  afp->do_gzip    = FALSE;
  afp->do_stdin   = FALSE;
#ifdef ESL_SSI_INCLUDED		/* SSI augmentation */
  afp->ssi = NULL;
#endif  

  n        = strlen(filename);
  ssifile  = NULL;

  if (strcmp(filename, "-") == 0)
    {
      afp->f         = stdin;
      afp->do_stdin  = TRUE; 
      afp->do_gzip   = FALSE;
      status         = esl_strdup("[STDIN]", -1, &(afp->fname));
      if (status != ESL_OK) 
	{ esl_msafile_Close(afp); return ESL_EMEM; }
    }
#ifdef ESL_POSIX_AUGMENTATION
  /* popen(), pclose() aren't portable to non-POSIX systems; 
   * disable this section in strict ANSI C mode.
   */
  /* tricky: if n= length of a string s, then
   * s+n-i repositions pointer s at the last i chars
   * of the string.
   */
  else if (strcmp(filename+n-3, ".gz") == 0)
    {
      char *cmd;

      /* Note that popen() will return "successfully"
       * if file doesn't exist, because gzip works fine
       * and prints an error! So we have to check for
       * existence of file ourself.
       */
      if (! esl_FileExists(filename))
	{ 
	  esl_msafile_Close(afp); 
	  return ESL_ENOTFOUND; 
	}
      if ((cmd = malloc(sizeof(char) * (n+1+strlen("gzip -dc ")))) == NULL)
	{ 
	  esl_msafile_Close(afp);
	  ESL_ERROR(ESL_EMEM, "malloc for cmd failed"); 
	}
      sprintf(cmd, "gzip -dc %s", filename);
      if ((afp->f = popen(cmd, "r")) == NULL)
	{ 
	  esl_msafile_Close(afp); 
	  return ESL_ENOTFOUND; 
	}
      status = esl_strdup(filename, n, &(afp->fname));
      afp->do_stdin = FALSE;
      afp->do_gzip  = TRUE;
      if (status != ESL_OK)
	{ esl_msafile_Close(afp); return ESL_EMEM; }
    }
#endif /*ESL_POSIX_AUGMENTATION*/
  else
    {
      char *envfile;

      /* When we open a file, it may be either in the current
       * directory, or in the directory indicated by the env
       * argument - and we construct an SSI filename accordingly.
       * (Whether or not we're SSI augmented, in fact, for simplicity.)
       */
      if ((afp->f = fopen(filename, "r")) != NULL)
	{
	  esl_FileNewSuffix(filename, "ssi", &ssifile);	/* FIXME: check return status */
	}
      else if (esl_FileEnvOpen(filename, env, &(afp->f), &envfile) == ESL_OK)
	{
	  esl_FileNewSuffix(envfile, "ssi", &ssifile);
	  free(envfile);
	}
      else 
	{ esl_msafile_Close(afp); return ESL_ENOTFOUND; }

      afp->do_stdin = FALSE;
      afp->do_gzip  = FALSE;
      status = esl_strdup(filename, n, &(afp->fname));
      if (status != ESL_OK)
	{ esl_msafile_Close(afp); return ESL_EMEM; }
    }

  /* If augmented by SSI indexing:
   * Open the SSI index file. If it doesn't exist, or
   * it's corrupt, or some error happens, afp->ssi stays NULL.
   */
#ifdef ESL_SSI_INCLUDED
  SSIOpen(ssifile, &(afp->ssi)); /* FIXME */
#endif
  if (ssifile != NULL) free (ssifile);

  /* Invoke autodetection if we haven't already been told what
   * to expect.
   */
  if (format == eslMSAFILE_UNKNOWN)
    {
      if (afp->do_stdin == TRUE || afp->do_gzip)
	{
	  esl_msafile_Close(afp);
	  ESL_ERROR(ESL_EINVAL, 
		    "Can't autodetect alignment file fmt in stdin, gzip pipe");
	}

      format = esl_msa_GuessFileFormat(afp);
      if (format == eslMSAFILE_UNKNOWN)
	{
	  esl_msafile_Close(afp);
	  return ESL_EFORMAT;
	}
    }

  afp->format     = format;
  if (ret_msafp != NULL) *ret_msafp = afp; else esl_msafile_Close(afp);
  return ESL_OK;
}

/* Function:  esl_msafile_Close()
 * Incept:    SRE, Sun Jan 23 08:18:39 2005 [St. Louis]
 *
 * Purpose:   Close an open <ESL_MSAFILE>.
 *
 * Xref:      squid's MSAFileClose().
 */
void
esl_msafile_Close(ESL_MSAFILE *afp)
{
  if (afp == NULL) return;

#ifdef ESL_POSIX_AUGMENTATION /* gzip functionality only on POSIX systems */
  if (afp->do_gzip)       pclose(afp->f);
#endif
  if (! afp->do_stdin)    fclose(afp->f);
  if (afp->fname != NULL) free(afp->fname);
  if (afp->buf  != NULL)  free(afp->buf);
#ifdef ESL_SSI_INCLUDED
  if (afp->ssi  != NULL)  SSIClose(afp->ssi); /* FIXME */
#endif
  free(afp);
}

/* msafile_getline():
 * load the next line of <afp> into <afp>::<buf>. 
 * Returns ESL_OK on success, ESL_EOF on normal eof.
 * Throws ESL_EMEM on alloc failure.
 */
int
msafile_getline(ESL_MSAFILE *afp)
{
  int status;
  status = esl_fgets(&(afp->buf), &(afp->buflen), afp->f);
  afp->linenumber++;
  return status;
}

/*-------------------- end of ESL_MSAFILE functions -------------------------*/




/******************************************************************************
 * Functions for i/o of Stockholm format                                      *
 *****************************************************************************/

/* Forward declarations of private Stockholm i/o functions
 */
static int is_blankline(char *s);
static int parse_gf(ESL_MSA *msa, char *buf);
static int parse_gs(ESL_MSA *msa, char *buf);
static int parse_gc(ESL_MSA *msa, char *buf);
static int parse_gr(ESL_MSA *msa, char *buf);
static int parse_comment(ESL_MSA *msa, char *buf);
static int parse_sequence(ESL_MSA *msa, char *buf);

/* Function:  esl_msa_ReadStockholm()
 * Incept:    SRE, Sun Jan 23 08:33:32 2005 [St. Louis]
 *
 * Purpose:   Parse the next alignment from an open Stockholm format alignment
 *            file <afp>, leaving the alignment in <ret_msa>.
 *
 * Returns:   <ESL_OK> on success, and the alignment is in <ret_msa>.
 *            <ESL_EOF> if there are no more alignments in <afp>. 
 *            <ESL_EFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem.
 *
 * Throws:    <ESL_EMEM> on allocation error.
 *
 * Xref:      squid's ReadStockholm(), 1999.
 */
int
esl_msa_ReadStockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa;
  char      *s;
  int        status;
  int        status2;

  if (ret_msa != NULL) *ret_msa = NULL;
  if (feof(afp->f)) return ESL_EOF;
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
  msa = esl_msa_Create(16, 0);
  if (msa == NULL) return ESL_EMEM;

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
    status = msafile_getline(afp);
    if (status != ESL_OK) /* normal EOF, or thrown EMEM */
      { esl_msa_Destroy(msa); return status; } 
  } while (is_blankline(afp->buf));

  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    { 
      sprintf(afp->errbuf, "missing \"# STOCKHOLM\" header");
      esl_msa_Destroy(msa); 
      return ESL_EFORMAT; 
    } 

  /* Read the alignment file one line at a time.
   */
  while ((status = msafile_getline(afp)) != ESL_OK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */

      if (*s == '#') {
	if      (strncmp(s, "#=GF", 4) == 0 && 
		 (status2 = parse_gf(msa, s)) != ESL_OK)
	  sprintf(afp->errbuf, "failed to parse #=GF line");

	else if (strncmp(s, "#=GS", 4) == 0 &&
		 (status2 = parse_gs(msa, s)) != ESL_OK)
	  sprintf(afp->errbuf, "failed to parse #=GS line");

	else if (strncmp(s, "#=GC", 4) == 0 &&
		 (status2 = parse_gc(msa, s)) != ESL_OK)
	  sprintf(afp->errbuf, "failed to parse #=GC line");

	else if (strncmp(s, "#=GR", 4) == 0 && 
		 (status2 = parse_gr(msa, s)) != ESL_OK)
	  sprintf(afp->errbuf, "failed to parse #=GR line");

	else if ((status2 = parse_comment(msa, s)) != ESL_OK)
	  sprintf(afp->errbuf, "failed to parse comment line");
      } 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n')                   continue;
      else if ((status2 = parse_sequence(msa, s)) != ESL_OK)
	sprintf(afp->errbuf, "failed to parse sequence line");

      if (status2 != ESL_OK) 
	{ esl_msa_Destroy(msa); return status2; }
    }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be ESL_OK.
   */ 
  if (status != ESL_OK)
    { 
      sprintf(afp->errbuf, "didn't find // at end of alignment %.128s",
	      msa->name == NULL ? "" : msa->name);
      esl_msa_Destroy(msa); 
      return ESL_EFORMAT; 
    } 
  
  /* Stockholm's complex, so give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != ESL_OK)
    { esl_msa_Destroy(msa); return ESL_EFORMAT; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return ESL_OK;
}

static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}

/* Format of a GF line:
 *    #=GF <tag> <text>
 * Returns ESL_OK on success; ESL_EFORMAT on parse failure.
 * Throws ESL_EMEM on allocation failure.
 */
static int
parse_gf(ESL_MSA *msa, char *buf)
{
  char *gf;
  char *tag;
  char *text;
  char *s;
  int   n;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gf,   NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, "\n\r",    &text, &n)   != ESL_OK) return ESL_EFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;

  if      (strcmp(tag, "ID") == 0)
    status = esl_strdup(text, n, &(msa->name));
  else if (strcmp(tag, "AC") == 0) 
    status = esl_strdup(text, n, &(msa->acc));
  else if (strcmp(tag, "DE") == 0) 
    status = esl_strdup(text, n, &(msa->desc));
  else if (strcmp(tag, "AU") == 0) 
    status = esl_strdup(text, n, &(msa->au));
  else if (strcmp(tag, "GA") == 0) 
    {				/* Pfam has GA1, GA2. Rfam just has GA1. */
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &text, NULL)) != ESL_OK) 
	return ESL_EFORMAT;
      msa->cutoff[eslMSA_GA1] = atof(text);
      msa->cutset[eslMSA_GA1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &text, NULL)) != ESL_OK) 
	{
	  msa->cutoff[eslMSA_GA2] = atof(text);
	  msa->cutset[eslMSA_GA2] = TRUE;
	}
      status = ESL_OK;
    }
  else if (strcmp(tag, "NC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &text, NULL)) != ESL_OK) 
	return ESL_EFORMAT;
      msa->cutoff[eslMSA_NC1] = atof(text);
      msa->cutset[eslMSA_NC1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &text, NULL)) != ESL_OK) 
	{
	  msa->cutoff[eslMSA_NC2] = atof(text);
	  msa->cutset[eslMSA_NC2] = TRUE;
	}
      status = ESL_OK;
    }
  else if (strcmp(tag, "TC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &text, NULL)) != ESL_OK) 
	return ESL_EFORMAT;
      msa->cutoff[eslMSA_TC1] = atof(text);
      msa->cutset[eslMSA_TC1] = TRUE;
      if ((esl_strtok(&s, "\t\n\r", &text, NULL)) != ESL_OK) 
	{
	  msa->cutoff[eslMSA_TC2] = atof(text);
	  msa->cutset[eslMSA_TC2] = TRUE;
	}
      status = ESL_OK;
    }
  else 				/* an unparsed #=GF: */
    status = add_gf(msa, tag, text);

  return status;
}


/* Format of a GS line:
 *    #=GS <seqname> <tag> <text>
 * Return <ESL_OK> on success; <ESL_EFORMAT> on parse error.
 * Throws <ESL_EMEM> on allocation error (trying to grow for a new
 *        name; <ESL_EINVAL> if we try to grow an ungrowable MSA.
 */
static int
parse_gs(ESL_MSA *msa, char *buf)
{
  char *gs;
  char *seqname;
  char *tag;
  char *text; 
  int   seqidx;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gs,      NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, "\n\r",    &text,    NULL) != ESL_OK) return ESL_EFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;
  
  /* GS usually follows another GS; guess lastidx+1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != ESL_OK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "WT") == 0)
    {
      msa->wgt[seqidx] = atof(text);
      msa->flags      |= eslMSA_HASWGTS;
      status           = ESL_OK;
    }
  else if (strcmp(tag, "AC") == 0)
    status = set_seq_accession(msa, seqidx, text);
  else if (strcmp(tag, "DE") == 0)
    status = set_seq_description(msa, seqidx, text);
  else				
    status = add_gs(msa, tag, seqidx, text);

  return status;
}



/* parse_gc():
 * Format of a GC line:
 *    #=GC <tag> <aligned text>
 */
static int 
parse_gc(ESL_MSA *msa, char *buf)
{
  char *gc;
  char *tag;
  char *text; 
  char *s;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gc,   NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text, &len) != ESL_OK) return ESL_EFORMAT;
  
  if (strcmp(tag, "SS_cons") == 0)
    status = esl_strcat(&(msa->ss_cons), -1, text, len);
  else if (strcmp(tag, "SA_cons") == 0)
    status = esl_strcat(&(msa->sa_cons), -1, text, len);
  else if (strcmp(tag, "RF") == 0)
    status = esl_strcat(&(msa->rf), -1, text, len);
  else
    status = append_gc(msa, tag, text);

  return status;
}

/* parse_gr():
 * Format of a GR line:
 *    #=GR <seqname> <featurename> <text>
 */
static int
parse_gr(ESL_MSA *msa, char *buf)
{
  char *gr;
  char *seqname;
  char *tag;
  char *text;
  int   seqidx;
  int   len;
  int   j;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gr,      NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != ESL_OK) return ESL_EFORMAT;

  /* GR usually follows sequence it refers to; guess msa->lastidx */
  status = get_seqidx(msa, seqname, msa->lastidx, &seqidx);
  if (status != ESL_OK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "SS") == 0) 
    {
      if (msa->ss == NULL)
	{
	  ESL_MALLOC(msa->ss,    sizeof(char *) * msa->sqalloc);
	  ESL_MALLOC(msa->sslen, sizeof(int)    * msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++)
	    {
	      msa->ss[j]    = NULL;
	      msa->sslen[j] = 0;
	    }
	}
      msa->sslen[seqidx] += len;
      status = esl_strcat(&(msa->ss[seqidx]), msa->sslen[seqidx], text, len);
    }
  else if (strcmp(tag, "SA") == 0)
    {
      if (msa->sa == NULL)
	{
	  ESL_MALLOC(msa->sa,    sizeof(char *) * msa->sqalloc);
	  ESL_MALLOC(msa->salen, sizeof(int)    * msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->sa[j]    = NULL;
	      msa->salen[j] = 0;
	    }
	}
      msa->salen[seqidx] += len;
      status = esl_strcat(&(msa->sa[seqidx]), msa->salen[seqidx], text, len);
    }
  else 
    status = append_gr(msa, tag, seqidx, text);

  return status;
}


/* parse_comment():
 * comments are simply stored verbatim, not parsed
 */
static int
parse_comment(ESL_MSA *msa, char *buf)
{
  char *s;
  char *comment;

  s = buf + 1;			               /* skip leading '#' */
  if (*s == '\n') { *s = '\0'; comment = s; }  /* deal with blank comment */
  else if (esl_strtok(&s, "\n\r", &comment, NULL)!= ESL_OK) return ESL_EFORMAT;
  return (add_comment(msa, comment));
}

/* parse_sequence():
 * Format of line is:
 *     <name>  <aligned text>
 */
static int
parse_sequence(ESL_MSA *msa, char *buf)
{
  char *s;
  char *seqname;
  char *text;
  int   seqidx;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != ESL_OK) return ESL_EFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != ESL_OK) return ESL_EFORMAT; 
  
  /* seq usually follows another seq; guess msa->lastidx +1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  msa->lastidx = seqidx;

  msa->sqlen[seqidx] += len;
  return (esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len));
}





/*-------------------- end of Stockholm format section ----------------------*/












/*****************************************************************
 * @LICENSE@
 *****************************************************************/
