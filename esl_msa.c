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

#include <easel.h>
#ifdef eslAUGMENT_KEYHASH
#include <esl_keyhash.h>
#endif
#include <esl_msa.h>


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

/* Function:  esl_msa_Create()
 * Incept:    SRE, Sun Jan 23 08:25:26 2005 [St. Louis]
 *
 * Purpose:   Creates and initializes an <ESL_MSA> object, and returns a
 *            pointer to it. If you know exactly the dimensions of the
 *            alignment, both <nseq> and <alen>, then <msa =
 *            esl_msa_Create(nseq, alen)> allocates the whole thing at
 *            once.  If you don't know the dimensions of the alignment
 *            (a typical situation if you're parsing an alignment
 *            file), then pass <alen>=0 (your parser must handle allocation
 *            of individual sequences), and pass an <nseq> that
 *            will be used as an initial allocation size; for example,
 *            <msa = esl_msa_Create(16, 0)>. This allocation can be
 *            expanded (by doubling) by calling <esl_msa_Expand()>, 
 *            for example: 
 *             <if (msa->nseq == msa->sqalloc) esl_msa_Expand(msa);>
 *
 *           A created <msa> can only be <_Expand()>'ed if <alen> is 0.
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or 0      
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *            Note that msa->nseq is initialized to 0, even though space
 *            is allocated.
 *           
 * Throws:    NULL on allocation failure.          
 *
 * Xref:      squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_Create(int nseq, int alen)
{
  ESL_MSA *msa;
  int  i;

  if ((msa = malloc(sizeof(ESL_MSA))) == NULL)
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  msa->aseq    = NULL;
  msa->sqname  = NULL;
  msa->wgt     = NULL;
  msa->alen    = alen;		/* if 0, then we're growable. */
  msa->nseq    = 0;
  msa->flags   = 0;
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

  msa->gs_tag         = NULL;
  msa->gs             = NULL;
  msa->ngs            = 0;

  msa->gc_tag         = NULL;
  msa->gc             = NULL;
  msa->ngc            = 0;

  msa->gr_tag         = NULL;
  msa->gr             = NULL;
  msa->ngr            = 0;

#ifdef eslKEYHASH_INCLUDED
  msa->index     = esl_keyhash_Create();
  msa->gs_idx    = NULL;
  msa->gc_idx    = NULL;
  msa->gr_idx    = NULL;
#endif /*keyhash augmentation*/

  /* Allocation, round 2.
   */
  if ((msa->aseq   = malloc(sizeof(char *) * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((msa->sqname = malloc(sizeof(char *) * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((msa->wgt    = malloc(sizeof(float)  * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((msa->sqlen  = malloc(sizeof(int)    * nseq)) == NULL) 
    { esl_msa_Destroy(msa); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }

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
	    ESL_ERROR_NULL(eslEMEM, "malloc failed"); 
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

  esl_Free2D((void **) msa->gs_tag,  msa->ngs);
  esl_Free3D((void ***)msa->gs,      msa->ngs, msa->nseq);
  esl_Free2D((void **) msa->gc_tag,  msa->ngc);
  esl_Free2D((void **) msa->gc,      msa->ngc);
  esl_Free2D((void **) msa->gr_tag,  msa->ngr);
  esl_Free3D((void ***)msa->gr,      msa->ngr, msa->nseq);

#ifdef eslKEYHASH_INCLUDED
  esl_keyhash_Destroy(msa->index);
  esl_keyhash_Destroy(msa->gs_idx);
  esl_keyhash_Destroy(msa->gc_idx);
  esl_keyhash_Destroy(msa->gr_idx);
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
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on reallocation failure; <msa> is undamaged,
 *            and the caller may attempt to recover from the error.
 *            
 *            <eslEINVAL> if <msa> is not growable: its <alen> field
 *            must be 0 to be growable.
 *
 * Xref:      squid's MSAExpand(), 1999.
 */
int
esl_msa_Expand(ESL_MSA *msa)
{
  int   old, new;		/* old & new allocation sizes */
  void *p;			/* tmp ptr to realloc'ed memory */
  int   i,j;

  if (msa->alen > 0) 
    ESL_ERROR(eslEINVAL, "that MSA is not growable");

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

  msa->sqalloc = new;
  return eslOK;
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
 * Returns:  <eslOK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws:   <eslEMEM> if we try to add a name and allocation fails.
 *           <eslEINVAL> if we try to add a name to a non-growable MSA.
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
    { *ret_idx = guess; return eslOK; }

  /* Else look it up - either brute force
   * or, if we're keyhash-augmented, by hashing.
   */
#ifdef eslKEYHASH_INCLUDED                  
  status = esl_gki_Store(msa->index, name, &seqidx);
  if (status == eslEDUP) return status;	/* already stored this name */
  if (status != eslOK)   return status; /* an error. */
#else
  for (seqidx = 0; seqidx < msa->nseq; seqidx++)
    if (strcmp(msa->sqname[seqidx], name) == 0) break;
  if (seqidx < msa->nseq) 
    { *ret_idx = seqidx; return eslOK; }
#endif

  /* If we reach here, then this is a new name that we're
   * adding.
   */
  if (seqidx >= msa->sqalloc &&  
     (status = esl_msa_Expand(msa)) != eslOK)
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
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
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
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
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
 * Returns:  <eslOK> on success.
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
 * Returns eslOK on success. 
 * Throws eslEMEM on allocation failure.
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
  if (status != eslOK) return status;
  status = esl_strdup(value, -1, &(msa->gf[msa->ngf]));
  if (status != eslOK) return status;
  msa->ngf++;

  return eslOK;
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
 * Returns:  <eslOK> on success
 * Throws:   <eslEMEM> on allocation failure
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
#ifdef eslKEYHASH_INCLUDED
      msa->gs_idx = esl_keyhash_Create();
      status = esl_gki_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK) return status;
      ESL_DASSERT1((tagidx == 0));
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
#ifdef eslKEYHASH_INCLUDED
      status = esl_gki_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
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
      if (status != eslOK) return status;
      msa->ngs++;
    }
  
  /* Store the annotation on the sequence.
   * If seq is unannotated, dup the value; if
   * seq already has a GS annotation, cat a \n, then cat the value.
   */
  if (msa->gs[tagidx][sqidx] == NULL)
    {
      status = esl_strdup(value, -1, &(msa->gs[tagidx][sqidx]));
      if (status != eslOK) return status;
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
  return eslOK;
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
 * Returns:  <eslOK> on success
 * 
 * Throws:   <eslEMEM> on allocation failure
 */
static int
append_gc(ESL_MSA *msa, char *tag, char *value)
{
  int   tagidx;
  int   status;
  void *p;

  /* Is this an unparsed tag name that we recognize?
   * If not, handle adding it to index, and reallocating
   * as needed.
   */
  if (msa->gc_tag == NULL)	/* first tag? init w/ malloc  */
    {
#ifdef eslKEYHASH_INCLUDED
      msa->gc_idx = esl_keyhash_Create();
      status = esl_gki_Store(msa->gc_idx, tag, &tagidx);      
      if (status != eslOK) return status;
      ESL_DASSERT1((tagidx == 0));
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
#ifdef eslKEYHASH_INCLUDED
      status = esl_gki_Store(msa->gc_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
#else
      for (tagidx = 0; tagidx < msa->ngc; tagidx++)
	if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocate, in block of one tag at a time
       */
      if (tagidx == msa->ngc)
	{
	  ESL_REALLOC(msa->gc_tag, p, (msa->ngc+1) * sizeof(char **));
	  ESL_REALLOC(msa->gc,     p, (msa->ngc+1) * sizeof(char **));
	  msa->gc[tagidx] = NULL;
	}
    }
  /* new tag? store it.
   */
  if (tagidx == msa->ngc) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gc_tag[tagidx]))) != eslOK)
	return status;
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
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
append_gr(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int tagidx;
  int i;
  int status;

  if (msa->gr_tag == NULL)	/* first tag? init w/ malloc  */
    {
#ifdef eslKEYHASH_INCLUDED
      msa->gr_idx = esl_keyhash_Create();
      status = esl_gki_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK) return status;
      ESL_DASSERT1((tagidx == 0));
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
#ifdef eslKEYHASH_INCLUDED
      status = esl_gki_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
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
      if ((status = esl_strdup(tag, -1, &(msa->gr_tag[tagidx]))) != eslOK)
	return status;
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
 * (<msa->alen> is set, for example). 
 * 
 * <errbuf> is a place to sprintf an informative message about the
 * reason for a parse error. The caller provides an <errbuf>
 * of at least 512 bytes.
 *
 * Returns:  <eslOK>, and errbuf is set to an empty string.
 *           
 * Throws:   <eslEFORMAT> if a problem is detected, and an
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
      return eslEFORMAT;
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
	  return eslEFORMAT;
	}

      /* either all weights must be set, or none of them */
      if ((msa->flags & eslMSA_HASWGTS) && msa->wgt[idx] == -1.0)
	{
	  sprintf(errbuf,
		  "MSA %.128s parse error: expected a weight for seq %.128s", 
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]);
	  return eslEFORMAT;
	}

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	{
	  sprintf(errbuf,
	  "MSA %.128s parse error: sequence %.128s: length %d, expected %d",
		  msa->name != NULL ? msa->name : "",
		  msa->sqname[idx], msa->sqlen[idx], msa->alen);
	  return eslEFORMAT;
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
	  return eslEFORMAT;
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
	  return eslEFORMAT;
	}
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC SS_cons markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->ss_cons), msa->alen);
      return eslEFORMAT;
    }

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC SA_cons markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->sa_cons), msa->alen);
      return eslEFORMAT;
    }

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    {
      sprintf(errbuf,
	      "MSA %.128s parse error: GC RF markup: len %d, expected %d",
	      msa->name != NULL ? msa->name : "",
	      strlen(msa->rf), msa->alen);
      return eslEFORMAT;
    }

  /* If no weights were set, set 'em all to 1.0 */
  if (!(msa->flags & eslMSA_HASWGTS))
    for (idx = 0; idx < msa->nseq; idx++)
      msa->wgt[idx] = 1.0;

  /* Clean up a little from the parser */
  if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
  if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
  if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }

  return eslOK;
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
 *           reading one alignment, or sequentially in the case of 
 *           multiple MSA databases (e.g. Stockholm format); returns
 *           the opened file pointer in <ret_msafp>.
 *          
 *           There are one or two special cases for <filename>. If
 *           <filename> is "-", then the alignment is read from
 *           <stdin>. If <filename> * ends in ".gz", then the file is
 *           assumed to be compressed * by gzip, and it is opened as a
 *           pipe from <gunzip -dc>. (Auto-decompression of gzipp'ed files
 *           is only available on POSIX-compliant systems w/ popen(), when 
 *           <HAVE_POPEN> is defined at compile-time.)
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
 * Returns:  <eslOK> on success, and <ret_msafp> is set to point at
 *           an open <ESL_MSAFILE>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           Returns <eslENOTFOUND> if <filename> cannot be opened,
 *           or <eslEFORMAT> if autodetection is attempted and 
 *           format cannot be determined.
 *           
 * Throws:   <eslEMEM> on allocation failure.
 *           <eslEINVAL> if format autodetection is attempted on 
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
    ESL_ERROR(eslEMEM, "malloc failed");
  afp->f          = NULL;
  afp->fname      = NULL;
  afp->linenumber = 0;
  afp->errbuf[0]  = '\0';
  afp->buf        = NULL;
  afp->buflen     = 0;
  afp->do_gzip    = FALSE;
  afp->do_stdin   = FALSE;
#ifdef eslSSI_INCLUDED		/* SSI augmentation */
  afp->ssi = NULL;
#endif  

  n        = strlen(filename);
  ssifile  = NULL;

  if (strcmp(filename, "-") == 0)
    {
      afp->f         = stdin;
      afp->do_stdin  = TRUE; 
      status         = esl_strdup("[STDIN]", -1, &(afp->fname));
      if (status != eslOK) 
	{ esl_msafile_Close(afp); return eslEMEM; }
    }
#ifdef HAVE_POPEN
  /* popen(), pclose() aren't portable to non-POSIX systems; 
   * disable this section in strict ANSI C mode.
   */
  /* tricky: if n= length of a string s, then
   * s+n-i repositions pointer s at the last i chars
   * of the string.
   */
  else if (n > 3 && strcmp(filename+n-3, ".gz") == 0)
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
	  return eslENOTFOUND; 
	}
      if ((cmd = malloc(sizeof(char) * (n+1+strlen("gzip -dc ")))) == NULL)
	{ 
	  esl_msafile_Close(afp);
	  ESL_ERROR(eslEMEM, "malloc for cmd failed"); 
	}
      sprintf(cmd, "gzip -dc %s", filename);
      if ((afp->f = popen(cmd, "r")) == NULL)
	{ 
	  esl_msafile_Close(afp); 
	  return eslENOTFOUND; 
	}
      status = esl_strdup(filename, n, &(afp->fname));
      afp->do_gzip  = TRUE;
      if (status != eslOK)
	{ esl_msafile_Close(afp); return eslEMEM; }
    }
#endif /*HAVE_POPEN*/
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
      else if (esl_FileEnvOpen(filename, env, &(afp->f), &envfile) == eslOK)
	{
	  esl_FileNewSuffix(envfile, "ssi", &ssifile);
	  free(envfile);
	}
      else 
	{ esl_msafile_Close(afp); return eslENOTFOUND; }

      afp->do_stdin = FALSE;
      afp->do_gzip  = FALSE;
      status = esl_strdup(filename, n, &(afp->fname));
      if (status != eslOK)
	{ esl_msafile_Close(afp); return eslEMEM; }
    }

  /* If augmented by SSI indexing:
   * Open the SSI index file. If it doesn't exist, or
   * it's corrupt, or some error happens, afp->ssi stays NULL.
   */
#ifdef eslSSI_INCLUDED
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
	  ESL_ERROR(eslEINVAL, 
		    "Can't autodetect alignment file fmt in stdin, gzip pipe");
	}

      if (esl_msa_GuessFileFormat(afp) != eslOK)
	{ esl_msafile_Close(afp); return eslEFORMAT; }
    }
  else 
    afp->format     = format;

  if (ret_msafp != NULL) *ret_msafp = afp; else esl_msafile_Close(afp);
  return eslOK;
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

#ifdef HAVE_POPEN /* gzip functionality */
  if (afp->do_gzip)       pclose(afp->f);
#endif
  if (! afp->do_stdin)    fclose(afp->f);
  if (afp->fname != NULL) free(afp->fname);
  if (afp->buf  != NULL)  free(afp->buf);
#ifdef eslSSI_INCLUDED
  if (afp->ssi  != NULL)  SSIClose(afp->ssi); /* FIXME */
#endif
  free(afp);
}

/* msafile_getline():
 * load the next line of <afp> into <afp->buf>. 
 * Returns eslOK on success, eslEOF on normal eof.
 * Throws eslEMEM on alloc failure.
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
 * Functions for general i/o of all formats                                   *
 *****************************************************************************/

/* Function:  esl_msa_Read()
 * Incept:    SRE, Fri Jan 28 08:10:49 2005 [St. Louis]
 *
 * Purpose:   Reads the next MSA from an open MSA file <afp>,
 *            and returns it via <ret_msa>. 
 *
 * Returns:   <eslOK> on success, and <ret_msa> points at the
 *            new MSA object.
 *            <eslEOF> if there are no more alignments in the file.
 *            <eslEFORMAT> if there is a parse error, and <afp->errbuf>
 *            is set to an informative message.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA *msa;
  int      status;

  *ret_msa = NULL;
  
  switch (afp->format) {
  case eslMSAFILE_STOCKHOLM: status = esl_msa_ReadStockholm(afp, &msa); break;
  default:
    ESL_ERROR(eslEINCONCEIVABLE, "no such format");
  }

  *ret_msa = msa;
  return status;
}

/* Function:  esl_msa_Write()
 * Incept:    SRE, Fri Jan 28 09:29:28 2005 [St. Louis]
 *
 * Purpose:   Writes an alignment <msa> to an open stream <fp>,
 *            in format specified by <fmt>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Write(FILE *fp, ESL_MSA *msa, int fmt)
{
  int status;
  switch (fmt) {
  case eslMSAFILE_STOCKHOLM: status = esl_msa_WriteStockholm(fp, msa); break;
  default: 
    ESL_ERROR(eslEINCONCEIVABLE, "no such format");
  } 
  return status;
}



/* Function:  esl_msa_GuessFileFormat()
 * Incept:    SRE, Fri Jan 28 07:29:00 2005 [St. Louis]
 *
 * Purpose:   Attempts to determine the format of an open alignment file
 *            <afp>, for which <afp->format> is <eslMSAFILE_UNKNOWN>. 
 *            If successful, sets <afp->format>.
 *            
 *            Currently a placeholder: it always guesses Stockholm!
 *
 * Returns:   <eslOK> on success, and sets <afp->format>. 
 *            <eslEFORMAT> if format can't be determined.
 *
 * Xref:      squid's MSAFileFormat()
 */
int
esl_msa_GuessFileFormat(ESL_MSAFILE *afp)
{
  /* Placeholder: FIXME: autodetection code goes here.
   */
  afp->format = eslMSAFILE_STOCKHOLM;
  return eslOK;
}



/*-------------------- end of general i/o functions -------------------------*/






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
static int actually_write_stockholm(FILE *fp, ESL_MSA *msa, int cpl);
static int maxwidth(char **s, int n);

/* Function:  esl_msa_ReadStockholm()
 * Incept:    SRE, Sun Jan 23 08:33:32 2005 [St. Louis]
 *
 * Purpose:   Parse the next alignment from an open Stockholm format alignment
 *            file <afp>, leaving the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            Returns <eslEOF> if there are no more alignments in <afp>,
 *            and <ret_msa> is set to NULL.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem, and <ret_msa> is
 *            set to NULL. 
 *
 * Throws:    <eslEMEM> on allocation error.
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
  if (feof(afp->f)) return eslEOF;
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
  msa = esl_msa_Create(16, 0);
  if (msa == NULL) return eslEMEM;

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
    status = msafile_getline(afp);
    if (status != eslOK) /* normal EOF, or thrown EMEM */
      { esl_msa_Destroy(msa); return status; } 
  } while (is_blankline(afp->buf));

  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    { 
      sprintf(afp->errbuf, "missing \"# STOCKHOLM\" header");
      esl_msa_Destroy(msa); 
      return eslEFORMAT; 
    } 

  /* Read the alignment file one line at a time.
   */
  while ((status = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */

      if (*s == '#') {

	if      (strncmp(s, "#=GF", 4) == 0)
	  {
	    if ((status2 = parse_gf(msa, s)) != eslOK)
	      {
		sprintf(afp->errbuf, "failed to parse #=GF line");
		esl_msa_Destroy(msa); 
		return status2; 
	      }
	  }

	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    if ((status2 = parse_gs(msa, s)) != eslOK)
	      {
		sprintf(afp->errbuf, "failed to parse #=GS line");
		esl_msa_Destroy(msa); 
		return status2; 
	      }
	  }

	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if  ((status2 = parse_gc(msa, s)) != eslOK)
	      {
		sprintf(afp->errbuf, "failed to parse #=GC line");
		esl_msa_Destroy(msa); 
		return status2; 
	      }
	  }

	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if ((status2 = parse_gr(msa, s)) != eslOK)
	      {
		sprintf(afp->errbuf, "failed to parse #=GR line");
		esl_msa_Destroy(msa); 
		return status2; 
	      }
	  }

	else if ((status2 = parse_comment(msa, s)) != eslOK)
	  {
	    sprintf(afp->errbuf, "failed to parse comment line");
	    esl_msa_Destroy(msa); 
	    return status2; 
	  }
      } 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n')                   continue;
      else if ((status2 = parse_sequence(msa, s)) != eslOK)
	{
	  sprintf(afp->errbuf, "failed to parse sequence line");
	  esl_msa_Destroy(msa); 
	  return status2; 
	}
    }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status != eslOK)
    { 
      sprintf(afp->errbuf, "didn't find // at end of alignment %.128s",
	      msa->name == NULL ? "" : msa->name);
      esl_msa_Destroy(msa); 
      return eslEFORMAT; 
    } 
  
  /* Stockholm's complex, so give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != eslOK)
    { esl_msa_Destroy(msa); return eslEFORMAT; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;
}

/* Function:  esl_msa_WriteStockholm()
 * Incept:    SRE, Fri Jan 28 09:24:02 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in multiblock format, with
 *            50 residues per line.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholm(), 1999.
 */
int
esl_msa_WriteStockholm(FILE *fp, ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, 50)); /* 50 char per block */
}

/* Function:  esl_msa_WriteStockholmOneBlock()
 * Incept:    SRE, Fri Jan 28 09:25:42 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in single block (Pfam) format.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholmOneBlock(), 1999.
 */
int
esl_msa_WriteStockholmOneBlock(FILE *fp, ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, msa->alen)); /* one big block */
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
 * Returns eslOK on success; eslEFORMAT on parse failure.
 * Throws eslEMEM on allocation failure.
 */
static int
parse_gf(ESL_MSA *msa, char *buf)
{
  char *gf;
  char *tag;
  char *text;
  char *tok;
  char *s;
  int   n;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gf,   NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text, &n)   != eslOK) return eslEFORMAT;
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
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_GA1] = atof(tok);
      msa->cutset[eslMSA_GA1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_GA2] = atof(tok);
	  msa->cutset[eslMSA_GA2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "NC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_NC1] = atof(tok);
      msa->cutset[eslMSA_NC1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_NC2] = atof(tok);
	  msa->cutset[eslMSA_NC2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "TC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok, NULL)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_TC1] = atof(tok);
      msa->cutset[eslMSA_TC1] = TRUE;
      if ((esl_strtok(&s, "\t\n\r", &tok, NULL)) == eslOK) 
	{
	  msa->cutoff[eslMSA_TC2] = atof(tok);
	  msa->cutset[eslMSA_TC2] = TRUE;
	}
      status = eslOK;
    }
  else 				/* an unparsed #=GF: */
    status = add_gf(msa, tag, text);

  return status;
}


/* Format of a GS line:
 *    #=GS <seqname> <tag> <text>
 * Return <eslOK> on success; <eslEFORMAT> on parse error.
 * Throws <eslEMEM> on allocation error (trying to grow for a new
 *        name; <eslEINVAL> if we try to grow an ungrowable MSA.
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
  if (esl_strtok(&s, " \t\n\r", &gs,      NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text,    NULL) != eslOK) return eslEFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;
  
  /* GS usually follows another GS; guess lastidx+1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "WT") == 0)
    {
      msa->wgt[seqidx] = atof(text);
      msa->flags      |= eslMSA_HASWGTS;
      status           = eslOK;
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
  if (esl_strtok(&s, " \t\n\r", &gc,   NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,  NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text, &len) != eslOK) return eslEFORMAT;
  
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
  if (esl_strtok(&s, " \t\n\r", &gr,      NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag,     NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != eslOK) return eslEFORMAT;

  /* GR usually follows sequence it refers to; guess msa->lastidx */
  status = get_seqidx(msa, seqname, msa->lastidx, &seqidx);
  if (status != eslOK) return status;
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
      status = esl_strcat(&(msa->ss[seqidx]), msa->sslen[seqidx], text, len);
      msa->sslen[seqidx] += len;
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
      status = esl_strcat(&(msa->sa[seqidx]), msa->salen[seqidx], text, len);
      msa->salen[seqidx] += len;
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
  else if (esl_strtok(&s, "\n\r", &comment, NULL)!= eslOK) return eslEFORMAT;
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
  if (esl_strtok(&s, " \t\n\r", &seqname, NULL) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &text,    &len) != eslOK) return eslEFORMAT; 
  
  /* seq usually follows another seq; guess msa->lastidx +1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
  msa->sqlen[seqidx] += len;
  return status;

}


/* actually_write_stockholm()
 * SRE, Fri May 21 17:39:22 1999 [St. Louis]
 *
 * Write an alignment in Stockholm format to an open file. This is the
 * function that actually does the work. The API's WriteStockholm()
 * and WriteStockholmOneBlock() are wrappers.
 *
 * Args:     fp    - file that's open for writing
 *           msa   - alignment to write        
 *           cpl   - characters to write per line in alignment block
 *
 * Returns:  eslOK on success.
 * 
 * Throws:   eslEMEM on allocation failure.
 */
static int
actually_write_stockholm(FILE *fp, ESL_MSA *msa, int cpl)
{
  int  i, j;
  int  maxname;		/* maximum name length     */
  int  maxgf;		/* max #=GF tag length     */
  int  maxgc;		/* max #=GC tag length     */
  int  maxgr; 		/* max #=GR tag length     */
  int  margin;        	/* total left margin width */
  int  gslen;		/* length of a #=GS tag    */
  char *buf;
  int  currpos;
  char *s, *tok;
  
  /* Figure out how much space we need for name + markup
   * to keep the alignment in register. Required by Stockholm
   * spec, even though our Stockholm parser doesn't care (Erik's does).
   *
   * The left margin of an alignment block can be composed of:
   * 
   * <seqname>                      max length: maxname + 1
   * #=GC <gc_tag>                  max length: 4 + 1 + maxgc + 1
   * #=GR <seqname> <gr_tag>        max length: 4 + 1 + maxname + 1 + maxgr + 1
   * 
   * <margin> is the max of these. It is the total length of the
   * left margin that we need to leave, inclusive of the last space.
   * 
   * Then when we output, we do:
   * name:  <leftmargin-1>
   * gc:    #=GC <leftmargin-6>
   * gr:    #=GR <maxname> <leftmargin-maxname-7>
   *
   * xref STL9/p17
   */
  maxname = maxwidth(msa->sqname, msa->nseq);
  
  maxgf   = maxwidth(msa->gf_tag, msa->ngf);
  if (maxgf < 2) maxgf = 2;

  maxgc   = maxwidth(msa->gc_tag, msa->ngc);
  if (msa->rf      !=NULL && maxgc < 2) maxgc = 2;
  if (msa->ss_cons !=NULL && maxgc < 7) maxgc = 7;
  if (msa->sa_cons !=NULL && maxgc < 7) maxgc = 7;

  maxgr   = maxwidth(msa->gr_tag, msa->ngr);
  if (msa->ss != NULL && maxgr < 2) maxgr = 2;
  if (msa->sa != NULL && maxgr < 2) maxgr = 2;

  margin = maxname + 1;
  if (maxgc > 0 && maxgc+6 > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  
  /* Allocate a tmp buffer to hold sequence chunks in
   */
  ESL_MALLOC(buf, sizeof(char) * (cpl+1));

  /* Magic Stockholm header
   */
  fprintf(fp, "# STOCKHOLM 1.0\n");

  /* Free text comments
   */
  for (i = 0;  i < msa->ncomment; i++)
    fprintf(fp, "# %s\n", msa->comment[i]);
  if (msa->ncomment > 0) fprintf(fp, "\n");

  /* GF section: per-file annotation
   */
  if (msa->name != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "ID", msa->name);
  if (msa->acc  != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AC", msa->acc);
  if (msa->desc != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "DE", msa->desc);
  if (msa->au   != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AU", msa->au);
  
  /* Thresholds are hacky. Pfam has two. Rfam has one.
   */
  if      (msa->cutset[eslMSA_GA1] && msa->cutset[eslMSA_GA2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1], msa->cutoff[eslMSA_GA2]);
  else if (msa->cutset[eslMSA_GA1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1]);

  if      (msa->cutset[eslMSA_NC1] && msa->cutset[eslMSA_NC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1], msa->cutoff[eslMSA_NC2]);
  else if (msa->cutset[eslMSA_NC1])
    fprintf(fp, "#=GF %-*s %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1]);

  if      (msa->cutset[eslMSA_TC1] && msa->cutset[eslMSA_TC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "TC", msa->cutoff[eslMSA_TC1], msa->cutoff[eslMSA_TC2]);
  else if (msa->cutset[eslMSA_TC1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "TC", msa->cutoff[eslMSA_TC1]);

  for (i = 0; i < msa->ngf; i++)
    fprintf(fp, "#=GF %-*s %s\n", maxgf, msa->gf_tag[i], msa->gf[i]); 
  fprintf(fp, "\n");


  /* GS section: per-sequence annotation
   */
  if (msa->flags & eslMSA_HASWGTS) 
    {
      for (i = 0; i < msa->nseq; i++) 
	fprintf(fp, "#=GS %-*s WT %.2f\n", 
		maxname, msa->sqname[i], msa->wgt[i]);		
      fprintf(fp, "\n");
    }

  if (msa->sqacc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqacc[i] != NULL)
	  fprintf(fp, "#=GS %-*s AC %s\n", 
		  maxname, msa->sqname[i], msa->sqacc[i]);
      fprintf(fp, "\n");
    }

  if (msa->sqdesc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqdesc[i] != NULL)
	  fprintf(fp, "#=GS %-*s DE %s\n", 
		  maxname, msa->sqname[i], msa->sqdesc[i]);
      fprintf(fp, "\n");
    }

  for (i = 0; i < msa->ngs; i++)
    {
      /* Multiannotated GS tags are possible; for example, 
       *     #=GS foo DR PDB; 1xxx;
       *     #=GS foo DR PDB; 2yyy;
       * These are stored, for example, as:
       *     msa->gs[0][0] = "PDB; 1xxx;\nPDB; 2yyy;"
       * and must be decomposed.
       */
      gslen = strlen(msa->gs_tag[i]);
      for (j = 0; j < msa->nseq; j++)
	if (msa->gs[i][j] != NULL)
	  {
	    s = msa->gs[i][j];
	    while (esl_strtok(&s, "\n", &tok, NULL) == eslOK)
	      fprintf(fp, "#=GS %-*s %-*s %s\n", 
		      maxname, msa->sqname[j],
		      gslen,   msa->gs_tag[i], 
		      tok);
	  }
      fprintf(fp, "\n");
    }

  /* Alignment section:
   * contains aligned sequence, #=GR annotation, and #=GC annotation
   */
  for (currpos = 0; currpos < msa->alen; currpos += cpl)
    {
      if (currpos > 0) fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
	  strncpy(buf, msa->aseq[i] + currpos, cpl);
	  buf[cpl] = '\0';	      
	  fprintf(fp, "%-*s %s\n", 
		  margin-1, msa->sqname[i], buf);

	  if (msa->ss != NULL && msa->ss[i] != NULL) {
	    strncpy(buf, msa->ss[i] + currpos, cpl);
	    buf[cpl] = '\0';	 
	    fprintf(fp, "#=GR %-*s %-*s %s\n", 
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SS",
		    buf);
	  }
	  if (msa->sa != NULL && msa->sa[i] != NULL) {
	    strncpy(buf, msa->sa[i] + currpos, cpl);
	    buf[cpl] = '\0';
	    fprintf(fp, "#=GR %-*s %-*s %s\n",
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SA",
		    buf);
	  }
	  for (j = 0; j < msa->ngr; j++)
	    if (msa->gr[j][i] != NULL) {
	      strncpy(buf, msa->gr[j][i] + currpos, cpl);
	      buf[cpl] = '\0';
	      fprintf(fp, "#=GR %-*s %-*s %s\n", 
		      maxname,          msa->sqname[i],
		      margin-maxname-7, msa->gr_tag[j],
		      buf);
	    }
	}
      if (msa->ss_cons != NULL) {
	strncpy(buf, msa->ss_cons + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SS_cons", buf);
      }

      if (msa->sa_cons != NULL) {
	strncpy(buf, msa->sa_cons + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SA_cons", buf);
      }

      if (msa->rf != NULL) {
	strncpy(buf, msa->rf + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "RF", buf);
      }
      for (j = 0; j < msa->ngc; j++) {
	strncpy(buf, msa->gc[j] + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, msa->gc_tag[j], buf);
      }
    }
  fprintf(fp, "//\n");
  free(buf);
  return eslOK;
}

/* maxwidth()
 * Return the length of the longest string in 
 * an array of strings.
 */
static int
maxwidth(char **s, int n)
{
  int max, i, len;
  
  max = 0;
  for (i = 0; i < n; i++)
    if (s[i] != NULL)
      {
	len = strlen(s[i]);
	if (len > max) max = len;
      }
  return max;
}


/*-------------------- end of Stockholm format section ----------------------*/



/******************************************************************************
 * Example and test driver
 *****************************************************************************/

#ifdef eslMSA_EXAMPLE
/* gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE msa.c easel.c 
 * time ./example SSU_rRNA_5
 * 
 * or add -DeslAUGMENT_KEYHASH, and
 * gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE -DeslAUGMENT_KEYHASH msa.c easel.c keyhash.c
 */
#include <stdio.h>

#include <easel.h>
#ifdef eslAUGMENT_KEYHASH
#include <esl_keyhash.h>
#endif
#include <esl_msa.h>

int
main(int argc, char **argv)
{
  char        *filename;
  int          fmt;
  ESL_MSAFILE *afp;
  ESL_MSA     *msa;
  int          status;
  int          nali;

  filename = argv[1];
  fmt      = eslMSAFILE_UNKNOWN;

  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND) 
    {
      fprintf(stderr, "Alignment file %s not readable\n", filename);
      exit(1);
    } 
  else if (status == eslEFORMAT) 
    {
      fprintf(stderr, "Couldn't determine format of alignment %s\n", filename);
      exit(1);
    } 
  else if (status != eslOK) 
    {
      fprintf(stderr, "Alignment file open failed with error %d\n", status);
      exit(1);
    }

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, msa->nseq, msa->alen);
      esl_msa_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
      esl_msa_Destroy(msa);
    }

  if (status == eslEFORMAT) 
    {
      fprintf(stderr, "Alignment parse error at line %d of file %s:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf);
      fprintf(stderr, "Offending line is: %s\n", afp->buf);
      esl_msafile_Close(afp);
      exit(1);
    } 
  else if (status != eslEOF) 
    {
      fprintf(stderr, "Alignment file read failed with error %d\n", status);
      esl_msafile_Close(afp);
      exit(1);
    }

  esl_msafile_Close(afp);
  exit(0);
}
#endif /*eslMSA_EXAMPLE*/
 
#ifdef eslMSA_TESTDRIVE
/* gcc -g -Wall -o test -I. -DeslMSA_TESTDRIVE msa.c easel.c 
 * ./test
 * 
 * or add -DeslAUGMENT_KEYHASH, and
 * gcc -g -Wall -o test -I. -DeslMSA_TESTDRIVE -DeslAUGMENT_KEYHASH msa.c easel.c keyhash.c
 */
#include <stdlib.h>
#include <stdio.h>

#include <easel.h>
#ifdef eslAUGMENT_KEYHASH
#include <esl_keyhash.h>
#endif
#include <esl_msa.h>

int
main(int argc, char **argv)
{
  char         filename[] = "tmpxxx.ali";
  int          fmt;
  FILE        *fp;
  ESL_MSAFILE *afp;
  ESL_MSA     *msa;
  int          status;

  /* Create a test alignment. 
   * Extensive format testing will rely on external example files;
   * this is just going to be a quickie test that nothing's grossly
   * wrong.
   */
  if ((fp = fopen(filename, "w")) == NULL) abort();
  fprintf(fp, "# STOCKHOLM 1.0\n");
  fprintf(fp, "seq1 ACDEFGHIKLMNPQRSTVWY\n");
  fprintf(fp, "seq2 ACDEFGHIKLMNPQRSTVWY\n");
  fprintf(fp, "\n");
  fprintf(fp, "seq1 ACDEFGHIKLMNPQRSTVWY\n");
  fprintf(fp, "seq2 ACDEFGHIKLMNPQRSTVWY\n");
  fprintf(fp, "//\n");
  fclose(fp);

  /* Read it back in
   */
  fmt = eslMSAFILE_UNKNOWN;
  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status != eslOK) abort();

  status = esl_msa_Read(afp, &msa);
  if (status != eslOK) abort();

  /* Check that it's ok.
   */
  if (msa->alen != 40) abort();
  if (msa->nseq != 2)  abort();
  if (strcmp(msa->aseq[0], "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY") != 0) abort();
  if (strcmp(msa->aseq[1], "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY") != 0) abort();


  /* Try to read one more; file should be empty, so we get EOF.
   */
  status = esl_msa_Read(afp, &msa); 
  if (status != eslEOF) abort();
  esl_msafile_Close(afp);
  exit(0);
}
#endif /*eslMSA_TESTDRIVE*/

/*-------------------- end of test drivers and examples ---------------------*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
