/* msa.c
 * Multiple sequence alignment file i/o.
 * 
 * SRE, Thu Jan 20 08:50:43 2005 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <easel/easel.h>
#include <easel/msa.h>



/* Function: esl_msa_Create()
 * Incept:   squid's MSACreate(), 1999.
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
 *                  msa = esl_msa_Create(10, 0);
 *                  if (msa->nseq == msa->sqalloc) esl_msa_Expand(msa);
 *                \end{cchunk}   
 *                
 *           A created <msa> can only be <Expand()>'ed if <alen> is 0
 *           (i.e. case 3). 
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or 0      
 *
 * Returns:  pointer to new MSA object, w/ all values initialized.
 *           Note that msa->nseq is initialized to 0, even though space
 *           is allocated.
 *           
 * Throws:   NULL on allocation failure.          
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
  msa->type    = eslUNKNOWN;	/* no alphabet type yet */
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
   * Stockholm tags require augmentation with the keyhash module.
   */
  msa->comment        = NULL;
  msa->ncomment       = 0;
  msa->alloc_ncomment = 0;
#ifdef ESL_KEYHASH_INCLUDED
  msa->gf_tag         = NULL;
  msa->gf             = NULL;
  msa->ngf            = 0;
  msa->alloc_ngf      = 0;

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
 * Incept:    squid's MSAFree(), 1999.
 *
 * Purpose:   Destroys <msa>.
 */
void
esl_msa_Destroy(ESL_MSA *msa)
{
  int i;

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

#ifdef ESL_KEYHASH_INCLUDED
  esl_Free2D((void **) msa->gf_tag,  msa->ngf);
  esl_Free2D((void **) msa->gf,      msa->ngf);
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
 * Incept:    squid's MSAExpand(), 1999.
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
 */
int
esl_msa_Expand(MSA *msa)
{
  int   old, new;		/* old & new allocation sizes */
  void *p;			/* tmp ptr to realloc'ed memory */
  int   i,j;

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



/* Function:  esl_msa_SetSeqAccession()
 * Incept:    squid's MSASetSeqAccession(), 1999.
 *
 * Purpose:   Sets the sequence accession field for sequence
 *            number <seqidx> in an alignment <msa>, by
 *            copying the string <acc>.
 *
 * Returns:   <ESL_OK> on success.
 * 
 * Throws:    <ESL_EMEM> on allocation failure.
 */
int
esl_msa_SetSeqAccession(MSA *msa, int seqidx, char *acc)
{
  int i;

  /* If this is the first acccession, we have to
   * initialize the whole optional array.
   */
  if (msa->sqacc == NULL) 
    {
      if ((msa->sqacc = malloc(sizeof(char *) * msa->sqalloc)) == NULL)
	ESL_ERROR(ESL_EMEM, "malloc failed");
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqacc[i] = NULL;
    }
  /* If we already had an accession, that's weird, but free it. 
   */
  if (msa->sqacc[seqidx] != NULL) free(msa->sqacc[seqidx]);
  
  msa->sqacc[seqidx] = esl_strdup(acc, -1);
}

/* Function: esl_msa_SetSeqDescription()
 * Incept:   squid's MSASetSeqDescription(), 1999.
 *
 * Purpose:  Set the sequence description field for sequence number
 *           <seqidx> in an alignment <msa> by copying the string <desc>.
 *
 * Returns:  <ESL_OK> on success.
 * 
 * Throws:   <ESL_EMEM> on allocation failure.
 */
int
esl_msa_SetSeqDescription(MSA *msa, int seqidx, char *desc)
{
  int i;

  if (msa->sqdesc == NULL) 
    {
      if ((msa->sqdesc = malloc(sizeof(char *) * msa->sqalloc)) == NULL)
	ESL_ERROR(ESL_EMEM, "malloc failed");
      for (i = 0; i < msa->sqalloc; i++)
	msa->sqdesc[i] = NULL;
  }
  if (msa->sqdesc[i] != NULL) free(msa->sqdesc[i]);
  msa->sqdesc[seqidx] = esl_strdup(desc, -1);
}


/* Function: esl_msa_AddComment()
 * Date:     squid's MSAAddComment(), 1999.
 *
 * Purpose:  Copy an unparsed comment line <s> into the MSA structure <msa>.
 *
 * Returns:  <ESL_OK>.
 * 
 * Throws:   <ESL_EMEM> on allocation failure.
 */
int
esl_msaAddComment(MSA *msa, char *s)
{
  void *p;
  int new;

  /* If this is our first recorded comment, we need to malloc().
   * Note the arbitrary starting size of 16 lines.
   */
  if (msa->comment == NULL) 
    {
      if ((msa->comment = malloc(sizeof(char *) * 16)) == NULL)
	ESL_ERROR(ESL_EMEM, "malloc failed");
      msa->alloc_ncomment = 16;
    }

  /* Check to see if reallocation for comment lines is needed;
   * if so, redouble.
   */
  if (msa->ncomment == msa->alloc_ncomment) 
    {
      new = 2*msa->alloc_ncomment;
      ESL_REALLOC(msa->comment, p, sizeof(char *) * new);
      msa->alloc_ncomment = new;
    }
  msa->comment[msa->ncomment] = esl_strdup(s, -1);
  msa->ncomment++;
  return ESL_OK;
}


/*****************************************************************
 * Deleted (relative to squid's msa):
 *   MSAAddGF()
 *   MSAAddGS()
 *   MSAAppendGC()
 *   MSAGetGC()
 *   MSAAppendGR()
 * Restore when keyhash/GKI augmentation is implemented.
 *****************************************************************/

/* Function: esl_msa_VerifyParse()
 * Date:     squid's MSAVerifyParse(), 1999.
 *
 * Purpose:  Last function called after a multiple alignment parser
 *           thinks it's done. Checks that parse was successful; makes sure
 *           required information is present; makes sure required
 *           information is consistent. Some fields that are
 *           only use during parsing may be freed (sqlen, for
 *           example), and some fields are finalized now (msa::alen
 *           is set, for example).
 *
 * Returns:  <ESL_OK>.
 *           
 * Throws:   <ESL_EFORMAT> if a problem is detected.
 */
int
esl_msa_VerifyParse(MSA *msa)
{
  int idx;

  if (msa->nseq == 0) 
    {
      esl_error(ESL_EFORMAT, __FILE__, __LINE__
		"MSA parse error: no sequences were found for alignment %s",
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
	  esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		    "MSA %s parse error: no sequence for %s",
		    msa->name != NULL ? msa->name : "", msa->sqname[idx]); 
	  return ESL_EFORMAT;
	}

      /* either all weights must be set, or none of them */
      if ((msa->flags & MSA_HASWGTS) && msa->wgt[idx] == -1.0)
	{
	  esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		    "MSA %s parse error: expected a weight for seq %s", 
		    msa->name != NULL ? msa->name : "", msa->sqname[idx]);
	  return ESL_EFORMAT;
	}

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	{
	  esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		    "MSA %s parse error: sequence %s: length %d, expected %d",
		    msa->name != NULL ? msa->name : "",
		    msa->sqname[idx], msa->sqlen[idx], msa->alen);
	  return ESL_EFORMAT;
	}

      /* if individual SS is present, it must have length right too */
      if (msa->ss != NULL &&
	  msa->ss[idx] != NULL && 
	  msa->sslen[idx] != msa->alen) 
	{
	  esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		    "MSA %s parse error: GR SS for %s: length %d, expected %d",
		    msa->name != NULL ? msa->name : "",
		    msa->sqname[idx], msa->sslen[idx], msa->alen);
	  return ESL_EFORMAT;
	}
				/* if SA is present, must have length right */
      if (msa->sa != NULL && 
	  msa->sa[idx] != NULL && 
	  msa->salen[idx] != msa->alen) 
	{
	  esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		    "MSA %s parse error: GR SA for %s: length %d, expected %d",
		    msa->name != NULL ? msa->name : "",
		    msa->sqname[idx], msa->salen[idx], msa->alen);
	  return ESL_EFORMAT;
	}
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    {
      esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		"MSA %s parse error: GC SS_cons markup: len %d, expected %d",
		msa->name != NULL ? msa->name : "",
		strlen(msa->ss_cons), msa->alen);
      return ESL_EFORMAT;
    }

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    {
      esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		"MSA %s parse error: GC SA_cons markup: len %d, expected %d",
		msa->name != NULL ? msa->name : "",
		strlen(msa->sa_cons), msa->alen);
      return ESL_EFORMAT;
    }

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    {
      esl_error(ESL_EFORMAT, __FILE__, __LINE__,
		"MSA %s parse error: GC RF markup: len %d, expected %d",
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

/*SRE_STOPPED_HERE*/


/* Function: MSAFileOpen()
 * Date:     SRE, Tue May 18 13:22:01 1999 [St. Louis]
 *
 * Purpose:  Open an alignment database file and prepare
 *           for reading one alignment, or sequentially
 *           in the (rare) case of multiple MSA databases
 *           (e.g. Stockholm format).
 *           
 * Args:     filename - name of file to open
 *                      if "-", read stdin
 *                      if it ends in ".gz", read from pipe to gunzip -dc
 *           format   - format of file (e.g. MSAFILE_STOCKHOLM)
 *           env      - environment variable for path (e.g. BLASTDB)
 *
 * Returns:  opened MSAFILE * on success.
 *           NULL on failure: 
 *             usually, because the file doesn't exist;
 *             for gzip'ed files, may also mean that gzip isn't in the path.
 */
MSAFILE *
MSAFileOpen(char *filename, int format, char *env)
{
  MSAFILE *afp;
  
  afp        = MallocOrDie(sizeof(MSAFILE));
  if (strcmp(filename, "-") == 0)
    {
      afp->f         = stdin;
      afp->do_stdin  = TRUE; 
      afp->do_gzip   = FALSE;
      afp->fname     = sre_strdup("[STDIN]", -1);
      afp->ssi       = NULL;	/* can't index stdin because we can't seek*/
    }
#ifndef SRE_STRICT_ANSI		
  /* popen(), pclose() aren't portable to non-POSIX systems; disable */
  else if (Strparse("^.*\\.gz$", filename, 0))
    {
      char cmd[256];

      /* Note that popen() will return "successfully"
       * if file doesn't exist, because gzip works fine
       * and prints an error! So we have to check for
       * existence of file ourself.
       */
      if (! FileExists(filename))
	Die("%s: file does not exist", filename);
      if (strlen(filename) + strlen("gzip -dc ") >= 256)
	Die("filename > 255 char in MSAFileOpen()"); 
      sprintf(cmd, "gzip -dc %s", filename);
      if ((afp->f = popen(cmd, "r")) == NULL)
	return NULL;

      afp->do_stdin = FALSE;
      afp->do_gzip  = TRUE;
      afp->fname    = sre_strdup(filename, -1);
      /* we can't index a .gz file, because we can't seek in a pipe afaik */
      afp->ssi      = NULL;	
    }
#endif /*SRE_STRICT_ANSI*/
  else
    {
      char *ssifile;
      char *dir;

      /* When we open a file, it may be either in the current
       * directory, or in the directory indicated by the env
       * argument - and we have to construct the SSI filename accordingly.
       */
      if ((afp->f = fopen(filename, "r")) != NULL)
	{
	  ssifile = MallocOrDie(sizeof(char) * (strlen(filename) + 5));
	  sprintf(ssifile, "%s.ssi", filename);
	}
      else if ((afp->f = EnvFileOpen(filename, env, &dir)) != NULL)
	{
	  char *full;
	  full = FileConcat(dir, filename);
	  ssifile = MallocOrDie(sizeof(char) * (strlen(full) + strlen(filename)  + 5));
	  sprintf(ssifile, "%s.ssi", full);
	  free(dir);
	}
      else return NULL;

      afp->do_stdin = FALSE;
      afp->do_gzip  = FALSE;
      afp->fname    = sre_strdup(filename, -1);
      afp->ssi      = NULL;

      /* Open the SSI index file. If it doesn't exist, or
       * it's corrupt, or some error happens, afp->ssi stays NULL.
       */
      SSIOpen(ssifile, &(afp->ssi));
      free(ssifile);
    }

  /* Invoke autodetection if we haven't already been told what
   * to expect.
   */
  if (format == MSAFILE_UNKNOWN)
    {
      if (afp->do_stdin == TRUE || afp->do_gzip)
	Die("Can't autodetect alignment file format from a stdin or gzip pipe");
      format = MSAFileFormat(afp);
      if (format == MSAFILE_UNKNOWN)
	Die("Can't determine format of multiple alignment file %s", afp->fname);
    }

  afp->format     = format;
  afp->linenumber = 0;
  afp->buf        = NULL;
  afp->buflen     = 0;

  return afp;
}


/* Function: MSAFilePositionByKey()
 *           MSAFilePositionByIndex()
 *           MSAFileRewind()
 * 
 * Date:     SRE, Tue Nov  9 19:02:54 1999 [St. Louis]
 *
 * Purpose:  Family of functions for repositioning in
 *           open MSA files; analogous to a similarly
 *           named function series in HMMER's hmmio.c.
 *
 * Args:     afp    - open alignment file
 *           offset - disk offset in bytes
 *           key    - key to look up in SSI indices 
 *           idx    - index of alignment.
 *
 * Returns:  0 on failure.
 *           1 on success.
 *           If called on a non-fseek()'able file (e.g. a gzip'ed
 *           or pipe'd alignment), returns 0 as a failure flag.
 */
int 
MSAFileRewind(MSAFILE *afp)
{
  if (afp->do_gzip || afp->do_stdin) return 0;
  rewind(afp->f);
  return 1;
}
int 
MSAFilePositionByKey(MSAFILE *afp, char *key)
{
  int       fh;			/* filehandle is ignored       */
  SSIOFFSET offset;		/* offset of the key alignment */

  if (afp->ssi == NULL) return 0;
  if (SSIGetOffsetByName(afp->ssi, key, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(afp->f, &offset) != 0) return 0;
  return 1;
}
int
MSAFilePositionByIndex(MSAFILE *afp, int idx)
{
  int       fh;			/* filehandled is passed but ignored */
  SSIOFFSET offset;		/* disk offset of desired alignment  */

  if (afp->ssi == NULL) return 0;
  if (SSIGetOffsetByNumber(afp->ssi, idx, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(afp->f, &offset) != 0) return 0;
  return 1;
}


/* Function: MSAFileRead()
 * Date:     SRE, Fri May 28 16:01:43 1999 [St. Louis]
 *
 * Purpose:  Read the next msa from an open alignment file.
 *           This is a wrapper around format-specific calls.
 *
 * Args:     afp     - open alignment file
 *
 * Returns:  next alignment, or NULL if out of alignments 
 */
MSA *
MSAFileRead(MSAFILE *afp)
{
  MSA *msa = NULL;

  switch (afp->format) {
  case MSAFILE_STOCKHOLM: msa = ReadStockholm(afp); break;
  case MSAFILE_MSF:       msa = ReadMSF(afp);       break;
  case MSAFILE_A2M:       msa = ReadA2M(afp);       break;
  case MSAFILE_CLUSTAL:   msa = ReadClustal(afp);   break;
  case MSAFILE_SELEX:     msa = ReadSELEX(afp);     break;
  case MSAFILE_PHYLIP:    msa = ReadPhylip(afp);    break;
  default:
    Die("MSAFILE corrupted: bad format index");
  }
  return msa;
}

/* Function: MSAFileClose()
 * Date:     SRE, Tue May 18 14:05:28 1999 [St. Louis]
 *
 * Purpose:  Close an open MSAFILE.
 *
 * Args:     afp  - ptr to an open MSAFILE.
 *
 * Returns:  void
 */
void
MSAFileClose(MSAFILE *afp)
{
#ifndef SRE_STRICT_ANSI	 /* gzip functionality only on POSIX systems */
  if (afp->do_gzip)    pclose(afp->f);
#endif
  if (! afp->do_stdin) fclose(afp->f);
  if (afp->buf  != NULL) free(afp->buf);
  if (afp->ssi  != NULL) SSIClose(afp->ssi);
  if (afp->fname != NULL) free(afp->fname);
  free(afp);
}

char *
MSAFileGetLine(MSAFILE *afp)
{
  char *s;
  if ((s = sre_fgets(&(afp->buf), &(afp->buflen), afp->f)) == NULL)
    return NULL;
  afp->linenumber++;
  return afp->buf;
}

void 
MSAFileWrite(FILE *fp, MSA *msa, int outfmt, int do_oneline)
{
  switch (outfmt) {
  case MSAFILE_A2M:       WriteA2M(fp, msa);     break;
  case MSAFILE_CLUSTAL:   WriteClustal(fp, msa); break;
  case MSAFILE_MSF:       WriteMSF(fp, msa);     break;
  case MSAFILE_PHYLIP:    WritePhylip(fp, msa);  break;
  case MSAFILE_SELEX:     WriteSELEX(fp, msa);   break;
  case MSAFILE_STOCKHOLM:
    if (do_oneline) WriteStockholmOneBlock(fp, msa);
    else            WriteStockholm(fp, msa);
    break;
  default:
    Die("can't write. no such alignment format %d\n", outfmt);
  }
}

/* Function: MSAGetSeqidx()
 * Date:     SRE, Wed May 19 15:08:25 1999 [St. Louis]
 *
 * Purpose:  From a sequence name, return seqidx appropriate
 *           for an MSA structure.
 *           
 *           1) try to guess the index. (pass -1 if you can't guess)
 *           2) Look up name in msa's hashtable.
 *           3) If it's a new name, store in msa's hashtable;
 *                                  expand allocs as needed;
 *                                  save sqname.
 *
 * Args:     msa   - alignment object
 *           name  - a sequence name
 *           guess - a guess at the right index, or -1 if no guess.
 *
 * Returns:  seqidx
 */
int
MSAGetSeqidx(MSA *msa, char *name, int guess)
{
  int seqidx;
				/* can we guess? */
  if (guess >= 0 && guess < msa->nseq && strcmp(name, msa->sqname[guess]) == 0) 
    return guess;
				/* else, a lookup in the index */
  if ((seqidx = GKIKeyIndex(msa->index, name)) >= 0)
    return seqidx;
				/* else, it's a new name */
  seqidx = GKIStoreKey(msa->index, name);
  if (seqidx >= msa->nseqalloc)  MSAExpand(msa);

  msa->sqname[seqidx] = sre_strdup(name, -1);
  msa->nseq++;
  return seqidx;
}






/* Function: MSAFileFormat()
 * Date:     SRE, Fri Jun 18 14:26:49 1999 [Sanger Centre]
 *
 * Purpose:  (Attempt to) determine the format of an alignment file.
 *           Since it rewinds the file pointer when it's done,
 *           cannot be used on a pipe or gzip'ed file. Works by
 *           calling SeqfileFormat() from sqio.c, then making sure
 *           that the format is indeed an alignment. If the format
 *           comes back as FASTA, it assumes that the format as A2M 
 *           (e.g. aligned FASTA).
 *
 * Args:     fname   - file to evaluate
 *
 * Returns:  format code; e.g. MSAFILE_STOCKHOLM
 */
int
MSAFileFormat(MSAFILE *afp)
{
  int fmt;

  fmt = SeqfileFormat(afp->f);

  if (fmt == SQFILE_FASTA) fmt = MSAFILE_A2M;

  if (fmt != MSAFILE_UNKNOWN && ! IsAlignmentFormat(fmt)) 
    Die("File %s does not appear to be an alignment file;\n\
rather, it appears to be an unaligned file in %s format.\n\
I'm expecting an alignment file in this context.\n",
	afp->fname,
	SeqfileFormat2String(fmt));
  return fmt;
}


/* Function: MSAMingap()
 * Date:     SRE, Mon Jun 28 18:57:54 1999 [on jury duty, St. Louis Civil Court]
 *
 * Purpose:  Remove all-gap columns from a multiple sequence alignment
 *           and its associated per-residue data.
 *
 * Args:     msa - the alignment
 *
 * Returns:  (void)
 */
void
MSAMingap(MSA *msa)
{
  int *useme;			/* array of TRUE/FALSE flags for which columns to keep */
  int apos;			/* position in original alignment */
  int idx;			/* sequence index */

  useme = MallocOrDie(sizeof(int) * msa->alen);
  for (apos = 0; apos < msa->alen; apos++)
    {
      for (idx = 0; idx < msa->nseq; idx++)
	if (! isgap(msa->aseq[idx][apos]))
	  break;
      if (idx == msa->nseq) useme[apos] = FALSE; else useme[apos] = TRUE;
    }
  MSAShorterAlignment(msa, useme);
  free(useme);
  return;
}

/* Function: MSANogap()
 * Date:     SRE, Wed Nov 17 09:59:51 1999 [St. Louis]
 *
 * Purpose:  Remove all columns from a multiple sequence alignment that
 *           contain any gaps -- used for filtering before phylogenetic
 *           analysis.
 *
 * Args:     msa - the alignment
 *
 * Returns:  (void). The alignment is modified, so if you want to keep
 *           the original for something, make a copy.
 */
void
MSANogap(MSA *msa)
{
  int *useme;			/* array of TRUE/FALSE flags for which columns to keep */
  int apos;			/* position in original alignment */
  int idx;			/* sequence index */

  useme = MallocOrDie(sizeof(int) * msa->alen);
  for (apos = 0; apos < msa->alen; apos++)
    {
      for (idx = 0; idx < msa->nseq; idx++)
	if (isgap(msa->aseq[idx][apos]))
	  break;
      if (idx == msa->nseq) useme[apos] = TRUE; else useme[apos] = FALSE;
    }
  MSAShorterAlignment(msa, useme);
  free(useme);
  return;
}


/* Function: MSAShorterAlignment()
 * Date:     SRE, Wed Nov 17 09:49:32 1999 [St. Louis]
 *
 * Purpose:  Given an array "useme" (0..alen-1) of TRUE/FALSE flags,
 *           where TRUE means "keep this column in the new alignment":
 *           Remove all columns annotated as "FALSE" in the useme
 *           array.
 *
 * Args:     msa   - the alignment. The alignment is changed, so
 *                   if you don't want the original screwed up, make
 *                   a copy of it first.
 *           useme - TRUE/FALSE flags for columns to keep: 0..alen-1
 *
 * Returns:  (void)
 */
void
MSAShorterAlignment(MSA *msa, int *useme)
{
  int apos;			/* position in original alignment */
  int mpos;			/* position in new alignment      */
  int idx;			/* sequence index */
  int i;			/* markup index */

  /* Since we're minimizing, we can overwrite, using already allocated
   * memory.
   */
  for (apos = 0, mpos = 0; apos < msa->alen; apos++)
    {
      if (useme[apos] == FALSE) continue;

			/* shift alignment and associated per-column+per-residue markup */
      if (mpos != apos)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    {
	      msa->aseq[idx][mpos] = msa->aseq[idx][apos];
	      if (msa->ss != NULL && msa->ss[idx] != NULL) msa->ss[idx][mpos] = msa->ss[idx][apos];
	      if (msa->sa != NULL && msa->sa[idx] != NULL) msa->sa[idx][mpos] = msa->sa[idx][apos];
	      
	      for (i = 0; i < msa->ngr; i++)
		if (msa->gr[i][idx] != NULL) msa->gr[i][idx][mpos] = msa->gr[i][idx][apos];
	    }
	  
	  if (msa->ss_cons != NULL) msa->ss_cons[mpos] = msa->ss_cons[apos];
	  if (msa->sa_cons != NULL) msa->sa_cons[mpos] = msa->sa_cons[apos];
	  if (msa->rf      != NULL) msa->rf[mpos]      = msa->rf[apos];

	  for (i = 0; i < msa->ngc; i++)
	    msa->gc[i][mpos] = msa->gc[i][apos];
	}
      mpos++;
    }
		
  msa->alen = mpos;		/* set new length */
				/* null terminate everything */
  for (idx = 0; idx < msa->nseq; idx++)
    {
      msa->aseq[idx][mpos] = '\0';
      if (msa->ss != NULL && msa->ss[idx] != NULL) msa->ss[idx][mpos] = '\0';
      if (msa->sa != NULL && msa->sa[idx] != NULL) msa->sa[idx][mpos] = '\0';
	      
      for (i = 0; i < msa->ngr; i++)
	if (msa->gr[i][idx] != NULL) msa->gr[i][idx][mpos] = '\0';
    }

  if (msa->ss_cons != NULL) msa->ss_cons[mpos] = '\0';
  if (msa->sa_cons != NULL) msa->sa_cons[mpos] = '\0';
  if (msa->rf != NULL)      msa->rf[mpos] = '\0';

  for (i = 0; i < msa->ngc; i++)
    msa->gc[i][mpos] = '\0';

  return;
}


/* Function: MSASmallerAlignment()
 * Date:     SRE, Wed Jun 30 09:56:08 1999 [St. Louis]
 *
 * Purpose:  Given an array "useme" of TRUE/FALSE flags for
 *           each sequence in an alignment, construct
 *           and return a new alignment containing only 
 *           those sequences that are flagged useme=TRUE.
 *           
 *           Used by routines such as MSAFilterAlignment()
 *           and MSASampleAlignment().
 *           
 * Limitations:
 *           Does not copy unparsed Stockholm markup.
 *
 *           Does not make assumptions about meaning of wgt;
 *           if you want the new wgt vector renormalized, do
 *           it yourself with FNorm(new->wgt, new->nseq). 
 *
 * Args:     msa     -- the original (larger) alignment
 *           useme   -- [0..nseq-1] array of TRUE/FALSE flags; TRUE means include 
 *                      this seq in new alignment
 *           ret_new -- RETURN: new alignment          
 *
 * Returns:  void
 *           ret_new is allocated here; free with MSAFree() 
 */
void
MSASmallerAlignment(MSA *msa, int *useme, MSA **ret_new)
{
  MSA *new;                     /* RETURN: new alignment */
  int nnew;			/* number of seqs in new msa (e.g. # of TRUEs) */
  int oidx, nidx;		/* old, new indices */
  int i;

  nnew = 0;
  for (oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx]) nnew++;
  if (nnew == 0) { *ret_new = NULL; return; }
  
  new  = MSAAlloc(nnew, 0);
  nidx = 0;
  for (oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx])
      {
	new->aseq[nidx]   = sre_strdup(msa->aseq[oidx],   msa->alen);
	new->sqname[nidx] = sre_strdup(msa->sqname[oidx], msa->alen);
	GKIStoreKey(new->index, msa->sqname[oidx]);
	new->wgt[nidx]    = msa->wgt[oidx];
	if (msa->sqacc != NULL)
	  MSASetSeqAccession(new, nidx, msa->sqacc[oidx]);
	if (msa->sqdesc != NULL)
	  MSASetSeqDescription(new, nidx, msa->sqdesc[oidx]);
	if (msa->ss != NULL && msa->ss[oidx] != NULL)
	  {
	    if (new->ss == NULL) new->ss = MallocOrDie(sizeof(char *) * new->nseq);
	    new->ss[nidx] = sre_strdup(msa->ss[oidx], -1);
	  }
	if (msa->sa != NULL && msa->sa[oidx] != NULL)
	  {
	    if (new->sa == NULL) new->sa = MallocOrDie(sizeof(char *) * new->nseq);
	    new->sa[nidx] = sre_strdup(msa->sa[oidx], -1);
	  }
	nidx++;
      }

  new->nseq    = nnew;
  new->alen    = msa->alen; 
  new->flags   = msa->flags;
  new->type    = msa->type;
  new->name    = sre_strdup(msa->name, -1);
  new->desc    = sre_strdup(msa->desc, -1);
  new->acc     = sre_strdup(msa->acc, -1);
  new->au      = sre_strdup(msa->au, -1);
  new->ss_cons = sre_strdup(msa->ss_cons, -1);
  new->sa_cons = sre_strdup(msa->sa_cons, -1);
  new->rf      = sre_strdup(msa->rf, -1);
  for (i = 0; i < MSA_MAXCUTOFFS; i++) {
    new->cutoff[i]        = msa->cutoff[i];
    new->cutoff_is_set[i] = msa->cutoff_is_set[i];
  }
  free(new->sqlen); new->sqlen = NULL;

  MSAMingap(new);
  *ret_new = new;
  return;
}


/*****************************************************************
 * Retrieval routines
 * 
 * Access to MSA structure data is possible through these routines.
 * I'm not doing this because of object oriented design, though
 * it might work in my favor someday.
 * I'm doing this because lots of MSA data is optional, and
 * checking through the chain of possible NULLs is a pain.
 *****************************************************************/

char *
MSAGetSeqAccession(MSA *msa, int idx)
{
  if (msa->sqacc != NULL && msa->sqacc[idx] != NULL)
    return msa->sqacc[idx];
  else
    return NULL;
}
char *
MSAGetSeqDescription(MSA *msa, int idx)
{
  if (msa->sqdesc != NULL && msa->sqdesc[idx] != NULL)
    return msa->sqdesc[idx];
  else
    return NULL;
}
char *
MSAGetSeqSS(MSA *msa, int idx)
{
  if (msa->ss != NULL && msa->ss[idx] != NULL)
    return msa->ss[idx];
  else
    return NULL;
}
char *
MSAGetSeqSA(MSA *msa, int idx)
{
  if (msa->sa != NULL && msa->sa[idx] != NULL)
    return msa->sa[idx];
  else
    return NULL;
}


/*****************************************************************
 * Information routines
 * 
 * Access information about the MSA.
 *****************************************************************/

/* Function: MSAAverageSequenceLength()
 * Date:     SRE, Sat Apr  6 09:41:34 2002 [St. Louis]
 *
 * Purpose:  Return the average length of the (unaligned) sequences
 *           in the MSA.
 *
 * Args:     msa  - the alignment
 *
 * Returns:  average length
 */
float
MSAAverageSequenceLength(MSA *msa)
{
  int   i;
  float avg;
  
  avg = 0.;
  for (i = 0; i < msa->nseq; i++) 
    avg += (float) DealignedLength(msa->aseq[i]);

  if (msa->nseq == 0) return 0.;
  else                return (avg / msa->nseq);
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
