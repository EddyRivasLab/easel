/* esl_dsqdata :: high performance sequence input
 * 
 * Contents:
 *   1. ESL_DSQDATA, high performance sequence data input
 *   2. ESL_DSQDATA_CHUNK, a chunk of input sequence data
 *   3. Loader and unpacker, the input threads
 */

#include "easel.h"
#include "esl_dsqdata.h"

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

static ESL_DSQDATA_CHUNK *dsqdata_chunk_Create(void);
static void               dsqdata_chunk_Destroy(ESL_DSQDATA_CHUNK *chu);

static void *dsqdata_loader_thread  (void *p);
static void *dsqdata_unpacker_thread(void *p);


/*****************************************************************
 * 1. ESL_DSQDATA: high performance sequence data input
 *****************************************************************/

/* Function:  esl_dsqdata_Open()
 * Synopsis:  Open a digital sequence database for reading
 * Incept:    SRE, Wed Jan 20 09:50:00 2016 [Amtrak 2150, NYP-BOS]
 *
 * Purpose:   Open digital sequence database <basename> for reading.
 *            Configure it for a specified number <nconsumers> of
 *            parallelized consumers. The consumers are one or more
 *            threads that are processing chunks of data in parallel.
 *            
 *            The file <basename> is a human-readable stub describing
 *            the database. The actual data are in three accompanying
 *            binary files: the index file <basename>.dsqi, the
 *            metadata file <basename>.dsqm, and the sequence file
 *            <basename>.dsqs.
 *            
 *            <byp_abc> provides a way to either tell <dsqdata> to
 *            expect a specific alphabet in the <basename> database
 *            (and return a normal failure on a mismatch), or, when
 *            the alphabet remains unknown, to figure out the alphabet
 *            in <basename> is and allocate and return a new alphabet.
 *            <byp_abc> uses a partial Easel "bypass" idiom for this:
 *            if <*byp_abc> is NULL, we allocate and return a new
 *            alphabet; if <*byp_abc> is a ptr to an existing
 *            alphabet, we use it for validation. That is,
 *                
 *                abc = NULL;
 *                esl_dsqdata_Open(&abc, basename...)
 *                // <abc> is now the alphabet of <basename>; you're responsible for Destroy'ing it
 *                
 *            or:
 *                abc = esl_alphabet_Create(eslAMINO);
 *                status = esl_dsqdata_Open(&abc, basename);
 *                // if status == eslEINCOMPAT, alphabet in basename doesn't match caller's expectation
 *
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEINCOMPAT> if caller provides a digital alphabet in
 *            <*byp_abc> and it doesn't match the database's alphabet.
 *            
 *            On any normal error, <*ret_dd> is returned but in an
 *            error state, and <dd->errbuf> is a user-directed error
 *            message that the caller can relay.
 *
 * Throws:    <eslEMEM> on allocation error.
 * 
 *            On any thrown exception, <*ret_dd> is returned NULL.
 */
int
esl_dsqdata_Open(ESL_ALPHABET **byp_abc, char *basename, int nconsumers, ESL_DSQDATA **ret_dd)
{
  ESL_DSQDATA *dd = NULL;
  int          status;
  
  ESL_DASSERT1(( nconsumers > 0 ));

  ESL_ALLOC(dd, sizeof(ESL_DSQDATA));
  dd->abc_r           = *byp_abc;        // This may be NULL, in which case we will create it later.
  dd->at_eof          = FALSE;
  dd->sfp             = NULL;
  dd->ifp             = NULL;
  dd->mfp             = NULL;
  dd->loader_outbox   = NULL;
  dd->unpacker_outbox = NULL;
  dd->recycling       = NULL;
  dd->errbuf[0]       = '\0';
  dd->nconsumers      = nconsumers;

  ESL_ALLOC( dd->basename, sizeof(char) * (strlen(basename) + 6)); // +5 for .dsqx; +1 for \0
  sprintf(dd->basename, "%s.dsqi", basename);
  dd->ifp = fopen(dd->basename, "rb");
  sprintf(dd->basename, "%s.dsqm", basename);
  dd->mfp = fopen(dd->basename, "rb");
  sprintf(dd->basename, "%s.dsqs", basename);
  dd->sfp = fopen(dd->basename, "rb");
  strcpy(dd->basename, basename);
  /* TODO: error checking
   *       include hash or random number to verify that files belong together
   *       include binary magic numbers for verification and byteswap detection
   *       if index file ever includes a header, load it here
   */
  

  pthread_mutex_init(&dd->loader_outbox_mutex,   NULL);                dd->lom_c = TRUE;
  pthread_mutex_init(&dd->unpacker_outbox_mutex, NULL);                dd->uom_c = TRUE;
  pthread_mutex_init(&dd->recycling_mutex,       NULL);                dd->rm_c  = TRUE;

  pthread_cond_init(&dd->loader_outbox_full_cv,     NULL);             dd->lof_c = TRUE;
  pthread_cond_init(&dd->loader_outbox_empty_cv,    NULL);             dd->loe_c = TRUE;
  pthread_cond_init(&dd->unpacker_outbox_full_cv,   NULL);             dd->uof_c = TRUE;
  pthread_cond_init(&dd->unpacker_outbox_empty_cv,  NULL);             dd->uoe_c = TRUE;
  pthread_cond_init(&dd->recycling_cv,              NULL);             dd->r_c   = TRUE;
  
  pthread_create(&dd->unpacker_t, NULL, dsqdata_unpacker_thread, dd);  dd->ut_c = TRUE;
  pthread_create(&dd->loader_t,   NULL, dsqdata_loader_thread,   dd);  dd->lt_c = TRUE;

  *ret_dd  = dd;
  *byp_abc = dd->abc_r;     // If caller provided <*byp_abc> this is a no-op, because we set abc_r = *byp_abc.
  return eslOK;

 ERROR:
  esl_dsqdata_Close(dd);
  return status;
}


/* Function:  esl_dsqdata_Read()
 * Synopsis:  Read next chunk of sequence data.
 * Incept:    SRE, Thu Jan 21 11:21:38 2016 [Harvard]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   <eslOK> on success. <*ret_chu> is a chunk of seq data.
 *            Caller needs to call esl_dsqdata_Recycle() on each chunk
 *            that it Read()'s.
 *             
 *            <eslEOF> if we've reached the end of the input file;
 *            <*ret_chu> is NULL.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_dsqdata_Read(ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK **ret_chu)
{
  ESL_DSQDATA_CHUNK *chu = NULL;

  /* The loader and unpacker have already done the work.  All that
   * _Read() needs to do is take a finished chunk from the unpacker's
   * outbox.  That finished chunk could be a final empty chunk, which
   * is the EOF signal.
   */
  
  /* If one reader has already processed eof, all subsequent Read() calls also return eslEOF */
  if (dd->at_eof) { *ret_chu = NULL; return eslEOF; }

  /* Get next chunk from unpacker. Wait if needed. */
  pthread_mutex_lock(&dd->unpacker_outbox_mutex);
  while (dd->unpacker_outbox == NULL)
    pthread_cond_wait(&dd->unpacker_outbox_full_cv, &dd->unpacker_outbox_mutex);
  chu = dd->unpacker_outbox;
  dd->unpacker_outbox = NULL;
  if (! chu->N) dd->at_eof = TRUE;        // The eof flag makes sure only one reader processes EOF chunk.
  pthread_mutex_unlock(&dd->unpacker_outbox_mutex);
  pthread_cond_signal (&dd->unpacker_outbox_empty_cv);
  
  /* If chunk has any data in it, go ahead and return it. */
  if (chu->N)
    {
      *ret_chu = chu;
      return eslOK;
    }
  /* Otherwise, an empty chunk is a signal that the loader and unpacker
   * are done. But the loader is responsible for freeing all the chunks
   * it allocated, so we have to get this chunk back to the loader, via
   * the recycling. (Alternatively, we could let the caller recycle 
   * the chunk on EOF, but letting the caller detect EOF on read and 
   * exit its loop, only recycling chunks inside the loop, is consistent
   * with all the rest of Easel's read idioms.
   */
  else
    {
      esl_dsqdata_Recycle(dd, chu);
      *ret_chu = NULL;
      return eslEOF;
    }
}


int
esl_dsqdata_Recycle(ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK *chu)
{
  pthread_mutex_lock(&dd->recycling_mutex);
  chu->nxt = dd->recycling;                    // Push chunk onto head of recycling stack
  dd->recycling = chu;
  pthread_mutex_unlock(&dd->recycling_mutex);
  pthread_cond_signal(&dd->recycling_cv);      // Tell loader that there's a chunk it can recycle
  return eslOK;
}



void
esl_dsqdata_Close(ESL_DSQDATA *dd)
{
  if (dd)
    {
      if (dd->lt_c)  pthread_join(dd->loader_t,   NULL);
      if (dd->ut_c)  pthread_join(dd->unpacker_t, NULL);
      if (dd->lof_c) pthread_cond_destroy(&dd->loader_outbox_full_cv);
      if (dd->loe_c) pthread_cond_destroy(&dd->loader_outbox_empty_cv);
      if (dd->uof_c) pthread_cond_destroy(&dd->unpacker_outbox_full_cv);
      if (dd->uoe_c) pthread_cond_destroy(&dd->unpacker_outbox_empty_cv);
      if (dd->r_c)   pthread_cond_destroy(&dd->recycling_cv);
      if (dd->lom_c) pthread_mutex_destroy(&dd->loader_outbox_mutex);
      if (dd->uom_c) pthread_mutex_destroy(&dd->unpacker_outbox_mutex);
      if (dd->rm_c)  pthread_mutex_destroy(&dd->recycling_mutex);

      if (dd->ifp) fclose(dd->ifp);
      if (dd->sfp) fclose(dd->sfp);
      if (dd->mfp) fclose(dd->mfp);

      if (dd->basename) free(dd->basename);

      /* Loader thread is responsible for freeing all chunks it created, even on error. */
      ESL_DASSERT1(( dd->loader_outbox   == NULL ));
      ESL_DASSERT1(( dd->unpacker_outbox == NULL ));
      ESL_DASSERT1(( dd->recycling       == NULL ));

      free(dd);
    }
}


/*****************************************************************
 * 2. ESL_DSQDATA_CHUNK: a chunk of input sequence data
 *****************************************************************/

static ESL_DSQDATA_CHUNK *
dsqdata_chunk_Create(void)
{
  ESL_DSQDATA_CHUNK *chu = NULL;
  int                U;               // max size of unpacked seq data, in bytes (smem allocation)
  int                status;

  ESL_ALLOC(chu, sizeof(ESL_DSQDATA_CHUNK));
  chu->i0       = 0;
  chu->N        = 0;
  chu->pn       = 0;
  chu->dsq      = NULL;
  chu->name     = NULL;
  chu->acc      = NULL;
  chu->desc     = NULL;
  chu->taxid    = NULL;
  chu->L        = NULL;
  chu->metadata = NULL;
  chu->smem     = NULL;
  chu->nxt      = NULL;

  /* dsq, name, acc, desc are arrays of pointers into smem, metadata.
   * taxid is cast to int, from the metadata.
   * L is figured out by the unpacker.
   * All of these are set by the unpacker.
   */
  ESL_ALLOC(chu->dsq,   eslDSQDATA_CHUNK_MAXSEQ * sizeof(ESL_DSQ *));   
  ESL_ALLOC(chu->name,  eslDSQDATA_CHUNK_MAXSEQ * sizeof(char *));
  ESL_ALLOC(chu->acc,   eslDSQDATA_CHUNK_MAXSEQ * sizeof(char *));
  ESL_ALLOC(chu->desc,  eslDSQDATA_CHUNK_MAXSEQ * sizeof(char *));
  ESL_ALLOC(chu->taxid, eslDSQDATA_CHUNK_MAXSEQ * sizeof(int));
  ESL_ALLOC(chu->L,     eslDSQDATA_CHUNK_MAXSEQ * sizeof(int64_t));

  /* On the <smem> allocation, and the <dsq> and <psq> pointers into it:
   *
   * _MAX (in uint32's) sets the maximum single fread() size: one load
   * of a new chunk of packed sequence, up to _MAX*4 bytes. <smem>
   * needs to be able to hold the fully unpacked sequence, because we
   * unpack in place. For protein sequence, each uint32 unpacks to at
   * most 6 residues (5-bit packing). We don't pack sentinels, so the
   * maximum unpacked size includes _MAXSEQ+1 sentinels... because we
   * concat the digital seqs so that the trailing sentinel of seq i is
   * the leading sentinel of seq i+1.
   *
   * The packed seq (max of P bytes) loads overlap with the unpacked
   * data (max of U bytes):
   *                   psq
   *                   v[    P bytes    ]
   * smem: 0........0........0..........0
   *       ^[         U bytes           ]
   *       ^dsq[0]  ^dsq[1]  ^dsq[2]
   *
   * and as long as we unpack psq left to right -- and as long as we
   * read the last uint32 before we write the last unpacked residues
   * to smem - we're guaranteed that the unpacking works without
   * overwriting any unpacked data.
   */
  U =  6 * eslDSQDATA_CHUNK_MAX + (eslDSQDATA_CHUNK_MAXSEQ + 1);
  ESL_ALLOC(chu->smem, sizeof(ESL_DSQ) * U);
  chu->psq = (uint32_t *) (chu->smem + U - 4*eslDSQDATA_CHUNK_MAX);
  

  /* We don't have any guarantees about the amount of metadata
   * associated with the N sequences, so <metadata> has to be a
   * reallocatable space. We make a lowball guess for the initial
   * alloc, on the off chance that the metadata size is small (names
   * only, no acc/desc): minimally, say 12 bytes of name, 3 \0's, and
   * 4 bytes for the taxid integer: call it 20.
   */
  chu->mdalloc = 20 * eslDSQDATA_CHUNK_MAXSEQ;
  ESL_ALLOC(chu->metadata, sizeof(char) * chu->mdalloc);

  return chu;
  
 ERROR:
  dsqdata_chunk_Destroy(chu);
  return NULL;
}


static void
dsqdata_chunk_Destroy(ESL_DSQDATA_CHUNK *chu)
{
  if (chu)
    {
      if (chu->metadata) free(chu->metadata);
      if (chu->smem)     free(chu->smem);
      if (chu->L)        free(chu->L);
      if (chu->taxid)    free(chu->taxid);
      if (chu->desc)     free(chu->desc);
      if (chu->acc)      free(chu->acc);
      if (chu->name)     free(chu->name);
      if (chu->dsq)      free(chu->dsq);
      free(chu);
    }
}


/*****************************************************************
 * 3. Loader and unpacker, the input threads
 *****************************************************************/

static void *
dsqdata_loader_thread(void *p)
{
  ESL_DSQDATA         *dd        = (ESL_DSQDATA *) p;
  ESL_DSQDATA_RECORD  *idx       = NULL;
  ESL_DSQDATA_CHUNK   *chu       = NULL;
  int                  nchunk    = 0;             // number of chunks we create, and need to destroy.
  int                  nidx      = 0;             // how many records in <idx>: usually MAXSEQ, until end
  int                  nload     = 0;             // how many sequences we load: >=1, <=nidx
  int                  ncarried  = 0;             // how many records carry over to next iteration: nidx-nload
  int                  nread     = 0;             // fread()'s return value
  int                  nmeta     = 0;             // how many bytes of metadata we want to read for this chunk
  int                  i0        = 0;             // absolute index of first record in <idx>, 0-offset
  int64_t              psq_last  = -1;            // psq_end for record i0-1
  int64_t              meta_last = -1;            // metadata_end for record i0-1
  int                  done      = FALSE;
  int                  status;
  
  ESL_ALLOC(idx, sizeof(ESL_DSQDATA_RECORD) * eslDSQDATA_CHUNK_MAXSEQ);

  while (! done)
    {

      /* Get a chunk - either by creating it, or recycling it.
       * We'll create up to <nconsumers>+2 of them.
       */
      if (nchunk < dd->nconsumers+2)
	{
	  chu = dsqdata_chunk_Create();  // TODO: check status
	  nchunk++;
	}
      else
	{
	  pthread_mutex_lock(&dd->recycling_mutex);
	  while (dd->recycling == NULL)
	    pthread_cond_wait(&dd->recycling_cv, &dd->recycling_mutex);
	  chu = dd->recycling;
	  dd->recycling = chu->nxt;                    // pop one off recycling stack
	  pthread_mutex_unlock(&dd->recycling_mutex);
	  pthread_cond_signal(&dd->recycling_cv);      // signal *after* unlocking mutex
	}
      
      /* Refill index. (The memmove is avoidable. Alt strategy: we could load in 2 frames)
       * The previous loop loaded packed sequence for <nload'> of the <nidx'> entries,
       * where the 's indicate the variable has carried over from prev iteration:
       *       |----- nload' ----||--- (ncarried) ---|
       *       |-------------- nidx' ----------------|
       * Now we're going to shift the remainder ncarried = nidx-nload to the left, then refill:
       *       |---- ncarried ----||--- (MAXSEQ-ncarried) ---|
       *       |-------------- MAXSEQ -----------------------|
       * while watching out for the terminal case where we run out of
       * data, loading less than (MAXSEQ-ncarried) records:
       *       |---- ncarried ----||--- nidx* ---|
       *       |------------- nidx --------------|
       * where the <nidx*> is what fread() returns to us.
       */
      i0      += nload;               // this chunk starts with seq #<i0>
      ncarried = (nidx - nload);
      memmove(idx, idx + nload, sizeof(ESL_DSQDATA_RECORD) * ncarried);
      nidx  = fread(idx + ncarried, sizeof(ESL_DSQDATA_RECORD), eslDSQDATA_CHUNK_MAXSEQ - ncarried, dd->ifp);
      nidx += ncarried;               // usually, this'll be MAXSEQ, unless we're near eof.
      
      if (nidx == 0) 
	{ // We're EOF. This chunk will be the empty EOF signal to unpacker, consumers.
	  chu->i0 = i0;
	  chu->N  = 0;
	  chu->pn = 0;
	  done    = TRUE;
	}
      else
	{
	  /* Figure out how many sequences we're going to load: <nload>
	   *  nload = max i : i <= MAXSEQ && idx[i].psq_end - psq_last <= CHUNK_MAX
	   */
	  ESL_DASSERT1(( idx[0].psq_end - psq_last <= eslDSQDATA_CHUNK_MAX ));
	  if (idx[nidx-1].psq_end - psq_last <= eslDSQDATA_CHUNK_MAX)
	    nload = nidx;
	  else
	    { // Binary search for nload = max_i idx[i-1].psq_end - lastend <= MAX
	      int righti = nidx;
	      int mid;
	      nload = 1;
	      while (righti - nload > 1)
		{
		  mid = nload + (righti - nload) / 2;
		  if (idx[mid-1].psq_end - psq_last <= eslDSQDATA_CHUNK_MAX) nload = mid;
		  else righti = mid;
		}                                                  
	    }
	  
	  /* Read packed sequence. */
	  chu->pn = idx[nload-1].psq_end - psq_last;
	  nread   = fread(chu->psq, sizeof(uint32_t), chu->pn, dd->sfp);
	  //printf("Read %d packed ints from seq file\n", nread);
	  ESL_DASSERT1(( nread == chu->pn )); // TODO: actually check this, it could fail
	      
	  /* Read metadata, reallocating if needed */
	  nmeta = idx[nload-1].metadata_end - meta_last;
	  if (nmeta > chu->mdalloc) {
	    ESL_REALLOC(chu->metadata, sizeof(char) * nmeta);   // should be realloc by doubling instead?
	    chu->mdalloc = nmeta;
	  }
	  nread  = fread(chu->metadata, sizeof(char), nmeta, dd->mfp);
	  ESL_DASSERT1(( nread == nmeta ));  // TODO: check this, fread() can fail

	  chu->i0   = i0;
	  chu->N    = nload;
	  psq_last  = idx[nload-1].psq_end;
	  meta_last = idx[nload-1].metadata_end;
	}

      /* Put the finished chunk into outbox;
       * unpacker will pick it up and unpack it.
       */
      pthread_mutex_lock(&dd->loader_outbox_mutex);  
      while (dd->loader_outbox != NULL) 
	pthread_cond_wait(&dd->loader_outbox_empty_cv, &dd->loader_outbox_mutex); 
      dd->loader_outbox = chu;   
      pthread_mutex_unlock(&dd->loader_outbox_mutex);
      pthread_cond_signal(&dd->loader_outbox_full_cv);
    }
  /* done == TRUE: we've sent the empty EOF chunk downstream, and now
   * we wait to get all our chunks back through the recycling, so we
   * can free them and exit cleanly. We counted them as they went out,
   * in <nchunk>, so we know how many need to come home.
   */

  while (nchunk)
    {
      pthread_mutex_lock(&dd->recycling_mutex);
      while (dd->recycling == NULL)                 // Readers may still be working, will Recycle() their chunks
	pthread_cond_wait(&dd->recycling_cv, &dd->recycling_mutex);
      while (dd->recycling != NULL) {               // Free entire stack, while we have the mutex locked.
	chu           = dd->recycling;   
	dd->recycling = chu->nxt;
	dsqdata_chunk_Destroy(chu);
	nchunk--;
      }
      pthread_mutex_unlock(&dd->recycling_mutex);
      /* Because the recycling is a stack, readers never have to wait
       * on a condition to Recycle(); the recycling, unlike the
       * outboxes, doesn't need to be empty.
       */
    }
  free(idx);
  pthread_exit(NULL);


 ERROR:   // TODO: collect the chunks in the error case, but what if threads itself are messed up?
  if (idx) free(idx);    
  pthread_exit(NULL);  // TODO: return an error status.
}



static void *
dsqdata_unpacker_thread(void *p)
{
  ESL_DSQDATA          *dd   = (ESL_DSQDATA *) p;
  ESL_DSQDATA_CHUNK    *chu  = NULL;
  int                   done = FALSE;
  ESL_DSQ              *dsq  = NULL;
  char                 *ptr;
  int                   i;                        // counter over seqs in chunk, 0..N-1
  int                   pos;			  // position in packed seq buffer, 0..pn-1
  int                   r;			  // counter over unpacked, concatenated digital seq residues
  uint32_t              v;			  // one packed integer
  int                   bitshift;

  while (! done)
    {
      /* Get a chunk from loader's outbox. Wait if necessary. */
      pthread_mutex_lock(&dd->loader_outbox_mutex);
      while (dd->loader_outbox == NULL) 
	pthread_cond_wait(&dd->loader_outbox_full_cv, &dd->loader_outbox_mutex); 
      chu = dd->loader_outbox;
      dd->loader_outbox  = NULL;
      pthread_mutex_unlock(&dd->loader_outbox_mutex);
      pthread_cond_signal(&dd->loader_outbox_empty_cv);

      /* Unpack the metadata.
       * If chunk is empty (N==0), it's the EOF signal - let it go straight out to a consumer.
       * (The first consumer that sees it will set the at_eof flag in <dd>, which all
       *  consumers check. So we only need the one empty EOF chunk to flow downstream.)
       */
      if (! chu->N)
	{
	  done = TRUE; // still need to pass the chunk along to a consumer.
	}
      else
	{
	  /* "Unpack" the metadata */
	  ptr = chu->metadata;
	  for (i = 0; i < chu->N; i++)
	    {
	      ESL_DASSERT1(( ptr < chu->metadata + chu->mdalloc ));
	      chu->name[i] = ptr;                           ptr = 1 + strchr(ptr, '\0');   ESL_DASSERT1(( ptr < chu->metadata + chu->mdalloc ));
	      chu->acc[i]  = ptr;                           ptr = 1 + strchr(ptr, '\0');   ESL_DASSERT1(( ptr < chu->metadata + chu->mdalloc ));
	      chu->desc[i] = ptr;                           ptr = 1 + strchr(ptr, '\0');   ESL_DASSERT1(( ptr < chu->metadata + chu->mdalloc ));
	      chu->taxid[i] = (int32_t) *((int32_t *) ptr); ptr += sizeof(int32_t);     
	    }
	  
	  /* Unpack sequence data */
	  i           = 0;
	  r           = 0;
	  dsq         = (ESL_DSQ *) chu->smem;
	  dsq[r]      = eslDSQ_SENTINEL;
	  chu->dsq[i] = &(dsq[r]);
	  r++;
	  for (pos = 0; pos < chu->pn; pos++)
	    {
	      v = chu->psq[pos];   // Must pick up the value, because of packed psq overlap w/ unpacked smem
	      if ( v & (1 << 31) ) // EOD bit set? Then this is last packed int for seq i
		{
		  for (bitshift = 25; bitshift >= 0 && ((v >> bitshift) & 31) != 31; bitshift -= 5)
		    dsq[r++] = (v >> bitshift) & 31;
		  chu->L[i] = &(dsq[r]) - chu->dsq[i] - 1;
		  i++;
		  if (i < chu->N) chu->dsq[i] = &(dsq[r]);
		  dsq[r++] = eslDSQ_SENTINEL;
		}
	      else
		{
		  dsq[r++] = (v >> 25) & 31;
		  dsq[r++] = (v >> 20) & 31;
		  dsq[r++] = (v >> 15) & 31;
		  dsq[r++] = (v >> 10) & 31;
		  dsq[r++] = (v >>  5) & 31;
		  dsq[r++] = (v >>  0) & 31;
		}
	    }
	  ESL_DASSERT1(( i == chu->N ));
	}

      /* Put unpacked chunk into the unpacker's outbox.
       * May need to wait for it to be empty/available.
       */
      pthread_mutex_lock(&dd->unpacker_outbox_mutex);
      while (dd->unpacker_outbox != NULL) 
	pthread_cond_wait(&dd->unpacker_outbox_empty_cv, &dd->unpacker_outbox_mutex); 
      dd->unpacker_outbox = chu;
      pthread_mutex_unlock(&dd->unpacker_outbox_mutex);
      pthread_cond_signal(&dd->unpacker_outbox_full_cv);      
    }

  pthread_exit(NULL);
}


/*****************************************************************
 * x. Writer
 *****************************************************************/
#ifdef eslDSQDATA_EXAMPLE2

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"


static int
pack5(ESL_DSQ *dsq, int n, uint32_t **ret_psq, int *ret_plen)
{
  uint32_t  *psq      = *ret_psq;   // might be NULL for initial alloc; might be reallocated for large plen
  int        plen     = 1 + n/6;    // +1 for sentinel. For DNA, with mixed pack size, will need to count more carefully
  int        i;                     // position in <dsq>
  int        pos      = 0;          // position in <psq>
  int        bitshift = 25;
  int        status;
  
  ESL_REALLOC(psq, sizeof(uint32_t) * ESL_MAX(1, plen));
  psq[pos] = (1 << 30);  // 5-bit pack bit is set on first packed int
  for (i = 1; i <= n; i++)
    {  
      psq[pos] |= (uint32_t) ( dsq[i] << bitshift );
      bitshift -= 5;
      if (bitshift < 0) 
	{
	  pos++;
	  psq[pos] = (1 << 30);  // 5-bit pack bit is set on next packed int
	  bitshift = 25;
	}
    }
  ESL_DASSERT1(( pos == plen-1 ));  // psq[pos] is the unfinished final packed int, which might become all sentinel

  psq[pos] |= (1 << 31);            // set the EOD bit on it
  while (bitshift >= 0)             // any remaining slots are set to sentinel 31, all 1's
    { 
      psq[pos] |= (uint32_t) (31 << bitshift);  
      bitshift -= 5;
    }

  *ret_psq  = psq;
  *ret_plen = plen;
  return eslOK;


 ERROR:
  *ret_psq  = NULL;
  *ret_plen = 0;
  return status;
}


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile_in> <binary seqfile_out>";
static char banner[] = "experimental: create binary database for esl_dsqdata";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *infile    = esl_opt_GetArg(go, 1);
  char           *basename  = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc       = esl_alphabet_Create(eslAMINO);
  int             format    = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp      = NULL;
  ESL_SQ         *sq        = esl_sq_CreateDigital(abc);
  char           *ifile     = NULL;
  char           *sfile     = NULL;
  char           *mfile     = NULL;
  ESL_DSQDATA_RECORD idx;
  uint32_t       *packsq    = NULL;
  FILE           *stubfp    = NULL;  // Human readable stub file, unparsed: <basename>
  FILE           *ifp       = NULL;  // Index file <basename>.dsqi
  FILE           *sfp       = NULL;  // Sequence file <basename>.dsqs
  FILE           *mfp       = NULL;  // Metadata file <basename>.dsqm
  int64_t         spos      = 0;     // current length of seq file (in uint32's)
  int64_t         mpos      = 0;     // current length of metadata file (in bytes)
  int             plen;
  int             n;
  int             status;

  esl_sprintf(&ifile, "%s.dsqi", basename);
  esl_sprintf(&sfile, "%s.dsqs", basename);
  esl_sprintf(&mfile, "%s.dsqm", basename);
  stubfp = fopen(basename, "w");
  ifp    = fopen(ifile, "wb");
  sfp    = fopen(sfile, "wb");
  mfp    = fopen(mfile, "wb");

  status = esl_sqfile_OpenDigital(abc, infile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      pack5(sq->dsq, sq->n, &packsq, &plen);
      fwrite(packsq, sizeof(uint32_t), plen, sfp);
      spos += plen;

      n = strlen(sq->name); fwrite(sq->name, sizeof(char), n+1, mfp);  mpos += n+1;
      n = strlen(sq->acc);  fwrite(sq->acc,  sizeof(char), n+1, mfp);  mpos += n+1;
      n = strlen(sq->desc); fwrite(sq->desc, sizeof(char), n+1, mfp);  mpos += n+1;
      fwrite( &(sq->tax_id), sizeof(int32_t), 1, mfp);                 mpos += sizeof(int32_t); 
      // TODO: byteswapping 
      
      idx.psq_end      = spos-1;  // could be -1, on 1st seq, if 1st seq L=0.
      idx.metadata_end = mpos-1; 
      fwrite(&idx, sizeof(ESL_DSQDATA_RECORD), 1, ifp);

      esl_sq_Reuse(sq);
    }

  fprintf(stubfp, "This is a test.\n");
  fprintf(stubfp, "If this were a real binary database...\n");

  fclose(stubfp);
  fclose(sfp);
  fclose(ifp);
  fclose(mfp);
  free(ifile);
  free(sfile);
  free(mfile);
  free(packsq);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*eslDSQDATA_EXAMPLE2*/


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef eslDSQDATA_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dsqdata.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name             type          default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <basename>";
static char banner[] = "example of using ESL_DSQDATA to read sequence db";

int
main(int argc, char **argv)
{
  ESL_GETOPTS       *go       = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET      *abc      = esl_alphabet_Create(eslAMINO);
  char              *basename = esl_opt_GetArg(go, 1);
  ESL_DSQDATA       *dd       = NULL;
  ESL_DSQDATA_CHUNK *chu      = NULL;
  int                i;
  int64_t            pos;
  int64_t            ct[128], total;
  int                x;
  int                ncpu     = 4;
  int                status;
  
  esl_dsqdata_Open(&abc, basename, ncpu, &dd);

  for (x = 0; x < 127; x++) ct[x] = 0;

  while ((status = esl_dsqdata_Read(dd, &chu)) != eslEOF)
    {
      for (i = 0; i < chu->N; i++)
      	for (pos = 1; pos <= chu->L[i]; pos++)
      	  ct[ chu->dsq[i][pos] ]++;

      esl_dsqdata_Recycle(dd, chu);
    }
  
  total = 0;
  for (x = 0; x < abc->Kp; x++) 
    {
      printf("%c  %" PRId64 "\n", abc->sym[x], ct[x]);
      total += ct[x];
    }
  printf("Total = %" PRId64 "\n", total);

  esl_alphabet_Destroy(abc);
  esl_dsqdata_Close(dd);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslDSQDATA_EXAMPLE*/
