#ifndef eslDSQDATA_INCLUDED
#define eslDSQDATA_INCLUDED

#include "easel.h"
#include "esl_alphabet.h"

#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

#define eslDSQDATA_CHUNK_MAXSEQ     4096
#define eslDSQDATA_CHUNK_MAX      262144      // maximum number of packed uint32't to load in one fread. 

typedef struct esl_dsqdata_record_s {
  int64_t  metadata_end;
  int64_t  psq_end;
} ESL_DSQDATA_RECORD;


typedef struct esl_dsqdata_chunk_s {
  int64_t   i0;
  int       N;
  ESL_DSQ **dsq;

  char    **name;         // Names, \0 terminated. <name>, <acc>, <desc> are ptrs into <metadata> buffer.
  char    **acc;          // Optional accessions, \0 terminated; "\0" if none. 
  char    **desc;         // Optional descriptions, \0 terminated; "\0" if none
  int32_t  *taxid;        // NCBI taxonomy identifiers. (-1 if none)

  int64_t  *L;            // Sequence lengths, in residues. The unpacker figures these out.

  char     *smem;         // Unpacked (dsq[]) and packed (psq) data ptrs share this allocation.
  uint32_t *psq;          // Pointer into smem; packed data fread()'s go here.
  int       pn;           // how many uint32's are loaded in <psq>

  char     *metadata;     // Raw fread() buffer of all name/acc/desc/taxid data.
  int       mdalloc;      // Current allocation size for <metadata> in bytes

  struct esl_dsqdata_chunk_s *nxt;
} ESL_DSQDATA_CHUNK;


/* ESL_DSQDATA:  high performance sequence database input
 *               This is the object that we open and read from.
 */              
typedef struct esl_dsqdata_s {
  char *basename;
  FILE *ifp;
  FILE *sfp;
  FILE *mfp;

  int   nconsumers;

  pthread_t          loader_t;
  ESL_DSQDATA_CHUNK *loader_outbox;
  pthread_mutex_t    loader_outbox_mutex;
  pthread_cond_t     loader_outbox_full_cv;  
  pthread_cond_t     loader_outbox_empty_cv; 

  pthread_t          unpacker_t;
  ESL_DSQDATA_CHUNK *unpacker_outbox;
  pthread_mutex_t    unpacker_outbox_mutex;
  pthread_cond_t     unpacker_outbox_full_cv;
  pthread_cond_t     unpacker_outbox_empty_cv;
  int                at_eof;                    // flag that goes up when we reach end of the input file; 
                                                // raising <at_eof> is inside the unpacker_outbox mutex

  ESL_DSQDATA_CHUNK *recycling;
  pthread_mutex_t    recycling_mutex;
  pthread_cond_t     recycling_cv;

  /* Unfortunately, pristine cleanup after errors requires that we use
   * separate booleans to keep track of which thread resources have
   * been created.
   */
  int lt_c;  int lom_c;  int lof_c;  int loe_c;
  int ut_c;  int uom_c;  int uof_c;  int uoe_c;
  int rm_c;  int r_c;

  ESL_ALPHABET *abc_r;                   // Copy of ptr to the alphabet the caller told us to read in.
  char          errbuf[eslERRBUFSIZE];   // User-directed error message in case of a failed open or read.
} ESL_DSQDATA;  
  
  
extern int  esl_dsqdata_Open   (ESL_ALPHABET **byp_abc, char *basename, int nconsumers, ESL_DSQDATA **ret_dd);
extern int  esl_dsqdata_Read   (ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK **ret_chu);
extern int  esl_dsqdata_Recycle(ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK *chu);
extern void esl_dsqdata_Close  (ESL_DSQDATA *dd);

#endif /*eslDSQDATA_INCLUDED*/
