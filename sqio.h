/* sqio.h
 * Sequence file i/o.
 * 
 * SVN $Id$
 */
#ifndef ESL_SQIO_INCLUDED
#define ESL_SQIO_INCLUDED

#include <stdio.h>

/* name, accession, description, and sequence itself are of unlimited
 * length, but they are initially allocated to something sensible, as
 * set below, in the hope that any given ESL_SQ object only has to
 * make one malloc() call for them. These lengths are inclusive of the
 * \0 NUL character (so ESL_SQ_NAMELEN of 32 means we expect names <=
 * 31 chars).
 * 
 * The reallocation rule is to double the allocation every time 
 * we need to reallocate. So, with the SEQCHUNK set to 256, for example,
 * sequences may be allocated with length 256, 512, 1024, 2048... etc.
 * This results in >50% utilization of our memory, which means
 * some wastage, but we have to compromise with efficiency in
 * sequence reading speed. The esl_sq_Squeeze() function provides
 * for optimizing memory after seqs have been read, at a cost of
 * a ~20% speed hit (for memcpy()'ing within the seq to perfected space.
 */
#define ESL_SQ_NAMECHUNK   32	/* allocation unit for name         */
#define ESL_SQ_ACCCHUNK    32	/* allocation unit for accession    */
#define ESL_SQ_DESCCHUNK  128	/* allocation unit for description  */
#define ESL_SQ_SEQCHUNK   256	/* allocation unit for seqs         */

/* fread() is apparently the fastest portable way to input from disk;
 * the READBUFSIZE is the fixed size of a block to bring in at one time,
 * for character-based parsers (like the FASTA parser).
 */
#define ESL_READBUFSIZE  4096

#define ESL_SQFORMAT_FASTA  1


typedef struct {
  FILE *fp;			/* Open file ptr                            */
  char *filename;		/* Name of file (for diagnostics)           */
  int   format;			/* Format of this file                      */
  int   do_gzip;		/* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;		/* TRUE if we're reading from stdin         */

  int  *inmap;			/* pointer to an input map, 0..127  */

  char  buf[ESL_READBUFSIZE];	/* buffer for fread() block input           */
  int   nc;			/* # of valid chars in buf (usually full, but less at EOF) */  
  int   pos;			/* current parsing position in the buffer   */
  int   linenumber;		/* What line of the file this is (1..N)     */

  /* SSI indexing eventually goes here, including rpl,bpl counting;
   * xref squid.
   */
} ESL_SQFILE;

/* Allocation of name, acc, desc is all in one malloc, of length
 * nalloc+aalloc+dalloc.
 */
typedef struct {
  char *name;	        /* name                     */
  int   nalloc;		/* allocated length of name */
  
  char *acc;	        /* accession                     */
  int   aalloc;	       	/* allocated length of accession */

  char *desc;		/* description */
  int   dalloc;		/* allocated length of description */

  char *seq;		/* growing sequence during a normal parse, or NULL */
  char *dsq;		/* growing digital seq during digital parse, or NULL */
  char *ss;		/* secondary structure annotation    */
  int   n;		/* current length of seq             */
  int   salloc;		/* current allocation length for seq */

  char *optmem;		/* optimized mem storage area; see esl_sq_Squeeze() */
} ESL_SQ;


extern ESL_SQ *esl_sq_Create(void);
extern int     esl_sq_Inflate(ESL_SQ *sq);
extern int     esl_sq_Reuse(ESL_SQ *sq);
extern int     esl_sq_Squeeze(ESL_SQ *sq);
extern void    esl_sq_Deflate(ESL_SQ *sq);
extern void    esl_sq_Destroy(ESL_SQ *sq);

extern int  esl_sqfile_OpenFASTA(char *seqfile, ESL_SQFILE **ret_sqfp);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);

extern int  esl_sio_ReadFASTA(ESL_SQFILE *sqfp, ESL_SQ *s);

#endif /*!ESL_SQIO_INCLUDED*/
