/* sqio.h
 * Sequence file i/o.
 * 
 * SVN $Id$
 */

#include <stdio.h>

#define SQ_MAXNAME     64	/* maximum length of sequence name/accession */
#define SQ_MAXDESC    128	/* maximum length of desc line               */
#define SQ_CHUNKSIZE 2048	/* defines allocation chunk size for seqs    */


#define ESL_SQFORMAT_FASTA  1


struct esl_sqfile_s {
  FILE *fp;			/* Open file ptr                            */
  char *filename;		/* Name of file (for diagnostics)           */
  int   format;			/* Format of this file                      */
  int   do_gzip;		/* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;		/* TRUE if we're reading from stdin         */

  int   linenumber;		/* What line of the file we're on (0..N-1) */

  char *buf;			/* dynamic esl_parse_fgets() buffer */
  int   buflen;			/* current alloc length of buf      */

  /* SSI indexing eventually goes here, including rpl,bpl counting;
   * xref squid.
   */
};
typedef struct esl_sqfile_s ESL_SQFILE;

struct esl_sequence_s {
  char  name[SQ_MAXNAME];	/* name        */
  char  acc[SQ_MAXNAME];	/* accession   */
  char  desc[SQ_MAXDESC];	/* description */

  char *seq;			/* growing sequence during a parse   */
  char *ss;			/* secondary structure annotation    */
  int   len;			/* current length of seq             */
  int   alloclen;		/* current allocation length for seq */
  int   allocchunk;		/* allocation chunksize (settable)   */
};
typedef struct esl_sequence_s ESL_SQ;



  
