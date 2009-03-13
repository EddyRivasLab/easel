/* Unaligned sequence file i/o.
 * 
 * SVN $Id$
 */
#ifndef ESL_SQIO_INCLUDED
#define ESL_SQIO_INCLUDED

#include <stdio.h>
#include "esl_sq.h"
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_MSA
#include "esl_msa.h"
#endif


/* ESL_SQFILE:
 * An open sequence file for reading.
 */
typedef struct esl_sqfile_s {
  FILE *fp;		      /* Open file ptr                            */
  char *filename;	      /* Name of file (for diagnostics)           */
  int   do_gzip;	      /* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;	      /* TRUE if we're reading from stdin         */
  char  errbuf[eslERRBUFSIZE];/* parse error mesg. Size must match msa.h  */

  /* all input first gets buffered in memory; this gives us enough
   * recall to use Guess*() functions even in nonrewindable streams
   */
  char    *mem;		      /* buffered input                           */
  int      allocm;	      /* <mem> size, multiples of eslREADBUFSIZE  */
  int      mn;		      /* number of chars in <mem> (up to allocm)  */
  int      mpos;	      /* pos of next <buf> to load from <mem>     */
  off_t    moff;	      /* disk offset to start of <mem>            */
  int      is_recording;      /* TRUE if we need to keep buffering more   */

  /* input is either character-based [fread()] or line-based (esl_fgets())*/
  char    *buf;		      /* buffer for fread() or fgets() input      */
  off_t    boff;	      /* disk offset to start of buffer           */
  int      balloc;	      /* allocated size of buf                    */
  int      nc;		      /* #chars in buf (usually full, less at EOF)*/ 
  int      bpos;	      /* current position in the buffer (0..nc-1) */
  int64_t  L;		      /* #residues seen so far in current seq     */
  int64_t  linenumber;	      /* What line of the file  (1..N; -1=unknown)*/
  off_t    bookmark_offset;   /* bookmark fwd position before reversing...*/
  int64_t  bookmark_linenum;  /* in both linenumber and disk offset       */

  /* In digital mode, we have an alphabet ptr                             */
  int   do_digital;	      /* TRUE if we're reading in digital mode    */
#if defined(eslAUGMENT_ALPHABET)  
  const ESL_ALPHABET *abc;
#else
  void               *abc;
#endif

  /* Format-specific configuration                                           */
  int   format;		      /* Format code of this file                    */
  int   is_linebased;	      /* TRUE for fgets() parsers; FALSE for fread() */
  int   eof_is_ok;	      /* TRUE if record can end on EOF               */
  int  (*parse_header)(struct esl_sqfile_s *, ESL_SQ *sq);
  int  (*parse_end)   (struct esl_sqfile_s *, ESL_SQ *sq); 
  ESL_DSQ inmap[128];	      /* an input map, 0..127                        */

  /* MSA augmentation confers reading MSA files as sequential seq files. */
#if defined(eslAUGMENT_MSA)
  ESL_MSAFILE *afp;	      /* open ESL_MSAFILE for reading           */
  ESL_MSA     *msa;	      /* preloaded alignment to draw seqs from  */
  int          idx;	      /* index of next seq to return, 0..nseq-1 */
#else
  void        *afp; 	      /* NULL */
  void        *msa;           /* NULL */
  int          idx;           /* 0    */
#endif /*eslAUGMENT_MSA*/

  /* SSI augmentation confers random access of records in a seq file        */
  char    *ssifile;	      /* path to expected SSI index file            */
  int      rpl;		      /* residues per line in file; -1=unset 0=inval*/
  int      bpl;		      /* bytes per line in file; -1=unset, 0=inval  */
  int      currpl;	      /* residues on current line (-1=unknown)      */
  int      curbpl;	      /* bytes on current line    (-1=unknown)      */
  int      prvrpl;	      /* residues on previous line                  */
  int      prvbpl;	      /* bytes on previous line                     */
#if defined(eslAUGMENT_SSI)
  ESL_SSI *ssi;		/* open ESL_SSI index, or NULL if none     */
#else
  void    *ssi;		/* NULL */
#endif /*eslAUGMENT_SSI*/
} ESL_SQFILE;


/* eslREADBUFSIZE is the fixed size of a block to bring in at one time,
 * in character-based (fread()) parsers (like the FASTA parser).
 */
#define eslREADBUFSIZE  4096


/* Unaligned file format codes
 * These codes are coordinated with the msa module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqio, for unaligned formats
 *   - >100  is reserved for msa, for aligned formats
 */
#define eslSQFILE_UNKNOWN 0
#define eslSQFILE_FASTA   1
#define eslSQFILE_EMBL    2	/* EMBL/Swissprot/TrEMBL */
#define eslSQFILE_GENBANK 3	/* Genbank */
#define eslSQFILE_DDBJ    4	/* DDBJ (currently passed to Genbank parser */
#define eslSQFILE_UNIPROT 5     /* Uniprot (passed to EMBL parser) */


extern int  esl_sqfile_Open(const char *seqfile, int fmt, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_GuessFileFormat(ESL_SQFILE *sqfp, int *ret_format);
extern int  esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);

#ifdef eslAUGMENT_ALPHABET
extern int  esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
extern int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);
#endif

extern int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp);
extern int   esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignoredchars);
extern int   esl_sqio_AcceptAs(ESL_SQFILE *sqfp, char *xchars, char readas);
extern int   esl_sqio_EncodeFormat(char *fmtstring);
extern char *esl_sqio_DecodeFormat(int fmt);
extern int   esl_sqio_IsAlignment(int fmt);

extern int   esl_sqio_Read      (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadInfo  (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
extern int   esl_sqio_Echo      (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);

#ifdef eslAUGMENT_SSI
extern int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
extern int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
extern int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);

extern int   esl_sqio_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

extern int   esl_sqio_Write(FILE *fp, ESL_SQ *s, int format);

#endif /*!ESL_SQIO_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
