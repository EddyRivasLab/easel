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
typedef struct {
  FILE *fp;		      /* Open file ptr                            */
  char *filename;	      /* Name of file (for diagnostics)           */
  int   do_gzip;	      /* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;	      /* TRUE if we're reading from stdin         */
  char  errbuf[eslERRBUFSIZE];/* parse error mesg. Size must match msa.h  */

  /* Our input buffer, whether using character-based parser [fread()]
   * or line-based parser (esl_fgets()).
   */
  char *buf;		      /* buffer for fread() or fgets() input      */
  off_t boff;		      /* disk offset to start of buffer           */
  int   balloc;		      /* allocated size of buf                    */
  int   nc;		      /* #chars in buf (usually full, less at EOF)*/ 
  int   pos;		      /* current parsing position in the buffer   */
  int   linenumber;	      /* What line of the file this is (1..N)     */

  /* Format-specific information
   */
  int   format;		      /* Format code of this file                    */
  int   is_linebased;	      /* TRUE for fgets() parsers; FALSE for fread() */
  int   addfirst;             /* TRUE to parse first line of seq record      */
  int   addend;	              /* TRUE to parse last line of seq record       */
  int   eof_is_ok;	      /* TRUE if record can end on EOF               */
  int  (*endTest)(char *);    /* ptr to function that tests if buffer is end */
  ESL_DSQ inmap[128];	      /* an input map, 0..127                        */

  /* If we have to GuessAlphabet(), we cache the first seq in the file */
  ESL_SQ *sq_cache;

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

  /* SSI augmentation confers random access of records in a seq file */
  char    *ssifile;	      /* path to expected SSI index file, even if nonexistent */
  int      rpl;		      /* residues per line; -1 if unset, 0 if inval */
  int      bpl;		      /* bytes per line; -1 if unset, 0 if inval    */
  int      lastrpl;	      /* tmp var used only when indexing            */
  int      lastbpl;	      /* ditto                                      */
#if defined(eslAUGMENT_SSI)
  ESL_SSI *ssi;		/* open ESL_SSI index, or NULL if none     */
#else
  void    *ssi;		/* NULL */
#endif /*eslAUGMENT_SSI*/
} ESL_SQFILE;

/* fread() is apparently the fastest portable way to input from disk;
 * the READBUFSIZE is the fixed size of a block to bring in at one time,
 * for character-based parsers (like the FASTA parser).
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


extern int  esl_sqfile_Open(char *seqfile, int fmt, char *env, ESL_SQFILE **ret_sqfp);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);
#ifdef eslAUGMENT_ALPHABET
extern int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);
#endif

extern int   esl_sqio_Read(ESL_SQFILE *sqfp, ESL_SQ *s);
extern int   esl_sqio_Write(FILE *fp, ESL_SQ *s, int format);
extern int   esl_sqio_Echo (FILE *ofp, ESL_SQFILE *sqfp);
extern int   esl_sqio_WhatFormat(FILE *fp);
extern int   esl_sqio_FormatCode(char *fmtstring);
extern char *esl_sqio_DescribeFormat(int fmt);
extern int   esl_sqio_IsAlignment(int fmt);
extern int   esl_sqio_Position(ESL_SQFILE *sqfp, off_t r_off);
extern int   esl_sqio_Rewind  (ESL_SQFILE *sqfp);

#ifdef eslAUGMENT_SSI
extern int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
extern int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
extern int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);
#endif /*eslAUGMENT_SSI*/
#endif /*!ESL_SQIO_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
