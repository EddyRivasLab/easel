/* esl_ssi.h
 * 
 * "simple sequence indices": fast record lookups in large files by keyword.
 * 
 * SVN $Id$
 * SRE, Thu Mar  2 15:54:51 2006 [St. Louis]
 */
#ifndef ESL_SSI_INCLUDED
#define ESL_SSI_INCLUDED


/* Limits.
 */
#define eslSSI_MAXFILES 32767	     /* 2^15-1 */
#define eslSSI_MAXKEYS  2147483647L  /* 2^31-1 */
#define eslSSI_MAXRAM   200	     /* >200MB indices trigger external sort */




/* ESL_OFFSET
 * Portable file position offsets compatible with arithmetic operations.
 * (The portable fpos_t type is opaque.)
 *
 * Use the union to save space, since the two offset types are
 * mutually exclusive, controlled by "mode"
 */
typedef struct {
  char mode;			/* eslOFFSET32, for example               */
  union {
    esl_uint32   i32;           /* an offset that fseek() can use         */
    esl_uint64   i64;           /* an offset that e.g. fseeko64() can use */
  } off;
} ESL_OFFSET;

#define eslOFFSET32    0
#define eslOFFSET64    1


/* ESL_SSI
 * Using an existing SSI index.
 */ 
typedef struct {
  FILE        *fp;              /* open SSI index file                 */
  esl_uint32   flags;           /* optional behavior flags             */
  esl_uint16   nfiles;          /* number of files = 16 bit int        */
  esl_uint32   nprimary;        /* number of primary keys              */
  esl_uint32   nsecondary;      /* number of secondary keys            */
  esl_uint32   flen;            /* length of filenames (inc '\0')      */
  esl_uint32   plen;            /* length of primary keys (inc '\0')   */
  esl_uint32   slen;            /* length of secondary keys (inc '\0') */
  esl_uint32   frecsize;        /* # bytes in a file record            */
  esl_uint32   precsize;        /* # bytes in a primary key record     */
  esl_uint32   srecsize;        /* # bytes in a secondary key record   */
  ESL_OFFSET   foffset;         /* disk offset, start of file records  */
  ESL_OFFSET   poffset;         /* disk offset, start of pri key recs  */
  ESL_OFFSET   soffset;         /* disk offset, start of sec key recs  */
  
  char imode;                   /* mode for index file offsets, 32 v. 64 bit    */
  char smode;                   /* mode for sequence file offsets, 32 v. 64 bit */

  /* File information:
   */
  char       **filename;        /* list of file names [0..nfiles-1]    */
  esl_uint32  *fileformat;      /* file formats                        */
  esl_uint32  *fileflags;       /* optional per-file behavior flags    */
  esl_uint32  *bpl;             /* bytes per line in file              */
  esl_uint32  *rpl;             /* residues per line in file           */
} ESL_SSI;

/* optional per-index behavior flags in ESL_SSIFILE's flags:
 */
#define eslSSI_USE64        1<<0	/* seq offsets are 64-bit        */
#define eslSSI_USE64_INDEX  1<<1	/* index file offsets are 64-bit */

/* optional per-file behavior flags in ESL_SSIFILE's fileflags:
 */
#define eslSSI_FAST_SUBSEQ  1<<0	/* can do subseq lookup in this file */









#endif /* ESL_SSI_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
