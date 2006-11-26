/* esl_ssi.h [created from esl_ssi.h.in by ./configure]
 * 
 * "simple sequence indices": 
 * fast sequence record lookup in large files by keywords, such
 * as names or accessions.
 * 
 * SVN $Id$
 * SRE, Thu Mar  2 15:54:51 2006 [St. Louis]
 */
#ifndef ESL_SSI_INCLUDED
#define ESL_SSI_INCLUDED

/* Conditional inclusion of non-ANSI C headers that make life easier.
 */
#ifdef HAVE_STDINT_H	
#include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

/* Limits.
 */
#define eslSSI_MAXFILES 32767	     /* 2^15-1 */
#define eslSSI_MAXKEYS  2147483647L  /* 2^31-1 */
#define eslSSI_MAXRAM   256	     /* >256MB indices trigger external sort */

#ifndef HAVE_FSEEKO
#define fseeko fseek
#define ftello ftell
#endif 


/* ESL_SSI
 * Using an existing SSI index file.
 */ 
typedef struct {
  FILE      *fp;              /* open SSI index file                 */
  uint32_t   flags;	      /* optional behavior flags             */
  uint16_t   nfiles;          /* number of files = 16 bit int        */
  uint32_t   nprimary;        /* number of primary keys              */
  uint32_t   nsecondary;      /* number of secondary keys            */
  uint32_t   flen;            /* length of filenames (inc '\0')      */
  uint32_t   plen;            /* length of primary keys (inc '\0')   */
  uint32_t   slen;            /* length of secondary keys (inc '\0') */
  uint32_t   frecsize;        /* # bytes in a file record            */
  uint32_t   precsize;        /* # bytes in a primary key record     */
  uint32_t   srecsize;        /* # bytes in a secondary key record   */
  off_t      foffset;         /* disk offset, start of file records  */
  off_t      poffset;         /* disk offset, start of pri key recs  */
  off_t      soffset;         /* disk offset, start of sec key recs  */
  
  char imode;                 /* mode for index file offsets, 32 v. 64 bit */
  char smode;                 /* mode for seq file offsets, 32 v. 64 bit */

  /* File information:
   */
  char     **filename;        /* list of file names [0..nfiles-1]    */
  uint32_t  *fileformat;      /* file formats                        */
  uint32_t  *fileflags;	      /* optional per-file behavior flags    */
  uint32_t  *bpl;             /* bytes per line in file              */
  uint32_t  *rpl;             /* residues per line in file           */
} ESL_SSI;

/* Flags for the <ssi->flags> bit vector. 
 * Can add flags, but don't change the existing ones; they're locked
 * into reverse compatibility with older SSI files.
 */
#define eslSSI_USE64        (1 << 0)  /* key offsets (in the indexed files) are 64-bit */
#define eslSSI_USE64_INDEX  (1 << 1)  /* index file itself is so large that its offsets are 64-bit */

/* Flags for the <ssi->fileflags> bit vectors.
 */
#define eslSSI_FASTSUBSEQ   (1<<0)    /* we can do fast subseq lookup calculations on this file */


/* ESL_NEWSSI
 * Used to create a new SSI index.
 */
typedef struct {		/* Primary key data: */
  char      *key;               /* key name          */
  uint16_t   fnum;		/* file number       */
  off_t      r_off;		/* record offset     */
  off_t      d_off;		/* data offset       */
  uint32_t   len;		/* sequence length   */
} ESL_PKEY;

typedef struct {		/* Secondary key data: */
  char        *key;             /* secondary key name  */
  char        *pkey;            /* primary key name    */ 
} ESL_SKEY;

typedef struct {
  int         external;	        /* TRUE if pkeys and skeys are on disk    */
  int         max_ram;	        /* threshold in MB to trigger extern sort */

  char      **filenames;
  uint32_t   *fileformat;
  uint32_t   *bpl;
  uint32_t   *rpl;
  uint32_t    flen;		/* length of longest filename, inc '\0' */
  uint16_t    nfiles;
  
  ESL_PKEY   *pkeys;
  uint32_t    plen;	        /* length of longest pkey, including '\0' */
  uint32_t    nprimary;
  char       *ptmpfile;		/* primary key tmpfile name, for extern sort */
  FILE       *ptmp;	        /* handle on open ptmpfile */

  ESL_SKEY   *skeys;
  uint32_t    slen;        	/* length of longest skey, including '\0' */
  uint32_t    nsecondary;
  char       *stmpfile;		/* secondary key tmpfile name, for extern sort */
  FILE       *stmp;	        /* handle on open ptmpfile */
} ESL_NEWSSI;


#define eslSSI_FCHUNK  16	/* chunk size for file name reallocation */
#define eslSSI_KCHUNK  128	/* and for key reallocation              */



/* 1. Using SSI indices
 */
extern int  esl_ssi_Open(char *filename, ESL_SSI **ret_ssi);
extern void esl_ssi_Close(ESL_SSI *ssi);
extern int  esl_ssi_FindName(ESL_SSI *ssi, char *key,
			     uint16_t *ret_fh, off_t *ret_offset);
extern int  esl_ssi_FindNumber(ESL_SSI *ssi, int nkey,
			       uint16_t *ret_fh, off_t *ret_offset);
extern int  esl_ssi_FindSubseq(ESL_SSI *ssi, char *key, long requested_start,
			       uint16_t *ret_fh, off_t *record_offset, off_t *data_offset, 
			       long *ret_actual_start);
extern int  esl_ssi_FileInfo(ESL_SSI *ssi, uint16_t fh,
			     char **ret_filename, int *ret_format);



/* 2. Creating SSI indices
 */
extern ESL_NEWSSI *esl_newssi_Create(void);
extern int  esl_newssi_AddFile(ESL_NEWSSI *ns, char *filename,
			       int fmt, uint16_t *ret_fh);
extern int  esl_newssi_SetSubseq(ESL_NEWSSI *ns, uint16_t fh,
				 int bpl, int rpl);
extern int  esl_newssi_AddKey(ESL_NEWSSI *ns, char *key, uint16_t fh, 
			      off_t r_off, off_t d_off, uint32_t L);
extern int  esl_newssi_AddAlias(ESL_NEWSSI *ns, char *alias, char *key);
extern int  esl_newssi_Write(FILE *fp, ESL_NEWSSI *ns);
extern void esl_newssi_Destroy(ESL_NEWSSI *ns);


/* 3. Binary file portability
 */
extern void     esl_byteswap(char *swap, int nbytes);
extern uint16_t esl_ntoh16(uint16_t netshort);
extern uint32_t esl_ntoh32(uint32_t netlong);
extern uint64_t esl_ntoh64(uint64_t net_int64);
extern uint16_t esl_hton16(uint16_t hostshort);
extern uint32_t esl_hton32(uint32_t hostlong);
extern uint64_t esl_hton64(uint64_t host_int64);
extern int      esl_fread_i16(FILE *fp, uint16_t *ret_result);
extern int      esl_fread_i32(FILE *fp, uint32_t *ret_result);
extern int      esl_fread_i64(FILE *fp, uint64_t *ret_result);
extern int      esl_fwrite_i16(FILE *fp, uint16_t n);
extern int      esl_fwrite_i32(FILE *fp, uint32_t n);
extern int      esl_fwrite_i64(FILE *fp, uint64_t n);
extern int	esl_fread_offset(FILE *fp, int mode, off_t *ret_offset);
extern int      esl_fwrite_offset(FILE *fp, off_t offset);


#endif /* ESL_SSI_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
