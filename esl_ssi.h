/* esl_ssi.h [created from esl_ssi.h.in by ./configure]
 * 
 * "simple sequence indices": fast record lookups in large files by keyword.
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
#define eslSSI_MAXRAM   200	     /* >200MB indices trigger external sort */

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
 * Creating a new ssi index.
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
  char       *ptmpfile;	        /* name of tmp file, for external sort mode */
  FILE       *ptmp;	        /* handle on open ptmpfile */

  ESL_SKEY   *skeys;
  uint32_t    slen;        	/* length of longest skey, including '\0' */
  uint32_t    nsecondary;
  char       *stmpfile;	        /* name of tmp file, for external sort mode */
  FILE       *stmp;	        /* handle on open ptmpfile */
} ESL_NEWSSI;



/* Binary file portability
 */
extern void     esl_byteswap(char *swap, int nbytes);
extern uint16_t esl_ntoh16(uint16_t netshort);
extern uint32_t esl_ntoh32(uint32_t netlong);
extern uint64_t esl_ntoh64(uint64_t net_int64);
extern uint16_t esl_hton16(uint16_t hostshort);
extern uint32_t esl_hton32(uint32_t hostlong);
extern uint64_t esl_hton64(uint64_t host_int64);
extern int      esl_fread_i16(FILE *fp, esl_uint16 *ret_result);
extern int      esl_fread_i32(FILE *fp, esl_uint32 *ret_result);
extern int      esl_fread_i64(FILE *fp, esl_uint64 *ret_result);
extern int      esl_fwrite_i16(FILE *fp, esl_uint16 n);
extern int      esl_fwrite_i32(FILE *fp, esl_uint32 n);
extern int      esl_fwrite_i64(FILE *fp, esl_uint64 n);
extern int	esl_fread_offset(FILE *fp, int mode, off_t *ret_offset);
extern int      esl_fwrite_offset(FILE *fp, off_t *offset);


#endif /* ESL_SSI_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
