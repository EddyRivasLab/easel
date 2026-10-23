/* sqfile: unaligned sequence file i/o.
 * 
 * To do: 
 *   :: Write tests for bad format detection, making sure that
 *      linenumber report is correct.
 */
#ifndef eslSQFILE_INCLUDED
#define eslSQFILE_INCLUDED
#include <esl_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"

/* Unaligned file format codes
 * These codes are coordinated with the msafile module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqfile,  for unaligned formats
 *   - >100  is reserved for msafile, for aligned formats
 */
#define eslSQFILE_UNKNOWN      0
#define eslSQFILE_FASTA        1   // FASTA format
#define eslSQFILE_EMBL         2   // EMBL DNA sequence
#define eslSQFILE_GENBANK      3   // Genbank DNA sequence
#define eslSQFILE_DDBJ         4   // DDBJ (currently identical to GenBank parser)
#define eslSQFILE_UNIPROT      5   // UniProt (currently identical to EMBL parser) 
                            // 6 was eslSQFILE_NCBI, NCBI blast db, v4, single file; now removed
#define eslSQFILE_DAEMON       7   // Farrar format, hmmpgmd queries: fasta + // terminator 
#define eslSQFILE_HMMPGMD      8   // Farrar hmmpgmd database format: fasta + # header 
#define eslSQFILE_FMINDEX      9   // Pressed FM-index format used in HMMER

/* eslSQFILE_READBUFSIZE is the fixed size of a block to bring in at one time,
 * in character-based (fread()) parsers (like the FASTA parser).
 */
#define eslSQFILE_READBUFSIZE  4096

/* eslSQFILE_BLOCK_MAXRES is the total # of residues in a block of protein seqs
 * or a window of DNA (with long_target on), specifically in _ReadBlock().
 */
#define eslSQFILE_BLOCK_MAXRES (1024 * 1024)


/* ESL_SQFILE:
 * An open sequence file for reading.
 */
typedef struct esl_sqfile_s {
  char *filename;	       // name of file that was opened, possibly including path (for diagnostics)
  FILE *fp;           	       // open file ptr
  char  errbuf[eslERRBUFSIZE]; // parse error mesg.  size must match msa.h

  int   do_digital;	       // TRUE if we're reading in digital mode
  int   do_gzip;	       // TRUE if we're reading from gzip -dc pipe
  int   do_stdin;	       // TRUE if we're reading from stdin
  int   do_buffer;             // TRUE if we're reading from a buffer
  const ESL_ALPHABET *abc;     // valid in digital mode; NULL in text mode 

  /* Format-specific configuration
   */
  int     format;	       // Format code of this file
  ESL_DSQ inmap[128];	       // an input map, 0..127
  int     is_linebased;	       // TRUE for fgets() parsers; FALSE for fread()
  int     eof_is_ok;	       // TRUE if record can end on EOF
  int    (*parse_header)(struct esl_sqfile_s *, ESL_SQ *sq);
  int    (*parse_end)   (struct esl_sqfile_s *, ESL_SQ *sq); 

  /* All input first gets buffered in memory; this gives us enough
   * recall to use Guess*() functions even in nonrewindable streams
   */
  char    *mem;		       // buffered input
  int      allocm;	       // <mem> size, multiples of eslSQFILE_READBUFSIZE
  int      mn;		       // number of chars in <mem> (up to allocm)
  int      mpos;	       // pos of next <buf> to load from <mem>
  off_t    moff;	       // disk offset to start of <mem>
  int      is_recording;       // TRUE if we need to keep buffering more

  /* Input is either character-based [fread()] or line-based [esl_fgets()]
   */
  char    *buf;		      // buffer for fread() or fgets() input
  off_t    boff;	      // disk offset to start of buffer
  int      balloc;	      // allocated size of buf
  int      nc;		      // #chars in buf (usually full, less at EOF)
  int      bpos;	      // current position in the buffer (0..nc-1)
  int64_t  L;		      // #residues seen so far in current seq
  int64_t  linenumber;	      // What line of the file  (1..N; -1=unknown)
  off_t    bookmark_offset;   // bookmark fwd position before reversing...
  int64_t  bookmark_linenum;  // in both linenumber and disk offset

  /* MSA files can be read as sequential seq files.
   */
  ESL_MSAFILE  *afp;	      // open ESL_MSAFILE for reading
  ESL_MSA      *msa;	      // preloaded alignment to draw seqs from
  int           idx;	      // index of next seq to return, 0..nseq-1

  /* SSI indexes allow fast random access of records in a seq file
   */
  char    *ssifile;          // path to expected SSI index file
  int      rpl;              // residues per line in file; -1=unset 0=inval
  int      bpl;              // bytes per line in file; -1=unset, 0=inval
  int      currpl;           // residues on current line (-1=unknown)
  int      curbpl;           // bytes on current line    (-1=unknown)
  int      prvrpl;           // residues on previous line
  int      prvbpl;           // bytes on previous line
  int      maxrpl;           // max residues on any seq line so far
  int      maxbpl;           // max bytes on any seq line so far
  int      maxgap;           // max nonresidues (e.g. spaces) on any seq line so far
  ESL_SSI *ssi;              // open ESL_SSI index, or NULL if none
} ESL_SQFILE;


/* ESL_SQCACHE:
 * A entire database cached into memory.
 */
typedef struct esl_sqcache_s {
  char               *filename;    /* Name of file (for diagnostics)              */
  int                 format;      /* Format code of this file                    */

  const ESL_ALPHABET *abc;         /* alphabet for database                       */

  uint32_t            seq_count;   /* number of sequences                         */
  uint64_t            res_count;   /* number of residues                          */
  uint32_t            max_seq;     /* longest sequence                            */

  ESL_SQ             *sq_list;     /* list of cached sequences [0 .. seq_count-1] */

  void               *residue_mem; /* memory holding the residues                 */
  void               *header_mem;  /* memory holding the header strings           */

  uint64_t            res_size;    /* size of residue memory allocation           */
  uint64_t            hdr_size;    /* size of header memory allocation            */
} ESL_SQCACHE;


extern int   esl_sqfile_Open(const char *filename, int fmt, const char *env, ESL_SQFILE **ret_sqfp);
extern int   esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
extern int   esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
extern void  esl_sqfile_Close(ESL_SQFILE *sqfp);

extern int   esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);

extern int   esl_sqfile_Read      (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqfile_ReadInfo  (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqfile_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
extern int   esl_sqfile_ReadBlock (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target);
extern int   esl_sqfile_Parse     (char *buffer, int size, ESL_SQ *sq, int format);

extern int   esl_sqfile_Write     (FILE *fp, ESL_SQ *s, int format, int update);
extern int   esl_sqfile_Echo      (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);
extern int   esl_sqfile_WriteFasta(FILE *fp, ESL_SQ *sq, int save_offsets);

extern const char *esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp);
extern int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp);
extern int   esl_sqfile_IsAlignment(int fmt);
extern int   esl_sqfile_EncodeFormat(char *fmtstring);
extern char *esl_sqfile_DecodeFormat(int fmt);
extern int   esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset);

extern int   esl_sqfile_Ignore  (ESL_SQFILE *sqfp, const char *ignored);
extern int   esl_sqfile_Accept  (ESL_SQFILE *sqfp, const char *accepted);
extern int   esl_sqfile_AcceptAs(ESL_SQFILE *sqfp, const char *accepted, char readas);

extern int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
extern int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
extern int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);
extern int   esl_sqfile_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqfile_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqfile_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);

extern int   esl_sqfile_Cache(const ESL_ALPHABET *abc, const char *seqfile, int fmt, const char *env, ESL_SQCACHE **ret_sqcache);
extern void  esl_sqfile_Free(ESL_SQCACHE *sqcache);

#endif /*eslSQFILE_INCLUDED*/

