/* Multiple sequence alignment file i/o
 */
#ifndef eslMSAFILE_INCLUDED
#define eslMSAFILE_INCLUDED

#include <stdio.h>

#include "esl_buffer.h"
#include "esl_msa.h"

#include "esl_alphabet.h"	/* AUGMENTATION: adds ability to use digital alphabet */
#include "esl_ssi.h"        	/* AUGMENTATION: adds ability to use SSI file indices */

/* Object: ESLX_MSAFILE
 * 
 * Defines an alignment file that's open for parsing.
 */
typedef struct {
  ESL_BUFFER         *bf;             /* input file/data being parsed                          */

  char               *line;	      /* line read from <bf> by <esl_msafile_GetLine()>        */
  esl_pos_t           n;	      /* length of line in bytes (line is not NUL-terminated)  */
  int64_t             linenumber;     /* input linenumber for diagnostics; -1 if we lose track */
  esl_pos_t           lineoffset;     /* offset of start of <line> in <bf> input               */

  int32_t             format;	      /* format of alignment file we're reading                */
  ESL_DSQ             inmap[128];     /* input map, 0..127                                     */
  const ESL_ALPHABET *abc;	      /* non-NULL if augmented and in digital mode             */
  ESL_SSI            *ssi;	      /* open SSI index; or NULL, if none or not augmented     */
  char                errmsg[eslERRBUFSIZE];   /* user-directed message for normal errors      */
} ESLX_MSAFILE;


/* Alignment file format codes.
 * Must coexist with sqio unaligned file format codes.
 * Rules:
 *     - 0 is an unknown/unassigned format 
 *     - <=100 reserved for unaligned formats
 *     - >100 reserved for aligned formats
 */
#define eslMSAFILE_UNKNOWN     0    /* unknown format                              */
#define eslMSAFILE_STOCKHOLM   101  /* Stockholm format, interleaved               */
#define eslMSAFILE_PFAM        102  /* Pfam/Rfam one-line-per-seq Stockholm format */
#define eslMSAFILE_A2M         103  /* UCSC SAM's fasta-like a2m format            */
#define eslMSAFILE_PSIBLAST    104  /* NCBI PSI-BLAST alignment format             */
#define eslMSAFILE_SELEX       105  /* old SELEX format (largely obsolete)         */
#define eslMSAFILE_AFA         106  /* aligned FASTA format                        */
#define eslMSAFILE_CLUSTAL     107  /* CLUSTAL format                              */
#define eslMSAFILE_CLUSTALLIKE 108  /* CLUSTAL-like formats (MUSCLE, PROBCONS)     */

/* 1. Opening/closing an ESLX_MSAFILE */
extern int   eslx_msafile_Open(ESL_ALPHABET **byp_abc, const char *msafile, int format, const char *env, ESLX_MSAFILE **ret_afp);
extern int   eslx_msafile_OpenMem(ESL_ALPHABET **byp_abc, char *p, esl_pos_t n, int format, ESLX_MSAFILE **ret_afp);
extern int   eslx_msafile_OpenBuffer(ESL_ALPHABET **byp_abc, ESL_BUFFER *bf, int format, ESLX_MSAFILE **ret_afp);
extern void  eslx_msafile_OpenFailure(ESLX_MSAFILE *afp, int status);
extern void  eslx_msafile_Close(ESLX_MSAFILE *afp);

/* 2. Utilities for different file formats */
extern int   eslx_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode);
extern int   eslx_msafile_EncodeFormat(char *fmtstring);
extern char *eslx_msafile_DecodeFormat(int fmt);

/* 3. Utilities for different alphabets */
#ifdef eslAUGMENT_ALPHABET
extern int eslx_msafile_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
#endif

/* 4. Reading an MSA from an ESLX_MSAFILE */
extern int  eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern void eslx_msafile_ReadFailure(ESLX_MSAFILE *afp, int status);

/* 5. Writing an MSA to a stream */
extern int eslx_msafile_Write(FILE *fp, ESL_MSA *msa, int fmt);

/* 6. Utilities for specific parsers */
extern int eslx_msafile_GetLine(ESLX_MSAFILE *afp, char **opt_p, esl_pos_t *opt_n);

#include "esl_msafile_a2m.h"
#include "esl_msafile_afa.h"
#include "esl_msafile_clustal.h"
#include "esl_msafile_psiblast.h"
#include "esl_msafile_selex.h"
#include "esl_msafile_stockholm.h"
#endif /*eslMSAFILE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
