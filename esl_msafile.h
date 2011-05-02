/* Multiple sequence alignment file i/o
 */
#ifndef eslMSAFILE_INCLUDED
#define eslMSAFILE_INCLUDED

#include <stdio.h>

#include "esl_buffer.h"
#include "esl_msa.h"

#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"	/* adds ability to use digital alphabet */
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"        	/* adds ability to use SSI file indices */
#endif

/* Object: ESLX_MSAFILE
 * 
 * Defines an alignment file that's open for parsing.
 */
typedef struct {
  ESL_BUFFER *bf;                     /* input file/data being parsed                          */
  int64_t     linenumber;             /* input linenumber for diagnostics; -1 if we lose track */
  int32_t     format;		      /* format of alignment file we're reading                */
  ESL_DSQ     inmap[128];	      /* input map, 0..127                                     */

#if defined(eslAUGMENT_ALPHABET)      /* AUGMENTATION (alphabet): digitized input              */
  const ESL_ALPHABET *abc;	      /* non-NULL if augmented and in digital mode             */
#else
  void               *abc;	      /* NULL if not digital mode, or unaugmented              */
#endif
#if defined(eslAUGMENT_SSI)	      /* AUGMENTATION: SSI indexing of an MSA db               */
  ESL_SSI *ssi;		              /* open SSI index file; or NULL, if none.                */
#else
  void    *ssi;
#endif
  char    errmsg[eslERRBUFSIZE];
} ESLX_MSAFILE;


/* Alignment file format codes.
 * Must coexist with sqio unaligned file format codes.
 * Rules:
 *     - 0 is an unknown/unassigned format 
 *     - <=100 reserved for unaligned formats
 *     - >100 reserved for aligned formats
 */
#define eslMSAFILE_UNKNOWN   0	  /* unknown format                              */
#define eslMSAFILE_STOCKHOLM 101  /* Stockholm format, interleaved               */
#define eslMSAFILE_PFAM      102  /* Pfam/Rfam one-line-per-seq Stockholm format */
#define eslMSAFILE_A2M       103  /* UCSC SAM's fasta-like a2m format            */
#define eslMSAFILE_PSIBLAST  104  /* NCBI PSI-BLAST alignment format             */
#define eslMSAFILE_SELEX     105  /* old SELEX format (largely obsolete)         */
#define eslMSAFILE_AFA       106  /* aligned FASTA format                        */
#define eslMSAFILE_CLUSTAL   107  /* CLUSTAL format                              */
#define eslMSAFILE_MUSCLE    108  /* MUSCLE format (essentially CLUSTAL)         */

extern int  eslx_msafile_Open(const char *msafile, int format, const char *env, ESLX_MSAFILE **ret_afp);
extern void eslx_msafile_OpenFailure(ESLX_MSAFILE *afp, int status);
extern void eslx_msafile_ReadFailure(ESLX_MSAFILE *afp, int status);
extern int  eslx_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode);
extern void eslx_msafile_Close(ESLX_MSAFILE *afp);


#endif /*eslMSAFILE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
