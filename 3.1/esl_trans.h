/* Translate DNA sequence into six frames, into individual ORFs.
 * 
 */

#ifndef eslTRANS_INCLUDED
#define eslTRANS_INCLUDED

#include "easel.h"
#include "esl_gencode.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_getopts.h"


/*****************************************************************
 * 1. A stateful structure, workstate_s, to support both ReadSeq and ReadWindow()
 *****************************************************************/

/* struct esl_trans_workstate_s
 *   keeps state in DNA sequence <sq>, allowing us to process a sequence
 *   either in a single gulp (using ReadSeq) or in overlapping windows 
 *   (using ReadWindow).
 *
 *   also contains one-time configuration information
 */
typedef struct esl_trans_workstate_s {
  /* stateful info (which may get updated with each new seq, strand, and/or window): */
  ESL_SQ *psq[3];     // Growing ORFs in each frame
  int8_t  in_orf[3];  // TRUE|FALSE: TRUE if we're growing an ORF in this frame
  int     apos;       // 1..L:  current nucleotide we're on (starting a codon) in <sq>
  int     frame;      // 0..2:  which frame <apos> is in
  int     codon;      // 0..63: Digitized codon for apos,apos+1,apos+2
  int     inval;      // 0..3:  how many apos increments we need to get past an ambiguous nucleotide
  int     is_revcomp; // TRUE|FALSE: TRUE if we're doing reverse complement strand
  int     orfcount;   // >=0:   How many ORFs we've processed so far
  ESL_SQ_BLOCK  *orf_block; // block of sequences to which to write ORFs
  
  /* one-time configuration information (from options) */
  int     do_watson;         // TRUE|FALSE:  TRUE if we translate the top strand
  int     do_crick;          // TRUE|FALSE:  TRUE if we translate the reverse complement strand
  int     using_initiators;  // TRUE|FALSE : TRUE if -m or -M, only valid initiators can start an ORF, and initiator codon always translates to Met
  int     minlen;            // >=0: minimum orf length that process_orf will deal with
  FILE   *outfp;             // default stdout: where to write output ORF data
  int     outformat;         // default eslSQFILE_FASTA: sqfile format to write ORFs in
} ESL_TRANS_WORKSTATE;

extern void esl_trans_WorkstateDestroy(ESL_TRANS_WORKSTATE *wrk);
extern ESL_TRANS_WORKSTATE * esl_trans_WorkstateCreate(ESL_GETOPTS *go, ESL_GENCODE *gcode);


/*****************************************************************
 *  2. Functions shared by the esl-translate's full or windowed reads style,
 *  as well as nhmmscant
 *****************************************************************/

extern int esl_trans_ProcessOrf(ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq);
extern void esl_trans_ProcessStart(ESL_GENCODE *gcode, ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq);
extern int esl_trans_ProcessPiece(ESL_GENCODE *gcode, ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq);
extern int esl_trans_ProcessEnd(ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq);


#endif /*eslTRANS_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
