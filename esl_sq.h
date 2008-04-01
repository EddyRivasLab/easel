/* A sequence.
 * 
 * SRE, Mon Mar 31 17:03:51 2008 [Janelia]
 * SVN $Id$
 */
#ifndef ESL_SQ_INCLUDED
#define ESL_SQ_INCLUDED

#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_MSA
#include "esl_msa.h"
#endif

/* ESL_SQ - a biosequence
 * 
 * Can be either in text mode <seq>, or digital mode <dsq>. 
 * One of them has to be NULL, and the other contains the data.
 * 
 * When in text mode, <ss> and <seq> can hold up to <n=salloc-1>
 * residues and a terminal '\0', and both are indexed <0..n-1>.
 * 
 * When in digital mode, <ss> and <dsq> can hold up to <n=salloc-2>
 * residues; both are indexed <1..n>, and positions 0 and n+1 are
 * sentinel bytes. The digital sequence <dsq> uses <eslDSQ_SENTINEL>
 * as its sentinels; as a hack, <ss> uses '\0' as sentinels.  This
 * means that <sq->ss+1> is a valid NUL-terminated C string, but
 * <sq->ss> itself would be a string of length 0 because of the
 * leading NUL sentinel.
 * 
 * To save on allocation calls, the structure is designed to be reused
 * for subsequent sequences, rather than free'd and reallocated -
 * thus, we keep track of the allocated sizes of all the strings.
 * 
 * Notes on when we need to reallocate:
 *    - In a text mode sequence (seq 0..n-1), byte salloc-1 is
 *      reserved for the NUL, so the sequence is full when
 *      n == salloc-1.
 *          
 *    - In a digital mode sequence (dsq 1..n), bytes 0 and salloc-1
 *      are reserved for sentinel bytes, so the reallocation condition
 *      is when n == salloc-2.
 *      
 * At least for now, the only way to set the <ss> structure annotation
 * field is by a CreateFrom(), by extraction from an MSA, or manually
 * (by directly allocating a sequence's <ss> field).
 */
typedef struct {
  /*::cexcerpt::sqio_sq::begin::*/
  char    *name;           /* name ("\0" if no name)                           */
  char    *acc;            /* optional accession ("\0" if no accession)        */
  char    *desc;           /* description ("\0" if no description)             */
  char    *seq;            /* sequence [0..n-1], or NULL if digital            */
  ESL_DSQ *dsq;            /* digitized sequence [1..n], or NULL if text       */
  char    *ss;             /* optional sec structure [0..n-1], [1..n], or NULL */
  int      n;              /* length of seq (or dsq) and ss                    */
  /*::cexcerpt::sqio_sq::end::*/

  int      nalloc;         /* allocated length of name                         */
  int      aalloc;         /* allocated length of accession                    */
  int      dalloc;         /* allocated length of description                  */
  int      salloc;         /* alloc for seq or dsq, and ss if present          */

  off_t    roff;	   /* record offset (start of record); -1 if none      */
  off_t    doff;	   /* data offset (start of sequence data); -1 if none */
#ifdef eslAUGMENT_ALPHABET
  const ESL_ALPHABET *abc; /* reference to the alphabet for <dsq>              */
#endif
} ESL_SQ;


/* These control default initial allocation sizes in an ESL_SQ.     */
#define eslSQ_NAMECHUNK   32	/* allocation unit for name         */
#define eslSQ_ACCCHUNK    32	/* allocation unit for accession    */
#define eslSQ_DESCCHUNK  128	/* allocation unit for description  */
#define eslSQ_SEQCHUNK   256	/* allocation unit for seqs         */


extern ESL_SQ *esl_sq_Create(void);
extern ESL_SQ *esl_sq_CreateFrom(const char *name, const char *seq,
				 const char *desc, const char *acc, const char *ss);
extern int     esl_sq_Grow  (ESL_SQ *sq, int *ret_nsafe);
extern int     esl_sq_GrowTo(ESL_SQ *sq, int  n);
extern int     esl_sq_Copy(const ESL_SQ *src, ESL_SQ *dst);
extern int     esl_sq_Reuse  (ESL_SQ *sq);
extern void    esl_sq_Destroy(ESL_SQ *sq);

extern int     esl_sq_SetName     (ESL_SQ *sq, char *name, ...);
extern int     esl_sq_SetAccession(ESL_SQ *sq, char *acc,  ...);
extern int     esl_sq_SetDesc     (ESL_SQ *sq, char *desc, ...);
extern int     esl_sq_CAddResidue (ESL_SQ *sq, char c);

#ifdef eslAUGMENT_ALPHABET
extern ESL_SQ *esl_sq_CreateDigital(const ESL_ALPHABET *abc);
extern ESL_SQ *esl_sq_CreateDigitalFrom(const ESL_ALPHABET *abc, const char *name, const ESL_DSQ *dsq, 
					int L, const char *desc, const char *acc,  const char *ss);
extern int     esl_sq_Digitize(const ESL_ALPHABET *abc, ESL_SQ *sq);
extern int     esl_sq_Textize(ESL_SQ *sq);
extern int     esl_sq_GuessAlphabet(ESL_SQ *sq, int *ret_type);
extern int     esl_sq_XAddResidue(ESL_SQ *sq, ESL_DSQ x);
#endif

#ifdef eslAUGMENT_MSA
extern int     esl_sq_GetFromMSA  (const ESL_MSA *msa, int which, ESL_SQ *sq);
extern int     esl_sq_FetchFromMSA(const ESL_MSA *msa, int which, ESL_SQ **ret_sq);
#endif

#endif /*!ESL_SQ_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
