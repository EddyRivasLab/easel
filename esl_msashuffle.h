/* Shuffling or bootstrapping multiple sequence alignments.
 * 
 * SRE, Tue Jan 22 09:18:09 2008 [Market Street Cafe, Leesburg]
 * SVN $Id$
 */
#ifndef ESL_MSASHUFFLE_INCLUDED
#define ESL_MSASHUFFLE_INCLUDED

#include "esl_random.h"
#include "esl_alphabet.h"

extern int esl_msashuffle_Shuffle  (ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *shuf);
extern int esl_msashuffle_Bootstrap(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *bootsample);

extern int esl_msashuffle_XQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, ESL_DSQ *x, ESL_DSQ *y, ESL_DSQ *xs, ESL_DSQ *ys);
extern int esl_msashuffle_CQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, char    *x, char    *y, char    *xs, char    *ys);

#endif /*ESL_MSASHUFFLE_INCLUDED*/
