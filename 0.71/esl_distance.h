/* esl_distance.h
 * Distances between aligned sequences, including both
 * probabilistic evolutionary models and ad hoc measures.
 * 
 * SVN $Id$
 * SRE, Fri Apr 28 06:41:13 2006 [New York]
 */
#ifndef ESL_DISTANCE_INCLUDED
#define ESL_DISTANCE_INCLUDED

#include <easel.h>		/* ESL_DSQ declaration      */
#ifdef eslAUGMENT_ALPHABET
#include <esl_alphabet.h>	/* ESL_ALPHABET declaration */
#endif
#ifdef eslAUGMENT_DMATRIX
#include <esl_dmatrix.h>	/* ESL_DMATRIX declaration  */
#endif


/* Basic pairwise distances for two aligned text sequences.
 */
extern int esl_dst_CPairId(char *asq1, char *asq2, 
			   double *ret_pid, int *ret_nid, int *ret_n);
extern int esl_dst_CJukesCantor(int K, char *as1, char *as2, 
				double *ret_distance, double *ret_variance);

/* Pairwise distance calculations on aligned digital sequences.
 */
#ifdef eslAUGMENT_ALPHABET
extern int esl_dst_XPairId(ESL_ALPHABET *abc, ESL_DSQ *ax1, ESL_DSQ *ax2, 
			   double *ret_pid, int *ret_nid, int *ret_n);
extern int esl_dst_XJukesCantor(ESL_ALPHABET *abc, ESL_DSQ *ax, ESL_DSQ *ay, 
				double *ret_distance, double *ret_variance);
#endif


/* Distance matrices on aligned text sequences
 */
#ifdef eslAUGMENT_DMATRIX
extern int esl_dst_CPairIdMx     (char **as, int N, ESL_DMATRIX **ret_S);
extern int esl_dst_CDiffMx       (char **as, int N, ESL_DMATRIX **ret_D);
extern int esl_dst_CJukesCantorMx(int K, char **as, int N, ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V);
#endif

/* Distance matrices on aligned digital sequences
 */
#if defined(eslAUGMENT_DMATRIX) && defined(eslAUGMENT_ALPHABET)
extern int esl_dst_XPairIdMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S);
extern int esl_dst_XDiffMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D);


extern int esl_dst_XJukesCantorMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int nseq, 
				  ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V);
#endif




#endif /*ESL_DISTANCE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
