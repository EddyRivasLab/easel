/* esl_distance.h
 * Distances between aligned sequences, including both
 * probabilistic evolutionary models and ad hoc measures.
 * 
 * SVN $Id$
 * SRE, Fri Apr 28 06:41:13 2006 [New York]
 */
#ifndef ESL_DISTANCE_INCLUDED
#define ESL_DISTANCE_INCLUDED

extern int esl_dst_CPairId(ESL_ALPHABET *abc, char *asq1, char *asq2, 
			   double *ret_pid, int *ret_nid, int *ret_n);
extern int esl_dst_XPairId(ESL_ALPHABET *abc, ESL_DSQ *ax1, ESL_DSQ *ax2, 
			   double *ret_pid, int *ret_nid, int *ret_n);
extern int esl_dst_CPairIdMx(ESL_ALPHABET *abc, char **as, int N, ESL_DMATRIX **ret_S);
extern int esl_dst_XPairIdMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S);
extern int esl_dst_CDiffMx(ESL_ALPHABET *abc, char **as, int N, ESL_DMATRIX **ret_D);
extern int esl_dst_XDiffMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D);
extern int esl_dst_CJukesCantor(ESL_ALPHABET *abc, char *as1, char *as2, 
				double *ret_distance, double *ret_variance);
extern int esl_dst_XJukesCantor(ESL_ALPHABET *abc, ESL_DSQ *ax, ESL_DSQ *ay, 
				double *ret_distance, double *ret_variance);
extern int esl_dst_CJukesCantorMx(ESL_ALPHABET *abc, char **aseq, int nseq, 
				  ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V);
extern int esl_dst_XJukesCantorMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int nseq, 
				  ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V);

#endif /*ESL_DISTANCE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
