/* esl_msacluster.h
 * Sequence clustering algorithms
 * 
 * SVN $Id$
 * SRE, Sun Nov  5 10:08:14 2006 [Janelia]
 */

#ifndef ESL_MSACLUSTER_INCLUDED
#define ESL_MSACLUSTER_INCLUDED

extern int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, int **ret_c, int *ret_nc);

#endif /*ESL_MSACLUSTER_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
