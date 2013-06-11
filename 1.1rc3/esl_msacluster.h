/* Clustering sequences in an MSA by % identity.
 * 
 * SRE, Sun Nov  5 10:08:14 2006 [Janelia]
 * SVN $Id$
 * SVN $URL$
 */
#ifndef eslMSACLUSTER_INCLUDED
#define eslMSACLUSTER_INCLUDED

extern int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, 
					int **opt_c, int **opt_nin, int *opt_nc);

#endif /*eslMSACLUSTER_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
