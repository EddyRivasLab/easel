/* Sequence clustering algorithms
 * 
 * SVN $Id$
 * SRE, Sun Nov  5 10:06:53 2006 [Janelia]
 */

#include <esl_config.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_msa.h>
#include <esl_vectorops.h>
#include <esl_distance.h>

#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
#include <ctype.h>
static double squid_distance(char *s1, char *s2);
static double squid_xdistance(ESL_ALPHABET *a, ESL_DSQ *x1, ESL_DSQ *x2);
#endif


/* Function:  esl_msacluster_SingleLinkage()
 * Synopsis:  Single linkage clustering by percent identity.
 * Incept:    SRE, Sun Nov  5 10:11:45 2006 [Janelia]
 *
 * Purpose:   Perform single link clustering of seqs in a
 *            sequence alignment. Any pair of sequences with
 *            percent identity $\geq$ <maxid> are linked (using
 *            the definition from the distance module).
 *            
 *            The resulting clustering is returned in the <ret_c>
 *            array, and the total number of clusters is in <ret_nc>.
 *            The <c[0..nseq-1]> array assigns a cluster index to each
 *            sequence. For example, <c[4] = 1> means that sequence 4
 *            is assigned to cluster 1. Cluster indices range from
 *            <0..nc-1>. 
 *
 *            Importantly, this algorithm runs in $O(N)$ memory, and
 *            produces one discrete clustering. Compare to
 *            <esl_tree_SingleLinkage()>, which requires an $O(N^2)$ 
 *            adjacency matrix, and produces a hierarchical clustering
 *            tree.
 *            
 *            The algorithm is $O(LN^2)$ time, for N sequences of
 *            length L. However, the worst case (no links at all) is
 *            unusual. More typically, time scales as about $LN \log
 *            N$. The best case scales as $LN$, when there is just one
 *            cluster in a completely connected graph.
 *            
 * Algorithm: 
 *           I don't know if this algorithm is published. I 
 *           haven't seen it in graph theory books, but that might
 *           be because it's so obvious that nobody's bothered.
 *           
 *           In brief, we're going to do a breadth-first search of the
 *           graph, and we're going to calculate links on the fly
 *           rather than precalculating them into a standard adjacency
 *           matrix.
 *           
 *           While working, we keep two stacks of maximum length N:
 *                a : list of vertices that are still unconnected.
 *                b : list of vertices that we've connected to 
 *                    in our current breadth level, but we haven't
 *                    yet tested for other connections to a.
 *           The current length (number of elements in) a and b are
 *           kept in na, nb.
 *                    
 *           We store our results in an array of length N:
 *                c : assigns each vertex to a component. for example
 *                    c[4] = 1 means that vertex 4 is in component 1.
 *                    nc is the number of components. Components
 *                    are numbered from 0 to nc-1. We return c and nc
 *                    to our caller.
 *                    
 *           The algorithm is:
 *           
 *           Initialisation: 
 *                a  <-- all the vertices
 *                na <-- N
 *                b  <-- empty set
 *                nb <-- 0
 *                nc <-- 0
 *                
 *           Then:
 *                while (a is not empty)
 *                  pop a vertex off a, push onto b
 *                  while (b is not empty)
 *                    pop vertex v off b
 *                    assign c[v] = nc
 *                    for each vertex w in a:
 *                       compare v,w. If w is linked to v, remove w
 *                       from a, push onto b.
 *                  nc++     
 *           q.e.d. 
 *
 * Args:      msa     - multiple alignment to cluster
 *            maxid   - pairwise identity threshold: cluster if $\geq$ <maxid>
 *            ret_c   - optRETURN: cluster assignments for each sequence
 *            ret_nc  - optRETURN: number of clusters        

 * Returns:   <eslOK> on success; the <ret_c[0..nseq-1]> array contains cluster
 *            indices <0..nc-1> assigned to each sequence, and <ret_nc> contains
 *            the number of clusters. The <ret_c> array is allocated here, and
 *            must be free'd by the caller. The input <msa> is unmodified.
 *            
 *            The caller may pass <NULL> for either <ret_c> or
 *            <ret_nc> if it is only interested in one of the two
 *            results.
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In these events, <ret_c> is set to <NULL>
 *            and <ret_nc> is set to 0, and the <msa> is unmodified.
 */
int
esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, int **ret_c, int *ret_nc)
{
  int  na, *a = NULL;           /* stack of available vertices */
  int  nb, *b = NULL;           /* stack of working vertices   */
  int *c = NULL;                /* array of results            */
  int  nc;			/* total number of components  */
  double pid;			/* percent pairwise identity   */
  int  v,w;			/* index of a working vertices */
  int  i;			/* loop counter */
  int  status;

  ESL_ALLOC(a, (sizeof(int) * msa->nseq));
  ESL_ALLOC(b, (sizeof(int) * msa->nseq));
  ESL_ALLOC(c, (sizeof(int) * msa->nseq));
  for (i = 0; i < msa->nseq; i++) a[i] = i;
  na = msa->nseq;
  nb = 0;
  nc = 0;

  while (na > 0)
    {
      v = a[na-1]; na--;	/* pop a vertex off a, */
      b[nb] = v;   nb++;	/* and push onto b     */
      while (nb > 0)
	{
	  v    = b[nb-1]; nb--;	/* pop vertex off b          */
	  c[v] = nc;		/* assign it to component nc */
	  for (i = na-1; i >= 0; i--)/* backwards, becase of deletion/swapping we do*/
	    {
#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
	      if (msa->flags & eslMSA_DIGITAL)
   	        pid = 1. - squid_xdistance(msa->abc, msa->ax[v], msa->ax[a[i]]);
	      else
		pid = 1. - squid_distance(msa->aseq[v], msa->aseq[a[i]]);
#else
	      if (msa->flags & eslMSA_DIGITAL)
		status = esl_dst_XPairId(msa->abc, msa->ax[v], msa->ax[a[i]], &pid, NULL, NULL);
	      else
		status = esl_dst_CPairId(msa->aseq[v], msa->aseq[a[i]], &pid, NULL, NULL);
	      if (status != eslOK) goto ERROR;
#endif
	      if (pid >= maxid) /* linked? */
		{			
		  w = a[i]; a[i] = a[na-1]; na--; /* delete w from a (note swap) */
		  b[nb] = w; nb++;                /* push w onto b */
		}
	    }
	}
      nc++;
    }

  /* Cleanup and return
   */
  free(a);
  free(b);
  if (ret_c  != NULL) *ret_c  = c; else free(c);
  if (ret_nc != NULL) *ret_nc = nc;
  return eslOK;

 ERROR:
  if (c != NULL) free(c);
  if (b != NULL) free(b);
  if (a != NULL) free(a);
  if (ret_c  != NULL) *ret_c  = NULL;
  if (ret_nc != NULL) *ret_nc = 0;
  return status;
}


/* 
 * When regression testing against squid, we have to replace
 * Easel's distance calculations with a simpler, (even) less robust 
 * calculation that squid did.
 */
#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
static double 
squid_distance(char *s1, char *s2)
{
  int diff  = 0;
  int valid = 0;

  for (; *s1 != '\0'; s1++, s2++)
    {
      if (!isalpha(*s1) || !isalpha(*s2)) continue;
      if (*s1 != *s2) diff++;
      valid++;
    }
  return (valid > 0 ? ((double) diff / (double) valid) : 0.0);
}
static double
squid_xdistance(ESL_ALPHABET *a, ESL_DSQ *x1, ESL_DSQ *x2)
{
  int diff  = 0;
  int valid = 0;

  for (; *x1 != eslDSQ_SENTINEL; x1++, x2++)
    {
      if (esl_abc_XIsGap(a, *x1) || esl_abc_XIsGap(a, *x2)) continue;
      if (*x1 != *x2) diff++;
      valid++;
    }
  return (valid > 0 ? ((double) diff / (double) valid) : 0.0);
}
#endif /* eslMSACLUSTER_REGRESSION || eslMSAWEIGHT_REGRESSION */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

