/* Find subset of sequences such that no pair is >t% identical (independent set)
 * Or find a pair of disjoint subsets X and Y such that no pair of sequences
 * one in X and one in Y are >t% identical (bipartite independent pair)
 *
 * Table of contents:
 *    1. Bipartite independent pair algorithms (Random, Cobalt, Blue)
 *    2. Independent set algorithms (Cobalt, Blue)
 *    3. Internal functions, interface to the clustering API
 *    4. Unit tests
 *    5. Test driver
 */
#include "esl_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cluster.h"
#include "esl_distance.h"
#include "esl_msa.h"

#include "esl_random.h"
#include "esl_msa_iset.h"
#include "esl_iset.h"

/* These functions are going to get defined in an internal regression
 * testing section further below:
 */
#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
#include <ctype.h>
static double squid_distance(char *s1, char *s2);
static double squid_xdistance(ESL_ALPHABET *a, ESL_DSQ *x1, ESL_DSQ *x2);
#endif

/* These functions will define linkage between a pair of text or
 *  digital aseq's:
 */
static int msacluster_clinkage(const void *v1, const void *v2, const void *p, int *ret_link);
static int msacluster_xlinkage(const void *v1, const void *v2, const void *p, int *ret_link);

/* In digital mode, we'll need to pass the clustering routine two parameters -
 * %id threshold and alphabet ptr - so make a structure that bundles them.
 */
struct msa_param_s {
  double        maxid;
  ESL_ALPHABET *abc;
};

/*****************************************************************
 * 1. Independent set algorithms
 *****************************************************************/

/* Function:  esl_msa_iset_Cobalt(), 
 * Synopsis:  Produces a independent set by a greedy algorithm with a random
 *            order
 *
 * Incept:    SNP, Oct 16 2020 
 *
 * Purpose:   Produce an independent set. For algorithm details, see
 *            description of esl_bi_iset_Cobalt in esl_iset.c.
 *
 * Args:      msa     - multiple alignment to find independent set within
 *            maxid   - pairwise identity threshold: no pair can be $\geq$ <maxid>
 *            opt_c   - optRETURN: set assignments for each sequence, [0..nseq-1]
 *            r       - source of randomness
 *
 * Returns:   <eslOK> on success; the <opt_c[0..nseq-1]> array contains
 *            set indices: 1 if sequence in iset, 0 if sequence not in iset.
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In either case, <opt_c> and <opt_nin> are set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_iset_Cobalt(const ESL_MSA *msa, double maxid,
           int **opt_c, int **opt_nin, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 2);
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  if (! (msa->flags & eslMSA_DIGITAL))
    status = esl_iset_Cobalt((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
               msacluster_clinkage, (void *) &maxid,
               workspace, assignment, r);
  else {
    param.maxid = maxid;
    param.abc   = msa->abc;
//    printf("calling esl_iset_Cobalt in else\n");
    status = esl_iset_Cobalt((void *) msa->ax, (size_t) msa->nseq, sizeof(ESL_DSQ *),
               msacluster_xlinkage, (void *) &param,
               workspace, assignment, r);
  }
  if (status != eslOK) goto ERROR;




  /* cleanup and return */
  free(workspace);
  if (opt_c  != NULL) *opt_c  = assignment; else free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (opt_c  != NULL) *opt_c  = NULL;
  return status;
}


/* Function:  esl_msa_iset_Blue(), 
 * Synopsis:  Produces an independent set by a multi-round election process
 *
 * Incept:    SNP, Oct 16 2020 
 *
 * Purpose:   Produce a bipartite independent pair. For algorithm details, see
 *            description of esl_iset_Blue in esl_iset.c.
 *
 * Args:      msa     - multiple alignment to find independent set within
 *            maxid   - pairwise identity threshold: no pair can be $\geq$ <maxid>
 *            opt_c   - optRETURN: set assignments for each sequence, [0..nseq-1]
 *            r       - source of randomness
 *
 * Returns:   <eslOK> on success; the <opt_c[0..nseq-1]> array contains
 *            set indices: 1 if sequence in iset, 0 if sequence not in iset.
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In either case, <opt_c> and <opt_nin> are set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_iset_Blue(const ESL_MSA *msa, double maxid,
			     int **opt_c, int **opt_nin, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 4);
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  if (! (msa->flags & eslMSA_DIGITAL))
    status = esl_iset_Blue((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
				       msacluster_clinkage, (void *) &maxid,
				       workspace, assignment, r);
  else {
    param.maxid = maxid;
    param.abc   = msa->abc;
//    printf("calling esl_iset_Cobalt in else\n");
    status = esl_iset_Blue((void *) msa->ax, (size_t) msa->nseq, sizeof(ESL_DSQ *),
				       msacluster_xlinkage, (void *) &param,
				       workspace, assignment, r);
  }
  if (status != eslOK) goto ERROR;




  /* cleanup and return */
  free(workspace);
  if (opt_c  != NULL) *opt_c  = assignment; else free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (opt_c  != NULL) *opt_c  = NULL;
  return status;
}


/*****************************************************************
 * 2. Bipartite independent pair algorithms
 *****************************************************************/

/* Function:  esl_msa_bi_iset_Random(), 
 * Synopsis:  Random biparitite indpendent pair algorithm
 *
 * Incept:    SNP, Oct 16 2020 
 *
 * Purpose:   Produce a bipartite independent pair, where one of the sets of 
 *            of the pair is chosen independently at random. For algorithm details, 
 *            see description of esl_bi_iset_Random in esl_iset.c.
 *
 * Args:      msa     - multiple alignment to find independent pair within
 *            maxid   - pairwise identity threshold: no pair can be $\geq$ <maxid>
 *            opt_c   - optRETURN: set assignments for each sequence, [0..nseq-1]
 *            t_prob  - each sequence is included in the random set independently
 *                      with probability t_prob
 *            r       - source of randomness
 *
 * Returns:   <eslOK> on success; the <opt_c[0..nseq-1]> array contains
 *            set indices: 
 *            0 - sequence not in bipartite independent pair
 *            1 - sequence in random set of bipartite independent pair
 *            2 - sequence in other set of bipartite independent pair  
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In either case, <opt_c> and <opt_nin> are set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_bi_iset_Random(const ESL_MSA *msa, double maxid,
           int **opt_c, int **opt_nin, ESL_RANDOMNESS *r, double t_prob)

{
  int   status;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  if (! (msa->flags & eslMSA_DIGITAL))
    status = esl_bi_iset_Random((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
               msacluster_clinkage, (void *) &maxid,
                assignment, r, t_prob);
  else {
    param.maxid = maxid;
    param.abc   = msa->abc;
//    printf("calling esl_iset_Cobalt in else\n");
    status = esl_bi_iset_Random((void *) msa->ax, (size_t) msa->nseq, sizeof(ESL_DSQ *),
               msacluster_xlinkage, (void *) &param,
               assignment, r, t_prob);
  }
  if (status != eslOK) goto ERROR;

  /* cleanup and return */
  if (opt_c  != NULL) *opt_c  = assignment; else free(assignment);
  return eslOK;

 ERROR:
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (opt_c  != NULL) *opt_c  = NULL;
  return status;
}

/* Function:  esl_msa_bi_iset_Cobalt(), 
 * Synopsis:  Produces a bipartite independent pair by a greedy algorithm with 
 *            a random order
 *
 * Incept:    SNP, Oct 16 2020 
 *
 * Purpose:   Produce a bipartite independent pair. For algorithm details, see
 *            description of esl_bi_iset_Cobalt in esl_iset.c.
 *
 * Args:      msa     - multiple alignment to find independent pair within
 *            maxid   - pairwise identity threshold: no pair can be $\geq$ <maxid>
 *            opt_c   - optRETURN: set assignments for each sequence, [0..nseq-1]
 *            r       - source of randomness
 *
 * Returns:   <eslOK> on success; the <opt_c[0..nseq-1]> array contains
 *            set indices: 
 *            0 - sequence not in bipartite independent pair
 *            1 - sequence in one set of bipartite independent pair
 *            2 - sequence in other set of bipartite independent pair
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In either case, <opt_c> and <opt_nin> are set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */


int
esl_msa_bi_iset_Cobalt(const ESL_MSA *msa, double maxid,
           int **opt_c, int **opt_nin, int *ret_larger, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  int   larger;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq*3);
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  if (! (msa->flags & eslMSA_DIGITAL))
    status = esl_bi_iset_Cobalt((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
               msacluster_clinkage, (void *) &maxid,
               workspace, assignment, &larger, r);
  else {
    param.maxid = maxid;
    param.abc   = msa->abc;
//    printf("calling esl_iset_Cobalt in else\n");
    status = esl_bi_iset_Cobalt((void *) msa->ax, (size_t) msa->nseq, sizeof(ESL_DSQ *),
               msacluster_xlinkage, (void *) &param,
               workspace, assignment, &larger, r);
  }
  if (status != eslOK) goto ERROR;

  /* cleanup and return */
  free(workspace);
  if (ret_larger != NULL) *ret_larger = larger;
  if (opt_c  != NULL) *opt_c  = assignment; else free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (ret_larger != NULL) *ret_larger = 0;
  if (opt_c  != NULL) *opt_c  = NULL;
  return status;
}

/* Function:  esl_msa_bi_iset_Blue(), 
 * Synopsis:  Produces a bipartite independent pair by a multi-round election 
 *            process
 *
 * Incept:    SNP, Oct 16 2020 
 *
 * Purpose:   Produce a bipartite independent pair. For algorithm details, see
 *            description of esl_bi_iset_Blue in esl_iset.c.
 *
 * Args:      msa     - multiple alignment to find independent pair within
 *            maxid   - pairwise identity threshold: no pair can be $\geq$ <maxid>
 *            opt_c   - optRETURN: set assignments for each sequence, [0..nseq-1]
 *            r       - source of randomness
 *
 * Returns:   <eslOK> on success; the <opt_c[0..nseq-1]> array contains
 *            set indices: 
 *            0 - sequence not in bipartite independent pair
 *            1 - sequence in one set of bipartite independent pair
 *            2 - sequence in other set of bipartite independent pair
 *
 * Throws:    <eslEMEM> on allocation failure, and <eslEINVAL> if a pairwise
 *            comparison is invalid (which means the MSA is corrupted, so it
 *            shouldn't happen). In either case, <opt_c> and <opt_nin> are set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_bi_iset_Blue(const ESL_MSA *msa, double maxid,
           int **opt_c, int **opt_nin, int *ret_larger, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  int   larger;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq*5);
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  if (! (msa->flags & eslMSA_DIGITAL))
    status = esl_bi_iset_Blue((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
               msacluster_clinkage, (void *) &maxid,
               workspace, assignment, &larger, r);
  else {
    param.maxid = maxid;
    param.abc   = msa->abc;
//    printf("calling esl_iset_Cobalt in else\n");
    status = esl_bi_iset_Blue((void *) msa->ax, (size_t) msa->nseq, sizeof(ESL_DSQ *),
               msacluster_xlinkage, (void *) &param,
               workspace, assignment, &larger, r);
  }
  if (status != eslOK) goto ERROR;

  /* cleanup and return */
  free(workspace);
  if (ret_larger != NULL) *ret_larger = larger;
  if (opt_c  != NULL) *opt_c  = assignment; else free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (ret_larger != NULL) *ret_larger = 0;
  if (opt_c  != NULL) *opt_c  = NULL;
  return status;
}



/*****************************************************************
 * 3. Internal functions, interface to the clustering API
 *****************************************************************/

/* Definition of %id linkage in text-mode aligned seqs (>= maxid): */
static int
msacluster_clinkage(const void *v1, const void *v2, const void *p, int *ret_link)
{
  char  *as1   = *(char **) v1;
  char  *as2   = *(char **) v2;
  double maxid = *(double *) p;
  double pid;
  int    status = eslOK;

#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
  pid = 1. - squid_distance(as1, as2);
#else
  if ((status = esl_dst_CPairId(as1, as2, &pid, NULL, NULL)) != eslOK) return status;
#endif

  *ret_link = (pid >= maxid ? TRUE : FALSE);
  return status;
}

/* Definition of % id linkage in digital aligned seqs (>= maxid) */
static int
msacluster_xlinkage(const void *v1, const void *v2, const void *p, int *ret_link)
{
  ESL_DSQ *ax1              = *(ESL_DSQ **) v1;
  ESL_DSQ *ax2              = *(ESL_DSQ **) v2;
  struct msa_param_s *param = (struct msa_param_s *) p;
  double   pid;
  int      status = eslOK;

#if defined(eslMSACLUSTER_REGRESSION) || defined(eslMSAWEIGHT_REGRESSION)
  pid = 1. - squid_xdistance(param->abc, ax1, ax2);
#else
  if ( (status = esl_dst_XPairId(param->abc, ax1, ax2, &pid, NULL, NULL)) != eslOK) return status;
#endif

  *ret_link = (pid >= param->maxid ? TRUE : FALSE);
  return status;
}









