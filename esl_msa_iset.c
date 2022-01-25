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
 *            shouldn't happen). In either case, <opt_c> is set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_iset_Cobalt(const ESL_MSA *msa, double maxid,
           int **opt_c, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;
  int allocated_assignment =0;
  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 2);
  if(opt_c != NULL){
    assignment = *opt_c;
  }
  else{
    ESL_ALLOC(assignment, sizeof(int) * msa->nseq);
    allocated_assignment = 1;
  }
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
  if(allocated_assignment){
    free(assignment);
  }
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (allocated_assignment) free(assignment);
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
 *            shouldn't happen). In either case, <opt_c> is set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_iset_Blue(const ESL_MSA *msa, double maxid,
			     int **opt_c, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;
  int allocated_assignment = 0;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 4);
  if(opt_c != NULL){
    assignment = *opt_c;
  }
  else{
    ESL_ALLOC(assignment, sizeof(int) * msa->nseq);
    allocated_assignment = 1;
  }

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
  if (allocated_assignment) free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (allocated_assignment) free(assignment);
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
 *            shouldn't happen). In either case, <opt_c> is set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_bi_iset_Random(const ESL_MSA *msa, double maxid,
           int **opt_c, ESL_RANDOMNESS *r, double t_prob)

{
  int   status;
  int  *assignment = NULL;
  int  *nin        = NULL;
  struct msa_param_s param;

  /* Allocations */
  int allocated_assignment =0;
  /* Allocations */
  if(opt_c != NULL){
    assignment = *opt_c;
  }
  else{
    ESL_ALLOC(assignment, sizeof(int) * msa->nseq);
    allocated_assignment = 1;
  }

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
  if (allocated_assignment) free(assignment);
  return eslOK;

 ERROR:
  if (allocated_assignment) free(assignment);
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
 *            shouldn't happen). In either case, <opt_c> is set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */


int
esl_msa_bi_iset_Cobalt(const ESL_MSA *msa, double maxid,
           int **opt_c, int *ret_larger, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  int   larger;
  struct msa_param_s param;

  /* Allocations */
  int allocated_assignment =0;
  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 3);
  if(opt_c != NULL){
    assignment = *opt_c;
  }
  else{
    ESL_ALLOC(assignment, sizeof(int) * msa->nseq);
    allocated_assignment = 1;
  }

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
  if (allocated_assignment) free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (allocated_assignment) free(assignment);
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
 *            shouldn't happen). In either case, <opt_c> is set to <NULL>
 *            and <opt_nc> is set to 0, and the <msa> is unmodified.
 */

int
esl_msa_bi_iset_Blue(const ESL_MSA *msa, double maxid,
           int **opt_c, int *ret_larger, ESL_RANDOMNESS *r)

{
  int   status;
  int  *workspace  = NULL;
  int  *assignment = NULL;
  int  *nin        = NULL;
  int   larger;
  struct msa_param_s param;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq*5);
  int allocated_assignment =0;
  if(opt_c != NULL){
    assignment = *opt_c;
  }
  else{
    ESL_ALLOC(assignment, sizeof(int) * msa->nseq);
    allocated_assignment = 1;
  }

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
  if (allocated_assignment) free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (allocated_assignment) free(assignment);
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



/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef eslMSA_ISET_TESTDRIVE

#include "esl_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_iset.h"

// Checking functions heavily cribbed from check_iset and check_msa_iset in esl_iset.c
/* Function: check_msa_iset()
 * Synopsis: Verify that a subset of vertices is an independent set
 * Incept:   NPC 1/24/22
 *
 * Purpose:  Given a subset of vertices, verify whether they form an 
 *           independent set 
 *
 * Args:     digital     - nonzero if source MSA was digital, zero otherwise
 *           base        - pointer to array of n fixed-size vertices in graph
 *           n           - number of vertices
 *           size        - size of each vertex element
 *           maxid       - maximum identity used in the partitioning
 *           abc         - alphabet used by the MSA. Only requiced if MSA digital
 *           assignments - array of 0/1s; 1 indicates a vertex is in the subset, 0 indicates
 *                         vertex not in the subset
 *
 * Returns:   <eslOK> if the subset is an independent set, eslFAIL if not
 *
 */

static int check_msa_iset(int digital, void *base, size_t n, double maxid, ESL_ALPHABET *abc, int *assignments)
{
   int i,j;
   int status;
   int do_link;
   struct msa_param_s param;

   if (!digital){ // ASCII MSA
    size_t size = sizeof(char *);
    for (i = 0; i < n; i++){
      for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==1) {
          if ((status = msacluster_clinkage( (char *) base + j*size, (char *) base + i*size, &maxid, &do_link)) != eslOK) goto ERROR;
	        if (do_link){
              return eslFAIL;
           }
        }
     }
    }
  }
  else{
    size_t size = sizeof(ESL_DSQ *);
    param.maxid = maxid;
    param.abc   = abc;
    for (i = 0; i < n; i++){
      for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==1) {
          if ((status = msacluster_xlinkage( (char *) base + j*size, (char *) base + i*size, &param, &do_link)) != eslOK) goto ERROR;
	        if (do_link){
              return eslFAIL;
           }
        }
     }
    }
  }
  return eslOK;

  ERROR:
  return eslFAIL; 
}

/* Function: check_msa_bi_iset()
 * Synopsis: Verify that a pair of subsets of vertices form a bipartite independent
 *           pair
 *
 * Incept:   NPC 1/24/22
 *  
 * Purpose:  Given a pair of disjoint subsets of vertices, verify that the pair is a 
 *           bipartite independent pair 
 *
 * Args:     digital     - nonzero if the source MSA was digital, 0 otherwise 
 *           base        - pointer to array of n fixed-size vertices in graph
 *           n           - number of vertices
 *           size        - size of each vertex element
 *           maxid       - maximum identity used in the partitioning
 *           abc         - alphabet used by the MSA. Only requiced if MSA digital
 *           assignments - array of 0/1/2s; 1 indicates the vertex is in one subset, 2 indicates a
 *                         the vertex in in the other subset, 0 indicates the vertex is in neither
 *
 * Returns:   <eslOK> if the pair forms a bipartite independent pair 
 *
 * Throws:   esl_fatal error if not a bipartite independent set pair
 */

static int check_msa_bi_iset( int digital, void *base, size_t n, double maxid, ESL_ALPHABET *abc, int *assignments)
{
   struct msa_param_s param;
   int i,j;
   int status;
   int do_link;
  if(!digital){
    size_t size = sizeof(char *);
    for (i = 0; i < n; i++){
      for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==2) {
           if ((status = msacluster_clinkage( (char *) base + j*size, (char *) base + i*size, &maxid, &do_link)) != eslOK) goto ERROR;
	         if (do_link){
              return(eslFAIL);
           }
        }
        else if (assignments[i]==2 && assignments[j]==1) {
          if ((status = msacluster_clinkage( (char *) base + j*size, (char *) base + i*size, &maxid, &do_link)) != eslOK) goto ERROR;
	        if (do_link){
              return(eslFAIL);
          }
        }
     }
   }
  }
  else{
    size_t size = sizeof(ESL_DSQ *);
    param.maxid = maxid;
    param.abc   = abc;
    for (i = 0; i < n; i++){
      for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==2) {
           if ((status = msacluster_xlinkage( (char *) base + j*size, (char *) base + i*size, &param, &do_link)) != eslOK) goto ERROR;
	         if (do_link){
              return(eslFAIL);
           }
        }
        else if (assignments[i]==2 && assignments[j]==1) {
          if ((status = msacluster_xlinkage( (char *) base + j*size, (char *) base + i*size, &param, &do_link)) != eslOK) goto ERROR;
	        if (do_link){
              return(eslFAIL);
          }
        }
     }
   }
  }
   return eslOK;

   ERROR:
   return status;

}

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for iset module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  int *assignments;
  ESL_MSA        *msa     = esl_msa_CreateFromString("\
# STOCKHOLM 1.0\n\
\n\
seq0  AAAAAAAAAA\n\
seq1  AAAAAAAAAA\n\
seq2  AAAAAAAAAC\n\
seq3  AAAAAAAADD\n\
seq4  AAAAAAAEEE\n\
seq5  AAAAAAFFFF\n\
seq6  AAAAAGGGGG\n\
seq7  AAAAHHHHHH\n\
seq8  AAAIIIIIII\n\
seq9  AAKKKKKKKK\n\
seq10 ALLLLLLLLL\n\
seq11 MMMMMMMMMM\n\
//",   eslMSAFILE_STOCKHOLM);
  int status;
  ESL_RANDOMNESS *r = NULL;
  r=esl_randomness_Create(0);
  ESL_ALLOC(assignments, 12 * sizeof(int));  //must be = # of sequences in alignment
  // Make digital copy of the msa
  ESL_MSA *msa2 = esl_msa_Create(12, 10);
  esl_msa_Copy(msa, msa2);
  esl_msa_Digitize(abc, msa2, NULL);
  int larger, larger_set;
  //Test 1: msa_iset_Cobalt on ASCII MSA
  if(esl_msa_iset_Cobalt(msa, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_iset(0, msa->aseq, 12, sizeof(char *), abc, assignments) != eslOK){
    return eslFAIL;
  }

  //Test2: msa_iset_Cobalt on Digital MSA
  if(esl_msa_iset_Cobalt(msa2, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_iset(1, msa2->ax, 12, sizeof(char *), abc, assignments) != eslOK){
    return eslFAIL;
  }

//Test 3: msa_iset_Blue on ASCII MSA
  if(esl_msa_iset_Blue(msa, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_iset(0, msa->aseq, 12, sizeof(char *), abc, assignments) != eslOK){
    return eslFAIL;
  }

  //Test4: msa_iset_Blue on Digital MSA
  if(esl_msa_iset_Blue(msa2, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_iset(1, msa2->ax, 12, sizeof(char *), abc, assignments) != eslOK){
    return eslFAIL;
  }

  //Test5: msa_bi_iset_Cobalt on ASCII MSA
  if(esl_msa_bi_iset_Cobalt(msa, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(0, msa->aseq, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }
  larger_set = 0;
  for(int i = 0; i < 12; i++){
    if(assignments[i] == 2){
      larger_set++;
    }
    if(assignments[i] ==1){
      larger_set--;
    }
  }

  if((larger_set < 0) && (larger == 2)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }
   if((larger_set > 0) && (larger == 1)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }

  //Test6: msa_bi_iset_Cobalt on Digital MSA
  if(esl_msa_bi_iset_Cobalt(msa2, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(1, msa2->ax, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }
  larger_set = 0;
  for(int i = 0; i < 12; i++){
    if(assignments[i] == 2){
      larger_set++;
    }
    if(assignments[i] ==1){
      larger_set--;
    }
  }

  if((larger_set < 0) && (larger == 2)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }
   if((larger_set > 0) && (larger == 1)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }

 //Test7: msa_bi_iset_Blue on ASCII MSA
  if(esl_msa_bi_iset_Blue(msa, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(0, msa->aseq, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }
  larger_set = 0;
  for(int i = 0; i < 12; i++){
    if(assignments[i] == 2){
      larger_set++;
    }
    if(assignments[i] ==1){
      larger_set--;
    }
  }

  if((larger_set < 0) && (larger == 2)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }
   if((larger_set > 0) && (larger == 1)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }

  //Test8: msa_bi_iset_Blue on Digital MSA
  if(esl_msa_bi_iset_Blue(msa2, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(1, msa2->ax, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }
  larger_set = 0;
  for(int i = 0; i < 12; i++){
    if(assignments[i] == 2){
      larger_set++;
    }
    if(assignments[i] ==1){
      larger_set--;
    }
  }

  if((larger_set < 0) && (larger == 2)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }
   if((larger_set > 0) && (larger == 1)){
    // check and return value disagree about which set is larger
    return eslFAIL;
  }
 //Test9: msa_bi_iset_Random on ASCII MSA
  if(esl_msa_bi_iset_Random(msa, 0.5, &assignments, r, 0.3) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(0, msa->aseq, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }


  //Test9: msa_bi_iset_Random on Digital MSA
  if(esl_msa_bi_iset_Random(msa2, 0.5, &assignments, r, 0.3) != eslOK){
    return eslFAIL;
  }
  if(check_msa_bi_iset(1, msa2->ax, 12, 0.5, abc, assignments)!= eslOK){
    return eslFAIL;
  }
  esl_msa_Destroy(msa);
  esl_msa_Destroy(msa2);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  free(assignments);
  return eslOK;

  ERROR:
  return eslFAIL;
}
#endif /* eslMSA_ISET_TESTDRIVE*/


/*****************************************************************
 * 6. Example.
 *****************************************************************/ 

#ifdef eslMSA_ISET_EXAMPLE
#include "esl_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_iset.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for iset module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  int *assignments;
  ESL_MSA        *msa     = esl_msa_CreateFromString("\
# STOCKHOLM 1.0\n\
\n\
seq0  AAAAAAAAAA\n\
seq1  AAAAAAAAAA\n\
seq2  AAAAAAAAAC\n\
seq3  AAAAAAAADD\n\
seq4  AAAAAAAEEE\n\
seq5  AAAAAAFFFF\n\
seq6  AAAAAGGGGG\n\
seq7  AAAAHHHHHH\n\
seq8  AAAIIIIIII\n\
seq9  AAKKKKKKKK\n\
seq10 ALLLLLLLLL\n\
seq11 MMMMMMMMMM\n\
//",   eslMSAFILE_STOCKHOLM);
  int status;
  ESL_RANDOMNESS *r = NULL;
  r=esl_randomness_Create(0);
  ESL_ALLOC(assignments, 12 * sizeof(int));  //must be = # of sequences in alignment
  // Make digital copy of the msa
  ESL_MSA *msa2 = esl_msa_Create(12, 10);
  esl_msa_Copy(msa, msa2);
  esl_msa_Digitize(abc, msa2, NULL);
  int larger, larger_set;

  //msa_iset_Cobalt on ASCII MSA
  if(esl_msa_iset_Cobalt(msa, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  
  //msa_iset_Cobalt on Digital MSA
  if(esl_msa_iset_Cobalt(msa2, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }

  //msa_iset_Blue on ASCII MSA
  if(esl_msa_iset_Blue(msa, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
 
  //msa_iset_Blue on Digital MSA
  if(esl_msa_iset_Blue(msa2, 0.5, &assignments, r) != eslOK){
    return eslFAIL;
  }
  
  //msa_bi_iset_Cobalt on ASCII MSA
  if(esl_msa_bi_iset_Cobalt(msa, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
 
  //msa_bi_iset_Cobalt on Digital MSA
  if(esl_msa_bi_iset_Cobalt(msa2, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }

 
 //msa_bi_iset_Blue on ASCII MSA
  if(esl_msa_bi_iset_Blue(msa, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
  
 
  //msa_bi_iset_Blue on Digital MSA
  if(esl_msa_bi_iset_Blue(msa2, 0.5, &assignments, &larger, r) != eslOK){
    return eslFAIL;
  }
 
  
 //msa_bi_iset_Random on ASCII MSA
  if(esl_msa_bi_iset_Random(msa, 0.5, &assignments, r, 0.3) != eslOK){
    return eslFAIL;
  }


  //msa_bi_iset_Random on Digital MSA
  if(esl_msa_bi_iset_Random(msa2, 0.5, &assignments, r, 0.3) != eslOK){
    return eslFAIL;
  }
 
  esl_msa_Destroy(msa);
  esl_msa_Destroy(msa2);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  free(assignments);
  return eslOK;

  ERROR:
  return eslFAIL;
}
#endif /* eslMSA_ISET_TESTDRIVE*/