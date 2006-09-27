/* esl_tree.c
 * Phylogenetic trees.
 * 
 * Contents:
 *   1. The ESL_TREE object.
 *   2. Tree comparison algorithms.
 *   3. Clustering algorithms for distance-based tree construction.
 *   4. Unit tests.
 *   5. Test driver.
 *   6. Example code.
 *   7. Copyright notice and license.
 * 
 * SVN $Id$
 * SRE, Tue May  2 14:08:42 2006 [St. Louis]
 */


#include <esl_config.h>

#include <math.h>

#include <easel.h>
#include <esl_tree.h>
#include <esl_dmatrix.h>
#include <esl_stack.h>
#include <esl_vectorops.h>

/*****************************************************************
 * 1. The ESL_TREE object.
 *****************************************************************/

/* Function:  esl_tree_Create()
 * Incept:    SRE, Tue May  2 14:10:17 2006 [St. Louis]
 *
 * Purpose:   Allocates an empty tree structure for <ntaxa> taxa,
 *            and return a ptr to it. <ntaxa> must be $\geq 2$.
 *
 * Args:      <ntaxa>   - number of taxa
 *
 * Returns:   pointer to the new <ESL_TREE> object; caller frees 
 *            this with <esl_tree_Destroy()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
ESL_TREE *
esl_tree_Create(int ntaxa)
{
  ESL_TREE *T = NULL;
  int       i;
  int       status;

  /* Contract verification  */
  ESL_DASSERT1((ntaxa >= 2));

  /* 1st allocation round  */
  ESL_ALLOC(T, sizeof(ESL_TREE));
  T->parent = NULL;
  T->left   = NULL;
  T->right  = NULL;
  T->ld     = NULL;
  T->rd     = NULL;
  
  /* 2nd allocation round */
  T->N    = ntaxa;
  ESL_ALLOC(T->parent, sizeof(int)    * (ntaxa-1));
  ESL_ALLOC(T->left,   sizeof(int)    * (ntaxa-1));
  ESL_ALLOC(T->right,  sizeof(int)    * (ntaxa-1));
  ESL_ALLOC(T->ld,     sizeof(double) * (ntaxa-1));
  ESL_ALLOC(T->rd,     sizeof(double) * (ntaxa-1));
  
  for (i = 0; i < ntaxa-1; i++)
    {
      T->parent[i] = 0;
      T->left[i  ] = 0;
      T->right[i]  = 0;
      T->ld[i]   = 0.;
      T->rd[i]   = 0.;
    }

  /* Optional info starts NULL
   */
  T->parent_of_otu = NULL;

  return T;
  
 FAILURE:
  esl_tree_Destroy(T);
  return NULL;
}

/* Function:  esl_tree_MapTaxaParents()
 * Synopsis:  Construct the lookup map for each taxon's parent node.
 * Incept:    SRE, Fri Sep 22 13:39:49 2006 [Janelia]
 *
 * Purpose:   Constructs the <T->parent_of_otu[]> map in the tree
 *            structure <T>, by an O(N) traversal of the tree.
 *            Upon return, <T->parent_of_otu[i]> is the index
 *            of the internal node that taxon <i> is a child of.
 *
 * Args:      T   - the tree structure to map
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on internal allocation error. In this case, the tree is 
 *            returned unchanged.
 *
 * Xref:      STL11/63
 */
int
esl_tree_MapTaxaParents(ESL_TREE *T)
{
  ESL_STACK *ns = NULL;
  int parent, child;
  int status;

  if (T->parent_of_otu != NULL) return eslOK; /* map already exists. */

  ESL_ALLOC(T->parent_of_otu, sizeof(int) * T->N);
#if (eslDEBUGLEVEL >= 1)
  esl_vec_ISet(T->parent_of_otu, T->N, -1);
#endif

  if ((ns = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto FAILURE; }
  if ((status = esl_stack_IPush(ns, 0)) != eslOK) goto FAILURE;	/* init: push root  */

  while ((status = esl_stack_IPop(ns, &parent)) == eslOK)
    {
      child = T->left[parent];
      if (child <= 0) T->parent_of_otu[child] = parent;
      else	      esl_stack_IPush(ns, child);

      child = T->right[parent];
      if (child <= 0) T->parent_of_otu[child] = parent;
      else	      esl_stack_IPush(ns, child);
    }
  esl_stack_Destroy(ns);

#if (eslDEBUGLEVEL >= 1)
  for (child = 0; child < N; child++) assert(T->parent_of_otu[child] >= 0);
#endif
  return eslOK;

 FAILURE:
  if (ns               != NULL) esl_stack_Destroy(ns);
  if (T->parent_of_otu != NULL) { free(T->parent_of_otu); T->parent_of_otu = NULL; }
  return status;
}
  


/* Function:  esl_tree_Destroy()
 * Incept:    SRE, Tue May  2 14:18:31 2006 [St. Louis]
 *
 * Purpose:   Frees an <ESL_TREE> object.
 */
void
esl_tree_Destroy(ESL_TREE *T)
{
  if (T == NULL) return;

  if (T->parent        != NULL) free(T->parent);
  if (T->left          != NULL) free(T->left);
  if (T->right         != NULL) free(T->right);
  if (T->ld            != NULL) free(T->ld);
  if (T->rd            != NULL) free(T->rd);
  if (T->parent_of_otu != NULL) free(T->parent_of_otu);
  free(T);
  return;
}
/*----------------- end, ESL_TREE object -----------------------*/


/*****************************************************************
 * 2. Tree comparison algorithms
 *****************************************************************/

/* Function:  esl_tree_Compare()
 * Incept:    SRE, Fri Sep 22 14:05:09 2006 [Janelia]
 *
 * Purpose:   Given two trees <T1> and <T2> for the same
 *            set of <N> taxa (represented in the trees by the same 
 *            indices, <0..N-1>), compare the topologies of the
 *            two trees.
 *            
 *            For comparing unrooted topologies, be sure that <T1> and
 *            <T2> both obey the unrooted tree convention that the
 *            "root" is placed on the branch to taxon 0. (That is,
 *            <T->left[0] = 0>.)
 *            
 * Returns:   <eslOK> if tree topologies are identical. <eslFAIL>
 *            if they aren't.           
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_tree_Compare(ESL_TREE *T1, ESL_TREE *T2)
{
  int *Mg;			/* the M(g) tree-mapping function    */
  int  g, child;		/* node indices for parent, children */
  int  a,b;
  int  status;

  /* We need taxon parent map in tree 2, but not tree 1.
   */
  if ((status = esl_tree_MapTaxaParents(T2)) != eslOK) goto FAILURE;

  /* We're going to use the tree mapping function M(g) [Goodman79]:
   * M[g] for node g in T1 is the index of the lowest node in T2
   * that contains the same children taxa as the subtree 
   * under g in T1.
   */
  ESL_ALLOC(Mg, sizeof(int) * (T1->N-1));

  /* We use the SDI algorithm [ZmasekEddy01] to construct M(g),
   * by postorder traversal of T1
   */
  for (g = T1->N-1; g >= 0; g--)
    {
      child = T1->left[g];
      if (child <= 0)  a = T2->parent_of_otu[child]; 
      else             a = T2->parent[Mg[child]];

      child = T1->right[g];
      if (child <= 0)  b = T2->parent_of_otu[child]; 
      else             b = T2->parent[Mg[child]];

      if (a != b) { free(Mg); return eslFAIL; } /* a shortcut in SDI: special case for exact tree comparison */
      Mg[g] = a;
    }

  free(Mg);
  return eslOK;

 FAILURE:
  if (Mg != NULL) free(Mg);
  return status;
}

/*----------------- end, tree comparison  -----------------------*/






/*****************************************************************
 * 3. Clustering algorithms for tree construction.
 *****************************************************************/

/* UPGMA, average-link, minimum-link, and maximum-link clustering
 * are all implemented by one algorithm, cluster_engine(). We define some flags
 * (within the scope of the tree module) to control the behavior,
 * as we call the algorithm engine from four different API functions.
 */
#define eslUPGMA            0
#define eslWPGMA            1
#define eslSINGLE_LINKAGE   2
#define eslCOMPLETE_LINKAGE 3

/* cluster_engine()
 * 
 * Implements four clustering algorithms for tree construction:
 * UPGMA, WPGMA, single-linkage, and maximum-linkage. These differ
 * only by the rule used to construct new distances after joining
 * two clusters i,j.
 * 
 * Input <D_original> is a symmetric distance matrix, for <D->n> taxa.
 * The diagonal is all 0's, and off-diagonals are $\geq 0$. <D->n>
 * must be at least two.
 * 
 * <mode> is one of <eslUPGMA>, <eslWPGMA>, <eslSINGLE_LINKAGE>, or
 * <eslCOMPLETE_LINKAGE>: a flag specifying which algorithm to use.
 * 
 * The output is a tree structure, returned in <ret_T>.
 * 
 * Returns <eslOK> on success.
 * 
 * Throws <eslEMEM> on allocation failure.
 */
static int
cluster_engine(ESL_DMATRIX *D_original, int mode, ESL_TREE **ret_T)
{
  ESL_DMATRIX *D = NULL;
  ESL_TREE    *T = NULL;
  double      *height = NULL;	/* height of internal nodes  [0..N-2]          */
  int         *idx    = NULL;	/* taxa or node index of row/col in D [0..N-1] */
  int         *nin    = NULL;	/* # of taxa in clade in row/col in D [0..N-1] */
  int          N;
  int          i = 0, j = 0;
  int          row,col;
  double       minD;
  int          status;

  /* Contract checks.
   */
  ESL_DASSERT1((D_original != NULL));               /* matrix exists      */
  ESL_DASSERT1((D_original->n == D_original->m));   /* D is NxN square    */
  ESL_DASSERT1((D_original->n >= 2));               /* >= 2 taxa          */
#if (eslDEBUGLEVEL >=1)
  for (i = 0; i < D_original->n; i++) {
    assert(D_original->mx[i][i] == 0.);	           /* self-self d = 0    */
    for (j = i+1; j < D_original->n; j++)	   /* D symmetric        */
      assert(D_original->mx[i][j] == D_original->mx[j][i]);
  }
#endif

  /* Allocations.
   * NxN copy of the distance matrix, which we'll iteratively whittle down to 2x2;
   * tree for N taxa;
   */
  if ((D = esl_dmatrix_Duplicate(D_original)) == NULL) return eslEMEM;
  if ((T = esl_tree_Create(D->n))             == NULL) return eslEMEM;
  ESL_ALLOC(idx,    sizeof(int)    *  D->n);
  ESL_ALLOC(nin,    sizeof(int)    *  D->n);
  ESL_ALLOC(height, sizeof(double) * (D->n-1));
  for (i = 0; i < D->n;   i++) idx[i]    = -i; /* assign taxa indices to row/col coords */
  for (i = 0; i < D->n;   i++) nin[i ]   = 1;  /* each cluster starts as 1  */
  for (i = 0; i < D->n-1; i++) height[i] = 0.; 


  for (N = D->n; N >= 2; N--)
    {
      /* Find minimum in our current N x N matrix
       */
      minD = HUGE_VAL;
      for (row = 0; row < N; row++)
	for (col = row+1; col < N; col++)
	  if (D->mx[row][col] < minD)
	    {
	      minD = D->mx[row][col];
	      i    = row;
	      j    = col;
	    }

      /* We're joining node at row/col i with node at row/col j.
       * Add node (index = N-2) to the tree at height minD/2.
       */
      T->left[N-2]  = idx[i];
      T->right[N-2] = idx[j];
      height[N-2]   = minD / 2.;

      /* Set the branch lengths
       */
      T->ld[N-2] = T->rd[N-2] = height[N-2];
      if (idx[i] > 0) T->ld[N-2] -= height[idx[i]];
      if (idx[j] > 0) T->ld[N-2] -= height[idx[j]];      
      
      /* If either node was an internal node, record parent in it.
       */
      if (idx[i] > 0)  T->parent[idx[i]] = N-2;
      if (idx[j] > 0)  T->parent[idx[j]] = N-2;
      
      /* Now, build a new matrix by merging row i+j and col i+j.
       *  1. move j to N-1 (unless it's already there)
       *  2. move i to N-2 (unless it's already there)
       */
      if (j != N-1)
	{
	  for (row = 0; row < N; row++)
	    ESL_SWAP(D->mx[row][N-1], D->mx[row][j], double);
	  for (col = 0; col < N; col++)
	    ESL_SWAP(D->mx[N-1][col], D->mx[j][col], double);
	  ESL_SWAP(idx[j],  idx[N-1],  int);
	  ESL_SWAP(nin[j], nin[N-1], int);
	}
      if (i != N-2)
	{
	  for (row = 0; row < N; row++)
	    ESL_SWAP(D->mx[row][N-2], D->mx[row][i], double);
	  for (col = 0; col < N; col++)
	    ESL_SWAP(D->mx[N-2][col], D->mx[i][col], double);
	  ESL_SWAP(idx[i], idx[N-2], int);
	  ESL_SWAP(nin[i], nin[N-2], int);
	}
      i = N-2;
      j = N-1;

      /* 3. merge i (now at N-2) with j (now at N-1) 
       *    according to the desired clustering rule.
       */
      for (col = 0; col < N; col++)
	{
	  switch (mode) {
	  case eslUPGMA: 
	    D->mx[i][col] = (nin[i] * D->mx[i][col] + nin[j] * D->mx[j][col]) / (double) (nin[i] + nin[j]);
	    break;
	  case eslWPGMA:            D->mx[i][col] = (D->mx[i][col] + D->mx[j][col]) / 2.;    break;
	  case eslSINGLE_LINKAGE:   D->mx[i][col] = ESL_MIN(D->mx[i][col], D->mx[j][col]);   break;
	  case eslCOMPLETE_LINKAGE: D->mx[i][col] = ESL_MAX(D->mx[i][col], D->mx[j][col]);   break;
	  default:                  ESL_FAIL(eslEINCONCEIVABLE, "no such strategy");
	  }
	  D->mx[col][i] = D->mx[i][col];
	}

      /* row/col i is now the new cluster, and it corresponds to node N-2
       * in the tree (remember, N is decrementing at each iteration).
       * row/col j (N-1) falls away when we go back to the start of the loop 
       * and decrement N. 
       */
      nin[i] += nin[j];
      idx[i]  = N-2;
    }  

  esl_dmatrix_Destroy(D);
  free(height);
  free(idx);
  free(nin);
  if (ret_T != NULL) *ret_T = T;
  return eslOK;

 FAILURE:
  if (D      != NULL) esl_dmatrix_Destroy(D);
  if (T      != NULL) esl_tree_Destroy(T);
  if (height != NULL) free(height);
  if (idx    != NULL) free(idx);
  if (nin    != NULL) free(nin);
  if (ret_T != NULL) *ret_T = NULL;
  return status;
}


/* Function:  esl_tree_UPGMA()
 * Incept:    SRE, Wed May  3 15:14:17 2006 [St. Louis]
 *
 * Purpose:   Given distance matrix <D>, use the UPGMA algorithm
 *            to construct a tree <T>.
 *
 * Returns:   <eslOK> on success; the tree is returned in <ret_T>,
 *            and must be freed by the caller with <esl_tree_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation problem, and <ret_T> is set <NULL>.
 */
int
esl_tree_UPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T)
{
  return cluster_engine(D, eslUPGMA, ret_T);
}

/* Function:  esl_tree_WPGMA()
 * Incept:    SRE, Wed May  3 15:47:13 2006 [St. Louis]
 *
 * Purpose:   Given distance matrix <D>, use the WPGMA algorithm
 *            to construct a tree <T>.
 *
 * Returns:   <eslOK> on success; the tree is returned in <ret_T>,
 *            and must be freed by the caller with <esl_tree_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation problem, and <ret_T> is set <NULL>.
 */
int
esl_tree_WPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T)
{
  return cluster_engine(D, eslWPGMA, ret_T);
}

/* Function:  esl_tree_SingleLinkage()
 * Incept:    SRE, Wed May  3 15:49:06 2006 [St. Louis]
 *
 * Purpose:   Given distance matrix <D>, construct a single-linkage
 *            (minimum distances) clustering tree <T>.
 *
 * Returns:   <eslOK> on success; the tree is returned in <ret_T>,
 *            and must be freed by the caller with <esl_tree_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation problem, and <ret_T> is set <NULL>.
 */
int
esl_tree_SingleLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T)
{
  return cluster_engine(D, eslSINGLE_LINKAGE, ret_T);
}

/* Function:  esl_tree_CompleteLinkage()
 * Incept:    SRE, Wed May  3 15:49:14 2006 [St. Louis]
 *
 * Purpose:   Given distance matrix <D>, construct a complete-linkage
 *            (maximum distances) clustering tree <T>.
 *
 * Returns:   <eslOK> on success; the tree is returned in <ret_T>,
 *            and must be freed by the caller with <esl_tree_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation problem, and <ret_T> is set <NULL>.
 */
int
esl_tree_CompleteLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T)
{
  return cluster_engine(D, eslCOMPLETE_LINKAGE, ret_T);
}
/*----------------- end, clustering algorithms  ----------------*/



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef eslTESTDRIVE_TREE


#endif /*eslTESTDRIVE_TREE*/
/*-------------------- end, unit tests  -------------------------*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef eslTESTDRIVE_TREE




#endif /*eslTESTDRIVE_TREE*/
/*-------------------- end, test driver  -------------------------*/




/*****************************************************************
 * 6. Example.
 *****************************************************************/
/*::cexcerpt::tree_example::begin::*/
/* To compile: gcc -g -Wall -o example -I. -DeslTREE_EXAMPLE esl_tree.c esl_dmatrix.c esl_msa.c easel.c -lm
 *         or: gcc -g -Wall -o example -I. -L. -DeslTREE_EXAMPLE esl_tree.c -leasel -lm
 *     To run: ./example <MSA file>
 */
#include <easel.h>
#include <esl_msa.h>
#include <esl_distance.h>
#include <esl_tree.h>

int main(int argc, char **argv)
{
  ESL_TREE    *tree;
  ESL_MSAFILE *afp;
  ESL_MSA     *msa;
  ESL_DMATRIX *D;

  esl_msafile_Open(argv[1], eslMSAFILE_UNKNOWN, NULL, &afp);
  esl_msa_Read(afp, &msa);
  esl_msafile_Close(afp);

  esl_dst_CDiffMx(msa->aseq, msa->nseq, &D);
  esl_tree_UPGMA(D, &tree);

  esl_tree_Destroy(tree);
  esl_msa_Destroy(msa);
  esl_dmatrix_Destroy(D);
  return eslOK;
}
/*::cexcerpt::tree_example::end::*/




/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
