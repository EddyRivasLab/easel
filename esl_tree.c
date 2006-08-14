/* esl_tree.c
 * Phylogenetic trees.
 * 
 * SVN $Id$
 * SRE, Tue May  2 14:08:42 2006 [St. Louis]
 */


#include <esl_config.h>

#include <math.h>

#include <easel.h>
#include <esl_tree.h>
#include <esl_dmatrix.h>

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
  return T;
  
 FAILURE:
  esl_tree_Destroy(T);
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

  if (T->parent != NULL) free(T->parent);
  if (T->left   != NULL) free(T->left);
  if (T->right  != NULL) free(T->right);
  if (T->ld     != NULL) free(T->ld);
  if (T->rd     != NULL) free(T->rd);
  free(T);
  return;
}


/*****************************************************************
 * 2. Clustering algorithms for tree construction.
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
int
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
  double       x;
  int          k;
  int          status;

  /* Contract checks.
   */
  ESL_ASSERT1((D_original != NULL));               /* matrix exists      */
  ESL_ASSERT1((D_original->n == D_original->m));   /* D is NxN square    */
  ESL_ASSERT1((D_original->n >= 2));               /* >= 2 taxa          */
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
  ESL_ALLOC(idx,    sizeof(int) *  D->n);
  ESL_ALLOC(nin,    sizeof(int) *  D->n);
  ESL_ALLOC(height, sizeof(int) * (D->n-1));
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
	    mx[i][col] = nin[i] * mx[i][col] + nin[j] * mx[j][col] / (double) (nin[i] + nin[j]);
	    break;
	  case eslWPGMA:            mx[i][col] = mx[i][col] + mx[j][col] / 2.;    break;
	  case eslSINGLE_LINKAGE:   mx[i][col] = ESL_MIN(mx[i][col], mx[j][col]); break;
	  case eslCOMPLETE_LINKAGE: mx[i][col] = ESL_MAX(mx[i][col], mx[j][col]); break;
	  default:
	  }
	  mx[col][i] = mx[i][col];
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


/*****************************************************************
 * Timing benchmark.
 *****************************************************************/

/* 
 */
#ifdef eslTREE_BENCHMARK

#include <easel.h>
#include <esl_tree.h>
#include <esl_distance.h>
#include <esl_getopts.h>

static ESL_OPTIONS options[] = {
  /* name    type        default env_var  range toggles req  incompat help                  docgroup */
 { "-h",     eslARG_NONE, FALSE,  NULL,   NULL,  NULL,  NULL,  NULL,  "show help and usage",       0 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[] = "Usage: ./tree-benchmark [-options]";

int
main(int argc, char **argv)
{
  int ntaxa;			
  
  



}
#endif /*eslTREE_BENCHMARK*/


/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
