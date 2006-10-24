/* esl_tree.c
 * Phylogenetic trees.
 * 
 * Contents:
 *   1. The ESL_TREE object.
 *   2. Newick format i/o
 *   3. Tree comparison algorithms.
 *   4. Clustering algorithms for distance-based tree construction.
 *   5. Generating simulated trees.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example code.
 *   9. Copyright notice and license.
 * 
 * SVN $Id$
 * SRE, Tue May  2 14:08:42 2006 [St. Louis]
 */


#include <esl_config.h>

#include <math.h>
#include <string.h>
#include <ctype.h>

#include <easel.h>
#include <esl_tree.h>
#include <esl_dmatrix.h>
#include <esl_stack.h>
#include <esl_vectorops.h>
#include <esl_random.h>

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
  T->taxonlabel    = NULL;
  T->nodelabel     = NULL;

  /* Tree output options default to PHYLIP style
   */
  T->show_unrooted            = FALSE;
  T->show_node_labels         = TRUE;
  T->show_root_branchlength   = FALSE;
  T->show_branchlengths       = TRUE;
  T->show_quoted_labels       = FALSE;
  T->show_numeric_taxonlabels = TRUE;

  T->nalloc = ntaxa;
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
      if (child <= 0) T->parent_of_otu[-child] = parent;
      else	      esl_stack_IPush(ns, child);

      child = T->right[parent];
      if (child <= 0) T->parent_of_otu[-child] = parent;
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
 * 2. Newick format i/o
 *****************************************************************/

/* newick_validate_unquoted():
 *   Returns <eslOK> if we can represent <label> as an unquoted label
 *   in Newick format. (Spaces are ok, but will be converted to
 *   _ on output.)
 */
static int 
newick_validate_unquoted(char *label)
{
  char *sptr;
  for (sptr = label; *sptr != '\0'; sptr++)
    {
      if (! isprint(*sptr))                  return eslFAIL;
      if (strchr("()[]':;,", *sptr) != NULL) return eslFAIL;
    }
  return eslOK;
}

/* newick_validate_quoted():
 *   Returns <eslOK> if we can represent <label> as a 
 *   quoted label in Newick format. (Single quotes will
 *   be converted to '' on output.)
 */
static int
newick_validate_quoted(char *label)
{
  char *sptr;
  for (sptr = label; *sptr != '\0'; sptr++)
    {
      if (! isprint(*sptr))                  return eslFAIL;
    }
  return eslOK;
} 

/* newick_write_unquoted():
 *   Prints <label> to <fp> as an unquoted Newick label.
 */
static int
newick_write_unquoted(FILE *fp, char *label)
{
  char *sptr;

  for (sptr = label; *sptr != '\0'; sptr++)
    {
      if (*sptr == ' ') fputc('_',   fp);
      else              fputc(*sptr, fp);
    }
  return eslOK;
}

/* newick_write_quoted():
 *   Prints <label> to <fp> as a quoted Newick label.
 */
static int
newick_write_quoted(FILE *fp, char *label)
{
  char *sptr;

  fputc('\'', fp);
  for (sptr = label; *sptr != '\0'; sptr++)
    {
      if (*sptr == '\'') fprintf(fp, "''");
      else               fputc(*sptr, fp);        
    }
  fputc('\'', fp);
  return eslOK;
}

/* newick_write_taxonlabel():
 *    Print the label for taxon <v> to stream <fp>.
 *    Tries to print label as an unquoted label, then
 *    as a quoted label, (then fails).
 *    If label isn't available, does nothing.
 *    If label contains invalid characters, throws <eslECORRUPT>.
 */
static int
newick_write_taxonlabel(FILE *fp, ESL_TREE *T, int v)
{
  if (T->taxonlabel == NULL || T->taxonlabel[v] == NULL)
    {
      if (T->show_numeric_taxonlabels)
	fprintf(fp, "%d", v);
      return eslOK;
    }

  if (! T->show_quoted_labels && newick_validate_unquoted(T->taxonlabel[v]) == eslOK)
    newick_write_unquoted(fp, T->taxonlabel[v]);
  else if (newick_validate_quoted(T->taxonlabel[v]) == eslOK)
    newick_write_quoted(fp, T->taxonlabel[v]);
  else
    ESL_ERROR(eslECORRUPT, "bad taxon label");

  return eslOK;
}

/* newick_write_nodelabel():
 *    Print the label for internal node <v> to stream <fp>.
 *    Tries to print label as an unquoted label, then
 *    as a quoted label. 
 *    If label isn't available, does nothing.
 *    If tree's options say not to print node labels, does nothing.
 *    If label contains invalid characters, throws <eslECORRUPT>.
 */
static int
newick_write_nodelabel(FILE *fp, ESL_TREE *T, int v)
{
  if (T->nodelabel    == NULL)      return eslOK;
  if (T->nodelabel[v] == NULL)      return eslOK;
  if (T->show_node_labels != TRUE)  return eslOK;
  
  if (! T->show_quoted_labels && newick_validate_unquoted(T->nodelabel[v]) == eslOK)
    newick_write_unquoted(fp, T->nodelabel[v]);
  else if (newick_validate_quoted(T->nodelabel[v]) == eslOK)
    newick_write_quoted(fp, T->nodelabel[v]);
  else
    ESL_ERROR(eslECORRUPT, "bad node label");

  return eslOK;
}

/* newick_write_branchlength()
 *    Writes the branch length *to* <v>.
 *    If <v> is negative, it's a leaf; if <v> is positive, it's an internal node.
 *    You can't pass the root node 0 to this. 0 always means taxon 0.
 *    There is no branch to the root node.
 */
int
newick_write_branchlength(FILE *fp, ESL_TREE *T, int v)
{
  double branchlength;

  if (! T->show_branchlengths)   return eslOK;
  if (T->parent_of_otu == NULL)  ESL_ERROR(eslECONTRACT, "T must have parent_of_otu");
  
  if (v <= 0)			/* leaf */
    {
      if      (T->left [T->parent_of_otu[-v]] == v) branchlength = T->ld[T->parent_of_otu[-v]];
      else if (T->right[T->parent_of_otu[-v]] == v) branchlength = T->rd[T->parent_of_otu[-v]]; 
      else    ESL_ERROR(eslECORRUPT, "Can't find branch length");
    }
  else				/* internal node */
    {
      if      (T->left [T->parent[v]] == v) branchlength = T->ld[T->parent[v]];
      else if (T->right[T->parent[v]] == v) branchlength = T->rd[T->parent[v]]; 
      else    ESL_ERROR(eslECORRUPT, "Can't find branch length");
    }

  fprintf(fp, ":%f", branchlength);
  return eslOK;
}

/* Function:  esl_tree_WriteNewick()
 * Incept:    SRE, Fri Oct  6 14:35:51 2006 [Janelia]
 *
 * Purpose:   Writes tree <T> to stream <fp> in Newick format.
 *  
 *            Certain options are set in <T> to control output style.
 *            If <T->show_unrooted> is <TRUE>, <T> is printed as an
 *            unrooted tree starting with a trifurcation, a la PHYLIP
 *            format (default=<FALSE>). If <T->show_node_labels> is
 *            <TRUE>, then labels are shown for internal nodes, if any
 *            are available (default=<TRUE>). If
 *            <T->show_branchlengths> is <TRUE>, then branch lengths
 *            are shown, as opposed to just printing a labeled
 *            topology (default=<TRUE>). If
 *            <T->show_root_branchlength> is also <TRUE>, then a 0.0
 *            branchlength is shown to the root node, a la Hein's
 *            TreeAlign Newick format (default=<FALSE>). If
 *            <T->show_quoted_labels> is <TRUE>, then all labels are
 *            shown in Newick's quoted format, as opposed to only
 *            using quoted labels where necessary (default=<FALSE>).
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINCONCEIVABLE> on internal error.
 *
 * Xref:      STL11/74
 */
int
esl_tree_WriteNewick(FILE *fp, ESL_TREE *T)
{
  ESL_STACK *vs = NULL;
  ESL_STACK *cs = NULL;
  int  v;
  char c;
  int  status;

  if ((vs = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto FAILURE; }
  if ((cs = esl_stack_CCreate()) == NULL) { status = eslEMEM; goto FAILURE; }
  
  if ((status = esl_tree_MapTaxaParents(T)) != eslOK) goto FAILURE;
  
  /* Initialization.
   * Push a trifurcation (swallowing the right internal node) if unrooted;
   * else push the first bifurcation.
   * 
   * When we push a trifurcation, the branch lengths will come out fine
   * on output, if the tree followed the correct convention of having
   * a T->rd[0] = 0.0.
   */
  fputc('(', fp);
  if (T->show_unrooted && T->right[0] > 0)
    {
      v = T->right[0];
      if ((status = esl_stack_CPush(cs, 'x'))         != eslOK) goto FAILURE;
      if ((status = esl_stack_IPush(vs, T->right[v])) != eslOK) goto FAILURE;
      if ((status = esl_stack_CPush(cs, ','))         != eslOK) goto FAILURE;
      if ((status = esl_stack_CPush(cs, 'x'))         != eslOK) goto FAILURE;
      if ((status = esl_stack_IPush(vs, T->left[v]))  != eslOK) goto FAILURE;
    }
  else 
    {
      if ((status = esl_stack_CPush(cs, 'x'))         != eslOK) goto FAILURE;
      if ((status = esl_stack_IPush(vs, T->right[0])) != eslOK) goto FAILURE;
    }
  if ((status = esl_stack_CPush(cs, ','))             != eslOK) goto FAILURE;
  if ((status = esl_stack_CPush(cs, 'x'))             != eslOK) goto FAILURE;
  if ((status = esl_stack_IPush(vs, T->left[0]))      != eslOK) goto FAILURE;


  /* Main iteration. Pop off stacks 'til they're empty.
   */
  while ((status = esl_stack_CPop(cs, &c)) == eslOK)
    {
      if (c == ',') { fputc(',', fp); continue; } /* comma doesn't have a v stacked with it */

      if ((status = esl_stack_IPop(vs, &v)) != eslOK) goto FAILURE;

      switch (c) {
      case 'x':			/* a subtree, which could be a node or a taxon: */
	if (v > 0)		/* internal node 1..N-2*/
	  {
	    fputc('(', fp);
	    if ((status = esl_stack_CPush(cs, ')'))         != eslOK) goto FAILURE;
	    if ((status = esl_stack_IPush(vs, v))           != eslOK) goto FAILURE;
	    if ((status = esl_stack_CPush(cs, 'x'))         != eslOK) goto FAILURE;
	    if ((status = esl_stack_IPush(vs, T->right[v])) != eslOK) goto FAILURE;
	    if ((status = esl_stack_CPush(cs, ','))         != eslOK) goto FAILURE;
	    if ((status = esl_stack_CPush(cs, 'x'))         != eslOK) goto FAILURE;
	    if ((status = esl_stack_IPush(vs, T->left[v]))  != eslOK) goto FAILURE;
	  }
	else			/* taxon -(N-1)..0 */
	  { 	    /* -v below to convert taxon code to 0..N-1 */
	    if ((status = newick_write_taxonlabel  (fp, T, -v)) != eslOK) goto FAILURE;
	    if ((status = newick_write_branchlength(fp, T,  v)) != eslOK) goto FAILURE;
	  }
	break;

      case ')':			/* closing an internal node. v > 0 is a node code. */
	fputc(')', fp);
	if ((status = newick_write_nodelabel   (fp, T, v)) != eslOK) goto FAILURE;
	if ((status = newick_write_branchlength(fp, T, v)) != eslOK) goto FAILURE;
	break;

      default:
	ESL_ERROR(eslEINCONCEIVABLE, "bad state code");
	break;
      }
    }
  
  /* Termination
   */
  fputc(')', fp);
  newick_write_nodelabel(fp, T, 0);
  if (T->show_branchlengths && T->show_root_branchlength) fprintf(fp, ":0.0");
  fputc(';', fp);
  fputc('\n', fp);

  esl_stack_Destroy(vs);
  esl_stack_Destroy(cs);
  return eslOK;

 FAILURE:
  if (vs != NULL) esl_stack_Destroy(vs);
  if (cs != NULL) esl_stack_Destroy(cs);
  return status;

}


/*-------------------- end, Newick i/o --------------------------*/


/*****************************************************************
 * 3. Tree comparison algorithms
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
  for (g = T1->N-2; g >= 0; g--)
    {
      child = T1->left[g];
      if (child <= 0)  a = T2->parent_of_otu[-child]; 
      else             a = T2->parent[Mg[child]];

      child = T1->right[g];
      if (child <= 0)  b = T2->parent_of_otu[-child]; 
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
 * 4. Clustering algorithms for tree construction.
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
      if (idx[j] > 0) T->rd[N-2] -= height[idx[j]];      
      
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
 * 5. Generating simulated trees
 *****************************************************************/

/* Function:  esl_tree_Simulate()
 * Synopsis:  Generate a random rooted ultrametric tree.
 * Incept:    SRE, Mon Oct  2 11:36:22 2006 [Janelia]
 *
 * Purpose:   Generate a random rooted ultrametric tree of <N> taxa,
 *            using the algorithm of Kuhner and Felsenstein (1996).
 *            
 *            The branch lengths are generated by choosing <N-1>
 *            exponentially distributed split times, with decreasing
 *            expectations of $\frac{1}{2},\frac{1}{3}..\frac{1}{N}$
 *            as the simulation proceeds from the root. Thus the
 *            total expected branch length on the tree is
 *            $\sum_{k=2}^{N} \frac{1}{k}$.
 *
 * Args:      r     - random number source
 *            N     - number of taxa (>= 2)
 *            ret_T - RETURN: sampled tree
 *
 * Returns:   <eslOK> on success, and the new tree is allocated
 *            here and returned via <ret_tree>; caller is 
 *            responsible for free'ing it.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case
 *            the <ret_T> is returned <NULL>.
 *
 * Xref:      STL11/65.
 */
int
esl_tree_Simulate(ESL_RANDOMNESS *r, int N, ESL_TREE **ret_T)
{
  ESL_TREE       *T          = NULL;
  int            *branchpapa = NULL;
  int            *branchside = NULL;
  int       nactive;
  double    d;
  int       node;
  int       bidx;	        	/* index of an active branch */
  int       status;

  ESL_DASSERT1(r != NULL);
  ESL_DASSERT1(N >= 2);

  /* Kuhner/Felsenstein uses a list of active branches,
   * which we implement by tracking the index of the parent
   * node (in <branchpapa>) and a 0/1 flag (in <branchside>)
   * for the branch to the left vs. right child.
   */
  if ((T = esl_tree_Create(N)) == NULL)  goto FAILURE;
  ESL_ALLOC(branchpapa, sizeof(int) * N);
  ESL_ALLOC(branchside, sizeof(int) * N);
  
  /* Initialize: add two branches from the root
   * onto the active list, and set internal node
   * counter to start at 1.
   */
  branchpapa[0] = 0;   branchside[0] = 0;
  branchpapa[1] = 0;   branchside[1] = 1;
  nactive = 2;
  node    = 1;			

  /* Algorithm proceeds by iterating:
   *    1. choose random time <d> from exponential(1/nactive)
   *    2. choose random active branch, <bidx>
   *    3. add new <node> to active branch at length d
   *    4. add d to all other active branches      
   *    5. delete the old parent branch from the active list,
   *       add the two new child branches to the active list
   */
  while (nactive < N)
    {
      d               = (double) nactive * -log(esl_rnd_UniformPositive(r));
      bidx            = esl_rnd_Choose(r, nactive);
      T->parent[node] = branchpapa[bidx];
      
      if (branchside[bidx] == 0) {
	T->left[branchpapa[bidx]]   = node;
	T->ld  [branchpapa[bidx]]  += d;
      } else {
	T->right[branchpapa[bidx]]  = node;
	T->rd   [branchpapa[bidx]] += d;
      }

      ESL_SWAP(branchpapa[bidx], branchpapa[nactive-1], int);
      ESL_SWAP(branchside[bidx], branchside[nactive-1], int);
      for (bidx = 0; bidx < nactive-1; bidx++) {
	if (branchside[bidx] == 0) T->ld[branchpapa[bidx]] += d;
	else                       T->rd[branchpapa[bidx]] += d;
      }
      
      /* delete the branch at nactive-1 that we just added to;
       * replace it with two new branches
       */
      branchpapa[nactive-1]  = node;  branchside[nactive-1] = 0;
      branchpapa[nactive]    = node;  branchside[nactive]   = 1;
      node++;
      nactive++;
    }

  /* Terminate by adding the N taxa to the N active branches.
   */
  d = (double) N * -log(esl_rnd_UniformPositive(r));
  for (bidx = 0; bidx < N; bidx++)
    {
      if (branchside[bidx] == 0) {
	T->left[branchpapa[bidx]]  =  -bidx; /* taxa indices stored as neg #'s */
	T->ld  [branchpapa[bidx]]  += d;
      } else {
	T->right[branchpapa[bidx]] =  -bidx;
	T->rd  [branchpapa[bidx]]  += d;
      }
    }

  *ret_T = T; 
  free(branchpapa);
  free(branchside);
  return eslOK;

 FAILURE:
  if (T          != NULL) esl_tree_Destroy(T);
  if (branchpapa != NULL) free(branchpapa);
  if (branchside != NULL) free(branchside);
  *ret_T = NULL;
  return status;
}


/* Function:  esl_tree_ToDistanceMatrix()
 * Synopsis:  Obtain a pairwise distance matrix from a tree.
 * Incept:    SRE, Fri Oct  6 13:50:37 2006 [Janelia]
 *
 * Purpose:   Given tree <T>, calculate a pairwise distance matrix
 *            and return it in <ret_D>.
 *            
 * Note:      Algorithm here is O(N^3). It can probably be improved.
 *            There ought to be a more efficient recursion that
 *            saves recalculating node-node distances inside the tree.
 *            All we do here is a brute force, upwards O(N) LCA 
 *            search for each of the N^2 taxon pairs. 
 *
 * Args:      T     - input tree 
 *            ret_D - RETURN: the new distance matrix    
 *
 * Returns:   <eslOK> on success, and <ret_D> points to the distance 
 *            matrix, which caller is responsible for free'ing with
 *            <esl_dmatrix_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case
 *            <ret_D> is returned <NULL>.
 *
 * Xref:      STL11/66.
 */
int
esl_tree_ToDistanceMatrix(ESL_TREE *T, ESL_DMATRIX **ret_D)
{
  ESL_DMATRIX *D = NULL;
  int i,j;			/* a pair of taxa {0..N-1}           */
  int a,b;			/* a pair of internal nodes {0..N-2} */
  int p;			/* a tmp parent index */
  double d;			/* ij distance */
  int status;

  D = esl_dmatrix_Create(T->N, T->N); /* creates a NxN square symmetric matrix; really only need triangular */
  if (D == NULL) { status = eslEMEM; goto FAILURE; }

  if ((status = esl_tree_MapTaxaParents(T)) != eslOK) goto FAILURE;

  for (i = 0; i < T->N; i++)
    {
      D->mx[i][i] = 0.;		/* by definition */
      for (j = i; j < T->N; j++)
	{
	  a  = T->parent_of_otu[i];
	  b  = T->parent_of_otu[j];
	  d  = (T->left[a] == -i) ? T->ld[a] : T->rd[a];
	  d += (T->left[b] == -j) ? T->ld[b] : T->rd[b];
	  while (a != b)	/* a brute force LCA algorithm */
	    {
	      if (a < b) ESL_SWAP(a, b, int);
	      p  = T->parent[a];
	      d += (T->left[p] == a) ? T->ld[p] : T->rd[p];
	      a  = p;
	    }

	  D->mx[i][j] = D->mx[j][i] = d;
	}
    }

  *ret_D = D;
  return eslOK;

 FAILURE:
  if (D != NULL) esl_dmatrix_Destroy(D);
  *ret_D = NULL;
  return status;
}



/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef eslTREE_TESTDRIVE

static int
verify_ultrametricity(ESL_TREE *T)
{
  double *d = NULL;		/* Distance from root for each OTU */
  int status;
  int i, child, parent;
  
  /* First, calculate distance from root to each taxon.
   * (This chunk of code might be useful to put on its own someday.)
   */
  ESL_ALLOC(d, sizeof(double) * T->N);
  esl_tree_MapTaxaParents(T);
  for (i = 0; i < T->N; i++)
    {
      d[i]   = 0.0;
      child  = i;
      parent = T->parent_of_otu[i];
      if       (T->left[parent]  == -i) d[i] += T->ld[parent];
      else if  (T->right[parent] == -i) d[i] += T->rd[parent];
      else     ESL_ERROR(eslEINCONCEIVABLE, "oops");

      while (parent != 0)	/* upwards to the root */
	{
	  child  = parent;
	  parent = T->parent[child];
	  if      (T->left[parent]  == child) d[i] += T->ld[parent];
	  else if (T->right[parent] == child) d[i] += T->rd[parent];
	  else    ESL_ERROR(eslEINCONCEIVABLE, "oops");
	}
    }

  /* In an ultrametric tree, all those distances must be equal.
   */
  status = eslOK;
  if (d[0] == 0.0)
    {
      for (i = 1; i < T->N; i++)
	if (d[i] != 0.0) 
	  { status = eslFAIL; break; }
    }
  else
    {
      for (i = 1; i < T->N; i++)
	if ((fabs(d[i] - d[0]) / d[0]) > 0.0001) 
	  { status = eslFAIL; break; }
    }

  free(d);
  return status;
  
 FAILURE:
  if (d != NULL) free(d);
  return status;
}

static int
utest_WriteNewick(ESL_RANDOMNESS *r, int ntaxa)
{
  ESL_TREE *T1 = NULL;

  if (esl_tree_Simulate(r, ntaxa, &T1) != eslOK)  abort();
  if (esl_tree_WriteNewick(stdout, T1) != eslOK)  abort();
  
  esl_tree_Destroy(T1);
  return eslOK;
}


static int
utest_UPGMA(ESL_RANDOMNESS *r, int ntaxa)
{
  ESL_TREE    *T1 = NULL;
  ESL_TREE    *T2 = NULL;
  ESL_DMATRIX *D  = NULL;

  if (esl_tree_Simulate(r, ntaxa, &T1)  != eslOK) abort();
  if (esl_tree_ToDistanceMatrix(T1, &D) != eslOK) abort();
  if (esl_tree_UPGMA(D, &T2)            != eslOK) abort();

  esl_dmatrix_Dump(stdout, D, NULL, NULL);
  esl_tree_WriteNewick(stdout, T1);
  esl_tree_WriteNewick(stdout, T2);

  if (verify_ultrametricity(T1)         != eslOK) abort();
  if (verify_ultrametricity(T2)         != eslOK) abort();
  if (esl_tree_Compare(T1, T2)          != 0)     abort();

  esl_tree_Destroy(T1);
  esl_tree_Destroy(T2);
  esl_dmatrix_Destroy(D);
  return eslOK;
}

#endif /*eslTREE_TESTDRIVE*/
/*-------------------- end, unit tests  -------------------------*/


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef eslTREE_TESTDRIVE

/* 
 * gcc -g -Wall -o test -L. -I. -DeslTREE_TESTDRIVE esl_tree.c -leasel -lm
 * gcc -g -Wall -o test -L. -I. -DeslTEST_THROWING -DeslTREE_TESTDRIVE esl_msa.c -leasel -lm
 * ./test
 */
#include <easel.h>
#include <esl_tree.h>
#include <esl_random.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r = NULL;
  int ntaxa;

  r     = esl_randomness_Create(42);
  ntaxa = 20;
  
  if ((utest_WriteNewick(r, ntaxa)) != eslOK)  abort();
  if ((utest_UPGMA(r, ntaxa))       != eslOK)  abort();

  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*eslTREE_TESTDRIVE*/
/*-------------------- end, test driver  -------------------------*/




/*****************************************************************
 * 8. Example.
 *****************************************************************/
#ifdef eslTREE_EXAMPLE
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
#endif /*eslTREE_EXAMPLE*/



/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
