/* esl_tree.c
 * Phylogenetic trees.
 * 
 * SVN $Id$
 * SRE, Tue May  2 13:54:30 2006 [St. Louis]
 */
#ifndef ESL_TREE_INCLUDED
#define ESL_TREE_INCLUDED

#include <esl_dmatrix.h>
#include <esl_random.h>

/* Object: ESL_TREE
 *
 * All trees are represented as rooted trees, starting from
 * node 0. For N taxa, there are N-1 internal nodes, numbered
 * 0..N-2. Taxa on leaves are numbered 0..N-1, and represented
 * in <parent>, <left>, <right> as negative numbers.
 * 
 */
typedef struct {
  int   N;		/* number of taxa */

  /* (Mandatory) information in the internal nodes of a rooted tree.
   * There are N-1 nodes, numbered 0..N-2, with the root at 0,
   * so each array below is indexed [0..N-2].
   * There is no ambiguity between taxon 0/root node 0, because 
   * a taxon can't be a parent, and the root node can't be a child.
   * For an unrooted tree, by convention, taxon 0 is the outgroup: T->left[0] = 0,
   * and T->rd[0] = 0.0.
   */
  int    *parent;	/* index of parent of node: values are 0..N-2; parent of root 0 = 0 */
  int    *left;		/* index of left child:  values are -(N-1)..0=taxa; 1..N-2=nodes */
  int    *right;	/* index of right child: values are -(N-1)..0=taxa; 1..N-2=nodes */
  double *ld;	        /* left branch length under node: values are >= 0 */
  double *rd;	        /* right branch length under node: values are >= 0 */

  /* Derived (optional) information, that we can reconstruct if
   * we need to from the mandatory info above.
   */
  int    *parent_of_otu;  /* for each taxon [0..N-1]: index of its parent node, 0..N-2. */

  /* Optional information
   */
  char  **taxonlabel;	  /* labels for taxa: [0..N-1] array of char strings */
  char  **nodelabel;	  /* labels for nodes: [0..N-2] array of char strings */

  /* Tree output options.
   */
  int   show_unrooted;	        /* TRUE to output 'root' as a trifurcation (a la PHYLIP) */
  int   show_node_labels;       /* TRUE to output labels for interior nodes */
  int   show_root_branchlength; /* TRUE to show 0.0 branch length to root node (a la TreeAlign) */
  int   show_branchlengths;	/* TRUE to output branch lengths */
  int   show_quoted_labels;	/* TRUE to output ALL labels as quoted labels */
  int   show_numeric_taxonlabels;/* TRUE to output taxa labels as their 0..N-1 indices if no other taxonlabel is present */

  /* Memory allocation information, when growing a tree (on input, for example)
   */
  int     nalloc;	/* current allocated # of taxa */

} ESL_TREE;



/* 1. The ESL_TREE object.
 */
extern ESL_TREE *esl_tree_Create(int ntaxa);
extern int       esl_tree_Grow(ESL_TREE *T);
extern int       esl_tree_MapTaxaParents(ESL_TREE *T);
extern void      esl_tree_Destroy(ESL_TREE *T);

/* 2. Newick format i/o
 */
extern int  esl_tree_WriteNewick(FILE *fp, ESL_TREE *T);
extern int  esl_tree_ReadNewick(FILE *fp, char *errbuf, ESL_TREE **ret_T);

/* 3. Tree comparison algorithms.
 */
extern int esl_tree_Compare(ESL_TREE *T1, ESL_TREE *T2);

/* 4. Clustering algorithms for distance-based tree construction.
 */
extern int esl_tree_UPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T);
extern int esl_tree_WPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T);
extern int esl_tree_SingleLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T);
extern int esl_tree_CompleteLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T);

/* 5. Generating simulated trees.
 */
extern int esl_tree_Simulate(ESL_RANDOMNESS *r, int N, ESL_TREE **ret_T);
extern int esl_tree_ToDistanceMatrix(ESL_TREE *T, ESL_DMATRIX **ret_D);


#endif /*!ESL_TREE_INCLUDED*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/

