# esl_tree - Phylogenetic trees

The `tree` module implements an object for representing
phylogenetic trees (`ESL_TREE`). It also implements four of the
most basic distance-based methods for tree inference and clustering
algorithms (UPGMA, WPGMA, single linkage clustering, and complete
linkage clustering).




## example

```c
/* To compile: gcc -g -Wall -o example -I. -DeslTREE_EXAMPLE esl_tree.c esl_dmatrix.c esl_msa.c easel.c -lm
 *         or: gcc -g -Wall -o example -I. -L. -DeslTREE_EXAMPLE esl_tree.c -leasel -lm
 *     To run: ./example <MSA file>
 */
#include "easel.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_distance.h"
#include "esl_tree.h"

int main(int argc, char **argv)
{
  ESL_TREE     *tree = NULL;
  ESL_MSAFILE  *afp  = NULL;
  ESL_MSA      *msa  = NULL;
  ESL_DMATRIX  *D    = NULL;

  esl_msafile_Open(NULL, argv[1], NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  esl_msafile_Read(afp, &msa);
  esl_msafile_Close(afp);

  esl_dst_CDiffMx(msa->aseq, msa->nseq, &D);
  esl_tree_UPGMA(D, &tree);

  esl_tree_Destroy(tree);
  esl_msa_Destroy(msa);
  esl_dmatrix_Destroy(D);
  return eslOK;
}
```
