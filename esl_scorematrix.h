/* Routines for manipulating sequence alignment score matrices.
 * 
 * SRE, Mon Apr  2 08:33:23 2007 [Janelia]
 * SVN $Id$
 */
#ifndef ESL_SCOREMATRIX_INCLUDED
#define ESL_SCOREMATRIX_INCLUDED

#include <esl_alphabet.h>
#include <esl_dmatrix.h>
#include <esl_random.h>

/* 
 * allocation is in one array in s[0].
 *
 * i,j can range from 0..Kp-1, including all characters valid in the alphabet.
 * Only values for 0..K-1 (canonical alphabet) are mandatory.
 */
typedef struct {
  int **s;			/* s[i][j] is the score of aligning residue i,j; i,j range 0..Kp-1 */
  int   K;			/* size of base alphabet (duplicate of S->abc_r->K) */
  int   Kp;			/* full size of s[][], including degeneracies (duplicate of S->abc_r->Kp) */

  /* bookkeeping for degenerate residues */
  char *isval;			/* 0..Kp-1: which residues of alphabet have valid scores in S. */
  const ESL_ALPHABET *abc_r;	/* reference to the alphabet: includes K, Kp, and sym order */

  /* bookkeeping that lets us output exactly the residue order we read in a matrix file */
  int   nc;			/* number of residues with scores (inclusive of *, if present) */
  char *outorder;		/* string 0..nc-1 giving order of residues in col/row labels   */
  int   has_stop;		/* TRUE if * is a residue */
  int   stopsc;			/* score for alignment to a * */
  int   stopstopsc;		/* score for a *-* alignment  */
} ESL_SCOREMATRIX;




/* 1. The ESL_SCOREMATRIX object. */
extern ESL_SCOREMATRIX *esl_scorematrix_Create(const ESL_ALPHABET *abc);
extern int              esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, const double lambda, const double t);
extern int              esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, const double lambda, const ESL_DMATRIX *P,
						     const double *fi, const double *fj);
extern int              esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2);
extern int              esl_scorematrix_Max(const ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_Min(const ESL_SCOREMATRIX *S);
extern void             esl_scorematrix_Destroy(ESL_SCOREMATRIX *S);

/* 2. Reading/writing score matrices. */
extern int  esl_scorematrix_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S);
extern int  esl_scorematrix_Write(FILE *fp, const ESL_SCOREMATRIX *S);

/* 3. Interpreting score matrices probabilistically. */
extern int esl_scorematrix_ObtainPij(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, const int lambda, ESL_DMATRIX *P);
extern int esl_scorematrix_SolveLambda(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, ESL_DMATRIX *P, double *ret_lambda);
extern int esl_scorematrix_ReverseEngineer(const ESL_SCOREMATRIX *S, ESL_DMATRIX *P, double *fi, double *fj, double *ret_lambda);


#endif /*ESL_SCOREMATRIX_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/ 



