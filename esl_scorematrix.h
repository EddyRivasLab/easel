/* Routines for manipulating sequence alignment score matrices.
 * 
 * SRE, Mon Apr  2 08:33:23 2007 [Janelia]
 * SVN $Id$
 */
#ifndef ESL_SCOREMATRIX_INCLUDED
#define ESL_SCOREMATRIX_INCLUDED

#include <esl_alphabet.h>

/* 
 * allocation is in one array in s[0].
 *
 * i,j can range from 0..Kp-1, including all characters valid in the alphabet.
 * Only values for 0..K-1 (canonical alphabet) are mandatory.
 */
typedef struct {
  int **s;			/* s[i][j] is the score of aligning residue i,j; i,j range 0..Kp-1 */
  char *isval;			/* 0..Kp-1: which residues of alphabet have valid scores. */
  ESL_ALPHABET *abc_r;		/* reference to the alphabet: includes K, Kp, and sym order */

  /* bookkeeping that lets us output exactly the residue order we read in a matrix file */
  int   nc;			/* number of residues with scores (inclusive of *, if present) */
  char *outorder;		/* string 0..nc-1 giving order of residues in col/row labels   */
  int   has_stop;		/* TRUE if * is a residue */
  int   stopsc;			/* score for alignment to a * */
  int   stopstopsc;		/* score for a *-* alignment  */
} ESL_SCOREMATRIX;




/* 1. The ESL_SCOREMATRIX object. */
extern ESL_SCOREMATRIX *esl_scorematrix_Create(ESL_ALPHABET *abc);
extern int              esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_Compare(ESL_SCOREMATRIX *S1, ESL_SCOREMATRIX *S2);
extern void             esl_scorematrix_Destroy(ESL_SCOREMATRIX *S);

/* 2. Reading/writing score matrices. */
extern int  esl_scorematrix_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S);
extern int  esl_scorematrix_Write(FILE *fp, ESL_SCOREMATRIX *S);

#endif /*ESL_SCOREMATRIX_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/ 



