/* Routines for manipulating sequence alignment score matrices,
 * such as the BLOSUM and PAM matrices.
 * 
 * Contents:
 *   1. The ESL_SCOREMATRIX object.
 *   2. Reading/writing score matrices.
 *   3. Interpreting score matrices probabilistically.
 *   4. Utility programs.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example program.
 *   8. License and copyright.
 * 
 * SRE, Mon Apr  2 08:25:05 2007 [Janelia]
 * SVN $Id$
 */

#include <esl_config.h>

#include <string.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_fileparser.h>
#include <esl_scorematrix.h>

/*****************************************************************
 * 1. The ESL_SCOREMATRIX object
 *****************************************************************/

/* Function:  esl_scorematrix_Create()
 * Incept:    SRE, Mon Apr  2 08:38:10 2007 [Janelia]
 *
 * Purpose:   Allocates a score matrix for alphabet <abc>, initializes
 *            all scores to zero.
 *
 * Args:      abc   - pointer to digital alphabet 
 *
 * Returns:   a pointer to the new object.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SCOREMATRIX *
esl_scorematrix_Create(ESL_ALPHABET *abc)
{
  int status;
  int i;
  ESL_SCOREMATRIX *S = NULL;

  ESL_ALLOC(S, sizeof(ESL_SCOREMATRIX));
  S->s          = NULL;
  S->isval      = NULL;
  S->abc_r      = abc;
  S->nc         = 0;
  S->outorder   = NULL;
  S->has_stop   = FALSE;
  S->stopsc     = 0;
  S->stopstopsc = 0;

  ESL_ALLOC(S->s, sizeof(int *) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->s[i] = NULL;
  ESL_ALLOC(S->isval, sizeof(char) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->isval[i] = FALSE;

  ESL_ALLOC(S->s[0], sizeof(int) * abc->Kp * abc->Kp);
  for (i = 1; i < abc->Kp; i++) S->s[i] = S->s[0] + abc->Kp * i;

  for (i = 0; i < abc->Kp*abc->Kp; i++) S->s[0][i] = 0;
  return S;

 ERROR:
  esl_scorematrix_Destroy(S);
  return NULL;
}

/* Function:  esl_scorematrix_SetBLOSUM62
 * Incept:    SRE, Tue Apr  3 13:22:03 2007 [Janelia]
 *
 * Purpose:   Set the 20x20 canonical residue scores in an 
 *            allocated amino acid score matrix <S> to BLOSUM62
 *            scores \citep{Henikoff92}.
 *
 * Returns:   <eslOK> on success, and the scores in <S> are set.
 */
int
esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S)
{
  int status;
  int x,y;
  static int blosum62[28][28] = {
    /*  A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    -    B    J    Z    O    U    X    ~  */
    {   4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2,   0,  -2,   0,  -1,   0,   0,   0,   0,  },
    {   0,   9,  -3,  -4,  -2,  -3,  -3,  -1,  -3,  -1,  -1,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -2,  -2,   0,  -3,   0,  -3,   0,   0,  -2,   0,  },
    {  -2,  -3,   6,   2,  -3,  -1,  -1,  -3,  -1,  -4,  -3,   1,  -1,   0,  -2,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,   0,  },
    {  -1,  -4,   2,   5,  -3,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,   0,  },
    {  -2,  -2,  -3,  -3,   6,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   1,   3,   0,  -3,   0,  -3,   0,   0,  -1,   0,  },
    {   0,  -3,  -1,  -2,  -3,   6,  -2,  -4,  -2,  -4,  -3,   0,  -2,  -2,  -2,   0,  -2,  -3,  -2,  -3,   0,  -1,   0,  -2,   0,   0,  -1,   0,  },
    {  -2,  -3,  -1,   0,  -1,  -2,   8,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2,   0,   0,   0,   0,   0,   0,  -1,   0,  },
    {  -1,  -1,  -3,  -3,   0,  -4,  -3,   4,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -3,  -1,   0,  -3,   0,  -3,   0,   0,  -1,   0,  },
    {  -1,  -3,  -1,   1,  -3,  -2,  -1,  -3,   5,  -2,  -1,   0,  -1,   1,   2,   0,  -1,  -2,  -3,  -2,   0,   0,   0,   1,   0,   0,  -1,   0,  },
    {  -1,  -1,  -4,  -3,   0,  -4,  -3,   2,  -2,   4,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1,   0,  -4,   0,  -3,   0,   0,  -1,   0,  },
    {  -1,  -1,  -3,  -2,   0,  -3,  -2,   1,  -1,   2,   5,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1,   0,  -3,   0,  -1,   0,   0,  -1,   0,  },
    {  -2,  -3,   1,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2,   0,   3,   0,   0,   0,   0,  -1,   0,  },
    {  -1,  -3,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   7,  -1,  -2,  -1,  -1,  -2,  -4,  -3,   0,  -2,   0,  -1,   0,   0,  -2,   0,  },
    {  -1,  -3,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   5,   1,   0,  -1,  -2,  -2,  -1,   0,   0,   0,   3,   0,   0,  -1,   0,  },
    {  -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   2,  -2,  -1,   0,  -2,   1,   5,  -1,  -1,  -3,  -3,  -2,   0,  -1,   0,   0,   0,   0,  -1,   0,  },
    {   1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -2,  -2,   0,  -1,   0,  -1,   0,   0,   0,   0,  },
    {   0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -2,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -3,  -1,   0,  -3,   0,  -2,   0,   0,  -1,   0,  },
    {  -3,  -2,  -4,  -3,   1,  -2,  -2,  -3,  -3,  -2,  -1,  -4,  -4,  -2,  -3,  -3,  -2,  -3,  11,   2,   0,  -4,   0,  -3,   0,   0,  -2,   0,  },
    {  -2,  -2,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   7,   0,  -3,   0,  -2,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {  -2,  -3,   4,   1,  -3,  -1,   0,  -3,   0,  -4,  -3,   3,  -2,   0,  -1,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {  -1,  -3,   1,   4,  -3,  -2,   0,  -3,   1,  -3,  -1,   0,  -1,   3,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,  -1,  -1,   0,   0,  -1,  -2,  -1,   0,  -1,   0,  -1,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
  };

  for (x = 0;           x < S->abc_r->K;  x++) S->isval[x] = TRUE;
  for (x = S->abc_r->K; x < S->abc_r->Kp; x++) S->isval[x] = FALSE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'B');   S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'Z');   S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'X');   S->isval[x] = TRUE;
    
  for (x = 0; x < S->abc_r->Kp; x++)
    for (y = 0; y < S->abc_r->Kp; y++)
      S->s[x][y] = blosum62[x][y];

  /* Bookkeeping necessary to be able to reproduce BLOSUM62 output format exactly, if we need to Write() */
  if ((status = esl_strdup("ARNDCQEGHILKMFPSTWYVBZX*", -1, &(S->outorder))) != eslOK) return status;
  S->nc         = strlen(S->outorder);
  S->has_stop   = TRUE;
  S->stopsc     = -4;
  S->stopstopsc = 1;

  return eslOK;
}


/* Function:  esl_scorematrix_Compare()
 * Incept:    SRE, Tue Apr  3 14:17:12 2007 [Janelia]
 *
 * Purpose:   Compares two score matrices; returns <eslOK> if they 
 *            are identical, <eslFAIL> if they differ.
 */
int
esl_scorematrix_Compare(ESL_SCOREMATRIX *S1, ESL_SCOREMATRIX *S2)
{
  int a,b;

  if (strcmp(S1->outorder, S2->outorder) != 0) return eslFAIL;
  if (S1->nc         != S2->nc)                return eslFAIL;
  if (S1->has_stop   != S2->has_stop)          return eslFAIL;
  if (S1->stopsc     != S2->stopsc)            return eslFAIL;
  if (S1->stopstopsc != S2->stopstopsc)        return eslFAIL;
  
  for (a = 0; a < S1->nc; a++)
    if (S1->isval[a] != S2->isval[a])          return eslFAIL;
  
  for (a = 0; a < S1->abc_r->Kp; a++)
    for (b = 0; b < S1->abc_r->Kp; b++)
      if (S1->s[a][b] != S2->s[a][b])          return eslFAIL;
  return eslOK;
}



/* Function:  esl_scorematrix_Destroy()
 * Incept:    SRE, Mon Apr  2 08:46:44 2007 [Janelia]
 *
 * Purpose:   Frees a score matrix.
 *
 * Returns:   (void).
 */
void
esl_scorematrix_Destroy(ESL_SCOREMATRIX *S)
{
  if (S == NULL) return;
  if (S->s != NULL) {
    if (S->s[0] != NULL) free(S->s[0]);
    free(S->s);
  }
  if (S->isval    != NULL) free(S->isval);
  if (S->outorder != NULL) free(S->outorder);
  return;
}


/*****************************************************************
 * 2. Reading/writing score matrices.
 *****************************************************************/

/* Function:  esl_scorematrix_Read()
 * Incept:    SRE, Mon Apr  2 08:26:40 2007 [Janelia]
 *
 * Purpose:   Given a pointer <efp> to an open file parser for a file
 *            containing a score matrix (such as a PAM or BLOSUM
 *            matrix), parse the file and create a new score matrix
 *            object. The scores are expected to be for the alphabet
 *            <abc>. 
 *            
 *            The score matrix file is in the format that BLAST or
 *            FASTA use. The first line is a header contains N
 *            single-letter codes for the residues. Each of N
 *            subsequent rows optionally contains a residue row label
 *            (in the same order as the columns), followed by N
 *            residue scores.  (Older matrix files do not contain the
 *            leading row label; newer ones do.) The residues may
 *            appear in any order. They must minimally include the
 *            canonical K residues (K=4 for DNA, K=20 for protein),
 *            and may also contain none, some, or all degeneracy
 *            codes. Any other residue code that is not in the Easel
 *            digital alphabet (including, in particular, the '*' code
 *            for a stop codon) is ignored by the parser.
 *
 * Returns:   <eslOK> on success, and <ret_S> points to a newly allocated 
 *            score matrix. 
 *
 *            Returns <eslEFORMAT> on parsing error; in which case, <ret_S> is
 *            returned <NULL>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_scorematrix_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S)
{
  int status;
  ESL_SCOREMATRIX *S     = NULL;
  int             *map   = NULL; /* maps col/row index to digital alphabet x */
  char            *tmp;
  int              lalloc;
  char            *tok;
  int              toklen;
  int              nc    = 0;
  int              c, x;
  int              row,col;

  /* Allocate the matrix
   */
  if ((S = esl_scorematrix_Create(abc)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Make sure we've got the comment character set properly in the fileparser.
   * Score matrices use #.
   */
  esl_fileparser_SetCommentChar(efp, '#');

  /* Look for the first non-blank, non-comment line in the file.  That line
   * gives us the single-letter codes in the order that the file's using.
   */
  if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "file appears to be empty");

  /* Read the characters: count them and store them in order in label[0..nc-1].
   */
  lalloc = 32;
  ESL_ALLOC(S->outorder, sizeof(char) * lalloc);
  while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) == eslOK)
    {
      if (toklen != 1) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Header can only contain single-char labels; %s is invalid", tok);
      S->outorder[nc++] = *tok;
      if (nc == lalloc) {
	lalloc *= 2;
	ESL_RALLOC(S->outorder, tmp, sizeof(char) * lalloc);
      }
    }
  if (status != eslEOL) ESL_XFAIL(status, efp->errbuf, "Unexpected failure of esl_fileparser_GetTokenOnLine()");
  S->nc = nc;
  
  /* Verify that these labels for the score matrix seem plausible, given our alphabet.
   * This sets S->isval array: which residues we have scores for.
   * It also sets the map[] array, which maps coord in label[] to x in alphabet.
   * It's possible to see a residue in the score matrix that's not in the alphabet (main example is '*', a stop codon)
   */
  ESL_ALLOC(map, sizeof(int) * nc);
  for (c = 0; c < nc; c++)
    {
      if (esl_abc_CIsValid(abc, S->outorder[c])) 
	{  
	  x = esl_abc_DigitizeSymbol(abc, S->outorder[c]);
	  map[c] = x;
	  S->isval[x] = TRUE;
	}
      else if (S->outorder[c] == '*')
	{
	  S->has_stop = TRUE;
	  map[c] = -1;
	}
      else
	ESL_XFAIL(eslEFORMAT, efp->errbuf, "Don't know how to deal with residue %c in matrix file", S->outorder[c]);
    }
  for (x = 0; x < abc->K; x++)
    if (! S->isval[x]) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected to see a column for residue %c", abc->sym[x]);


  /* Read nc rows, one at a time;
   * on each row, read nc+1 or nc tokens, of which nc are scores (may lead with a label or not)
   */
  for (row = 0; row < nc; row++)
    {
      if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of lines in file");
      for (col = 0; col < nc; col++)
	{
	  if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of fields on line");
	  if (col == 0 && *tok == S->outorder[row]) { col--; continue; } /* skip leading label */

	  if (map[row] >= 0 && map[col] >= 0)
	    S->s[map[row]][map[col]] = atoi(tok);
	  else if (map[row] == -1 && map[col] == -1) /* stop/stop alignment */
	    S->stopstopsc = atoi(tok);
	  else 
	    S->stopsc = atoi(tok); /* this'll reset the stop score 2*nc-1 times, wastefully, and assuming they're all identical */
	}
      if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslEOL)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many fields on line");
    }
  if ((status = esl_fileparser_NextLine(efp)) != eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many lines in file");
  

  free(map);
  *ret_S = S;
  return eslOK;

 ERROR:
  esl_scorematrix_Destroy(S);
  if (map != NULL) free(map);
  *ret_S = NULL;
  return status;
}

/* Function:  esl_scorematrix_Write()
 * Incept:    SRE, Tue Apr  3 13:55:10 2007 [Janelia]
 *
 * Purpose:   Writes a score matrix <S> to an open stream <fp>, in 
 *            format compatible with BLAST, FASTA, and other common
 *            sequence alignment software.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_scorematrix_Write(FILE *fp, ESL_SCOREMATRIX *S)
{
  int a,b;			
  int x,y;
  int nc = 0;
  
  /* Total paranoia: we have two redundant ways to determine the
   * number of residues in this matrix, and they should match:
   * S->nc, or the sum of the isval[] flags + has_stop.
   */
  if (S->has_stop) nc++;
  for (x = 0; x < S->abc_r->Kp; x++)
    if (S->isval[x]) nc++;
  if (nc != S->nc) ESL_EXCEPTION(eslEINVAL, "nc's don't match. matrix is corrupt");

  /* The header line, with column labels for residues */
  fprintf(fp, "  ");
  for (a = 0; a < nc; a++) fprintf(fp, "  %c ", S->outorder[a]);
  fprintf(fp, "\n");
  
  /* The data. Watch out for those pesky *'s, which aren't in the Easel digital alphabet (yet)
   */
  for (a = 0; a < nc; a++)
    {
      fprintf(fp, "%c ", S->outorder[a]);
      for (b = 0; b < nc; b++)
	{
	  if (S->outorder[a] != '*' && S->outorder[b] != '*') 
	    {
	      x = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[a]);
	      y = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[b]);
	      fprintf(fp, "%3d ", S->s[x][y]);
	    } 
	  else if (S->outorder[a] != '*' || S->outorder[b] != '*')
	    fprintf(fp, "%3d ", S->stopsc);
	  else
	    fprintf(fp, "%3d ", S->stopstopsc);
	}
      fprintf(fp, "\n");
    }
  
  return eslOK;
}

/*****************************************************************
 * 3. Interpreting score matrices probabilistically.
 *****************************************************************/ 








/*****************************************************************
 * 4. Utilities
 *****************************************************************/ 

/* Reformat a score matrix file, canonical residues only, into
 * Easel internal digital alphabet order, suitable for making 
 * a static data structure.
 */
#ifdef eslSCOREMATRIX_UTILITY1
/* 
    gcc -g -Wall -o utility -I. -L. -DeslSCOREMATRIX_UTILITY1 esl_scorematrix.c -leasel -lm
    ./utility BLOSUM62
*/
#include <easel.h>
#include <esl_alphabet.h>
#include <esl_scorematrix.h>
#include <esl_fileparser.h>

int
main(int argc, char **argv)
{
  char *infile = argv[1];
  ESL_ALPHABET    *abc;
  ESL_FILEPARSER  *efp;
  ESL_SCOREMATRIX *S;
  int x,y;

  abc = esl_alphabet_Create(eslAMINO);

  if (esl_fileparser_Open(infile, &efp)  != eslOK) esl_fatal("Failed to open %s\n", infile);
  if (esl_scorematrix_Read(efp, abc, &S) != eslOK) esl_fatal("parse failed: %s", efp->errbuf);

  for (x = 0; x < abc->Kp; x++) {
    printf("{ ");
    for (y = 0; y < abc->Kp; y++)
      printf("%3d, ", S->s[x][y]);
    printf(" },\n");
  }
  
  esl_scorematrix_Destroy(S);
  esl_fileparser_Close(efp);
  esl_alphabet_Destroy(abc);
  return eslOK;
}
#endif /*eslSCOREMATRIX_UTILITY1*/


/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/

#ifdef eslSCOREMATRIX_TESTDRIVE

static void
utest_ReadWrite(void)
{
  char tmpfile[16]     = "esltmpXXXXXX";
  FILE            *fp  = NULL;
  ESL_ALPHABET    *abc = NULL;
  ESL_SCOREMATRIX *S   = NULL;
  ESL_SCOREMATRIX *S2  = NULL;
  ESL_FILEPARSER  *efp = NULL;
  
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("allocation of alphabet failed");
  if ((S   = esl_scorematrix_Create(abc))   == NULL)  esl_fatal("allocation of scorematrix failed");
  if (esl_scorematrix_SetBLOSUM62(S)        != eslOK) esl_fatal("failed to set mx to BLOSUM62");
  if (esl_tmpfile_named(tmpfile, &fp)       != eslOK) esl_fatal("failed to open tmp file");
  if (esl_scorematrix_Write(fp, S)          != eslOK) esl_fatal("failed to write BLOSUM62 matrix");
  fclose(fp);

  if (esl_fileparser_Open(tmpfile, &efp)    != eslOK) esl_fatal("failed to open tmpfile containing BLOSUM62 matrix");
  if (esl_scorematrix_Read(efp, abc, &S2)   != eslOK) esl_fatal("failed to read tmpfile containing BLOSUM62 matrix");
  if (esl_scorematrix_Compare(S, S2)        != eslOK) esl_fatal("the two test matrices aren't identical");
  
  /*  remove(tmpfile); */
  esl_fileparser_Close(efp);
  esl_scorematrix_Destroy(S2);
  esl_scorematrix_Destroy(S);
  esl_alphabet_Destroy(abc);
  return;
}

#endif /*eslSCOREMATRIX_TESTDRIVE*/


/*****************************************************************
 * 6. Test driver.
 *****************************************************************/

/* 
    gcc -g -Wall -I. -L. -o test -DeslSCOREMATRIX_TESTDRIVE esl_scorematrix.c -leasel -lm
    ./test
*/
#ifdef eslSCOREMATRIX_TESTDRIVE
#include <easel.h>
#include <esl_scorematrix.h>

int 
main(int argc, char **argv)
{
  utest_ReadWrite();

  return 0;
}
#endif /*eslSCOREMATRIX_TESTDRIVE*/

/*****************************************************************
 * 7. Example program
 *****************************************************************/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/ 
