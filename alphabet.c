/* alphabet.c
 * Implements the standard digitized alphabets for biosequences.
 * 
 * SVN $Id$
 * SRE, Tue Dec  7 13:49:43 2004
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/alphabet.h>

static void set_symmap(ESL_ALPHABET *a);
static void init_degenerate(ESL_ALPHABET *a);
static void set_degenerate(ESL_ALPHABET *a, char iupac, char *syms);


/* Function: esl_alphabet_Set()
 * 
 * Purpose:  Set the alphabet globals, given an alphabet <type>
 *           of <eslAMINO>, <eslDNA>, or <eslRNA>.
 *           
 * Returns:  <ESL_OK>
 * 
 * Throws:   <ESL_EINVAL> if <type> is not one of <eslAMINO>, <eslDNA>, 
 *           or <eslRNA>.
 */
int
esl_alphabet_Set(ESL_ALPHABET *a, int type)
{
  switch(type) { 
  case eslAMINO:  esl_alphabet_SetAmino(a); break;
  case eslDNA:    esl_alphabet_SetDNA(a);   break;
  case eslRNA:    esl_alphabet_SetRNA(a);   break;
  default:        ESL_ERROR(ESL_EINVAL, "No support for alphabets other than DNA, RNA, protein");
  }
  return ESL_OK;
}

/* Function:  esl_alphabet_SetDNA()
 * Incept:    SRE, Tue Dec  7 15:22:11 2004 [St. Louis]
 *
 * Purpose:   Sets alphabet <a> to the 4/16 DNA code.
 *
 * Returns:   <ESL_OK>
 */
int
esl_alphabet_SetDNA(ESL_ALPHABET *a)
{
  
  a->type = ESL_DNA;
  a->K    = 4;
  a->Kall = 17;
  strcpy(a->sym, "ACGTURYMKSWHBVDXN");

  set_symmap(a);

  init_degenerate(a);
  set_degenerate(a, 'U', "T");
  set_degenerate(a, 'R', "AG");
  set_degenerate(a, 'Y', "CT");
  set_degenerate(a, 'M', "AC");
  set_degenerate(a, 'K', "GT");
  set_degenerate(a, 'S', "CG");
  set_degenerate(a, 'W', "AT");
  set_degenerate(a, 'H', "ACT");
  set_degenerate(a, 'B', "CGT");
  set_degenerate(a, 'V', "ACG");
  set_degenerate(a, 'D', "AGT");  
  set_degenerate(a, 'X', "ACGT");
  set_degenerate(a, 'N', "ACGT");

  return ESL_OK;
}

/* Function:  esl_alphabet_SetRNA()
 * Incept:    SRE, Tue Dec  7 15:36:46 2004 [St. Louis]
 *
 * Purpose:   Sets alphabet <a> to the 4/16 RNA code.
 *
 * Returns:   <ESL_OK>
 */
int
esl_alphabet_SetRNA(ESL_ALPHABET *a)
{
  a->type = ESL_RNA;
  a->K    = 4;
  a->Kall = 17;
  strcpy(a->sym, "ACGUTRYMKSWHBVDXN");

  set_symmap(a);

  init_degenerate(a);
  set_degenerate(a, 'T', "U");
  set_degenerate(a, 'R', "AG");
  set_degenerate(a, 'Y', "CT");
  set_degenerate(a, 'M', "AC");
  set_degenerate(a, 'K', "GT");
  set_degenerate(a, 'S', "CG");
  set_degenerate(a, 'W', "AT");
  set_degenerate(a, 'H', "ACT");
  set_degenerate(a, 'B', "CGT");
  set_degenerate(a, 'V', "ACG");
  set_degenerate(a, 'D', "AGT");  
  set_degenerate(a, 'X', "ACGT");
  set_degenerate(a, 'N', "ACGT");

  return ESL_OK;
}


/* Function:  esl_alphabet_SetAmino()
 * Incept:    SRE, Tue Dec  7 15:01:55 2004 [St. Louis]
 *
 * Purpose:   Sets alphabet <a> to the 20/23 amino acid code.
 *
 * Returns:   <ESL_OK>
 */
int
esl_alphabet_SetAmino(ESL_ALPHABET *a)
{
  int i,x,y;

  a->type = ESL_AMINO;
  a->K    = 20;
  a->Kall = 24;
  strcpy(a->sym, "ACDEFGHIKLMNPQRSTVWYBZUX");

  set_symmap(a);
  
  init_degenerate(a);
  set_degenerate(a, 'U', "S"); /* treat Sec as Ser */
  set_degenerate(a, 'B', "ND");
  set_degenerate(a, 'Z', "QE");
  set_degenerate(a, 'X', "ACDEFGHIKLMNPQRSTVWY");  

  return ESL_OK;
}

/* Function: esl_alph_DigitizeSequence()
 * 
 * Purpose:  Internal representation of a sequence in Easel is
 *           as an unsigned char array. 1..L are the indices 
 *           of seq symbols in Alphabet[]. 0,L+1 are sentinel
 *           bytes, set to be Alphabet_iupac -- i.e. one more
 *           than the maximum allowed index.  
 *           
 *           Assumes that 'N' or 'X', the fully degenerate characters
 *           for DNA/RNA or protein, respectively, is the last
 *           character in the allowed alphabet.
 *           
 * Args:     seq - sequence to be digitized (0..L-1)
 *           L   - length of sequence      
 *           
 * Return:   digitized sequence, dsq.
 *           dsq is allocated here and must be free'd by caller.
 */
unsigned char *
DigitizeSequence(char *seq, int L)
{
  unsigned char *dsq;
  int i;

  /* TODO: make this malloc more robust. */
  if ((dsq = malloc (sizeof(unsigned char) * (L+2))) == NULL) return NULL;
  dsq[0] = dsq[L+1] = (unsigned char) Alphabet_iupac;
  for (i = 1; i <= L; i++) 
    dsq[i] = SymbolIndex(seq[i-1]);
  return dsq;
}


/* set_symmap()
 * 
 * Builds the symbol map, which maps characters in a seq
 * to digitized indices 0..Kall-1 in a dsq, case-insensitively.
 * a->sym and a->Kall must be set in the alphabet <a>.
 */
static void
set_symmap(ESL_ALPHABET *a)
{
  int i,x;

  for (i = 0; i < 128; i++) a->symmap[i] = -1;
  for (i = 0; i < a->Kall; i++) {
    x = a->sym[i];          a->symmap[x] = i;
    x = tolower(a->sym[i]); a->symmap[x] = i;
  }
}

/* init_degenerate()
 * 
 * Initialize the degen and ndegen fields of an alphabet
 * to all zeros, prior to starting to set them appropriately;
 * and for nondegenerate characters, set each to its unique
 * single 1's.
 */
static void
init_degenerate(ESL_ALPHABET *a)
{
  int x;

  for (x = 0; x < ESL_MAXCODE; x++)
    {
      a->ndegen[x] = 0;
      for (y = 0; y < ESL_MAXABET; y++)
	a->degen[x][y] = 0;
    }
  for (x = 0; x < a->K; x++) {
    a->degen[x][x] = 1;
    a->ndegen[x]   = 1;
  }
}

/* set_degenerate()
 * 
 * Convenience function for setting up degen and ndegen fields of an
 * alphabet.  Given a degenerate symbol <iupac>, and a string of
 * unique characters that it represents, fill in the alphabet <a>
 * appropriately.
 */
static void 
set_degenerate(ESL_ALPHABET *a, char iupac, char *syms)
{
  int x,y;

  x = strchr(a->sym, iupac)-a->sym;
  a->ndegen[x] = strlen(syms);

  while (*syms) {
    y = strchr(a->sym, *syms)-a->sym;
    a->degen[x][y] = 1;
    syms++;
  }
}




