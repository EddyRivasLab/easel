/* alphabet.c
 * Implements the standard digitized alphabets for biosequences.
 * (See alphabet.tex for API and other documentation.)
 * 
 * SVN $Id$
 * SRE, Tue Dec  7 13:49:43 2004
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>

static ESL_ALPHABET *create_dna(void);
static ESL_ALPHABET *create_rna(void);
static ESL_ALPHABET *create_amino(void);


/* Function:  esl_alphabet_Create()
 * Incept:    SRE, Mon Dec 20 10:21:54 2004 [Zaragoza]
 *
 * Purpose:   Creates one of the three standard bio alphabets:
 *            <eslDNA>, <eslRNA>, or <eslAMINO>, and returns
 *            a pointer to it.
 *
 * Args:      type  - <eslDNA>, <eslRNA>, or <eslAMINO>. 
 *
 * Returns:   pointer to the new alphabet.
 *
 * Throws:    NULL if any allocation or initialization fails.
 */
ESL_ALPHABET *
esl_alphabet_Create(int type)
{
  ESL_ALPHABET *a;

  switch(type) { 
  case eslAMINO:  a = create_amino(); break;
  case eslDNA:    a = create_dna();   break;
  case eslRNA:    a = create_rna();   break;
  default:    
    ESL_ERROR_NULL(eslEINVAL,
		   "Standard alphabets include only DNA, RNA, protein.");
  }
  return a;
}


/* Function:  esl_alphabet_CreateCustom()
 * Incept:    SRE, Mon Dec 20 09:18:28 2004 [Zaragoza]
 *
 * Purpose:   Creates a customized biosequence alphabet,
 *            and returns a ptr to it. The alphabet type is set 
 *            to eslNONSTANDARD.
 *            
 *            <alphabet> is the internal alphabet string;
 *            <K> is the size of the base alphabet;
 *            <Kp> is the total size of the alphabet string. 
 *            
 *            In the alphabet string, residues 0..K-1 are the base alphabet; 
 *            residue K is the canonical gap (indel) symbol; 
 *            residue Kp-1 is an "any" symbol (such as N or X); 
 *            and residues K+1..Kp-1 are additional degeneracy symbols.
 *            The gap and the "any" symbol are mandatory even for
 *            nonstandard alphabets, so $Kp >= K+2$.
 *            
 * Args:      alphabet - internal alphabet; example "ACGT-RYMKSWHBVDN"
 *            K        - base size; example 4
 *            Kp       - total size; example 16
 *
 * Returns:   pointer to new <ESL_ALPHABET> structure
 *
 * Throws:    NULL if any allocation or initialization fails.
 */
ESL_ALPHABET *
esl_alphabet_CreateCustom(char *alphabet, int K, int Kp)
{
  ESL_ALPHABET *a;
  int           c,x,y;

  /* Argument checks.
   */
  if (strlen(alphabet) != Kp) 
    ESL_ERROR_NULL(eslEINVAL, "alphabet length != Kp");
  if (Kp < K+2)               
    ESL_ERROR_NULL(eslEINVAL, "Kp too small in alphabet"); 

  /* Allocation/init, level 1.
   */
  if ((a = malloc(sizeof(ESL_ALPHABET))) == NULL) 
    ESL_ERROR_NULL(eslEMEM, "malloc failed");
  a->sym    = NULL;
  a->degen  = NULL;
  a->ndegen = NULL;
  
  /* Allocation/init, level 2.
   */
  if ((a->sym    = malloc(sizeof(char)   * (Kp+1))) == NULL) 
    { esl_alphabet_Destroy(a); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((a->degen  = malloc(sizeof(char *) * Kp))     == NULL)
    { esl_alphabet_Destroy(a); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((a->ndegen = malloc(sizeof(int)    * Kp))     == NULL)
    { esl_alphabet_Destroy(a); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  a->degen[0] = NULL;

  /* Allocation/init, level 3.
   */
  if ((a->degen[0] = malloc(sizeof(char) * (Kp*K))) == NULL) 
    { esl_alphabet_Destroy(a); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  for (x = 1; x < Kp; x++)
    a->degen[x] = a->degen[0]+(K*x);

  /* Initialize the internal alphabet: 
   */
  a->type = eslNONSTANDARD;
  a->K    = K;
  a->Kp   = Kp;
  strcpy(a->sym, alphabet);

  /* Initialize the input map;
   * note duplication of similar work in sqio.c, which we should
   * someday abstract away and standardize.
   */
  for (c = 0; c < 128; c++)   a->inmap[c]               = ESL_ILLEGAL_CHAR;
  for (x = 0; x < a->Kp; x++) a->inmap[(int) a->sym[x]] = x;  

  /* Initialize the degeneracy map:
   *  Base alphabet (first K syms) are automatically
   *  mapped uniquely; last character (Kp-1) is assumed to be
   *  the "any" character; other degen chars (K+1..Kp-2) are 
   *  unset; gap character is unused.
   */
  for (x = 0; x < a->Kp; x++)  	/* clear everything */
    {
      a->ndegen[x] = 0;
      for (y = 0; y < a->K; y++) a->degen[x][y] = 0;
    }
  for (x = 0; x < a->K; x++) 	/* base alphabet */
    {
      a->ndegen[x]   = 1;
      a->degen[x][x] = 1;
    }
                                /* "any" character */
  a->ndegen[Kp-1]  = K;
  for (x = 0; x < a->K; x++) a->degen[Kp-1][x] = 1;

  /* Successful return.
   * (Caller still must set degeneracies and synonyms.)
   */
  return a;
}



/* create_dna(): 
 * creates and returns a standard DNA alphabet.
 */
static ESL_ALPHABET *
create_dna(void)
{
  ESL_ALPHABET *a;

  /* Create the fundamental alphabet.
   */
  if ((a = esl_alphabet_CreateCustom("ACGT-RYMKSWHBVDN", 4, 16)) == NULL) return NULL;
  a->type = eslDNA;
  
  /* Add desired synonyms in the input map.
   */
  esl_alphabet_SetSynonym(a, 'U', 'T');	    /* read U as a T */
  esl_alphabet_SetSynonym(a, 'X', 'N');	    /* read X as an N (many seq maskers use X) */
  esl_alphabet_SetSynonym(a, '_', '-');     /* allow _ as a gap too */
  esl_alphabet_SetSynonym(a, '.', '-');     /* allow . as a gap too */
  esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input */

  /* Define IUBMB degenerate symbols other than the N.
   */
  esl_alphabet_SetDegeneracy(a, 'R', "AG");
  esl_alphabet_SetDegeneracy(a, 'Y', "CT");
  esl_alphabet_SetDegeneracy(a, 'M', "AC");
  esl_alphabet_SetDegeneracy(a, 'K', "GT");
  esl_alphabet_SetDegeneracy(a, 'S', "CG");
  esl_alphabet_SetDegeneracy(a, 'W', "AT");
  esl_alphabet_SetDegeneracy(a, 'H', "ACT");
  esl_alphabet_SetDegeneracy(a, 'B', "CGT");
  esl_alphabet_SetDegeneracy(a, 'V', "ACG");
  esl_alphabet_SetDegeneracy(a, 'D', "AGT");  

  return a;
}


/* create_rna(): 
 * Creates a standard RNA alphabet.
 */
static ESL_ALPHABET *
create_rna(void)
{
  ESL_ALPHABET *a;

  /* Create the fundamental alphabet
   */
  if ((a = esl_alphabet_CreateCustom("ACGU-RYMKSWHBVDN", 4, 16)) == NULL) return NULL;
  a->type = eslRNA;
  
  /* Add desired synonyms in the input map.
   */
  esl_alphabet_SetSynonym(a, 'T', 'U');	    /* read T as a U */
  esl_alphabet_SetSynonym(a, 'X', 'N');	    /* read X as an N (many seq maskers use X) */
  esl_alphabet_SetSynonym(a, '_', '-');     /* allow _ as a gap too */
  esl_alphabet_SetSynonym(a, '.', '-');     /* allow . as a gap too */
  esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input */

  
  /* Define IUBMB degenerate symbols other than the N.
   */
  esl_alphabet_SetDegeneracy(a, 'R', "AG");
  esl_alphabet_SetDegeneracy(a, 'Y', "CU");
  esl_alphabet_SetDegeneracy(a, 'M', "AC");
  esl_alphabet_SetDegeneracy(a, 'K', "GU");
  esl_alphabet_SetDegeneracy(a, 'S', "CG");
  esl_alphabet_SetDegeneracy(a, 'W', "AU");
  esl_alphabet_SetDegeneracy(a, 'H', "ACU");
  esl_alphabet_SetDegeneracy(a, 'B', "CGU");
  esl_alphabet_SetDegeneracy(a, 'V', "ACG");
  esl_alphabet_SetDegeneracy(a, 'D', "AGU");  

  return a;
}


/* create_amino():
 * Creates a new standard amino acid alphabet.
 */
static ESL_ALPHABET *
create_amino(void)
{
  ESL_ALPHABET *a;

  /* Create the fundamental alphabet
   */
  if ((a = esl_alphabet_CreateCustom("ACDEFGHIKLMNPQRSTVWY-BZX", 20, 24)) == NULL) return NULL;
  a->type = eslAMINO;
  
  /* Add desired synonyms in the input map.
   */
  esl_alphabet_SetSynonym(a, 'U', 'S');	    /* read SelCys U as a serine S */
  esl_alphabet_SetSynonym(a, '_', '-');     /* allow _ as a gap too */
  esl_alphabet_SetSynonym(a, '.', '-');     /* allow . as a gap too */
  esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input */
  
  /* Define IUPAC degenerate symbols other than the X.
   */
  esl_alphabet_SetDegeneracy(a, 'B', "ND");
  esl_alphabet_SetDegeneracy(a, 'Z', "QE");

  return a;
}



/* Function:  esl_alphabet_Destroy()
 * Incept:    SRE, Mon Dec 20 10:27:23 2004 [Zaragoza]
 *
 * Purpose:   Free's an <ESL_ALPHABET> structure.
 *
 * Args:      a  - the <ESL_ALPHABET> to free.
 *
 * Returns:   (void).
 */
void
esl_alphabet_Destroy(ESL_ALPHABET *a)
{
  if (a == NULL) return;

  if (a->sym      != NULL) free(a->sym);
  if (a->ndegen   != NULL) free(a->ndegen);
  if (a->degen    != NULL) 
    {
      if (a->degen[0] != NULL) free(a->degen[0]);
      free(a->degen);
    }
  free(a);
}




/* Function:  esl_alphabet_SetSynonym()
 * Incept:    SRE, Mon Dec 20 10:40:33 2004 [Zaragoza]
 *
 * Purpose:   Maps an additional input alphabetic symbol <sym> to 
 *            an internal alphabet symbol <c>; for example,
 *            we might map T to U for an RNA alphabet, so that we
 *            allow for reading input DNA sequences.
 *
 * Args:      sym   - symbol in the input alphabet; 'T' for example
 *            c     - symbol in the internal alphabet; 'U' for example
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <c> is not in the internal alphabet.
 */
int
esl_alphabet_SetSynonym(ESL_ALPHABET *a, char sym, char c)
{
  char *sp;
  int   x;

  if ((sp = strchr(a->sym, c)) == NULL) 
    ESL_ERROR(eslEINVAL, "symbol not in the alphabet");
  x = sp - a->sym;
  a->inmap[(int) sym] = x;
  return eslOK;
}

/* Function:  esl_alphabet_SetCaseInsensitive()
 * Incept:    SRE, Mon Dec 20 15:31:12 2004 [Zaragoza]
 *
 * Purpose:   Given an alphabet <a>, with all synonyms set,
 *            make the input map case-insensitive: for every
 *            letter that is mapped in either lower or upper
 *            case, map the other case to the same internal
 *            residue.
 *
 * Args:      a  - alphabet to make inmap case-insensitive.
 *                 
 * Returns:   <eslOK> on success.                
 */
int
esl_alphabet_SetCaseInsensitive(ESL_ALPHABET *a)
{
  int lc, uc;

  for (lc = 'a'; lc <= 'z'; lc++)
    {
      uc = lc; toupper(uc);
      if (a->inmap[lc] >= 0 && a->inmap[uc] < 0) a->inmap[uc] = a->inmap[lc];
      if (a->inmap[uc] >= 0 && a->inmap[lc] < 0) a->inmap[lc] = a->inmap[uc];
    }
  return eslOK;
}

/* Function:  esl_alphabet_SetDegeneracy()
 * Incept:    SRE, Mon Dec 20 15:42:23 2004 [Zaragoza]
 *
 * Purpose:   Given an alphabet under construction, 
 *            define the degenerate character <c> to mean
 *            any of the characters in the string <ds>.
 *
 * Args:      a   - an alphabet under construction.
 *            c   - degenerate character code; example: 'R'
 *            ds  - string of base characters for c; example: "AG"
 *
 * Returns:   <eslOK> on success.
 */
int
esl_alphabet_SetDegeneracy(ESL_ALPHABET *a, char c, char *ds)
{
  char *sp;
  int   x,y;

  if ((sp = strchr(a->sym, c)) == NULL)
    ESL_ERROR(eslEINVAL, "no such degenerate character");
  x = sp - a->sym;
  
  while (*ds != '\0') {
    if ((sp = strchr(a->sym, *ds)) == NULL) 
      ESL_ERROR(eslEINVAL, "no such base character");
    y = sp - a->sym;

    a->degen[x][y] = 1;
    a->ndegen[x]++;
    ds++;
  }

  return eslOK;
}


/* Function: esl_abc_CreateDigitalSequence()
 * 
 * Purpose:  Given an alphabet <a> and a sequence <seq> of length <L> 
 *           residues, allocate and create a digitized sequence
 *           and return it through <ret_dsq>. Caller must free the dsq.
 *           
 * Args:     a       - internal alphabet
 *           seq     - sequence to be digitized (0..L-1, alphabetic)
 *           L       - length of sequence      
 *           ret_dsq - RETURN: the new digital sequence
 *           
 * Returns:  <eslOK> on success, and digitized sequence is passed 
 *           back to the caller via <ret_dsq>; caller is responsible
 *           for freeing this with <free(dsq)>.
 *           Returns <eslEINVAL> if <seq> contains one or more characters
 *           that are not recognized in the alphabet <a>.
 *           
 * Throws:   <eslEMEM> if allocation fails.          
 */
int
esl_abc_CreateDigitalSequence(ESL_ALPHABET *a, char *seq, int L,char **ret_dsq)
{
  char *dsq;

  *ret_dsq = NULL;
  if ((dsq = malloc (sizeof(char) * (L+2))) == NULL) 
    ESL_ERROR(eslEMEM, "malloc failed");
  if (esl_abc_DigitizeSequence(a, seq, L, dsq) != eslOK)
    ESL_ERROR(eslEINVAL, "failed to digitize sequence");
  *ret_dsq = dsq;
  return eslOK;
}

/* Function:  esl_abc_DigitizeSequence()
 * Incept:    SRE, Mon Dec 20 16:40:31 2004 [Zaragoza]
 *
 * Purpose:   Given an allocated <dsq> of length <L> (that is, array of
 *            <L>+2 chars), reuse it, digitizing up to <L> characters of <seq>,
 *            according to the alphabet <a>. <seq> may be of any length, but
 *            <dsq> will not contain more than <L> characters of it.
 *            
 *            Usually, both <dsq> and <seq> have the same number of
 *            residues <L>, but the API allows <dsq> to be a window on
 *            a longer sequence.
 *            
 * Args:      a      - internal alphabet
 *            seq    - alphabetic input sequence
 *            L      - allocated length of dsq; max len of <seq> to digitize.
 *            dsq    - allocated space for digital sequence, 1..L
 *
 * Returns:   <eslOK> on success, and <dsq> contains newly digitized <seq>.
 *            <eslEINVAL> if any character of <seq> is not in the input map 
 *            of the alphabet <a>.
 */
int
esl_abc_DigitizeSequence(ESL_ALPHABET *a, char *seq, int L, char *dsq)
{
  int i;
  int x;

  dsq[0] = eslSENTINEL;
  for (i = 1; i <= L && seq[i-1] != '\0'; i++) 
    { 
      x = a->inmap[(int) seq[i-1]];
      if (x < 0)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "Symbol %c is not in the input alphabet", seq[i-1]);
	  return eslEINVAL;
	}
      else dsq[i] = x;
    }

  dsq[i] = eslSENTINEL;
  return eslOK;
}


/* Function:  esl_abc_AvgIScore()
 * Incept:    SRE, Tue Dec 21 10:53:57 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet
 *            <a>, and an array of integer scores <sc> for the residues
 *            in the base alphabet, calculate and return the 
 *            average score (rounded to nearest integer).
 */
int
esl_abc_AvgIScore(ESL_ALPHABET *a, char x, int *sc)
{
  float result = 0.;
  int i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += (float) sc[i];
  result /= (float) a->ndegen[(int) x];
  if (result < 0) return (int) (result - 0.5);
  else            return (int) (result + 0.5);
}

/* Function:  esl_abc_AvgFScore()
 * Incept:    SRE, Tue Dec 21 10:56:52 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet
 *            <a>, and an array of float scores <sc> for the residues
 *            in the base alphabet, calculate and return the 
 *            average score.
 */
float
esl_abc_AvgFScore(ESL_ALPHABET *a, char x, float *sc)
{
  float result = 0.;
  int   i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += sc[i];
  result /= (float) a->ndegen[(int) x];
  return result;
}
  
/* Function:  esl_abc_AvgDScore()
 * Incept:    SRE, Tue Dec 21 10:57:30 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet
 *            <a>, and an array of double scores <sc> for the residues
 *            in the base alphabet, calculate and return the 
 *            average score.
 */
double
esl_abc_AvgDScore(ESL_ALPHABET *a, char x, double *sc)
{
  double result = 0.;
  int    i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += sc[i];
  result /= (double) a->ndegen[(int) x];
  return result;
}


/* Function:  esl_abc_ExpectIScore()
 * Incept:    SRE, Tue Dec 21 11:02:46 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet <a>, an
 *            array of integer scores <sc> for the residues in the base
 *            alphabet, and background frequencies <p> for the
 *            occurrence frequencies of the residues in the base
 *            alphabet, calculate and return the expected score
 *            (weighted by the occurrence frequencies).
 */
int
esl_abc_ExpectIScore(ESL_ALPHABET *a, char x, int *sc, float *p)
{
  float  result = 0.;
  float  denom  = 0.;
  int    i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) { 
      result += (float) sc[i] * p[i];
      denom  += p[i];
    }
  result /= denom;
  if (result < 0) return (int) (result - 0.5);
  else            return (int) (result + 0.5);
}

/* Function:  esl_abc_ExpectFScore()
 * Incept:    SRE, Tue Dec 21 11:04:48 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet <a>, an
 *            array of float scores <sc> for the residues in the base
 *            alphabet, and background frequencies <p> for the
 *            occurrence frequencies of the residues in the base
 *            alphabet, calculate and return the expected score
 *            (weighted by the occurrence frequencies).
 */
float
esl_abc_ExpectFScore(ESL_ALPHABET *a, char x, float *sc, float *p)
{
  float  result = 0.;
  float  denom  = 0.;
  int    i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) { 
      result += sc[i] * p[i];
      denom  += p[i];
    }
  result /= denom;
  return result;
}

/* Function:  esl_abc_ExpectDScore()
 * Incept:    SRE, Tue Dec 21 11:06:10 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet <a>, an
 *            array of double scores <sc> for the residues in the base
 *            alphabet, and background frequencies <p> for the
 *            occurrence frequencies of the residues in the base
 *            alphabet, calculate and return the expected score
 *            (weighted by the occurrence frequencies).
 */
double
esl_abc_ExpectDScore(ESL_ALPHABET *a, char x, double *sc, double *p)
{
  double result = 0.;
  double denom  = 0.;
  int    i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) { 
      result += sc[i] * p[i];
      denom  += p[i];
    }
  result /= denom;
  return result;
}




/*****************************************************************
 * Test driver:
 * gcc -g -Wall -I. -o alphabet_test -DESL_ALPHABET_TESTDRIVE\
 *    alphabet.c easel.c 
 *****************************************************************
 */ 
#ifdef ESL_ALPHABET_TESTDRIVE

static void basic_examples(void);
static void degeneracy_integer_scores(void);
static void degeneracy_float_scores(void);
static void degeneracy_double_scores(void);

int
main(void)
{
  basic_examples();
  degeneracy_integer_scores();
  degeneracy_float_scores();
  degeneracy_double_scores();

  return eslOK;
}

static void
basic_examples(void)
{
  ESL_ALPHABET  *a1, *a2;
  char           dnaseq[] = "GARYTCN";
  char           aaseq[]  = "EFILQZU";
  int            L;
  char          *dsq, *dsq2;
  int            i;

  /* Example 1. 
   * Create a DNA alphabet; digitize a DNA sequence.
   */
  a1 = esl_alphabet_Create(eslDNA);
  L  = strlen(dnaseq);
  esl_abc_CreateDigitalSequence(a1, dnaseq, L, &dsq);
  esl_alphabet_Destroy(a1);

  /* Example 2. 
   * Create an RNA alphabet; digitize the same DNA sequence;
   * make sure it is equal to the dsq above (so T=U were
   * correctly synonymous on input).
   */
  a2 = esl_alphabet_Create(eslRNA);
  esl_abc_CreateDigitalSequence(a2, dnaseq, L, &dsq2);
  for (i = 1; i <= L; i++)
    if (dsq[i] != dsq2[i]) abort();
  esl_alphabet_Destroy(a2);

  /* Example 3.
   * Create an amino alphabet; digitize a protein sequence, 
   * while reusing memory already allocated in dsq.
   */
  a1 = esl_alphabet_Create(eslAMINO);
  esl_abc_DigitizeSequence(a1, aaseq, L, dsq);
  
  /* Example 4.
   * Create a custom alphabet almost the same as the amino
   * acid alphabet; digitize the same protein seq, reusing
   * memory in dsq2; check that seqs are identical.
   */
  a2 = esl_alphabet_CreateCustom("ACDEFGHIKLMNPQRSTVWY-BZX", 20, 24);
  esl_alphabet_SetSynonym(a2, 'U', 'S');     /* read selenocys U as serine S */
  esl_alphabet_SetCaseInsensitive(a2);       /* allow lower case input */
  esl_alphabet_SetDegeneracy(a2, 'Z', "QE");

  esl_abc_DigitizeSequence(a2, aaseq, L, dsq2);
  for (i = 1; i <= L; i++)
    if (dsq[i] != dsq2[i]) abort();

  /* clean up.
   */
  esl_alphabet_Destroy(a1);
  esl_alphabet_Destroy(a2);
  free(dsq);
  free(dsq2);
}


static void
degeneracy_integer_scores(void)
{
  ESL_ALPHABET *a;
  int           x;
  float         p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  int           sc[] = { -1,  -6,   6,   1};
  int           val;

  a     = esl_alphabet_Create(eslDNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_AvgIScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (val != 0) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_ExpectIScore(a, x, sc, p);  
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (val != -2) abort();

  esl_alphabet_Destroy(a);
  return;
}

static void
degeneracy_float_scores(void)
{
  ESL_ALPHABET *a;
  int           x;
  float         p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  float         sc[] = { -1., -6.,  6., 1.};
  float         val;

  a     = esl_alphabet_Create(eslRNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_AvgFScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (fabs(val - 0.) > 0.0001) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_ExpectFScore(a, x, sc, p);  
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (fabs(val + 2.) > 0.0001) abort();

  esl_alphabet_Destroy(a);
  return;
}

static void
degeneracy_double_scores(void)
{
  ESL_ALPHABET *a;
  int           x;
  double        p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  double        sc[] = { -1., -6.,  6., 1.};
  double        val;

  a     = esl_alphabet_Create(eslRNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_AvgDScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (fabs(val - 0.) > 0.0001) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_ExpectDScore(a, x, sc, p); 
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (fabs(val + 2.) > 0.0001) abort();

  esl_alphabet_Destroy(a);
  return;
}

#endif /*ESL_ALPHABET_TESTDRIVE*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
