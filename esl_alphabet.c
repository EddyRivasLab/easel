/* esl_alphabet.c
 * Implements the standard digitized alphabets for biosequences.
 * 
 *    1. ESL_ALPHABET object  - a digital alphabet
 *    2. Digitized sequences (ESL_DSQ *)
 *    3. Other routines in the API
 *    4. Example code
 *    5. Test driver
 * 
 * SVN $Id$
 * SRE, Tue Dec  7 13:49:43 2004
 */
#include <esl_config.h>

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>


/*****************************************************************
 * 1. The ESL_ALPHABET object
 *****************************************************************/ 

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
  int           status;
  ESL_ALPHABET *a;

  switch(type) { 
  case eslAMINO:  a = create_amino(); break;
  case eslDNA:    a = create_dna();   break;
  case eslRNA:    a = create_rna();   break;
  default:    
    ESL_FAIL(eslEINVAL, "Standard alphabets include only DNA, RNA, protein.");
  }
  return a;

 FAILURE:
  return NULL;
}


/* Function:  esl_alphabet_CreateCustom()
 * Incept:    SRE, Mon Dec 20 09:18:28 2004 [Zaragoza]
 *
 * Purpose:   Creates a customized biosequence alphabet,
 *            and returns a ptr to it. The alphabet type is set 
 *            to <eslNONSTANDARD>.
 *            
 *            <alphabet> is the internal alphabet string;
 *            <K> is the size of the base alphabet;
 *            <Kp> is the total size of the alphabet string. 
 *            
 *            In the alphabet string, residues <0..K-1> are the base alphabet; 
 *            residue <K> is the canonical gap (indel) symbol; 
 *            residues <K+1..Kp-3> are additional degeneracy symbols (possibly 0 of them);
 *            residue <Kp-2> is an "any" symbol (such as N or X); 
 *            and residue <Kp-1> is a "missing data" gap symbol.
 *            
 *            The two gap symbols and the "any" symbol are mandatory even for
 *            nonstandard alphabets, so $Kp \geq K+3$.
 *            
 * Args:      alphabet - internal alphabet; example "ACGT-RYMKSWHBVDN~"
 *            K        - base size; example 4
 *            Kp       - total size, including gap, degeneracies; example 17
 *
 * Returns:   pointer to new <ESL_ALPHABET> structure
 *
 * Throws:    <NULL> if any allocation or initialization fails.
 */
ESL_ALPHABET *
esl_alphabet_CreateCustom(char *alphabet, int K, int Kp)
{
  ESL_ALPHABET *a;
  int           c,x,y;
  int           status;

  /* Argument checks.
   */
  if (strlen(alphabet) != Kp) ESL_FAIL(eslEINVAL, "alphabet length != Kp");
  if (Kp < K+3)               ESL_FAIL(eslEINVAL, "Kp too small in alphabet"); 

  /* Allocation/init, level 1.
   */
  ESL_ALLOC(a, sizeof(ESL_ALPHABET));
  a->sym    = NULL;
  a->degen  = NULL;
  a->ndegen = NULL;
  
  /* Allocation/init, level 2.
   */
  ESL_ALLOC(a->sym,    sizeof(char)   * (Kp+1));
  ESL_ALLOC(a->degen,  sizeof(char *) * Kp);
  ESL_ALLOC(a->ndegen, sizeof(int)    * Kp);
  a->degen[0] = NULL;

  /* Allocation/init, level 3.
   */
  ESL_ALLOC(a->degen[0], sizeof(char) * (Kp*K));
  for (x = 1; x < Kp; x++)
    a->degen[x] = a->degen[0]+(K*x);

  /* Initialize the internal alphabet: 
   */
  a->type = eslNONSTANDARD;
  a->K    = K;
  a->Kp   = Kp;
  strcpy(a->sym, alphabet);

  /* Initialize the input map, mapping ASCII seq chars to digital codes,
   * and eslDSQ_ILLEGAL for everything else.
   */
  for (c = 0; c < 128; c++)   a->inmap[c]               = eslDSQ_ILLEGAL;
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
  a->ndegen[Kp-2]  = K;
  for (x = 0; x < a->K; x++) a->degen[Kp-2][x] = 1;

  return a;

 FAILURE:
  esl_alphabet_Destroy(a);
  return NULL;
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
  if ((a = esl_alphabet_CreateCustom("ACGT-RYMKSWHBVDN~", 4, 17)) == NULL) return NULL;
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
  if ((a = esl_alphabet_CreateCustom("ACGU-RYMKSWHBVDN~", 4, 17)) == NULL) return NULL;
  a->type = eslRNA;
  
  /* Add desired synonyms in the input map.
   */
  esl_alphabet_SetSynonym(a, 'T', 'U');	    /* read T as a U */
  esl_alphabet_SetSynonym(a, 'X', 'N');	    /* read X as an N (many seq maskers use X) */
  esl_alphabet_SetSynonym(a, '_', '-');     /* allow _ as a gap too */
  esl_alphabet_SetSynonym(a, '.', '-');     /* allow . as a gap too */
  esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input */

  /* Define degenerate symbols.
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

  /* Create the internal alphabet
   */
  if ((a = esl_alphabet_CreateCustom("ACDEFGHIKLMNPQRSTVWY-BJZOUX~", 20, 28)) == NULL) return NULL;
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
  esl_alphabet_SetDegeneracy(a, 'J', "IL");
  esl_alphabet_SetDegeneracy(a, 'Z', "QE");

  /* Define unusual residues as one-to-one degeneracies.
   */
  esl_alphabet_SetDegeneracy(a, 'U', "C"); /* selenocysteine is scored as cysteine */
  esl_alphabet_SetDegeneracy(a, 'O', "K"); /* pyrrolysine is scored as lysine      */

  return a;
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
  char    *sp;
  ESL_DSQ  x;

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
 * 
 * Throws:    <eslECORRUPT> if any lower/uppercase symbol pairs
 *            are already both mapped to different symbols.
 */
int
esl_alphabet_SetCaseInsensitive(ESL_ALPHABET *a)
{
  int lc, uc;

  for (lc = 'a'; lc <= 'z'; lc++)
    {
      uc = toupper(lc);

      if      (esl_abc_CIsValid(a, lc) && ! esl_abc_CIsValid(a, uc)) a->inmap[uc] = a->inmap[lc];
      else if (esl_abc_CIsValid(a, uc) && ! esl_abc_CIsValid(a, lc)) a->inmap[lc] = a->inmap[uc];
      else if (esl_abc_CIsValid(a, lc) && esl_abc_CIsValid(a, uc) && a->inmap[uc] != a->inmap[lc])
	ESL_ERROR(eslECORRUPT, "symbols %c and %c map differently already (%c vs. %c)",
		  lc, uc, a->inmap[lc], a->inmap[uc]);
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
  char   *sp;
  ESL_DSQ x,y;

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



/*****************************************************************
 * 2. Digitized sequences: char *dsq[1..L], with sentinels at 0,L+1
 *****************************************************************/ 

/* Function: esl_abc_Digitize()
 * Incept:   SRE, Sun Aug 27 11:18:56 2006 [Leesburg]
 * 
 * Purpose:  Given an alphabet <a> and an ASCII sequence <seq> of
 *           length <L> residues, digitize the sequence and put it
 *           in <dsq>. Caller provides space in <dsq> allocated for
 *           at least <L+2> <ESL_DSQ> residues.
 *           
 * Args:     a       - internal alphabet
 *           seq     - sequence to be digitized (0..L-1, alphabetic);
 *                     does not need to be NUL-terminated.
 *           L       - length of sequence      
 *           dsq     - RETURN: the new digital sequence (caller allocates,
 *                     at least <(L+2) * sizeof(ESL_DSQ)>).
 *           
 * Returns:  <eslOK> on success.
 *           Returns <eslEINVAL> if <seq> contains one or more characters
 *           that are not recognized in the alphabet <a>. (This is classed
 *           as a normal error, because the <seq> may be user input.)
 *           If this happens, the digital sequence <dsq> is still valid upon
 *           return; invalid ASCII characters are replaced by ambiguities
 *           (X or N).
 */
int
esl_abc_Digitize(ESL_ALPHABET *a, char *seq, int L, ESL_DSQ *dsq)
{
  int     status;
  int     i;
  ESL_DSQ x;

  status = eslOK;
  dsq[0] = eslDSQ_SENTINEL;
  for (i = 1; i <= L && seq[i-1] != '\0'; i++) 
    { 
      x = a->inmap[(int) seq[i-1]];
      if (esl_abc_XIsValid(a, x))
	dsq[i] = x;
      else
	{
	  status = eslEINVAL;
	  dsq[i] = esl_abc_XGetUnknown(a);
	}
    }
  dsq[L+1] = eslDSQ_SENTINEL;
  return status;
}

/* Function:  esl_abc_Textize()
 * Incept:    SRE, Sun Aug 27 11:14:58 2006 [Leesburg]
 *
 * Purpose:   Make an ASCII sequence <seq> by converting a digital
 *            sequence <dsq> of length <L> back to text, according to
 *            the digital alphabet <a>. 
 *            
 *            Caller provides space in <seq> allocated for at least
 *            <L+1> bytes (<(L+1) * sizeof(char)>).
 *
 * Args:      a   - internal alphabet
 *            dsq - digital sequence to be converted (1..L)
 *            L   - length of dsq
 *            seq - RETURN: the new text sequence (caller allocated
 *                  space, at least <(L+1) * sizeof(char)>).
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslECORRUPT> if something's wrong with the <dsq>, like
 *            a premature sentinel byte or an invalid residue. 
 */
int
esl_abc_Textize(ESL_ALPHABET *a, ESL_DSQ *dsq, int L, char *seq)
{
  int i;
  
  for (i = 0; i < L; i++)
    {
      if (! esl_abc_XIsValid(a, dsq[i+1]))
	ESL_ERROR(eslECORRUPT, "bad code in dsq");
      seq[i] = a->sym[dsq[i+1]];
    }
  seq[i] = '\0';
  return eslOK;
}


/* Function:  esl_abc_TextizeN()
 * Incept:    SRE, Tue Sep  5 09:28:38 2006 [Janelia] STL11/54.
 *
 * Purpose:   Similar in semantics to <strncpy()>, this procedure takes
 *            a window of <L> residues in a digitized sequence
 *            starting at the residue pointed to by <dptr>,
 *            converts them to ASCII text representation, and 
 *            copies them into the buffer <buf>.
 *            
 *            <buf> must be at least <L> residues long; <L+1>, if the
 *            caller needs to NUL-terminate it.
 *            
 *            If a sentinel byte is encountered in the digitized
 *            sequence before <L> residues have been copied, <buf> is
 *            NUL-terminated there. Otherwise, like <strncpy()>, <buf>
 *            will not be NUL-terminated.
 *            
 *            Note that because digital sequences are indexed <1..N>,
 *            not <0..N-1>, the caller must be careful about
 *            off-by-one errors in <dptr>. For example, to copy from
 *            the first residue of a digital sequence <dsq>, you must
 *            pass <dptr=dsq+1>, not <dptr=dsq>. The text in <buf>
 *            on the other hand is a normal C string indexed <0..L-1>.
 *
 * Args:      a     - reference to an internal alphabet
 *            dptr  - ptr to starting residue in a digital sequence
 *            L     - number of residues to convert and copy
 *            buf   - text buffer to store the <L> converted residues in
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslECORRUPT> if something's wrong with <dsq>, like a bogus
 *            symbol.
 */
int
esl_abc_TextizeN(ESL_ALPHABET *a, ESL_DSQ *dptr, int L, char *buf)
{
  int i;

  for (i = 0; i < L; i++)
    {
      if (dptr[i] == eslDSQ_SENTINEL) 
	{ 
	  buf[i] = '\0';
	  return eslOK;
	}

      if (! esl_abc_XIsValid(a, dptr[i]))
	ESL_ERROR(eslECORRUPT, "bad code in dsq");

      buf[i] = a->sym[dptr[i]];
    }
  return eslOK;
}


/* Function:  esl_abc_dsqdup()
 * Incept:    SRE, Tue Aug 29 13:51:05 2006 [Janelia]
 *
 * Purpose:   Like <esl_strdup()>, but for digitized sequences:
 *            make a duplicate of <dsq> and leave it in <ret_dup>.
 *            Caller can pass the string length <L> if it's known, saving
 *            some overhead; else pass <-1> and the length will be
 *            determined for you.
 *            
 *            Tolerates <dsq> being <NULL>; in which case, returns
 *            <eslOK> with <*ret_dup> set to <NULL>.
 *
 * Args:      dsq     - digital sequence to duplicate (w/ sentinels at 0,L+1)
 *            L       - length of dsq in residues, if known; -1 if unknown
 *            ret_dup - RETURN: allocated duplicate of <dsq>, which caller will
 *                      free.
 *
 * Returns:   <eslOK> on success, and leaves a pointer in <ret_dup>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      STL11/48
 */
int 
esl_abc_dsqdup(ESL_DSQ *dsq, int L, ESL_DSQ **ret_dup)
{
  int      status;
  ESL_DSQ *new = NULL;

  *ret_dup = NULL;
  if (dsq == NULL) return eslOK;
  if (L < 0) L = esl_abc_dsqlen(dsq);

  ESL_ALLOC(new, sizeof(ESL_DSQ) * (L+2));
  memcpy(new, dsq, sizeof(ESL_DSQ) * (L+2));
  
  *ret_dup = new;
  return eslOK;

 FAILURE:
  if (new     != NULL) free(new);
  *ret_dup = NULL;
  return status;
}


/* Function:  esl_abc_dsqcat()
 * Incept:    SRE, Tue Aug 29 14:01:59 2006 [Janelia]
 *
 * Purpose: Like <esl_strcat()>, except specialized for digitizing a
 *            biosequence text string and appending it to a growing
 *            digital sequence. The growing digital sequence is <dsq>,
 *            currently of length <L> residues; we append <s> to it,
 *            of length <n> symbols, after digitization.  Upon return,
 *            <dsq> has been reallocated and <L> is set to the new
 *            length (which is why both must be passed by reference).
 *            
 *            Note that the final <L> is not necessarily the initial
 *            <L> plus <n>, because the text string <s> may contain
 *            symbols that are defined to be ignored
 *            (<eslDSQ_IGNORED>) in the input map of this alphabet.
 *            (The final <L> is guaranteed to be $\leq$ <L+n> though.>
 *            
 *            If the initial <L> is unknown, pass <-1>, and it will be
 *            determined by counting the residues in <dsq>.
 *            
 *            Similarly, if <n> is unknown, pass <-1> and it will be
 *            determined by counting the symbols in <s>
 *            
 *            <dsq> may be <NULL>, in which case this call is
 *            equivalent to an allocation and digitization just of
 *            <s>.
 *            
 *            <s> may also be <NULL>, in which case <dsq> is
 *            unmodified; <L> would be set to the correct length of
 *            <dsq> if it was passed as <-1> (unknown).
 *            
 * Args:      abc  - digital alphabet to use
 *            dsq  - reference to the current digital seq to append to 
 *                   (with sentinel bytes at 0,L+1); may be <NULL>. 
 *                   Upon return, this will probably have 
 *                   been reallocated, and it will contain the original
 *                   <dsq> with <s> digitized and appended.
 *            L    - reference to the current length of <dsq> in residues;
 *                   may be <-1> if unknown. Upon return, <L> is set to
 *                   the new length of <dsq>, after <s> is appended.
 *            s    - NUL-terminated ASCII text sequence to append. May
 *                   contain ignored text characters (flagged with
 *                   <eslDSQ_IGNORED> in the input map of alphabet <abc>).  
 *            n    - Length of <s> in characters, if known; or <-1> if 
 *                   unknown.
 *
 * Returns:   <eslOK> on success; <dsq> contains the result of digitizing
 *            and appending <s> to the original <dsq>; and <L> contains
 *            the new length of the <dsq> result in residues.
 *            
 *            If any of the characters in <s> are illegal in the alphabet
 *            <abc>, these characters are digitized as unknown residues, 
 *            and the function returns <eslEINVAL>. The caller might want
 *            to call <esl_abc_ValidateSeq()> on <s> if it wants to figure
 *            out where digitization goes awry and get a more informative
 *            error report. This is a normal error, because the string <s>
 *            might be user input.
 *
 * Throws:    <eslEMEM> on allocation or reallocation failure;
 *
 * Xref:      STL11/48.
 */
int
esl_abc_dsqcat(ESL_ALPHABET *a, ESL_DSQ **dsq, int *L, char *s, int n)
{
  int     status;
  void   *p;
  int     newL;
  int     xpos, cpos;
  ESL_DSQ x;

  if (*L < 0) newL = ((*dsq == NULL) ? 0 : esl_abc_dsqlen(*dsq));
  else        newL = *L;

  if (n < 0)  n = ((s == NULL) ? 0 : strlen(s));

  /* below handles weird case of empty s (including empty dsq and empty s):
   * just hand dsq and its length right back to the caller.
   */
  if (n == 0) { *L = newL; return eslOK; } 

  if (*dsq == NULL) {		/* an entirely new dsq must be allocated *and* initialized with left sentinel. */
    ESL_ALLOC(*dsq, sizeof(ESL_DSQ)     * (n+2));
    (*dsq)[0] = eslDSQ_SENTINEL;
  } else			/* else, existing dsq is just reallocated; left sentinel already in place. */
    ESL_RALLOC(*dsq, p, sizeof(ESL_DSQ) * (newL+n+2)); /* most we'll need */

  /* Watch these coords. Start in the 0..n-1 text string at 0;
   * start in the 1..L dsq at L+1, overwriting its terminal 
   * sentinel byte.
   */
  status = eslOK;
  for (xpos = newL+1, cpos = 0; s[cpos] != '\0'; cpos++)
    {
      x = a->inmap[(int) s[cpos]];
      if (esl_abc_XIsValid(a, x))
	(*dsq)[xpos++] = x;
      else if (x == eslDSQ_IGNORED)
	;
      else 
	{
	  (*dsq)[xpos++] = esl_abc_XGetUnknown(a);
	  status = eslEINVAL;
	}
    }
  (*dsq)[xpos] = eslDSQ_SENTINEL;
  *L = xpos-1;
  return status;

 FAILURE:
  *L = newL;
  return status;
}


/* Function:  esl_abc_dsqlen()
 * Incept:    SRE, Tue Aug 29 13:49:02 2006 [Janelia]
 *
 * Purpose:   Returns the length of digitized sequence <dsq>
 *            in residues. The <dsq> must be properly terminated
 *            by a sentinel byte (<eslDSQ_SENTINEL>). 
 */
int 
esl_abc_dsqlen(ESL_DSQ *dsq)
{
  int n = 0;
  while (dsq[n+1] != eslDSQ_SENTINEL) n++;
  return n;
}


/*****************************************************************
 * 3. Other routines
 *****************************************************************/ 

/* Function:  esl_abc_IAvgScore()
 * Incept:    SRE, Tue Dec 21 10:53:57 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet
 *            <a>, and an array of integer scores <sc> for the residues
 *            in the base alphabet, calculate and return the 
 *            average score (rounded to nearest integer).
 *            
 *            <esl_abc_FAvgScore()> and <esl_abc_DAvgScore()> do the
 *            same, but for float and double scores instead of integers
 *            (for real-valued scores, no rounding is done).
 */
int
esl_abc_IAvgScore(ESL_ALPHABET *a, ESL_DSQ x, int *sc)
{
  float result = 0.;
  int i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += (float) sc[i];
  result /= (float) a->ndegen[(int) x];
  if (result < 0) return (int) (result - 0.5);
  else            return (int) (result + 0.5);
}
float
esl_abc_FAvgScore(ESL_ALPHABET *a, ESL_DSQ x, float *sc)
{
  float result = 0.;
  int   i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += sc[i];
  result /= (float) a->ndegen[(int) x];
  return result;
}
double
esl_abc_DAvgScore(ESL_ALPHABET *a, ESL_DSQ x, double *sc)
{
  double result = 0.;
  int    i;

  for (i = 0; i < a->K; i++)
    if (a->degen[(int) x][i]) result += sc[i];
  result /= (double) a->ndegen[(int) x];
  return result;
}


/* Function:  esl_abc_IExpectScore()
 * Incept:    SRE, Tue Dec 21 11:02:46 2004 [Zaragoza]
 *
 * Purpose:   Given a (degenerate) residue code <x> in alphabet <a>, an
 *            array of integer scores <sc> for the residues in the base
 *            alphabet, and background frequencies <p> for the
 *            occurrence frequencies of the residues in the base
 *            alphabet, calculate and return the expected score
 *            (weighted by the occurrence frequencies).
 *            
 *            <esl_abc_FExpectScore()> and <esl_abc_DExpectScore()> do the
 *            same, but for float and double scores instead of integers
 *            (for real-valued scores, no rounding is done).
 */
int
esl_abc_IExpectScore(ESL_ALPHABET *a, ESL_DSQ x, int *sc, float *p)
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
float
esl_abc_FExpectScore(ESL_ALPHABET *a, ESL_DSQ x, float *sc, float *p)
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
double
esl_abc_DExpectScore(ESL_ALPHABET *a, ESL_DSQ x, double *sc, double *p)
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

/* Function:  esl_abc_IAvgScVec()
 * Incept:    SRE, Thu Apr  6 12:12:25 2006 [AA890 enroute to Boston]
 *
 * Purpose:   Given an alphabet <a> and a score vector <sc> of length
 *            <a->Kp> that contains integer scores for the base
 *            alphabet (<0..a->K-1>), fill out the rest of the score 
 *            vector, calculating average scores for 
 *            degenerate residues using <esl_abc_IAvgScore()>.
 *            
 *            The score, if any, for a gap character is not touched by
 *            this function.
 *            
 *            <esl_abc_FAvgScVec()> and <esl_abc_DAvgScVec()> do
 *            the same, but for score vectors of floats or doubles,
 *            respectively.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_abc_IAvgScVec(ESL_ALPHABET *a, int *sc)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_IAvgScore(a, x, sc);
  return eslOK;
}
int
esl_abc_FAvgScVec(ESL_ALPHABET *a, float *sc)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_FAvgScore(a, x, sc);
  return eslOK;
}
int
esl_abc_DAvgScVec(ESL_ALPHABET *a, double *sc)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_DAvgScore(a, x, sc);
  return eslOK;
}

/* Function:  esl_abc_IExpectScVec()
 * Incept:    SRE, Thu Apr  6 12:23:52 2006 [AA 890 enroute to Boston]
 *
 * Purpose:   Given an alphabet <a>, a score vector <sc> of length
 *            <a->Kp> that contains integer scores for the base
 *            alphabet (<0..a->K-1>), and residue occurrence probabilities
 *            <p[0..a->K-1]>; fill out the rest of the score 
 *            vector, calculating expected scores for 
 *            degenerate residues using <esl_abc_IExpectScore()>.
 *            
 *            The score, if any, for a gap character is not touched by
 *            this function.
 *            
 *            <esl_abc_FExpectScVec()> and <esl_abc_DExpectScVec()> do
 *            the same, but for score vectors of floats or doubles,
 *            respectively. The probabilities <p> are floats for the
 *            integer and float versions, and doubles for the double
 *            version.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_abc_IExpectScVec(ESL_ALPHABET *a, int *sc, float *p)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_IExpectScore(a, x, sc, p);
  return eslOK;
}
int
esl_abc_FExpectScVec(ESL_ALPHABET *a, float *sc, float *p)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_FExpectScore(a, x, sc, p);
  return eslOK;
}
int
esl_abc_DExpectScVec(ESL_ALPHABET *a, double *sc, double *p)
{
  ESL_DSQ x;
  for (x = a->Kp+1; x <= a->Kp; x++)
    sc[x] = esl_abc_DExpectScore(a, x, sc, p);
  return eslOK;
}


/* Function:  esl_abc_FCount()
 * Incept:    SRE, Wed Apr 12 17:16:35 2006 [St. Louis]
 *
 * Purpose:   Count a possibly degenerate digital symbol <x> (0..Kp-1)
 *            into a counts array <ct> for base symbols (0..K-1).
 *            Assign the symbol a weight of <wt> (often just 1.0).
 *            The count weight of a degenerate symbol is divided equally
 *            across the possible base symbols. 
 *            
 *            <x> can be a gap; if so, <ct> must be allocated 0..K,
 *            not 0..K-1.
 *            
 *            <esl_abc_DCount()> does the same, but for double-precision
 *            count vectors and weights.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_abc_FCount(ESL_ALPHABET *abc, float *ct, ESL_DSQ x, float wt)
{
  ESL_DSQ y;

  if (x < abc->K) 
    ct[x] += wt;
  else
    for (y = 0; y < abc->K; y++) {
      if (abc->degen[x][y])
	ct[x] += wt / (float) abc->ndegen[x];
    }
  return eslOK;
}
int
esl_abc_DCount(ESL_ALPHABET *abc, double *ct, ESL_DSQ x, double wt)
{
  ESL_DSQ y;

  if (x < abc->K) 
    ct[x] += wt;
  else
    for (y = 0; y < abc->K; y++) {
      if (abc->degen[x][y])
	ct[x] += wt / (double) abc->ndegen[x];
    }
  return eslOK;
}

/* Function:  esl_abc_Type()
 * Incept:    SRE, Wed Apr 12 12:23:24 2006 [St. Louis]
 *
 * Purpose:   For diagnostics and other output: given an internal
 *            alphabet code <type> (<eslRNA>, for example), return
 *            ptr to an internal string ("RNA", for example). 
 */
char *
esl_abc_Type(int type)
{
  switch (type) {
  case eslUNKNOWN:     return "unknown";
  case eslRNA:         return "RNA";
  case eslDNA:         return "DNA";
  case eslAMINO:       return "protein";
  case eslNONSTANDARD: return "nonstandard/custom";
  default:             return "BOGUS";
  }
}


/* Function:  esl_abc_ValidateSeq()
 * Incept:    SRE, Sat Aug 26 17:40:00 2006 [St. Louis]
 *
 * Purpose:   Check that sequence <seq> of length <L> can be digitized
 *            without error; all its symbols are valid in alphabet
 *            <a>. If so, return <eslOK>. If not, return <eslEINVAL>.
 *            
 *            <errbuf> is either passed as <NULL>, or a pointer to an
 *            error string buffer allocated by the caller for
 *            <eslERRBUFSIZE> characters. If <errbuf> is non-NULL, and
 *            the sequence is invalid, an error message is placed in
 *            <errbuf>.
 *
 * Args:      a      - digital alphabet
 *            seq    - sequence to validate [0..L-1]; NUL-termination unnecessary
 *            L      - length of <seq>
 *            errbuf - NULL, or ptr to <eslERRBUFSIZE> chars of allocated space 
 *                     for an error message.
 *
 * Returns:   <eslOK> if <seq> is valid; <eslEINVAL> if not.
 *
 * Throws:    (no abnormal error conditions).
 */
int
esl_abc_ValidateSeq(ESL_ALPHABET *a, char *seq, int L, char *errbuf)
{
  int i;
  int firstpos = -1;
  int nbad     = 0;

  for (i = 0; i < L; i++)
    {
      if (! esl_abc_CIsValid(a, seq[i]))
	{
	  if (firstpos == -1) firstpos = i;
	  nbad++;
	}
    }

  if (nbad == 0)  return eslOK;

  if (errbuf != NULL) 
    sprintf(errbuf, "%d bad chars (including bad %c at pos %d)", 
	    nbad, seq[firstpos], firstpos);
  return eslEINVAL;
}




/*****************************************************************
 * 4. Examples
 *****************************************************************/ 

/*   gcc -g -Wall -o example -I. -DeslALPHABET_EXAMPLE esl_alphabet.c easel.c
 */
#ifdef eslALPHABET_EXAMPLE
/*::cexcerpt::alphabet_example::begin::*/
#include <easel.h>
#include <esl_alphabet.h>
int main(void)
{
  ESL_ALPHABET *a;
  char          dnaseq[] = "GARYTC";
  int           L        = 6;
  ESL_DSQ      *dsq;
  
  a = esl_alphabet_Create(eslDNA);

  if ((dsq = malloc(sizeof(ESL_DSQ * (L+2)))) == NULL)
    esl_fatal("malloc failed");
    
  if (esl_abc_Digitize(a, dnaseq, L, dsq) != eslOK) 
    esl_fatal("failed to digitize the sequence");

  free(dsq);
  esl_alphabet_Destroy(a);
  return 0;
}
/*::cexcerpt::alphabet_example::end::*/
#endif /*eslALPHABET_EXAMPLE*/


/*   gcc -g -Wall -o example -I. -DeslALPHABET_EXAMPLE2 esl_alphabet.c easel.c
 */
#ifdef eslALPHABET_EXAMPLE2
/*::cexcerpt::alphabet_example2::begin::*/
#include <easel.h>
#include <esl_alphabet.h>
int main(void)
{ 
  ESL_ALPHABET *a;

  /* 1. Create the base alphabet structure. */
  a = esl_alphabet_CreateCustom("ACDEFGHIKLMNOPQRSTUVWY-BJZX~", 22, 28);

  /* 2. Set your equivalences in the input map.  */
  esl_alphabet_SetEquiv(a, '.', '-');     /* allow . as a gap character too */

  /* 3. After all synonyms are set, (optionally) make map case-insensitive. */
  esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input too */

  /* 4. Define your optional degeneracy codes in the alphabet, one at a time.
   *    The 'any' character X was automatically set up.  */
  esl_alphabet_SetDegeneracy(a, 'B', "DN"); /* read B as {D|N} */
  esl_alphabet_SetDegeneracy(a, 'J', "IL"); /* read B as {I|L} */
  esl_alphabet_SetDegeneracy(a, 'Z', "QE"); /* read Z as {Q|E} */

  /* 5. (do your stuff) */

  /* 6. Remember to free it when you're done with it. */
  esl_alphabet_Destroy(a);
  return 0;
}
/*::cexcerpt::alphabet_example2::end::*/
#endif /*eslALPHABET_EXAMPLE2*/


/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

/* gcc -g -Wall -I. -L. -o test -DeslALPHABET_TESTDRIVE esl_alphabet.c easel.c -lm
 */
#ifdef eslALPHABET_TESTDRIVE
static void basic_examples(void);
static void degeneracy_integer_scores(void);
static void degeneracy_float_scores(void);
static void degeneracy_double_scores(void);

#include <easel.h>
#include <esl_alphabet.h>

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
  char           dnaseq[] = "GARYtcN";
  char           aaseq[]  = "EFILqzU";
  int            L;
  ESL_DSQ       *dsq, *dsq2;
  int            i;

  /* Example 1. 
   * Create a DNA alphabet; digitize a DNA sequence.
   */
  if ((a1 = esl_alphabet_Create(eslDNA)) == NULL)      abort();
  L  = strlen(dnaseq);
  if ((dsq = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL) abort();
  if (esl_abc_Digitize(a1, dnaseq, L, dsq) != eslOK)   abort();
  if (esl_abc_dsqlen(dsq) != L)                        abort();
  esl_alphabet_Destroy(a1);

  /* Example 2. 
   * Create an RNA alphabet; digitize the same DNA sequence;
   * make sure it is equal to the dsq above (so T=U were
   * correctly synonymous on input).
   */
  if ((a2 = esl_alphabet_Create(eslRNA)) == NULL)       abort();
  if ((dsq2 = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL) abort();
  if (esl_abc_Digitize(a2, dnaseq, L, dsq2) != eslOK)   abort();
  for (i = 1; i <= L; i++)
    if (dsq[i] != dsq2[i]) abort();
  esl_alphabet_Destroy(a2);

  /* Example 3.
   * Create an amino alphabet; digitize a protein sequence, 
   * while reusing memory already allocated in dsq.
   */
  if ((a1 = esl_alphabet_Create(eslAMINO)) == NULL)     abort();
  if (esl_abc_Digitize(a1, aaseq, L, dsq) != eslOK)     abort();
  
  /* Example 4.
   * Create a custom alphabet almost the same as the amino
   * acid alphabet; digitize the same protein seq, reusing
   * memory in dsq2; check that seqs are identical.
   */
  if ((a2 = esl_alphabet_CreateCustom("ACDEFGHIKLMNPQRSTVWY-BJZOUX~", 20, 28)) == NULL) abort();
  if (esl_alphabet_SetSynonym(a2, 'U', 'S') != eslOK)     abort();  /* read selenocys U as serine S */
  if (esl_alphabet_SetCaseInsensitive(a2)   != eslOK)     abort();  /* allow lower case input */
  if (esl_alphabet_SetDegeneracy(a2, 'Z', "QE") != eslOK) abort();

  if (esl_abc_Digitize(a2, aaseq, L, dsq2) != eslOK)      abort();
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
  ESL_DSQ       x;
  float         p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  int           sc[] = { -1,  -6,   6,   1};
  int           val;

  a     = esl_alphabet_Create(eslDNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_IAvgScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (val != 0) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_IExpectScore(a, x, sc, p);  
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (val != -2) abort();

  esl_alphabet_Destroy(a);
  return;
}

static void
degeneracy_float_scores(void)
{
  ESL_ALPHABET *a;
  ESL_DSQ       x;
  float         p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  float         sc[] = { -1., -6.,  6., 1.};
  float         val;

  a     = esl_alphabet_Create(eslRNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_FAvgScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (fabs(val - 0.) > 0.0001) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_FExpectScore(a, x, sc, p);  
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (fabs(val + 2.) > 0.0001) abort();

  esl_alphabet_Destroy(a);
  return;
}

static void
degeneracy_double_scores(void)
{
  ESL_ALPHABET *a;
  ESL_DSQ       x;
  double        p[]  = {0.4, 0.1, 0.1, 0.4}; /* A/T biased background */
  double        sc[] = { -1., -6.,  6., 1.};
  double        val;

  a     = esl_alphabet_Create(eslRNA);  

  x     = esl_abc_DigitizeSymbol(a, 'N'); /* any: A/C/G/T */
  val   = esl_abc_DAvgScore(a, x, sc); 
  /* average of -1,-6,6,1 = 0 */
  if (fabs(val - 0.) > 0.0001) abort();

  x     = esl_abc_DigitizeSymbol(a, 'M');     /* M = A/C */
  val   = esl_abc_DExpectScore(a, x, sc, p); 
  /* expectation of -1,-6 given p = 0.4,0.1 = -2 */
  if (fabs(val + 2.) > 0.0001) abort();

  esl_alphabet_Destroy(a);
  return;
}

#endif /*eslALPHABET_TESTDRIVE*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/

