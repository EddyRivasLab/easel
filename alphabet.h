/* alphabet.h
 * SRE, Tue Nov 23 19:44:01 2004 [St. Louis]
 * 
 * The biological symbol alphabets, and how they're standardly digitized
 * in Easel. 
 * 
 * An application would typically define one ESL_ALPHABET structure globally,
 * and initialize it once.
 * 
 * The alphabet can be ESL_AMINO, ESL_DNA, or ESL_RNA. 
 * Their symbol alphabets follow these conventions (and parts of the code will
 * assume that these conventions are obeyed):
 *    1. The first K syms (sym[0..K-1]) are the standard alphabet, in alphabetical order.
 *    2. The remaining symbols are the IUPAC degenerate symbols.
 *    3. The final symbol (indexed Kall-1) is the fully degenerate residue, N or X.
 * Thus:
 *   DNA:     "ACGTRYMKSWHBVDXN"           K=4  Kall=16
 *   RNA:     "ACGURYMKSWHBVDXN"           K=4  Kall=16
 *   Amino:   "ACDEFGHIKLMNPQRSTVWYBZUX"   K=20 Kall=24
 *
 *  The DNA and RNA alphabets follow published IUBMB recommendations
 *  ["Nomenclature for incompletely specified bases in nucleic acid
 *  sequences", Eur. J. Biochem. 150:1-5 (1985);
 *  http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html], with the addition
 *  of X as a synonym for N (acquiescing to the BLAST filter standard of 
 *  using X's to mask residues), and the use of U in RNA sequences in
 *  place of T.
 *  
 *  The one-letter code for amino acids follows section 3AA-21 of the
 *  IUPAC recommendations ("Nomenclature and symbolism for amino acids
 *  and peptides", Eur. J. Biochem. 138:9-37, 1985;
 *  http://www.chem.qmul.ac.uk/iupac/AminoAcid/]; augmented by U for
 *  selenocysteine, as recommended by the JCBN/NC-IUBMB Newsletter, 1999
 *  [http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html], though
 *  it is not really a "degenerate" residue. Since we must map it onto one
 *  of the 20-letter code, we map it onto serine (S).
 */
#define ESL_MAXABET 20          /* maximum normal alphabet size  (4 or 20)             */
#define ESL_MAXCODE 24	 	/* maximum degenerate (IUPAC) alphabet size (16 or 24) */

#define eslNOTSETYET 0
#define eslDNA       2		/* compatibility with squid's kRNA, HMMER's hmmNUCLEIC   */
#define eslAMINO     3		/* compatibility with squid's kAmino, HMMER's hmmAMINO   */
#define eslRNA       4

struct esl_alphabet_s {
  int   type;		          /* eslDNA, eslRNA, or eslAMINO                       */
  int   K;		          /* uniq alphabet size: 4 or 20                       */
  int   Kall;		          /* total size of alphabet + IUPAC degen; 16 or 24    */
  char  sym[MAXCODE+1];           /* "ACGTRYMKSWHBVDXN", for instance                  */ 
  int   symmap[128];              /* symmap['A'] = 0, etc: dsq[] index for a symbol    */
  char  degen[MAXCODE][MAXABET];  /* 1/0 arrays, for whether IUPAC code includes a residue */
  int   ndegen[MAXCODE];
};
typedef struct esl_alphabet_s ESL_ALPHABET;


extern void           SetAlphabet(int type);
extern unsigned char  SymbolIndex(char sym);
extern unsigned char *DigitizeSequence(char *seq, int L);
