/* alphabet.h [abc]
 * Digital representation of biosequence symbols in Easel.
 * SRE, Tue Nov 23 19:44:01 2004 [St. Louis]
 */
#ifndef ESL_ALPHABET_INCLUDED
#define ESL_ALPHABET_INCLUDED

/* Flags for alphabet types.
 */
#define eslUNKNOWN     0 /* 0=unknown is easel-wide convention; don't change */
#define eslRNA         1
#define eslDNA         2		
#define eslAMINO       3		
#define eslNONSTANDARD 4

/* Value of the sentinel bytes (0,L+1) in a dsq.
 */
#define eslSENTINEL    127

struct esl_alphabet_s {
  int    type;		 /* eslDNA, eslRNA, esl_AMINO, or eslNONSTANDARD    */
  int    K;		 /* uniq alphabet size: 4 or 20                     */
  int    Kp;		 /* total size of alphabet + gap + IUPAC degen      */
  char  *sym;            /* "ACGT-RYMKSWHBVDN", for instance    [0..Kp-1]   */ 
  int    inmap[128];     /* inmap['A'] = 0, etc: dsq[] index for a symbol   */
  char **degen;          /* 1/0, which syms inc which res [0..Kp-1][0..K-1] */
  int   *ndegen;	 /* # of degenerate residues per code  [0..Kp-1]    */
};
typedef struct esl_alphabet_s ESL_ALPHABET;

extern ESL_ALPHABET *esl_alphabet_Create(int type);
extern ESL_ALPHABET *esl_alphabet_CreateCustom(char *alphabet, int K, int Kp);
extern void    esl_alphabet_Destroy(ESL_ALPHABET *a);
extern int     esl_alphabet_SetSynonym(ESL_ALPHABET *a, char sym, char c);
extern int     esl_alphabet_SetCaseInsensitive(ESL_ALPHABET *a);
extern int     esl_alphabet_SetDegeneracy(ESL_ALPHABET *a, char c, char *ds);

extern int esl_abc_CreateDigitalSequence(ESL_ALPHABET *a, char *seq,
					 int L, char **ret_dsq);
extern int esl_abc_DigitizeSequence(ESL_ALPHABET *a, char *seq, 
				    int L, char *dsq);

extern int    esl_abc_AvgIScore(ESL_ALPHABET *a, char x, int *sc);
extern float  esl_abc_AvgFScore(ESL_ALPHABET *a, char x, float *sc);
extern double esl_abc_AvgDScore(ESL_ALPHABET *a, char x, double *sc);
extern int    esl_abc_ExpectIScore(ESL_ALPHABET *a, char x, int *sc, float *p);
extern float  esl_abc_ExpectFScore(ESL_ALPHABET *a, char x, float *sc,
				   float *p);
extern double esl_abc_ExpectDScore(ESL_ALPHABET *a, char x, double *sc,
				   double *p);

#define esl_abc_DigitizeSymbol(a, c) ((char) a->inmap[c])
#define esl_abc_IsDegenerate(a, x)   ((x) > a->K && (x) < a->Kp)
#define esl_abc_IsBasic(a, x)        ((x) >= 0   && (x) < a->K)
#define esl_abc_IsGap(a, x)          ((x) == a->K)

#endif /*!ESL_ALPHABET_INCLUDED*/
