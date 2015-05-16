/* Genetic code tables for translation, whether canonical or non.
 */
#ifndef eslGENCODE_INCLUDED
#define eslGENCODE_INCLUDED

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"

typedef struct {
  int     trans_table;       // NCBI trans_table number, or -1
  char    desc[128];         // Description

  ESL_DSQ basic[64];         // Basic code table. aacode[0..63; pos1^16 + pos2^4 + pos3] = residue code for amino acid, 0..19. No degeneracies.
  int8_t  is_initiator[64];  // TRUE for allowed initiator codons; FALSE if not

  ESL_ALPHABET *nt_abc;      // A reference to nucleic alphabet that caller is maintaining elsewhere
  ESL_ALPHABET *aa_abc;      // A reference to amino alphabet that caller is maintaining 
} ESL_GENCODE;

extern ESL_GENCODE *esl_gencode_Create(ESL_ALPHABET *nt_abc, ESL_ALPHABET *aa_abc);
extern void         esl_gencode_Destroy(ESL_GENCODE *gcode);
extern int          esl_gencode_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *nucleic_abc, ESL_ALPHABET *amino_abc, ESL_GENCODE **ret_gcode);

extern int   esl_gencode_TranslateCodon(ESL_GENCODE *gcode, ESL_DSQ *dsq);

extern char *esl_gencode_DecodeDigicodon(ESL_GENCODE *gcode, int digicodon, char *codon);

#endif	/*eslGENCODE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
