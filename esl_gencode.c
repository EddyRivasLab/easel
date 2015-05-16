/* Genetic code tables for translation, whether canonical or non.
 * 
 * Table of contents:
 *   1. NCBI genetic code tables, in Easel digital form
 *   2. ESL_GENCODE object
 *   3. Parsing NCBI genetic code data 
 *   4. Debugging/development utilities
 *   5. Examples
 *   6. Copyright and license
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"
#include "esl_regexp.h"
#include "esl_gencode.h"


/*****************************************************************
 * 1. Genetic code tables, from NCBI
 *****************************************************************/

/* 
 * From: http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
 * Digitized by the example driver.
 * 
 * The NCBI page has useful information about these code tables, references and caveats.
 */

static const ESL_GENCODE ncbi_trans_tables[] = {
  { 1, "standard",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },
  
  { 2, "vertebrate mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 27, 15, 27, 15, 10,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   *   S   *   S   M   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 3, "yeast mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 10,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14, 16, 16, 16, 16,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   M   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   T   T   T   T   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 4, "mold, protozoan, coelenterate mitochondrial; Mycoplasma/Spiroplasma",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 5, "invertebrate mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 10,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   S   S   S   S   M   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },
  
  { 6, "ciliate, dasycladacean, Hexamita nuclear",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 13, 19, 13, 19, 15, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   Q   Y   Q   Y   S   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },

  { 9, "echinoderm and flatworm mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    { 11, 11,  8, 11, 16, 16, 16, 16, 15, 15, 15, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   N   N   K   N   T   T   T   T   S   S   S   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 10, "euplotid nuclear",
   /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15,  1,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   C   C   W   C   L   F   L   F */
    NULL, NULL },

  { 11, "bacterial, archaeal, and plant plastid",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },

  { 12, "alternative yeast", 
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9, 15,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   S   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },

  { 13, "ascidian mitochondrial", 
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16,  5, 15,  5, 15, 10,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   G   S   G   S   M   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 14, "alternative flatworm mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    { 11, 11,  8, 11, 16, 16, 16, 16, 15, 15, 15, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 19, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   N   N   K   N   T   T   T   T   S   S   S   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   Y   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 16, "chlorophycean mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19,  9, 19, 15, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   L   Y   S   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },

  { 21, "trematode mitochondrial", 
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    { 11, 11,  8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 10,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   N   N   K   N   T   T   T   T   S   S   S   S   M   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 22, "Scenedesmus obliquus mitochondrial",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19,  9, 19, 27, 15, 15, 15, 27,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   L   Y   *   S   S   S   *   C   W   C   L   F   L   F */
    NULL, NULL },

  { 23, "Thraustochytrium mitochondrial", 
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27,  1, 18,  1, 27,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   *   C   W   C   *   F   L   F */
    NULL, NULL },

  { 24, "Pterobranchia mitochondrial", 
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 15, 15,  8, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   S   S   K   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   W   C   W   C   L   F   L   F */
    NULL, NULL },

  { 25, "Candidate Division SR1 and Gracilibacteria",
  /* AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT */
    {  8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7, 13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9,  3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15,  5,  1, 18,  1,  9,  4,  9,  4 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
  /*   K   N   K   N   T   T   T   T   R   S   R   S   I   I   M   I   Q   H   Q   H   P   P   P   P   R   R   R   R   L   L   L   L   E   D   E   D   A   A   A   A   G   G   G   G   V   V   V   V   *   Y   *   Y   S   S   S   S   G   C   W   C   L   F   L   F */
    NULL, NULL },
};


/*****************************************************************
 * 2. The ESL_GENCODE genetic code object
 *****************************************************************/


ESL_GENCODE *
esl_gencode_Create(ESL_ALPHABET *nt_abc, ESL_ALPHABET *aa_abc)
{
  ESL_GENCODE *gcode = NULL;
  int digicodon;
  int status;

  ESL_ALLOC(gcode, sizeof(ESL_GENCODE));

  /* Default = standard code (NCBI trans table 1) */
  ESL_DASSERT1((  ncbi_trans_tables[0].trans_table == 1 ));
  gcode->trans_table =  ncbi_trans_tables[0].trans_table;
  strcpy(gcode->desc,   ncbi_trans_tables[0].desc);
  for (digicodon = 0; digicodon < 64; digicodon++)
    {
      gcode->basic[digicodon]        = ncbi_trans_tables[0].basic[digicodon];
      gcode->is_initiator[digicodon] = ncbi_trans_tables[0].is_initiator[digicodon];
    }
  gcode->nt_abc = nt_abc;  // Keep a reference to the nucleic alphabet; caller remains responsible for it
  gcode->aa_abc = aa_abc;  //  ditto for amino alphabet
  return gcode;

 ERROR:
  esl_gencode_Destroy(gcode);
  return NULL;
}


void
esl_gencode_Destroy(ESL_GENCODE *gcode)
{
  if (gcode)
    free(gcode);
}


/*****************************************************************
 * 3. Parsing NCBI genetic code data
 *****************************************************************/


/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *            Example of an NCBI genetic code datafile:
 * 
 *            AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
 *            Starts = ---M---------------M---------------M----------------------------
 *            Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
 *            Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
 *            Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 *
 *            This function is, and must remain independent of the order of
 *            residues in the amino and nucleic alphabets. This allows us to
 *            convert NCBI genetic code text files to digitized Easel translation
 *            tables even if we someday change the Easel alpahabet(s). Once 
 *            digitized, Easel encodings of the genetic code are dependent on
 *            the <eslAMINO> and <eslNUCLEIC> alphabets they were created with.
 *
 * Args:      
 *
 * Returns:   
 *            
 *            On a parse error, returns <eslEFORMAT>, and an informative message is
 *            left in <efp->errbuf>.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
esl_gencode_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *nucleic_abc, ESL_ALPHABET *amino_abc, ESL_GENCODE **ret_gcode)
{
  ESL_GENCODE *gcode = esl_gencode_Create(nucleic_abc, amino_abc);
  ESL_REGEXP  *mach  = esl_regexp_Create();
  int   start, end, s, e;
  char  aas[65];
  char  mline[65];
  char  base1[65];
  char  base2[65];
  char  base3[65];
  int   aa_seen[20];
  int   stop_seen;
  int   codon_seen[64];
  int   x, codon, pos;
  int   status;

  ESL_DASSERT1(( nucleic_abc->K == 4  ));  // We're going to hardcode ncodons = 64, so "trust but verify"
  ESL_DASSERT1(( amino_abc->K   == 20 ));

  if ((status = esl_fileparser_NextLine(efp))                                           != eslOK)  {  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "File empty or truncated? No AAs line found");  else goto ERROR; }
  if ((status = esl_regexp_Match(mach, "^\\s*[Aa][Aa]s\\s*=\\s*(\\S+)\\s*$", efp->buf)) != eslOK)  {  if (status == eslEOD) ESL_XFAIL(eslEFORMAT, efp->errbuf, "First data line doesn't start with 'AAs ='");  else goto ERROR; }
  if ((status = esl_regexp_SubmatchCoords(mach, efp->buf, 1, &start, &end))             != eslOK)  goto ERROR;  
  if (end - start + 1 != 64) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected 64 char of AAs data");
  strncpy(aas, efp->buf+start, 64); 
  aas[64] = '\0';
  
  if ((status = esl_fileparser_NextLine(efp))                                           != eslOK)  {  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "File empty or truncated? No Starts line found");   else goto ERROR; }
  if ((status = esl_regexp_Match(mach, "^\\s*[Ss]tarts\\s*=\\s*(\\S+)\\s*$", efp->buf)) != eslOK)  {  if (status == eslEOD) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Second data line doesn't start with 'Starts ='");  else goto ERROR; }
  if ((status = esl_regexp_SubmatchCoords(mach, efp->buf, 1, &s, &e))                   != eslOK)  goto ERROR;
  if (e - s + 1 != 64)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected 64 char of Starts data");
  if (s != start)       ESL_XFAIL(eslEFORMAT, efp->errbuf, "Starts data is not aligned with AAs data above it");
  strncpy(mline, efp->buf+start, 64); 
  mline[64] = '\0';

  if ((status = esl_fileparser_NextLine(efp))                                           != eslOK)  {  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "File empty or truncated? No Base1 line found");  else goto ERROR; }
  if ((status = esl_regexp_Match(mach, "^\\s*[Bb]ase1\\s*=\\s*(\\S+)\\s*$", efp->buf))  != eslOK)  {  if (status == eslEOD) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Third data line doesn't start with 'Base1 ='");  else goto ERROR; }
  if ((status = esl_regexp_SubmatchCoords(mach, efp->buf, 1, &s, &e))                   != eslOK)  goto ERROR;
  if (e - s + 1 != 64)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected 64 char of Base1 data");
  if (s != start)       ESL_XFAIL(eslEFORMAT, efp->errbuf, "Base1 data is not aligned with data above it");
  strncpy(base1, efp->buf+start, 64); 
  base1[64] = '\0';

  if ((status = esl_fileparser_NextLine(efp))                                           != eslOK)  {  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "File empty or truncated? No Base2 line found");  else goto ERROR; }
  if ((status = esl_regexp_Match(mach, "^\\s*[Bb]ase2\\s*=\\s*(\\S+)\\s*$", efp->buf))  != eslOK)  {  if (status == eslEOD) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Fourth data line doesn't start with 'Base2 ='"); else goto ERROR; }
  if ((status = esl_regexp_SubmatchCoords(mach, efp->buf, 1, &s, &e))                   != eslOK)  goto ERROR;
  if (e - s + 1 != 64)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected 64 char of Base2 data");
  if (s != start)       ESL_XFAIL(eslEFORMAT, efp->errbuf, "Base2 data is not aligned with data above it");
  strncpy(base2, efp->buf+start, 64); 
  base2[64] = '\0';

  if ((status = esl_fileparser_NextLine(efp))                                           != eslOK)  {  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "File empty or truncated? No Base3 line found"); else goto ERROR; }
  if ((status = esl_regexp_Match(mach, "^\\s*[Bb]ase3\\s*=\\s*(\\S+)\\s*$", efp->buf))  != eslOK)  {  if (status == eslEOD) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Fifth data line doesn't start with 'Base3 ='"); else goto ERROR; }
  if ((status = esl_regexp_SubmatchCoords(mach, efp->buf, 1, &s, &e))                   != eslOK)  goto ERROR;
  if (e - s + 1 != 64)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected 64 char of Base3 data");
  if (s != start)       ESL_XFAIL(eslEFORMAT, efp->errbuf, "Base3 data is not aligned with data above it");
  strncpy(base3, efp->buf+start, 64); 
  base3[64] = '\0';

  stop_seen = FALSE;
  for (    x = 0;     x < 20;     x++)    aa_seen[x] = FALSE;
  for (codon = 0; codon < 64; codon++) codon_seen[x] = FALSE;
  
  for (pos = 0; pos < 64; pos++)
    {
      if (! esl_abc_CIsValid(amino_abc,   aas[pos])   || ! (esl_abc_CIsCanonical(amino_abc, aas[pos]) || esl_abc_CIsNonresidue(amino_abc, aas[pos])))  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Character %c on AAs line is not an amino acid or a * (stop)", aas[pos]);
      if (! esl_abc_CIsValid(nucleic_abc, base1[pos]) || ! esl_abc_CIsCanonical(nucleic_abc, base1[pos]))                                              ESL_XFAIL(eslEFORMAT, efp->errbuf, "Character %c on Base1 line is not a nucleotide", base1[pos]);
      if (! esl_abc_CIsValid(nucleic_abc, base2[pos]) || ! esl_abc_CIsCanonical(nucleic_abc, base2[pos]))                                              ESL_XFAIL(eslEFORMAT, efp->errbuf, "Character %c on Base2 line is not a nucleotide", base2[pos]);
      if (! esl_abc_CIsValid(nucleic_abc, base3[pos]) || ! esl_abc_CIsCanonical(nucleic_abc, base3[pos]))                                              ESL_XFAIL(eslEFORMAT, efp->errbuf, "Character %c on Base3 line is not a nucleotide", base3[pos]);
      if ( mline[pos] != '-' && mline[pos] != 'm' && mline[pos] != 'M')                                                                                ESL_XFAIL(eslEFORMAT, efp->errbuf, "Character %c on Starts line is neither a - or an M", mline[pos]);

      codon = 16 * esl_abc_DigitizeSymbol(nucleic_abc, base1[pos]) +
	       4 * esl_abc_DigitizeSymbol(nucleic_abc, base2[pos]) +
                   esl_abc_DigitizeSymbol(nucleic_abc, base3[pos]);
      x    = esl_abc_DigitizeSymbol(amino_abc, aas[pos]);

      ESL_DASSERT1(( codon >= 0 && codon < 64 ));
      ESL_DASSERT1(( x >= 0 && (x < 20 || x == esl_abc_XGetNonresidue(amino_abc))));

      if (x < 20) aa_seen[x]++; else stop_seen++;
      codon_seen[codon]++;
      
      gcode->basic[codon]        = x;
      gcode->is_initiator[codon] = ( mline[pos] == '-' ? FALSE : TRUE );   // We already checked above that it's one of "-mM"
    }

  /* A genetic code must provide a translation for all 64 codons, and
   * all 20 amino acids to be encoded. (No organism is yet known to
   * encode fewer than 20 amino acids [Kawahara-Kobayashi et al, NAR
   * 40:10576, 2012].) The code must include at least one stop codon.
   */
  if (! stop_seen)           ESL_XFAIL(eslEFORMAT, efp->errbuf, "No stop codon found in that genetic code");
  for (codon = 0; codon < 64; codon++)
    if (! codon_seen[codon]) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Data for fewer than 64 codons was found");
  for (x = 0; x < 20; x++)
    if (aa_seen[x] == 0)     ESL_XFAIL(eslEFORMAT, efp->errbuf, "No codon for residue %c found", amino_abc->sym[x]);

  
  esl_regexp_Destroy(mach);
  gcode->trans_table = -1;          // It was initialized to 1, the NCBI standard table; reset
  gcode->desc[0]     = '\0';        // Was initialized to desc of NCBI table 1; blank it
  gcode->nt_abc      = nucleic_abc;
  gcode->aa_abc      = amino_abc;
  *ret_gcode         = gcode;
  return eslOK;

 ERROR:
  if (gcode) esl_gencode_Destroy(gcode);
  if (mach)  esl_regexp_Destroy(mach);
  *ret_gcode = NULL;
  return status;
}


/*****************************************************************
 * x. Other API tools
 *****************************************************************/


/* Translate the digital codon dsq[0..2] */
int
esl_gencode_TranslateCodon(ESL_GENCODE *gcode, ESL_DSQ *dsq)
{
  ESL_DSQ x, y, z;
  int     codon;
  int     aa = -1;

  if (esl_abc_XIsCanonical(gcode->nt_abc, dsq[0]) && esl_abc_XIsCanonical(gcode->nt_abc, dsq[1]) && esl_abc_XIsCanonical(gcode->nt_abc, dsq[2]))
    {
      codon = 16*dsq[0] + 4*dsq[1] + dsq[2];
      return gcode->basic[codon];
    }

  for (x = 0; x < 20; x++)
    {
      if (! gcode->nt_abc->degen[dsq[0]][x]) continue;
      for (y = 0; y < 20; y++)
	{
	  if (! gcode->nt_abc->degen[dsq[1]][y]) continue;
	  for (z = 0; z < 20; z++)
	    {
	      if (! gcode->nt_abc->degen[dsq[2]][z]) continue;
	      /* xyz is one possible basic codon included in the dsq[3] degeneracy */
	      codon = x * 16 + y * 4 + z;
	      if      (aa == -1) aa = gcode->basic[codon];
	      else if (aa != gcode->basic[codon]) return esl_abc_XGetUnknown(gcode->aa_abc);
	    }
	}
    }
  return aa;
}


/*****************************************************************
 * 4. Debugging/development utilities
 *****************************************************************/ 

/* <codon> must be allocated for at least 4 chars:
 * the codon + \0. The returned ptr is <codon> itself,
 * which allows caller to use <esl_gencode_DecodeDigicodon()>
 * directly as an arg in printf(), etc.
 */
char *
esl_gencode_DecodeDigicodon(ESL_GENCODE *gcode, int digicodon, char *codon)
{
  codon[0] = gcode->nt_abc->sym[ digicodon / 16 ];
  codon[1] = gcode->nt_abc->sym[ (digicodon % 16) / 4 ];
  codon[2] = gcode->nt_abc->sym[ digicodon % 4 ];
  codon[3] = '\0';
  return codon;
}



/****************************************************************
 * 5. Example
 ****************************************************************/

#ifdef eslGENCODE_EXAMPLE
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"
#include "esl_gencode.h"

#include <stdio.h>

int
main(int argc, char **argv)
{
  char           *codefile = argv[1];
  ESL_FILEPARSER *efp      = NULL;
  ESL_GENCODE    *gcode    = NULL;
  ESL_ALPHABET   *nt_abc   = esl_alphabet_Create(eslDNA);
  ESL_ALPHABET   *aa_abc   = esl_alphabet_Create(eslAMINO);
  int  digicodon;
  char codon[4];
  int  status;

  if (esl_fileparser_Open(codefile, /*env=*/NULL, &efp) != eslOK) esl_fatal("Failed to open code file %s", codefile);
  esl_fileparser_SetCommentChar(efp, '#');

  status = esl_gencode_Read(efp, nt_abc, aa_abc, &gcode);
  if      (status == eslEFORMAT) esl_fatal("Failed to parse genetic code datafile %s\n  %s\n", codefile, efp->errbuf);
  else if (status != eslOK)      esl_fatal("Unexpected failure parsing genetic code datafile %s : code %d\n", codefile, status);

  printf("/* ");
  for (digicodon = 0; digicodon < 64; digicodon++)
    printf("%3s ", esl_gencode_DecodeDigicodon(gcode, digicodon, codon));
  printf("*/\n");

  printf("  {");
  for (digicodon = 0; digicodon < 64; digicodon++)
    printf("%3d%c", gcode->basic[digicodon], (digicodon < 63 ? ',' : ' '));
  printf("},\n");

  printf("  {");
  for (digicodon = 0; digicodon < 64; digicodon++)
    printf("%3d%c", gcode->is_initiator[digicodon], (digicodon < 63 ? ',' : ' '));
  printf("},\n");

  printf("/* ");
  for (digicodon = 0; digicodon < 64; digicodon++)
    printf("  %c ", gcode->aa_abc->sym [gcode->basic[digicodon]]);
  printf("*/\n");

  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_gencode_Destroy(gcode);
  esl_fileparser_Close(efp);
}
#endif /*eslGENCODE_EXAMPLE*/


/****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 ****************************************************************/
