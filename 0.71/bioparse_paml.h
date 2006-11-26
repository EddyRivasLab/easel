/* bioparse_paml.h
 * 
 * Parse datafiles from PAML:
 *   "Phylogenetic Analysis by Maximum Likelihood"
 *   Ziheng Yang
 *   http://abacus.gene.ucl.ac.uk/software/paml.html
 * 
 * SRE, Tue Jul 13 13:20:08 2004 [St. Louis]
 * SVN $Id$
 */

#ifndef ESL_BIOPARSE_PAML_INCLUDED
#define ESL_BIOPARSE_PAML_INCLUDED

#include <stdio.h>
#include <easel/dmatrix.h>

extern int esl_bio_ParsePAMLRateMatrix(FILE *fp, ESL_DMATRIX **ret_s, double **ret_pi);

#endif /*ESL_BIOPARSE_PAML_INCLUDED*/
