/* bioparse_paml.c
 * 
 * Parsers for datafiles from PAML:
 *   "Phylogenetic Analysis by Maximum Likelihood"
 *   Ziheng Yang
 *   http://abacus.gene.ucl.ac.uk/software/paml.html
 * 
 * SRE, Tue Jul 13 13:23:21 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <easel/easel.h>
#include <easel/dmatrix.h>
#include <easel/parse.h>
#include <easel/bioparse_paml.h>

/* Function:  esl_bio_ParsePAMLRateData()
 * Incept:    SRE, Fri Jul  9 09:27:24 2004 [St. Louis]
 *
 * Purpose:   Parse a rate matrix in PAML format, such as the
 *            Whelan and Goldman WAG matrix.
 *
 *            Format: First 190 numbers are a lower-triangular
 *            matrix of amino acid exchangeabilities s_ij.
 *            Next 20 numbers are the amino acid frequencies
 *            \pi. Remainder of the datafile is ignored.
 *            
 *            The alphabet order in the matrix and the frequency
 *            vector is assumed to be "ARNDCQEGHILKMFPSTWYV" 
 *            (alphabetical by three-letter code); this is
 *            transformed to Easel's "ACDEFGHIKLMNPQRSTVWY"
 *            (alphabetical by one-letter code) in the s_ij
 *            and \pi_i that are returned.
 *            
 * Args:      fp     - open datafile for reading.
 *            ret_s  - RETURN: s_ij matrix of amino acid exchangeabilities,
 *                     lower triangular only (only j<i valid),
 *                     in Easel order "ACDEFGHIKLMNPQRSTVWY".
 *                     Allocated here; caller must free.
 *                     Pass NULL if not wanted.
 *            ret_pi - RETURN: \pi_i vector of amino acid frequencies,
 *                     in Easel order "ACDEFGHIKLMNPQRSTVWY". 
 *                     Allocated here; caller must free.
 *                     Pass NULL if not wanted.
 *
 * Returns:   ESL_OK on success;
 *              s, pi are allocated here, caller is responsible for freeing.
 *            
 *            on failure:
 *            ESL_EMEM: memory allocation failure.
 *            ESL_EOF:  premature end of file, parse failed.
 *
 * Xref:      STL8/p.56.
 */
int
esl_bio_ParsePAMLRateData(FILE *fp, ESL_DMATRIX **ret_s, double **ret_pi)
{
  ESL_FILEPARSER *efp = NULL;
  ESL_DMATRIX    *s   = NULL;
  double         *pi  = NULL;
  char           *tok;
  int             status;
  int             i,j;
  char           *pamlalpha = "ARNDCQEGHILKMFPSTWYV";
  char           *eslalpha  = "ACDEFGHIKLMNPQRSTVWY";
  int             perm[20];
  char           *sptr;  

  status = esl_fileparse_create(fp, &efp);
  if (status != ESL_OK) goto FAILURE;

  status = esl_fileparse_set_commentchar(efp, '#');
  if (status != ESL_OK) goto FAILURE;

  if ((s    = esl_dmatrix_Create(20,20))   == NULL) { status = ESL_EMEM; goto FAILURE; }
  if ((pi   = malloc(sizeof(double) * 20)) == NULL) { status = ESL_EMEM; goto FAILURE; }

  /* constructs the alphabet permutation we need.
   * perm[i] -> original row/column i goes to row/column perm[i]
   */
   for (i = 0; i < 20; i++)
     perm[i] = (int) (strchr(eslalpha, pamlalpha[i]) - eslalpha);

   /* Read the s_ij matrix data in, permuting as we go.
    */
   esl_dmx_SetZero(s);
   for (i = 1; i < 20; i++)
    for (j = 0; j < i; j++)
      {
	if ((status = esl_fileparse_token(efp, &tok, NULL)) != ESL_OK) goto FAILURE;	
	s->mx[perm[i]][perm[j]] = atof(tok);
	s->mx[perm[j]][perm[i]] = s->mx[perm[i]][perm[j]];
      }

   /* Read the pi_i vector in, permuting as we read.
    */
  for (i = 0; i < 20; i++)
    {
      if ((status = esl_fileparse_token(efp, &tok, NULL)) != ESL_OK) goto FAILURE;	
      pi[perm[i]] = atof(tok);
    }

  esl_fileparse_free(efp);
  if (ret_s  != NULL) *ret_s  = s;  else esl_dmx_Free(s);
  if (ret_pi != NULL) *ret_pi = pi; else free(pi);
  return ESL_OK;

 FAILURE:
  if (efp != NULL) esl_fileparse_free(efp);
  if (s   != NULL) esl_dmx_Free(s);
  if (pi  != NULL) free(pi);
  return status;
}
