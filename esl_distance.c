/* esl_distance.c
 * 
 * Distances between aligned sequence pairs, including both
 * probabilistic evolutionary models and ad hoc measures;
 * including both digital sequences (ESL_ALPHABET) and 
 * "analog" char sequences; and functions for calculating
 * complete NxN distance matrices from input alignments.
 * 
 *   
 * SRE, Mon Apr 17 20:05:43 2006 [St. Louis]
 * SVN $Id$
 */
#include <esl_config.h>

#include <ctype.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>

/* Function:  esl_dst_CPairId()
 * Incept:    SRE, Mon Apr 17 20:06:07 2006 [St. Louis]
 *
 * Purpose:   Calculates the pairwise fractional identity between two
 *            aligned character strings <asq1> and <asq2>.  Fractional
 *            identity is defined as <#idents / MIN(len1, len2)>,
 *            (where <len1> and <len2> are the number of residues in
 *            the two sequences, not counting gaps).  Optionally
 *            return this distance in <ret_pid>; optionally return the
 *            number of identities counted in <ret_nid>; and
 *            optionally return the denominator <MIN(len1,len2)> in
 *            <ret_n>.
 *            
 *            If <abc> is non-NULL, residues and identities are
 *            counted according to alphabet <abc>. This correctly 
 *            counts synonyms (such as T/U in nucleic acid alignments)
 *            and degeneracies (in DNA, N/T is counted as 0.25 of an
 *            identity; Y/T as 0.5 of an identity; R/T as 0; Y/Y as
 *            0.5; N/N as 0.25; and so on), as well as being 
 *            case-insensitive. 
 *            
 *            If <abc> is NULL, the comparison rule is simpler and not
 *            bio-alphabet-aware.  Any nonalphabetic character is
 *            assumed to be a gap symbol. Alphabetic symbols are
 *            compared for identity case-insensitively, but otherwise
 *            literally. Aligned pairs involving degeneracies and
 *            synonyms will generally be counted erroneously, but for
 *            a simple application with consistent, unambiguous sequences,
 *            this may not be of concern.
 *
 *            There are many ways to calculate pairwise identity,
 *            because there are a variety of choices for the
 *            denominator. Easel primarily uses percent identity
 *            calculations in ad hoc sequence weighting of multiple
 *            sequence alignments during profile HMM or profile SCFG
 *            construction (for example). We are therefore more
 *            concerned here about robustness to what real multiple
 *            alignments can throw at us, as opposed to correct
 *            phylogenetic distance inference. Multiple alignments
 *            often contain short sequence fragments, and we have to
 *            deal with cases where two short fragments may have
 *            little overlap (or none at all). The more
 *            phylogenetically motivated calculation of pairwise
 *            identity, <idents/(idents+mismat)> (the starting point
 *            for a Jukes-Cantor distance) is not robust enough,
 *            because alignments with few aligned residues (either
 *            because they are highly gappy, or they are partially
 *            overlapping fragments) might receive artifactually high
 *            identities. Two other ad hoc definitions,
 *            <idents/(AVG|MAX)(len1,len2)>, both have the
 *            disadvantage that alignments of fragments to longer
 *            sequences would have artifactually low identities.
 *            
 *            In the unusual case where <MIN(len1,len2)=0> -- that is,
 *            one of the sequences is completely blank -- the percent identity
 *            (0/0) is defined as 0. This is for robustness against
 *            length 0 sequences, which do arise in real applications.
 *
 * Args:      abc          - NULL, or the bioalphabet
 *            asq1         - aligned character string 1
 *            asq2         - aligned character string 2
 *            ret_pid      - optRETURN: pairwise identity, 0<=x<=1
 *            ret_nid      - optRETURN: # of identities
 *            ret_n        - optRETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <ret_pid>, <ret_nid>, <ret_n>
 *            contain the answers (for whichever were passed non-NULL). 
 *
 * Throws:    <eslECORRUPT> if either string contains an
 *            illegal non-sequence character;  
 *            <eslEINVAL> if the strings are different lengths
 *            (not aligned).
 */
int
esl_dst_CPairId(ESL_ALPHABET *abc, char *asq1, char *asq2, 
		double *ret_pid, int *ret_nid, int *ret_n)
{
  int     status;
  char    x1,x2;
  int     idents;               /* total identical positions  */
  int     len1, len2;           /* lengths of seqs            */
  int     i;                    /* position in aligned seqs   */

  idents = len1 = len2 = 0;
  if (abc == NULL)
    {
      for (i = 0; asq1[i] != '\0' && asq2[i] != '\0'; i++) 
	{
	  if (isalpha(asq1[i])) len1++;
	  if (isalpha(asq2[i])) len2++;
	  if (isalpha(asq1[i]) && isalpha(asq2[i])
	      && toupper(asq1[i]) == toupper(asq2[i])) idents++;
	}
  else
    {
      for (i = 0; asq1[i] != '\0' && asq2[i] != '\0'; i++) 
	{
	  x1 = esl_abc_DigitizeSymbol(abc, asq1[i]);
	  x2 = esl_abc_DigitizeSymbol(abc, asq2[i]);

	  if (x1 == eslILLEGAL_CHAR || x2 == eslILLEGAL_CHAR)
	    ESL_FAIL(eslEINVAL, "illegal 

	  if (esl_abc_XIsBasic(abc, x1)) len1++;
	  if (esl_abc_XIsBasic(abc, x2)) len2++;

	  if (esl_abc_XIsBasic(x1) && esl_abc_XIsBasic(x2) && x1 == x2)
	    idents++;
	}
    }
  if (len2 < len1) len1 = len2;

  if (asq1[i] != '\0' || asq2[i] != '\0') 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");

  if (ret_distance != NULL)  *ret_distance = ( len1==0 ? 0. : (double) idents / (double) len1 );
  if (ret_nid      != NULL)  *ret_nid      = idents;
  if (ret_n        != NULL)  *ret_n        = len1;
  return eslOK;

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = 0.;
  if (ret_nid      != NULL)  *ret_nid      = 0;
  if (ret_n        != NULL)  *ret_n        = 0;
  return status;
}


/* Function:  esl_dst_XPairId()
 * Incept:    SRE, Tue Apr 18 09:24:05 2006 [St. Louis]
 *
 * Purpose:   Digital version of <esl_dst_PairId()>: <adsq1> and
 *            <adsq2> are digitized aligned sequences, in alphabet
 *            <abc>. Otherwise, same as <esl_dst_PairId()>.
 *            
 * Args:      adsq1        - aligned digital seq 1
 *            adsq2        - aligned digital seq 2
 *            ret_distance - optRETURN: pairwise identity, 0<=x<=1
 *            ret_nid      - optRETURN: # of identities
 *            ret_n        - optRETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <ret_distance>, <ret_nid>, <ret_n>
 *            contain the answers (for whichever were passed non-NULL). 
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 *            <eslEDIVZERO> if MIN(len1,len2) is 0. *            
 */
int
esl_dst_XPairId(ESL_ALPHABET *abc, char *adsq1, char *adsq2, 
		double *ret_distance, int *ret_nid, int *ret_n)
{
  int     status;
  int     idents;               /* total identical positions  */
  int     len1, len2;           /* lengths of seqs            */
  int     i;                    /* position in aligned seqs   */

  idents = len1 = len2 = 0;
  for (i = 1; adsq1[i] != eslSENTINEL && adsq2[i] != eslSENTINEL; i++) 
    {
      if (esl_abc_XIsBasic(abc, adsq1[i])) len1++;
      if (esl_abc_XIsBasic(abc, adsq2[i])) len2++;

      if (esl_abc_XIsBasic(asdq1[i]) && esl_abc_XIsBasic(asdq2[i])
	  && adsq1[i] == adsq2[i])
	idents++;
    }
  if (len2 < len1) len1 = len2;

  if (adsq1[i] != eslSENTINEL || adsq2[i] != eslSENTINEL) 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");

  if (ret_distance != NULL)  *ret_distance = (double) idents / (double) len1;
  if (ret_nid      != NULL)  *ret_nid      = idents;
  if (ret_n        != NULL)  *ret_n        = len1;
  return eslOK;

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = 0.;
  if (ret_nid      != NULL)  *ret_nid      = 0;
  if (ret_n        != NULL)  *ret_n        = 0;
  return status;
}


/* jukescantor()
 * 
 * The generalized Jukes/Cantor distance calculation.
 * 
 * 
 */
static int
jukescantor(int n1, int n2, int alphabet_size, double *ret_distance, double *ret_variance)
{
  int    status;
  double D, K, N;
  double x;
  double distance, variance;

  if (n1+n2 == 0) { status = eslEDIVZERO; goto FAILURE; }

  K = (double) alphabet_size;
  D = (double) n2 / (double) (n1+n2);
  N = (double) (n1+n2);

  x = 1. - D * K/(K-1.);
  if (x <= 0.) 
    {
      distance = HUGE_VAL;
      variance = HUGE_VAL;
    }
  else
    {
      distance =   -log(x) * K/(K-1);
      variance =  exp( 2.*K*distance/(K-1) ) * D * (1.-D) / N;
    }
  if (*ret_distance != NULL)  *ret_distance = distance;
  if (*ret_variance != NULL)  *ret_variance = variance;
  return eslOK;


 FAILURE:
  if (*ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (*ret_variance != NULL)  *ret_variance = HUGE_VAL;
  return status;
}

/* Function:  esl_dst_CJukesCantor()
 * Incept:    SRE, Tue Apr 18 14:00:37 2006 [St. Louis]
 *
 * Purpose:   Calculate the generalized Jukes-Cantor distance between two
 *            aligned character strings <as1> and <as2>, in substitutions/site, 
 *            using alphabet <abc> to evaluate identities and differences.
 *            The maximum likelihood estimate for the distance is returned in
 *            <ret_distance>. The large-sample variance for the distance
 *            estimate is optionally returned in <ret_variance>.
 *            
 *            Only aligned pairs of unambiguous residues (in <abc>)
 *            are counted towards identities (n1) and substitutions
 *            (n2) (including synonyms such as U/T, for a nucleic acid
 *            alphabet). Pairs that involve a gap symbol or degeneracy
 *            are ignored. The fractional difference <D> is <n1/n1+n2>.
 *            
 *            The alphabet <abc> is required: we to know the alphabet
 *            size <abc->K> to calculate a generalized Jukes-Cantor
 *            distance.
 *            
 *            A Jukes-Cantor model assumes that all positions are
 *            substituted at the same rate $\alpha$. It implies
 *            equiprobable stationary probabilities.
 *            The calculation is:
 *            
 *            $d = -\frac{K-1}{K} \log \left( 1 - \frac{K}{K-1} D \right)$
 *                                                
 *            where $D$ is the fractional difference, and $K$ is the
 *            alphabet size. The variance is:
 *            
 *            $\sigma^2 = e^\frac{2Kd}{K-1} \frac{D(1-D)}{N}
 *            
 *            where $N$ is the total number of columns counted, <n1+n2>.
 *
 * Args:      abc          - bioalphabet to use for comparisons
 *            as1          - 1st aligned seq, 0..L-1, \0-terminated
 *            as2          - 2nd aligned seq, 0..L-1, \0-terminated 
 *            ret_distance - RETURN: ML estimate of distance d
 *            ret_variance - optRETURN: large-sample variance of d
 *
 * Returns:   <eslOK> on success.
 * 
 *            Infinite distances are quite possible, in which case
 *            distance and variance are both <HUGE_VAL>. Caller
 *            has to deal with this case as it sees fit.
 *
 * Throws:    <eslEINVAL> if the two strings aren't the same length (and
 *            thus can't have been properly aligned).
 *            <eslEDIVZERO> if no aligned residues were counted.
 */
int
esl_dst_CJukesCantor(ESL_ALPHABET *abc, char *as1, char *as2, 
		     double *ret_distance, double *ret_variance)
{
  int     status;
  char    x1,x2;		/* digitized residues */
  int     n1, n2;               /* number of observed identities, substitutions */
  int     i;                    /* position in aligned seqs   */

  /* 1. Count identities, mismatches.
   */
  n1 = n2 = 0;
  for (i = 0; as1[i] != '\0' && as2[i] != '\0'; i++) 
    {
      x1 = esl_abc_DigitizeSymbol(abc, as1[i]);
      x2 = esl_abc_DigitizeSymbol(abc, as2[i]);
      if (esl_abc_XIsBasic(x1) && esl_abc_XIsBasic(x2))
	{
	  if (x1 == x2) n1++;
	  else          n2++;
	}
    }
  if (as1[i] != '\0' || as2[i] != '\0') 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");
  
  return jukescantor(n1, n2, abc->K, ret_distance, ret_variance);

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = -HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = -HUGE_VAL;
  return status;
}


/* Function:  esl_dst_XJukesCantor()
 * Incept:    SRE, Tue Apr 18 15:26:51 2006 [St. Louis]
 *
 * Purpose:   Calculate the generalized Jukes-Cantor distance between two
 *            aligned digital strings <ax> and <ay>, in substitutions/site, 
 *            using alphabet <abc> to evaluate identities and differences.
 *            The maximum likelihood estimate for the distance is returned in
 *            <ret_distance>. The large-sample variance for the distance
 *            estimate is optionally returned in <ret_variance>.
 *            
 *            Identical to <esl_dst_CJukesCantor()>, except that it takes
 *            digital sequences instead of character strings.
 *
 * Args:      abc          - bioalphabet to use for comparisons
 *            ax           - 1st digital aligned seq
 *            ay           - 2nd digital aligned seq
 *            ret_distance - RETURN: ML estimate of distance d
 *            ret_variance - optRETURN: large-sample variance of d
 *
 * Returns:   <eslOK> on success. As in <esl_dst_CJukesCantor()>, the
 *            distance and variance may be infinite, in which case they
 *            are returned as <HUGE_VAL>.
 *
 * Throws:    <eslEINVAL> if the two strings aren't the same length (and
 *            thus can't have been properly aligned).
 *            <eslEDIVZERO> if no aligned residues were counted.
 *            On either failure, the distance and variance are set
 *            to <HUGE_VAL>.
 */
int
esl_dst_XJukesCantor(ESL_ALPHABET *abc, char *ax, char *ay, 
		     double *ret_distance, double *ret_variance)
{
  int     status;
  int     n1, n2;               /* number of observed identities, substitutions */
  int     i;                    /* position in aligned seqs   */

  n1 = n2 = 0;
  for (i = 1; ax[i] != eslSENTINEL && ay != eslSENTINEL; i++) 
    {
      if (esl_abc_XIsBasic(ax[i]) && esl_abc_XIsBasic(ay[i]))
	{
	  if (ax[i] == ay[i]) n1++;
	  else                n2++;
	}
    }
  if (ax[i] != '\0' || ay[i] != '\0') 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");
  
  return jukescantor(n1, n2, abc->K, ret_distance, ret_variance);

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
  return status;
}



/* Function:  esl_dst_CJukesCantorMx()
 * Incept:    SRE, Tue Apr 18 16:00:16 2006 [St. Louis]
 *
 * Purpose:   Given a multiple sequence alignment <aseq>, consisting of
 *            <nseq> aligned character sequences in bioalphabet <abc>,
 *            calculate a symmetric Jukes/Cantor pairwise distance
 *            matrix for all sequence pairs, and return the distance
 *            matrix in <ret_D>. Optionally also return the
 *            large-sample variances for those ML distance estimates
 *            in <ret_V>.
 *            
 *            Infinite distances (and variances) are possible; they
 *            are represented as <HUGE_VAL> in <D> and <V>. Caller must
 *            be prepared to deal with them as appropriate.
 *
 * Args:      abc
 *            aseq
 *            nseq
 *            ret_D
 *            ret_V
 *
 * Returns:   <eslOK> on success. <D> (and optionally <V>) contain the
 *            distance matrix (and variances); caller frees these with
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if any pair of sequences have differing lengths
 *            (and thus cannot have been properly aligned). 
 *            <eslEDIVZERO> 
 *
 * Xref:      
 */


int
esl_dst_CJukesCantorMx(ESL_ALPHABET *abc, char **aseq, int nseq, 
		       ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V)
{
  int          status;
  ESL_DMATRIX *D = NULL;
  ESL_DMATRIX *V = NULL;
  int          i,j;

  D = esl_dmatrix_Create(nseq, nseq);   if (D == NULL) goto FAILURE;
  V = esl_dmatrix_Create(nseq, nseq);   if (V == NULL) goto FAILURE;
  for (i = 0; i < nseq; i++)
    {
      D->mx[i][i] = 0.;
      V->mx[i][i] = 0.;
      for (j = i+1; j < nseq; j++)
	{
	  status = esl_dst_CJukesCantor(abc, aseq[i], aseq[j], 
					&(D->mx[i][j]), &(V->mx[i][j]));
	  if (status != eslOK) 
	    {
	      esl_error(status, "J/C calculation failed at seqs %d,%d", i,j);
	      goto FAILURE;
	    }
	  D->mx[j][i] = D->mx[i][j];
	  V->mx[j][i] = V->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D;  else esl_dmatrix_Destroy(D);
  if (ret_V != NULL) *ret_D = V;  else esl_dmatrix_Destroy(V);
  return eslOK;

 FAILURE:
  if (D     != NULL) esl_dmatrix_Destroy(D);
  if (V     != NULL) esl_dmatrix_Destroy(V);
  if (ret_D != NULL) *ret_D = NULL;
  if (ret_V != NULL) *ret_V = NULL;
  return status;
}



 /*****************************************************************
  * @LICENSE@
  *****************************************************************/

