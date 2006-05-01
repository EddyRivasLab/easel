/* esl_distance.c
 * 
 * Distances between aligned sequences, including both
 * probabilistic evolutionary models and ad hoc measures;
 * including both digital sequences (ESL_ALPHABET) and 
 * "analog" char sequences; and functions for calculating
 * complete NxN distance matrices from input alignments.
 * 
 * SVN $Id$
 * SRE, Mon Apr 17 20:05:43 2006 [St. Louis]
 */
#include <esl_config.h>

#include <ctype.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_dmatrix.h>


/* Function:  esl_dst_CPairId()
 * Incept:    SRE, Mon Apr 17 20:06:07 2006 [St. Louis]
 *
 * Purpose:   Calculates pairwise fractional identity between two
 *            aligned character strings <asq1> and <asq2>.  Fractional
 *            identity is defined as <#idents / MIN(len1, len2)>,
 *            where <len1> and <len2> are the number of residues in
 *            the two sequences, not counting gaps.  Optionally
 *            return this distance in <ret_pid>; optionally return the
 *            number of identities counted in <ret_nid>; and
 *            optionally return the denominator <MIN(len1,len2)> in
 *            <ret_n>.
 *            
 *            If alphabet <abc> is non-NULL, residues and identities
 *            are counted in a 'bio-aware' mode using <abc>. This
 *            correctly counts synonyms (such as T/U in nucleic acid
 *            alignments) and degeneracies (in DNA, N/T is counted as
 *            0.25 of an identity; Y/T as 0.5 of an identity; R/T as
 *            0; Y/Y as 0.5; N/N as 0.25; and so on), as well as being
 *            case-insensitive.
 *            
 *            If <abc> is NULL, the comparison rule is simpler and not
 *            bio-alphabet-aware.  Any nonalphabetic character is
 *            assumed to be a gap symbol. Alphabetic symbols are
 *            compared for identity literally, albeit
 *            case-insensitively. Thus aligned pairs involving
 *            degeneracies and synonyms will generally be counted
 *            erroneously, but for a simple application this may not
 *            be of concern.
 *
 *            There are many ways to calculate pairwise identity
 *            because there are a variety of choices for the
 *            denominator. We primarily use percent identity
 *            calculations in ad hoc sequence weighting of multiple
 *            sequence alignments during profile HMM or profile SCFG
 *            construction. We are therefore more concerned here about
 *            robustness to what real multiple alignments can throw at
 *            us, as opposed to correct phylogenetic distance
 *            inference. Multiple alignments often contain short
 *            sequence fragments, and we have to deal with cases where
 *            two short fragments may have little overlap (or none at
 *            all). The more phylogenetically 'correct' calculation of
 *            pairwise identity, <idents/(idents+mismat)> -- the
 *            starting point for a Jukes-Cantor distance -- is not
 *            robust enough, because alignments with few aligned
 *            residues (either because they are highly gappy, or they
 *            are partially overlapping fragments) might receive
 *            artifactually high identities. Two other ad hoc
 *            definitions, <idents/(AVG|MAX)(len1,len2)>, both have
 *            the disadvantage that alignments of fragments to longer
 *            sequences would have artifactually low identities.
 *            
 *            In the unusual case where <MIN(len1,len2)=0> -- that is,
 *            one of the sequences is completely gaps -- the percent
 *            identity (0/0) is defined as 0. The calculation is then
 *            robust against length 0 sequences, which do arise in
 *            real applications.
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
	      && toupper(asq1[i]) == toupper(asq2[i])) 
	    idents++;
	}
    }
  else
    {
      for (i = 0; asq1[i] != '\0' && asq2[i] != '\0'; i++) 
	{
	  x1 = esl_abc_DigitizeSymbol(abc, asq1[i]);
	  x2 = esl_abc_DigitizeSymbol(abc, asq2[i]);

	  if (x1 == eslILLEGAL_CHAR) ESL_FAIL(eslECORRUPT, "illegal character %c in seq 1", asq1[i]);
	  if (x2 == eslILLEGAL_CHAR) ESL_FAIL(eslECORRUPT, "illegal character %c in seq 2", asq2[i]);

	  if (esl_abc_XIsBasic(abc, x1)) len1++;
	  if (esl_abc_XIsBasic(abc, x2)) len2++;

	  if (esl_abc_XIsBasic(abc, x1) && esl_abc_XIsBasic(abc, x2) && x1 == x2)
	    idents++;
	}
    }
  if (len2 < len1) len1 = len2;

  if (asq1[i] != '\0' || asq2[i] != '\0') 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");

  if (ret_pid  != NULL)  *ret_pid = ( len1==0 ? 0. : (double) idents / (double) len1 );
  if (ret_nid  != NULL)  *ret_nid = idents;
  if (ret_n    != NULL)  *ret_n   = len1;
  return eslOK;

 FAILURE:
  if (ret_pid  != NULL)  *ret_pid = 0.;
  if (ret_nid  != NULL)  *ret_nid = 0;
  if (ret_n    != NULL)  *ret_n   = 0;
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
 *            ret_pid      - optRETURN: pairwise identity, 0<=x<=1
 *            ret_nid      - optRETURN: # of identities
 *            ret_n        - optRETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <ret_distance>, <ret_nid>, <ret_n>
 *            contain the answers (for whichever were passed non-NULL). 
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 */
int
esl_dst_XPairId(ESL_ALPHABET *abc, char *ax1, char *ax2, 
		double *ret_distance, int *ret_nid, int *ret_n)
{
  int     status;
  int     idents;               /* total identical positions  */
  int     len1, len2;           /* lengths of seqs            */
  int     i;                    /* position in aligned seqs   */

  idents = len1 = len2 = 0;
  for (i = 1; ax1[i] != eslSENTINEL && ax2[i] != eslSENTINEL; i++) 
    {
      if (esl_abc_XIsBasic(abc, ax1[i])) len1++;
      if (esl_abc_XIsBasic(abc, ax2[i])) len2++;

      if (esl_abc_XIsBasic(abc, ax1[i]) && esl_abc_XIsBasic(abc, ax2[i])
	  && ax1[i] == ax2[i])
	idents++;
    }
  if (len2 < len1) len1 = len2;

  if (ax1[i] != eslSENTINEL || ax2[i] != eslSENTINEL) 
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


/* Function:  esl_dst_CPairIdMx()
 * Incept:    SRE, Thu Apr 27 08:46:08 2006 [New York]
 *
 * Purpose:   Given a multiple sequence alignment <as>, consisting
 *            of <N> aligned character strings, and optionally
 *            a bioalphabet <abc> (or NULL); calculate 
 *            a symmetric pairwise identity matrix by $N(N-1)/2$
 *            calls to <esl_dst_CPairId()>, and return it in 
 *            <ret_D>.
 *
 * Returns:   <eslOK> on success, and <ret_S> contains the identity
 *            matrix; caller is obligated to free <S> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslECORRUPT> if a seq has an illegal nonsequence char
 *            according to <abc>; <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CPairIdMx(ESL_ALPHABET *abc, char **as, int N, ESL_DMATRIX **ret_S)
{
  int status;
  ESL_DMATRIX *S = NULL;
  int i,j;

  if (( S = esl_dmatrix_Create(N,N) ) == NULL) goto FAILURE;
  
  for (i = 0; i < N; i++)
    {
      S->mx[i][i] = 1.;
      for (j = i+1; j < N; j++)
	{
	  status = esl_dst_CPairId(abc, as[i], as[j], &(S->mx[i][j]), NULL, NULL);
	  if (status != eslOK)
	    ESL_FAIL(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
	  S->mx[j][i] =  S->mx[i][j];
	}
    }

  *ret_S = S;
  return eslOK;

 FAILURE:
  if (S != NULL) esl_dmatrix_Destroy(S);
  *ret_S = NULL;
  return status;
}

/* Function:  esl_dst_XPairIdMx()
 * Incept:    SRE, Thu Apr 27 09:08:11 2006 [New York]
 *
 * Purpose:   Given a digitized multiple sequence alignment <ax>, consisting
 *            of <N> aligned digital sequences in alphabet <abc>; calculate
 *            a symmetric pairwise identity matrix by $N(N-1)/2$
 *            calls to <esl_dst_XPairId()>, and return it in <ret_S>.
 *
 * Returns:   <eslOK> on success, and <ret_S> contains the distance
 *            matrix; caller is obligated to free <S> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_S> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XPairIdMx(ESL_ALPHABET *abc, char **ax, int N, ESL_DMATRIX **ret_S)
{
  int status;
  ESL_DMATRIX *S = NULL;
  int i,j;

  if (( S = esl_dmatrix_Create(N,N) ) == NULL) goto FAILURE;
  
  for (i = 0; i < N; i++)
    {
      S->mx[i][i] = 1.;
      for (j = i+1; j < N; j++)
	{
	  status = esl_dst_XPairId(abc, ax[i], ax[j], &(S->mx[i][j]), NULL, NULL);
	  if (status != eslOK)
	    ESL_FAIL(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
	  S->mx[j][i] =  S->mx[i][j];
	}
    }
  *ret_S = S;
  return eslOK;

 FAILURE:
  if (S != NULL) esl_dmatrix_Destroy(S);
  *ret_S = NULL;
  return status;
}


/* Function:  esl_dst_CDiffMx()
 * Incept:    SRE, Fri Apr 28 06:27:20 2006 [New York]
 *
 * Purpose:   Same as <esl_dst_CPairIdMx()>, but calculates
 *            the fractional difference <d=1-s> instead of the
 *            fractional identity <s> for each pair.
 *
 * Returns:   <eslOK> on success, and <ret_D> contains the
 *            difference matrix. Caller free's <D> with 
 *            <esl_dmatrix_Destroy()>.
 *
 * Throws:    <eslECORRUPT> if one or more seqs have illegal
 *            nonsequence chars according to <abc>.
 *            <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CDiffMx(ESL_ALPHABET *abc, char **as, int N, ESL_DMATRIX **ret_D)
{
  int status;
  ESL_DMATRIX *D = NULL;
  int i,j;

  status = esl_dst_CPairIdMx(abc, as, N, &D);
  if (status != eslOK) goto FAILURE;

  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      for (j = i+1; j < N; j++) 
	{
	  D->mx[i][j] = 1. - D->mx[i][j];
	  D->mx[j][i] = D->mx[i][j];
	}
    }
  *ret_D = D;
  return eslOK;

 FAILURE:
  if (D != NULL) esl_dmatrix_Destroy(D);
  *ret_D = NULL;
  return status;

}

/* Function:  esl_dst_XDiffMx()
 * Incept:    SRE, Fri Apr 28 06:37:29 2006 [New York]
 *
 * Purpose:   Same as <esl_dst_XPairIdMx()>, but calculates fractional
 *            difference <1-s> instead of fractional identity <s> for
 *            each pair.
 *
 * Returns:   <eslOK> on success, and <ret_D> contains the difference
 *            matrix; caller is obligated to free <D> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XDiffMx(ESL_ALPHABET *abc, char **ax, int N, ESL_DMATRIX **ret_D)
{
  int status;
  ESL_DMATRIX *D = NULL;
  int i,j;

  status = esl_dst_XPairIdMx(abc, ax, N, &D);
  if (status != eslOK) goto FAILURE;

  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      for (j = i+1; j < N; j++) 
	{
	  D->mx[i][j] = 1. - D->mx[i][j];
	  D->mx[j][i] = D->mx[i][j];
	}
    }
  *ret_D = D;
  return eslOK;

 FAILURE:
  if (D != NULL) esl_dmatrix_Destroy(D);
  *ret_D = NULL;
  return status;
}


/* jukescantor()
 * 
 * The generalized Jukes/Cantor distance calculation.
 * Given <n1> identities and <n2> differences, for a
 * base alphabet size of <alphabet_size> (4 or 20);
 * calculate J/C distance in substitutions/site and
 * return it in <ret_distance>; calculate large-sample
 * variance and return it in <ret_variance>.
 *
 * Returns <eslEDIVZERO> if there are no data (<n1+n2=0>).
 */
static int
jukescantor(int n1, int n2, int alphabet_size, double *ret_distance, double *ret_variance)
{
  int    status;
  double D, K, N;
  double x;
  double distance, variance;

  ESL_DASSERT1( (n1 >= 0) );
  ESL_DASSERT1( (n2 >= 0) );
  ESL_DASSERT1( (alphabet_size >= 0) );

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
  if (ret_distance != NULL)  *ret_distance = distance;
  if (ret_variance != NULL)  *ret_variance = variance;
  return eslOK;

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
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
 *            The alphabet <abc> is required: we must know the alphabet
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
 *            Infinite distances are possible, in which case distance
 *            and variance are both <HUGE_VAL>. Caller has to deal
 *            with this case as it sees fit, perhaps by enforcing
 *            an arbitrary maximum distance.
 *
 * Throws:    <eslEINVAL> if the two strings aren't the same length (and
 *            thus can't have been properly aligned).
 *            <eslEDIVZERO> if no aligned residues were counted.
 *            On either failure, distance and variance are both returned
 *            as <HUGE_VAL>.
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
      if (esl_abc_XIsBasic(abc, x1) && esl_abc_XIsBasic(abc, x2))
	{
	  if (x1 == x2) n1++;
	  else          n2++;
	}
    }
  if (as1[i] != '\0' || as2[i] != '\0') 
    ESL_FAIL(eslEINVAL, "strings not same length, not aligned");
  
  return jukescantor(n1, n2, abc->K, ret_distance, ret_variance);

 FAILURE:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
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
  for (i = 1; ax[i] != eslSENTINEL && ay[i] != eslSENTINEL; i++) 
    {
      if (esl_abc_XIsBasic(abc, ax[i]) && esl_abc_XIsBasic(abc, ay[i]))
	{
	  if (ax[i] == ay[i]) n1++;
	  else                n2++;
	}
    }
  if (ax[i] != eslSENTINEL || ay[i] != eslSENTINEL) 
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
 * Args:      abc    - bioalphabet for <aseq>
 *            aseq   - aligned sequences [0.nseq-1][0..L-1]
 *            nseq   - number of aseqs
 *            ret_D  - RETURN: [0..nseq-1]x[0..nseq-1] symmetric distance mx
 *            ret_V  - optRETURN: matrix of variances.
 *
 * Returns:   <eslOK> on success. <D> (and optionally <V>) contain the
 *            distance matrix (and variances); caller frees these with
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if any pair of sequences have differing lengths
 *            (and thus cannot have been properly aligned). 
 *            <eslEDIVZERO> if some pair of sequences had no aligned
 *            residues. On failure, <D> and <V> are both returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CJukesCantorMx(ESL_ALPHABET *abc, char **aseq, int nseq, 
		       ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V)
{
  int          status;
  ESL_DMATRIX *D = NULL;
  ESL_DMATRIX *V = NULL;
  int          i,j;

  if (( D = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto FAILURE;
  if (( V = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto FAILURE;

  for (i = 0; i < nseq; i++)
    {
      D->mx[i][i] = 0.;
      V->mx[i][i] = 0.;
      for (j = i+1; j < nseq; j++)
	{
	  status = esl_dst_CJukesCantor(abc, aseq[i], aseq[j], 
					&(D->mx[i][j]), &(V->mx[i][j]));
	  if (status != eslOK) 
	    ESL_FAIL(status, "J/C calculation failed at seqs %d,%d", i,j);

	  D->mx[j][i] = D->mx[i][j];
	  V->mx[j][i] = V->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D;  else esl_dmatrix_Destroy(D);
  if (ret_V != NULL) *ret_V = V;  else esl_dmatrix_Destroy(V);
  return eslOK;

 FAILURE:
  if (D     != NULL) esl_dmatrix_Destroy(D);
  if (V     != NULL) esl_dmatrix_Destroy(V);
  if (ret_D != NULL) *ret_D = NULL;
  if (ret_V != NULL) *ret_V = NULL;
  return status;
}


/* Function:  esl_dst_XJukesCantorMx()
 * Incept:    SRE, Thu Apr 27 08:38:08 2006 [New York City]
 *
 * Purpose:   Given a digitized multiple sequence alignment <ax>,
 *            consisting of <nseq> aligned digital sequences in
 *            bioalphabet <abc>, calculate a symmetric Jukes/Cantor
 *            pairwise distance matrix for all sequence pairs, and
 *            return the distance matrix in <ret_D>. Optionally also
 *            return the large-sample variances for those ML distance
 *            estimates in <ret_V>.
 *            
 *            Infinite distances (and variances) are possible; they
 *            are represented as <HUGE_VAL> in <D> and <V>. Caller must
 *            be prepared to deal with them as appropriate.
 *
 * Args:      abc    - bioalphabet for <aseq>
 *            ax     - aligned digital sequences [0.nseq-1][1..L]
 *            nseq   - number of aseqs
 *            ret_D  - RETURN: [0..nseq-1]x[0..nseq-1] symmetric distance mx
 *            ret_V  - optRETURN: matrix of variances.
 *
 * Returns:   <eslOK> on success. <D> (and optionally <V>) contain the
 *            distance matrix (and variances); caller frees these with
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if any pair of sequences have differing lengths
 *            (and thus cannot have been properly aligned). 
 *            <eslEDIVZERO> if some pair of sequences had no aligned
 *            residues. On failure, <D> and <V> are both returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XJukesCantorMx(ESL_ALPHABET *abc, char **ax, int nseq, 
		       ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V)
{
  int          status;
  ESL_DMATRIX *D = NULL;
  ESL_DMATRIX *V = NULL;
  int          i,j;

  if (( D = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto FAILURE;
  if (( V = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto FAILURE;

  for (i = 0; i < nseq; i++)
    {
      D->mx[i][i] = 0.;
      V->mx[i][i] = 0.;
      for (j = i+1; j < nseq; j++)
	{
	  status = esl_dst_XJukesCantor(abc, ax[i], ax[j], 
					&(D->mx[i][j]), &(V->mx[i][j]));
	  if (status != eslOK) 
	    ESL_FAIL(status, "J/C calculation failed at digital aseqs %d,%d", i,j);

	  D->mx[j][i] = D->mx[i][j];
	  V->mx[j][i] = V->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D;  else esl_dmatrix_Destroy(D);
  if (ret_V != NULL) *ret_V = V;  else esl_dmatrix_Destroy(V);
  return eslOK;

 FAILURE:
  if (D     != NULL) esl_dmatrix_Destroy(D);
  if (V     != NULL) esl_dmatrix_Destroy(V);
  if (ret_D != NULL) *ret_D = NULL;
  if (ret_V != NULL) *ret_V = NULL;
  return status;
}

/*****************************************************************
 * Example.
 *****************************************************************/ 
/* 
   gcc -g -Wall -o example -I. -DeslDISTANCE_EXAMPLE esl_distance.c\
        esl_alphabet.c esl_dmatrix.c esl_msa.c easel.c -lm
 */
#ifdef eslDISTANCE_EXAMPLE
/*::cexcerpt::distance_example::begin::*/
#include <easel.h>
#include <esl_distance.h>
#include <esl_alphabet.h>
#include <esl_msa.h>
int main(int argc, char **argv)
{
  char        *filename;
  int          fmt;
  ESL_MSAFILE *afp; 
  ESL_MSA     *msa;
  ESL_DMATRIX *P;
  int          status;
  int          i,j;
  double       min, avg, max;

  filename = argv[1];
  fmt      = eslMSAFILE_UNKNOWN;
  
  if ( (status = esl_msafile_Open(filename, fmt, NULL, &afp)) != eslOK) 
    esl_fatal("Alignment file open failed with error %d\n", status);
  if ( (status = esl_msa_Read(afp, &msa)) != eslOK)
    esl_fatal("Failed to read alignment from %s\n", filename);

  esl_dst_CPairIdMx(NULL, msa->aseq, msa->nseq, &P);

  min = 1.0;
  max = 0.0;
  avg = 0.0;
  for (i = 0; i < msa->nseq; i++)
    for (j = i+1; j < msa->nseq; j++)
      {
	avg += P->mx[i][j];
	if (P->mx[i][j] < min) min = P->mx[i][j];
	if (P->mx[i][j] > max) max = P->mx[i][j];
      }
  avg /= (double) (msa->nseq * (msa->nseq-1) / 2);

  printf("Average pairwise %% id:  %.1f%%\n", avg * 100.);
  printf("Minimum pairwise %% id:  %.1f%%\n", min * 100.);
  printf("Maximum pairwise %% id:  %.1f%%\n", max * 100.);

  esl_dmatrix_Destroy(P);
  esl_msa_Destroy(msa);
  esl_msafile_Close(afp);
  return 0;
}
/*::cexcerpt::distance_example::end::*/
#endif /*eslDISTANCE_EXAMPLE*/



/*****************************************************************
 * Test driver.
 *****************************************************************/ 
/* 
   gcc -g -Wall -o testdriver -I. -L. -DeslDISTANCE_TESTDRIVE esl_distance.c -leasel -lm
 */
#ifdef eslDISTANCE_TESTDRIVE
#include <easel.h>
#include <esl_getopts.h>
#include <esl_distance.h>
#include <esl_random.h>
#include <esl_alphabet.h>
#include <esl_vectorops.h>

static ESL_OPTIONS options[] = {
  /* name        type       def   env  range toggles reqs incomp help                       docgroup*/
  { "-h",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",            0},
  { "-N",     eslARG_INT,     "6", NULL,"n>2", NULL, NULL, NULL, "number of iid seqs in alignment",0},
  { "-L",     eslARG_INT,    "50", NULL,"n>0", NULL, NULL, NULL, "length of seqs in alignment",    0},
  { "--seed", eslARG_INT,    "42", NULL,"n>0", NULL, NULL, NULL, "random # seed",                  0},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[] = "Usage: ./testdrive-distance [-options]";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go = NULL;
  ESL_ALPHABET *abc = NULL;
  ESL_RANDOMNESS *r = NULL;
  ESL_DMATRIX *S = NULL,
              *D = NULL,
              *V = NULL;
  int status;
  int N,L;			/* # of seqs, and their aligned lengths */
  int show_help;
  int i;
  int j;
  int seed;
  char **as = NULL;		/* aligned character seqs (random, iid) */
  char **ax = NULL;		/* digitized alignment                  */
  double p[4];			/* ACGT probabilities */
  double pid;
  int    nid;
  int    nres;
  double distance, variance;

  /* Process command line
   */
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  esl_opt_GetBooleanOption(go, "-h", &show_help);
  if (show_help) {
    puts(usage); 
    puts("\n  where options are:");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2=indentation; 80=width */
    return 0;
  }
  esl_opt_GetIntegerOption(go, "-L",     &L);
  esl_opt_GetIntegerOption(go, "-N",     &N);
  esl_opt_GetIntegerOption(go, "--seed", &seed);
  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return 1;
  }
  esl_getopts_Destroy(go);

  /* Create a random DNA alignment and digitize it.
   */
  r   = esl_randomness_Create(seed);
  abc = esl_alphabet_Create(eslDNA);

  esl_vec_DSet(p, 4, 0.25);
  ESL_ALLOC(as, sizeof(char *) * N);
  ESL_ALLOC(ax, sizeof(char *) * N);
  for (i = 0; i < N; i++) as[i] = ax[i] = NULL;
  for (i = 0; i < N; i++) esl_rnd_IID(r, "ACGT", p, 4, L, &(as[i]));
  for (i = 0; i < N; i++) esl_dsq_Create(abc, as[i], L, &(ax[i]));

  /* CPairId() tests: self-comparisons give 1
   */
  esl_dst_CPairId(abc, as[0], as[0], &pid, &nid, &nres);
  if (pid  != 1.0) ESL_FAIL(eslFAIL, "pid != 1.0 in self comparison, CPairId");
  if (nid  != L)   ESL_FAIL(eslFAIL, "nid != L in self comparison, CPairId");
  if (nres != L)   ESL_FAIL(eslFAIL, "nres != L in self comparison, CPairId");
  esl_dst_CPairId(NULL, as[0], as[0], &pid, &nid, &nres);  /* w/o alphabet */
  if (pid  != 1.0 || nid != L || nres != L) ESL_FAIL(eslFAIL, "problem w/ no-alphabet CPairId");
  esl_dst_CPairId(abc, as[0], as[0], &pid, NULL, NULL); 
  if (pid != 1.0) ESL_FAIL(eslFAIL, "CPairId() optional args test");

  /* XPairId() self tests
   */
  esl_dst_XPairId(abc, ax[0], ax[0], &pid, &nid, &nres);
  if (pid  != 1.0) ESL_FAIL(eslFAIL, "pid != 1.0 in self comparison, XPairId");
  if (nid  != L)   ESL_FAIL(eslFAIL, "nid != L in self comparison, XPairId");
  if (nres != L)   ESL_FAIL(eslFAIL, "nres != L in self comparison, XPairId");
  esl_dst_XPairId(abc, ax[0], ax[0], &pid, NULL, NULL); 
  if (pid != 1.0) ESL_FAIL(eslFAIL, "XPairId() optional args test");

  /* CPairIdMx() tests: on random seqs, we expect avg off diag to be 
   * about 25%, and diag to be all 100%
   */
  esl_dst_CPairIdMx(abc, as, N, &S);
  for (i = 0; i < N; i++) 
    if (S->mx[i][i] != 1.0) ESL_FAIL(eslFAIL, "CPairIdMx: non-1.0 diagonal");
  pid = 0.;
  for (i = 0; i < N; i++)
    for (j = i+1; j < N; j++)
      pid += S->mx[i][j];
  pid /= (double) (N * (N-1) / 2);
  if (pid < 0.15 || pid > 0.35) ESL_FAIL(eslFAIL, "CPairIdMx: offdiags aren't ~25%%");
  esl_dmatrix_Destroy(S); S=NULL;

  /* XPairIdMx() tests
   */
  esl_dst_XPairIdMx(abc, ax, N, &S);
  for (i = 0; i < N; i++) 
    if (S->mx[i][i] != 1.0) ESL_FAIL(eslFAIL, "XPairIdMx: non-1.0 diagonal");
  pid = 0.;
  for (i = 0; i < N; i++)
    for (j = i+1; j < N; j++)
      pid += S->mx[i][j];
  pid /= (double) (N * (N-1) / 2);
  if (pid < 0.15 || pid > 0.35) ESL_FAIL(eslFAIL, "XPairIdMx: offdiags aren't ~25%%");
  esl_dmatrix_Destroy(S); S=NULL;

  /* CDiffMx() tests
   */
  esl_dst_CDiffMx(abc, as, N, &S);
  for (i = 0; i < N; i++) 
    if (S->mx[i][i] != 0.0) ESL_FAIL(eslFAIL, "CDiffMx: non-0 diagonal");
  pid = 0.;
  for (i = 0; i < N; i++)
    for (j = i+1; j < N; j++)
      pid += S->mx[i][j];
  pid /= (double) (N * (N-1) / 2);
  if (pid < 0.65 || pid > 0.85) ESL_FAIL(eslFAIL, "CDiffMx: offdiags aren't ~75%%");
  esl_dmatrix_Destroy(S); S=NULL;
  
  /* XDiffMx() tests
   */
  esl_dst_XDiffMx(abc, ax, N, &S);
  for (i = 0; i < N; i++) 
    if (S->mx[i][i] != 0.0) ESL_FAIL(eslFAIL, "XDiffMx: non-0 diagonal");
  pid = 0.;
  for (i = 0; i < N; i++)
    for (j = i+1; j < N; j++)
      pid += S->mx[i][j];
  pid /= (double) (N * (N-1) / 2);
  if (pid < 0.65 || pid > 0.85) ESL_FAIL(eslFAIL, "XDiffMx: offdiags aren't ~75%%");
  esl_dmatrix_Destroy(S); S=NULL;
      

  /* CJukesCantor() self-comparison tests
   */
  esl_dst_CJukesCantor(abc, as[0], as[0], &distance, &variance);
  if (distance  != 0.0) ESL_FAIL(eslFAIL, "distance != 0.0 in self comparison, CJukesCantor()");
  esl_dst_CJukesCantor(abc, as[0], as[0], &distance, NULL);
  if (distance != 0.0) ESL_FAIL(eslFAIL, "CJukesCantor() optional args test");

  /* XJukesCantor() self-comparison tests
   */
  esl_dst_XJukesCantor(abc, ax[0], ax[0], &distance, &variance);
  if (distance  != 0.0) ESL_FAIL(eslFAIL, "distance != 0.0 in self comparison, XJukesCantor()");
  esl_dst_XJukesCantor(abc, ax[0], ax[0], &distance, NULL);
  if (distance != 0.0) ESL_FAIL(eslFAIL, "XJukesCantor() optional args test");

  /* Jukes-Cantor distance matrices - just test for crashing.
   */
  esl_dst_CJukesCantorMx(abc, as, N, &D, &V);
  esl_dmatrix_Destroy(D); D= NULL;
  esl_dmatrix_Destroy(V); V= NULL;
  esl_dst_CJukesCantorMx(abc, as, N, &D, NULL);
  esl_dmatrix_Destroy(D); D= NULL;
  esl_dst_XJukesCantorMx(abc, ax, N, &D, &V);
  esl_dmatrix_Destroy(D); D= NULL;
  esl_dmatrix_Destroy(V); V= NULL;
  esl_dst_XJukesCantorMx(abc, ax, N, &D, NULL);
  esl_dmatrix_Destroy(D); D= NULL;


  status = eslOK;
  /* ok to flow through below: cleanup only */
FAILURE:
  if (r   != NULL) esl_randomness_Destroy(r);
  if (abc != NULL) esl_alphabet_Destroy(abc);
  if (S   != NULL) esl_dmatrix_Destroy(S);
  if (D   != NULL) esl_dmatrix_Destroy(D);
  if (V   != NULL) esl_dmatrix_Destroy(V);
  esl_Free2D((void **) as, N);
  esl_Free2D((void **) ax, N);
  return status;
}
#endif /*eslDISTANCE_TESTDRIVE*/

 /*****************************************************************
  * @LICENSE@
  *****************************************************************/

