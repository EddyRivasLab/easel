/* esl_wuss : RNA secondary structure markup annotation format
 *
 * Full (Infernal output) WUSS notation:
 *    <>  : base pairs in stem-loops, i.e. <<<___>>>
 *    ()  : base pairs in helices enclosing multifurcation of all <> stems
 *    []  : base pairs in helices enclosing multifurc of at least one () helix
 *    {}  : base pairs in helices enclosing even deeper multifurcations
 *    
 *    _   : (i.e. underscore) nucleotide in a hairpin loop
 *    -   : (i.e. dash) nucleotide in a bulge or interior loop
 *    ,   : nucleotide in a multifurcation loop (mnemonic: "stem1, stem2," )  
 *    :   : nucleotide in external single strand
 *    
 *    Aa  : (and Bb, Cc, etc) pseudoknotted base pairs, upper case on left, lower case on right.
 *
 *  and in alignments of a seq to an RNA structure profile, as in Infernal:
 *    .   : insertion relative to a known consensus structure
 *    ~   : nucleotide is unaligned to a structure profile, because of local structure alignment
 *
 * "Simple" (input) WUSS:
 *  ..<<..AA..>>..aa..  dots and brackets, plus Aa etc for PKs.
 *
 *
 * Contents:
 *    1. Functions that allow annotation to be a string or esl_mem.
 *    2. Functions that assume annotation is a string.
 *    3. Interconversion with obsolete Konings-Hogeweg convention.
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 *    
 */
#include <esl_config.h>

#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_stack.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

/*****************************************************************
 * 1. Functions that allow annotation to be a string or esl_mem.
 *****************************************************************/

/* Function:  esl_wuss_IsSimple()
 * Synopsis:  Returns TRUE if string looks like simple WUSS annotation.
 * Incept:    SRE, Fri 02 Jan 2026
 *
 * Purpose:   Return TRUE if input string or buf <ss> of length <n> is
 *            consistent with simple WUSS secondary structure
 *            annotation format. Else, return FALSE.
 *
 *            If <ss> is an Easel parsing buffer (no \0 terminator), length
 *            <n> must be provided. If <ss> is a string (\0-terminated)
 *            and its length is unknown, <n> can be passed as -1.
 *
 *            The string is considered consistent if it is composed
 *            only of '.', one choice of open/close parenthesis-like
 *            characters "<>" or "()" or "[]" or "{}" (using more than
 *            one choice is not simple format), or alphabetic
 *            characters [A-Za-z] (for pseudoknots). Moreover, the
 *            pairs indicated by the annotation must be validly
 *            balanced: the number of open-pair characters must equal
 *            the number of close-pair characters.
 *
 *            Because this test will tend to fail quickly on a full
 *            WUSS annotation string (often on the first character),
 *            there is negligible cost in evaluating the ternary
 *            simple|full|not WUSS at all by cascading if's:
 * 
 *               if      (esl_wuss_IsSimple(ss)) {...}
 *               else if (esl_wuss_IsFull(ss)    {...}
 *               else                            {...}
 */
int
esl_wuss_IsSimple(const char *ss, int n)
{
  int has_angle = 0;
  int has_paren = 0;
  int has_brack = 0;
  int has_curly = 0;
  int nopen     = 0;
  int nclose    = 0;
  int nuppers[26];
  int nlowers[26];
  int i;

  esl_vec_ISet(nuppers, 26, 0);
  esl_vec_ISet(nlowers, 26, 0);

  if (n < 0) n = strlen(ss);

  for (i = 0; i < n; i++)
    {
      if      (strchr("<>", ss[i])) has_angle = 1;
      else if (strchr("()", ss[i])) has_paren = 1;
      else if (strchr("[]", ss[i])) has_brack = 1;
      else if (strchr("{}", ss[i])) has_curly = 1;
      else if (ss[i] != '.' && !isalpha(ss[i])) return FALSE;

      if       (strchr("<([{", ss[i])) nopen++;
      else if  (strchr(">)]}", ss[i])) nclose++;
      else if  (isupper(ss[i]))        nuppers[ss[i]-'A']++;
      else if  (islower(ss[i]))        nlowers[ss[i]-'a']++;
    }

  if (has_angle + has_paren + has_brack + has_curly > 1) return FALSE;
  if (nopen != nclose) return FALSE;
  for (i = 0; i < 26; i++) if (nuppers[i] != nlowers[i]) return FALSE;

  return TRUE;
}

/* Function:  esl_wuss_IsFull()
 * Synopsis:  Returns TRUE if string looks like full WUSS annotation.
 * Incept:    SRE, Fri 02 Jan 2026
 *
 * Purpose:   Return TRUE if input string <ss> is consistent with full
 *            WUSS secondary structure annotation format. Else, return
 *            FALSE.
 *
 *            If <ss> is an Easel parsing buffer (no \0 terminator), length
 *            <n> must be provided. If <ss> is a string (\0-terminated)
 *            and its length is unknown, <n> can be passed as -1.
 *
 *
 *            The string is considered consistent if it is composed
 *            only of pairing chars "<>()[]{}", unpaired chars "_-,:",
 *            gap pad chars ".~", or pseudoknot pairing chars
 *            [A-Za-z].  Moreover, annotated pairs must be validly
 *            balanced: the number of open-pair characters must equal
 *            the number of close-pair characters, for each type of
 *            open char "<([{" and [A-Z], and each type of close char
 *            ">)]}" and [a-z].
 *
 *            Simple WUSS annotation is valid as full WUSS, so if you're
 *            aiming to distinguish simple from full, test simple first:
 *            e.g.:
 *
 *               if      (esl_wuss_IsSimple(ss)) {...}
 *               else if (esl_wuss_IsFull(ss)    {...}
 *               else                            {...}
 */
int
esl_wuss_IsFull(const char *ss, int n)
{
  int nl_angle = 0;
  int nl_paren = 0;
  int nl_brack = 0;
  int nl_curly = 0;
  int nr_angle = 0;
  int nr_paren = 0;
  int nr_brack = 0;
  int nr_curly = 0;
  int nuppers[26];
  int nlowers[26];
  int i;

  esl_vec_ISet(nuppers, 26, 0);
  esl_vec_ISet(nlowers, 26, 0);

  if (n < 0) n = strlen(ss);
  
  for (i = 0; i < n; i++)
    {
      switch (ss[i]) {
      case '<': nl_angle++; break;
      case '(': nl_paren++; break;
      case '[': nl_brack++; break;
      case '{': nl_curly++; break;
      case '>': nr_angle++; break;
      case ')': nr_paren++; break;
      case ']': nr_brack++; break;
      case '}': nr_curly++; break;
      default:
        if      (isupper(ss[i])) nuppers[ss[i]-'A']++;
        else if (islower(ss[i])) nlowers[ss[i]-'a']++;
        else if (!strchr("_-,:.~", ss[i])) return FALSE;
      }
    }

  if (nl_angle != nr_angle ||
      nl_paren != nr_paren ||
      nl_brack != nr_brack ||
      nl_curly != nr_curly)
    return FALSE;
  for (i = 0; i < 26; i++)
    if (nuppers[i] != nlowers[i])
      return FALSE;

  return TRUE;
}



/* Function:  esl_wuss2ct()
 * Incept:    SRE, Tue Feb 15 08:44:54 2005 [St. Louis]
 *
 * Purpose:   Given a secondary structure string <ss>, <0..len-1>,
 *            in WUSS notation, convert it to a CT array, <1..len>,
 *            in <ct>. Caller provides a <ct> allocated for at least 
 *            <len+1> ints. <ct[i]> is the position that residue i
 *            base pairs to, or 0 if i is unpaired. <ct[0]> is undefined
 *            (but if you care: it is set to 0).
 *            
 *            WUSS notation is interpreted loosely here, as input
 *            WUSS.  Any matching bracket pair or upper/lower case
 *            alphabetic pair is interpreted as a base pair; any other
 *            WUSS annotation is interpreted as unpaired.
 *            
 * Returns:   <eslOK> on success.
 *            <eslESYNTAX> if the WUSS string isn't valid.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 */
int 
esl_wuss2ct(const char *ss, int len, int *ct)
{
  ESL_STACK *pda[27];     /* 1 secondary structure + up to 26 levels of pk's */
  int        i;
  int        pos, pair;
  int        status;      /* success or failure return status */

 /* Initialization: always initialize the main pda (0);
  * we'll init the pk pda's on demand.
  */
  for (i = 1; i <= 26; i++) pda[i] = NULL;
  if ((pda[0] = esl_stack_ICreate()) == NULL) goto FINISH;

  for (pos = 0; pos <= len; pos++) ct[pos] = 0;

  for (pos = 1; pos <= len; pos++)
    {
      if (!isprint((int) ss[pos-1]))  /* armor against garbage */
	{ status = eslESYNTAX; goto FINISH; }

      /* left side of a pair: push position onto stack 0 (pos = 1..L) */
      else if (ss[pos-1] == '<' ||
	       ss[pos-1] == '(' ||
	       ss[pos-1] == '[' ||
	       ss[pos-1] == '{')
	{
	  if ((status = esl_stack_IPush(pda[0], pos)) != eslOK) goto FINISH;
	}
      
      /* right side of a pair; resolve pair; check for agreement */
      else if (ss[pos-1] == '>' || 
	       ss[pos-1] == ')' ||
	       ss[pos-1] == ']' ||
	       ss[pos-1] == '}')
        {
          if (esl_stack_IPop(pda[0], &pair) == eslEOD)
            { status = eslESYNTAX; goto FINISH;} /* no closing bracket */
          else if ((ss[pair-1] == '<' && ss[pos-1] != '>') ||
		   (ss[pair-1] == '(' && ss[pos-1] != ')') ||
		   (ss[pair-1] == '[' && ss[pos-1] != ']') ||
		   (ss[pair-1] == '{' && ss[pos-1] != '}'))
	    { status = eslESYNTAX; goto FINISH; }  /* brackets don't match */
	  else
	    {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
                                /* same stuff for pseudoknots */
      else if (isupper((int) ss[pos-1])) 
	{
	  /* Create the PK stacks on demand.
	   */
	  i = ss[pos-1] - 'A' + 1;
	  if (pda[i] == NULL) 
	    if ((pda[i] = esl_stack_ICreate()) == NULL) 
	      { status = eslEMEM; goto FINISH; }

	  if ((status = esl_stack_IPush(pda[i], pos)) != eslOK) goto FINISH;
	}
      else if (islower((int) ss[pos-1])) 
	{
	  i = ss[pos-1] - 'a' + 1;
	  if (pda[i] == NULL || 
	      esl_stack_IPop(pda[i], &pair) == eslEOD)
            { status = eslESYNTAX; goto FINISH;}
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
	}
      else if (strchr(":,_-.~", ss[pos-1]) == NULL)
	{ status = eslESYNTAX; goto FINISH; } /* bogus character */
    }
  status = eslOK;

 FINISH:
  for (i = 0; i <= 26; i++)
    if (pda[i] != NULL) 
      { /* nothing should be left on stacks */
	if (esl_stack_ObjectCount(pda[i]) != 0)
	  status = eslESYNTAX;
	esl_stack_Destroy(pda[i]);
      }
  return status;
}


/* Function:  esl_ct2wuss()
 * Incept:    SRE, Wed Feb 16 11:22:53 2005 [St. Louis]
 *            pseudoknot support by ER, Sat Aug 18 13:22:03 EDT 2012
 *
 * Purpose:   Convert a CT array <ct> for <n> residues (1..n) to a WUSS
 *            format string <ss>. <ss> must be allocated for at least
 *            n+1 chars (+1 for the terminal '\0').
 *
 *            Pseudoknots are annotated as AA...aa, BB...bb, ...,
 *            ZZ..zz, with upper case letters for open-pair and lower
 *            case for close-pair.
 *
 *            Because pseudoknotted structures can be resolved
 *            different ways into what's called "main structure" and
 *            what's called pseudoknot, if you convert a WUSS
 *            annotation for a pseudoknotted structure to CT and then
 *            back to WUSS, you may not get the same WUSS string, even
 *            though the set of base pairs is identical.
 *            
 *            The WUSS annotation generated here is for an
 *            unaligned/ungapped single sequence or consensus
 *            structure of length <n>. If the CT was made from a WUSS
 *            annotation on an Infernal alignment of alen <n> that
 *            contains `.` and `~` characters for insertion padding
 *            and local structural alignment padding, the WUSS
 *            generated here will not have any of that padding. There
 *            are no `.` or `~` characters in the generated WUSS.
 *
 *            Attempting to convert a <ct> that requires more letters
 *            than [A-Z] will return an <eslEINVAL> error.
 *
 *            Attempting to convert a <ct> that involves triplet interactions
 *            will also return an <eslEINVAL> error.
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> if the CT implies a structure that can't be rendered in WUSS.
 * 
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on "impossible" internal failures.
 */
int
esl_ct2wuss(const int *ct, int n, char *ss)
{
  int        rb[26];                // array that delimits the right bound of a pseudoknot character
  ESL_STACK *pda    = NULL;         // stack for "main" secondary structure
  ESL_STACK *auxpk  = NULL;	    // aux stack for pseudoknot
  ESL_STACK *auxss  = NULL;	    // aux stack for single stranded
  int       *cct    = NULL;         // copy of ct vector
  int        nfaces;                // number of faces in a cWW structure
  int        minface;               // max depth of faces in a cWW structure
  int        leftbound, rightbound; // left and right bound to find basepairs belonging to a given pseudoknot
  int        xpk = 0;               // number of pseudoknot characters used
  int        npairs = 0;            // total number of basepairs
  int        npairs_reached = 0;    // number of basepairs found so far
  int        found_partner;         // true if we've found left partner of a given base in stack pda
  int        i,j,k;                 // sequence indices
  int        x;                     // index for pseudoknot characters
  int        status = eslEMEM;	    // exit status 'til proven otherwise

  /* total number of basepairs */
  for (j = 1; j <= n; j ++) { if (ct[j] > 0 && j < ct[j]) npairs ++; }
  
  /* Copy of ct; if a pseudoknotted structure, cct will be modified later.
   */
  ESL_ALLOC(cct, sizeof(int)*(n+1));
  esl_vec_ICopy(ct, (n+1), cct);
  
  /* Initialize rightbounds for all 26 pseudoknot indices */
  for (x = 0; x < 26; x ++) rb[x] = -1;

  /* init ss[] to single stranded */
  for (j = 0; j < n; j ++) { ss[j] = ':'; }  
  ss[n] = '\0'; 
 
  /* Initialization*/
  if ((pda   = esl_stack_ICreate()) == NULL) goto FINISH;
  if ((auxpk = esl_stack_ICreate()) == NULL) goto FINISH;
  if ((auxss = esl_stack_ICreate()) == NULL) goto FINISH;
  
  for (j = 1; j <= n; j++)
    {
      if (cct[j] == 0)	/* unpaired: push j. */
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else if (cct[j] > j) /* left side of a bp: push j. */
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else   /* right side of a bp; main routine: find the left partner */
	{
	  found_partner = FALSE;
	  /* Pop back until we find the left partner of j;
	   * In case this is not a nested structure, finding
	   * the left partner of j will require to put bases 
	   * aside into stack auxpk.
	   *
	   * After we find the left partner of j,
	   * store single stranded residues in auxss;
	   * keep track of #faces and the maximum face depth.
	   */
	  nfaces  = 0;
	  minface = -1;
	 
	  while (esl_stack_ObjectCount(pda)) 
	    {
	      if (esl_stack_IPop(pda, &i) != eslOK) goto FINISH;
	      
	      if (i < 0) 		/* a face counter */
		{
		  nfaces++;
		  if (i < minface) minface = i;
		}

	      else if (cct[i] == j)  /* we found the i,j pair. */
		{
		  found_partner = TRUE;
		  npairs_reached ++;	
		  /* Now we know i,j pair; and we know how many faces are
		   * above them; and we know the max depth of those faces.
		   * That's enough to label the pair in WUSS notation.
		   * if nfaces == 0, minface is -1; <> a closing bp of a hairpin.
		   * if nfaces == 1, inherit minface, we're continuing a stem.
		   * if nfaces > 1, bump minface in depth; we're closing a bifurc.
		   */
		  if (nfaces > 1 && minface > -4) minface--;
		  switch (minface) {
		  case -1: ss[i-1] = '<'; ss[j-1] = '>'; break;
		  case -2: ss[i-1] = '('; ss[j-1] = ')'; break;
		  case -3: ss[i-1] = '['; ss[j-1] = ']'; break;
		  case -4: ss[i-1] = '{'; ss[j-1] = '}'; break;
		  default:
		    esl_stack_Destroy(pda); esl_stack_Destroy(auxpk); esl_stack_Destroy(auxss); free(cct); 
		    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such face code");
		  }
		  if (esl_stack_IPush(pda, minface) != eslOK) goto FINISH;
		  
		  /* Now, aux contains all the unpaired residues we need to label,
		   * according to the # of faces "above" them:
		   *  nfaces = 0: hairpin loop
		   *  nfaces = 1: bulge or interior loop
		   *  nfaces > 1: multifurc
		   */
		  while (esl_stack_IPop(auxss, &i) == eslOK)
		    {
		      switch (nfaces) {
			
		      case 0:  ss[i-1] = '_'; break;
		      case 1:  ss[i-1] = '-'; break;
		      default: ss[i-1] = ','; break; /* nfaces > 1 */
		      }
		    }
		  break;
		}
	      
	      else if (cct[i] == 0) 
		{
		  /* add to auxss only if originally single stranded */
		  if (ct[i] == 0) { if (esl_stack_IPush(auxss, i) != eslOK) goto FINISH; }
		}

	      else /* cct[i]>0, != j: i is paired, but not to j: pseudoknot! */
		{
		  /* i is in the way to find j's left partner. 
		   * Move i to stack auxpk; resolve pseudoknot(s) after we've found partner for j.
		   */ 
		  if (esl_stack_IPush(auxpk, i) != eslOK) goto FINISH;
		}
	    } 
	  
	  if (!found_partner) { status = eslEINVAL; goto ERROR; }  // Cannot find left partner ct[j] of base j. Likely a triplet.
	} /* finished finding the left partner of j */
      
      /* After we've found the left partner of j, resolve pks found along the way.
       * Then, remove the pseudoknotted based from cct so we can find the rest of the structure.
       */
      if (esl_stack_ObjectCount(auxpk)) {

	/* init for first pseudoknot */
	leftbound  = cct[j];
	rightbound = leftbound + 1;
	xpk        = -1;            /* start with 'A' if possible again */

	while (esl_stack_IPop(auxpk, &i) == eslOK) {

	  for (k = rightbound-1; k > leftbound; k --) 
	    {
	      if      (cct[k] == 0)          { continue; } 
	      else if (cct[k] >  rightbound) { continue; } 
	      else if (cct[k] == i)          { break; }                  /* i continues the given pseudoknot */
	      else                           { k = leftbound; break; }   /* a new pseudoknot */		    		
	    }
	  
	  if (k == leftbound) /* a new pseudoknot */
	    {
	      // npk ++;  // if you wanted to count PK's, here's where to do it
	      xpk ++;
	      /* figure out if we can use this alphabet index, or bump it up if necessary */
	      while (i < rb[xpk]) { xpk ++; }
	      
	      leftbound  = (rightbound < cct[i])? rightbound : cct[j];
	      rightbound = cct[i];
	    }
	      
	  npairs_reached ++;
	  if (xpk+(int)('a') <= (int)('z')) {

	    /* update the rightbound of this pk index if necessary */
	    if (cct[i] > rb[xpk]) rb[xpk] = cct[i];
	    
	    /* Add pk indices for this basepair */
	    ss[i-1]      = (char)(xpk+(int)('A'));
	    ss[cct[i]-1] = (char)(xpk+(int)('a'));
	    
	    /* remove pseudoknotted pair from cct */
	    cct[i]     = 0;
	    cct[ct[i]] = 0;
	  }
	  else { status = eslEINVAL; goto ERROR; } // Ran out of letters for PK's.
	} 	
      } /* while there is something in auxpk stack */
  } /* finished loop over j: end position on seq, 1..n*/ 
  if (npairs != npairs_reached) ESL_XEXCEPTION(eslEINCONCEIVABLE, "only found %d out of %d pairs.", npairs_reached, npairs);

  status = eslOK;
 ERROR:
 FINISH:
  esl_stack_Destroy(pda);
  esl_stack_Destroy(auxpk);
  esl_stack_Destroy(auxss);
  free(cct);
  return status;
}

/* Function:  esl_ct2simplewuss()
 * Incept:    ER, Wed Aug 22 13:31:54 EDT 2012 [Janelia]
 *
 * Purpose:   Convert a CT array <ct> for <n> residues (1..n) to a simple WUSS
 *            format string <ss>. <ss> must be allocated for at least
 *            n+1 chars (+1 for the terminal NUL). 
 *
 *            This function can be used with the <ct> of a secondary
 *            structure including arbitrary pseudoknots, or for the 
 *            <ct> of a tertiary structure (say cWH, tWH, cSS,... H bonds). 
 *
 *            The string <ss> has basepairs annotated as <>, Aa, Bb, ..., Zz;
 *            unpaired bases are annotated as '.'.
 *
 *            Attempting to convert a <ct> that requires more letters
 *            than [A-Z] will return an <eslEINVAL> error.
 *
 *            Attempting to convert a <ct> that involves triplet interactions
 *            will return an <eslEINVAL> error.
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> if structure in CT cannot be represented in WUSS format.
 * 
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal failure.
 */
int
esl_ct2simplewuss(const int *ct, int n, char *ss)
{
  int        rb[26];                // array that delimits the right bound of a pseudoknot character
  ESL_STACK *pda    = NULL;         // stack for "main" secondary structure
  ESL_STACK *auxpk  = NULL;	    // aux stack for pseudoknot
  int       *cct    = NULL;         // copy of ct vector
  int        leftbound, rightbound; // left and right bound to find basepairs belonging to a given pseudoknot
  int        xpk = 0;               // number of pseudoknot characters used
  int        npairs = 0;            // total number of basepairs
  int        npairs_reached = 0;    // number of basepairs found so far
  int        found_partner;         // true if we've found left partner of a given base in stack pda
  int        i,j,k;                 // sequence indices
  int        x;                     // index for pseudoknot characters
  int        status;

  /* total number of basepairs */
  for (j = 1; j <= n; j ++) { if (ct[j] > 0 && j < ct[j]) npairs ++; }
  
  /* Copy of ct; if a pseudoknotted structure, cct will be modified later.
   */
  ESL_ALLOC(cct, sizeof(int)*(n+1));
  esl_vec_ICopy(ct, (n+1), cct);
  
  /* Initialize rightbounds for all 26 pseudoknot indices */
  for (x = 0; x < 26; x ++) rb[x] = -1;

  /* init ss[] to single stranded */
  for (j = 0; j < n; j ++) { ss[j] = '.'; }  
  ss[n] = '\0'; 
 
  /* Initialization*/
  if ((pda   = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto ERROR; }
  if ((auxpk = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto ERROR; }
  
  for (j = 1; j <= n; j++)
    {
      if (cct[j] == 0)	/* unpaired: push j. */
	{
	  if ((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	}
      else if (cct[j] > j) /* left side of a bp: push j. */
	{
	  if ((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	}
      else   /* right side of a bp; main routine: fingh the left partner */
	{
	  found_partner = FALSE;

	  /* Pop back until we find the left partner of j;
	   * In case this is not a nested structure, finding
	   * the left partner of j will require to put bases 
	   * aside into stack auxpk.
	   */	 
	  while (esl_stack_ObjectCount(pda)) 
	    {
	      if ((status = esl_stack_IPop(pda, &i)) != eslOK) goto ERROR;
	      
	      if (cct[i] == j)  /* we found the i,j pair. */
		{
		  found_partner = TRUE;
		  npairs_reached ++;	

		  ss[i-1] = '<';
		  ss[j-1] = '>';
		  break;
		}
	      
	      else if (cct[i] == 0) 
		{
		  if (ct[i] == 0) ss[i-1] = '.';
		}

	      else /* cct[i]>0, != j: i is paired, but not to j: pseudoknot! */
		{
		  /* i is in the way to find j's left partner. 
		   * Move i to stack auxpk; resolve pseudoknot(s) after we've found partern for j.
		   */ 
		  if ((status = esl_stack_IPush(auxpk, i)) != eslOK) goto ERROR;
		}
	    } 
	  if (!found_partner) { status = eslEINVAL; goto ERROR; } // Cannot find left partner (ct[j] of base j. Likely a triplet.
	} /* finished finding the left partner of j */
      
      /* After we've found the left partner of j, resolve pks found along the way.
       * Then, remove the pseudoknotted based from cct so we can find the rest of the structure.
       */
      if (esl_stack_ObjectCount(auxpk)) {

	/* init for first pseudoknot */
	leftbound  = cct[j];
	rightbound = leftbound + 1;
	xpk        = -1;            /* start with 'A' if possible again */

	while (esl_stack_IPop(auxpk, &i) == eslOK) {

	  for (k = rightbound-1; k > leftbound; k --) 
	    {
	      if      (cct[k] == 0)          { continue; } 
	      else if (cct[k] >  rightbound) { continue; } 
	      else if (cct[k] == i)          { break; }                  /* i continues the given pseudoknot */
	      else                           { k = leftbound; break; }   /* a new pseudoknot */		    		
	    }
	  
	  if (k == leftbound) /* a new pseudoknot */
	    {
	      // npk ++;   // if you want to count PK's, here's where to do it
	      xpk ++;
	      /* figure out if we can use this alphabet index, or bump it up if necessary */
	      while (i < rb[xpk]) { xpk ++; }
	      
	      leftbound  = (rightbound < cct[i])? rightbound : cct[j];
	      rightbound = cct[i];
	    }
	      
	  npairs_reached ++;
	  if (xpk+(int)('a') <= (int)('z')) {

	    /* update the rightbound of this pk index if necessary */
	    if (cct[i] > rb[xpk]) rb[xpk] = cct[i];
	    
	    /* Add pk indices for this basepair */
	    ss[i-1]      = (char)(xpk+(int)('A'));
	    ss[cct[i]-1] = (char)(xpk+(int)('a'));
	    
	    /* remove pseudoknotted pair from cct */
	    cct[i]     = 0;
	    cct[ct[i]] = 0;
	  }
	  else { status = eslEINVAL; goto ERROR; } // Don't have enough letters to describe all different pseudoknots.
	} 	
      } /* while there is something in auxpk stack */
    } /* finished loop over j: end position on seq, 1..n*/ 
  if (npairs != npairs_reached) ESL_XEXCEPTION(eslEINCONCEIVABLE, "only found %d out of %d pairs.", npairs_reached, npairs);

  status = eslOK;
 ERROR:
  esl_stack_Destroy(pda);
  esl_stack_Destroy(auxpk);
  free(cct);
  return status;
}
/*------------ end, annotation as string or esl_mem -------------*/


/*****************************************************************
 * 2. Functions that assume annotation is a string.
 *
 *    These also allow overwriting: input <ss> and output <ss>
 *    can be the same.
 *****************************************************************/

/* Function:  esl_wuss_full()
 * Incept:    SRE, Mon Feb 28 09:44:40 2005 [St. Louis]
 *
 * Purpose:   Given a simple ("input") WUSS format annotation string <oldss>,
 *            convert it to full ("output") WUSS format in <newss>.
 *            <newss> must be allocated by the caller to be at least as 
 *            long as <oldss>. <oldss> and <newss> can be the same,
 *            to convert a secondary structure string in place.
 *            
 *            Any pseudoknot annotation on <oldss> is preserved exactly
 *            on <newss>. 
 *
 * Returns:   <eslOK> on success.
 *            <eslSYNTAX> if <oldss> isn't in valid WUSS format.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal error that can't happen.
 */
int
esl_wuss_full(const char *oldss, char *newss)
{
  char *tmp = NULL;
  int  *ct  = NULL;
  int   n;
  int   i;
  int   status;

  /* The assignment of base pairs to the main nested structure vs. to
   * pseudoknots, and the A-Z label on each pseudoknotted stem, are
   * arbitrary choices. So if we do the obvious thing of converting to
   * a CT array then using ct2wuss(), for pseudoknotted structures, we
   * might get a full WUSS that's been arbitrarily rearranged relative
   * to the input WUSS. That's undesirable.
   *
   * To work around this, we remove pseudoknots... then to CT, then
   * ct2wuss()... and then put the pseudoknots back exactly where they
   * were.
   */
  n = strlen(oldss);
  ESL_ALLOC(ct,  sizeof(int)  * (n+1));
  ESL_ALLOC(tmp, sizeof(char) * (n+1));
  
  if ((status = esl_wuss_nopseudo(oldss, tmp)) != eslOK) goto ERROR; // tmp = nonpseudoknotted oldss
  if ((status = esl_wuss2ct(tmp, n, ct))       != eslOK) goto ERROR; // ct  = oldss in ct format, no pks 
  if ((status = esl_ct2wuss(ct, n, tmp))       != eslOK) goto ERROR; // now tmp is a full WUSS string, no pks. (eslEINVAL can't happen here)
  
  for (i = 0; i < n; i++)
    newss[i] = ( isalpha(oldss[i]) ? oldss[i] : tmp[i] );  // old pk annotation | new WUSS

  free(ct);
  free(tmp);
  return eslOK;

 ERROR:
  free(ct);
  free(tmp);
  return status;
}



/* Function:  esl_wuss_nopseudo()
 * Incept:    SRE, Tue Feb 15 11:02:43 2005 [St. Louis]
 *
 * Purpose:   Given a WUSS format annotation string <ss1>,
 *            removes all pseudoknot annotation to create a new 
 *            WUSS string <ss2> that contains only a "canonical"
 *            (nonpseudoknotted) structure. <ss2> must be allocated to
 *            be at least as large as <ss1>. <ss1> and <ss2>
 *            may be the same, in which case the conversion is
 *            done in place.
 *
 *            Pseudoknot annotation in <ss1> is simply replaced by <.>
 *            in <ss2>. The resulting <ss2> WUSS string is therefore
 *            in valid simplified format, but may not be valid full
 *            format WUSS.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_wuss_nopseudo(const char *ss1, char *ss2)
{
  while (*ss1 != '\0') 
    {
      if (isalpha(*ss1)) *ss2 = '.';
      else *ss2 = *ss1;
      ss1++;
      ss2++;
    }
  *ss2 = '\0';
  return eslOK;
}


/* Function:  esl_wuss_reverse()
 * Synopsis:  "Reverse complement" a WUSS annotation
 * Incept:    SRE, Wed Feb 10 12:46:51 2016 [JB251 BOS-MCO]
 *
 * Purpose:   If we need to reverse complement a structure-annotated RNA
 *            sequence, we need to "reverse complement" the WUSS
 *            annotation string. Reverse complement the annotation string
 *            <ss> into caller-provided space <new>. To revcomp an annotation 
 *            in place, use <esl_wuss_reverse(ss, ss)>.
 *            
 *            Old SELEX files use a different structure annotation
 *            format, with angle brackets pointing the opposite
 *            direction: \ccode{><} for a base pair. As a convenient
 *            side effect, <esl_wuss_reverse()> will also reverse
 *            complement SELEX annotation lines.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_wuss_reverse(const char *ss, char *new)
{
  int i, n;

  /* first, "complement" the annotation */
  for (i = 0; ss[i] != '\0'; i++)
    {
      if      (isupper(ss[i])) new[i] = tolower(ss[i]);
      else if (islower(ss[i])) new[i] = toupper(ss[i]);
      else {
	switch (ss[i]) {
	case '<': new[i] = '>';   break;
	case '>': new[i] = '<';   break;
	case '(': new[i] = ')';   break;
	case ')': new[i] = '(';   break;
	case '[': new[i] = ']';   break;
	case ']': new[i] = '[';   break;
	case '{': new[i] = '}';   break;
	case '}': new[i] = '{';   break;
	default:  new[i] = ss[i]; break;
	}
      }
    }
  n = i;
  /* Then, reverse it in place. */
  for (i = 0; i < n/2; i++)
    ESL_SWAP(new[i], new[n-i-1], char);

  return eslOK;
}
/*------------------ end, WUSS as string only -------------------*/


/***************************************************************** 
 * 3. Interconversion with obsolete Konings-Hogeweg convention.
 *
 *    "KH" annotation has the angle brackets facing the other way; for
 *    example, ">>>>....<<<<" is a stem loop. COVE, the predecessor to
 *    Infernal, used this convention in ~1993-2002.
 *****************************************************************/

/* Function:  esl_wuss2kh()
 * Incept:    SRE, Tue Feb 15 10:05:35 2005 [St. Louis]
 *
 * Purpose:   Converts a secondary structure string <ss> in 
 *            WUSS notation back to old KHS format in <kh>.
 *            <kh> must be allocated for at least as much
 *            space as <ss>. <kh> may be the same as <ss>,
 *            in which case the conversion is done in-place.
 *
 * Note:      Left bp chars  are converted to >   (left base of base pairs)
 *            Right bp chars are converted to <   (right base of base pairs)
 *            Characters _-,:~ are converted to . (unpaired bases)
 *            Character  .     is untouched       (unpaired)
 *            Everything else is untouched, including any pseudoknot notation.
 * 
 * Returns:   <eslOK> on success.
 */
int
esl_wuss2kh(const char *ss, char *kh)
{
  while (*ss != '\0')
    {
      if       (*ss == '<') *kh = '>';
      else if  (*ss == '(') *kh = '>';
      else if  (*ss == '[') *kh = '>';
      else if  (*ss == '{') *kh = '>';
      else if  (*ss == '>') *kh = '<';
      else if  (*ss == ')') *kh = '<';
      else if  (*ss == ']') *kh = '<';
      else if  (*ss == '}') *kh = '<';
      else if  (*ss == '_') *kh = '.';
      else if  (*ss == '-') *kh = '.';
      else if  (*ss == ',') *kh = '.';
      else if  (*ss == ':') *kh = '.';
      else if  (*ss == '~') *kh = '.';
      else *kh = *ss;
      ss++;
      kh++;
    }
  *kh = '\0';
  return eslOK;
}


/* Function:  esl_kh2wuss()
 * Incept:    SRE, Tue Feb 15 10:10:40 2005 [St. Louis]
 *
 * Purpose:   Converts an old format secondary structure string <kh>
 *            to shorthand WUSS format <ss>. <ss> must be allocated at least
 *            as large as <kh>. <ss> can be identical to <kh>, in which
 *            case the conversion is done in-place.
 *
 * Note:      Character > is converted to <  (left base of base pairs)
 *            Character < is converted to >  (right base of base pairs)
 *            A space is converted to .      (just in case)      
 *
 * Returns:   <eslOK> on success.
 */
int
esl_kh2wuss(const char *kh, char *ss)
{
  while (*kh != '\0')
    {
      if      (*kh == '>') *ss = '<';
      else if (*kh == '<') *ss = '>';
      else if (*kh == ' ') *ss = '.';
      else *ss = *kh;
      kh++;
      ss++;
    }
  *ss = '\0';
  return eslOK;
}
/*----------------- end, KH interconversion ---------------------*/



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef eslWUSS_TESTDRIVE

/* utest_basic()
 *
 * Tests for most functions and their main uses (as opposed to testing
 * pathological edge cases).
 *
 * The provided <ss> can be simple or full WUSS, but when represented
 * as full WUSS, it must be distinguishable from simple WUSS. (This
 * would usually be true; even a hairpin "<<<<....>>>>"
 * vs. "<<<<____>>>>" is distinguishable.)
 */
static void
utest_basic(const char *ss)
{
  char msg[]      = "esl_wuss::utest_basic() failed";
  int  n          = strlen(ss);
  int  nbp_nested = 0;
  int  nbp_pk     = 0;
  int  nbp        = 0;
  int  *ct        = NULL;
  int  *ct2       = NULL;
  char *ss2       = NULL;
  char *ss3       = NULL;
  int   i;
  int   status;

  ESL_ALLOC(ct,  sizeof(int) * (n+1));
  ESL_ALLOC(ct2, sizeof(int) * (n+1));
  ESL_ALLOC(ss2, sizeof(int) * (n+1));
  ESL_ALLOC(ss3, sizeof(int) * (n+1));

  /* <ct> is CT array derived from the input <ss>.
   * We won't change <ct> again, having created it - but we will reuse <ct2>,<ss2>,<ss3>
   */
  if (esl_wuss2ct(ss,  n, ct) != eslOK) esl_fatal(msg);         

  for (i = 0; i < n; i++)                                        // count bps in <ss>
    if      (strchr("{[(<", ss[i]) != NULL) nbp_nested++;
    else if (isupper(ss[i]))                nbp_pk++;
  
  for (i = 1; i <= n; i++)                                       // count bps in corresponding <ct>, check that they match.
    if (ct[i] > i) nbp++;
  if (nbp != nbp_nested + nbp_pk) esl_fatal(msg);

  /* Converting ct => full wuss => ct gives us the same ct
   */
  if (esl_ct2wuss(ct, n, ss2)         != eslOK) esl_fatal(msg);  
  if (esl_wuss_IsFull(ss2, n)         != TRUE)  esl_fatal(msg);
  if (esl_wuss_IsSimple(ss2, n)       != FALSE) esl_fatal(msg);
  if (esl_wuss2ct(ss2, n, ct2)        != eslOK) esl_fatal(msg);
  if (esl_vec_ICompare(ct, ct2, n+1)  != eslOK) esl_fatal(msg);
   
  /* Converting ct => simple wuss => ct also gives the same ct
   */
  if (esl_ct2simplewuss(ct, n, ss2)   != eslOK) esl_fatal(msg);  
  if (esl_wuss_IsSimple(ss2, n)       != TRUE)  esl_fatal(msg);
  if (esl_wuss_IsFull(ss2, n)         != TRUE)  esl_fatal(msg);  // this evaluates TRUE because simple WUSS is a subset of full
  if (esl_wuss2ct(ss2, n, ct2)        != eslOK) esl_fatal(msg);
  if (esl_vec_ICompare(ct, ct2, n+1)  != eslOK) esl_fatal(msg);

  /* Converting simple wuss => full wuss => ct also gives same ct
   * <ss2> is simple WUSS from above.
   */
  if (esl_wuss_full(ss2, ss3)         != eslOK) esl_fatal(msg);
  if (esl_wuss_IsFull(ss3, n)         != TRUE)  esl_fatal(msg);  
  if (esl_wuss_IsSimple(ss3, n)       != FALSE) esl_fatal(msg); // FALSE depends on <ss> being sufficiently complicated
  if (esl_wuss2ct(ss3, n, ct2)        != eslOK) esl_fatal(msg);
  if (esl_vec_ICompare(ct, ct2, n+1)  != eslOK) esl_fatal(msg);
  
  /* Test esl_wuss_nopseudo()
   * simple wuss => ct => simple wuss is guaranteed to give the same WUSS if no PKs (though *not* when there are PKs)
   * <ss2> is still a simple WUSS from above. (<ss> itself may not have been)
   * Don't count bp's in <ct2> and expect it to equal nbp_nested; we got to <ss2> by a wuss->ct->wuss conversion of a PK'd structure
   */
  if (esl_wuss_nopseudo(ss2, ss2)     != eslOK) esl_fatal(msg);   // strip PKs in place
  if (esl_wuss2ct(ss2, n, ct2)        != eslOK) esl_fatal(msg);   // ct2 is now the no-PK version
  if (esl_ct2simplewuss(ct2, n, ss3)  != eslOK) esl_fatal(msg);   // make ss3 from it...
  if (strcmp(ss2, ss3)                != 0)     esl_fatal(msg);   // ss2 == ss3

  /* Test esl_wuss_reverse()
   * Reverse complement of a reverse complement of either simple or full WUSS is the original WUSS.
   */
  if (esl_wuss_reverse(ss,  ss2) != eslOK) esl_fatal(msg);
  if (esl_wuss_reverse(ss2, ss3) != eslOK) esl_fatal(msg);
  if (strcmp(ss, ss3)            != 0)     esl_fatal(msg);

  /* Test the conversions in/out of obsolete KH-style annotation */
  if (esl_wuss2kh(ss,  ss2)           != eslOK) esl_fatal(msg);
  if (esl_kh2wuss(ss2, ss2)           != eslOK) esl_fatal(msg);
  if (esl_wuss2ct(ss2, n, ct2)        != eslOK) esl_fatal(msg);
  if (esl_vec_ICompare(ct, ct2, n+1)  != eslOK) esl_fatal(msg);
 
  free(ct);
  free(ct2);
  free(ss2);
  free(ss3);
  return;

 ERROR:
  esl_fatal(msg);
}
#endif //eslWUSS_TESTDRIVE
/*--------------------- end, unit tests -------------------------*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef eslWUSS_TESTDRIVE

int
main(int argc, char **argv)
{
  /* Some example structure annotations, for real examples with pseudoknots.
   *   ss_examples[0] : Twister consensus [Roth..Breaker, Nat Chem Bio 2014]
   *              [1] : theta ribozyme    [Riccitelli..Luptak, 2014]
   *              [2] : S4 autoregulation [Peselis & Serganov, 2014]
   *              [3] : E. coli RNAse P.  [J Brown figure 10.3.00]
   *
   * Some . and ~'s artificially added to the RNase P to test aligned WUSS.             
   */
  char *ss_examples[] = { "<<<<...<<<<AA<<<..BBB..aa.>>>>>>>..bbb>>>>",
                          "<<<<<<..BBBBBB<<<.AA....>>>>>>>>>aa..........bbbbbb",
                          "...<<<<<<<.<<<<<........AAAA...AAAAA.....BBBB.>>>>>>>>>>>>..........................aaaaa....aaaa...........bbbb.",
                          "\
.....{{{{{{{{{{{{{{{{{{,<<<<<<<<<<<<<-<<<<<____~~~~~>>>>>>>>>->>>>>>>>.....\
>,,,,AAA-AAAAA[[[[---BBBB-[[[[[<<<<<_____>>>>><<<<____>>>->(\
(---(((((,,,,,,,,,,,,<<<<<--<<<<<<<<____>>>>>->>>>>>-->>,,,,\
,,,<<<<<<_______>>>>>><<<<<<<<<____>>>->>>>>->,,)))--))))]]]\
]]]]]],,,<<<<------<<<<<<----<<<<<_bbbb>>>>>>>>>>>----->>>>,\
,,,,,<<<<<<<<____>>>>>>>>,,,,,,,,,,}}}}}}}----------aaaaaaaa\
-}-}}}}}}}}}}::::" };
  int nexamples = sizeof(ss_examples) / sizeof(char *);
  int i;

  fprintf(stderr, "## %s\n", argv[0]);

  for (i = 0; i < nexamples; i++)
    utest_basic(ss_examples[i]);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*eslWUSS_TESTDRIVE*/
/*---------------------- end, test driver -----------------------*/


/***************************************************************** 
 * 6. Example
 *****************************************************************/
#ifdef eslWUSS_EXAMPLE

int
main(void)
{
  /* theta ribozyme drz-Mtgn-1, 51nt, 17bp, 9/17 nested
   * [Riccitelli..Luptak, 2014; Kienbeck..Sigel, bioRxiv 2023] 
   *        5'-GGUAGCACACCUAUGCGUUCCCGUCGCGCUACUGAUUUAGACUAAAUAGGU-3'
   */
  char  ss[] = "<<<<<<..AAAAAA<<<.BB....>>>>>>>>>bb..........aaaaaa";
  int   n;
  int  *ct;   // "CT" array: ct[i=1..n] = 0 if i unpaired, j if i:j basepair. ct[0] unused (=0)
  char *ss2;

  /* ss annotations are length n: ss[0..n-1].
   * A CT array is length n+1: ct[0,1..n], with ct[0] unused.
   *
   * Some esl_wuss functions are used in file parsers, where <ss>
   * might be an esl_mem buffer (no \0 termination); these functions
   * take the length <n> as an argument. Other esl_wuss functions
   * assume that <ss> is a string, with \0 termination.
   */
  n  = strlen(ss);
  if ((ct  = malloc(sizeof(int) * (n+1))) == NULL) esl_fatal("malloc failed");
  if ((ss2 = malloc(sizeof(int) * n))     == NULL) esl_fatal("malloc failed");

  /* This example is "simple" WUSS, a subset of "full" WUSS.
   *
   * We sometimes also refer to simple WUSS as "input" WUSS because
   * it's easier for a person to write, and full WUSS as "output" WUSS
   * because it's the output of Infernal.
   *
   * Simple WUSS is valid as full WUSS, but the reverse is generally
   * not true.
   */
  if ( esl_wuss_IsSimple(ss, n) != TRUE ) esl_fatal("yes it is");
  if ( esl_wuss_IsFull  (ss, n) != TRUE ) esl_fatal("yes it is");

  /* A "CT" array ct[i=1..n] gives the position j that i is base
   * paired to, or 0 if i is unpaired; ct[0] is unused.
   *
   * The name "ct" comes from "Connectivity Table" format, originally
   * used (I believe) by Michael Zuker in MFOLD. CT format has other
   * fields [https://rna.urmc.rochester.edu/Text/File_Formats.html].
   * Easel just needs the basepairing information.
   */
  if (esl_wuss2ct(ss, n, ct) != eslOK)  esl_fatal("conversion to CT failed");

  /* A CT array can be converted to simple or full WUSS.
   *
   * Pseudoknot annotation is not unique: the choice of what to call
   * "main structure" vs "pseudoknot" is arbitrary, and so is the A-Z
   * labelling of each PK.  Converting a pseudoknotted WUSS to CT and
   * back to WUSS may give a different WUSS.
   */
  printf("\n");
  printf("Original SS:           %s\n", ss);

  if (esl_ct2wuss      (ct, n, ss2) != eslOK) esl_fatal("ct2wuss conversion failed");
  printf("full WUSS via CT:      %s\n", ss2);
  
  if (esl_ct2simplewuss(ct, n, ss2) != eslOK) esl_fatal("ct2simplewuss conversion failed");
  printf("simple WUSS via CT:    %s\n", ss2);

  /* Simple WUSS can be converted to full.
   *
   * This function *does* guarantee that the pseudoknot annotation will be unchanged.
   */
  if (esl_wuss_full (ss, ss2) != eslOK) esl_fatal("esl_wuss_full conversion failed");
  printf("full WUSS from simple: %s\n", ss2);
  
  free(ct);
  free(ss2);
  return 0;
}
#endif /*eslWUSS_EXAMPLE*/
/*----------------------- end, example --------------------------*/

