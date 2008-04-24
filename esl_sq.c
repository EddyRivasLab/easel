/* A sequence.
 * 
 * Contents:
 *   1. Text version of the ESL_SQ object.
 *   2. Digitized version of the ESL_SQ object.     [with eslAUGMENT_ALPHABET]
 *   3. Other functions that operate on sequences.
 *   4. Internal functions.
 *   x. Unit tests.
 *   x. Test driver.
 *   x. Example.
 *   x. License and copyright
 * 
 * SRE, Mon Mar 31 17:18:59 2008 [Janelia]
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"	/* alphabet aug adds digital sequences */
#endif 
#ifdef eslAUGMENT_MSA
#include "esl_msa.h"		/* msa aug adds ability to extract sq from an MSA  */
#endif
#include "esl_sq.h"

/* Shared parts of text/digital creation functions (defined in "internal functions" section) */
static ESL_SQ *sq_create(int do_digital);
static ESL_SQ *sq_create_from(const char *name, const char *desc, const char *acc);


/*****************************************************************
 *# 1. Text version of the <ESL_SQ> object.
 *****************************************************************/

/* Function:  esl_sq_Create()
 * Incept:    SRE, Thu Dec 23 11:57:00 2004 [Zaragoza]
 *
 * Purpose:   Creates an empty <ESL_SQ> sequence object, in text mode, with
 *            internal fields allocated to reasonable initial sizes. 
 *            
 * Args:      (void)
 *
 * Returns:   a pointer to the new <ESL_SQ>. Caller frees this with
 *            <esl_sq_Destroy()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
ESL_SQ *
esl_sq_Create(void)
{
  return sq_create(FALSE);
}

/* Function:  esl_sq_CreateFrom()
 * Incept:    SRE, Wed Mar 22 09:17:04 2006 [St. Louis]
 *
 * Purpose:   Create a new <ESL_SQ> object in text mode from elemental data.
 *            This provides an interface between non-Easel code
 *            and Easel's object.
 *            
 *            Makes copies of all data. Caller is still
 *            responsible for memory of name, seq, etc.
 *            
 *            <ss> is an optional alphabetic secondary structure 
 *            annotation string. If provided, its length must match
 *            the length of <seq>.
 *            
 * Args:      name    -  name of the sequence (NUL-terminated)
 *            seq     -  the sequence (alphabetic; NUL-terminated)
 *            desc    -  optional: description line (or NULL)
 *            acc     -  optional: accession (or NULL)
 *            ss      -  optional: secondary structure annotation (or NULL)
 *
 * Returns:   a pointer to the new object. Free with
 *            <esl_sq_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SQ *
esl_sq_CreateFrom(const char *name, const char *seq, const char *desc, const char *acc, const char *ss)
{
  ESL_SQ *sq = NULL;
  int     n  = strlen(seq);
  int     status;

  if ((sq     = sq_create_from(name, desc, acc)) == NULL)  goto ERROR;
  if ((status = esl_strdup(seq, n, &(sq->seq)))  != eslOK) goto ERROR;

  if (ss != NULL) 
    {
      if (strlen(ss) != n) ESL_XEXCEPTION(eslEINVAL, "ss, seq lengths mismatch");
      if ((status = esl_strdup(ss, n, &(sq->ss))) != eslOK) goto ERROR;
    } 
  else sq->ss = NULL;

  sq->n      = n;
  sq->salloc = n+1;
  return sq;

 ERROR:
  esl_sq_Destroy(sq);
  return NULL;
}

/* Function:  esl_sq_Grow()
 * Incept:    SRE, Wed Jan 10 08:26:23 2007 [Janelia]
 *
 * Purpose:   Assure that the sequence <sq> can hold at least
 *            one more residue, whether in digital or text mode.
 *            Reallocate if necessary. Optionally returns the number
 *            of residues that can be added before the next call
 *            to <esl_sq_Grow()> in <opt_nsafe>.
 *            
 *            The terminal <NUL> or sentinel count as a 'residue' for
 *            allocation purposes: that is, you may need to call
 *            <esl_sq_Grow()> before terminating a new sequence.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. In this case, the
 *            original <sq> is untouched, and <*opt_nsafe> is returned
 *            as 0.
 *
 * Xref:      STL11/125.
 */
int
esl_sq_Grow(ESL_SQ *sq, int *opt_nsafe)
{
  void *tmp;
  int   new;
  int   nsafe;
  int   status;

  if (sq->seq != NULL)  nsafe = sq->salloc - sq->n;         /* text */
  else                  nsafe = (sq->salloc-1) - sq->n;     /* digital: -1 because 0 is a sentinel       */

  if (nsafe < 1)
    {  /* reallocate by doubling (shouldn't need more, but if we do, keep doubling) */
      new = sq->salloc;
      do { nsafe += new; new*=2; } while (nsafe < 1);
      
      if (sq->seq != NULL) ESL_RALLOC(sq->seq, tmp, new * sizeof(char));	/* text    */
      else                 ESL_RALLOC(sq->dsq, tmp, new * sizeof(ESL_DSQ));	/* digital */
      if (sq->ss != NULL)  ESL_RALLOC(sq->ss,  tmp, new * sizeof(char));
      sq->salloc = new;
    }
  if (opt_nsafe != NULL) *opt_nsafe = nsafe;
  return eslOK;

 ERROR:
  if (opt_nsafe != NULL) *opt_nsafe = 0;
  return status;
}

/* Function:  esl_sq_GrowTo()
 * Synopsis:  Grows an <ESL_SQ> to hold a seq of at least <n> residues.
 * Incept:    SRE, Fri Jan 18 11:06:50 2008 [UA5233 Westchester-Dulles]
 *
 * Purpose:   Assure that the appropriate (text or digital) sequence
 *            field in <sq> can hold up to a total of <n> residues,
 *            reallocating as needed.
 *            
 *            If reallocated, the allocation will be $\geq (n+1)$ for
 *            text mode (the +1 is for the terminal NUL byte), $\geq
 *            (n+2)$ for digital mode (+2 for sentinel bytes at each
 *            end). That is, you don't need to take these extra bytes into
 *            account in your <n>; <n> is the number of residues, not
 *            bytes.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sq_GrowTo(ESL_SQ *sq, int n)
{
  void *tmp;
  int   status;

  if (sq->seq != NULL)		/* text mode */
    {
      if (n+1 > sq->salloc) {
	ESL_RALLOC(sq->seq, tmp, (n+1) * sizeof(char));
	if (sq->ss != NULL) ESL_RALLOC(sq->ss, tmp, (n+1) * sizeof(char));
	sq->salloc = n+1;
      }
    }
  else				/* digital mode */
    {
      if (n+2 > sq->salloc) {
	ESL_RALLOC(sq->dsq, tmp, (n+2) * sizeof(ESL_DSQ));
	if (sq->ss != NULL) ESL_RALLOC(sq->ss, tmp, (n+2) * sizeof(char));
	sq->salloc = n+2;
      }
    }
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_sq_Copy()
 * Synopsis:  Make a copy of an <ESL_SQ>.
 * Incept:    SRE, Sun Feb 24 17:59:24 2008 [UA5315 to St. Louis]
 *
 * Purpose:   Copies a source sequence object <src> into 
 *            destination sequence object <dst>.
 *            
 *            The two objects don't have to be matched as far as
 *            text/digital mode go; if mismatched, appropriate
 *            text/digital conversion will be done.
 *            
 *            The destination sequence <dst> is reallocated internally
 *            as necessary to hold a copy of <src>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 * 
 * Note:      Note the shenanigans involved in copying ss; you have
 *            to pay attention to the ss being a 0..n-1 string in text
 *            mode versus a 1..n string in digital mode.
 */
int
esl_sq_Copy(const ESL_SQ *src, ESL_SQ *dst)
{
  int status;

  /* If <src> has structure annotation and <dst> does not, initialize an allocation in <dst> */
  if (src->ss != NULL) ESL_ALLOC(dst->ss, sizeof(char) * dst->salloc);

  if ((status = esl_sq_SetName     (dst, src->name)) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(dst, src->acc))  != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (dst, src->desc)) != eslOK) goto ERROR;
  if ((status = esl_sq_GrowTo      (dst, src->n))    != eslOK) goto ERROR;

  if (src->seq != NULL && dst->seq != NULL) /* text to text */
    {
      strcpy(dst->seq, src->seq);
      if (src->ss != NULL) strcpy(dst->ss, src->ss);
    }
#ifdef eslAUGMENT_ALPHABET
  else if (src->seq != NULL && dst->dsq != NULL) /* text to digital */
    {
      if ((status = esl_abc_Digitize(dst->abc, src->seq, dst->dsq)) != eslOK) goto ERROR;      
      if (src->ss != NULL) {
	strcpy(dst->ss+1, src->ss);
	dst->ss[0] = '\0';
      }
    }
  else if (src->dsq != NULL && dst->seq != NULL) /* digital to text */
    {
      if ((status = esl_abc_Textize(src->abc, src->dsq, src->n, dst->seq)) != eslOK) goto ERROR;
      if (src->ss != NULL) strcpy(dst->ss, src->ss+1);
    }
  else 				/* digital to digital */
    {
      if (src->abc->type != dst->abc->type) 
	ESL_XEXCEPTION(eslEINCOMPAT, "seq objects involved in Copy differ in digital alphabet");
      if ((status = esl_abc_dsqcpy(src->dsq, src->n, dst->dsq)) != eslOK) goto ERROR;
      if (src->ss != NULL) {
	strcpy(dst->ss+1, src->ss+1);
	dst->ss[0] = '\0';
      }
    }
#endif
  
  dst->n     = src->n;
  dst->roff  = src->roff;
  dst->doff  = src->doff;
  /* don't copy allocations (nalloc, etc); dst knows its own memory */
  return eslOK;

 ERROR:
  esl_sq_Reuse(dst);
  return status;
}

/* Function:  esl_sq_Reuse()
 * Incept:    SRE, Thu Dec 23 12:23:51 2004 [Zaragoza]
 *
 * Purpose:   Given a sequence object <sq> already in use;
 *            reinitialize all its data, so a new seq
 *            may be read into it. This allows sequential sequence
 *            input without a lot of wasted allocation/free cycling.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_sq_Reuse(ESL_SQ *sq)
{
  sq->name[0] = '\0';
  sq->acc[0]  = '\0';
  sq->desc[0] = '\0';
  if (sq->seq != NULL) sq->seq[0] = '\0';
  if (sq->dsq != NULL) sq->dsq[0] = sq->dsq[1] = eslDSQ_SENTINEL;
  if (sq->ss  != NULL) {
    if (sq->seq != NULL) sq->ss[0] = '\0';
    else                 sq->ss[0] = sq->ss[1] = '\0'; /* in digital mode, ss string is 1..n; 0 is a dummy \0 byte*/
  }
  sq->n    = 0;
  sq->doff = -1;
  sq->roff = -1;
  return eslOK;
}

/* Function:  esl_sq_Destroy()
 * Incept:    SRE, Thu Dec 23 12:28:07 2004 [Zaragoza]
 *
 * Purpose:   Free a Create()'d <sq>.
 */
void
esl_sq_Destroy(ESL_SQ *sq)
{
  if (sq == NULL) return;

  if (sq->name   != NULL) free(sq->name);  
  if (sq->acc    != NULL) free(sq->acc);   
  if (sq->desc   != NULL) free(sq->desc);  
  if (sq->seq    != NULL) free(sq->seq);   
  if (sq->dsq    != NULL) free(sq->dsq);   
  if (sq->ss     != NULL) free(sq->ss);    
  free(sq);
  return;
}
/*--------------- end of ESL_SQ object functions ----------------*/





/*****************************************************************
 *# 2. Digitized version of the <ESL_SQ> object.
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  esl_sq_CreateDigital()
 * Incept:    SRE, Tue Jan  9 16:42:35 2007 [Janelia]
 *
 * Purpose:   Same as <esl_sq_Create()>, except the returned <sq> is
 *            configured for a digital sequence using internal
 *            alphabet <abc>, rather than a text sequence. Creates an
 *            empty digital <ESL_SQ> sequence object, with internal
 *            fields allocated to reasonable initial sizes.
 *            
 * Args:      abc      - pointer to internal alphabet
 * 
 * Returns:   a pointer to the new <ESL_SQ>. Caller frees this with
 *            <esl_sq_Destroy()>.
 * 
 * Throws:    <NULL> if an allocation fails.
 *
 * Xref:      STL11/124
 */
ESL_SQ *
esl_sq_CreateDigital(const ESL_ALPHABET *abc)
{
  ESL_SQ *s;
  if ((s = sq_create(TRUE)) == NULL) return NULL;
  s->abc    = abc;
  return s;
}

/* Function:  esl_sq_CreateDigitalFrom()
 * Incept:    EPN, Fri Aug 24 13:38:56 2007
 *
 * Purpose:   Create a new <ESL_SQ> object from elemental data.
 *            Same as <esl_sq_CreateFrom> except takes a digital <ESL_DSQ *dsq>
 *            instead of a text <char *seq> as the sequence to copy.
 *            
 *            Makes copies of all data. Caller is still
 *            responsible for memory of name, seq, etc.
 *            
 *            <ss> is an optional alphabetic secondary structure
 *            annotation string, <0..L-1>. If provided, its length
 *            must match the length of <seq>. (Note that although the
 *            argument <ss> is provided as a standard <0..L-1> C
 *            string, <ss> is stored internally as a <1..L> string in
 *            a digital sequence object, so that both the digital
 *            sequence and its alphabetic secondary structure
 *            annotation are indexed the same.)
 *            
 *            The object is growable; you can use <esl_sq_Reuse()>
 *            on it.
 *
 * Args:      name    -  name of the sequence
 *            dsq     -  digital sequence <1..L>
 *            L       -  length of digitized sequence in residues (or -1 if unknown)
 *            desc    -  optional: description line (or NULL)
 *            acc     -  optional: accession (or NULL)
 *            ss      -  optional: secondary structure annotation (or NULL)
 *
 * Returns:   a pointer to the new object. Free with
 *            <esl_sq_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SQ *
esl_sq_CreateDigitalFrom(const ESL_ALPHABET *abc, const char *name, const ESL_DSQ *dsq, int L,
			 const char *desc, const char *acc, const char *ss)
{
  ESL_SQ *sq = NULL;
  int     status;

  if((sq = sq_create_from(name, desc, acc)) == NULL) goto ERROR;
  sq->n = (L == -1) ? esl_abc_dsqlen(dsq) : L;
  if ((status = esl_abc_dsqdup(dsq, sq->n, &(sq->dsq))) != eslOK) goto ERROR;

  if (ss != NULL)
    {
      if (strlen(ss) != sq->n) ESL_XEXCEPTION(eslEINVAL, "ss, seq lengths mismatch");
      ESL_ALLOC(sq->ss, sizeof(char) * (sq->n+2));
      sq->ss[0] = '\0';
      strcpy(sq->ss+1, ss);
    }

  sq->salloc = sq->n+2;
  sq->abc   =  abc;
  return sq;

 ERROR:
  esl_sq_Destroy(sq);
  return NULL;
}


/* Function:  esl_sq_Digitize()
 * Incept:    EPN, Mon Feb 12 11:09:06 2007
 *
 * Purpose:   Given a sequence <sq> in text mode, convert it to
 *            digital mode, using alphabet <abc>.
 *            
 *            Internally, the <dsq> digital sequence field is filled,
 *            the <seq> text alignment field is destroyed and free'd,
 *            a copy of the alphabet pointer is kept in the sq's
 *            <abc> reference.
 *
 * Args:      abc    - digital alphabet
 *            sq     - sequence to digitize
 *
 * Returns:   <eslOK> on success;
 *            <eslEINVAL> if the sequence contains invalid characters
 *            that can't be digitized. If this happens, the <sq> is returned
 *            unaltered - left in text mode, with <seq> as it was. (This is
 *            a normal error, because <sq->seq> may be user input that we 
 *            haven't validated yet.)
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, state of <sq> may be 
 *            wedged, and it should only be destroyed, not used.
 */
int
esl_sq_Digitize(const ESL_ALPHABET *abc, ESL_SQ *sq)
{
  int status;

  /* Contract checks */
  if (sq->dsq   != NULL) return eslOK;
  if (sq->seq   == NULL) ESL_EXCEPTION(eslEINVAL, "sq has no text sequence");

  /* Validate before we convert, so we leave <seq> untouched if it's bad. */
  if (esl_abc_ValidateSeq(abc, sq->seq, sq->n, NULL) != eslOK) return eslEINVAL;

  /* Allocate dsq, ss properly; these are our last failure points. */
  if (sq->salloc < sq->n+2) {	/* it's possible (though unlikely) for salloc to be 1 residue too small */
    sq->salloc = sq->n+2;
    if (sq->ss != NULL) {
      void *tmp;
      ESL_RALLOC(sq->ss, tmp, sizeof(char) * sq->salloc);
    }
  }
  ESL_ALLOC(sq->dsq, (sq->salloc) * sizeof(ESL_DSQ));

  /* Now convert. */
  if ((status = esl_abc_Digitize(abc, sq->seq, sq->dsq)) != eslOK) goto ERROR;
  if (sq->ss != NULL) {
    memmove(sq->ss+1, sq->ss, sq->n+1);
    sq->ss[0] = '\0';
  }
  free(sq->seq);
  sq->seq = NULL;
  sq->abc =  abc;
  return eslOK;

 ERROR:
  if (sq->dsq != NULL) free(sq->dsq);
  return status;
}

/* Function:  esl_sq_Textize()
 * Incept:    EPN, Mon Feb 12 11:15:06 2007
 *
 * Purpose:   Given a sequence <sq> in digital mode, convert it
 *            to text mode.
 *            
 *            Internally, the <seq> text alignment field is filled, the
 *            <dsq> digital alignment field is destroyed and free'd, the
 *            sq's <abc> digital alphabet reference is nullified.
 *            
 * Args:      sq   - sequence object to convert to text
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslECORRUPT> if the digitized alignment string contains 
 *                          invalid characters.
 */
int
esl_sq_Textize(ESL_SQ *sq)
{
  int status;

  /* Contract checks */
  if (sq->seq  != NULL) return eslOK;
  if (sq->dsq  == NULL) ESL_EXCEPTION(eslEINVAL, "sq has no digital sequence");
  if (sq->abc  == NULL) ESL_EXCEPTION(eslEINVAL, "sq has no digital alphabet");

  /* Allocate. salloc is guaranteed big enough, if it was big enough for digital. */
  ESL_ALLOC(sq->seq, sq->salloc * sizeof(char));
  
  /* Convert. */
  if ((status = esl_abc_Textize(sq->abc, sq->dsq, sq->n, sq->seq)) != eslOK) goto ERROR;
  if (sq->ss != NULL) 
    memmove(sq->ss, sq->ss+1, sq->n+1);	/* slide back to 0..n-1; +1 includes terminal \0 */

  free(sq->dsq);
  sq->dsq = NULL;
  sq->abc = NULL;           /* nullify reference (caller still owns real abc) */
  return eslOK;

 ERROR:
  if (sq->seq != NULL) free(sq->seq);
  return status;
}

/* Function:  esl_sq_GuessAlphabet()
 * Synopsis:  Guess alphabet type of a sequence.
 * Incept:    SRE, Wed May 16 11:03:44 2007 [Janelia]
 *
 * Purpose:   Guess the alphabet type of biosequence <sq>, and store the
 *            guess in <*ret_type>.
 *            
 *            All 26 letters are valid in the amino alphabet (even O
 *            and J now), so the DNA alphabet is necessarily a subset.
 *            Therefore most protein sequences can be identified
 *            unambiguously (because they use letters that only occur
 *            in amino acid sequence), but DNA sequences cannot be.
 *            
 *            The sequence must contain more than 10 residues, or it
 *            is called <eslUNKNOWN>.
 *            
 *            Specifically, this routine calls the sequence <eslDNA>
 *            if it consists only of the residues ACGTN and all four
 *            of ACGT occur. (And analogously for <eslRNA>, ACGU$+$N.)
 *            It calls the sequence <eslAMINO> either if it contains
 *            an amino-specific letter (EFIJLOPQZ), or if it contains
 *            at least 15 of the 20 canonical amino acids and consists
 *            only of canonical amino acids or X.

 *            Thus DNA sequences containing IUPAC degeneracies other
 *            than N are called <eslUNKNOWN>, rather than hazarding a
 *            guess. It may be possible to improve on this in the
 *            future by using residue occurrence frequencies.
 *            
 *            Note that a sequence of "ATATATA..." will be called
 *            <eslUNKNOWN>, whereas a sequence "ACGTACGTACGT..."
 *            (which could conceivably be "ala-cys-gly-thr...") will
 *            be called <eslDNA>. Peptides of simple mono and di-amino
 *            acid compositions are known, but I have not (yet) seen a
 *            peptide consisting only of all four residues ACGT.
 *            
 *            The routine is designed to be conservative, calling
 *            <eslUNKNOWN> rather than making errors. In a test on the
 *            Oct 2006 version of the NCBI nonredundant databases,
 *            this routine called 0 <eslDNA> and 5694 <eslUNKNOWN> on
 *            4.0M protein sequences (99.9% classification, no errors)
 *            and 0 <eslAMINO> and 155756 <eslUNKNOWN> in 4.4M DNA
 *            sequences (96% classification, no errors). (Actually,
 *            one DNA call was made in the protein database. That
 *            entry was indeed a DNA contaminant, and it has since
 *            been removed by NCBI.)
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set to
 *            <eslAMINO>, <eslRNA>, or <eslDNA>.
 *
 *            Returns <eslEAMBIGUOUS> if unable to determine the
 *            alphabet type; in this case, <*ret_type> is set to 
 *            <eslUNKNOWN>.
 *
 * Xref:      J1/62; 2007/0517-easel-guess-alphabet
 */
int
esl_sq_GuessAlphabet(ESL_SQ *sq, int *ret_type)
{
  int ct[26];
  int x, i;
  int n = 0;

  for (x = 0; x < 26; x++) ct[x] = 0;
  for (i = 0; i < sq->n; i++) {
    x = toupper(sq->seq[i]) - 'A';
    if (x < 0 || x > 26) continue;
    ct[x]++;
    n++;
    if (n > 10000) break;	/* we oughta know by now! */
  }
  return esl_abc_GuessAlphabet(ct, ret_type);
}
#endif /*eslAUGMENT_ALPHABET*/

/*---------- end of digitized ESL_SQ object functions -----------*/



/*****************************************************************
 *# 3. Other functions that operate on sequences.
 *****************************************************************/

/* Function:  esl_sq_SetName()
 * Synopsis:  Format and set a name of a sequence.
 * Incept:    SRE, Thu Jan 11 08:42:53 2007 [Janelia]
 *
 * Purpose:   Set the name of the sequence <sq> to <name>, reallocating
 *            as needed. <name> can be a <printf()>-style format with
 *            arguments; for example, <esl_sq_SetName(sq, "random%d", i)>.
 * 
 *            A copy of <name> is made, so if caller had <name> allocated, 
 *            it is still responsible for freeing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      STL11/125
 */
int
esl_sq_SetName(ESL_SQ *sq, char *name, ...)
{
  va_list argp;
  va_list argp2;
  int   n;
  void *tmp;
  int   status;

  va_start(argp, name);
  va_copy(argp2, argp);
  if ((n = vsnprintf(sq->name, sq->nalloc, name, argp)) > sq->nalloc)
    {
      ESL_RALLOC(sq->name, tmp, sizeof(char) * (n+1)); 
      sq->nalloc = n+1;
      vsnprintf(sq->name, sq->nalloc, name, argp2);
    }
  va_end(argp);
  va_end(argp2);
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_sq_SetAccession()
 * Synopsis:  Format and set the accession field in a sequence.
 * Incept:    SRE, Fri Jan 18 09:48:54 2008 [Westchester airport]
 *
 * Purpose:   Set the accession of the sequence <sq> to <acc>, reallocating
 *            as needed. <acc> can be a <printf()>-style format with
 *            arguments; for example, <esl_sq_SetAccession(sq, "ACC%06d", i)>.
 * 
 *            A copy of <acc> is made, so if caller had <acc> allocated, 
 *            it is still responsible for freeing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      STL11/125
 */
int
esl_sq_SetAccession(ESL_SQ *sq, char *acc, ...)
{
  va_list argp, argp2;
  int     n;
  void   *tmp;
  int     status;

  va_start(argp, acc);
  va_copy(argp2, argp);
  if ((n = vsnprintf(sq->acc, sq->aalloc, acc, argp)) > sq->aalloc)
    {
      ESL_RALLOC(sq->acc, tmp, sizeof(char) * (n+1)); 
      sq->aalloc = n+1;
      vsnprintf(sq->acc, sq->aalloc, acc, argp2);
    }
  va_end(argp);
  va_end(argp2);
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_sq_SetDesc()
 * Synopsis:  Format and set the description field in a sequence.
 * Incept:    SRE, Fri Jan 18 09:46:14 2008 [Westchester airport]
 *
 * Purpose:   Set the description of the sequence <sq> to <desc>, reallocating
 *            as needed. <desc> can be a <printf()>-style format with
 *            arguments; for example, <esl_sq_SetDesc(sq, "random sequence %d", i)>.
 * 
 *            A copy of <desc> is made, so if caller had <desc> allocated, 
 *            it is still responsible for freeing it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      STL11/125
 */
int
esl_sq_SetDesc(ESL_SQ *sq, char *desc, ...)
{
  va_list argp, argp2;
  int     n;
  void   *tmp;
  int     status;

  va_start(argp, desc);
  va_copy(argp2, argp);
  if ((n = vsnprintf(sq->desc, sq->dalloc, desc, argp)) > sq->dalloc)
    {
      ESL_RALLOC(sq->desc, tmp, sizeof(char) * (n+1)); 
      sq->dalloc = n+1;
      vsnprintf(sq->desc, sq->dalloc, desc, argp2);
    }
  va_end(argp);  
  va_end(argp2);
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_sq_CAddResidue()
 * Synopsis:  Add one residue (or terminal NUL) to a text seq.
 * Incept:    SRE, Wed Jan 10 07:58:20 2007 [Janelia]
 *
 * Purpose:   Add one residue <c> onto a growing text mode sequence <sq>,
 *            and deal with any necessary reallocation.
 *
 *            The sequence in <sq> is not <NUL>-terminated. To 
 *            finish and NUL-terminate <sq>, call 
 *            <esl_sq_CAddResidue(sq, 0)>.
 *            
 * Note:      Not the most efficient routine, but convenient in some
 *            routines. Parsers (where speed is at a premium) typically
 *            use an addseq() kind of function instead.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on re-allocation failure.
 *
 * Xref:      STL11/125.
 */
int
esl_sq_CAddResidue(ESL_SQ *sq, char c)
{
  if (esl_sq_Grow(sq, NULL) != eslOK) return eslEMEM;
  sq->seq[sq->n] = c;
  if (c != '\0') sq->n++;
  return eslOK;
}

#ifdef eslAUGMENT_ALPHABET
/* Function:  esl_sq_XAddResidue()
 * Synopsis:  Add one residue (or terminal sentinel) to digital seq.
 * Incept:    SRE, Wed Jan 10 08:23:23 2007 [Janelia]
 *
 * Purpose:   Like <esl_sq_CAddResidue()>, except for digital mode
 *            sequence: add a digital residue <x> onto a growing
 *            digital sequence <sq>. 
 *            
 *            The digital sequence in <sq> must be explicitly
 *            terminated when you're done adding to it; call
 *            <esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on re-allocation failure.
 *
 * Xref:      STL11/125.
 */
int
esl_sq_XAddResidue(ESL_SQ *sq, ESL_DSQ x)
{
  if (esl_sq_Grow(sq, NULL) != eslOK) return eslEMEM;
  sq->dsq[sq->n+1] = x;
  if (x != eslDSQ_SENTINEL) sq->n++;
  return eslOK;
}
#endif /* eslAUGMENT_ALPHABET */

#ifdef eslAUGMENT_MSA

/* Function:  esl_sq_GetFromMSA()
 * Synopsis:  Get a single sequence from an MSA.
 * Incept:    SRE, Tue Apr  1 11:13:28 2008 [Janelia]
 *
 * Purpose:   Retrieve sequence number <which> (<0..msa->nseq-1>) from
 *            <msa> and store it in the <sq> that the caller allocated
 *            and provided. This version (as opposed to
 *            <esl_sq_FetchFromMSA()>, below) allows caller to reuse
 *            the same <sq> container for retrieving sequences one at
 *            a time from an MSA.
 *            
 *            The retrieved sequence <sq> must be in the same mode as
 *            the source <msa>, text versus digital.
 * 
 *            The retrieved sequence is dealigned. For a text mode
 *            sequence, gap characters to be removed are assumed to be
 *            <-_.>. For a digital mode sequence, gap characters are
 *            defined by the digital alphabet.
 *
 * Returns:   <eslOK> on success, and the retrieved sequence is in <sq>.
 *            Some of the internals of <sq> may have been reallocated if
 *            necessary. 
 *            
 *            Returns <eslEOD> if there is no sequence number <which>.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslEINVAL> if <sq> is in a different text/digital mode than
 *            <msa>.
 */
int
esl_sq_GetFromMSA(const ESL_MSA *msa, int which, ESL_SQ *sq)
{
  char   *gapchars = "-_.";	/* hardcoded for now */
  char   *acc      = NULL;
  char   *desc     = NULL;
  char   *ss       = NULL;
  int     status;

  if (which >= msa->nseq || which < 0) return eslEOD;
  if ( (msa->flags & eslMSA_DIGITAL) && sq->dsq == NULL) ESL_XEXCEPTION(eslEINVAL, "msa is digital, sq is not");
  if (!(msa->flags & eslMSA_DIGITAL) && sq->seq == NULL) ESL_XEXCEPTION(eslEINVAL, "msa is text, sq is not");

  /* watch out for optional msa annotations being totally NULL */
  if (msa->sqacc  != NULL) acc  = msa->sqacc[which]; 
  if (msa->sqdesc != NULL) desc = msa->sqdesc[which];
  if (msa->ss     != NULL) ss   = msa->ss[which]; 

  if ((status = esl_sq_SetName     (sq, msa->sqname[which])) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(sq, acc))                != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (sq, desc))               != eslOK) goto ERROR;
  if ((status = esl_sq_GrowTo      (sq, msa->alen))          != eslOK) goto ERROR; /* can't be more than alen residues */

  if (! msa->flags & eslMSA_DIGITAL) /* text mode to text mode */
    {
      strcpy(sq->seq, msa->aseq[which]);
      if (ss != NULL) strcpy(sq->ss, msa->ss[which]);
      esl_strdealign(sq->ss,  sq->seq, gapchars, NULL);
      esl_strdealign(sq->seq, sq->seq, gapchars, &(sq->n));      
    }
#ifdef eslAUGMENT_ALPHABET
  else
    {
      esl_abc_dsqcpy(msa->ax[which], msa->alen, sq->dsq);
      if (ss != NULL) { strcpy(sq->ss+1, ss); sq->ss[0] = '\0'; }
      esl_abc_CDealign(sq->abc, sq->ss+1, sq->dsq, NULL);
      esl_abc_XDealign(sq->abc, sq->dsq,  sq->dsq, &(sq->n));      
    }
#endif /*eslAUGMENT_ALPHABET*/
  
  sq->roff = -1;
  sq->doff = -1;
  return eslOK;

 ERROR:
  return status;
}




/* Function:  esl_sq_FetchFromMSA()
 * Synopsis:  Fetch a single sequence from an MSA.
 * Incept:    SRE, Sun Mar 30 13:39:06 2008 [Janelia]
 *
 * Purpose:   Retrieve sequence number <which> (<0..msa->nseq-1>) from <msa>, in a newly
 *            allocated sequence object; return a pointer to this object
 *            in <ret_sq>.
 * 
 *            The retrieved sequence is in the same mode as the source
 *            <msa>, text versus digital.
 * 
 *            The retrieved sequence is dealigned. For a text mode
 *            sequence, gap characters to be removed are assumed to be
 *            <-_.>. For a digital mode sequence, gap characters are
 *            defined by the digital alphabet.
 *
 * Returns:   <eslOK> on success, and a pointer to the newly fetched
 *            sequence is in <*ret_sq>, which caller is responsible for freeing.
 *            
 *            Returns <eslEOD> if there is no sequence number <which>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sq_FetchFromMSA(const ESL_MSA *msa, int which, ESL_SQ **ret_sq)
{
  ESL_SQ *sq       = NULL;
  char   *acc      = NULL;
  char   *desc     = NULL;
  char   *ss       = NULL;
  char   *gapchars = "-_.";	/* hardcoded for now */

  if (which >= msa->nseq || which < 0) return eslEOD;

  /* watch out for optional msa annotations being totally NULL */
  if (msa->sqacc  != NULL) acc  = msa->sqacc[which]; 
  if (msa->sqdesc != NULL) desc = msa->sqdesc[which];
  if (msa->ss     != NULL) ss   = msa->ss[which]; 

  if (! (msa->flags & eslMSA_DIGITAL)) /* text mode MSA to text mode sequence */
    {
      if ((sq = esl_sq_CreateFrom(msa->sqname[which], msa->aseq[which], desc, acc, ss)) == NULL) goto ERROR;
      if (sq->ss != NULL) esl_strdealign(sq->ss,  sq->seq, gapchars, NULL);
      esl_strdealign(sq->seq, sq->seq, gapchars, &(sq->n));
    }
#ifdef eslAUGMENT_ALPHABET
  else				/* digital mode MSA to digital mode sequence */
    {
      if ((sq = esl_sq_CreateDigitalFrom(msa->abc, msa->sqname[which], msa->ax[which], msa->alen, desc, acc, ss)) == NULL) goto ERROR; 
      if (sq->ss != NULL) esl_abc_CDealign(sq->abc, sq->ss+1, sq->dsq, NULL);
      esl_abc_XDealign(sq->abc, sq->dsq,  sq->dsq, &(sq->n));
    }
#endif

  *ret_sq = sq;
  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  *ret_sq = NULL;
  return eslEMEM;
}
#endif /*eslAUGMENT_MSA*/
/*---------- end, other functions in the API --------------------*/





/*****************************************************************
 * 4. Internal functions
 *****************************************************************/

/* Create and CreateDigital() are almost identical, so
 * their shared guts are here:
 */
static ESL_SQ *
sq_create(int do_digital)
{
  ESL_SQ *sq = NULL;
  int status;

  ESL_ALLOC(sq, sizeof(ESL_SQ));

  sq->name     = NULL;
  sq->acc      = NULL;
  sq->desc     = NULL;
  sq->seq      = NULL;
  sq->dsq      = NULL;	
  sq->ss       = NULL;		/* Note that ss is optional - it will only be allocated if needed */

  sq->nalloc   = eslSQ_NAMECHUNK;	
  sq->aalloc   = eslSQ_ACCCHUNK;
  sq->dalloc   = eslSQ_DESCCHUNK;
  sq->salloc   = eslSQ_SEQCHUNK; 

  ESL_ALLOC(sq->name, sizeof(char) * sq->nalloc);
  ESL_ALLOC(sq->acc,  sizeof(char) * sq->aalloc);
  ESL_ALLOC(sq->desc, sizeof(char) * sq->dalloc);
  if (do_digital) ESL_ALLOC(sq->dsq,  sizeof(ESL_DSQ) * sq->salloc);
  else            ESL_ALLOC(sq->seq,  sizeof(char)    * sq->salloc);

  esl_sq_Reuse(sq);	/* initialization of sq->n, offsets, and strings */
  return sq;

 ERROR:
  esl_sq_Destroy(sq);
  return NULL;
}  

/* CreateFrom and CreateDigitalFrom() are almost identical, so
 * their shared guts are here:
 */
static ESL_SQ *
sq_create_from(const char *name, const char *desc, const char *acc)
{
  ESL_SQ *sq = NULL;
  int     n;
  int     status;

  ESL_ALLOC(sq, sizeof(ESL_SQ));
  sq->name   = NULL;
  sq->acc    = NULL;
  sq->desc   = NULL;
  sq->seq    = NULL;
  sq->dsq    = NULL;
  sq->ss     = NULL;
  
  if (name != NULL)
    {
      n = strlen(name)+1;
      ESL_ALLOC(sq->name, sizeof(char) * n);
      strcpy(sq->name, name);
      sq->nalloc = n;
    }
  else 
    {
      sq->nalloc = eslSQ_NAMECHUNK;
      ESL_ALLOC(sq->name, sizeof(char) * sq->nalloc);
      sq->name[0] = '\0';
    }
  
  if (desc != NULL) 
    {
      n = strlen(desc)+1;
      ESL_ALLOC(sq->desc, sizeof(char) * n);
      strcpy(sq->desc, desc);
      sq->dalloc = n;
    } 
  else 
    {
      sq->dalloc   = eslSQ_DESCCHUNK;
      ESL_ALLOC(sq->desc, sizeof(char) * sq->dalloc);    
      sq->desc[0] = '\0';
    }

  if (acc != NULL) 
    {
      n = strlen(acc)+1;
      ESL_ALLOC(sq->acc, sizeof(char) * n);
      strcpy(sq->acc, acc);
      sq->aalloc = n;
    } 
  else 
    {
      sq->aalloc   = eslSQ_ACCCHUNK;
      ESL_ALLOC(sq->acc,  sizeof(char) * sq->aalloc);
      sq->acc[0] = '\0';
    }

  sq->doff = -1;
  sq->roff = -1;
  return sq;

 ERROR:
  esl_sq_Destroy(sq);
  return NULL;
}
/*----------------- end, internal functions ---------------------*/


/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/

/*--------------------- end, unit tests -------------------------*/

/*****************************************************************
 * 6. Test driver.
 *****************************************************************/


/*-------------------- end, test driver -------------------------*/

/*****************************************************************
 * 7. Example.
 *****************************************************************/




/*------------------ end, example driver ------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
