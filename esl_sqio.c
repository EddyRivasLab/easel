/* Unaligned sequence file i/o.
 * 
 * Sections:
 *    1. The ESL_SQ object API.
 *    2. The ESL_SQFILE object API.
 *    3. Digitized sequences. (Alphabet augmentation required.)
 *    4. The sequence input/output API.
 *    5. Internal functions.
 *    6. Test and example code.
 * 
 * Shares remote evolutionary homology with Don Gilbert's seminal,
 * public domain ReadSeq package. Last common ancestor was 1991 or so.
 * Vestiges of that history still remain in the design. Thanks Don!
 * 
 * BUG:  SRE, Tue Dec 11 14:52:15 2007: 
 * esl_sqio_Read() in digital mode will *not* throw an illegal character
 * exception; for example, put an L in a DNA seq and see it return eslOK.
 * 
 * SRE, Thu Feb 17 17:45:51 2005
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
#include "esl_msa.h"		/* msa aug adds ability to read MSAs as unaligned seqs  */
#endif
#include "esl_sqio.h"

/* Shared parts of text/digital creation functions */
static ESL_SQ *sq_create(int do_digital);
static ESL_SQ *sq_create_from(char *name, char *desc, char *acc, char *ss);

/* Generic functions for line-based parsers.
 */
static int is_blankline(char *s);
static int loadline(ESL_SQFILE *sqfp);
static int addseq(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int generic_readseq(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int set_name(ESL_SQ *sq, char *s, char *delim);
static int set_accession(ESL_SQ *sq, char *s, char *delim);
static int append_description(ESL_SQ *sq, char *s, char *delim);

/* EMBL format; also Uniprot, TrEMBL
 */
static void config_embl(ESL_SQFILE *sqfp);
static int  read_embl(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_embl(char *buf);

/* Genbank format; also DDBJ
 */
static void config_genbank(ESL_SQFILE *sqfp);
static int  read_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_genbank(char *buf);

/* FASTA format: uses faster character-based fread() i/o.
 */
static void config_fasta(ESL_SQFILE *sqfp);
static int  read_fasta(ESL_SQFILE *sqfp, ESL_SQ *s);
static int  write_fasta(FILE *fp, ESL_SQ *s);
static int  check_buffers(FILE *fp, off_t *boff, char *buf, int *nc, int *pos, 
			 char **s, ESL_DSQ **dsq, int i, int *slen);
#ifdef eslAUGMENT_ALPHABET
static int  write_digital_fasta(FILE *fp, ESL_SQ *s);
#endif

/* Optional MSA<->sqio interoperability */
#ifdef eslAUGMENT_MSA
static int convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa);
#endif




/*****************************************************************
 * 1. Routines for dealing with the ESL_SQ object.
 *****************************************************************/ 

/* Create and CreateDigital() (see below) are almost identical, so
 * their shared guts are here:
 */
static ESL_SQ *
sq_create(int do_digital)
{
  int status;
  ESL_SQ *sq = NULL;

  ESL_ALLOC(sq, sizeof(ESL_SQ));

  sq->name     = NULL;
  sq->acc      = NULL;
  sq->desc     = NULL;
  sq->seq      = NULL;
  sq->ss       = NULL;	/* secondary structure input currently unimplemented */
  sq->dsq      = NULL;	
  sq->optmem   = NULL;	/* this stays NULL unless we Squeeze() the structure */
  sq->flags    = 0;
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

/* Function:  esl_sq_Create()
 * Incept:    SRE, Thu Dec 23 11:57:00 2004 [Zaragoza]
 *
 * Purpose:   Creates an empty <ESL_SQ> sequence object, with
 *            internal fields allocated to reasonable initial sizes. 
 *            
 * Args:      (void)
 *
 * Returns:   a pointer to the new <ESL_SQ>. Caller frees this with
 *            <esl_sq_Destroy()>.
 *
 * Throws:    NULL if allocation fails.
 */
ESL_SQ *
esl_sq_Create(void)
{
  return sq_create(FALSE);
}

/* CreateFrom and CreateDigitalFrom() (see below) are almost identical, so
 * their shared guts are here:
 */
static ESL_SQ *
sq_create_from(char *name, char *desc, char *acc, char *ss)
{
  int status;
  ESL_SQ *sq = NULL;
  int  n;

  if (name == NULL) ESL_XEXCEPTION(eslEINVAL, "must provide seq name");

  ESL_ALLOC(sq, sizeof(ESL_SQ));
  sq->name   = NULL;
  sq->acc    = NULL;
  sq->desc   = NULL;
  sq->seq    = NULL;
  sq->dsq    = NULL;
  sq->ss     = NULL;
  sq->optmem = NULL;
  sq->flags  = 0;
  
  n = strlen(name)+1;
  ESL_ALLOC(sq->name, sizeof(char) * n);
  strcpy(sq->name, name);
  sq->nalloc = n;
  
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

  if (ss != NULL) 
    {
      n = strlen(ss)+1;
      if (n != sq->salloc) ESL_XEXCEPTION(eslEINVAL, "ss, seq lengths mismatch");
      ESL_ALLOC(sq->ss, sizeof(char) * n);
      strcpy(sq->ss, ss);
    } 

  sq->doff = -1;
  sq->roff = -1;
  return sq;

 ERROR:
  esl_sq_Destroy(sq);
  return NULL;
}

/* Function:  esl_sq_CreateFrom()
 * Incept:    SRE, Wed Mar 22 09:17:04 2006 [St. Louis]
 *
 * Purpose:   Create a new <ESL_SQ> object from elemental data.
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
 *            The object is growable; you can use <esl_sq_Reuse()>
 *            on it.
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

  if ((sq     = sq_create_from(name, desc, acc, ss)) == NULL)  goto ERROR;
  if ((status = esl_strdup(seq, n, &(sq->seq)))      != eslOK) goto ERROR;
  sq->n      = n;
  sq->salloc = n;
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
 *            to <esl_sq_Grow()> in <ret_nsafe>.
 *            
 *            The terminal <NUL> or sentinel count as a 'residue'
 *            for this: that is, you may need to call <esl_sq_Grow()>
 *            before terminating a new sequence.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. In this case, the
 *            original <sq> is untouched, and <*ret_nsafe> is returned
 *            as 0.
 *
 * Xref:      STL11/125.
 */
int
esl_sq_Grow(ESL_SQ *sq, int *ret_nsafe)
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
      do { nsafe += new; new  *=2; } while (nsafe < 1);
      
      if (sq->seq != NULL) ESL_RALLOC(sq->seq, tmp, new * sizeof(ESL_DSQ));	/* text */
      else                 ESL_RALLOC(sq->dsq, tmp, new * sizeof(ESL_DSQ));	/* digital */
      sq->salloc = new;
    }

  if (ret_nsafe != NULL) *ret_nsafe = nsafe;
  return eslOK;

 ERROR:
  if (ret_nsafe != NULL) *ret_nsafe = 0;
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
	sq->salloc = n+1;
      }
    }
  else				/* digital mode */
    {
      if (n+2 > sq->salloc) {
	ESL_RALLOC(sq->dsq, tmp, (n+2) * sizeof(ESL_DSQ));
      }
    }
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_sq_Copy()
 * Synopsis:  Make a copy of an <ESL_SQ>
 * Incept:    SRE, Sun Feb 24 17:59:24 2008 [UA5315 to St. Louis]
 *
 * Purpose:   Copies a source sequence object <src> into 
 *            destination sequence object <dst>.
 *            
 *            The two objects don't have to be matched as far as
 *            text/digital mode; if mismatched, appropriate
 *            text/digital conversion will be done.
 *            
 *            The destination sequence <sq> must be allocated, with
 *            enough room to hold the <src> sequence (see
 *            <esl_sq_GrowTo()>). However, <dst> does not need to be
 *            allocated for other info, such as the name or
 *            description. These fields will be reallocated if
 *            necessary.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sq_Copy(const ESL_SQ *src, ESL_SQ *dst)
{
  int status;

  if ((status = esl_sq_SetName     (dst, src->name)) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(dst, src->acc))  != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (dst, src->desc)) != eslOK) goto ERROR;
  if ((status = esl_sq_GrowTo      (dst, src->n))    != eslOK) goto ERROR;

  if (! (src->flags & eslSQ_DIGITAL) && ! (dst->flags & eslSQ_DIGITAL))
    strcpy(dst->seq, src->seq);
#ifdef eslAUGMENT_ALPHABET
  else if (! (src->flags & eslSQ_DIGITAL) && (dst->flags & eslSQ_DIGITAL))
    {
      if ((status = esl_abc_Digitize(dst->abc, src->seq, dst->dsq)) != eslOK) goto ERROR;
    }
  else if ((src->flags & eslSQ_DIGITAL) &&  ! (dst->flags & eslSQ_DIGITAL))
    {
      if ((status = esl_abc_Textize(src->abc, src->dsq, src->n, dst->seq)) != eslOK) goto ERROR;
    }
  else 
    {
      if (src->abc->type != dst->abc->type) 
	ESL_XEXCEPTION(eslEINCOMPAT, "seq objects involved in Copy differ in digital alphabet");
      if ((status = esl_abc_dsqcpy(src->dsq, src->n, dst->dsq)) != eslOK) goto ERROR;
    }
#endif
      
  dst->n     = src->n;
  dst->roff  = src->roff;
  dst->doff  = src->doff;
  /* don't copy flags, because all it holds is digital status, and dst retains preexisting status. */
  /* don't copy allocations (nalloc, etc); dst knows its own memory */
  /* and don't copy optmem; dst is assumed to be growable here */
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
  if (sq->dsq != NULL) sq->dsq[0] = '\0';
  if (sq->ss  != NULL) sq->ss[0]  = '\0';
  sq->n    = 0;
  sq->doff = -1;
  sq->roff = -1;
  return eslOK;
}


/* Function:  esl_sq_Squeeze()
 * Incept:    SRE, Sat Dec 25 04:30:47 2004 [Zaragoza]
 *
 * Purpose:   Given a <sq>, optimize its memory usage.
 *            The <sq> may never again be used for sequence
 *            input, because its dynamic buffers are 
 *            destroyed by this call.
 *
 *            When a sequence is input, data spaces are
 *            dynamically allocated to allow unlimited
 *            lengths. This results in somewhat inefficient 
 *            memory usage (up to 50\%). If an application
 *            is reading through a sequence database one
 *            seq at a time, this is acceptable, but if
 *            an app needs to read in a lot of seqs at
 *            once, it may care about optimizing memory.
 *            This function perfectly reallocates and
 *            copies the internal data, and free's the
 *            dynamic input buffers. 
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if reallocation fails.
 */
int
esl_sq_Squeeze(ESL_SQ *sq)
{
  int      status;
  int      nlen, alen, dlen, len;
  char    *name, *acc, *desc, *seq, *ss;
  ESL_DSQ *dsq;

  nlen = strlen(sq->name);
  alen = strlen(sq->acc);
  dlen = strlen(sq->desc);

  len = nlen + alen + dlen + 3; 
  if (sq->seq != NULL) len += sq->n+1;
  if (sq->ss  != NULL) len += sq->n+1;

  ESL_ALLOC(sq->optmem, sizeof(char) * len);
  
  len  = 0;
  name = sq->optmem+len; memcpy(name, sq->name, nlen+1);  len+=nlen+1;
  acc  = sq->optmem+len; memcpy(acc,  sq->acc,  alen+1);  len+=alen+1;
  desc = sq->optmem+len; memcpy(desc, sq->desc, dlen+1);  len+=dlen+1;
  if (sq->seq != NULL) 
    { seq  = sq->optmem+len; memcpy(seq, sq->seq,  sq->n+1); len+=sq->n+1;}
  if (sq->ss  != NULL)
    { ss   = sq->optmem+len; memcpy(ss,  sq->ss,   sq->n+1); len+=sq->n+1;}

  if (sq->dsq != NULL) {
    ESL_ALLOC(dsq, sizeof(ESL_DSQ) * sq->n+2);
    memcpy(dsq, sq->dsq, sizeof(ESL_DSQ) * sq->n+2);
  }

  free(sq->name); sq->nalloc = 0; sq->name = name;
  free(sq->acc);  sq->aalloc = 0; sq->acc  = acc;
  free(sq->desc); sq->dalloc = 0; sq->desc = desc;
  if (sq->seq != NULL) { free(sq->seq); sq->salloc = 0; sq->seq  = seq; }
  if (sq->dsq != NULL) { free(sq->dsq); sq->salloc = 0; sq->dsq  = dsq; }
  if (sq->ss  != NULL) { free(sq->ss);  sq->ss  = ss; }
  return eslOK;

 ERROR:
  return status;
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

  if (sq->optmem != NULL)
    { 
      free(sq->optmem); 
    }
  else
    {
      if (sq->name   != NULL) free(sq->name);  
      if (sq->acc    != NULL) free(sq->acc);   
      if (sq->desc   != NULL) free(sq->desc);  
      if (sq->seq    != NULL) free(sq->seq);   
      if (sq->dsq    != NULL) free(sq->dsq);   
      if (sq->ss     != NULL) free(sq->ss);    
    }
  free(sq);
  return;
}
/*----------------- end of ESL_SQ object functions -----------------*/


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
  int   n;
  void *tmp;
  int   status;

  va_start(argp, name);
  if ((n = vsnprintf(sq->name, sq->nalloc, name, argp)) > sq->nalloc)
    {
      ESL_RALLOC(sq->name, tmp, sizeof(char) * n); 
      sq->nalloc = n;
      vsnprintf(sq->name, sq->nalloc, name, argp);
    }
  va_end(argp);
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_sq_SetAccession()
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
  va_list argp;
  int     n;
  void   *tmp;
  int     status;

  va_start(argp, acc);
  if ((n = vsnprintf(sq->acc, sq->aalloc, acc, argp)) > sq->aalloc)
    {
      ESL_RALLOC(sq->acc, tmp, sizeof(char) * n); 
      sq->aalloc = n;
      vsnprintf(sq->acc, sq->aalloc, acc, argp);
    }
  va_end(argp);
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_sq_SetDesc()
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
  va_list argp;
  int     n;
  void   *tmp;
  int     status;

  va_start(argp, desc);
  if ((n = vsnprintf(sq->desc, sq->dalloc, desc, argp)) > sq->dalloc)
    {
      ESL_RALLOC(sq->desc, tmp, sizeof(char) * n); 
      sq->nalloc = n;
      vsnprintf(sq->desc, sq->dalloc, desc, argp);
    }
  va_end(argp);
  return eslOK;

 ERROR:
  return status;
}



/* Function:  esl_sq_CAddResidue()
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
 *            use addseq() instead.
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


/*****************************************************************
 * Section 2. Routines for dealing with the ESL_SQFILE object.
 *****************************************************************/ 
/* Function:  esl_sqfile_Open()
 * Incept:    SRE, Thu Feb 17 08:22:16 2005 [St. Louis]
 *
 * Purpose:   Open a sequence file <filename> for sequential reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The format of the file is asserted to be <format> (for
 *            example, <eslSQFILE_FASTA>).
 *            If <format> is <eslSQFILE_UNKNOWN> then format
 *            autodetection is invoked. 
 *            
 *            There are two special cases for <filename>. If
 *            <filename> is "-", the sequence data are read from a
 *            <STDIN> pipe. If <filename> ends in ".gz", the file is assumed
 *            to be compressed with <gzip>, and it is opened by a pipe
 *            from <gzip -dc>; this only works on POSIX-compliant
 *            systems that have pipes (specifically, the POSIX.2
 *            popen() call); this code is included only if
 *            <HAVE_POPEN> is defined at compile time. To use either
 *            of these abilities, <format> must be defined, not unknown;
 *            format autodetection requires a two-pass parse on a rewindable
 *            stream, but pipes are not rewindable.
 *            
 *            If <env> is non-NULL, it is the name of an environment
 *            variable that contains a colon-delimited list of
 *            directories in which we may find this <filename>.
 *            For example, if we had 
 *            <setenv BLASTDB /nfs/db/blast-db:/nfs/db/genomes/>
 *            in the environment, a database search application
 *            could pass "BLASTDB" as <env>.
 *            
 * Returns:   <eslOK> on success, and <*ret_sqfp> points to a new
 *            open <ESL_SQFILE>. Caller deallocates this object with
 *            <esl_sqfile_Close()>. 
 *            
 *            Returns <eslENOTFOUND> if <filename> can't be opened.
 *            Returns <eslEFORMAT> if the file is empty, or if
 *            autodetection is attempted and the format can't be
 *            determined.  Returns <eslEINVAL> if autodetection is
 *            attempted on a stdin or gunzip pipe.  On any of these error
 *            conditions, <*ret_sqfp> is returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
int
esl_sqfile_Open(char *filename, int format, char *env, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp    = NULL;
  char       *envfile = NULL;
  int         status;		/* return status from an ESL call */
  int         n;
  ESL_DSQ     x;

  /* Allocate and initialize the structure to base values;
   * only format is set correctly, though.
   */
  ESL_ALLOC(sqfp, sizeof(ESL_SQFILE));
  *ret_sqfp        = NULL;
  sqfp->fp         = NULL;
  sqfp->filename   = NULL;
  sqfp->ssifile    = NULL;
  sqfp->format     = format;
  sqfp->do_gzip    = FALSE;
  sqfp->do_stdin   = FALSE;
  sqfp->errbuf[0]  = '\0';
  sqfp->buf        = NULL;
  sqfp->boff       = 0;
  sqfp->balloc     = 0;
  sqfp->nc         = 0;
  sqfp->pos        = 0;
  sqfp->linenumber = 1;
  sqfp->rpl        = -1;	/* -1=unset */
  sqfp->bpl        = -1;	/* -1=unset */
  sqfp->lastrpl    = -1;	/* -1=unset */
  sqfp->lastbpl    = 0;		/* -1=unset */
  sqfp->sq_cache   = NULL;	/* only used if we GuessAlphabet() later */
#ifdef eslAUGMENT_MSA
  sqfp->afp        = NULL;
  sqfp->msa        = NULL;
#endif /*eslAUGMENT_MSA*/

  /* Open the file. It may either be in the current directory,
   * or in a directory indicated by the <env> argument. We have
   * to construct the SSI filename accordingly. For normal operation
   * (no pipes from stdin, gzip), this section opens the sqfp->fp,
   * stores the filename in sqfp->filename, and sets sqfp->ssifile to
   * the name of the SSI file that we should look for for this seq file.
   * 
   * stdin special case is handled here. fp is stdin pipe; filename
   * is [STDIN]; ssifile left NULL.
   */
  if (strcmp(filename, "-") == 0) /* stdin */
    {
      if ((status = esl_strdup("[STDIN]", -1, &(sqfp->filename))) != eslOK) goto ERROR;
      sqfp->fp       = stdin;
      sqfp->do_stdin = TRUE;
    }
  else
    {
      n = strlen(filename);  

      /* Check the current working directory first.
       */
      if ((sqfp->fp = fopen(filename, "r")) != NULL)
	{
	  if ((status = esl_FileNewSuffix(filename, "ssi", &(sqfp->ssifile))) != eslOK) goto ERROR;
	  if ((status = esl_strdup(filename, n, &(sqfp->filename)))           != eslOK) goto ERROR;
	}
      /* then the env variable.
       */
      else if (env != NULL &&
	       esl_FileEnvOpen(filename, env, &(sqfp->fp), &envfile)== eslOK)
	{
	  if ((status = esl_FileNewSuffix(envfile, "ssi", &(sqfp->ssifile))) != eslOK) goto ERROR;
	  if ((status = esl_strdup(envfile, -1, &(sqfp->filename)))          != eslOK) goto ERROR;
	  free(envfile); envfile = NULL;
	}
      else
	{ status = eslENOTFOUND; goto ERROR;}
    }


  /* Deal with the .gz special case.
   * 
   * To popen(), "success" means it found and executed gzip -dc.  If
   * gzip -dc doesn't find our file, popen() still blithely returns
   * success, so we have to be sure the file exists. Fortunately, we
   * already know that, because we fopen()'ed it as a normal file in
   * the section above.
   * 
   * For a .gz, close the fp we've got, and reopen it as a pipe from
   * gzip -dc w/ popen(). (But if HAVE_POPEN isn't defined, then a .gz
   * file is treated as a normal file.)
   *
   * After this section, fp, filename, ssifile, do_gzip, and do_stdin are
   * all set correctly in the sqfile object.
   */                           
#ifdef HAVE_POPEN
  n = strlen(sqfp->filename);
  if (n > 3 && strcmp(sqfp->filename+n-3, ".gz") == 0) 
    {
      char *cmd;

      fclose(sqfp->fp);

      ESL_ALLOC(cmd, sizeof(char) * (n+1+strlen("gzip -dc ")));
      sprintf(cmd, "gzip -dc %s", sqfp->filename);
      if ((sqfp->fp = popen(cmd, "r")) == NULL)	{ status = eslENOTFOUND; goto ERROR; }
      if ((status = esl_strdup(sqfp->filename, n, &(sqfp->filename))) != eslOK) goto ERROR;
      sqfp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN*/


  /* Init the input map. 
   */
  for (x = 0; x < 128; x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;

  /* If we don't know the format yet, autodetect it now.
   */
  if (sqfp->format == eslSQFILE_UNKNOWN)
    {
      if (sqfp->do_stdin || sqfp->do_gzip)   { status = eslEINVAL;  goto ERROR; }

      sqfp->format = esl_sqio_WhatFormat(sqfp->fp);

      if (sqfp->format == eslSQFILE_UNKNOWN) { status = eslEFORMAT; goto ERROR; }
    }


  /* Configure the <sqfp> for this specific format.
   */
  switch (sqfp->format) {
  case eslSQFILE_EMBL:     config_embl(sqfp);    break;
  case eslSQFILE_UNIPROT:  config_embl(sqfp);    break;
  case eslSQFILE_GENBANK:  config_genbank(sqfp); break;
  case eslSQFILE_DDBJ:     config_genbank(sqfp); break;
  case eslSQFILE_FASTA:    config_fasta(sqfp);   break;

#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM:
    sqfp->linenumber = 0;	/* line-oriented input */
    sqfp->is_linebased = TRUE;
    sqfp->addfirst = FALSE;	/* no-op for msa's */
    sqfp->addend   = FALSE;	/* no-op for msa's */
    sqfp->eof_is_ok= FALSE;	/* no-op for msa's */
    sqfp->endTest  = NULL;	/* no-op for msa's */
    ESL_ALLOC(sqfp->afp, sizeof(ESL_MSAFILE));
    sqfp->afp->f          = sqfp->fp;
    sqfp->afp->fname      = sqfp->filename;
    sqfp->afp->linenumber = sqfp->linenumber;
    sqfp->afp->errbuf[0]  = '\0';
    sqfp->afp->buf        = NULL;
    sqfp->afp->buflen     = 0;
    sqfp->afp->do_gzip    = sqfp->do_gzip;
    sqfp->afp->do_stdin   = sqfp->do_stdin;
    sqfp->afp->format     = sqfp->format;
    sqfp->afp->do_digital = FALSE;
    sqfp->afp->abc        = NULL;
    sqfp->afp->ssi        = NULL;
    sqfp->afp->msa_cache  = NULL;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  /* Preload the first line or chunk of the file into buf.
   */
  if (! esl_sqio_IsAlignment(sqfp->format))
    {
      if (sqfp->is_linebased)
	{
	  sqfp->linenumber = 0;
	  status = loadline(sqfp);
	  if (status == eslEOF)     { status = eslEFORMAT; goto ERROR; }
	  else if (status != eslOK) { goto ERROR; }
	}
      else
	{
	  sqfp->linenumber = 1;
	  sqfp->balloc = eslREADBUFSIZE;
	  ESL_ALLOC(sqfp->buf, sizeof(char) * sqfp->balloc);
	  sqfp->boff = 0;
	  sqfp->nc   = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
	  if (ferror(sqfp->fp)) { status = eslEFORMAT; goto ERROR; }
	}
    }

  if (envfile != NULL) free(envfile);
  *ret_sqfp = sqfp;
  return eslOK;

 ERROR:
  if (envfile != NULL) free(envfile);
  esl_sqfile_Close(sqfp); 
  *ret_sqfp = NULL;
  return status;
}


/* Function:  esl_sqfile_Close()
 * Incept:    SRE, Thu Dec 23 13:19:43 2004 [St. Louis]
 *
 * Purpose:   Closes an open <sqfp>.
 *
 * Returns:   (void).
 */
void
esl_sqfile_Close(ESL_SQFILE *sqfp)
{
  if (sqfp == NULL) return;

#ifdef HAVE_POPEN
  if (sqfp->do_gzip)          pclose(sqfp->fp);
#endif
  if (! sqfp->do_stdin && sqfp->fp != NULL) fclose(sqfp->fp);
  if (sqfp->filename != NULL) free(sqfp->filename);
  if (sqfp->ssifile  != NULL) free(sqfp->ssifile);
  if (sqfp->buf      != NULL) free(sqfp->buf);
  if (sqfp->sq_cache != NULL) esl_sq_Destroy(sqfp->sq_cache);

#ifdef eslAUGMENT_MSA
  if (sqfp->afp      != NULL) 
    { /* Because we copied info from the seqfile object to
       * create the msafile object, we can't just close the 
       * msafile, or we'd end up w/ double fclose()/free()'s.
       */
      if (sqfp->afp->buf != NULL) free(sqfp->afp->buf);
      free(sqfp->afp);
    }
  if (sqfp->msa      != NULL) esl_msa_Destroy(sqfp->msa);
#endif /*eslAUGMENT_MSA*/

  free(sqfp);
  return;
}

/*------------------- ESL_SQFILE open/close -----------------------*/

/*****************************************************************
 * Section 3. Digitized sequences (ALPHABET augmentation required)
 *****************************************************************/ 

#ifdef eslAUGMENT_ALPHABET
/* Function:  esl_sq_CreateDigital()
 * Incept:    SRE, Tue Jan  9 16:42:35 2007 [Janelia]
 *
 * Purpose:   Same as <esl_sq_Create()>, except the returned sq is configured
 *            for a digital sequence using internal alphabet <abc>, rather than
 *            a text sequence. Creates an empty digital <ESL_SQ> sequence 
 *            object, with internal fields allocated to reasonable initial sizes.
 *            Additionally, the <eslSQ_DIGITAL> flag is raised. 
 * 
 *            The first byte of the digital sequence 
 *            (<s->dsq>, internally) is initialized to a sentinel,
 *            but the terminal sentinel byte is the caller's
 *            responsibility. Sequence generation or i/o routines
 *            will generally handle this.
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
  if ((s = sq_create(TRUE)) != NULL)
    {
      s->abc    = abc;
      s->dsq[0] = eslDSQ_SENTINEL;
      s->flags |= eslSQ_DIGITAL;
    }
  return s;
}

/* Function:  esl_sq_CreateDigitalFrom()
 * Incept:    EPN, Fri Aug 24 13:38:56 2007
 *
 * Purpose:   Create a new <ESL_SQ> object from elemental data;
 *            Same as <esl_sq_CreateFrom> except takes a digital <ESL_DSQ *dsq>
 *            instead of a text <char *seq> as the sequence to copy.
 *            Originally written to ease reverse complementing
 *            database sequences for Infernal searches.
 *            
 *            Makes copies of all data. Caller is still
 *            responsible for memory of name, seq, etc.
 *            
 *            <ss> is an optional alphabetic secondary structure 
 *            annotation string. If provided, its length must match
 *            the length of <seq>.
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

  if((sq = sq_create_from(name, desc, acc, ss)) == NULL) goto ERROR;
  sq->n = (L == -1) ? esl_abc_dsqlen(dsq) : L;
  if ((status = esl_abc_dsqdup(dsq, sq->n, &(sq->dsq))) != eslOK) goto ERROR;
  sq->salloc = sq->n;
  sq->abc   =  abc;
  sq->flags |= eslSQ_DIGITAL;
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
 *            <abc> reference, and the <eslSEQ_DIGITAL> flag is raised
 *            in <flags>.
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

  /* Contract checks
   */
  if (sq->seq   == NULL)           ESL_EXCEPTION(eslEINVAL, "sq has no text sequence");
  if (sq->dsq   != NULL)           ESL_EXCEPTION(eslEINVAL, "sq already has digital sequence");
  if (sq->flags & eslSQ_DIGITAL)   ESL_EXCEPTION(eslEINVAL, "sq is flagged as digital");

  /* Validate before we convert. Then we can leave the <seq> untouched if
   * the sequence contains invalid characters.
   */
  if (esl_abc_ValidateSeq(abc, sq->seq, sq->n, NULL) != eslOK) 
      return eslEINVAL;

  /* Convert, free seq.
   */
  ESL_ALLOC(sq->dsq, (sq->n+2) * sizeof(ESL_DSQ));
  status = esl_abc_Digitize(abc, sq->seq, sq->dsq);
  if (status != eslOK) goto ERROR;
  free(sq->seq);
  sq->seq = NULL;

  sq->abc   =  abc;
  sq->flags |= eslSQ_DIGITAL;
  return eslOK;

 ERROR:
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
 *            sq's <abc> digital alphabet reference is nullified, and 
 *            the <eslSQ_DIGITAL> flag is dropped in <flags>.
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

  /* Contract checks
   */
  if (sq->dsq   == NULL)               ESL_EXCEPTION(eslEINVAL, "sq has no digital sequence");
  if (sq->seq   != NULL)               ESL_EXCEPTION(eslEINVAL, "sq already has text sequence");
  if (! (sq->flags & eslSQ_DIGITAL))   ESL_EXCEPTION(eslEINVAL, "sq is not flagged as digital");
  if (sq->abc   == NULL)               ESL_EXCEPTION(eslEINVAL, "sq has no digital alphabet");

  /* Convert, free dsq. 
   */
  ESL_ALLOC(sq->seq, (sq->n+2) * sizeof(char));
  status = esl_abc_Textize(sq->abc, sq->dsq, sq->n, sq->seq);
  if (status != eslOK) goto ERROR;
  free(sq->dsq);
  sq->dsq = NULL;
  
  sq->abc    = NULL;           /* nullify reference (caller still owns real abc) */
  sq->flags &= ~eslSQ_DIGITAL; /* drop the flag */
  return eslOK;

 ERROR:
  return status;
}
#endif /* eslAUGMENT_ALPHABET */


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

/* Function:  esl_sqfile_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open <ESL_SQFILE>
 * Incept:    SRE, Sun Feb 24 17:14:55 2008 [UA5315 to St. Louis]
 *
 * Purpose:   After opening <sqfp>, attempt to guess what alphabet
 *            its sequences are in, by inspecting the first sequence
 *            in the file, and return this alphabet type in <*ret_type>.
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set to <eslDNA>,
 *            <eslRNA>, or <eslAMINO>.
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to 
 *            <eslUNKNOWN> if the first sequence (or alignment)
 *            in the file contains no more than ten residues total,
 *            or if its alphabet cannot be guessed (i.e. it contains
 *            IUPAC degeneracy codes, but no amino acid specific
 *            residues).
 *            
 *            Returns <eslEFORMAT> if a parse error is encountered in
 *            trying to read the sequence file. <sqfp->errbuf> is set
 *            to a useful error message if this occurs,
 *            <sqfp->linenumber> is the line on which the error
 *            occurred, and <*ret_type> is set to <eslUNKNOWN>.
 *            
 *            Returns <eslENODATA> and sets <*ret_type> to <eslUNKNOWN>
 *            if the file appears to be empty.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslEINCONCEIVABLE> on unimaginable internal errors.
 */
int
esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)
{
  ESL_SQ *sq = NULL;
  int     status;

  /* Special case: for MSA files, we already have first MSA cached. */
  if (esl_sqio_IsAlignment(sqfp->format)) return esl_msa_GuessAlphabet(sqfp->msa, ret_type);

  /* Special case: already something cached; GuessAlphabet() was already called? */
  if (sqfp->sq_cache != NULL) return esl_sq_GuessAlphabet(sqfp->sq_cache, ret_type);

  /* Read and cache the first seq, and guess alphabet based on that.
   * This is risky - the first seq might be short/atypical, and fool us about
   * the rest of the file.
   */
  if ((sq = esl_sq_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEOF) { status = eslENODATA; goto ERROR; }
  else if (status != eslOK)  goto ERROR; 
  
  sqfp->sq_cache = sq;
  return esl_sq_GuessAlphabet(sqfp->sq_cache, ret_type);

 ERROR:
  esl_sq_Destroy(sq);
  sqfp->sq_cache = NULL;
  *ret_type      = eslUNKNOWN;
  return status;
}

/*-------------------- end of digital sequence functions --------------------*/


/*****************************************************************
 * Section 4. Sequence i/o API
 *****************************************************************/ 

/* Function:  esl_sqio_Read()
 * Incept:    SRE, Thu Feb 17 14:24:21 2005 [St. Louis]
 *
 * Purpose:   Reads the next sequence from open sequence file <sqfp> into 
 *            <sq>. Caller provides an allocated and initialized <s>, which
 *            will be internally reallocated if its space is insufficient.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character; the line number that the parse
 *            error occurs on is in <sqfp->linenumber>, and an informative
 *            error message is placed in <sqfp->errbuf>. 
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_sqio_Read(ESL_SQFILE *sqfp, ESL_SQ *s)
{
  int status;

  /* Special case: we already have a sequence cached. */
  if (sqfp->sq_cache != NULL) 
    {
      status = esl_sq_Copy(sqfp->sq_cache, s);
      esl_sq_Destroy(sqfp->sq_cache);
      sqfp->sq_cache = NULL;
      return status;
    }

  switch (sqfp->format) {
  case eslSQFILE_FASTA:    
    status = read_fasta(sqfp, s);   break;
    
  case eslSQFILE_EMBL:     
  case eslSQFILE_UNIPROT:
    status = read_embl(sqfp, s);    break;

  case eslSQFILE_GENBANK:  
  case eslSQFILE_DDBJ:
    status = read_genbank(sqfp, s); break;
    
#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM:
    if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
      {				/* load a new alignment */
	esl_msa_Destroy(sqfp->msa);
	status = esl_msa_Read(sqfp->afp, &(sqfp->msa));
	if (status == eslEFORMAT)
	  { /* oops, a parse error; upload the error info from afp to sqfp */
	    sqfp->linenumber = sqfp->afp->linenumber;
	    strcpy(sqfp->errbuf, sqfp->afp->errbuf); /* errbufs same size! */ 
	    return eslEFORMAT;
	  }
	if (status != eslOK) return status;
	sqfp->idx = 0;
      }
    /* grab next seq from alignment */
    /* this is inefficient; it goes via a temporarily allocated copy of the sequence */
    status = esl_sq_FetchFromMSA(sqfp->msa, sqfp->idx, &tmpsq);
    esl_sq_GrowTo(s, tmpsq->n);
    esl_sq_Copy(tmpsq, s);
    esl_sq_Destroy(tmpsq);

    sqfp->idx++;
    status = eslOK;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  return status;
}


/* Function:  esl_sqio_Write()
 * Incept:    SRE, Fri Feb 25 16:10:32 2005 [St. Louis]
 *
 * Purpose:   Write sequence <s> to an open FILE <fp> in 
 *            file format <format>. 
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sqio_Write(FILE *fp, ESL_SQ *s, int format)
{
  int status;
#ifdef eslAUGMENT_MSA
  ESL_MSA *msa;
#endif

  if (s->seq != NULL) 		/* text mode */
    {
      switch (format) {
      case eslSQFILE_FASTA: status = write_fasta(fp, s); break;

#ifdef eslAUGMENT_MSA
      case eslMSAFILE_STOCKHOLM:
      case eslMSAFILE_PFAM:
	/* For writing single sequences in "alignment" format,
	 * we convert the SQ object to an MSA, then write using
	 * the MSA API.
	 */
	if ((status = convert_sq_to_msa(s, &msa)) != eslOK) return status;
	status = esl_msa_Write(fp, msa, format);
	esl_msa_Destroy(msa);
	break;
#endif /* msa augmentation */

      default: 
	ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
      }
    }
  else				/* digital mode */
    {
#if defined (eslAUGMENT_ALPHABET)
      switch (format) {
      case eslSQFILE_FASTA: status = write_digital_fasta(fp, s); break;
      default:              ESL_EXCEPTION(eslEINCONCEIVABLE, "only supporting fasta for digital output currently");
      }
#else
      ESL_EXCEPTION(eslEINCONCEIVABLE, "whoops, how did I get a digital sequence?");
#endif
    }
  return status;
}

/* Function:  esl_sqio_WhatFormat()
 * Incept:    SRE, Mon Jun 20 19:07:44 2005 [St. Louis]
 *
 * Purpose:   Determine the format of a (rewindable) open file <fp>;
 *            return the appropriate code, or <eslSQFILE_UNKNOWN> 
 *            if the file format can't be determined.
 *            
 *            Rewinds the <fp> before returning.
 *
 * Returns:   File format code, such as <eslSQFILE_FASTA>.
 */
int
esl_sqio_WhatFormat(FILE *fp)
{
  char buf[eslREADBUFSIZE];
  int fmt;

  /* get first nonblank line */
  do {
    if (fgets(buf, eslREADBUFSIZE, fp) == NULL) 
      { rewind(fp); return eslSQFILE_UNKNOWN; }
  } while (is_blankline(buf));

  /* formats that can be determined from the first line:
   */
  if      (*buf == '>')                                       fmt = eslSQFILE_FASTA;
  else if (strncmp(buf, "ID   ", 5)    == 0)                  fmt = eslSQFILE_EMBL;
  else if (strncmp(buf, "LOCUS   ", 8) == 0)                  fmt = eslSQFILE_GENBANK;
  else if (strstr(buf, "Genetic Sequence Data Bank") != NULL) fmt = eslSQFILE_GENBANK;
#ifdef eslAUGMENT_MSA
  else if (strncmp(buf, "# STOCKHOLM", 11) == 0)              fmt = eslMSAFILE_STOCKHOLM;
#endif
  else                                                        fmt = eslSQFILE_UNKNOWN;

  rewind(fp);
  return fmt;
}

/* Function:  esl_sqio_FormatCode()
 * Incept:    SRE, Sun Feb 27 09:18:36 2005 [St. Louis]
 *
 * Purpose:   Given <fmtstring>, return format code.  For example, if
 *            <fmtstring> is "fasta", returns <eslSQFILE_FASTA>. Returns 
 *            <eslSQFILE_UNKNOWN> if <fmtstring> doesn't exactly match a 
 *            known format.
 *            
 *            The match is aggressively case insensitive: the <fmtstring>
 *            is converted to all upper case. (We would use strcasecmp(),
 *            but that isn't ANSI C.)
 *            
 *            When augmented by msa, then alignment file formats
 *            are recognized in addition to unaligned file formats.
 */
int
esl_sqio_FormatCode(char *fmtstring)
{
  if (strcasecmp(fmtstring, "fasta")     == 0) return eslSQFILE_FASTA;
  if (strcasecmp(fmtstring, "embl")      == 0) return eslSQFILE_EMBL;
  if (strcasecmp(fmtstring, "genbank")   == 0) return eslSQFILE_GENBANK;
  if (strcasecmp(fmtstring, "ddbj")      == 0) return eslSQFILE_DDBJ;
  if (strcasecmp(fmtstring, "uniprot")   == 0) return eslSQFILE_UNIPROT;
#ifdef eslAUGMENT_MSA
  if (strcasecmp(fmtstring, "stockholm") == 0) return eslMSAFILE_STOCKHOLM;
  if (strcasecmp(fmtstring, "pfam")      == 0) return eslMSAFILE_PFAM;
#endif
  return eslSQFILE_UNKNOWN;
}


/* Function:  esl_sqio_DescribeFormat()
 * Synopsis:  Returns descriptive string for file format code.
 * Incept:    SRE, Sun Feb 27 09:24:04 2005 [St. Louis]
 *
 * Purpose:   Given a format code <fmt>, returns a string label for
 *            that format. For example, if <fmt> is <eslSQFILE_FASTA>,
 *            returns "FASTA". 
 *            
 *            When augmented by msa, then alignment file format codes
 *            are recognized in addition to unaligned file format codes.
 */
char *
esl_sqio_DescribeFormat(int fmt)
{
  switch (fmt) {
  case eslSQFILE_UNKNOWN:    return "unknown";
  case eslSQFILE_FASTA:      return "FASTA";
  case eslSQFILE_EMBL:       return "EMBL";
  case eslSQFILE_GENBANK:    return "Genbank";
  case eslSQFILE_DDBJ:       return "DDBJ";
  case eslSQFILE_UNIPROT:    return "Uniprot";
#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM: return "Stockholm";
  case eslMSAFILE_PFAM:      return "Pfam";
#endif
  default: esl_fatal("no such format code");
  }
  /*NOTREACHED*/
  return NULL;
}

/* Function:  esl_sqio_IsAlignment()
 * Incept:    SRE, Sun Feb 27 09:36:23 2005 [St. Louis]
 *
 * Purpose:   Returns TRUE if <fmt> is an alignment file
 *            format code; else returns FALSE.
 *            
 *            This function only checks the convention
 *            that <fmt> codes $<$100 are unaligned formats,
 *            and $\geq$100 are aligned formats. It does
 *            not check that <fmt> is a recognized format
 *            code.
 *
 */
int
esl_sqio_IsAlignment(int fmt)
{
  if (fmt >= 100) return TRUE;
  else            return FALSE;
}


/* Function:  esl_sqio_Position()
 * Incept:    SRE, Tue Mar 28 13:21:47 2006 [St. Louis]
 *
 * Purpose:   Reposition an open <sqfp> to offset <r_off>, which
 *            must be the offset to the start of a sequence record.
 *            
 *            Only normal sequence files can be positioned; not
 *            a standard input stream, gunzip stream, or a multiple
 *            alignment file interface.
 *            
 *            After <esl_sqio_Position()> is called, 
 *            <sqfp->linenumber> is relative to that start position.
 *            
 *            See the SSI module for manipulating offsets and indices.
 *
 * Returns:   <eslOK>     on success;
 *            <eslEINVAL> if the <sqfp> is not positionable;
 *            <eslEOF>    if no data can be read from this position.
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails;
 *            <eslEMEM> on (re-)allocation failure.
 */
int
esl_sqio_Position(ESL_SQFILE *sqfp, off_t r_off)
{
  int status;

  if (sqfp->do_stdin || sqfp->do_gzip  ||  
      esl_sqio_IsAlignment(sqfp->format)) return eslEINVAL;

  if (fseeko(sqfp->fp, r_off, SEEK_SET) != 0)
    ESL_EXCEPTION(eslESYS, "fseeko() failed");

  if (sqfp->is_linebased)
    {
      sqfp->linenumber = 0;
      if ((status = loadline(sqfp)) != eslOK) return status;
    }
  else
    {
      sqfp->linenumber = 1;
      sqfp->boff = r_off;
      sqfp->nc   = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      if (ferror(sqfp->fp)) { return eslESYS; }
    }
  return eslOK;
}

/* Function:  esl_sqio_Rewind()
 * Incept:    SRE, Tue Mar 28 14:10:56 2006 [St. Louis]
 *
 * Purpose:   Rewind an open <sqfp> to its beginning.   
 *
 *            Only normal sequence files can be positioned; not
 *            a standard input stream, gunzip stream, or a multiple
 *            alignment file interface.
 *
 * Returns:   <eslOK>     on success;
 *            <eslEINVAL> if the <sqfp> is not positionable;
 *            <eslEOF>    if no data can be read from this position.
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails;
 *            <eslEMEM> on (re-)allocation failure.
 */
int
esl_sqio_Rewind(ESL_SQFILE *sqfp)
{
  return esl_sqio_Position(sqfp, 0);
}



/*--------------------- end of i/o API ----------------------------*/





/*****************************************************************
 * Section 5. Internal routines for line-oriented parsers
 *****************************************************************/ 

static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}


/* Fetch a line into the sqfile's buffer.
 * 
 * On return:
 *    sqfp->buf        the next line in the file
 *    sqfp->nc         the length of buf in bytes, inclusive of \n
 *    sqfp->boff       disk offset to start of buf
 *    sqfp->pos        initialized to 0 (start of buf)
 *    sqfp->linenumber what line buf is (1..N in file)
 *    sqfp->balloc     current buffer allocation might have been increased
 * 
 * Returns <eslOK>  on success;
 *         <eslEOF> if no more data is left in file (no errmsg)
 * Throws  <eslEMEM> on realloc failure            
 */
static int
loadline(ESL_SQFILE *sqfp)
{
  int status;

  sqfp->boff = ftello(sqfp->fp);
  status = esl_fgets(&(sqfp->buf), &(sqfp->balloc), sqfp->fp);
  if (status != eslOK)  return status;

  sqfp->nc = strlen(sqfp->buf);
  sqfp->pos = 0;
  sqfp->linenumber++;
  return status;   
}

/* addseq():
 * 
 * <sqfp->buf> is a sequence data line.
 * Add residues from it to the growing sequence in <sq>.
 *
 * Uses the <sqfp->inmap> to decide whether to skip a character
 * (eslDSQ_IGNORED), report an error (eslDSQ_ILLEGAL), or store it.
 * 
 * On return:
 *   sq->seq     now includes residues from sqfp->buf; not nul-terminated yet;
 *                 generic_readseq() is responsible for the eventual \0.
 *   sq->n       has increased by # of residues on this line
 *   sq->salloc  may have increased, if sq->seq was reallocated
 *   
 *   sqfp->pos     points to \0 at end of buf.
 *   sqfp->rpl     might have been init'd to prev lastrpl, or invalidated (0)
 *   sqfp->bpl     might have been init'd to prev lastbpl, or invalidated (0)
 *   sqfp->lastrpl contains # of residues read from this line
 *   sqfp->lastbpl contains # of bytes on this line (incl of \n).
 *
 * Returns <eslOK> on success;
 *         <eslEFORMAT> on detecting an illegal character. (msg recorded)
 * Throws  <eslEMEM> on realloc failure.                   
 */
static int
addseq(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int   status;
  void *tmp;
  int   n0;
  int   symbol;

  /* First, some bookkeeping, based on the *previous* seq line we read.
   * Each time we add a sequence line, we're potentially responsible for 
   * updating rpl and bpl, for sequence indexing. (xref stl10/128)
   */
  if (sqfp->rpl != 0 && sqfp->lastrpl != -1) {
    if      (sqfp->rpl     == -1) 	 sqfp->rpl = sqfp->lastrpl; /* init */
    else if (sqfp->lastrpl != sqfp->rpl) sqfp->rpl = 0;	            /* inval*/
  }
  if (sqfp->bpl != 0 && sqfp->lastbpl != -1) {
    if      (sqfp->bpl     == -1)        sqfp->bpl = sqfp->lastbpl; /* init  */
    else if (sqfp->lastbpl != sqfp->bpl) sqfp->bpl = 0;             /* inval */
  }

  /* Now, add the line.
   */
  n0 = sq->n;
  while ((symbol = sqfp->buf[sqfp->pos]) != '\0')
    {
      if (esl_inmap_IsValid(sqfp->inmap, symbol) == eslOK)
	{
#ifdef eslAUGMENT_ALPHABET
	  if (sq->flags & eslSQ_DIGITAL)
	    sq->dsq[++sq->n] = esl_abc_DigitizeSymbol(sq->abc, symbol);
	  else
#endif
          sq->seq[sq->n++] = sqfp->inmap[symbol];
	  sqfp->pos++;
	}
      else if (! isascii(symbol))
	{
	  sprintf(sqfp->errbuf, "Non-ASCII char %d in sequence", symbol);
	  return eslEFORMAT;
	}
      else if (sqfp->inmap[symbol] == eslDSQ_ILLEGAL)
	{
	  sprintf(sqfp->errbuf, "Illegal %c in sequence", symbol);
	  return eslEFORMAT;
	}
      else if (sqfp->inmap[symbol] == eslDSQ_IGNORED)
	{
	  sqfp->pos++;
	}
      else 
	{			/* inmap[] shouldn't have any other value */
	  sprintf(sqfp->errbuf, "Internal inmap corruption");
	  return eslECORRUPT;
	}

      /* Realloc seq as needed. Careful, dsq runs 1..n, seq 0..n-1 */
#ifdef eslAUGMENT_ALPHABET
	  if (sq->flags & eslSQ_DIGITAL && sq->n == (sq->salloc-1))
	    {
	      ESL_RALLOC(sq->dsq, tmp, sizeof(ESL_DSQ) * sq->salloc * 2);
	      sq->salloc *= 2; /* doubling */
	    }
	  else
#endif
      if (sq->n == sq->salloc)
	{
	  ESL_RALLOC(sq->seq, tmp, sizeof(char) * sq->salloc * 2);
	  sq->salloc *= 2; /* doubling */
	}
    }

  sqfp->lastrpl = sq->n - n0;	/* remember # of residues on this line. */
  sqfp->lastbpl = sqfp->nc;     /* remember # of bytes on this line.    */
  return eslOK;			/* eslOK; eslEFORMAT, eslEMEM           */

 ERROR:
  return status;
}

/* generic_readseq():
 *
 * The <sqfp> is positioned at the beginning of the sequence data
 * in a record. If <sqfp->addfirst> is TRUE, the sequence data
 * includes the *current* line in <sqfp->buf>. Else, the first line
 * of sequence data is the *next* line in the file.
 * 
 * Reads sequence data until the format's endTest() returns TRUE
 * on the last line; or on EOF, if the format can also end a record on 
 * an EOF.
 *
 * On return:
 *   sq->seq    contains the sequence of this record, NUL-terminated
 *   sq->n      the number of residues in seq
 *   sq->doff   disk offset to the start of the sequence data
 *   sq->salloc might be increased, if seq was reallocated
 *
 *   sqfp->buf  contains the last line that was read. If <sqfp->addend>
 *              is TRUE, this line was read as sequence data for this
 *              record; if <sqfp->addend> is false, it is to be the first
 *              line of the next sequence record in the file.
 *   sqfp->nc   number of bytes in buf, incl of \n
 *   sqfp->pos  is undefined (it depends on whether buf was a seq line or not)
 *   sqfp->rpl  is either initialized (>0), inval (0), or left unset (-1)
 *   sqfp->bpl  (ditto)             
 *   sqfp->lastrpl  contains # of residues in the final seq data line.
 *   sqfp->lastbpl  # of bytes on the final seq data line.
 * 
 * Returns <eslOK> on success;
 *         <eslEOF> on abnormal (premature) EOF;     (no mesg)
 *         <eslEFORMAT> on any other format problem. (mesg recorded)
 * Throws <eslEMEM> on allocation failure.           
 */
static int
generic_readseq(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int status;
  int done   = 0;

  sqfp->lastrpl = sqfp->lastbpl = -1; /* init */

  if (sqfp->addfirst) {
    sq->doff = sqfp->boff;
    status = addseq(sqfp, sq);
    if (status != eslOK) return status;
  }
  
  do {
    status = loadline(sqfp);
    if      (status == eslEOF && sqfp->eof_is_ok) done = TRUE;
    else if (status != eslOK) return status;
    
    done |= sqfp->endTest(sqfp->buf);

    if (!done || sqfp->addend) {
      status = addseq(sqfp, sq);
      if (status != eslOK) return status;
    }
  } while (!done);

#ifdef eslAUGMENT_ALPHABET
  if (sq->flags & eslSQ_DIGITAL)
    sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else
#endif
  sq->seq[sq->n] = '\0';

  
  return status;	/* ESL_OK; ESL_EOF; ESL_EMEM */
}


/* set_name(), set_accession(), append_description();
 * 
 * Given a ptr <s> into a line buffer;
 * strtok it using <delim> to extract a sequence name/acc/desc;
 * copy (or append) that into <sq>, reallocating as needed.
 * 
 * sq->name is set; it may have been reallocated, in which
 * case sq->nalloc is increased; the buffer that <s> pointed
 * to is modified w/ a \0 by the strtok(). (Analogously for
 * acc, desc).
 *
 * Returns eslOK on success.
 * Returns eslEFORMAT if the strtok fails; (no mesg)
 * Throws  eslEMEM if a realloc fails.   
 */
static int
set_name(ESL_SQ *sq, char *s, char *delim) 
{
  void *tmp;
  char *tok;
  int   toklen;
  int   status;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  if (toklen >= sq->nalloc) {
    ESL_RALLOC(sq->name, tmp, sizeof(char) * (toklen+eslSQ_NAMECHUNK));
    sq->nalloc = toklen + eslSQ_NAMECHUNK;
  }
  strcpy(sq->name, tok);
  return eslOK;

 ERROR:
  return status;
}
static int
set_accession(ESL_SQ *sq, char *s, char *delim) 
{
  void *tmp;
  char *tok;
  int   toklen;
  int   status;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  if (toklen >= sq->aalloc) {
    ESL_RALLOC(sq->acc, tmp, sizeof(char) * (toklen+eslSQ_ACCCHUNK));
    sq->aalloc = toklen + eslSQ_ACCCHUNK;
  }
  strcpy(sq->acc, tok);
  return eslOK;

 ERROR:
  return status;
}
static int
append_description(ESL_SQ *sq, char *s, char *delim) 
{
  void *tmp;
  char *tok;
  int   toklen;
  int   status;
  int   dlen;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  dlen = strlen(sq->desc);

  if (dlen + toklen + 1 >= sq->dalloc) { /* +1 for \n */
    ESL_RALLOC(sq->desc, tmp, sizeof(char) * (toklen+dlen+eslSQ_DESCCHUNK));
    sq->dalloc = dlen + toklen + eslSQ_DESCCHUNK;
  }

  if (dlen > 0) sq->desc[dlen] = '\n';
  strcpy(sq->desc + dlen + 1, tok);
  return eslOK;

 ERROR:
  return status;
}
/*------------------- line-oriented parsers -----------------------*/


/*****************************************************************
 * EMBL format (including Uniprot and TrEMBL
 *****************************************************************/ 
/* EMBL and Uniprot protein sequence database format.
 * See: http://us.expasy.org/sprot/userman.html
 */

static void
config_embl(ESL_SQFILE *sqfp)
{
  int x;

  sqfp->is_linebased = TRUE;
  sqfp->addfirst     = FALSE;
  sqfp->addend       = FALSE;
  sqfp->eof_is_ok    = FALSE;	/* records end with // */
  sqfp->endTest      = &end_embl;

  /* The input map.
   */
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
}

/* read_embl()
 * 
 * Called by esl_sqio_Read() as the EMBL-specific parser;
 * <sqfp> is an opened <ESL_SQFILE>;
 * <sq> is an allocated and initialized <ESL_SQ>.
 * 
 * Returns <eslOK> on success and <sq> contains the input sequence.
 * 
 * Returns <eslEOF> on normal end: no sequence was read and we're
 * out of data in the file.  
 * 
 * Returns <eslEFORMAT> on a format problem, including illegal char in
 * the sequence; line number that the parse error occurs on is in
 * <sqfp->linenumber>, and an informative error message is placed in
 * <sqfp->errbuf>.
 * 
 * Throws <eslEMEM> on an allocation failure;
 *        <eslEINCONCEIVABLE> on internal error.
 */
static int
read_embl(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int   status;

  /* Find first line:
   * "Each entry must begin with an identification line (ID)..."
   * "The two-character line-type code that begins each line is always
   *  followed by three blanks..."
   */
  if (feof(sqfp->fp))  return eslEOF;
  while (is_blankline(sqfp->buf)) {
    if ((status = loadline(sqfp)) == eslEOF) return eslEOF; /* normal */
    else if (status != eslOK) return status; /* abnormal */
  } 

  /* ID line is defined as:
   *     ID   ENTRY_NAME DATA_CLASS; MOLECULE_TYPE; SEQUENCE_LENGTH.
   *  and we're only after the ENTRY_NAME.
   */
  if (strncmp(sqfp->buf, "ID   ", 5) != 0) {
    sprintf(sqfp->errbuf, "Failed to find ID line");
    return eslEFORMAT;
  }
  status = set_name(sq, sqfp->buf+5, " ");
  if (status != eslOK) {
    sprintf(sqfp->errbuf, "Failed to parse name on ID line"); 
    return status;
  }
  sq->roff = sqfp->boff;	/* record the offset of the ID line */
  
  /* Look for SQ line; parsing optional info as we go.
   */
  do {
    if ((status = loadline(sqfp)) != eslOK) {
      sprintf(sqfp->errbuf, "Failed to find SQ line");
      return eslEFORMAT;
    }

    /* "The format of the AC line is:
     *    AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]
     *  Researchers who wish to cite entries in their publications
     *  should always cite the first accession number. This is
     *  commonly referred to as the 'primary accession
     *  number'."
     */
    if (strncmp(sqfp->buf, "AC   ", 5) == 0)
      {
	status = set_accession(sq, sqfp->buf+5, ";");
	if (status != eslOK) {
	  sprintf(sqfp->errbuf, "Failed to parse accession on AC line");
	  return status;
	}
      }

    /* "The format of the DE line is:
     *    DE   Description.
     * ...In cases where more than one DE line is required, the text is
     * only divided between words and only the last DE line is
     * terminated by a period."
     */
    if (strncmp(sqfp->buf, "DE   ", 5) == 0)
      {
	status = append_description(sq, sqfp->buf+5, "\n");
	if (status != eslOK) {
	  sprintf(sqfp->errbuf, "Failed to parse description on DE line");
	  return status;
	}
      }

    /* "The format of the SQ line is:
     *  SQ   SEQUENCE XXXX AA; XXXXX MW; XXXXXXXXXXXXXXXX CRC64;"
     */
  } while (strncmp(sqfp->buf, "SQ   ", 5) != 0);
  
  /* Read the sequence
   */
  status = generic_readseq(sqfp, sq);
  if (status == eslEOF) { /* premature EOF becomes an EFORMAT error */
    sprintf(sqfp->errbuf, "Premature EOF; no // found at end of seq record");
    return eslEFORMAT;
  }
  else if (status != eslOK) return status; /* throw all other errors */
  
  /* Load next line */
  status = loadline(sqfp);
  if (status == eslEOF) return eslOK;	/* defer EOF report 'til next read */
  else if (status != eslOK) return status;

  return eslOK;
}

static int
end_embl(char *buf)
{
  if (strncmp(buf, "//", 2) == 0) return 1;
  return 0;
}
/*---------------------- EMBL format ---------------------------------*/



/*****************************************************************
 * Genbank format 
 *****************************************************************/ 
/* NCBI Genbank sequence database format.
 * See Genbank release notes; for example,
 * ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
 */

static void
config_genbank(ESL_SQFILE *sqfp)
{
  int x;

  sqfp->is_linebased = TRUE;
  sqfp->addfirst     = FALSE;
  sqfp->addend       = FALSE;
  sqfp->eof_is_ok    = FALSE;	/* records end with //  */
  sqfp->endTest      = &end_genbank;

  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  for (x = '0'; x <= '9'; x++) sqfp->inmap[x] = eslDSQ_IGNORED;
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
} 

static int
read_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int   status;

  /* Find LOCUS line, allowing for ignoration of a file header.
   */
  if (feof(sqfp->fp))  return eslEOF;
  while (strncmp(sqfp->buf, "LOCUS   ", 8) != 0) {
    if ((status = loadline(sqfp)) == eslEOF) return eslEOF; /* normal */
    else if (status != eslOK) return status; /* abnormal */
  } 
  status = set_name(sq, sqfp->buf+12, " ");
  if (status != eslOK) {
    sprintf(sqfp->errbuf, "Failed to parse name on LOCUS line"); 
    return status;
  }
  sq->roff = sqfp->boff;	/* record the disk offset to the LOCUS line */
  
  /* Look for ORIGIN line, parsing optional info as we go.
   */
  do {
    if ((status = loadline(sqfp)) != eslOK) {
      sprintf(sqfp->errbuf, "Failed to find ORIGIN line");
      return eslEFORMAT;
    }

    /* Optional VERSION line is parsed as "accession".
     */
    if (strncmp(sqfp->buf, "VERSION   ", 10) == 0)
      {
	status = set_accession(sq, sqfp->buf+12, " ");
	if (status != eslOK) {
	  sprintf(sqfp->errbuf, "Failed to parse 'accession' on VERSION line");
	  return status;
	}
      }

    /* Optional DEFINITION Line is parsed as "description".
     */
    if (strncmp(sqfp->buf, "DEFINITION ", 11) == 0)
      {
	status = append_description(sq, sqfp->buf+12, "\n");
	if (status != eslOK) {
	  sprintf(sqfp->errbuf, "Failed to parse desc on DEFINITION line");
	  return status;
	}
      }
  } while (strncmp(sqfp->buf, "ORIGIN", 6) != 0);
  
  /* Read the sequence
   */
  status = generic_readseq(sqfp, sq);
  if (status == eslEOF) { /* premature EOF becomes an EFORMAT error */
    sprintf(sqfp->errbuf, "Premature EOF; no // found at end of seq record");
    return eslEFORMAT;
  }
  else if (status != eslOK) return status; /* throw all other errors */
  
  /* Load next line */
  status = loadline(sqfp);
  if (status == eslEOF) return eslOK;	/* defer EOF report 'til next read */
  else if (status != eslOK) return status;

  return eslOK;
}

static int
end_genbank(char *buf)
{
  if (strncmp(buf, "//", 2) == 0) return 1;
  return 0;
}
/*----------------- end Genbank format -------------------------------*/




/*****************************************************************
 * FASTA format i/o
 *****************************************************************/

static void
config_fasta(ESL_SQFILE *sqfp)
{
  int x;

  sqfp->is_linebased = FALSE;
  sqfp->addfirst     = FALSE;	/* no-op in a fread() parser */
  sqfp->addend       = FALSE;	/* ditto */
  sqfp->eof_is_ok    = TRUE;	/* unused, but fasta can indeed end w/ eof. */
  sqfp->endTest      = NULL;	/* unused in a fread() parser */
  
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
  /* \n is special - fasta reader detects it as an eol */
}

/* read_fasta()
 * SRE, Thu Dec 23 13:57:59 2004 [Zaragoza]
 *
 * Purpose:   Given an open <sqfp> for a FASTA file; read the next 
 *            sequence into <s>. Caller is responsible for creating
 *            <s> initially; but it will be reallocated here if its space is 
 *            insufficient.
 *            
 *            <sqfp->pos> is at the first byte in the file (which
 *            must be a $>$ if it's FASTA format); or at a '$>$' 
 *            for a subsequent sequence 2..N; or at EOF, byte B+1
 *            in a B-byte file, in which case we return <eslEOF>.
 *            One of these conditions is guaranteed if you only call 
 *            <read_fasta()> on an open FASTA file, but
 *            operations that muck with the internals of a <sqfp> 
 *            (such as indexing/lookup) have to be careful of this
 *            requirement.
 *            
 *            The file must be a UNIX or DOS/Windows textfile, obeying
 *            EOL conventions of \verb+\n+ or \verb+\r\n+. Mac files pre 
 *            OS9 that use an EOL convention of \verb+\r+ will fail.
 *
 * Args:      sqfp   - open ESL_SQFILE for reading a FASTA-format datafile
 *            s      - allocated ESL_SQ object         
 *
 * Returns:   <eslOK> on success; the newly read sequence info
 *               is stored in <s>.
 *            <eslEOF> when there is no sequence left in the file;
 *               (including first attempt to read an empty file).
 *            <eslEFORMAT> if there's a problem with the format,
 *               such as an illegal character. The linenumber that 
 *               the error occurred at is in <sqfp->linenumber>, 
 *               and the line itself is <sqfp->buf>, which an 
 *               application can use to format a useful error
 *               message.
 *
 * Throws:    <eslEMEM> on an allocation failure, either in 
 *            name/description lengths or in the sequence data
 *            itself.
 *            <eslECORRUPT> if the inmap is corrupted somehow.
 *
 * Xref:      STL8/p148; 2004/1224-fileread-speed
 *            Design goals: improve speed over SQUID's ReadSeq(); remove
 *            name, description length limitations by dynamically allocating
 *            as needed; impose no line length limitations.
 *            
 *            Redesigned to use fread() and character-based parsing (with
 *            a finite automaton) instead of fgets() equivalents and line-based
 *            parsing (with strtok() equivalents); and to use an inmap[]
 *            to validate sequence characters instead of checking against
 *            a string of valid chars. Approximate 4x speedup relative to
 *            SQUID ReadSeq(). 
 *            
 *            The trickiest bit of the implementation is the
 *            interaction between dynamic allocation of name, desc,
 *            seq while processing one char at a time in an FSA; at
 *            each char, have to make sure you *both* have a place to
 *            put it, and that it's loaded in the input buffer, but
 *            you don't want to be making two checks per char. Thus the
 *            "nsafe" idiom in this code, which calls check_buffers()
 *            to check both the input buffer and the storage area, and
 *            returns the number of chars we can safely read without
 *            having to check either one again. This results in an idiom
 *            of two nested while() loops, which necessitates using
 *            goto's to break out of that nesting.
 *            
 *            In the end, this code is about 4x faster than squid's ReadSeq(),
 *            and is actually almost 2x faster than test code that simply
 *            fread()'s the file and counts '>' characters, which is puzzling,
 *            but indicates we're running pretty efficiently. The test
 *            data are in 1224-fileread-speed. 
 */
static int
read_fasta(ESL_SQFILE *sqfp, ESL_SQ *s)
{
  int   npos = 0;	/* position in stored name               */
  int   dpos = 0;	/* position in stored description        */
  int   spos = 0;	/* position in stored sequence           */
  int   c;		/* character we're processing            */
  char *buf;		/* ptr to the input buffer               */
  int   nc;		/* number of chars in the buffer         */
  int   pos;		/* position in the buffer                */
  unsigned char *inmap; /* ptr to the input map                  */
  char *seq;            /* ptr to the growing sequence           */
  char **seq_addr;      /* address of seq* or NULL if seq = NULL */
  ESL_DSQ *dsq;         /* ptr to growing digitized seq          */
  ESL_DSQ **dsq_addr;   /* address of dsq* or NULL if dsq = NULL */
  int   state;		/* state of our FSA parser               */
  int   nsafe;          /* #bytes we can safely move in both input/storage */
  int   at_linestart;	/* TRUE when we're at first char of a data line    */

  /* If we have no more data, return EOF; we're done. (a normal return)
   */
  if (sqfp->nc == 0 && feof(sqfp->fp)) return eslEOF;

  buf   = sqfp->buf;	/* These are some optimizations, avoid dereferences */
  nc    = sqfp->nc;
  pos   = sqfp->pos;
  inmap = sqfp->inmap;
  dsq   = NULL;
#ifdef eslAUGMENT_ALPHABET
  if (s->flags & eslSQ_DIGITAL)
    { dsq  = s->dsq; seq = NULL; } 
  else
#endif
  seq   = s->seq;
  /* We parse one char at a time with simple state machine; states are indexed:
   *   0 = START     (on the >)    accepts >, moves to 1
   *   1 = NAMESPACE (>^name)      accepts space, stays;    else moves to 2
   *   2 = NAME                    accepts nonspace, stays; else moves to 3
   *   3 = DESCSPACE (name^desc)   accepts space and stays; else moves to 4
   *   4 = DESC                    accepts \n to move to 5; else stays.
   *   5 = SEQ                     accepts !> and stays;    else (on >) 6
   *   6 = END                           
   */
  state = 0;
  while (state != 6) {
    switch (state) {
    case 0:     /* START. Accept >, move to state 1. Skip blank lines.*/
      if      (isspace(buf[pos])) { pos++; }
      else if (buf[pos] == '>')   
	{ 
	  s->roff = sqfp->boff + pos;
	  pos++; 
	  state = 1;
	} 
      else 
	{ 
	  sprintf(sqfp->errbuf, "No > name/descline found."); 
	  return eslEFORMAT; 
	}
      break;

    case 1:     /* NAMESPACE. Switch to NAME state when we see a non-space. */
      c = buf[pos];
      if   (c != ' ' && c != '\t') state = 2; else pos++; 
      break;

    case 2:     /* NAME. Accept/store non-whitespace. Else, move to state 3. */
      while ((nsafe = check_buffers(sqfp->fp, &(sqfp->boff), buf, &nc, &pos, 
				    &(s->name), NULL, npos, &(s->nalloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if   (!isspace(c)) { s->name[npos++] = c;  pos++;     }
	  else { state = 3; goto NAMEDONE; }
	}
      }
      if (nsafe == -1) ESL_EXCEPTION(eslEMEM, "realloc failed");
      sprintf(sqfp->errbuf, "File ended within a seq name."); 
      return eslEFORMAT; /* we ran out of data while still in a name. */
    NAMEDONE:
      break;

    case 3:   /* DESCSPACE.  Accepts non-newline space and stays; else 4. */
      c = buf[pos]; 
      if (c != ' ' && c != '\t') state = 4; else pos++;
      break;
	  
    case 4:   /* DESC. Accepts and stores up to \n; accepts \n & moves to 5. */
      while ((nsafe = check_buffers(sqfp->fp, &(sqfp->boff), buf, &nc, &pos,
				    &(s->desc), NULL, dpos, &(s->dalloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if      (c != '\n' && c != '\r') 
	    { s->desc[dpos++] = c;  pos++; }
	  else if (c == '\r') 
	    { pos++; } /* ignore \r part of DOS \r\n EOL */
	  else                
	    { state = 5; pos++; sqfp->linenumber++; goto DESCDONE; }
	}
      }
      if (nsafe == -1) ESL_EXCEPTION(eslEMEM, "realloc failed");
      sprintf(sqfp->errbuf, "File ended within a description line."); 
      return eslEFORMAT;	/* ran out of data while still in desc */
    DESCDONE:
      break;

    case 5:   /* SEQ. Accept/process according to inmap; on '>',  move to 6. */
      s->doff = sqfp->boff + pos;
      sqfp->lastrpl = sqfp->lastbpl = -1;
      at_linestart  = TRUE;

      if(seq == NULL) seq_addr = NULL; else seq_addr = &seq;
      if(dsq == NULL) dsq_addr = NULL; else dsq_addr = &dsq;

      while ((nsafe = check_buffers(sqfp->fp, &(sqfp->boff), buf, &nc, &pos, 
				    seq_addr, dsq_addr, spos, &(s->salloc))) > 0) {
	while (nsafe--) {
	  /* At start of every new data line, do bookkeeping for rpl, bpl
	   * based on *previous* line lengths */
	  if (at_linestart) {
	    if (sqfp->rpl != 0 && sqfp->lastrpl != -1) {
	      if      (sqfp->rpl     == -1) 	 sqfp->rpl = sqfp->lastrpl; /* init */
	      else if (sqfp->lastrpl != sqfp->rpl) sqfp->rpl = 0;	            /* inval*/
	    }
	    if (sqfp->bpl != 0 && sqfp->lastbpl != -1) {
	      if      (sqfp->bpl     == -1)        sqfp->bpl = sqfp->lastbpl; /* init  */
	      else if (sqfp->lastbpl != sqfp->bpl) sqfp->bpl = 0;             /* inval */
	    }
	    at_linestart = FALSE;
	    sqfp->lastrpl = 0;
	    sqfp->lastbpl = 0;
	  }

	  /* bookkeeping complete, now deal with the character.
	   */
	  c = buf[pos];
	  if (esl_inmap_IsValid(sqfp->inmap, c))
	    {
#ifdef eslAUGMENT_ALPHABET
	      /* 02.27.07 */
	      if (s->flags & eslSQ_DIGITAL)
		dsq[++spos] = esl_abc_DigitizeSymbol(s->abc, c);
	      else
#endif
		seq[spos++] = c;
	      pos++;
	      sqfp->lastrpl++;
	    }
	  else if (! isascii(c))
	    {
	      sprintf(sqfp->errbuf, "Non-ASCII char %d found in sequence.", c); 
	      return eslEFORMAT;
	    }
	  else if (c == '>')               
	    goto FINISH; 
	  else if (c == '\n') 	/* end of a seq line. */
	    { 
	      pos++; 
	      sqfp->linenumber++; 
	      at_linestart = TRUE;
	    }
	  else if (inmap[c] == eslDSQ_ILLEGAL) 
	    {
	      sprintf(sqfp->errbuf, "Illegal char %c found in sequence.", c); 
	      return eslEFORMAT;
	    }
	  else if (inmap[c] == eslDSQ_IGNORED) 
	    { pos++; } 
	  else		
	    {
	      sprintf(sqfp->errbuf, "Internal corruption of an inmap"); 
	      return eslECORRUPT;
	    }

	  sqfp->lastbpl++;
	}
      }
      if (nsafe == -1) ESL_EXCEPTION(eslEMEM, "realloc failed");
      state = 6;
      break;
    } /* end of switch over FSA states */

    if (pos == nc) {		/* reload the buffer when it empties */
      sqfp->boff = ftello(sqfp->fp);
      nc  = fread(buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      pos = 0;
    }
    
  } /* end of while loop waiting to reach state 6 */

 FINISH:
  /* check_buffers() was careful to leave at least one
   * free byte (or two in case of dsq) on the storage 
   * strings, for the NUL (or sentinel), so we don't 
   * need to check/reallocate them again.
   */
  s->name[npos] = '\0';
  s->desc[dpos] = '\0';
#ifdef eslAUGMENT_ALPHABET
  if (s->flags & eslSQ_DIGITAL)
    dsq[spos+1] = eslDSQ_SENTINEL;
  else
#endif
  seq[spos]  = '\0';

  /* Reset the data that we took copies of.
   */
#ifdef eslAUGMENT_ALPHABET
  if (s->flags & eslSQ_DIGITAL)
    s->dsq = dsq;
  else
#endif
    s->seq        = seq;
  s->n          = spos;
  sqfp->pos     = pos;
  sqfp->nc      = nc;

  return eslOK;
}

/* write_fasta():
 * SRE, Fri Feb 25 16:18:45 2005
 *
 * Write a sequence <s> in FASTA format to the open stream <fp>.
 *
 * Returns <eslOK> on success.
 */
static int
write_fasta(FILE *fp, ESL_SQ *s)
{
  char  buf[61];
  int   pos;

  fprintf(fp, ">%s %s %s\n", s->name, s->acc, s->desc);
  buf[60] = '\0';
  for (pos = 0; pos < s->n; pos += 60)
    {
      strncpy(buf, s->seq+pos, 60);
      fprintf(fp, "%s\n", buf);
    }
  return eslOK;
}

#ifdef eslAUGMENT_ALPHABET
/* write_digital_fasta():
 * SRE, Tue Jan  9 16:26:52 2007 [Janelia] [Keb' Mo', Suitcase]
 * 
 * Write a digital sequence <s> in FASTA format to open stream <fp>.
 * 
 * Returns <eslOK> on success.
 */
static int
write_digital_fasta(FILE *fp, ESL_SQ *s)
{
  char buf[61];
  int  pos;
  
  fprintf(fp, ">%s", s->name);
  if (s->acc[0]  != 0) fprintf(fp, " %s", s->acc);
  if (s->desc[0] != 0) fprintf(fp, " %s", s->desc);
  fputc('\n', fp);

  buf[60] = '\0';
  for (pos = 1; pos <= s->n; pos+=60)
    {
      esl_abc_TextizeN(s->abc, s->dsq+pos, 60, buf);
      fputs(buf, fp);
      fputc('\n', fp);
    }
  return eslOK;
}
#endif /*eslAUGMENT_ALPHABET*/


/* check_buffers()
 * 
 * Given the input fread() buffer, and the storage location for where
 * we're putting the name, acc, desc, seq. If we've run out of input buffer, 
 * fread() a new block. If we're out of storage space, reallocate
 * it by doubling. Return the minimum number of bytes we can safely
 * process before we either run out of input buffer or we run out of storage
 * space, or -1 if a realloc fails.
 *
 * EPN: Added digitized sequence flexibility in ugly way, to reallocate
 *      a ESL_DSQ * instead of char *, pass in NULL for 'char **s' and
 *      a valid ESL_DSQ * address for 'ESL_DSQ **dsq'. Added contract
 *      check to ensure exactly one of these is NULL. 
 *
 * This supports an idiom of
 *     while (nsafe = check_buffers()) {
 *       while (nsafe--) {
 *         process buf[pos];
 *         pos++;
 *       }        
 *     }
 *     if (nsafe == -1) ESL_EXCEPTION(eslEMEM, "realloc failed");
 * 
 * which avoids having to check our buffers every character.
 * 
 */
static int
check_buffers(FILE *fp, off_t *boff, char *buf, int *nc, int *pos, 
	      char **s, ESL_DSQ **dsq, int i, int *slen)
{
  int inlen, savelen;

  /* Contract check. */
  if (dsq == NULL && s == NULL) ESL_EXCEPTION(eslEINVAL, "both s and dsq NULL, exactly 1 should be NULL.\n");
  if (dsq != NULL && s != NULL) ESL_EXCEPTION(eslEINVAL, "both s and dsq non-NULL, exactly 1 should be NULL\n");

  inlen = *nc - *pos;  	/* can read this many bytes before reloading buffer */
  if (inlen == 0)	/* if we're at the end, reload now. */
    {
      *boff = ftello(fp);
      *nc   = fread(buf, sizeof(char), eslREADBUFSIZE, fp);
      *pos  = 0;
      inlen = *nc - *pos;	/* (if this is still 0, we're at EOF.) */
    }

  /* Note the -1 on savelen, which leaves room for a NUL terminator.
   */
  if(dsq == NULL) savelen = *slen - i - 1;	/* can save this many before realloc'ing */
  else            savelen = *slen - i - 2;	/* can save this many before realloc'ing */
  if (savelen == 0)		/* if we need to reallocate now...       */
    {				/* then double our space. */
      savelen = *slen;		
      *slen  += *slen;		
      /* Exactly 1 of s and dsq is NULL, part of the contract */
      if(s != NULL)
	{
	  *s = realloc(*s, sizeof(char) * *slen);
	  if (*s == NULL) return -1;
	}
      else /* dsq must be non-NULL */
	{
	  *dsq = realloc(*dsq, sizeof(ESL_DSQ) * *slen);
	  if (*dsq == NULL) return -1;
	}
    }

  /* Now, return the minimum safe bytecount.
   */
  if (savelen < inlen) return savelen;  else return inlen;
}
/*------------------- end of FASTA i/o ---------------------------*/	       




/*****************************************************************
 * Functions specific to sqio <-> msa interoperation; 
 * require augmentation w/ msa module.
 *****************************************************************/

#ifdef eslAUGMENT_MSA
/* Function:  esl_sq_FetchFromMSA()
 * Synopsis:  Fetch a single sequence from an MSA.
 * Incept:    SRE, Sun Mar 30 13:39:06 2008 [Janelia]
 *
 * Purpose:   Retrieve sequence number <which> from <msa>, in a newly
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
esl_sq_FetchFromMSA(ESL_MSA *msa, int which, ESL_SQ *ret_sq)
{
  ESL_SQ *sq       = NULL;
  char    gapchars = "-_.";	/* hardcoded for now */

  if (which >= msa->nseq || which < 0) return eslEOD;

  if (! (msa->flags & eslMSA_DIGITAL))
    {
      if ((sq = esl_sq_CreateFrom(msa->sqname[which], msa->aseq[which],
				  msa->sqdesc[which], msa->sqacc[which], msa->ss[which])) == NULL) goto ERROR;
      esl_strdealign(sq->ss,  sq->seq, gapchars, NULL);
      esl_strdealign(sq->seq, sq->seq, gapchars, &(sq->n));
    }
#ifdef eslAUGMENT_ALPHABET
  else
    {
      if ((sq = esl_sq_CreateDigitalFrom(msa->abc, msa->sqname[which], msa->ax[which], msa->alen, 
					 msa->sqdesc[which], msa->sqacc[which], msa->ss[which])) == NULL) goto ERROR;
      esl_sq_CDealign(sq->abc, sq->ss,  sq->dsq, NULL);
      esl_sq_XDealign(sq->abc, sq->dsq, sq->dsq, &(sq->n));
    }
#endif

  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  return eslEMEM;
}


/* extract_sq_from_msa():
 * Move sequence <idx> from the <msa> into <s>, and dealign
 * it and any associated per-residue annotation.
 * 
 * This is destructive - the pointers are redirected so that the <s>
 * structure now points to data that the <msa> previously maintained,
 * the <msa> data is NULL'ed, and any previous storage in <s> is free'd. 
 * The <esl_msa_Destroy()> function checks for NULL'ed fields before 
 * freeing them, so this wholesale pillaging of the <msa> is safe, so
 * long as the caller has no intention of using it for anything else.
 * 
 * Limitation: hardcodes the gapstring "-_."
 */
static int
extract_sq_from_msa(ESL_MSA *msa, int idx, ESL_SQ *s)
{
  int n;
  
  /* Name.
   */
  n = strlen(msa->sqname[idx]);
  if (s->name != NULL) free(s->name);
  s->name   = msa->sqname[idx];
  s->nalloc = n;
  msa->sqname[idx] = NULL;

  /* Accession.
   */
  if (msa->sqacc != NULL && msa->sqacc[idx] != NULL)
    {
      n = strlen(msa->sqacc[idx]);
      if (s->acc != NULL) free(s->acc);
      s->acc    = msa->sqacc[idx];
      s->aalloc = n;
      msa->sqacc[idx] = NULL;
    }
  
  /* Description.
   */
  if (msa->sqdesc != NULL && msa->sqdesc[idx] != NULL)
    {
      n = strlen(msa->sqdesc[idx]);
      if (s->desc != NULL) free(s->desc);
      s->desc   = msa->sqdesc[idx];
      s->dalloc = n;
      msa->sqdesc[idx] = NULL;
    }

  /* Sequence... still aligned, for now
   */
  if (s->seq != NULL) free(s->seq);
  s->seq         = msa->aseq[idx];
  s->salloc      = msa->alen;
  msa->aseq[idx] = NULL;
  
  /* Structure... still aligned, for now
   */
  if (msa->ss != NULL && msa->ss[idx] != NULL)
    {
      if (s->ss != NULL) free(s->ss);
      s->ss        = msa->ss[idx];
      msa->ss[idx] = NULL;
    }

  /* Digital seq (dsq) is UNIMPLEMENTED, untouched here;
   * and optmem is untouched.
   */

  /* Dealign the ss and the seq.
   * ASSUMES that the gap characters are -_.
   */
  esl_strdealign(s->ss,  s->seq, "-_.", NULL);
  esl_strdealign(s->seq, s->seq, "-_.", &(sq->n));

  return eslOK;
}

/* convert_sq_to_msa()
 * SRE, Fri Feb 25 16:06:18 2005
 * 
 * Given a <sq>, create and return an "MSA" through <ret_msa>, which
 * contains only the single unaligned sequence. <sq> is 
 * not affected in any way. This is only to convert from the SQ
 * object to an MSA object for the purpose of writing SQ in an MSA format
 * file format.
 * 
 * Returns <eslOK> on success, and <*ret_msa> points to
 * a new "alignment".
 * 
 * Throws <eslEMEM> on allocation error, and <*ret_msa> is NULL.
 */
static int
convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa)
{
  ESL_MSA *msa;
  int      status;

  if ((msa = esl_msa_Create(1, sq->n)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((status = esl_strdup(sq->name, -1, &(msa->sqname[0]))) != eslOK) goto ERROR;
  
  if (*sq->acc != '\0')
    {
      ESL_ALLOC(msa->sqacc, sizeof(char *) * 1);
      if ((status = esl_strdup(sq->acc, -1, &(msa->sqacc[0]))) != eslOK) goto ERROR;
    }

  if (*sq->desc != '\0')
    {
      ESL_ALLOC(msa->sqdesc, sizeof(char *) * 1);
      if ((status = esl_strdup(sq->desc, -1, &(msa->sqdesc[0]))) != eslOK) goto ERROR;
    }

  strcpy(msa->aseq[0], sq->seq);
  
  if (sq->ss != NULL)
    {
      ESL_ALLOC(msa->ss, sizeof(char *) * 1);
      if ((status = esl_strdup(sq->ss, -1, &(msa->ss[0]))) != eslOK) goto ERROR;
    }
  
  msa->alen = sq->n;
  msa->nseq = 1;

  *ret_msa = msa;
  return eslOK;

 ERROR:
  esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}

#endif /*eslAUGMENT_MSA*/
/*---------- end of msa <-> sqio module interop -----------------*/

/*****************************************************************
 * Section 6. Test and example code.
 *****************************************************************/ 





/*****************************************************************
 * Test driver:
 * gcc -g -Wall -I. -L. -o testdrive -DeslSQIO_TESTDRIVE esl_sqio.c -leasel -lm
 * ./testdrive
 *****************************************************************/
#ifdef eslSQIO_TESTDRIVE
/*::cexcerpt::sqio_test::begin::*/
#include <stdlib.h>
#include <stdio.h>
#include <easel.h>
#include <esl_sqio.h>

int
main(void)
{
  ESL_SQ     *sq;
  ESL_SQFILE *sqfp;
  FILE       *fp;
  char        seq1[] = "GAATTC";
  char        seq2[] = "AAGCTT";
  char        tmpfile[32] = "esltmpXXXXXX";
  int         n;
  int         status;
  char       *textseq = NULL; /* textized digitized seq */

#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET   *abc  = NULL;
#endif

  /* Create a FASTA file containing two sequences.
   */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("failed to open tmpfile");
  fprintf(fp, ">seq1 seq1's description goes here\n");
  fprintf(fp, "%s\n", seq1);
  fprintf(fp, ">seq2 seq2's description goes here\n");
  fprintf(fp, "%s\n", seq2);	  
  fclose(fp);

  /* Example of the API for opening and reading 
   * seqs from a FASTA file.
   */
  if (esl_sqfile_Open(tmpfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) abort();
  sq = esl_sq_Create();

  n=0;
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      if (n==0 && strcmp(sq->seq, seq1) != 0) abort();
      if (n==1 && strcmp(sq->seq, seq2) != 0) abort();

      n++;
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %d, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);

  /* Potentially repeat using digital sequences. */
#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL) 
    esl_fatal("alphabet creation failed");

  if (esl_sqfile_Open(tmpfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) abort();
  sq = esl_sq_CreateDigital(abc);

  n=0;
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      ESL_ALLOC(textseq, sizeof(char) * (sq->n+2));
      esl_abc_Textize(abc, sq->dsq, sq->n, textseq);
      if (n==0)
	if(strcmp(textseq, seq1) != 0) abort();
      if (n==1)
	if(strcmp(textseq, seq2) != 0) abort();
      n++;
      esl_sq_Reuse(sq);
      free(textseq);
    }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %d, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
#endif

  remove(tmpfile);
  return 0;

 ERROR:
  if (textseq     != NULL)  free(textseq);
  return status;
}
/*::cexcerpt::sqio_test::end::*/
#endif /*eslSQIO_TESTDRIVE*/

/*****************************************************************
 * Benchmark driver (testing sqio speed)
 * gcc -O2 -I. -L. -o benchmark -DeslSQIO_BENCHMARK esl_sqio.c -leasel
 * ./benchmark <seqfile>
 *****************************************************************/
#ifdef eslSQIO_BENCHMARK
#include <stdlib.h>
#include <stdio.h>
#include <easel.h>
#include <esl_sqio.h>
#include <esl_stopwatch.h>

int
main(int argc, char **argv)
{
  ESL_STOPWATCH *w;
  ESL_SQ        *sq;
  ESL_SQFILE    *sqfp;
  FILE          *fp;
  char          *filename;
  int            n;
  int            format = eslSQFILE_FASTA;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET   *abc  = NULL;
#endif

  filename = argv[1];

  w = esl_stopwatch_Create();
  sq = esl_sq_Create();
  if (esl_sqfile_Open(filename, format, NULL, &sqfp) != eslOK) abort();

  n=0;
  esl_stopwatch_Start(w);
  while (esl_sqio_Read(sqfp, sq) == eslOK)
    {
      n++;
      esl_sq_Reuse(sq);
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, NULL);
  printf("Read %d sequences.\n", n);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_stopwatch_Destroy(w);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL) 
    esl_fatal("alphabet creation failed");

  /* EPN digital mode: repeat all, only diff is esl_sq_CreateDigital() call */
  w = esl_stopwatch_Create(); /* This would be unnec if there's an equivalent 
			       * func to StopWatchZero() in Squid. */
  sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_Open(filename, format, NULL, &sqfp) != eslOK) abort();

  n=0;
  esl_stopwatch_Start(w);
  while (esl_sqio_Read(sqfp, sq) == eslOK)
    {
      n++;
      esl_sq_Reuse(sq);
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, NULL);
  printf("Read %d sequences in digital mode.\n", n);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
#endif

  return 0;
}
#endif /*eslSQIO_BENCHMARK*/

/*****************************************************************
 * Example main()
 *****************************************************************/
#ifdef eslSQIO_EXAMPLE
/*::cexcerpt::sqio_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslSQIO_EXAMPLE esl_sqio.c easel.c
 * run:     ./example <FASTA file>
 */
#include <easel.h>
#include <esl_sqio.h>

int
main(int argc, char **argv)
{
  ESL_SQ     *sq;
  ESL_SQFILE *sqfp;
  int         status;
  int         format = eslSQFILE_UNKNOWN;
  char       *seqfile = argv[1];
  int         type;

  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  sq = esl_sq_Create();
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
  {
    /* use the sequence for whatever you want */
    printf("Read %12s: length %6d\n", sq->name, sq->n);
    esl_sq_Reuse(sq);
  }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %d, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  return 0;
}
/*::cexcerpt::sqio_example::end::*/
#endif /*eslSQIO_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
