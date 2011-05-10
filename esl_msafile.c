/* Multiple sequence alignment file i/o 
 * 
 * Table of contents:
 *    1. Opening/closing an ESLX_MSAFILE.
 *    2. Guessing file formats.
 *    3. Guessing alphabet.
 *    4. Reading an MSA from an ESLX_MSAFILE.
 *    5. Writing an MSA to a stream.
 *    x. Copyright and license.
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_mem.h"
#include "esl_msafile.h"
#include "esl_alphabet.h"

/*****************************************************************
 *# 1. Opening/closing an ESLX_MSAFILE
 *****************************************************************/
static int msafile_SetInmap(ESLX_MSAFILE *afp); 

/* Function:  eslx_msafile_Open()
 * Synopsis:  Open a multiple sequence alignment file for input.
 *
 * Purpose:   Open a multiple sequence alignment file <msafile> for input.
 *
 *            Caller asserts that <msafile> is in format code
 *            <format>, such as <eslMSAFILE_STOCKHOLM>,
 *            <eslMSAFILE_AFA>, <eslMSAFILE_CLUSTAL>.  If <format> is
 *            <eslMSAFILE_UNKNOWN>, format autodetection is performed.
 *            
 *            Alignments may be input in either text or digital mode,
 *            depending on the setting of the passed-by-reference
 *            alphabet pointer <byp_abc>. If caller passes NULL for
 *            the <byp_abc> argument, input is in text mode. If caller
 *            provides a valid non-NULL <byp_abc> pointer but
 *            <*byp_abc> is NULL (that is, caller has declared
 *            <ESL_ALPHABET *abc = NULL> and passed <&abc> as an
 *            argument), then we attempt to guess the digital alphabet
 *            using <eslx_msafile_GuessAlphabet()>, based on the first
 *            alignment in the input. In this case, the new alphabet
 *            is allocated here and returned to the caller. If caller
 *            provides a digital alphabet (that is, <ESL_ALPHABET *abc
 *            = esl_alphabet_Create...()> and passed <&abc>), that's
 *            the alphabet we use.
 *            
 *            Optionally, caller can provide in <env> the name of an
 *            environment variable ("PFAMDB", perhaps), in which the
 *            routine can find a colon-delimited list of directories.
 *            Then, if <msafile> is not found in the current working
 *            directory, we look for it in these directories, in the
 *            order they're listed. Otherwise caller sets <env> to
 *            <NULL>.
 *
 *            <msafile> is usually the name of a file. Alignments may
 *            also be read from standard input, or from
 *            gzip-compressed files.  If <msafile> is ``-'', alignment
 *            input is taken from the standard input stream. If
 *            <msafile> ends in ``.gz'', alignment input is read
 *            through a pipe from <gzip -dc>.
 *
 * Args:      byp_abc   - digital alphabet to use, or NULL for text mode
 *                        if <*byp_abc> is NULL, guess the digital alphabet,
 *                        create it, and return it in <*byp_abc>.
 *                        If <*byp_abc> is a digital alphabet, use it.
 *            msafile   - name of alignment input to open;
 *                        if "-", read standard input;
 *                        if "*.gz", read through a <gzip -dc> pipe.
 *            format    - format code, such as <eslMSAFILE_STOCKHOLM>;
 *                        or <eslMSAFILE_UNKNOWN> to autodetect format.
 *            env       - <NULL>, or the name of an environment variable
 *                        containing colon-delimited list of directories
 *                        in which to search for <msafile> (e.g. "PFAMDB").
 *            *ret_afp  - RETURN: open MSA input stream.
 *
 * Returns:   <eslOK> on success, and <*ret_afp> is the newly opened msa file.
 *
 *            <eslENOTFOUND> if <msafile> doesn't exist or can't be
 *            opened for reading; or (in the case of a <.gz> file) if
 *            a <gzip> executable doesn't exist in user's <PATH> or
 *            can't be executed.
 *
 *            <eslFAIL> in the case of a <.gz> file and the <gzip -dc> command
 *            fails on it.
 *            
 *            <eslENODATA> if the <msafile> is empty (contains no
 *            alignment), but we tried to read one to guess the
 *            digital alphabet.
 *            
 *            <eslEFORMAT> if we tried to autodetect the file format
 *            (caller provided <format=eslMSAFILE_UNKNOWN>), and
 *            failed.
 *            
 *            <eslEAMBIGUOUS> if we tried to autodetect the alphabet
 *            (based on the first alignment in the file), but its
 *            alphabet could not be reliably guessed.
 *
 *            On any of these normal errors, <*ret_afp> is returned in
 *            an error state, containing a user-directed error message
 *            in <afp->errmsg> and (if relevant) the full path to
 *            <msafile> that we attempted to open in
 *            <afp->bf->filename>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> on a system call failure, such as <fread()>.
 *            <eslEINVAL> if we tried to use <stdin> but the <stdin> stream was
 *            invalid (in an error state, <NULL>, or at <EOF>).
 *            On thrown exceptions, <*ret_afp> is <NULL>.
 */
int
eslx_msafile_Open(ESL_ALPHABET **byp_abc, const char *msafile, int format, const char *env, ESLX_MSAFILE **ret_afp)
{
  ESLX_MSAFILE *afp       = NULL;
  ESL_ALPHABET *abc       = NULL;
  int           alphatype = eslUNKNOWN;
  int           status;
  int           x;

  ESL_ALLOC(afp, sizeof(ESLX_MSAFILE));
  afp->bf         = NULL;
  afp->linenumber = 0;
  afp->abc        = NULL;
  afp->ssi        = NULL;

  if ((status = esl_buffer_Open(msafile, env, &(afp->bf))) != eslOK)
    ESL_XFAIL(status, afp->errmsg, "%s", afp->bf->errmsg); /* ENOTFOUND; FAIL are normal here */

  /* Determine the format */
  if (format == eslMSAFILE_UNKNOWN && (status = eslx_msafile_GuessFileFormat(afp->bf, &format)) != eslOK)
    ESL_XFAIL(status, afp->errmsg, "couldn't determine alignment input format"); /* EFORMAT is normal failure */
  afp->format = format;

  /* Set up a text-mode inmap. (We may soon switch to digital mode, but if we're
   * guessing the alphabet, we'll read the first MSA in text mode.)
   */
  if (! afp->abc) {
    for (x = 0; x < 128; x++) 
      afp->inmap[x] = (isgraph(x) ? x : eslDSQ_ILLEGAL); 
    afp->inmap[0] = '?';
  }
  msafile_SetInmap(afp); /* this does any remaining format-specific inmap configuration */

  /* Determine the alphabet; set <abc>. (<abc> == NULL means text mode.)  */
#ifdef eslAUGMENT_ALPHABET
  if (byp_abc && *byp_abc)	/* Digital mode, and caller provided the alphabet */
    { 
      abc       = *byp_abc;
      alphatype = abc->type;
    } 
  else if (byp_abc)		/* Digital mode, and caller wants us to guess and create an alphabet */
    {
      if ((status = eslx_msafile_GuessAlphabet(afp, &alphatype)) != eslOK) goto ERROR; /* includes normal EAMBIGUOUS, EFORMAT, ENODATA */
      if ( (abc = esl_alphabet_Create(alphatype)) == NULL) { status = eslEMEM; goto ERROR; }
    }    
#endif
  if (abc && ! byp_abc) ESL_EXCEPTION(eslEINCONCEIVABLE, "Your version of Easel does not include digital alphabet code."); 
  /* ^^^^^^^^^^^^^^^^^  this test interacts tricksily with the #ifdef above */

  /* Set up the inmap for digital mode. */
#ifdef eslAUGMENT_ALPHABET
  if (abc) {
    for (x = 0; x < 128; x++)
      afp->inmap[x] = abc->inmap[x]; 
    afp->inmap[0] = esl_abc_XGetUnknown(abc);
    msafile_SetInmap(afp); 
  }
#endif

  afp->abc = abc;
  if (esl_byp_IsReturned(byp_abc)) *byp_abc = abc;
  *ret_afp = afp;
  return eslOK;

 ERROR:  /* on normal errors, afp is returned in an error state */
  if (abc && ! esl_byp_IsProvided(byp_abc)) { esl_alphabet_Destroy(abc); }
  if (esl_byp_IsReturned(byp_abc)) *byp_abc = NULL;
  if (status == eslENOTFOUND || status == eslFAIL || status == eslEFORMAT || status == eslENODATA || eslEAMBIGUOUS) {
    afp->abc = NULL;
    *ret_afp = afp;
  } else {
    if (afp) eslx_msafile_Close(afp);
    *ret_afp = NULL;
  }
  return status;
}

/* Function:  eslx_msafile_OpenFailure()
 * Synopsis:  Report diagnostics of normal error in opening MSA file, and exit.
 *
 * Purpose:   Report user-directed diagnostics of a normal error in opening
 *            an MSA input. Print a message to <stderr>, then exit. 
 */
void
eslx_msafile_OpenFailure(ESLX_MSAFILE *afp, int status)
{
  fprintf(stderr, "alignment input open failed:\n   %s\n", afp->errmsg); 
  
  eslx_msafile_Close(afp);
  exit(status);
}

/* Function:  esl_msafile_Close()
 * Synopsis:  Close an open <ESLX_MSAFILE>.
 */
void
eslx_msafile_Close(ESLX_MSAFILE *afp)
{
  if (afp) {
    if (afp->bf)        esl_buffer_Close(afp->bf);
    if (afp->ssi)       esl_ssi_Close(afp->ssi);
    if (afp->msa_cache) esl_msa_Destroy(afp->msa_cache);
    free(afp);
  }
}

/* msafile_SetInmap()
 *
 * Finish the configuration of the input map of a newly opened
 * <ESLX_MSAFILE>. The input map is already configured for the
 * appropriate alphabet (including text mode, in which case it's set
 * to map all <isgraph()> characters to themselves), and <inmap[0]> is
 * set to an appropriate "unknown" character to replace any invalid
 * input with.
 *            
 * Now do any remaining format-specific initialization: Any
 * characters that need to be ignored in input are set to
 * <eslDSQ_IGNORE>. Any whitespace characters that need to be
 * recognized (as gap characters, for example) are mapped. Any
 * additional gap characters besides the usual "_-.~" in digital
 * alphabets are mapped.
 *            
 * (In fact this is just a dispatcher to format-specific functions;
 * the above documentation tells you what all of those functions are
 * doing, for their format.)
 */
static int
msafile_SetInmap(ESLX_MSAFILE *afp)
{
  int status;

  switch (afp->format) {
  case eslMSAFILE_AFA:          status = esl_msafile_afa_SetInmap(    afp); break;
  case eslMSAFILE_CLUSTAL:      status = esl_msafile_clustal_SetInmap(afp); break;
  case eslMSAFILE_CLUSTALLIKE:  status = esl_msafile_clustal_SetInmap(afp); break;
  default:                      ESL_EXCEPTION(eslEINCONCEIVABLE, "no such msa file format");
  }
  return status;
}

/*---------------- end, open/close, text mode -------------------*/


/*****************************************************************
 *# 2. Guessing file format.
 *****************************************************************/


/* Function:  eslx_msafile_GuessFileFormat()
 * Synopsis:  Guess the MSA file format of an open buffer.
 *
 * Purpose:   Peek into an open buffer, and try to determine what
 *            alignment file format (if any) its input is in. If a
 *            format can be determined, return <eslOK> and set
 *            <*ret_fmtcode> to the format code.  If not, return
 *            <eslEFORMAT> and set <*ret_fmtcode> to
 *            <eslMSAFILE_UNKNOWN.  In either case, the buffer <bf> is
 *            restored to its original position upon return.
 *
 *            If the <bf> corresponds to an open file with a file
 *            name, we attempt to use the suffix as a clue. Suffix
 *            rules for alignment files are as follows:
 *                 | Stockholm     |  .sto .sth .stk |
 *                 | Aligned FASTA |  .afa .afasta   |
 *                 | CLUSTAL       |  .aln           |
 *                 | Pfam          |  .pfam          |
 *                 | MSF           |  .msf           |
 *                 | A2M           |  .a2m           | 
 *                 | SELEX         |  .slx .selex    |   
 *
 * Args:      bf          - the open buffer to read input from
 *            ret_fmtcode - RETURN: format code that's determined
 *
 * Returns:   <eslOK> on success, and <*ret_fmtcode> contains the format code.
 *            <eslEFORMAT> if format can't be guessed, and <*ret_fmtcode> contains
 *            <eslMSAFILE_UNKNOWN>.
 *
 * Throws:    (no abnormal error conditions)
 */
int
eslx_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode)
{
  esl_pos_t  initial_offset;
  char      *p;
  esl_pos_t  n;
  int        fmt_bysuffix    = eslMSAFILE_UNKNOWN;
  int        fmt_byfirstline = eslMSAFILE_UNKNOWN;
  int        status;

  /* As we start, save parser status:
   *   remember the offset where we started (usually 0, but not necessarily)
   *   set an anchor to be sure that this offset stays in the buffer's memory
   */
  initial_offset = esl_buffer_GetOffset(bf);
  esl_buffer_SetAnchor(bf, initial_offset);

  /* First we try to guess based on the filename suffix.
   * (if there's a filename, and if it has a suffix, anyway.)
   */
  if (bf->filename)
    {
      esl_file_Extension(bf->filename, 0, &p, &n);
      if (esl_memstrcmp(p, n, ".gz"))
	esl_file_Extension(bf->filename, 3, &p, &n);

      if (p)
	{
	  if      (esl_memstrcmp(p, n, ".sto"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".stk"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".sth"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".pfam"))   fmt_bysuffix = eslMSAFILE_PFAM;
	  else if (esl_memstrcmp(p, n, ".a2m"))    fmt_bysuffix = eslMSAFILE_A2M;
	  else if (esl_memstrcmp(p, n, ".psi"))    fmt_bysuffix = eslMSAFILE_PSIBLAST;
	  else if (esl_memstrcmp(p, n, ".slx"))    fmt_bysuffix = eslMSAFILE_SELEX;
	  else if (esl_memstrcmp(p, n, ".selex"))  fmt_bysuffix = eslMSAFILE_SELEX;
	  else if (esl_memstrcmp(p, n, ".afa"))    fmt_bysuffix = eslMSAFILE_AFA;
	  else if (esl_memstrcmp(p, n, ".afasta")) fmt_bysuffix = eslMSAFILE_AFA;
	  else if (esl_memstrcmp(p, n, ".aln"))    fmt_bysuffix = eslMSAFILE_CLUSTAL;
	}
    }

  /* Second, we peek at the first non-blank line of the file.
   * Multiple sequence alignment files are generally identifiable by 
   * a token on this line.
   */
  /* Skip blank lines, get first non-blank one */
  do   { 
    status = esl_buffer_GetLine(bf, &p, &n);
  } while (status == eslOK && esl_memspn(p, n, " \t\r\n") == n);

  if      (esl_memstrpfx(p, n, "# STOCKHOLM"))                      fmt_byfirstline = eslMSAFILE_STOCKHOLM;
  else if (esl_memstrpfx(p, n, ">"))                                fmt_byfirstline = eslMSAFILE_AFA;
  else if (esl_memstrpfx(p, n, "CLUSTAL"))                          fmt_byfirstline = eslMSAFILE_CLUSTAL;
  else if (esl_memstrcontains(p, n, "multiple sequence alignment")) fmt_byfirstline = eslMSAFILE_CLUSTALLIKE;

  if      (fmt_byfirstline != eslMSAFILE_UNKNOWN && fmt_bysuffix == eslMSAFILE_UNKNOWN)  *ret_fmtcode = fmt_byfirstline;
  else if (fmt_byfirstline == eslMSAFILE_UNKNOWN && fmt_bysuffix != eslMSAFILE_UNKNOWN)  *ret_fmtcode = fmt_bysuffix;     
  else if (fmt_byfirstline == fmt_bysuffix)                                              *ret_fmtcode = fmt_byfirstline;   /* including MSAFILE_UNKNOWN */
  else if (fmt_byfirstline == eslMSAFILE_CLUSTALLIKE)                                    *ret_fmtcode = fmt_byfirstline;   /* MUSCLE files can have any suffix */
  else    *ret_fmtcode = fmt_bysuffix;

  /* As we return, restore parser status:
   *  put it back where it was when we started;
   *  clear the anchor that made sure that position stayed in memory.
   */
  esl_buffer_SetOffset  (bf, initial_offset);
  esl_buffer_RaiseAnchor(bf, initial_offset);
  return ((*ret_fmtcode == eslMSAFILE_UNKNOWN) ? eslEFORMAT: eslOK);
}


/* Function:  eslx_msafile_EncodeFormat()
 * Synopsis:  Convert text string to an MSA file format code.
 *
 * Purpose:   Given a text string, match it case-insensitively
 *            against a list of possible formats, and return the
 *            appropriate MSA file format code. For example,
 *            <esl_msa_EncodeFormat("Stockholm")> returns
 *            <eslMSAFILE_STOCKHOLM>.
 *            
 *            If the format is unrecognized, return
 *            <eslMSAFILE_UNKNOWN>.
 *            
 * Note:      Keep in sync with <esl_sqio_EncodeFormat()>, 
 *            which decodes all possible sequence file formats,
 *            both unaligned and aligned.           
 */
int
eslx_msafile_EncodeFormat(char *fmtstring)
{
  if (strcasecmp(fmtstring, "stockholm")   == 0) return eslMSAFILE_STOCKHOLM;
  if (strcasecmp(fmtstring, "pfam")        == 0) return eslMSAFILE_PFAM;
  if (strcasecmp(fmtstring, "a2m")         == 0) return eslMSAFILE_A2M;
  if (strcasecmp(fmtstring, "psiblast")    == 0) return eslMSAFILE_PSIBLAST;
  if (strcasecmp(fmtstring, "selex")       == 0) return eslMSAFILE_SELEX;
  if (strcasecmp(fmtstring, "afa")         == 0) return eslMSAFILE_AFA;
  if (strcasecmp(fmtstring, "clustal")     == 0) return eslMSAFILE_CLUSTAL;
  if (strcasecmp(fmtstring, "clustallike") == 0) return eslMSAFILE_CLUSTALLIKE;
  return eslMSAFILE_UNKNOWN;
}


/* Function:  eslx_msafile_DecodeFormat()
 * Synopsis:  Convert internal file format code to text string.
 *
 * Purpose:   Given an internal file format code <fmt> 
 *            (<eslMSAFILE_STOCKHOLM>, for example), returns
 *            a string suitable for printing ("Stockholm",
 *            for example).
 *            
 * Returns:   a pointer to a static description string.
 * 
 * Throws:    If code isn't valid, throws an <eslEINCONCEIVABLE> exception 
 *            internally, and returns <NULL>. 
 *            
 * Note:      Keep in sync with <esl_sqio_DecodeFormat()>.
 */
char *
eslx_msafile_DecodeFormat(int fmt)
{
  switch (fmt) {
  case eslMSAFILE_UNKNOWN:     return "unknown";
  case eslMSAFILE_STOCKHOLM:   return "Stockholm";
  case eslMSAFILE_PFAM:        return "Pfam";
  case eslMSAFILE_A2M:         return "UCSC A2M";
  case eslMSAFILE_PSIBLAST:    return "PSI-BLAST";
  case eslMSAFILE_SELEX:       return "SELEX";
  case eslMSAFILE_AFA:         return "aligned FASTA";
  case eslMSAFILE_CLUSTAL:     return "Clustal";
  case eslMSAFILE_CLUSTALLIKE: return "Clustal-like";
  default:                     break;
  }
  esl_exception(eslEINVAL, __FILE__, __LINE__, "no such msa format code %d\n", fmt);
  return NULL;
}

/*---------------- end, file format utilities -------------------*/


/*****************************************************************
 *# 3. Guessing alphabet
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  eslx_msafile_GuessAlphabet()
 * Synopsis:  Guess what kind of sequences the MSA file contains.
 *
 * Purpose:   Guess the alphabet of the sequences in the open
 *            <ESL_MSAFILE> <afp> -- <eslDNA>, <eslRNA>, or <eslAMINO> --
 *            based on the composition of the next MSA in the
 *            file. Usually this would be the first MSA, because we
 *            would call <esl_msafile_GuessAlphabet()> immediately
 *            after opening a new MSA file.
 *            
 * Returns:   Returns <eslOK> on success, and <*ret_type> is set
 *            to <eslDNA>, <eslRNA>, or <eslAMINO>. Now <afp> 
 *            might hold the next MSA in its cache, where the next 
 *            <esl_msafile*Read()> call will retrieve it.
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to
 *            <eslUNKNOWN> if the first alignment in the file contains
 *            no more than ten residues total, or if its alphabet
 *            cannot be reliably guessed (it contains IUPAC degeneracy
 *            codes, but no amino acid specific residues). <afp>
 *            holds the text-mode alignment in cache.
 * 
 *            Returns <eslEFORMAT> if a parse error is encountered
 *            in trying to read the alignment file. <afp->errmsg>
 *            is set to a useful error message if this occurs; 
 *            <*ret_type> is set to <eslUNKNOWN> and no alignment 
 *            is cached.
 *
 *            Returns <eslENODATA> if the file is empty and no
 *            alignment was found. <*ret_type> is set to <eslUNKNOWN>,
 *            <afp->errmsg> is blank, and the cache is empty.
 *
 * Throws:    <eslEMEM> - an allocation error.
 *            <eslESYS> - a system call such as <fread()> failed.
 *            <eslEINCONCEIVABLE> - "impossible" corruption, internal bug
 *
 * Xref:      SRE:J1/62.
 */
int
eslx_msafile_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type)
{
  ESL_MSA *msa = NULL;
  int      status;

  /* If <afp> is already digital mode, we already know the type (so why
   * are we being called?) If do_digital is TRUE, msafp->abc is
   * non-NULL: 
   */
  if (afp->abc) { *ret_type = afp->abc->type; return eslOK; } /* that was easy */

  /* If there's already an MSA cached, we've already called
   * GuessAlphabet(); don't read another one, or we'll overwrite the
   * first.
   */
  if (afp->msa_cache) return esl_msa_GuessAlphabet(afp->msa_cache, ret_type);

  /* Read and cache the first alignment in input.  */
  status = eslx_msafile_Read(afp, &msa);
  if      (status == eslEOF)  return eslENODATA;
  else if (status != eslOK)   return status;
  afp->msa_cache = msa;

  /* And over to msa_GuessAlphabet() for the decision. */
  return esl_msa_GuessAlphabet(msa, ret_type);
}
#endif /*eslAUGMENT_ALPHABET*/
/*----------- end, utilities for alphabets ----------------------*/



/*****************************************************************
 *# 4. Reading MSAs from input
 *****************************************************************/


/* Function:  eslx_msafile_Read()
 * Synopsis:  Read next MSA from input.
 *
 * Purpose:   Reads the next MSA from open MSA input <afp>, and return it in 
 *            <*ret_msa>.
 *
 * Args:      afp      - open alignment input stream
 8            *ret_msa - RETURN: alignment
 *
 * Returns:   <eslOK> on success. 
 *
 *            <eslEFORMAT> on a parse error, and <afp->errmsg> is set
 *            to a user-directed error message; <*ret_msa> is <NULL>.
 *
 *            If no alignment is found at all, returns <eslEOF>,
 *            and <afp->errmsg> is blank; <*ret_msa> is <NULL>.
 *
 *            On normal error, <afp> and the return status code may be
 *            passed to <eslx_msafile_ReadFailure()> to print diagnostics
 *            to <stderr> (including input source information and line
 *            number) and exit.
 *
 * Throws:    <eslEMEM> - an allocation failed.
 *            <eslESYS> - a system call such as fread() failed
 *            <eslEINCONCEIVABLE> - "impossible" corruption 
 */
int
eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA *msa    = NULL;
  int      status = eslOK;
  
  switch (afp->format) {
  case eslMSAFILE_AFA:          status = esl_msafile_afa_Read    (afp, &msa); break;
  case eslMSAFILE_CLUSTAL:      status = esl_msafile_clustal_Read(afp, &msa); break;
  case eslMSAFILE_CLUSTALLIKE:  status = esl_msafile_clustal_Read(afp, &msa); break;
  default:                      ESL_EXCEPTION(eslEINCONCEIVABLE, "no such msa file format");
  }
  
  *ret_msa = msa;
  return status;
}

/* Function:  eslx_msafile_ReadFailure()
 * Synopsis:  Report diagnostics of a normal error in parsing MSA file, and exit.
 *
 * Purpose:   Report user-directed diagnostics of a normal error from
 *            parsing an MSA file.  Output the error message to
 *            <stderr>, along with information about what we were
 *            parsing (filename, if it was a file) and where we were
 *            in the input (linenumber, if we know it). This
 *            information is all available in <afp>. Then close <afp>
 *            and exit with the <status> provided by the caller.
 *
 * Args:      afp    - open ESLX_MSAFILE, containing information about
 *                     the error and the input source.
 *            status - exit status; generally eslEFORMAT. 
 *
 * Returns:   no return. Exits here with <status>.
 */
void
eslx_msafile_ReadFailure(ESLX_MSAFILE *afp, int status)
{
  switch (status) {
  case eslEFORMAT:  fprintf(stderr, "alignment input parse error: %s\n", afp->errmsg);           break;
  case eslEOF:      fprintf(stderr, "alignment input appears empty?\n");                         break;
  default:          fprintf(stderr, "alignment input read error; unexpected code %d\n", status); break;
  }
  
  switch (afp->bf->mode_is) {
  case eslBUFFER_STREAM:   fprintf(stderr, "   while reading from an input stream (not a file)\n");   break;
  case eslBUFFER_CMDPIPE:  fprintf(stderr, "   while reading through a pipe (not a file)\n");         break;
  case eslBUFFER_FILE:     
  case eslBUFFER_ALLFILE:
  case eslBUFFER_MMAP:     fprintf(stderr, "   while reading file %s\n", afp->bf->filename);          break;
  case eslBUFFER_STRING:   fprintf(stderr, "   while reading from a provided string (not a file)\n"); break;
  default:                 break; 
  }

  if (afp->linenumber > 0) fprintf(stderr, "   at or near line %" PRIu64 "\n", afp->linenumber);
  else                     fprintf(stderr, "   at or near byte %" PRIu64 "\n", esl_buffer_GetOffset(afp->bf));

  eslx_msafile_Close(afp);
  exit(status);
}


/* Function:  eslx_msafile_Decache()
 * Synopsis:  Retrieve already-read MSA from cache.
 *
 * Purpose:   Retrieve a cached MSA from <afp>, instead of reading the
 *            next MSA from input. Return it in <*ret_msa>.
 *            
 *            <eslx_msafile_GuessAlphabet()> reads one MSA from input
 *            to examine its alphabet, then caches that alignment for
 *            the next <Read()> call. All <esl_msafile*_Read()>
 *            functions first check the <afp>'s cache, before reading
 *            a new MSA.
 *
 * Args:      afp     - open alignment input source
 *            ret_msa - RETURN: msa retrieved from cache
 *
 * Returns:   <eslOK> on success; <*ret_msa> contains the cached <msa>,
 *            in text/digital mode appropriate to the <afp>.
 *            
 *            <eslENODATA> if no MSA is cached, and <*ret_afp> is <NULL.>
 *            
 *            <eslEINVAL> if we try to digitize the <msa> (as requested
 *            by <afp>), but one or more sequences contains invalid
 *            characters that can't be digitized. If this happens,
 *            the MSA is left unaltered in the cache, and <afp->errmsg>
 *            is set to a user-directed error message. Now <*ret_msa>
 *            is <NULL>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslECORRUPT> if we try to convert a digital MSA to text
 *            but the digital format is invalid/corrupt.
 */
int
eslx_msafile_Decache(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  int status;

  *ret_msa = NULL;
  if (! afp->msa_cache) return eslENODATA;
  
#ifdef eslAUGMENT_ALPHABET
  if      (afp->abc  && !(afp->msa_cache->flags & eslMSA_DIGITAL)) {
    if ((status = esl_msa_Digitize(afp->abc, afp->msa_cache, afp->errmsg)) != eslOK) return status; 
  }
  else if (! afp->abc && (afp->msa_cache->flags & eslMSA_DIGITAL)) {
    if ((status = esl_msa_Textize(afp->msa_cache)) != eslOK) return status;
  }
#endif
  *ret_msa = afp->msa_cache;
  afp->msa_cache = NULL;
  return eslOK;
}
/*------------ end, reading MSA from ESLX_MSAFILE ---------------*/




/*****************************************************************
 * 5. Writing an MSA to a stream.
 *****************************************************************/

/* Function:  eslx_msafile_Write()
 * Synopsis:  Write an MSA to a stream.
 *
 * Purpose:   Writes alignment <msa> to open stream <fp> in format <fmt>.
 * 
 *            In general, the <msa> is unchanged, but there are some
 *            exceptions. For example, writing an alignment in A2M format
 *            will alter alignment data (marking missing data
 *            symbols on heuristically defined sequence fragments) and
 *            create an <\#=RF> annotation line, if an <msa->rf>
 *            annotation line isn't already present in the <msa>.
 *
 * Args:      fp   - open stream (such as <stdout>)  
 *            msa  - alignment to write
 *            fmt  - format code (such as <eslMSAFILE_STOCKHOLM>)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
eslx_msafile_Write(FILE *fp, ESL_MSA *msa, int fmt)
{
  int status;

  switch (fmt) {
  case eslMSAFILE_CLUSTAL:     status = esl_msafile_clustal_Write(fp, msa, eslMSAFILE_CLUSTAL);     break;
  case eslMSAFILE_CLUSTALLIKE: status = esl_msafile_clustal_Write(fp, msa, eslMSAFILE_CLUSTALLIKE); break;
  default:                     ESL_EXCEPTION(eslEINCONCEIVABLE, "no such msa file format");
  }

  return status;
}

/*-------------- end, writing MSA to stream ---------------------*/


/*****************************************************************
 *# 6. Utilities for specific parsers
 *****************************************************************/

/* Function:  eslx_msafile_ReadLine()
 * Synopsis:  Read next line of input alignment file.
 *
 * Purpose:   Read next line of input <afp>, into its internal
 *            data fields: <afp->line> points to the start of the
 *            line in <afp->bf>, <afp->nline> is its length in
 *            bytes, <afp->lineoffset> is the offset in the input
 *            to the start of the line, and <afp->linenumber> is
 *            the linenumber from <1..N> for <N> total lines in the
 *            input.
 *
 * Args:      <afp> : an open alignment file input
 *
 * Returns:   <eslOK> on success.
 *            
 *            <eslEOF> at EOF. Now <afp->line> is <NULL>, <afp->nline>
 *            is <0>, and <afp->lineoffset> is <0>. <afp->linenumber>
 *            is the total number of lines in the input.
 *            
 * Throws:    <eslEMEM> if an allocation fails.
 *            <eslESYS> if a system call fails, such as fread().
 *            <eslEINCONCEIVABLE> on internal code errors.
 */
int
eslx_msafile_ReadLine(ESLX_MSAFILE *afp)
{
  int status;

  afp->lineoffset = esl_buffer_GetOffset(afp->bf);
  if ((status = esl_buffer_GetLine(afp->bf, &(afp->line), &(afp->nline))) != eslOK) goto ERROR;
  if (afp->linenumber != -1) afp->linenumber++;

 ERROR:
  afp->line       = NULL;
  afp->nline      = 0;
  afp->lineoffset = 0;
  /* leave linenumber alone. on EOF, it is the number of lines in the file, and that might be useful. */
  return status;
}



/*****************************************************************
 * x. Examples.
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET
#ifdef eslMSAFILE_EXAMPLE
/*::cexcerpt::msafile_example::begin::*/
#include <stdio.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

static ESL_OPTIONS options[] = {
  /* name            type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",        0 },
  { "--informat", eslARG_STRING,  NULL,  NULL, NULL,  NULL,  NULL, NULL, "specify the input MSA file is in format <s>", 0 }, 
  { "--dna",      eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",      eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",    eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "example of digital MSA reading using the msafile module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *msafile   = esl_opt_GetArg(go, 1);
  int           fmt       = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET *abc       = NULL;
  ESLX_MSAFILE *afp       = NULL;
  ESL_MSA      *msa       = NULL;
  int           nali      = 0;
  int           status;

  /* If you know the alphabet you want, create it */
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  /* If you know the MSA file format, set it */
  if (esl_opt_IsOn(go, "--informat") &&
      (fmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  /* Open in digital mode. If fmt is unknown, guess format; if abc unknown, guess alphabet */
  if ( (status = eslx_msafile_Open(&abc, msafile, fmt, NULL, &afp)) != eslOK)
    eslx_msafile_OpenFailure(afp, status);
  
  /* Now the MSA's that you read are digital data in msa->ax[] */
  while ((status = eslx_msafile_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, (int) msa->nseq, (int) msa->alen);
      eslx_msafile_Write(stdout, msa, eslMSAFILE_CLUSTAL);
      esl_msa_Destroy(msa);
    }
  if (nali == 0 || status != eslEOF) eslx_msafile_ReadFailure(afp, status);

  esl_alphabet_Destroy(abc);
  eslx_msafile_Close(afp);
  esl_getopts_Destroy(go);
  exit(0);
}
/*::cexcerpt::msafile_example::end::*/
#endif /*eslMSAFILE_EXAMPLE*/
#endif /*eslAUGMENT_ALPHABET*/

/*------------------------ end of examples -----------------------*/
      

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

