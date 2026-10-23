/* Unaligned sequence file i/o.
 * 
 * Contents:
 *    1. Opening an <ESL_SQFILE>, in text or digital mode.
 *    2. Guessing the alphabet
 *    3. Sequence reading (sequential)
 *    4. Sequence writing.
 *    5. Miscellaneous functions
 *    6. Sequence/subsequence fetching, random access [with <ssi>]
 *    7. Sequence database caching.
 *    8. Internal functions (all formats)
 *    9. Internal functions (FASTA)
 *   10. Internal functions (EMBL|UniProt|TrEMBL)
 *   11. Internal functions (GenBank)
 *   12. Internal functions (daemon)
 *   13. Internal functions (hmmpgmd)
 *   14. Benchmark driver.
 *   15. Unit tests.
 *   16. Test driver.
 *   17. Examples.
 * 
 * This module shares remote evolutionary homology with Don Gilbert's
 * seminal, public domain ReadSeq package - though the last common
 * ancestor was circa 1991 and few recognizable vestiges
 * remain. Thanks Don!
 *
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>		// POSIX strcasecmp() 
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dsq.h"
#include "esl_mem.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_sq.h"
#include "esl_ssi.h"

static void sqfile_init(ESL_SQFILE *sqfp);
static int  sqfile_open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
static int  sqfile_guess_format(ESL_SQFILE *sqfp, int *ret_fmt);

/* Internal functions shared by parsers. */
static int  loadmem  (ESL_SQFILE *sqfp);
static int  loadbuf  (ESL_SQFILE *sqfp);
static int  nextchar (ESL_SQFILE *sqfp, char *ret_c);
static int  seebuf   (ESL_SQFILE *sqfp, int64_t maxn, int64_t *opt_nres, int64_t *opt_endpos);
static void addbuf   (ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nres);
static void skipbuf  (ESL_SQFILE *sqfp, int64_t nskip);
static int  read_nres(ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nskip, int64_t nres, int64_t *opt_actual_nres);
static int  skip_whitespace(ESL_SQFILE *sqfp);
static int  convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa);

/* EMBL format; also UniProt, TrEMBL */
static void config_embl(ESL_SQFILE *sqfp);
static void inmap_embl (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_embl(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_embl   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* GenBank format; also DDBJ */
static void config_genbank(ESL_SQFILE *sqfp);
static void inmap_genbank (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_genbank   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* FASTA format */
static void config_fasta(ESL_SQFILE *sqfp);
static void inmap_fasta (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_fasta   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* daemon format */
static void config_daemon(ESL_SQFILE *sqfp);
static void inmap_daemon (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  end_daemon   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* HMMPGMD format */
static int  fileheader_hmmpgmd(ESL_SQFILE *sqfp);


/*****************************************************************
 * 1. Opening an <ESL_SQFILE>, in text or digital mode.
 *****************************************************************/ 

/* Function:  esl_sqfile_Open()
 * Synopsis:  Open a sequence file <filename> for reading. 
 *
 * Purpose:   Open a sequence file <filename> for reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The format of the file is asserted to be <format> (for
 *            example, <eslSQFILE_FASTA>). If <format> is
 *            <eslSQFILE_UNKNOWN> then we attempt to autodetect the
 *            file format.
 *            
 *            If <env> is non-NULL, it is the name of an environment
 *            variable that contains a colon-delimited list of
 *            directories in which we may find this <filename>.
 *            For example, if we had 
 *            <setenv SEQDB /databases/uniprot:/databases/genomes/>
 *            in the environment, a database search application
 *            could pass "SEQDB" as <env>.
 *
 *            There are two special cases for <filename>. If
 *            <filename> is "-", the sequence data are read from a
 *            <STDIN> pipe. If <filename> ends in ".gz", the file is
 *            assumed to be compressed with <gzip>, and it is opened
 *            by a pipe from <gzip -dc>. Reading gzip files only works
 *            on POSIX-compliant systems that have pipes
 *            (specifically, the POSIX.2 popen() call). 
 *            
 * Returns:   <eslOK> on success, and <*ret_sqfp> points to a new
 *            open <ESL_SQFILE>. Caller deallocates this object with
 *            <esl_sqfile_Close()>. 
 *            
 *            Returns <eslENOTFOUND> if <filename> can't be opened for
 *            any reason.
 *            
 *            Returns <eslEFORMAT> if the file is empty, or if format
 *            autodetection was attempted and the format can't be
 *            determined.
 *
 *            On any error condition, <*ret_sqfp> is returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqfile_Open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
{
  return sqfile_open(filename, format, env, ret_sqfp);
}



/* Function:  esl_sqfile_OpenDigital()
 * Synopsis:  Open an <ESL_SQFILE> for digital input.
 *
 * Purpose:   Same as <esl_sqfile_Open()>, but we will expect all
 *            sequence input to conform to the digital alphabet <abc>.
 *            
 *            Normally, after opening the sequence file in digital
 *            mode, you'd read sequence into a digital <ESL_SQ>.
 *            However, you don't actually have to. The state of the
 *            <ESL_SQ> controls how the input is stored; the state of
 *            the <ESL_SQFILE> controls how the input is validated.
 *            
 * Returns:   <eslOK> on success, and <*ret_sqfp> points to a new
 *            open <ESL_SQFILE>.
 *            
 *            Returns <eslENOTFOUND> if <filename> can't be opened.
 *            Returns <eslEFORMAT> if the file is empty, or if
 *            autodetection is attempted and the format can't be
 *            determined.  On any error conditions, <*ret_sqfp> is
 *            returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
{
  int status;

  if ((status = sqfile_open(filename, format, env, ret_sqfp)) != eslOK) goto ERROR;
  if ((status = esl_sqfile_SetDigital(*ret_sqfp, abc))        != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  esl_sqfile_Close(*ret_sqfp);
  *ret_sqfp = NULL;
  return status;
}


/* Function:  esl_sqfile_SetDigital()
 * Synopsis:  Set an open <ESL_SQFILE> to read in digital mode.
 *
 * Purpose:   Given an <ESL_SQFILE> that's already been opened,
 *            configure it to expect subsequent input to conform
 *            to the digital alphabet <abc>.
 *            
 *            Calling <esl_sqfile_Open(); esl_sqfile_SetDigital()> is
 *            equivalent to <esl_sqfile_OpenDigital()>. The two-step
 *            version is useful when you need a
 *            <esl_sqfile_GuessAlphabet()> call in between, guessing
 *            the file's alphabet in text mode before you set it to
 *            digital mode.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
{
  if (! esl_sqio_IsAlignment(sqfp->format))
    {
      switch (sqfp->format) {
      case eslSQFILE_EMBL:       inmap_embl   (sqfp, abc->inmap); break;
      case eslSQFILE_UNIPROT:    inmap_embl   (sqfp, abc->inmap); break;
      case eslSQFILE_GENBANK:    inmap_genbank(sqfp, abc->inmap); break;
      case eslSQFILE_DDBJ:       inmap_genbank(sqfp, abc->inmap); break;
      case eslSQFILE_FASTA:      inmap_fasta  (sqfp, abc->inmap); break;
      case eslSQFILE_DAEMON:     inmap_daemon (sqfp, abc->inmap); break;
      case eslSQFILE_HMMPGMD:    inmap_fasta  (sqfp, abc->inmap); break;
      default:                   ESL_EXCEPTION(eslEINCONCEIVABLE, "no such seqfile format");
      }
    }
  else
    esl_msafile_SetDigital(sqfp->afp, abc);

  sqfp->do_digital = TRUE;
  sqfp->abc        = abc;
  return eslOK;
}


/* Function:  esl_sqfile_Close()
 * Synopsis:  Close a sequence file.
 *
 * Purpose:   Closes an open <sqfp>.
 *
 * Returns:   (void).
 */
void
esl_sqfile_Close(ESL_SQFILE *sqfp)
{
  if (sqfp)
    {
      if (! sqfp->do_gzip && ! sqfp->do_stdin && sqfp->fp) fclose(sqfp->fp);
#ifdef HAVE_POPEN
      if (sqfp->do_gzip && sqfp->fp)  pclose(sqfp->fp);
#endif
      if (sqfp->balloc > 0) free(sqfp->buf);
      if (sqfp->ssi)        esl_ssi_Close(sqfp->ssi);
      if (sqfp->afp)        esl_msafile_Close(sqfp->afp);
      if (sqfp->msa)        esl_msa_Destroy(sqfp->msa);

      free(sqfp->filename);
      free(sqfp->mem);
      free(sqfp->ssifile);
      free(sqfp);
    }
  return;
}

/* sqfile_init()
 * Initialize a newly opened ESL_SQFILE.
 * Just to pull a verbose block of routine code out of sqfile_open().
 * Matches order of elements in ESL_SQFILE as declared in esl_sqio.h.
 */
static void
sqfile_init(ESL_SQFILE *sqfp)
{
  sqfp->filename         = NULL;
  sqfp->fp               = NULL;
  sqfp->errbuf[0]        = '\0';

  sqfp->do_digital       = FALSE;
  sqfp->do_gzip          = FALSE;
  sqfp->do_stdin         = FALSE;
  sqfp->do_buffer        = FALSE;
  sqfp->abc              = NULL;

  sqfp->format           = eslSQFILE_UNKNOWN;
  sqfp->is_linebased     = FALSE;
  sqfp->eof_is_ok        = FALSE;
  sqfp->parse_header     = NULL;
  sqfp->parse_end        = NULL;

  sqfp->mem              = NULL;
  sqfp->allocm           = 0;
  sqfp->mn               = 0;
  sqfp->mpos             = 0;
  sqfp->moff             = -1;
  sqfp->is_recording     = FALSE;

  sqfp->buf              = NULL;
  sqfp->boff             = 0;
  sqfp->balloc           = 0;
  sqfp->nc               = 0;
  sqfp->bpos             = 0;
  sqfp->L                = 0;
  sqfp->linenumber       = 1;
  sqfp->bookmark_offset  = 0;
  sqfp->bookmark_linenum = 0;

  sqfp->afp              = NULL;
  sqfp->msa              = NULL;
  sqfp->idx              = -1;

  sqfp->ssifile          = NULL;
  sqfp->rpl              = -1; 
  sqfp->bpl              = -1; 
  sqfp->prvrpl           = -1; 
  sqfp->prvbpl           = -1; 
  sqfp->currpl           = -1;
  sqfp->curbpl           = -1;
  sqfp->maxrpl           = -1;
  sqfp->maxbpl           = -1;
  sqfp->maxgap           = -1;
  sqfp->ssi              = NULL;
}

/* sqfile_open():
 * This is the routine that actually opens an ESL_SQFILE.
 * esl_sqfile_Open() and esl_sqfile_OpenDigital() are
 * small wrappers around it.
 */
static int
sqfile_open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp   = NULL;
  char       *cmd    = NULL;
  int         n;
  int         status;

  /* Allocate and initialize the ESL_SQFILE structure
   */
  ESL_ALLOC(sqfp, sizeof(ESL_SQFILE));
  sqfile_init(sqfp);
  sqfp->format = format;  // <format> itself is never used again; only sqfp->format from now on.

  /* Find path to <filename> (watching out for "-" meaning stdin),
   * checking <env> variable for a colon-delimited list of dirs if
   * <filename> isn't in current working dir.
   *
   * Store path to file we're going to open in <sqfp->filename>.
   */
  if (strcmp(filename, "-") == 0 || esl_FileExists(filename))
    { if ((status = esl_strdup(filename, -1, &(sqfp->filename)))             != eslOK) goto ERROR; }
  else
    { if ((status = esl_FileEnvOpen(filename, env, NULL, &(sqfp->filename))) != eslOK) goto ERROR; }
  n = strlen(sqfp->filename);

  /* Open the file as a stream, and if <format> hasn't been asserted,
   * peek into it to guess the file format.
   * 
   * If we've already been told it's an MSA format, skip all of this:
   * the MSA parser will handle opening the stream instead, including
   * stdin and .gz files.
   */
  if (! esl_sqio_IsAlignment(sqfp->format))
    {
      if (strcmp(sqfp->filename, "-") == 0)  /* stdin special case */
        {
          sqfp->fp       = stdin;
          sqfp->do_stdin = TRUE;
        }
      else if ((sqfp->fp = fopen(sqfp->filename, "r")) == NULL)
        { status = eslENOTFOUND; goto ERROR; }

#ifdef HAVE_POPEN
      /* Deal with the .gz special case: to popen(), "success" only means
       * it found and executed gzip -dc.  If gzip -dc doesn't find our
       * file, popen() still blithely returns success, so we have to be
       * sure the file exists. That's why we fopen()'ed it above, only to
       * close it and popen() it here.
       */     
      if (n > 3 && strcmp(sqfp->filename+n-3, ".gz") == 0)
        {
          fclose(sqfp->fp);
          sqfp->fp = NULL;
          if ((status = esl_sprintf(&cmd, "gzip -dc %s", sqfp->filename)) != eslOK) goto ERROR; 
          if ((sqfp->fp = popen(cmd, "r")) == NULL)   // FIXME: this is a security problem; see plan
            { status = eslENOTFOUND; goto ERROR; }
          sqfp->do_gzip = TRUE;
        }
#endif

      /* If we don't know the format yet, try to autodetect it now. */
      if (sqfp->format == eslSQFILE_UNKNOWN)
        {
          status = sqfile_guess_format(sqfp, &(sqfp->format));
          if (status != eslOK && status != eslEFORMAT) goto ERROR; // <format> might still be eslSQFILE_UNKNOWN: MSA files for example
        }

      /* If the format is still unknown, it may be an MSA file.  The
       * msafile module is capable of autodetecting format even in a .gz
       * or stdin pipe, but the stuff above has already read from these
       * nonrewindable sources, trying to guess an unaligned format.  We
       * could open a second .gz pipe, but that's ugly; and in any case,
       * we can't rewind stdin. Eventually, this will get resolved, by
       * having sqio open an ESL_BUFFER, then doing an
       * esl_msafile_OpenBuffer() if we need to hand control to the
       * msafile module. For now, sqio is already documented to be
       * unable to autodetect MSA file formats in stdin or .gz pipes,
       * so leave it that way.
       */
      if (sqfp->format == eslSQFILE_UNKNOWN && (sqfp->do_gzip || sqfp->do_stdin)) 
        { status = eslEFORMAT; goto ERROR; }
    }

  /* If format is definitely an MSA, open it through the msafile interface.
   * Or, if format is still unknown, try to open the file as an MSA file,
   * using msafile autodetection. 
   */
  if (sqfp->format == eslSQFILE_UNKNOWN || esl_sqio_IsAlignment(sqfp->format))
    {
      status = esl_msafile_Open(NULL, sqfp->filename, NULL, sqfp->format, NULL, &(sqfp->afp));
      if (status != eslOK) { status = eslEFORMAT; goto ERROR; } // This was our last attempt. Failure to open == failure to detect format
      sqfp->format = sqfp->afp->format;
    }
  if (sqfp->format == eslSQFILE_UNKNOWN) { status = eslEFORMAT; goto ERROR; }

  /* Configure the <sqfp>'s parser and inmaps for this format. */
  if ( ! esl_sqio_IsAlignment(sqfp->format)) 
    {
      switch (sqfp->format) {
      case eslSQFILE_EMBL:     config_embl(sqfp);    inmap_embl(sqfp,    NULL);   break;
      case eslSQFILE_UNIPROT:  config_embl(sqfp);    inmap_embl(sqfp,    NULL);   break;
      case eslSQFILE_GENBANK:  config_genbank(sqfp); inmap_genbank(sqfp, NULL);   break;
      case eslSQFILE_DDBJ:     config_genbank(sqfp); inmap_genbank(sqfp, NULL);   break;
      case eslSQFILE_FASTA:    config_fasta(sqfp);   inmap_fasta(sqfp,   NULL);   break;
      case eslSQFILE_DAEMON:   config_daemon(sqfp);  inmap_daemon(sqfp,  NULL);   break;
      case eslSQFILE_HMMPGMD:  config_fasta(sqfp);   inmap_fasta(sqfp,   NULL);   break;
      default: status = eslEFORMAT; goto ERROR;
      }

      /* Preload the first line or chunk of file. */
      status = loadbuf(sqfp);
      if      (status == eslEOF) { status = eslEFORMAT; goto ERROR; }
      else if (status != eslOK)  { goto ERROR; }

      /* hmmpgmd is a special case: we need to skip first line before parsing it.
       * generalize that a little: this could be a section for parsing a file header,
       * and leaving the buf positioned at the first char of the first record
       * (just as expected if there's no file header)
       */
      switch (sqfp->format) {
      case eslSQFILE_HMMPGMD:   status = fileheader_hmmpgmd(sqfp); break;
      default:                  status = eslOK;                    break;
      }
      if (status != eslOK) goto ERROR;
    }
  else
    {
      sqfp->is_linebased = TRUE;
      sqfp->eof_is_ok    = FALSE; // no-op for msa's 
      sqfp->parse_header = NULL;  // no-op for msa's 
      sqfp->parse_end    = NULL;  // no-op for msa's 
    }

  free(cmd);
  *ret_sqfp = sqfp;
  return eslOK;

 ERROR:
  free(cmd);
  esl_sqfile_Close(sqfp); 
  *ret_sqfp = NULL;
  return status;
}


/* sqfile_guess_format()
 * Guess the format of a newly opened <ESL_SQFILE>.
 *
 * Try to guess the sequence file format of <sqfp>, and return the
 * format code in <*ret_fmt>.
 *            
 * First we attempt to guess based on the <filename>'s suffix. <*.fa>,
 * <*.faa>, <*.fna> are assumed to be in FASTA format; <*.gb> is
 * assumed to be in GenBank format.
 *            
 * If that fails, we attempt to guess based on peeking at the first
 * nonblank line of <filename>. If the line starts with $>$, we assume
 * FASTA format; if the line starts with <ID>, we assume EMBL format;
 * if the line starts with <LOCUS> or it contains the string <Genetic
 * Sequence Data Bank> we assume GenBank format.
 *            
 * If that fails too, return an <eslEFORMAT> error, and
 * <*ret_fmt> is set to <eslSQFILE_UNKNOWN>.
 *            
 * Returns:   <eslOK> on success, and <*ret_fmt> contains
 *            a valid sequence file format code, such as 
 *            <eslSQFILE_FASTA>.
 *            
 *            <eslEFORMAT> if we opened <filename> but it contains no
 *            nonblank lines, or if we peeked at the first nonblank
 *            line and still couldn't guess the format; <*ret_fmt> is
 *            then <eslSQFILE_UNKNOWN>.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 */
static int
sqfile_guess_format(ESL_SQFILE *sqfp, int *ret_fmt)
{
  int   n             = strlen(sqfp->filename);
  const char *sfx     = NULL;
  int   is_gzip       = FALSE;
  int   nheader_fasta = 0;
  int   nsfx;
  int   status;

  /* Is <filename> gzip'ed? Look at suffix. */
  if (n > 3 && strcmp(sqfp->filename+n-3, ".gz") == 0) is_gzip = TRUE;

  /* Locate the suffix that might indicate format (ignore .gz) */
  for (nsfx = 0, sfx = sqfp->filename + n - 1 - (is_gzip ? 3 : 0);
       sfx != sqfp->filename && *sfx != '.' && *sfx != eslDIRSLASH; 
       sfx--) 
    nsfx++;

  /* now sfx points either to filename (we didn't find a suffix),
   * or to a directory separator '/' (also means we didn't find a suffix),
   * or to the . of the suffix,
   * and nsfx is the suffix length exclusive of the . 
   */
  
  /* Attempt to guess file format based on file name suffix.
   * Note that esl_memstrcmp() returns TRUE|FALSE for match|no match.
   */
  if (nsfx)
    {
      if      (esl_memstrcmp(sfx, nsfx+1, ".fa"))     { *ret_fmt = eslSQFILE_FASTA;      return eslOK; }
      else if (esl_memstrcmp(sfx, nsfx+1, ".faa"))    { *ret_fmt = eslSQFILE_FASTA;      return eslOK; }
      else if (esl_memstrcmp(sfx, nsfx+1, ".fna"))    { *ret_fmt = eslSQFILE_FASTA;      return eslOK; }
      else if (esl_memstrcmp(sfx, nsfx+1, ".gb"))     { *ret_fmt = eslSQFILE_GENBANK;    return eslOK; }
      else if (esl_memstrcmp(sfx, nsfx+1, ".afa"))    { *ret_fmt = eslSQFILE_UNKNOWN;    return eslEFORMAT; }  // special case: don't open aligned FASTA MSA as unaligned FASTA.
      else if (esl_memstrcmp(sfx, nsfx+1, ".afasta")) { *ret_fmt = eslSQFILE_UNKNOWN;    return eslEFORMAT; }  // ditto, .afasta
      else if (esl_memstrcmp(sfx, nsfx+1, ".a2m"))    { *ret_fmt = eslSQFILE_UNKNOWN;    return eslEFORMAT; }  // ditto, .a2m
      else if (esl_memstrcmp(sfx, nsfx+1, ".a3m"))    { *ret_fmt = eslSQFILE_UNKNOWN;    return eslEFORMAT; }  // ditto, .a3m
    }
    
  /* If that didn't work, we'll have a peek at the stream; 
   * turn recording on, and set for line based input while we do this.
   */
  if (sqfp->is_recording == -1) ESL_XEXCEPTION(eslEINVAL, "sq file already too advanced");
  sqfp->is_recording = TRUE;
  sqfp->is_linebased = TRUE;
  status = loadbuf(sqfp);
  if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, sqfp->errbuf, "No data found in file");
  else if (status != eslOK)  goto ERROR;
  // now sqfp->buf is a line of the file 

  /* get first nonblank line (and noncomment, in the case of FASTA with metadata*/
  while (esl_str_IsBlank(sqfp->buf) || sqfp->buf[0] == ';')
    {
      if (sqfp->buf[0] == ';') nheader_fasta++;
      status = loadbuf(sqfp);
      if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, sqfp->errbuf, "No non-blank, non-comment data found in file");
      else if (status != eslOK)  goto ERROR;
    } 

  /* formats that can be determined from the first (nonblank, noncomment) line: */
  if      (*(sqfp->buf) == '>')                                                       *ret_fmt = eslSQFILE_FASTA;
  else if (!nheader_fasta && strncmp(sqfp->buf, "ID   ", 5)    == 0)                  *ret_fmt = eslSQFILE_EMBL;
  else if (!nheader_fasta && strncmp(sqfp->buf, "LOCUS   ", 8) == 0)                  *ret_fmt = eslSQFILE_GENBANK;
  else if (!nheader_fasta && strstr(sqfp->buf, "Genetic Sequence Data Bank") != NULL) *ret_fmt = eslSQFILE_GENBANK;
  else *ret_fmt = eslSQFILE_UNKNOWN;

  /* Whether we've been able to determine the format or not, reset
   * sqfp. Turn off recording and linebased mode. Leave sqfp->mem
   * alone, and set sqfp->mpos back to the start of the file's data:
   * after we return, the reader will treat sqfp->mem as preloaded
   * data.
   */
  sqfp->mpos         = 0;
  sqfp->is_recording = FALSE;
  sqfp->is_linebased = FALSE;
  free(sqfp->buf);
  sqfp->buf    = NULL;
  sqfp->balloc = 0;
  return (*ret_fmt == eslSQFILE_UNKNOWN) ? eslEFORMAT : eslOK;

 ERROR:
  *ret_fmt           = eslSQFILE_UNKNOWN;
  sqfp->mpos         = 0;
  sqfp->is_recording = FALSE;
  sqfp->is_linebased = FALSE;
  if (sqfp->buf) { free(sqfp->buf); sqfp->balloc = 0; }
  return status;
}
/*------------------ end, opening an ESL_SQFILE ----------------------*/




/*****************************************************************
 * 2. Guessing the alphabet
 *****************************************************************/

/* Function:  esl_sqfile_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open <ESL_SQFILE>.
 *
 * Purpose:   After opening <sqfp>, attempt to guess what alphabet its
 *            sequences are in, by inspecting up to 4000 residues of
 *            the first sequence in the file, and return this alphabet
 *            type in <*ret_type>.
 *
 *            The limit of a 4000 residue window is so we don't go
 *            reading a whole 100M+ chromosome for this.
 *
 *            Because it only looks at the first sequence in the file,
 *            you can occasionally have a frustrating situation where
 *            the first sequence is some short ambiguous thing and the
 *            second sequence is obviously DNA or protein, and we'll
 *            still return <eslENOALPHABET>. Our programs allow
 *            callers to assert the alphabet with options like
 *            --dna|--rna|--amino so we don't spuriously fail for this
 *            reason on high-throughput work where the input files
 *            have known alphabets.
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set to <eslDNA>,
 *            <eslRNA>, or <eslAMINO>.
 *            
 *            Returns <eslENOALPHABET> if the alphabet can't be
 *            reliably guessed.
 *            
 *            Returns <eslEFORMAT> if a parse error is encountered.
 *            Call <esl_sqfile_GetErrorBuf()> to get a ptr to a
 *            user-directed error message describing the problem,
 *            including the line number on which it was found.
 *            
 *            Returns <eslENODATA> if the file appears to be empty.
 *
 *            On any error, <*ret_type> is <eslUNKNOWN>.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslEINCONCEIVABLE> on unimaginable internal errors.
 */
int
esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)
{
  ESL_SQ *sq = NULL;
  int     status;

  /* Special case: for MSA files, hand this off to msafile_GuessAlphabet. */
  if (esl_sqio_IsAlignment(sqfp->format)) return esl_msafile_GuessAlphabet(sqfp->afp, ret_type);

  /* set the sqfp to record; we'll rewind afterwards and use the recording */
  sqfp->is_recording = TRUE;

  if ((sq = esl_sq_Create()) == NULL) { status = eslEMEM; goto ERROR; }

  status = esl_sqio_ReadWindow(sqfp, 0, 4000, sq);   // ReadWindow() so we don't read a whole 100M+ chromosome for this
  if      (status == eslEOF) { status = eslENODATA; goto ERROR; }
  else if ((status != eslOK) && (status != eslEOD)) goto ERROR;

  if ((status = esl_sq_GuessAlphabet(sq, ret_type)) != eslOK) goto ERROR;

  /* reset the sqfp, so subsequence reads will use the recording next */
  sqfp->mpos         = 0;
  sqfp->linenumber   = 1;
  sqfp->is_recording = FALSE;
  if ((status = loadbuf(sqfp)) != eslOK) ESL_EXCEPTION(status, "buffer load failed, but shouldn't have"); // only affects line-oriented parser: reload buf from mem
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  *ret_type = eslUNKNOWN;
  return status; 
}
/*---------------- end, guessing the alphabet -------------------*/



/*****************************************************************
 * 3. Sequence reading (sequential)
 *****************************************************************/ 

/* Function:  esl_sqio_Read()
 * Synopsis:  Read the next sequence from a file.
 *
 * Purpose:   Reads the next sequence from open sequence file <sqfp> into 
 *            <sq>. Caller provides an allocated and initialized <sq>, which
 *            will be internally reallocated if its space is insufficient.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sq>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if a parse error is encountered.
 *            Call <esl_sqfile_GetErrorBuf()> to get a ptr to a
 *            user-directed error message describing the problem,
 *            including the line number on which it was found.
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_sqio_Read(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int64_t epos;
  int64_t n;
  int     status;

  if (esl_sqio_IsAlignment(sqfp->format))
  {
    ESL_SQ *tmpsq = NULL;
    if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
      { // do we need to load a new alignment?
        esl_msa_Destroy(sqfp->msa);
        status = esl_msafile_Read(sqfp->afp, &(sqfp->msa));
        if (status == eslEFORMAT)
          { // oops, a parse error; upload the error info from afp to sqfp 
            sqfp->linenumber = sqfp->afp->linenumber;
            strcpy(sqfp->errbuf, sqfp->afp->errmsg); // errbufs must be same size 
            return eslEFORMAT;
          }
        if (status != eslOK) return status;
        sqfp->idx = 0;
      }
      
    /* grab next seq from alignment */
    /* this is inefficient; it goes via a temporarily allocated copy of the sequence */
    if ((status = esl_sq_FetchFromMSA(sqfp->msa, sqfp->idx, &tmpsq)) != eslOK) return status;
    esl_sq_GrowTo(sq, tmpsq->n);
    esl_sq_Copy(tmpsq, sq);
    esl_sq_Destroy(tmpsq);
    sqfp->idx++;

    sq->start = 1;
    sq->end   = sq->n;
    sq->C     = 0;
    sq->W     = sq->n;
    sq->L     = sq->n;
    return eslOK;
  }

  /* Main case: read next seq from sqfp's stream */
  if (sqfp->nc == 0) return eslEOF;
  if ((status = sqfp->parse_header(sqfp, sq)) != eslOK) return status; // EMEM, EOF, EFORMAT 

  do {
    if ((status = seebuf(sqfp, -1, &n, &epos)) == eslEFORMAT) return status;
    if (esl_sq_GrowTo(sq, sq->n + n) != eslOK) return eslEMEM;
    addbuf(sqfp, sq, n);
    sqfp->L   += n;
    sq->eoff   = sqfp->boff + epos - 1;
    if (status == eslEOD)     break;
  } while ((status = loadbuf(sqfp)) == eslOK);
    
  if (status == eslEOF)
    {
      if (! sqfp->eof_is_ok) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Unexpected EOF; file truncated?"); 
      if ((status = sqfp->parse_end(sqfp, sq)) != eslOK) return status;
    }
  else if (status == eslEOD)
    {
      sqfp->bpos = epos;
      if ((status = sqfp->parse_end(sqfp, sq)) != eslOK) return status;
    }
  else if (status != eslOK) return status;

  if (sq->dsq != NULL) sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else                 sq->seq[sq->n] = '\0';
  sq->start = 1;
  sq->end   = sq->n;
  sq->C     = 0;
  sq->W     = sq->n;
  sq->L     = sq->n;
  return eslOK;
}


/* Function:  esl_sqio_ReadInfo()
 * Synopsis:  Read sequence info, but not the sequence itself.
 *
 * Purpose:   Read the next sequence from open sequence file <sqfp>,
 *            but don't store the sequence (or secondary structure).
 *            Upon successful return, <s> holds all the available 
 *            information about the sequence -- its name, accession,
 *            description, and overall length <sq->L>. 
 *            
 *            This is useful for indexing sequence files, where
 *            individual sequences might be ginormous, and we'd rather
 *            avoid reading complete seqs into memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sqio_ReadInfo(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int64_t epos;
  int64_t n;
  int     status;

  if (esl_sqio_IsAlignment(sqfp->format))
    {
      ESL_SQ *tmpsq = NULL;
      if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
        { // do we need to load a new alignment?
          esl_msa_Destroy(sqfp->msa);
          status = esl_msafile_Read(sqfp->afp, &(sqfp->msa));
          if (status == eslEFORMAT)
            { // oops, a parse error; upload the error info from afp to sqfp 
              sqfp->linenumber = sqfp->afp->linenumber;
              strcpy(sqfp->errbuf, sqfp->afp->errmsg); // errbufs must be same size! 
              return eslEFORMAT;
            }
          if (status != eslOK) return status;
          sqfp->idx = 0;
        }
      
      /* grab next seq from alignment */
      /* this is inefficient; it goes via a temporarily allocated copy of the sequence */
      if ((status = esl_sq_FetchFromMSA(sqfp->msa, sqfp->idx, &tmpsq)) != eslOK) return status;
      if (tmpsq->dsq != NULL) tmpsq->dsq[1] = eslDSQ_SENTINEL;
      else                    tmpsq->seq[0] = '\0';
      esl_sq_Copy(tmpsq, sq);
      esl_sq_Destroy(tmpsq);
      sqfp->idx++;

      if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL;
      else                 sq->seq[0] = '\0';
      if (sq->ss  != NULL) { free(sq->ss); sq->ss = NULL; }

      sq->n     = 0;
      sq->start = 0;
      sq->end   = 0;
      sq->C     = 0;
      sq->W     = 0;
      return eslOK;
    }

  if (sqfp->nc == 0) return eslEOF;
  if ((status = sqfp->parse_header(sqfp, sq)) != eslOK) return status; /* EOF, EFORMAT */

  sqfp->L       = 0;
  do {
    status = seebuf(sqfp, -1, &n, &epos);
    sqfp->L += n;
    sq->eoff = sqfp->boff + epos - 1;
    if (status == eslEFORMAT) return status;
    if (status == eslEOD)     break;
  } while ((status = loadbuf(sqfp)) == eslOK);
    
  if      (status == eslEOF) 
    {
      if (! sqfp->eof_is_ok) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Unexpected EOF; file truncated?"); 
    }
  else if (status == eslEOD)
    {
      sqfp->bpos = epos;
      if ((status = sqfp->parse_end(sqfp, sq)) != eslOK) return status;
    }
  else if (status != eslOK) return status;
  sq->L = sqfp->L;

  /* Set coord system for an info-only ESL_SQ  */
  if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL;
  else                 sq->seq[0] = '\0';
  if (sq->ss  != NULL) { free(sq->ss); sq->ss = NULL; }
  sq->n     = 0;
  sq->start = 0;
  sq->end   = 0;
  sq->C     = 0;
  sq->W     = 0;
  return eslOK;
}



/* Function:  esl_sqio_ReadWindow()
 * Synopsis:  Read next window of sequence.
 *
 * Purpose:   Read a next window of <W> residues from open file <sqfp>,
 *            keeping <C> residues from the previous window as
 *            context, and keeping previous annotation in the <sq>
 *            as before. 
 *            
 *            If this is the first window of a new sequence record,
 *            <C> is ignored (there's no previous context yet), and
 *            the annotation fields of the <sq> (name, accession, and
 *            description) are initialized by reading the sequence
 *            record's header. This is the only time the annotation
 *            fields are initialized.
 *            
 *            On return, <sq->dsq[]> contains the window and its
 *            context; residues <1..sq->C> are the previous context,
 *            and residues <sq->C+1..sq->n> are the new window.  The
 *            start and end coordinates of the whole <dsq[1..n]>
 *            (including context) in the original source sequence are
 *            <sq->start..sq->end>. (Or, for text mode sequences,
 *            <sq->seq[0..sq->C-1,sq->C..sq->n-1]>, while <start> and
 *            <end> coords are still <1..L>.)
 *
 *            When a sequence record is completed and no more data
 *            remain, <eslEOD> is returned, with an ``info'' <sq>
 *            structure (containing the annotation and the total
 *            sequence length <L>, but no sequence). (The total
 *            sequence length <L> is unknown in <sq> until this
 *            <eslEOD> return.)
 *            
 *            The caller may then do one of two things before calling
 *            <esl_sq_ReadWindow()> again; it can reset the sequence
 *            with <esl_sq_Reuse()> to continue reading the next
 *            sequence in the file, or it can set a negative <W> as a
 *            signal to read windows from the reverse complement
 *            (Crick) strand. Reverse complement reading only works
 *            for nucleic acid sequence. 
 *            
 *            If you read the reverse complement strand, you must read
 *            the whole thing, calling <esl_sqio_ReadWindow()> with
 *            negative <W> windows until <eslEOD> is returned again
 *            with an empty (info-only) <sq> structure. When that
 *            <EOD> is reached, the <sqfp> is repositioned at the
 *            start of the next sequence record; the caller should now
 *            <Reuse()> the <sq>, and the next <esl_sqio_ReadWindow()>
 *            call must have a positive <W>, corresponding to starting
 *            to read the Watson strand of the next sequence.
 *
 *            Note that the <ReadWindow()> interface is designed for
 *            an idiom of sequential reading of complete sequences in
 *            overlapping windows, possibly on both strands; if you
 *            want more freedom to move around in the sequence
 *            grabbing windows in another order, you can use the
 *            <FetchSubseq()> interface.
 *
 *            Reading the reverse complement strand requires file
 *            repositioning, so it will not work on non-repositionable
 *            streams like gzipped files or a stdin pipe. Moreover,
 *            for reverse complement input to be efficient, the
 *            sequence file should have consistent line lengths, 
 *            suitable for SSI's fast subsequence indexing.
 *            
 * Returns:   <eslOK> on success; <sq> now contains next window of
 *            sequence, with at least 1 new residue. The number
 *            of new residues is <sq->W>; <sq->C> residues are 
 *            saved from the previous window. Caller may now
 *            process residues <sq->dsq[sq->C+1]..sq->dsq[sq->n]>.
 *            
 *            <eslEOD> if no new residues were read for this sequence
 *            and strand, and <sq> now contains an empty info-only
 *            structure (annotation and <L> are valid). Before calling
 *            <esl_sqio_ReadWindow()> again, caller will either want
 *            to make <W> negative (to start reading the Crick strand
 *            of the current sequence), or it will want to reset the
 *            <sq> (with <esl_sq_Reuse()>) to go on the next sequence.
 *            
 *            <eslEOF> if we've already returned <eslEOD> before to
 *            signal the end of the previous seq record, and moreover,
 *            there's no more sequence records in the file.
 *            
 *            <eslEINVAL> if an invalid residue is found in the
 *            sequence, or if you attempt to take the reverse
 *            complement of a sequence that can't be reverse
 *            complemented.
 *
 * Throws:    <eslESYNTAX> if you try to read a reverse window before
 *            you've read forward strand.
 *            
 *            <eslECORRUPT> if something goes awry internally in the
 *            coordinate system.
 *            
 *            <eslEMEM> on allocation error.
 */
int
esl_sqio_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq)
{
  ESL_SQ *tmpsq = NULL;
  int     actual_start;
  int64_t nres;
  int64_t line;
  off_t   offset;
  int     status;

  if (esl_sqio_IsAlignment(sqfp->format))
    {
      // special: if we're initializing a revcomp window read, back sqfp->idx up one
      if (W < 0 && sq->start == 0) sqfp->idx--;

      if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
        { // need new alignment? 
          esl_msa_Destroy(sqfp->msa);
          status = esl_msafile_Read(sqfp->afp, &(sqfp->msa));
          if (status == eslEFORMAT)
            { // oops, a parse error; upload the error info from afp to sqfp 
              sqfp->linenumber = sqfp->afp->linenumber;
              strcpy(sqfp->errbuf, sqfp->afp->errmsg); // errbufs must be same size
              return eslEFORMAT;
            }
          else if (status != eslOK) goto ERROR;
          sqfp->idx = 0;
        }
      
      /* grab appropriate seq from alignment into tmpsq */
      if ((status = esl_sq_FetchFromMSA(sqfp->msa, sqfp->idx, &tmpsq)) != eslOK) goto ERROR;

      /* by default, tmpsq is text-mode, convert it to digital if that's what sq is*/
      if (sq->seq == NULL &&
          (( status = esl_sq_Digitize(sq->abc, tmpsq)) != eslOK)) 
        goto ERROR;

    /* Figure out tmpsq coords we'll put in sq */
    if (W > 0)
      { /* forward strand */
        sq->C     = ESL_MIN(sq->n, C);
        sq->start = sq->end - sq->C + 1;
        sq->end   = ESL_MIN(tmpsq->L, sq->end + W);
        sq->n     = sq->end - sq->start + 1;
        sq->W     = sq->n - sq->C;
      }
    else
      { /* reverse strand */
        if (sq->L == -1) ESL_XEXCEPTION(eslESYNTAX, "Can't read reverse complement until you've read forward strand");

        sq->C     = ESL_MIN(sq->n, sq->end + C - 1);
        sq->end   = (sq->start == 0 ? sq->L : sq->end + sq->C - 1);
        sq->start = ESL_MAX(1, sq->end + W - sq->C - 1);
        sq->n     = sq->end - sq->start + 1;
        sq->W     = sq->n - sq->C;
      }

    if (sq->W == 0) /* no new sequence? that's the EOD case */
      {
        sq->start = 0;
        sq->end   = 0;
        sq->C     = 0;
        sq->W     = 0;
        sq->n     = 0;
        sq->L     = tmpsq->L;
        if      (sq->dsq) sq->dsq[1] = eslDSQ_SENTINEL;
        else if (sq->seq) sq->seq[0] = '\0';

        sqfp->idx++;
        esl_sq_Destroy(tmpsq);
        return eslEOD;
      }

    /* Copy the sequence frag.  */
    if (tmpsq->ss != NULL && sq->ss == NULL) ESL_ALLOC(sq->ss, sizeof(char) * (sq->salloc)); /* this *must* be for salloc  */
    esl_sq_GrowTo(sq, sq->n);
    if (tmpsq->seq != NULL)
      { /* text mode */
        memcpy(sq->seq, tmpsq->seq + sq->start - 1, sizeof(char) * sq->n);
        sq->seq[sq->n] = '\0';
        if (tmpsq->ss != NULL) {
          memcpy(sq->ss, tmpsq->ss + sq->start - 1, sizeof(char) * sq->n);
          sq->ss[sq->n] = '\0';
        }
      }
    else
      {
        memcpy(sq->dsq + 1, tmpsq->dsq + sq->start, sizeof(ESL_DSQ) * sq->n);
        sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
        if (tmpsq->ss != NULL) {
          memcpy(sq->ss + 1, tmpsq->ss + sq->start, sizeof(char) * sq->n);
          sq->ss[sq->n+1] = '\0';
        }
      }
    if (W < 0 && esl_sq_ReverseComplement(sq) != eslOK)
      ESL_XFAIL(eslEINVAL, sqfp->errbuf, "Can't reverse complement that sequence window");

    /* Copy annotation */
    if ((status = esl_sq_SetName     (sq, tmpsq->name))   != eslOK) goto ERROR;
    if ((status = esl_sq_SetSource   (sq, tmpsq->name))   != eslOK) goto ERROR;
    if ((status = esl_sq_SetAccession(sq, tmpsq->acc))    != eslOK) goto ERROR;
    if ((status = esl_sq_SetDesc     (sq, tmpsq->desc))   != eslOK) goto ERROR;
    sq->roff = -1;
    sq->doff = -1;
    sq->eoff = -1;
    sq->hoff = -1;

    esl_sq_Destroy(tmpsq);
    return eslOK;
  }

  /* Now for the normal case: we're reading a normal unaligned seq file, not an alignment. */
  /* Negative W indicates reverse complement direction */
  if (W < 0)
  {
    if (sq->L == -1) ESL_EXCEPTION(eslESYNTAX, "Can't read reverse complement until you've read forward strand");

    if (sq->end == 1 || sq->L == 0)
      { /* last end == 1 means last window was the final one on reverse strand,
         * so we're EOD; jump back to last forward position.
         *
         * Also check for the unusual case of sq->L == 0, a completely empty sequence:
         * in that case, immediately return eslEOD.
         */
        if (sqfp->bookmark_offset > 0) 
          {
            if (esl_sqfile_Position(sqfp, sqfp->bookmark_offset) != eslOK)
              ESL_EXCEPTION(eslECORRUPT, "Failed to reposition seq file at last forward bookmark");
            sqfp->linenumber = sqfp->bookmark_linenum;
          } 
        else 
          sqfp->nc = 0; /* signals EOF */

        sqfp->bookmark_offset  = 0;
        sqfp->bookmark_linenum = 0;

        sq->start      = 0;
        sq->end        = 0;
        sq->C          = 0;
        sq->W          = 0;
        sq->n          = 0;
        /* sq->L stays as it is */
        if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL;
        else                 sq->seq[0] = '\0';
        return eslEOD;
      }

    /* If sq->start == 0, we haven't read any reverse windows yet;
     * init reading from sq->L
     */
    W = -W;
    if (sq->start == 0)
      {
        sq->start        = ESL_MAX(1, (sq->L - W + 1));
        sq->end          = sq->L;
        sq->C            = 0;
        sq->W            = sq->end - sq->start + 1;
        sqfp->curbpl     = -1;
        sqfp->currpl     = -1;
        sqfp->prvbpl     = -1;
        sqfp->prvrpl     = -1;
        sqfp->linenumber = -1;
        sqfp->L          = -1;
      }
    else
      { /* Else, we're continuing to next window; prv was <end>..<start> */
        sq->C     = ESL_MIN(C, sq->L - sq->end + 1);    /* based on prev window's end */
        sq->end   = sq->end + sq->C - 1;                /* also based on prev end     */
        sq->start = ESL_MAX(1, (sq->end - W - sq->C + 1));
        sq->W     = sq->end - sq->start + 1 - sq->C;
      }

    /* Now position for a subseq fetch of <start..end> on fwd strand, using SSI offset calc  */
    ESL_DASSERT1(( sq->doff != 0 ));
    if (sqfp->bpl <= 0 || sqfp->rpl <= 0)      /* no help; brute force resolution. */
      {
        offset       = sq->doff;
        actual_start = 1;
      }
    else if (sqfp->bpl == sqfp->rpl+1)         /* residue resolution */
      {
        line = (sq->start-1) / sqfp->rpl; /* data line #0.. that <end> is on */
        offset       = sq->doff + line * sqfp->bpl + (sq->start-1)%sqfp->rpl;
        actual_start = sq->start;
      }
    else /* line resolution */
      {
        line         = (sq->start-1) / sqfp->rpl; /* data line #0.. that <end> is on */
        offset       = sq->doff + line * sqfp->bpl;
        actual_start = 1 + line * sqfp->rpl;
      }

    if (esl_sqfile_Position(sqfp, offset) != eslOK)
      ESL_EXCEPTION(eslECORRUPT, "Failed to reposition seq file for reverse window read");

    /* grab the subseq and rev comp it */
    if ((status = esl_sq_GrowTo(sq, sq->C+sq->W)) != eslOK) return status;
    sq->n = 0;
    status = read_nres(sqfp, sq, (sq->start - actual_start), (sq->end - sq->start + 1), &nres);

    if (status != eslOK || nres < (sq->end - sq->start + 1))
      ESL_EXCEPTION(eslECORRUPT, "Failed to extract %d..%d", sq->start, sq->end);

    status = esl_sq_ReverseComplement(sq);
    if      (status    == eslEINVAL) ESL_FAIL(eslEINVAL, sqfp->errbuf, "can't reverse complement that seq - it's not DNA/RNA");
    else if (status    != eslOK)     return status;

    return eslOK;
  }

  /* Else, we're reading the forward strand */
  else
  { /* sq->start == 0 means we haven't read any windows on this sequence yet...
     * it's a new record, and we need to initialize with the header and
     * the first window. This is the only case that we're allowed to return
     * EOF from.
     */
    if (sq->start == 0)
      {
        if (sqfp->nc == 0) return eslEOF;
        if ((status = sqfp->parse_header(sqfp, sq)) != eslOK) return status; /* EOF, EFORMAT */
        sq->start     = 1;
        sq->C         = 0;  // no context in first window
        sq->L         = -1; // won't be known 'til EOD.
        sqfp->L       = 0;  // init to 0, so we can count residues as we go
        esl_sq_SetSource(sq, sq->name);
        /* the <sqfp->buf> is now positioned at the start of seq data */
        /* sqfp->linenumber is ok where it is */
        /* the header_*() routines initialized rpl,bpl bookkeeping at start of seq line,
         * and also sq->doff,roff.
         */
      }
    else
      { /* else we're reading a window other than first; slide context over. */
        sq->C = ESL_MIN(C, sq->n);

        /* if the case where the window is smaller than the context and the
         * context is not full, it is not necessary to move the context part
         * of the sequence that has been read in.
         */
        if (sq->C >= C) {
          /* now handle the case where the context is full */
          if (sq->seq != NULL) memmove(sq->seq,   sq->seq + sq->n - sq->C,     sq->C);
          else                 memmove(sq->dsq+1, sq->dsq + sq->n - sq->C + 1, sq->C);
          sq->start = sqfp->L - sq->C + 1;
          sq->n = C;
        }
      }

    if ((status = esl_sq_GrowTo(sq, C+W)) != eslOK) return status; // EMEM
    status = read_nres(sqfp, sq, 0, W, &nres);
    sqfp->L += nres;

    if (status == eslEOD)
      { // Forward strand is done. 0 residues were read. Return eslEOD and an empty (info) <sq>.
        if ((status = sqfp->parse_end(sqfp, sq)) != eslOK) return status;
        sq->start      = 0;
        sq->end        = 0;
        sq->C          = 0;
        sq->W          = 0;
        sq->L          = sqfp->L;
        sq->n          = 0;

        if (sqfp->nc > 0) 
          sqfp->bookmark_offset  = sqfp->boff+sqfp->bpos; // remember where the next seq starts
        else {
          sqfp->bookmark_offset  = 0;                     // signals for EOF, no more seqs
          sqfp->bookmark_linenum = 0;
        }

        if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL; // erase the saved context 
        else                 sq->seq[0] = '\0';
        return eslEOD;
      }
    else if (status == eslOK)
      { // Forward strand is still in progress. <= W residues were read. Return eslOK.
        sq->end        = sq->start + sq->C + nres - 1;
        sq->W          = nres;
        return eslOK;
      }
    else return status; // EFORMAT,EMEM
  }
  /*NOTREACHED*/
  return eslOK;

 ERROR:
  if (tmpsq) esl_sq_Destroy(tmpsq);
  return status;
}

/* Function:  esl_sqio_ReadBlock()
 * Synopsis:  Read the next block of sequences from a file.
 *
 * Purpose:   Reads a block of sequences from open sequence file <sqfp> into 
 *            <sqBlock>.
 *            
 *            If <long_target> is FALSE, the sequences are expected to
 *            be protein: individual sequences won't be long, so we
 *            can read complete sequences. If <max_sequences> is > 0,
 *            read <max_sequences> sequences, up to at most
 *            eslSQFILE_BLOCK_MAXRES residues. <max_init_window> is
 *            irrelevant if <long_target> is false.
 *
 *            If <long_target> is TRUE, the sequences are expected to
 *            be DNA.  Because individual DNA sequences can exceed
 *            eslSQFILE_BLOCK_MAXRES, we use ReadWindow() to read
 *            chunks of sequence no larger than <max_residues>, and
 *            must allow for the possibility that a request will be
 *            made to continue reading a partly-read sequence. This
 *            case also respects the <max_sequences> limit.
 * 
 *            If <long_target> is true and <max_init_window> is TRUE,
 *            the first window read from each sequence (of length L)
 *            is always min(L, <max_residues>). If <max_init_window>
 *            is FALSE, then the length of the first window read from
 *            each sequence is calculated differently as
 *            max(<max_residues> - <size>, <max_residues> * .05);
 *            where <size> is total number of residues already
 *            existing in the block. <max_init_window> == TRUE mode
 *            was added to ensure that the window boundaries read are
 *            not dependent on the order of the sequence in the file,
 *            thus ensuring reproducibility if (for example) a user
 *            extracts one sequence from a file and reruns a program
 *            on it (and all else remains equal).
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sqBlock>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if a parse error is encountered.
 *            Call <esl_sqfile_GetErrorBuf()> to get a ptr to a
 *            user-directed error message describing the problem,
 *            including the line number on which it was found.
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_sqio_ReadBlock(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target)
{
  ESL_SQ *tmpsq  = NULL;
  int     i      = 0;
  int     size   = 0;
  int     status = eslOK;

  sqBlock->count = 0;
  if (max_sequences < 1 || max_sequences > sqBlock->listSize)
    max_sequences = sqBlock->listSize;

  if (!long_target)
    {  /* in these cases, an individual sequence won't ever be really long,
          so just read in a sequence at a time  */
      for (i = 0; i < max_sequences && size < eslSQFILE_BLOCK_MAXRES; ++i)
        {
          status = esl_sqio_Read(sqfp, sqBlock->list + i);
          if (status != eslOK) break;
          size += sqBlock->list[i].n;
          ++sqBlock->count;
        }
    }
  else
    { /* DNA, not an alignment.  Might be really long sequences */

      if (max_residues < 1)
        max_residues = eslSQFILE_BLOCK_MAXRES;

      tmpsq = esl_sq_CreateDigital(sqBlock->list->abc);
      // if complete flag is set to FALSE, then the prior block must have ended with a window that was a possibly
      // incomplete part of it's full sequence. Read another overlapping window.
      if (! sqBlock->complete )
        {
          //overloading C as indicator of how big C should be for this window reading action
          status = esl_sqio_ReadWindow(sqfp, sqBlock->list->C, max_residues, sqBlock->list);
          if (status == eslOK)
            {
              sqBlock->count = i = 1;
              size = sqBlock->list->n - sqBlock->list->C;
              sqBlock->list->L = sqfp->L;
              if (size == max_residues)
                { // Filled the block with a single very long window.
                  sqBlock->complete = FALSE; // default value, unless overridden below
                  status = skip_whitespace(sqfp);
                  if ( status != eslOK ) { // either EOD or end of buffer (EOF) was reached before the next character was seen
                    sqBlock->complete = TRUE;
                    status = eslOK;
                  }
                  if (tmpsq) esl_sq_Destroy(tmpsq);
                  return status;
                }
              else
                {
                  // Burn off EOD (see notes for similar entry ~25 lines below), then go fetch the next sequence
                  esl_sq_Reuse(tmpsq);
                  tmpsq->start =  sqBlock->list->start ;
                  tmpsq->C = 0;
                  status = esl_sqio_ReadWindow(sqfp, 0, max_residues, tmpsq);
                  if (status != eslEOD) {
                    if (tmpsq != NULL) esl_sq_Destroy(tmpsq);
                    return status; //surprising
                  }
                  //sqBlock->list->L = tmpsq->L;
                }
            }
          else if (status == eslEOD)
            { // turns out there isn't any more of the sequence to read, after all
            }
          else
            {
              if(tmpsq) esl_sq_Destroy(tmpsq);
              return status;
            }
        } // otherwise, just start at the beginning


      for (  ; i < max_sequences && size < max_residues; ++i)
        {
          /* restricted request_size is used to ensure that all blocks are pretty close to the
           * same size. Without it, we may either naively keep asking for max_residue windows,
           * which can result in a window with ~2*max_residues ... or we can end up with absurdly
           * short fragments at the end of blocks
           */
          int request_size = (max_init_window) ? max_residues : ESL_MAX(max_residues-size, max_residues * .05);

          esl_sq_Reuse(tmpsq);
          esl_sq_Reuse(sqBlock->list + i);

          status = esl_sqio_ReadWindow(sqfp, 0, request_size, tmpsq); 
          esl_sq_Copy(tmpsq, sqBlock->list +i);
          if (status != eslOK && status != eslEOD) break; /* end of sequences (eslEOF), or we read an empty seq (eslEOD) or error (other)  */

          size += sqBlock->list[i].n - sqBlock->list[i].C;
          sqBlock->list[i].L = sqfp->L;
          ++(sqBlock->count);

          if (size >= max_residues)
            {
              // a full window worth of sequence has been read; did we reach the end of the final sequence in the block?
              sqBlock->complete = FALSE; // default value, unless overridden below
              status = skip_whitespace(sqfp);
              if ( status != eslOK ) { // either EOD or end of buffer (EOF) was reached before the next character was seen
                sqBlock->complete = TRUE;
                status = eslOK;
              }

              if(tmpsq != NULL) esl_sq_Destroy(tmpsq);
              return status;
            }
          else if (status == eslEOD)
            {
              /* We've read an empty sequence of length 0, rare, but
               * possible, and we need to be able to handle it
               * gracefully. Ensure L is 0, set status to eslOK and move
               * on, we've already incremented sqBlock->count by 1
               * above. This means our block may contain zero-length
               * sequences when we return (that is, we still add these
               * seqs onto the block instead of skipping them altogether).
               */
              sqBlock->list[i].L = 0; /* actually, this should already be 0... */
              status = eslOK;
            }
          else {
            /* Sequence finished, but haven't yet reached max_residues. Need to burn off the EOD value
             * that will be returned by the next ReadWindow call. Can just use a tmp sq, after setting
             *  a couple values ReadWindow needs to see for correct processing.
             */
            esl_sq_Reuse(tmpsq);
            tmpsq->start =  sqBlock->list[i].start ;
            tmpsq->C = 0;
            status = esl_sqio_ReadWindow(sqfp, 0, max_residues, tmpsq);

            if (status != eslEOD) {
              if (tmpsq) esl_sq_Destroy(tmpsq);
              return status; //surprising
            }
            //sqBlock->list[i].L = tmpsq->L;
            status = eslOK;
          }
        }
    }
  
  /* EOF will be returned only in the case were no sequences were read */
  if (status == eslEOF && i > 0) status = eslOK;
  sqBlock->complete = TRUE;
  if (tmpsq) esl_sq_Destroy(tmpsq);
  return status;
}

/* Function:  esl_sqio_Parse()
 * Synopsis:  Parse a sequence already read into a buffer.
 *
 * Purpose:   Parse the buffer <buf> for a sequence <s> of type
 *            <format>.  The buffer must contain the entire sequence.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM>  on allocation error.
 *            <eslEFORMAT>  on parsing error.
 *            <eslEINVAL> on unsupported format.
 */
int
esl_sqio_Parse(char *buf, int size, ESL_SQ *sq, int format)
{
  ESL_SQFILE sqfp;
  int64_t    epos;
  int64_t    n;
  int        status;

  /* fill in a dummy esl_sqfile structure used to parse buf */
  sqfile_init(&sqfp);
  sqfp.mem          = buf;
  sqfp.mn           = size;
  sqfp.do_buffer    = TRUE;

  /* Configure the <sqfp>'s parser and inmaps for this format. */
  switch (format) {
  case eslSQFILE_EMBL:     
  case eslSQFILE_UNIPROT:  
    config_embl(&sqfp);    
    inmap_embl(&sqfp, NULL);
    break;
  case eslSQFILE_GENBANK:  
  case eslSQFILE_DDBJ:     
    config_genbank(&sqfp); 
    inmap_genbank(&sqfp, NULL);
    break;
  case eslSQFILE_FASTA:    
    config_fasta(&sqfp);   
    inmap_fasta(&sqfp, NULL);
    break;
  case eslSQFILE_DAEMON:    
    config_daemon(&sqfp);   
    inmap_daemon(&sqfp, NULL);
    break;
  default:
    ESL_EXCEPTION(eslEINVAL, "can't parse that format");
  }

  /* Main case: read next seq from sqfp's stream */
  if ((status = sqfp.parse_header(&sqfp, sq)) != eslOK) goto ERROR; /* EOF, EFORMAT */

  do {
    if ((status = seebuf(&sqfp, -1, &n, &epos)) == eslEFORMAT)      goto ERROR;     // EOD is ok! Don't overwrite <status> now, until we reach the status==eslEOD check below. 
    if ( esl_sq_GrowTo(sq, sq->n + n) != eslOK) { status = eslEMEM; goto ERROR; }   // (though overwriting <status> as we fail is ok)
    addbuf(&sqfp, sq, n);
    sqfp.L   += n;
    sq->eoff   = sqfp.boff + epos - 1;
    if (status == eslEOD)     break;
  } while ((status = loadbuf(&sqfp)) == eslOK);
    
  if      (status == eslEOF)
    {
      if (! sqfp.eof_is_ok) ESL_XFAIL(eslEFORMAT, sqfp.errbuf, "Unexpected EOF; file truncated?"); 
      if ((status = sqfp.parse_end(&sqfp, sq)) != eslOK) goto ERROR;
    }
  else if (status == eslEOD)
    {
      sqfp.bpos = epos;
      if ((status = sqfp.parse_end(&sqfp, sq)) != eslOK) goto ERROR;
    }
  else if (status != eslOK) goto ERROR;

  if (sq->dsq != NULL) sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else                 sq->seq[sq->n] = '\0';
  sq->start = 1;
  sq->end   = sq->n;
  sq->C     = 0;
  sq->W     = sq->n;
  sq->L     = sq->n;

  if (sqfp.balloc > 0) free(sqfp.buf);
  return eslOK;

 ERROR:
  if (sqfp.balloc > 0) free(sqfp.buf); // clang static analyzer can report a false positive leak here; 
  return status;                           // it can't see that we only allocated sqfp.buf if balloc > 0
}
/*-------------- end, sequence reading (sequential) -------------*/



/*****************************************************************
 * 4. Sequence writing
 *****************************************************************/

/* Function:  esl_sqio_Write()
 * Synopsis:  Write a sequence to a file.
 *
 * Purpose:   Write sequence <s> to an open FILE <fp> in file format
 *            <format>.
 * 
 *            If <update> is <TRUE>, set the offsets for sequence <s>.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEWRITE> on system write error, such as filled disk.
 */
int
esl_sqio_Write(FILE *fp, ESL_SQ *s, int format, int update)
{
  ESL_MSA *msa;
  int status;

  if (esl_sqio_IsAlignment(format))
    {
      if ((status = convert_sq_to_msa(s, &msa)) != eslOK) return status;
      status = esl_msafile_Write(fp, msa, format);
      esl_msa_Destroy(msa);
      return status;
    }

  switch (format) {
  case eslSQFILE_FASTA:   
  case eslSQFILE_HMMPGMD:
    status = esl_sqio_WriteFasta(fp, s, update); 
    break;
  default: 
    ESL_EXCEPTION(eslEINCONCEIVABLE, "can't write that format");
  }
  return status;
}

/* Function:  esl_sqio_Echo()
 * Synopsis:  Echo a sequence's record onto output stream.
 *
 * Purpose:   Given a complete <sq> that we have read by some means
 *            from an open <sqfp>; echo that sequence's record
 *            onto the output stream <ofp>. 
 *
 *            This allows records to be regurgitated exactly as they
 *            appear, rather than writing the subset of information
 *            stored in an <ESL_SQ>. <easel sfetch> in the miniapps uses
 *            this, for example.
 *            
 *            Because this relies on repositioning the <sqfp>, it
 *            cannot be called on non-positionable streams (stdin or
 *            gzipped files). Because it relies on the sequence lying
 *            in a contiguous sequence of bytes in the file, it cannot
 *            be called on a sequence in a multiple alignment file.
 *            Trying to do so throws an <eslEINVAL> exception.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL>   if <sqfp> isn't a repositionable sequence file.
 *            <eslECORRUPT> if we run out of data, probably from bad offsets
 *            <eslEMEM>     on allocation failure.
 *            <eslESYS>     on system call failures.
 */
int
esl_sqio_Echo(ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)
{
  int64_t save_linenumber;
  int     save_currpl;
  int     save_curbpl;
  int     save_prvrpl;
  int     save_prvbpl;
  int64_t save_L;
  int     n;
  int     nwritten;
  int     status;

  if (sqfp->do_stdin)                     ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from standard input");
  if (sqfp->do_gzip)                      ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from a gzipped file");
  if (esl_sqio_IsAlignment(sqfp->format)) ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from an alignment file");
  if (sq->roff == -1 || sq->eoff == -1)   ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence without disk offset info");

  save_linenumber = sqfp->linenumber;
  save_currpl     = sqfp->currpl;
  save_curbpl     = sqfp->curbpl;
  save_prvrpl     = sqfp->prvrpl;
  save_prvbpl     = sqfp->prvbpl;
  save_L          = sqfp->L;

  status = esl_sqfile_Position(sqfp, sq->roff);
  if      (status == eslEOF) ESL_EXCEPTION(eslECORRUPT, "repositioning failed; bad offset?");
  else if (status != eslOK)  return status;

  while (sqfp->boff + sqfp->nc <= sq->eoff)
    {
      if ((int) fwrite(sqfp->buf, sizeof(char), sqfp->nc, ofp) != sqfp->nc) ESL_EXCEPTION(eslESYS, "fwrite() failed");
      if (loadbuf(sqfp) != eslOK)  ESL_EXCEPTION(eslECORRUPT, "repositioning failed; bad offset?");
    } 
  n =  sq->eoff - sqfp->boff + 1;
  nwritten = fwrite(sqfp->buf, sizeof(char), n, ofp);
  if (nwritten != n) ESL_EXCEPTION(eslESYS, "fwrite() failed");

  status = esl_sqfile_Position(sqfp, sq->roff);
  if      (status == eslEOF) ESL_EXCEPTION(eslECORRUPT, "repositioning failed; bad offset?");
  else if (status != eslOK)  return status;

  sqfp->linenumber = save_linenumber;
  sqfp->currpl     = save_currpl;
  sqfp->curbpl     = save_curbpl;
  sqfp->prvrpl     = save_prvrpl;
  sqfp->prvbpl     = save_prvbpl;
  sqfp->L          = save_L;
  return eslOK;  
}



/* Function:  esl_sqio_WriteFasta()
 * Synopsis:  Write a sequence in FASTA format
 *
 * Purpose:   Write sequence <sq> in FASTA format to the open stream <fp>.
 * 
 *            If <save_offsets> is TRUE, then store record, data, and end 
 *            offsets in <sq>; this ability is used by unit tests.
 *
 * Returns:   <eslOK> on success.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEWRITE> on system write error.
 */
int
esl_sqio_WriteFasta(FILE *fp, ESL_SQ *sq, int save_offsets)
{
  char     buf[61];
  int64_t  pos;

  if (save_offsets) sq->roff = ftello(fp);
  if (fprintf(fp, ">%s", sq->name)                     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "fasta seq write failed");
  if (sq->acc[0]  != 0 && fprintf(fp, " %s", sq->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "fasta seq write failed");
  if (sq->desc[0] != 0 && fprintf(fp, " %s", sq->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "fasta seq write failed");
  if (save_offsets) sq->hoff = ftello(fp);
  if (fputc('\n', fp)                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "fasta seq write failed");

  buf[60] = '\0';
  if (save_offsets) sq->doff = ftello(fp);
  for (pos = 0; pos < sq->n; pos += 60)
  {
      if (sq->dsq != NULL) esl_dsq_TextizeN(sq->abc, sq->dsq+pos+1, 60, buf);
      else                 strncpy(buf, sq->seq+pos, 60);
      if (fprintf(fp, "%s\n", buf) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "fasta seq write failed");
  }
  if (save_offsets) sq->eoff = ftello(fp) - 1;
  return eslOK;
}
/*------------------ end, sequence writing  ---------------------*/



/*****************************************************************
 * 5. Miscellaneous functions
 *****************************************************************/ 

/* Function:  esl_sqfile_GetErrorBuf()
 * Synopsis:  Return ptr to the error buffer
 *
 * Purpose:   Returns the pointer to the error buffer.
 *            Each parser is responsible for formatting
 *            a zero terminated string describing the
 *            error condition.
 *
 * Returns:   A pointer to the error message.
 */
const char *
esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp)
{
  return sqfp->errbuf;
}


/* Function:  esl_sqfile_IsRewindable()
 * Synopsis:  Return <TRUE> if <sqfp> can be rewound.
 *
 * Purpose:   Returns <TRUE> if <sqfp> can be rewound (positioned 
 *            to an offset of zero), in order to read it a second
 *            time.
 */
int
esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp)
{
  if (sqfp->do_gzip || sqfp->do_stdin) return FALSE;
  else                                 return TRUE;
}

/* Function:  esl_sqio_IsAlignment()
 * Synopsis:  Return TRUE for alignment file format codes.
 *
 * Purpose:   Returns TRUE if <fmt> is an alignment file
 *            format code; else returns FALSE.
 *            
 *            This function only checks the convention
 *            that <fmt> codes $<$100 are unaligned formats,
 *            and $\geq$100 are aligned formats. It does
 *            not check that <fmt> is a recognized format
 *            code.
 */
int
esl_sqio_IsAlignment(int fmt)
{
  return (fmt >= 100 ? TRUE : FALSE);
}


/* Function:  esl_sqio_EncodeFormat()
 * Synopsis:  Convert a string to an internal format code.
 *
 * Purpose:   Given <fmtstring>, return format code.  For example, if
 *            <fmtstring> is "fasta", returns <eslSQFILE_FASTA>. Returns 
 *            <eslSQFILE_UNKNOWN> if <fmtstring> doesn't exactly match a 
 *            known format.
 *            
 *            Matching is case insensitive; fasta, FASTA, and FastA
 *            all return <eslSQFILE_FASTA>, for example.
 */
int
esl_sqio_EncodeFormat(char *fmtstring)
{
  if (strcasecmp(fmtstring, "fasta")     == 0) return eslSQFILE_FASTA;
  if (strcasecmp(fmtstring, "embl")      == 0) return eslSQFILE_EMBL;
  if (strcasecmp(fmtstring, "genbank")   == 0) return eslSQFILE_GENBANK;
  if (strcasecmp(fmtstring, "ddbj")      == 0) return eslSQFILE_DDBJ;
  if (strcasecmp(fmtstring, "uniprot")   == 0) return eslSQFILE_UNIPROT;
  if (strcasecmp(fmtstring, "daemon")    == 0) return eslSQFILE_DAEMON;
  if (strcasecmp(fmtstring, "hmmpgmd")   == 0) return eslSQFILE_HMMPGMD;
  if (strcasecmp(fmtstring, "fmindex")   == 0) return eslSQFILE_FMINDEX;
  return esl_msafile_EncodeFormat(fmtstring);
}

/* Function:  esl_sqio_DecodeFormat()
 * Synopsis:  Returns descriptive string for file format code.
 *
 * Purpose:   Given a format code <fmt>, returns a string label for
 *            that format. For example, if <fmt> is <eslSQFILE_FASTA>,
 *            returns "FASTA". 
 */
char *
esl_sqio_DecodeFormat(int fmt)
{
  if (esl_sqio_IsAlignment(fmt)) return esl_msafile_DecodeFormat(fmt);

  switch (fmt) {
  case eslSQFILE_UNKNOWN:    return "unknown";
  case eslSQFILE_FASTA:      return "FASTA";
  case eslSQFILE_EMBL:       return "EMBL";
  case eslSQFILE_GENBANK:    return "GenBank";
  case eslSQFILE_DDBJ:       return "DDBJ";
  case eslSQFILE_UNIPROT:    return "UniProt";
  case eslSQFILE_DAEMON:     return "daemon";
  case eslSQFILE_HMMPGMD:    return "hmmpgmd";
  case eslSQFILE_FMINDEX:    return "fmindex";
  default:                   break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__,  "no such sqio format code %d", fmt);
  return NULL;
}


/* Function:  esl_sqfile_Position()
 * Synopsis:  Reposition an open sequence file to an offset.
 *
 * Purpose:   Reposition an open <sqfp> to offset <offset>.
 *            <offset> would usually be the first byte of a
 *            desired sequence record.
 *            
 *            Only normal sequence files can be positioned to a
 *            nonzero offset. If <sqfp> corresponds to a standard
 *            input stream or gzip -dc stream, it may not be
 *            repositioned. If <sqfp> corresponds to a multiple
 *            sequence alignment file, the only legal <offset>
 *            is 0, to rewind the file to the beginning and 
 *            be able to read the entire thing again.
 *            
 *            After <esl_sqfile_Position()> is called on a nonzero
 *            <offset>, <sqfp->linenumber> and other bookkeeping
 *            information is unknown.  If caller knows it, it should
 *            set it explicitly.
 *            
 *            See the SSI module for manipulating offsets and indices.
 *
 * Returns:   <eslOK>     on success;
 *            <eslEOF>    if no data can be read from this position.
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails.
 *            <eslEMEM> on (re-)allocation failure.
 *            <eslEINVAL> if the <sqfp> is not positionable, or <offset> is invalid.
 *            <eslENOTFOUND> if in trying to rewind an alignment file  
 *              by closing and reopening it, the open fails.
 *
 *            On errors, the state of <sqfp> is indeterminate, and
 *            it should not be used again.
 */
int
esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset)
{
  int status;

  if (sqfp->do_stdin)                  ESL_EXCEPTION(eslEINVAL, "can't Position() in standard input");
  if (sqfp->do_gzip)                   ESL_EXCEPTION(eslEINVAL, "can't Position() in a gzipped file");
  if (offset < 0)                      ESL_EXCEPTION(eslEINVAL, "bad offset");
  if (offset > 0 && sqfp->afp != NULL) ESL_EXCEPTION(eslEINVAL, "can't use esl_sqfile_Position() w/ nonzero offset on MSA file");

  if (esl_sqio_IsAlignment(sqfp->format)) 
    { /* msa file: close and reopen. maybe sometime we'll have esl_msafile_Rewind() */
      /* we have already verified that offset==0 for MSA file */
      esl_msafile_Close(sqfp->afp);
      if (sqfp->msa) esl_msa_Destroy(sqfp->msa);
      sqfp->afp = NULL;
      sqfp->msa = NULL;
      sqfp->idx = 0;
      
      /* we know we successfully opened it the first time, so a
         failure to reopen is an exception, not a user-reportable
         normal error. ENOTFOUND is the only normal error;
         EFORMAT error can't occur because we know the format and
         don't use autodetection.
       */
      status = esl_msafile_Open(NULL, sqfp->filename, NULL, sqfp->format, NULL, &(sqfp->afp));
      if      (status == eslENOTFOUND) ESL_EXCEPTION(eslENOTFOUND, "failed to reopen alignment file");
      else if (status != eslOK)        return status;
    }
  else /* normal case: unaligned sequence file */
    {
      if (fseeko(sqfp->fp, offset, SEEK_SET) != 0) ESL_EXCEPTION(eslESYS, "fseeko() failed");

      sqfp->currpl     = -1;
      sqfp->curbpl     = -1;
      sqfp->prvrpl     = -1;
      sqfp->prvbpl     = -1;
      sqfp->linenumber = (offset == 0) ? 1 : -1;  // -1 is "unknown" 
      sqfp->L          = -1;
      sqfp->mpos       = sqfp->mn;               // this forces loadbuf to load new data 
      if ((status = loadbuf(sqfp)) != eslOK) return status;
    }
  return eslOK;
}

/* Function:  esl_sqio_Ignore()
 * Synopsis:  Sets the input map to ignore one or more input characters.
 *
 * Purpose:   Set the input map of the open <sqfp> to allow
 *            the characters in the string <ignored> to appear
 *            in input sequences. These characters will be ignored.
 *
 *            For example, an application might want to ignore '*'
 *            characters in its sequence input, because some translated
 *            peptide files use '*' to indicate stop codons.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignored)
{
  int i;
  for (i = 0; ignored[i] != '\0'; i++)
    sqfp->inmap[(int) ignored[i]] = eslDSQ_IGNORED;
  return eslOK;
}

/* Function:  esl_sqio_Accept()
 * Synopsis:  Add more characters to a text-mode sqfile's inmap.
 * Incept:    SRE, Sun 01 Jun 2025 [B6 920 SAN-BOS]
 *
 * Purpose:   Add each character in the <accepted> string to the input
 *            map for <sqfp>. By default, text mode input from
 *            unaligned seqfile formats allows a-zA-Z and * as input
 *            sequence residue characters.
 *
 *            This should be called immediately after opening the
 *            <sqfp>, before the first call to esl_sqio_Read().
 *
 *            This only works for a text mode <sqfp>. You can't
 *            extend a digital alphabet that's already defined (not
 *            without making a new custom alphabet). You can, however,
 *            have a digital alphabet accept additional input chars
 *            and map them to the digital alphabet's "unknown residue"
 *            code: see <esl_sqio_AcceptAs()>.
 *
 * Args:      sqfp     - recently opened seqfile, in text mode
 *            accepted - string of chars to add to the inmap
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEINVAL> if called on a digital-mode <sqfp>. You can't
 *            extend a digital alphabet that's already defined (not
 *            without making a new custom alphabet).
 *
 * Xref:      [2025/0531-easel-reformat]
 */
int
esl_sqio_Accept(ESL_SQFILE *sqfp, const char *accepted)
{
  int i;
  
  if (sqfp->do_digital)
    ESL_EXCEPTION(eslEINVAL, "bad esl_sqio_Accept(): a digital alphabet can't accept additional chars");

  for (i = 0; accepted[i] != '\0'; i++)
    sqfp->inmap[(int) accepted[i]] = accepted[i];

  return eslOK;
}

/* Function:  esl_sqio_AcceptAs()
 * Synopsis:  Map a list of additional characters.
 *
 * Purpose:   Set the input map of the open <sqfp> to allow the 
 *            characters in the string <accepted> to appear in 
 *            input sequences. These characters will all be 
 *            mapped to the character <readas> (or, for digital
 *            sequence input, to the digitized representation 
 *            of the text character <readas> in the <sqfp>'s
 *            digital alphabet).
 *            
 *            For example, an application might want to read
 *            '*' as 'X' when reading translated peptide files
 *            that use '*' to indicate a stop codon.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_sqio_AcceptAs(ESL_SQFILE *sqfp, const char *accepted, char readas)
{
  int i;
  
  if (sqfp->do_digital) {
    for (i = 0; accepted[i] != '\0'; i++)
      sqfp->inmap[(int) accepted[i]] = esl_abc_DigitizeSymbol(sqfp->abc, readas);
  } else {
    for (i = 0; accepted[i] != '\0'; i++)
      sqfp->inmap[(int) accepted[i]] = readas;
  }
  return eslOK;
}
/*--------------- end, miscellaneous functions ------------------*/



/*****************************************************************
 * 6. Sequence/subsequence fetching, random access [with <ssi>]
 *****************************************************************/

/* Function:  esl_sqfile_OpenSSI()
 * Synopsis:  Opens an SSI index associated with a sequence file.
 *
 * Purpose:   Opens an SSI index file associated with the already open
 *            sequence file <sqfp>. If successful, the necessary
 *            information about the open SSI file is stored internally
 *            in <sqfp>.
 *            
 *            The SSI index file name is determined in one of two
 *            ways, depending on whether a non-<NULL> <ssifile_hint>
 *            is provided.
 *            
 *            If <ssifile_hint> is <NULL>, the default for
 *            constructing the SSI filename from the sequence
 *            filename, by using exactly the same path (if any) for
 *            the sequence filename, and appending the suffix <.ssi>.
 *            For example, the SSI index for <foo> is <foo.ssi>, for
 *            <./foo.fa> is <./foo.fa.ssi>, and for
 *            </my/path/to/foo.1.fa> is </my/path/to/foo.1.fa.ssi>.
 *            
 *            If <ssifile_hint> is <non-NULL>, this exact fully
 *            qualified path is used as the SSI file name.
 *
 * Returns:   <eslOK> on success, and <sqfp->ssi> is now internally
 *            valid.
 *            
 *            <eslENOTFOUND> if no SSI index file is found;
 *            <eslEFORMAT> if it's found, but appears to be in incorrect format;
 *            <eslERANGE> if the SSI file uses 64-bit offsets but we're on
 *            a system that doesn't support 64-bit file offsets.
 *
 * Throws:    <eslEINVAL> if the open sequence file <sqfp> doesn't
 *            correspond to a normal sequence flatfile -- we can't
 *            random access in .gz compressed files, standard input,
 *            or multiple alignment files that we're reading
 *            sequentially.
 *            
 *            Throws <eslEMEM> on allocation error.
 */
int
esl_sqfile_OpenSSI(ESL_SQFILE *sqfp, const char *ssifile_hint)
{
  int status;

  if (sqfp->do_gzip)     ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for a .gz compressed seq file");
  if (sqfp->do_stdin)    ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for standard input");
  if (sqfp->afp != NULL) ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for sequential input from an MSA");

  if (ssifile_hint == NULL) {
    if ((status = esl_strdup(sqfp->filename, -1, &(sqfp->ssifile)))           != eslOK) return status;
    if ((status = esl_strcat(&(sqfp->ssifile), -1, ".ssi", 4))                != eslOK) return status;
  } else {
    if ((status = esl_strdup(ssifile_hint, -1, &(sqfp->ssifile)))             != eslOK) return status;
  }

  return esl_ssi_Open(sqfp->ssifile, &(sqfp->ssi));
}



/* Function:  esl_sqfile_PositionByKey()
 * Synopsis:  Use SSI to reposition seq file to a particular sequence.
 *
 * Purpose:   Reposition <sqfp> so that the next sequence we read will
 *            be the one named (or accessioned) <key>.
 *            
 *            <sqfp->linenumber> is reset to be relative to the start
 *            of the record named <key>, rather than the start of the
 *            file.
 *
 * Returns:   <eslOK> on success, and the file <sqfp> is repositioned
 *            so that the next <esl_sqio_Read()> call will read the
 *            sequence named <key>.
 *            
 *            Returns <eslENOTFOUND> if <key> isn't found in the
 *            index; in this case, the position of <sqfp> in the file
 *            is unchanged.
 *            
 *            Returns <eslEFORMAT> if something goes wrong trying to
 *            read the index, almost certainly indicating a format
 *            problem in the SSI file.
 *            
 *            Returns <eslEOF> if, after repositioning, we fail to
 *            load the next line or buffer from the sequence file;
 *            this probably also indicates a format problem in the SSI
 *            file.
 * 
 * Throws:    <eslEMEM>   on allocation error;
 *            <eslEINVAL> if there's no open SSI index in <sqfp>;
 *            <eslESYS>   if the <fseek()> fails.
 *            
 *            In all these cases, the state of <sqfp> becomes
 *            undefined, and the caller should not use it again.
 */
int
esl_sqfile_PositionByKey(ESL_SQFILE *sqfp, const char *key)
{
  uint16_t fh;
  off_t    offset;
  int      status;

  if (sqfp->ssi == NULL) ESL_EXCEPTION(eslEINVAL,"Need an open SSI index to call esl_sqfile_PositionByKey()");
  if ((status = esl_ssi_FindName(sqfp->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  return esl_sqfile_Position(sqfp, offset);
}


/* Function:  esl_sqfile_PositionByNumber()
 * Synopsis:  Use SSI to reposition by sequence number
 *
 * Purpose:   Reposition <sqfp> so that the next sequence we 
 *            read will be the <which>'th sequence, where <which>
 *            is <0..sqfp->ssi->nprimary-1>. 
 *            
 *            <sqfp->linenumber> is reset to be relative to the start
 *            of the record named <key>, rather than the start of the
 *            file.
 *
 * Returns:   <eslOK> on success, and the file <sqfp> is repositioned.
 *            
 *            Returns <eslENOTFOUND> if there is no sequence number
 *            <which> in the index; in this case, the position of
 *            <sqfp> in the file is unchanged.
 *            
 *            Returns <eslEFORMAT> if something goes wrong trying to
 *            read the index, almost certainly indicating a format
 *            problem in the SSI file.
 *            
 *            Returns <eslEOF> if, after repositioning, we fail to
 *            load the next line or buffer from the sequence file;
 *            this probably also indicates a format problem in the SSI
 *            file.
 * 
 * Throws:    <eslEMEM>   on allocation error;
 *            <eslEINVAL> if there's no open SSI index in <sqfp>;
 *            <eslESYS>   if the <fseek()> fails.
 *            
 *            In all these cases, the state of <sqfp> becomes
 *            undefined, and the caller should not use it again.
 */
int
esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which)
{
  uint16_t fh;
  off_t    offset;
  int      status;

  if (sqfp->ssi == NULL) ESL_EXCEPTION(eslEINVAL,"Need open SSI index to call esl_sqfile_PositionByNumber()");
  if ((status = esl_ssi_FindNumber(sqfp->ssi, which, &fh, &offset, NULL, NULL, NULL)) != eslOK) return status;
  return esl_sqfile_Position(sqfp, offset);
}


/* Function:  esl_sqio_Fetch()
 * Synopsis:  Fetch a complete sequence, using SSI indexing.
 *
 * Purpose:   Fetch a sequence named (or accessioned) <key> from
 *            the repositionable, open sequence file <sqfp>.
 *            The open <sqfp> must have an open SSI index.
 *            The sequence is returned in <sq>.
 *
 * Returns:   <eslOK> on soccess.
 *            <eslEINVAL> if no SSI index is present, or if <sqfp> can't
 *            be repositioned.
 *            <eslENOTFOUND> if <source> isn't found in the file.
 *            <eslEFORMAT> if either the index file or the sequence file
 *            can't be parsed, because of unexpected format issues.
 *       
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sqio_Fetch(ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
{
  int status;

  if (sqfp->ssi == NULL) ESL_FAIL(eslEINVAL, sqfp->errbuf, "No SSI index for %s; can't fetch sequences", sqfp->filename);
  if ((status = esl_sqfile_PositionByKey(sqfp, key)) != eslOK) return status;
  if ((status = esl_sqio_Read(sqfp, sq))             != eslOK) return status;
  return eslOK;
}
  
/* Function:  esl_sqio_FetchInfo()
 * Synopsis:  Fetch a sequence's info, using SSI indexing.
 *
 * Purpose:   Fetch a sequence named (or accessioned) <key> from
 *            the repositionable, open sequence file <sqfp>, reading
 *            all info except the sequence (and secondary structure).
 *            The open <sqfp> must have an open SSI index.
 *            The sequence info is returned in <sq>.
 *
 * Returns:   <eslOK> on soccess.
 *            <eslEINVAL> if no SSI index is present, or if <sqfp> can't
 *            be repositioned.
 *            <eslENOTFOUND> if <source> isn't found in the file.
 *            <eslEFORMAT> if either the index file or the sequence file
 *            can't be parsed, because of unexpected format issues.
 *       
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sqio_FetchInfo(ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq)
{
  int status;

  if (sqfp->ssi == NULL) ESL_FAIL(eslEINVAL, sqfp->errbuf, "No SSI index for %s; can't fetch subsequences", sqfp->filename);
  if ((status = esl_sqfile_PositionByKey(sqfp, key)) != eslOK) return status;
  if ((status = esl_sqio_ReadInfo(sqfp, sq))         != eslOK) return status;
  return eslOK;
}
  

/* Function:  esl_sqio_FetchSubseq()
 * Synopsis:  Fetch a subsequence, using SSI indexing.
 *
 * Purpose:   Fetch subsequence <start..end> from a sequence named (or
 *            accessioned) <source>, in the repositionable, open sequence file <sqfp>.
 *            The open <sqfp> must have an SSI index. Put the
 *            subsequence in <sq>. 
 *            
 *            As a special case, if <end> is 0, the subsequence is
 *            fetched all the way to the end, so you don't need to
 *            look up the sequence length <L> to fetch a suffix.
 *            
 *            The caller may want to rename/reaccession/reannotate the
 *            subsequence.  Upon successful return, <sq->name> is set
 *            to <source/start-end>, and <sq->source> is set to
 *            <source> The accession and description <sq->acc> and
 *            <sq->desc> are set to the accession and description of
 *            the source sequence.
 *            
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> if no SSI index is present, or if <sqfp> can't
 *            be repositioned.
 *            <eslENOTFOUND> if <source> isn't found in the file.
 *            <eslEFORMAT> if either the index file or the sequence file
 *            can't be parsed, because of unexpected format issues.
 *            <eslERANGE> if the <start..end> coords don't lie entirely
 *            within the <source> sequence.
 *
 * Throws:    <eslEMEM> on allocation errors.
 */
int
esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq)
{
  uint16_t fh; /* SSI file handle */
  off_t    r_off, d_off;
  int64_t  L;
  int64_t  actual_start;
  int64_t  nskip;
  int64_t  nres;
  int64_t  n;
  int      status;

  if (sqfp->ssi == NULL) ESL_FAIL(eslEINVAL, sqfp->errbuf, "No SSI index for %s; can't fetch subsequences", sqfp->filename);

  /* Find sequence info in the index */
  status = esl_ssi_FindSubseq(sqfp->ssi, source, start, &fh, &r_off, &d_off, &L, &actual_start);
  if      (status == eslENOTFOUND) ESL_FAIL(status, sqfp->errbuf, "Didn't find sequence %s in the index", source);
  else if (status == eslEFORMAT)   ESL_FAIL(status, sqfp->errbuf, "Failure reading SSI index; corrupt or bad format");
  else if (status == eslERANGE)    ESL_FAIL(status, sqfp->errbuf, "Requested start %" PRIi64 " isn't in the sequence %s", start, source);
  else if (status != eslOK)        ESL_FAIL(status, sqfp->errbuf, "Unexpected failure in finding subseq offset");

  /* The special case of end=0, asking for suffix fetch */
  if (end == 0) end = L;

  /* Validate coords if we can */
  if (start > end)       ESL_FAIL(eslERANGE, sqfp->errbuf, "Subsequence start %" PRIi64 " is greater than end %" PRIi64 "\n", start, end);
  if (L > 0 && end > L)  ESL_FAIL(eslERANGE, sqfp->errbuf, "Subsequence end %" PRIi64 " is greater than length %" PRIi64 "\n", end, L);

  /* Position the file at the record header; read the header info */
  status = esl_sqfile_Position(sqfp, r_off);
  if      (status == eslEOF)    ESL_FAIL(status, sqfp->errbuf, "Position appears to be off the end of the file");
  else if (status == eslEINVAL) ESL_FAIL(status, sqfp->errbuf, "Sequence file is not repositionable");
  else if (status != eslOK)     ESL_FAIL(status, sqfp->errbuf, "Failure in positioning sequence file");
  if ((status = sqfp->parse_header(sqfp, sq)) != eslOK) return status;

  /* Position the file close to the subseq: either at the start of the line
   * where the subseq starts, or exactly at the residue.
   */
  if (d_off != 0) 
    {
      status = esl_sqfile_Position(sqfp, d_off);
      if      (status == eslEOF)    ESL_FAIL(eslERANGE, sqfp->errbuf, "Position appears to be off the end of the file");
      else if (status == eslEINVAL) ESL_FAIL(status,    sqfp->errbuf, "Sequence file is not repositionable");
      else if (status != eslOK)     ESL_FAIL(status,    sqfp->errbuf, "Failure in positioning sequence file");
    }
  /* even if we didn't have a data offset, we're positioned at the
   * start of the sequence anyway, because we parsed the full header 
   */
  nskip = start - actual_start; /* how many residues do we still need to skip to reach start       */
  nres  = end - start + 1;      /* how many residues do we need to read as subseq                  */

  if ((status = esl_sq_GrowTo(sq, nres)) != eslOK) return status;
  status = read_nres(sqfp, sq, nskip, nres, &n);
  if (status != eslOK || n < nres) ESL_EXCEPTION(eslEINCONCEIVABLE, "Failed to fetch subsequence residues -- corrupt coords?");

  /* Set the coords */
  sq->start = start;
  sq->end   = end;
  sq->C     = 0;
  sq->W     = sq->n;
  sq->L     = (L > 0 ? L : -1);
  esl_sq_FormatName(sq, "%s/%" PRId64 "-%" PRId64, source, start, end);
  esl_sq_SetSource (sq, source);
  return eslOK;
}  
/*------------- end, random sequence access with SSI -------------------*/


/*****************************************************************
 * 7. Sequence database caching.
 *****************************************************************/ 

/* Function:  esl_sqfile_Cache()
 * Synopsis:  Read a database into memory.
 *
 * Purpose:   Read an entire database into memory building a cached
 *            structure <ESL_SQCACHE>.  The cache structure has basic
 *            information about the database, ie number of sequences
 *            number of residues, etc.
 *
 *            All sequences <ESL_SQ> are in a memory array <sq_list>.
 *            The number of elements in the list is <seq_count>.  The
 *            header pointers, ie name, acc and desc are pointers into
 *            the <header_mem> buffer.  All digitized sequences are pointers
 *            into the <residue_mem> buffer.
 *
 * Returns:   <eslOK> on success.
 *            
 *            Returns <eslEFORMAT> if a parse error is encountered in
 *            trying to read the sequence file.
 *            
 *            Returns <eslENODATA> if the file appears to be empty.
 *
 * Throws:    <eslEMEM> on allocation error;
 */
int  
esl_sqfile_Cache(const ESL_ALPHABET *abc, const char *seqfile, int fmt, const char *env, ESL_SQCACHE **ret_sqcache)
{
  int          status;

  int          n;

  uint32_t     len;
  uint32_t     max;
  uint32_t     count;

  uint64_t     res_size = 1;
  uint64_t     hdr_size = 1;

  ESL_SQFILE  *sqfp    = NULL;

  ESL_SQ      *c       = NULL;
  ESL_SQ      *sq      = NULL;
  ESL_SQCACHE *cache   = NULL;

  ESL_DSQ     *res_ptr = NULL;
  char        *hdr_ptr = NULL;

  /* open the database */
  status = esl_sqfile_OpenDigital(abc, seqfile, fmt, env, &sqfp);
  if (status != eslOK) return status;

  /* if we can't rewind the database, stop now.  */
  if (!esl_sqfile_IsRewindable(sqfp)) return eslFAIL;

  /* loop through the database reading all the sequences */
  max = 0;
  count = 0;
  sq  = esl_sq_CreateDigital(abc);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    ++count;

    res_size += sq->n + 1;
    if (sq->n > max) max = sq->n;

    len = strlen(sq->name);
    if (len > 0) ++len;
    hdr_size += len;

    len = strlen(sq->acc);
    if (len > 0) ++len;
    hdr_size += len;

    len = strlen(sq->desc);
    if (len > 0) ++len;
    hdr_size += len;

    esl_sq_Reuse(sq);
  }
  if (status != eslEOF) goto ERROR;
  
  /* now that the database information is know, allocate the memory to
   * hold the data.  different memory blocks will be used to hold the
   * redisues and header info.  the idea is that since the header info,
   * ie, name, acc, etc is used infrenquently (only when there is a hit)
   * if some pages need to be swapped out, hopefully it will be the
   * header pages first leaving the sequences in memory.
   */
  ESL_ALLOC(cache, sizeof(ESL_SQCACHE));

  cache->filename    = NULL;
  cache->sq_list     = NULL;
  cache->residue_mem = NULL;
  cache->header_mem  = NULL;

  cache->abc         = abc;
  cache->format      = fmt;
  cache->seq_count   = count;
  cache->res_count   = res_size;
  cache->max_seq     = max;

  cache->res_size    = res_size + 2;
  cache->hdr_size    = hdr_size;

  ESL_ALLOC(cache->filename, strlen(seqfile) + 1);
  strcpy(cache->filename, seqfile);

  ESL_ALLOC(cache->sq_list, sizeof(ESL_SQ) * (count + 1));

  /* different memory blocks will be used to hold the residues and header
   * info.  the idea is that since the header info, ie, name, acc, etc.
   * is used infrenquently (only when there is a hit) if some pages need
   * to be swapped out, hopefully it will be the header pages first
   * leaving the sequences in memory.
   */
  ESL_ALLOC(cache->residue_mem, res_size + 2);
  ESL_ALLOC(cache->header_mem, hdr_size);

  hdr_ptr  = cache->header_mem;
  *hdr_ptr = 0;

  res_ptr  = cache->residue_mem;
  *res_ptr = eslDSQ_SENTINEL;

  /* loop through the database filling in the cache */
  n = 0;
  esl_sqfile_Position(sqfp, 0);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    c = cache->sq_list + n;

    /* if header fields have been defined, copy them to the cache */
    c->name = hdr_ptr;
    if (sq->name[0] != 0) {
      c->name = hdr_ptr + 1;
      strcpy(c->name, sq->name);
      hdr_ptr += strlen(sq->name) + 1;
    }

    c->acc = hdr_ptr;
    if (sq->acc[0] != 0) {
      c->acc = hdr_ptr + 1;
      strcpy(c->acc, sq->acc);
      hdr_ptr += strlen(sq->acc) + 1;
    }

    c->desc = hdr_ptr;
    if (sq->desc[0] != 0) {
      c->desc = hdr_ptr + 1;
      strcpy(c->desc, sq->desc);
      hdr_ptr += strlen(sq->desc) + 1;
    }

    c->tax_id = sq->tax_id;
    c->seq    = NULL;
    c->ss     = NULL;
    c->nxr    = 0;
    c->xr_tag = NULL;
    c->xr     = NULL;

    /* copy the digitized sequence */
    memcpy(res_ptr + 1, sq->dsq + 1, sq->n + 1);
    c->dsq   = res_ptr;
    c->n     = sq->n;
    res_ptr += sq->n + 1;

    /* Coordinate info */
    c->start = sq->start;
    c->end = sq->end;
    c->C = sq->C;
    c->W = sq->W;
    c->L = sq->L;

    c->source = NULL;

    /* allocated lengths */
    c->nalloc   = -1;
    c->aalloc   = -1;
    c->dalloc   = -1;
    c->salloc   = -1;
    c->srcalloc = -1;

    /* Disk offset bookkeeping */
    c->idx  = n;
    c->roff = sq->roff;
    c->hoff = sq->hoff;
    c->doff = sq->doff;
    c->eoff = sq->eoff;

    c->abc = abc;

    esl_sq_Reuse(sq);
    ++n;
  }
  if (status != eslEOF) goto ERROR;

  /* add on last empty sequence */
  c = cache->sq_list + count;
  *(res_ptr + 1) = eslDSQ_SENTINEL;

  c->name     = hdr_ptr;
  c->acc      = hdr_ptr;
  c->desc     = hdr_ptr;

  c->tax_id   = -1;
  c->seq      = NULL;
  c->ss       = NULL;
  c->nxr      = 0;
  c->xr_tag   = NULL;
  c->xr       = NULL;

  c->dsq      = res_ptr;
  c->n        = 0;

  c->start    = 0;
  c->end      = 0;
  c->C        = 0;
  c->W        = 0;
  c->L        = -1;

  c->source   = NULL;

  c->nalloc   = -1;
  c->aalloc   = -1;
  c->dalloc   = -1;
  c->salloc   = -1;
  c->srcalloc = -1;

  c->idx      = count;
  c->roff     = -1;
  c->hoff     = -1;
  c->doff     = -1;
  c->eoff     = -1;

  c->abc      = NULL;
 
  if (sq != NULL) esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);

  *ret_sqcache = cache;

  return eslOK;

ERROR:
  if (sq != NULL) esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);

  esl_sqfile_Free(cache);

  return status;
}

/* Function:  esl_sqfile_Free()
 * Synopsis:  Free a cached database <ESL_SQCACHE>.
 *
 * Purpose:   Frees all the memory used to cache the sequence database.
 *
 * Returns:   none.
 */
void  
esl_sqfile_Free(ESL_SQCACHE *sqcache)
{
  if (sqcache == NULL) return;

  if (sqcache->filename    != NULL) free(sqcache->filename);
  if (sqcache->sq_list     != NULL) free(sqcache->sq_list);
  if (sqcache->residue_mem != NULL) free(sqcache->residue_mem);
  if (sqcache->header_mem  != NULL) free(sqcache->header_mem);

  sqcache->abc         = NULL;
  sqcache->filename    = NULL;
  sqcache->sq_list     = NULL;
  sqcache->residue_mem = NULL;
  sqcache->header_mem  = NULL;

  free(sqcache);
}
/*---------------- end, sequence database caching ---------------*/





/*****************************************************************
 *  8. Internal functions (all formats)
 *****************************************************************/

/* loadmem() 
 *
 * Load the next block of data from stream into mem buffer,
 * either concatenating to previous buffer (if we're recording) or
 * overwriting (if not). 
 * 
 * This block is loaded at sqfp->mem + sqfp->mpos.
 * 
 * Upon return:
 * sqfp->mem     now contains up to eslSQFILE_READBUFSIZE more chars
 * sqfp->mpos    is position of first byte in newly read block
 * sqfp->allocm  may have increased by eslSQFILE_READBUFSIZE, if we concatenated
 * sqfp->mn      is # of chars in <mem>; <mn-1> is pos of last byte in new block
 * 
 * Returns <eslEOF> (and mpos == mn) if no new data can be read;
 * Returns <eslOK>  (and mpos < mn) if new data is read. 
 * Throws <eslEMEM> on allocation error.
 */
static int
loadmem(ESL_SQFILE *sqfp)
{
  void *tmp;
  int   n = 0;
  int   status;

  if (sqfp->do_buffer)
    {
      sqfp->mpos = 0;
      sqfp->mn   = 0;
    }
  else if (sqfp->is_recording == TRUE)
    {
      if (sqfp->mem == NULL) sqfp->moff = ftello(sqfp->fp);        /* first time init of the offset */
      ESL_RALLOC(sqfp->mem, tmp, sizeof(char) * (sqfp->allocm + eslSQFILE_READBUFSIZE));
      sqfp->allocm += eslSQFILE_READBUFSIZE;
      n = fread(sqfp->mem + sqfp->mpos, sizeof(char), eslSQFILE_READBUFSIZE, sqfp->fp);
      sqfp->mn += n;
    }
  else
    {
      if (sqfp->mem == NULL) {
        ESL_ALLOC(sqfp->mem, sizeof(char) * eslSQFILE_READBUFSIZE);
        sqfp->allocm = eslSQFILE_READBUFSIZE;
      }
      sqfp->is_recording = -1;/* no more recording is possible now */
      sqfp->mpos = 0;
      sqfp->moff = ftello(sqfp->fp);
      n = fread(sqfp->mem, sizeof(char), eslSQFILE_READBUFSIZE, sqfp->fp); /* see note [1] below */
      sqfp->mn   = n;
    }
  return (n == 0 ? eslEOF : eslOK);

 ERROR:
  return status;
}

/* [1] Be alert for a possible problem above in that fread().
 *     Farrar had inserted an alternative case as follows:
 *     "If we are reading from stdin, buffered read cannot be used
 *      because it will block until EOF or the buffer is full, ie
 *      eslSQFILE_READBUFSIZE characters have been read.  Usually this would
 *      not be a problem, unless stdin is from a pipe.  In that case
 *      if the sequence is less than eslSQFILE_READBUFSIZE we would block.
 *
 *      NOTE:  any changes to the IO stream sqfp->fp, such as fseek, 
 *      might not have any affect on the file descriptor for the stream.
 *
 *   if (sqfp->do_stdin) {
 *     n = read(fileno(sqfp->fp), sqfp->mem, eslSQFILE_READBUFSIZE);
 *   } else {
 *   ...
 *  
 * but that's a bug, because you can't mix read and fread;
 * the i17-stdin.pl test fails, in particular.
 */




/* loadbuf()
 * Set sqfp->buf to contain next line of data, or point to next block.
 * This might just mean working with previously buffered memory in <sqfp->mem>
 * or might require reading new data from <sqfp->fp>.
 *
 * Reset sqfp->boff to be the position of the start of the block/line.
 * Reset sqfp->bpos to 0.
 * Reset sqfp->nc to the number of chars (bytes) in the new block/line.
 * Returns eslOK on success; eslEOF if there's no more data in the file.
 * (sqfp->nc == 0 is the same as eslEOF: no data in the new buffer.)
 * Can throw an <eslEMEM> error.
 */
static int
loadbuf(ESL_SQFILE *sqfp)
{
  void *tmp;
  char *nlp;
  int   n;
  int   status = eslOK;

  if (! sqfp->is_linebased)
    {
      if (sqfp->mpos >= sqfp->mn) {
        if ((status = loadmem(sqfp)) == eslEMEM) return status;
      }
      sqfp->buf    = sqfp->mem  + sqfp->mpos;
      sqfp->boff   = sqfp->moff + sqfp->mpos;
      sqfp->balloc = 0;
      sqfp->bpos   = 0;
      sqfp->nc     = sqfp->mn - sqfp->mpos;
      sqfp->mpos  += sqfp->nc;
    }
  else
    { /* Copy next line from <mem> into <buf>. Might require new load(s) into <mem>. */
      if (sqfp->mpos >= sqfp->mn) {
        if ((status = loadmem(sqfp)) == eslEMEM) return status;
      }
      sqfp->boff = sqfp->moff + sqfp->mpos;      
      sqfp->nc   = 0;
      nlp        = memchr(sqfp->mem + sqfp->mpos, '\n', sqfp->mn - sqfp->mpos);
      while (nlp == NULL) 
        {
          n = sqfp->mn - sqfp->mpos;
          while (sqfp->nc + n + 1 > sqfp->balloc) { /* +1: it'll hold the terminal \0 */
            ESL_RALLOC(sqfp->buf, tmp, sizeof(char) * (sqfp->balloc + eslSQFILE_READBUFSIZE));
            sqfp->balloc += eslSQFILE_READBUFSIZE;
          }
          memcpy(sqfp->buf + sqfp->nc, sqfp->mem + sqfp->mpos, n);
          sqfp->mpos += n;
          sqfp->nc   += n;
          status = loadmem(sqfp);
          if      (status == eslEOF) { break; }
          else if (status != eslOK)  return status;
          nlp = memchr(sqfp->mem + sqfp->mpos, '\n', sqfp->mn - sqfp->mpos);
        }
      if (status != eslEOF) {
        n = nlp - (sqfp->mem + sqfp->mpos) + 1; /* inclusive of \n */
        if (sqfp->nc + n + 1 > sqfp->balloc) {
          ESL_RALLOC(sqfp->buf, tmp, sizeof(char) * (sqfp->balloc + eslSQFILE_READBUFSIZE));
          sqfp->balloc += eslSQFILE_READBUFSIZE;
        }
        memcpy(sqfp->buf + sqfp->nc, sqfp->mem + sqfp->mpos, n);
        sqfp->mpos += n;
        sqfp->nc   += n;
      }
      sqfp->bpos  = 0;
      sqfp->buf[sqfp->nc] = '\0';
    }
  return (sqfp->nc == 0 ? eslEOF : eslOK);

 ERROR:
  return status;
}

/* nextchar()
 * 
 * Load next char from sqfp->buf into <*ret_c> and sets sqfp->bpos to
 * its position; usually this is c = sqfp->buf[++sqfp->bpos], but
 * we will refill the buffer w/ fresh fread() when needed, in which
 * case c =  sqfp->buf[0] and sqfp->bpos = 0.
 * 
 * Returns <eslOK> on success.
 * Return  <eslEOF> if we ran out of data in <sqfp>.
 * May throw an <eslEMEM> error.
 */
static int
nextchar(ESL_SQFILE *sqfp, char *ret_c)
{
  int status;

  sqfp->bpos++;
  if (sqfp->nc == sqfp->bpos && (status = loadbuf(sqfp)) != eslOK) return status;
  *ret_c = sqfp->buf[sqfp->bpos];
  return eslOK;
}

/* seebuf()
 * 
 * Examine and validate the current buffer <sqfp->buf> from its
 * current position <sqfp->bpos> until either the buffer ends (we run
 * out of characters) or the sequence data ends (we see whatever
 * character indicates EOD in this format) or we've seen <maxn>
 * residues. If <maxn> is passed as -1, parse the entire buffer,
 * without a residue limit.
 * 
 * There are three possible outcomes:
 *   <eslOK>:      The buffer is all residues that belong to the current
 *                 seq we're parsing (or chars we can ignore), at least
 *                 up to the <maxn> residue limit (if present).
 *   <eslEOD>:     Part of the buffer may be residues, but the current sequence
 *                 ends in this buffer (before <maxn> was reached).
 *   <eslEFORMAT>: Somewhere before we reached the end of the buffer or
 *                 the sequence record, we saw an illegal character.
 * 
 * On <eslOK>:
 *    *opt_nres    is the number of residues in the buffer (up to <maxn>)
 *    *opt_endpos  is sqfp->nc (off the end of the buffer by one)
 *    The caller will want to deal with the buffer, then load the next one.
 *    
 * On <eslEOD>: same as OK, except:
 *    *opt_endpos  is where sqfp->bpos *would* be at when we saw the EOD
 *                 signal (the next '>', in FASTA files) had we been parsing residues
 *    Therefore on EOD, the caller will want to deal with the <*opt_nres>
 *    residues in this buffer, then reposition the buffer by
 *    <sqfp->bpos = *opt_epos> (without reloading the buffer), so
 *    the next read will pick up there.
 *    
 * On <eslEFORMAT>:
 *    sqfp->errbuf  contains informative message about the format error.
 *    
 * seebuf() also handles linenumber and SSI bookkeeping in
 * <sqfp>. Every newline character seen increments <linenumber> (thus,
 * on EFORMAT return, linenumber is set to the line on which the bad
 * char occurred). <curbpl>,<currpl>,<prvbpl>,<prvrpl> keep track of # of bytes,
 * residues on the current,prev line; they keep state across calls to seebuf().
 * <bpl>,<rpl> are tracking whether there's a constant number of
 * bytes/residues per line; these are either -1 for "not set yet", 0
 * for "no, not constant", or a number > 0. Because of this bookkeeping, it's important
 * to make sure that <seebuf()> never counts the same byte twice (hence
 * the need for the <maxn> limit, which ReadWindow() uses.)
 */
static int
seebuf(ESL_SQFILE *sqfp, int64_t maxn, int64_t *opt_nres, int64_t *opt_endpos)
{
  int     bpos;
  int64_t nres  = 0;
  int64_t nres2 = 0;/* an optimization for determining lastrpl from nres, without incrementing lastrpl on every char */
  int     sym;
  ESL_DSQ x;
  int     lasteol;
  int     status  = eslOK;

  lasteol = sqfp->bpos - 1;
  if (maxn == -1) maxn = sqfp->nc; /* makes for a more efficient test. nc is a guaranteed upper bound on nres */

  for (bpos = sqfp->bpos; nres < maxn && bpos < sqfp->nc; bpos++)
    {
      sym = sqfp->buf[bpos];
      //printf ("nres: %d, bpos: %d  (%d)\n", nres, bpos, sym);
      if (!isascii(sym)) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": non-ASCII character %c in sequence", sqfp->linenumber, sym); 
      x   = sqfp->inmap[sym];

      if      (x <= 127) nres++;
      else if (x == eslDSQ_EOL)
        {
          if (sqfp->curbpl != -1) sqfp->curbpl += bpos - lasteol;
          if (sqfp->currpl != -1) sqfp->currpl += nres - nres2;
          nres2        += nres - nres2;

          /* If bpl/rpl are set, then nonterminal seq lines must match
           * them exactly, and terminal seq lines must fit within
           * them. We may not know bpl/rpl yet when we see a terminal
           * seq line (first sequence record(s) may have only one line
           * of sequence), so track maxima over all seq lines.
           */
          if (sqfp->curbpl                 > sqfp->maxbpl) sqfp->maxbpl = sqfp->curbpl;
          if (sqfp->currpl                 > sqfp->maxrpl) sqfp->maxrpl = sqfp->currpl;
          if (sqfp->curbpl - sqfp->currpl > sqfp->maxgap) sqfp->maxgap = sqfp->curbpl - sqfp->currpl;

          if (sqfp->rpl != 0 && sqfp->prvrpl != -1) {
            if      (sqfp->rpl    == -1)         sqfp->rpl = sqfp->prvrpl; /* init  */
            else if (sqfp->prvrpl != sqfp->rpl) sqfp->rpl = 0;             /* inval */
          }
          if (sqfp->bpl != 0 && sqfp->prvbpl != -1) {
            if      (sqfp->bpl    == -1)         sqfp->bpl = sqfp->prvbpl; /* init  */
            else if (sqfp->prvbpl != sqfp->bpl) sqfp->bpl = 0;             /* inval */
          }

          if (sqfp->rpl > 0 && sqfp->maxrpl > sqfp->rpl) sqfp->rpl = 0;   /* inval: some line is longer than rpl */
          if (sqfp->bpl > 0 && sqfp->maxbpl > sqfp->bpl) sqfp->bpl = 0;   /* inval: some line is longer than bpl */
          if (sqfp->bpl > 0 && sqfp->bpl == sqfp->rpl+1 && sqfp->maxgap > 1)  sqfp->bpl = 0;   /* inval: bpl==rpl+1 implies residue resolution, but some line has extra nonresidue bytes */

          sqfp->prvbpl  = sqfp->curbpl;
          sqfp->prvrpl  = sqfp->currpl;
          sqfp->curbpl  = 0;
          sqfp->currpl  = 0;
          lasteol       = bpos;
          if (sqfp->linenumber != -1) sqfp->linenumber++;
        }
      else if (x == eslDSQ_ILLEGAL) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": illegal character %c", sqfp->linenumber, sym);
      else if (x == eslDSQ_EOD)     { status = eslEOD; break; }
      else if (x != eslDSQ_IGNORED) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "inmap corruption?");
    }

  if (sqfp->curbpl != -1) sqfp->curbpl += bpos - lasteol - 1;
  if (sqfp->currpl != -1) sqfp->currpl += nres - nres2;

  /* Watch out for a final line with no EOL.
   * This may also be a partial line on a buffer boundary, but that's ok;
   * curbpl/currpl can be underestimated for this line and we won't invalidate bpl/rpl indexing.
   * The +1 in maxgap is for the EOL that we didn't see (yet, or if EOF, ever).
   */
  if (sqfp->curbpl != -1) {
    if (sqfp->curbpl                     > sqfp->maxbpl) sqfp->maxbpl = sqfp->curbpl;
    if (sqfp->currpl                     > sqfp->maxrpl) sqfp->maxrpl = sqfp->currpl;
    if (sqfp->curbpl - sqfp->currpl + 1  > sqfp->maxgap) sqfp->maxgap = sqfp->curbpl - sqfp->currpl + 1;
  }
  if (sqfp->rpl > 0 && sqfp->maxrpl > sqfp->rpl) sqfp->rpl = 0;
  if (sqfp->bpl > 0 && sqfp->maxbpl > sqfp->bpl) sqfp->bpl = 0;
  if (sqfp->bpl > 0 && sqfp->bpl == sqfp->rpl+1 && sqfp->maxgap > 1) sqfp->bpl = 0;

  if (opt_nres   != NULL) *opt_nres   = nres;
  if (opt_endpos != NULL) *opt_endpos = bpos;
  return status;
}

/* addbuf() 
 * Add <nres> residues from the current buffer <sqfp->buf> to <sq>.
 * This is designed to work when we're constructing a complete
 * sequence (add the whole buffer); when we're adding a suffix
 * of the buffer (<sqfp->bpos> is skipped ahead already);
 * or when we're adding a prefix of the buffer (terminating a subseq
 * or window load).
 * 
 * The caller must know that there are at least <nres> residues in
 * this buffer, and that all the characters are valid in the
 * format and alphabet, via a previous call to <seebuf()>. 
 * 
 * The caller also must have already allocated <sq> to hold at least
 * <nres> more residues.
 * 
 * On input:
 *   sqfp->buf[]  contains an fread() buffer
 *   sqfp->bpos   is set to where we're going to start parsing residues
 *   sqfp->nc     is the length of <buf>
 *   
 * On return:
 *   sqfp->buf[]  still contains the same buffer (no new freads here)
 *   sqfp->bpos   is set after the last residue we parsed 
 *   sq->seq/dsq  now holds <nres> new residues
 *   sq->n        is incremented by <nres>
 */
static void
addbuf(ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nres)
{
  ESL_DSQ x;

  if (sq->dsq != NULL) 
    {
      while (nres) {
        x  = sq->abc->inmap[(int) sqfp->buf[sqfp->bpos++]];
        if (x <= 127) { nres--; sq->dsq[++sq->n] = x; }
      } /* we skipped IGNORED, EOL. EOD, ILLEGAL don't occur; seebuf() already checked  */
    } 
  else
    {
      while (nres) {
        x   = sqfp->inmap[(int) sqfp->buf[sqfp->bpos++]];
        if (x <= 127) { nres--; sq->seq[sq->n++] = x; }
      }
    }
}

/* skipbuf() 
 * Like addbuf(), but we skip <nskip> residues instead of
 * reading them.
 */
static void
skipbuf(ESL_SQFILE *sqfp, int64_t nskip)
{
  ESL_DSQ x;

  while (nskip) {
    x  = sqfp->inmap[(int) sqfp->buf[sqfp->bpos++]];
    if (x <= 127) nskip--;/* skip IGNORED, EOL. */
  }
}


/* skip_whitespace()
 * Like skipbuf(), but instead of skipping a fixed number of
 * residues, skip forward until one of three conditions is met:
 *
 * (1) end of the sequence record (a character indicating
 *     the beginning of a new sequence); set sqfp->bpos
 *     to the beginning of the new record, and return eslEOD;
 * (2) a non-whitespace character in the current sequence is
 *     reached that does not indicate the end of a sequence
 *     record; set sqfp->bpos to that character's position,
 *     and return eslOK;
 * (3) end of file;  return eslEOF.
 *
 */
static int
skip_whitespace(ESL_SQFILE *sqfp)
{
  int c;
  ESL_DSQ x;

  if (sqfp->nc == 0) return eslEOF;

  /* if at end of buffer, reload it */
  if (sqfp->bpos == sqfp->nc)
    if (loadbuf(sqfp) == eslEOF) return eslEOF;

  c = (int) sqfp->buf[sqfp->bpos];
  x  = sqfp->inmap[c];

  while ( isspace(c) ) {
    sqfp->bpos++;

    /* if at end of buffer, reload it */
    if (sqfp->bpos == sqfp->nc)
      if ( loadbuf(sqfp) == eslEOF) return eslEOF;

    c = (int) sqfp->buf[sqfp->bpos];
    x  = sqfp->inmap[c];
  }
  if (x == eslDSQ_EOD) return eslEOD;
  return eslOK;
}



/* read_nres()
 * Read the next <nres> residues from <sqfp> after skipping <nskip> residues, then stop.
 * 
 * Returns <eslOK> and <0 < *ret_actual_nres <= nres> if it succeeded, and
 *                 there's more residues in the current seq record.
 * Returns <eslEOD> and <*ret_actual_nres == 0> if no more residues are
 *                 seen in the sequence record. 
 * 
 * Even on <eslEOD>, the <dsq/seq> is appropriately terminated here,
 * and <sq->n> is left the way it was (no new residues added - but there
 * may have been saved context C from a previous window).
 *
 * Returns <eslEFORMAT> on any parsing problem, and <sqfp->errbuf> is set.
 *
 * On <eslOK>, sqfp->bpos is positioned on the next character past the last residue we store;
 * on <eslEOD>, sqfp->bpos is positioned for reading the next sequence.
 * 
 * FetchSubseq() uses this with <nskip>, <nres>, and expects an
 * <eslOK> with <*opt_actual_nres = nres>. On <EOD>, or if fewer than
 * <nres> residues are obtained, the coords must've been screwed up,
 * because we didn't read the whole subseq we asked for.
 *
 * ReadWindow() on forward strand uses this with <nskip=0>, <nres=W>.
 * The last window might normally return <eslEOD> with
 * <*ret_actual_nres == 0>, and now <sqfp->bpos> is positioned at the
 * start of the next sequence on <EOD>, and at the next residue on
 * <OK>.
 * 
 * ReadWindow() in reverse complement acts like a subseq fetch.
 * 
 */
static int
read_nres(ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nskip, int64_t nres, int64_t *opt_actual_nres)
{
  int64_t n;
  int64_t epos;
  int64_t actual_nres = 0;
  int     status      = eslOK;
  
  status = seebuf(sqfp, nskip+nres, &n, &epos);
  while (status == eslOK && nskip - n > 0) {
    nskip   -= n;
    if ((status = loadbuf(sqfp)) == eslEOF) break;
    status = seebuf(sqfp, nskip+nres, &n, &epos);
  }
  
  if         (status == eslEOF) { 
    if (! sqfp->eof_is_ok) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Premature EOF before end of seq record");
    if (nskip > 0)         ESL_EXCEPTION(eslECORRUPT, "premature EOD while trying to skip residues"); 
    n = 0;
  } else if  (status == eslEOD) { 
    if (n < nskip)         ESL_EXCEPTION(eslECORRUPT, "premature EOD while trying to skip residues"); 
  } else if  (status != eslOK) 
    return status;

  skipbuf(sqfp, nskip); 
  n -= nskip; 

  while (status == eslOK && nres - n > 0) 
    {
      addbuf(sqfp, sq, n);
      actual_nres += n;
      nres        -= n;
      if ((status = loadbuf(sqfp)) == eslEOF) break;
      status = seebuf(sqfp, nres, &n, &epos);
    }

  if (status == eslEOF) { 
    if (! sqfp->eof_is_ok) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Premature EOF before end of seq record");
    n = 0;
  } else if  (status == eslEFORMAT) { return status; }

  n = ESL_MIN(nres, n); 
  addbuf(sqfp, sq, n);   /* bpos now at last residue + 1 if OK/EOD, 0 if EOF  */
  actual_nres += n;

  if (sq->dsq != NULL) sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else                 sq->seq[sq->n]   = '\0';
  
  if (status == eslEOD) { 
    sqfp->bpos = epos; 
  }

  if (opt_actual_nres != NULL) *opt_actual_nres = actual_nres;
  return (actual_nres == 0 ? eslEOD : eslOK);
}

/* convert_sq_to_msa()
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
  int      x;        /* counter for extra-residue markups */
  int      status;

  if (sq->dsq != NULL) 
    { 
      if ((msa = esl_msa_CreateDigital(sq->abc, 1, sq->n)) == NULL) { status = eslEMEM; goto ERROR; }
    } 
  else if ((msa = esl_msa_Create(1, sq->n)) == NULL) { status = eslEMEM; goto ERROR; }

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

  if (sq->dsq != NULL) esl_dsq_Copy(sq->dsq, sq->n, msa->ax[0]);
  else                 strcpy(msa->aseq[0], sq->seq);
  
  if (sq->ss != NULL)
    {
      ESL_ALLOC(msa->ss, sizeof(char *) * 1);

      if (sq->dsq != NULL) {	/* sq->ss is 1..L in digital mode; but msa->ss is always 0..L-1 */
	if      ((status = esl_strdup(sq->ss+1, -1, &(msa->ss[0]))) != eslOK) goto ERROR;
      } else if ((status = esl_strdup(sq->ss,   -1, &(msa->ss[0]))) != eslOK) goto ERROR;     	
    }

  if (sq->nxr > 0) {
    msa->ngr = sq->nxr;
    ESL_ALLOC(msa->gr,     sizeof(char **) * msa->ngr);    
    ESL_ALLOC(msa->gr_tag, sizeof(char  *) * msa->ngr);

    for (x = 0; x < msa->ngr; x ++) {
      ESL_ALLOC(msa->gr[x],     sizeof(char *));  
      ESL_ALLOC(msa->gr_tag[x], sizeof(char));
   
      if (sq->dsq != NULL) {	/* sq->xr is 1..L in digital mode; but msa->gr is always 0..L-1 */
        if      ((status = esl_strdup(sq->xr[x]+1, -1, &(msa->gr[x][0]))) != eslOK) goto ERROR;
      } else if ((status = esl_strdup(sq->xr[x],   -1, &(msa->gr[x][0]))) != eslOK) goto ERROR;     	

      if ((status = esl_strdup(sq->xr_tag[x], -1, &(msa->gr_tag[x]))) != eslOK) goto ERROR;     	  
    }
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
/*-------- end, internal functions (all formats)  ---------------*/



/*****************************************************************
 * 9. Internal functions (FASTA)
 *****************************************************************/

static void
config_fasta(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased = FALSE;
  sqfp->eof_is_ok    = TRUE;
  sqfp->parse_header = &header_fasta;
  sqfp->parse_end    = &end_fasta;
}

static void
inmap_fasta(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
    sqfp->inmap['-']  = eslDSQ_ILLEGAL;                      // don't let the abc_inmap map the gap char; this is an ungapped file format
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  sqfp->inmap['*']  = '*';         /* accept * as a nonresidue/stop codon character */
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;/* DOS eol compatibility */
  sqfp->inmap['\n'] = eslDSQ_EOL;
  sqfp->inmap['>']  = eslDSQ_EOD;
  /* \n is special - fasta reader detects it as an eol */
}


/* header_fasta()
 * 
 * sqfp->buf[sqfp->bpos] is sitting at the start of a FASTA record, or
 * at a space before it (in which case we'll advance, skipping whitespace,
 * until a > is reached).
 * Parse the header line, storing name and description in <sq>.
 * 
 * On success, returns <eslOK> and:
 *    sq->name contains sequence name (and may have been reallocated, changing sq->nalloc)
 *    sq->desc contains description line (and may have been reallocated, changing sq->dalloc)
 *    sq->roff has been set to the record offset
 *    sq->doff has been set to the data offset (start of sequence line)
 *    sqfp->buf[sqfp->bpos] is sitting at the start of the seq line.
 *    sqfp->currpl,curbpl set to 0, to start bookkeeping data line lengths 
 *
 * If no more seqs are found in the file, returns <eslEOF>.
 * On parse failure, return <eslEFORMAT>, leaves as mesg in sqfp->errbuf.
 *    
 * May also throw <eslEMEM> on allocation errors.
 */
static int
header_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  char  c;
  int   status = eslOK;
  void *tmp;
  int   pos;

  /* make sure there are characters in the buffer */
  if (sqfp->nc == sqfp->bpos && (status = loadbuf(sqfp)) != eslOK) return status;

  c =  sqfp->buf[sqfp->bpos];
  while (status == eslOK && isspace(c)) status = nextchar(sqfp, &c); /* skip space (including \n) */

  if (status == eslEOF) return eslEOF;

  if (status == eslOK && c == '>') {    /* accept the > */
    sq->roff = sqfp->boff + sqfp->bpos; /* store SSI record offset */
    status = nextchar(sqfp, &c);
  } else if (c != '>') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": unexpected char %c; expected FASTA to start with >", sqfp->linenumber, c);
  
  while (status == eslOK && (c == '\t' || c == ' ')) status = nextchar(sqfp, &c); /* skip space */

  /* Store the name (space delimited) */
  pos = 0;
  while (status == eslOK && ! isspace(c))
    {
      sq->name[pos++] = c;
      if (pos == sq->nalloc-1) { ESL_RALLOC(sq->name, tmp, sq->nalloc*2); sq->nalloc*=2; }
      status = nextchar(sqfp, &c); 
    }
  if (pos == 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": no FASTA name found", sqfp->linenumber);
  sq->name[pos] = '\0';
  
  while (status == eslOK &&  (c == '\t' || c == ' ')) status = nextchar(sqfp, &c);   /* skip space */

  /* Store the description (end-of-line delimited) */
  /* Patched to deal with NCBI NR desclines: delimit by ctrl-A (0x01) too. [SRE:H1/82] */
  pos = 0;
  while (status == eslOK && c != '\n' && c != '\r' && c != 1)
    {
      sq->desc[pos++] = c;
      if (pos == sq->dalloc-1) { ESL_RALLOC(sq->desc, tmp, sq->dalloc*2); sq->dalloc*= 2; }
      status = nextchar(sqfp, &c); 
    }
  sq->desc[pos] = '\0';

  /* Because of the NCBI NR patch, c might be0x01 ctrl-A now; skip to eol. 
   * (TODO: I'm worried about the efficiency of this nextchar() stuff. Revisit.)
   */
  while (status == eslOK && c != '\n' && c != '\r') 
    status = nextchar(sqfp, &c);
  sq->hoff = sqfp->boff + sqfp->bpos;
  
  while (status == eslOK && (c == '\n' || c == '\r')) status = nextchar(sqfp, &c); /* skip past eol (DOS \r\n, MAC \r, UNIX \n */
  if (status != eslOK && status != eslEOF) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Unexpected failure in parsing FASTA name/description line");
  /* Edge case: if the last sequence in the file is L=0, no residues, we are EOF now, not OK; but we'll return OK because we parsed the header line */

  sq->doff = sqfp->boff + sqfp->bpos;
  sqfp->prvrpl = sqfp->prvbpl = -1;
  sqfp->currpl = sqfp->curbpl = 0;
  sqfp->linenumber++;
  return eslOK;

 ERROR:
  return status;/* eslEMEM, from failed realloc */
}

static int 
end_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  if (sqfp->bpos < sqfp->nc) {
    if (sqfp->buf[sqfp->bpos] != '>') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Whoops, FASTA reader is corrupted");
    sq->eoff = sqfp->boff + sqfp->bpos - 1; /* this puts eoff at the last \n */
  } /* else, EOF, and we don't have to do anything. */
  return eslOK;
}
/*-------------- end, internal functions (FASTA) ----------------*/



/*****************************************************************
 * 10. Internal functions (EMBL|UniProt|TrEMBL)
 *****************************************************************/ 
/* EMBL and UniProt protein sequence database format.
 *   See: http://us.expasy.org/sprot/userman.html
 *   and: http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html#3
 * We use the same parser for both formats, so we have to be 
 * careful to only parse the conserved intersection of these two
 * very similar formats.
 */
static void
config_embl(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased      = TRUE;
  sqfp->eof_is_ok         = FALSE;/* records end with // */
  sqfp->parse_header      = &header_embl;
  sqfp->parse_end         = &end_embl;
}

static void
inmap_embl(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
    sqfp->inmap['-']  = eslDSQ_ILLEGAL;                      // don't let the abc_inmap map the gap char; this is an ungapped file format
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  for (x = '0'; x <= '9'; x++)
    sqfp->inmap[x] = eslDSQ_IGNORED;    /* EMBL DNA sequence format puts coordinates after each line */
  sqfp->inmap['*']  = '*';              /* accept * as a nonresidue/stop codon character */
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;  /* DOS eol compatibility */
  sqfp->inmap['/']  = eslDSQ_EOD;
}

/* header_embl()
 * 
 * See: http://us.expasy.org/sprot/userman.html
 * And: http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html#3
 * Our parser must work on the highest common denominator of EMBL DNA
 * and UniProt protein sequence files.
 *
 * sqfp->buf is the first (ID) line of the entry, or a blank line before
 * it (in which case we'll scan forwards skipping blank lines to find 
 * the ID line).
 * 
 * On success, returns <eslOK> and:
 *   sq->name  contains sequence name (and may have been reallocated, changing sq->nalloc)
 *   sq->acc   contains seq accession (and may have been reallocated, changing sq->aalloc)
 *   sq->desc  contains description line (and may have been reallocated, changing sq->dalloc)
 *   sq->roff  has been set to the record offset
 *   sq->doff  has been set to the data offset (start of sequence line)
 *   sqfp->buf is the first seq line.
 * 
 * If no more seqs are found in the file, returns <eslEOF>.
 * On parse failure, returns <eslEFORMAT>, leaves as mesg in sqfp->errbuf.
 * 
 * May also throw <eslEMEM> on allocation errors.
 */
static int
header_embl(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  char *s;
  char *tok;
  int   status;

  /* Find first line:
   * "Each entry must begin with an identification line (ID)..."
   * "The two-character line-type code that begins each line is always
   *  followed by three blanks..."
   */
  if (sqfp->nc == 0) return eslEOF;
  while (esl_str_IsBlank(sqfp->buf)) {
    if ((status = loadbuf(sqfp)) == eslEOF) return eslEOF; /* normal */
    else if (status != eslOK) return status; /* abnormal */
  } 

  /* ID line is defined as:
   *     ID   ENTRY_NAME DATA_CLASS; MOLECULE_TYPE; SEQUENCE_LENGTH.
   * We're only after the ENTRY_NAME.
   * Examples:
   *  ID   SNRPA_DROME    STANDARD;      PRT;   216 AA.
   *  ID   SNRPA_DROME             Reviewed;         216 AA.
   *  ID   X06347; SV 1; linear; mRNA; STD; HUM; 1209 BP.
   */
  if (strncmp(sqfp->buf, "ID   ", 5) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to find ID line", sqfp->linenumber);
  
  s = sqfp->buf+5;
  if ( esl_strtok(&s, " ;", &tok) != eslOK)
    ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to parse name on ID line", sqfp->linenumber);
  if ((status = esl_sq_SetName(sq, tok)) != eslOK) return status;
  sq->roff = sqfp->boff;/* record the offset of the ID line */
  
  /* Look for SQ line; parsing optional info as we go.
   */
  do {
    if ( loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to find SQ line", sqfp->linenumber);

    /* "The format of the AC line is:
     *    AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]
     *  Researchers who wish to cite entries in their publications
     *  should always cite the first accession number. This is
     *  commonly referred to as the 'primary accession
     *  number'."
     *  
     *  Examples:
     *   AC   P43332; Q9W4D7;
     *   AC   X06347;
     *   
     *  Note that Easel only stores primary accessions.
     *  Because there can be more than one accession line, we check to 
     *  see if the accession is already set before storing a line.
     */
    if (strncmp(sqfp->buf, "AC   ", 5) == 0 && sq->acc[0] == '\0')
      {
        s = sqfp->buf+5;
        if ( esl_strtok(&s, ";", &tok) != eslOK)
          ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to parse accession on AC line", sqfp->linenumber);
        if ((status = esl_sq_SetAccession(sq, tok)) != eslOK) return status;
      }

    /* "The format of the DE line is:
     *    DE   Description.
     * ...In cases where more than one DE line is required, the text is
     * only divided between words and only the last DE line is
     * terminated by a period."
     * 
     * Examples:
     *   DE   U1 small nuclear ribonucleoprotein A (U1 snRNP protein A) (U1-A) (Sex
     *   DE   determination protein snf).
     *   
     *   DE   Human mRNA for U1 small nuclear RNP-specific A protein
     *
     *   DE   RecName: Full=U1 small nuclear ribonucleoprotein A;
     *   DE            Short=U1 snRNP protein A;
     *   DE            Short=U1-A;
     *   DE   AltName: Full=Sex determination protein snf;
     *
     * We'll make no attempt to parse the structured UniProt description header,
     * for the moment.
     */
    if (strncmp(sqfp->buf, "DE   ", 5) == 0)
      {
        s = sqfp->buf+5;
        esl_strchop(s, sqfp->nc-5);
        if ((status = esl_sq_AppendDesc(sq, s)) != eslOK)
          ESL_FAIL(status, sqfp->errbuf, "Line %" PRId64 ": failed to parse description on DE line", sqfp->linenumber);
      }

    /* UniProt: "The format of the SQ line is:
     *  SQ   SEQUENCE XXXX AA; XXXXX MW; XXXXXXXXXXXXXXXX CRC64;"
     * EMBL:    "The SQ (SeQuence header) line marks the beginning of 
     *           the sequence data and Gives a summary of its content. 
     *           An example is:
     *  SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;"
     *  
     * We don't parse this line; we just look for it as the last line
     * before the sequence starts.
     */
  } while (strncmp(sqfp->buf, "SQ   ", 5) != 0);

  if (loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find any sequence");
  sq->hoff = sqfp->boff - 1;
  sq->doff = sqfp->boff;
  return eslOK;
}

static int
end_embl(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int status;

  if (strncmp(sqfp->buf, "//", 2) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", sqfp->linenumber);
  sq->eoff = sqfp->boff + sqfp->nc - 1;
  status = loadbuf(sqfp);
  if      (status == eslEOF) return eslOK; /* ok, actually. */
  else if (status == eslOK)  return eslOK;
  else                       return status;
}
/*-------- end, internal functions (EMBL|UniProt|TrEMBL)---------*/



/*****************************************************************
 * 11. Internal functions (GenBank) 
 *****************************************************************/ 
/* NCBI GenBank sequence database format.
 * See GenBank release notes; for example,
 * ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
 */

static void
config_genbank(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased      = TRUE;
  sqfp->eof_is_ok         = FALSE;/* records end with //  */
  sqfp->parse_header      = &header_genbank;
  sqfp->parse_end         = &end_genbank;
}

static void
inmap_genbank(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
    sqfp->inmap['-']  = eslDSQ_ILLEGAL;                      // don't let the abc_inmap map the gap char; this is an ungapped file format
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  for (x = '0'; x <= '9'; x++)
    sqfp->inmap[x] = eslDSQ_IGNORED;
  sqfp->inmap['*']  = '*';         /* accept * as a nonresidue/stop codon character */
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;/* DOS eol compatibility */
  sqfp->inmap['/']  = eslDSQ_EOD;
}

/* header_genbank()
 * 
 * sqfp->buf is the first (LOCUS) line of the entry, or a line before
 * it (in which case we'll scan forwards to find the LOCUS line - even
 * skipping non-blank lines, because there are sometimes headers at
 * the start of GenBank files).
 * 
 * On success, returns <eslOK> and:
 *   sq->name  contains sequence name (and may have been reallocated, changing sq->nalloc)
 *   sq->acc   contains seq accession (and may have been reallocated, changing sq->aalloc)
 *   sq->desc  contains description line (and may have been reallocated, changing sq->dalloc)
 *   sq->roff  has been set to the record offset
 *   sq->doff  has been set to the data offset (start of sequence line)
 *   sqfp->buf is the first seq line.
 * 
 * If no more seqs are found in the file, returns <eslEOF>.
 * On parse failure, returns <eslEFORMAT>, leaves as mesg in sqfp->errbuf.
 * 
 * May also throw <eslEMEM> on allocation errors.
 */
static int
header_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  char *s;
  char *tok;
  int   status;

  /* Find LOCUS line, ignoring a file header.  */
  if (sqfp->nc == 0) return eslEOF;
  while (strncmp(sqfp->buf, "LOCUS   ", 8) != 0) {
    if ((status = loadbuf(sqfp)) == eslEOF) return eslEOF; /* normal   */
    else if (status != eslOK) return status;                /* abnormal */
  } 
  
  s = sqfp->buf+12;
  if ( esl_strtok(&s, " ", &tok) != eslOK)
    ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to parse name on LOCUS line", sqfp->linenumber);
  if ((status = esl_sq_SetName(sq, tok)) != eslOK) return status;
  sq->roff = sqfp->boff;/* record the disk offset to the LOCUS line */
  
  /* Look for ORIGIN line, parsing optional info as we go. */
  do {
    if ( loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find ORIGIN line");

    /* Optional VERSION line is parsed as "accession". */
    if (strncmp(sqfp->buf, "VERSION   ", 10) == 0)
    {
      s = sqfp->buf+12;
      if ( esl_strtok(&s, " \t\n", &tok) != eslOK)
        ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": failed to parse VERSION line", sqfp->linenumber);
      if ((status = esl_sq_SetAccession(sq, tok)) != eslOK) return status;
    }

    /* Optional DEFINITION Line is parsed as "description". */
    if (strncmp(sqfp->buf, "DEFINITION ", 11) == 0)
    {
      s = sqfp->buf+12;
      esl_strchop(s, sqfp->nc-12);
      if ((status = esl_sq_AppendDesc(sq, s)) != eslOK)
        ESL_FAIL(status, sqfp->errbuf, "Line %" PRId64 ": failed to parse desc on DEFINITION line", sqfp->linenumber);
    }
  } while (strncmp(sqfp->buf, "ORIGIN", 6) != 0);

  if (loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find any sequence");
  sq->hoff = sqfp->boff - 1;
  sq->doff = sqfp->boff;
  return eslOK;
}
  
static int
end_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int status;

  if (strncmp(sqfp->buf, "//", 2) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", sqfp->linenumber);
  sq->eoff = sqfp->boff + sqfp->nc - 1;
  status = loadbuf(sqfp);
  if      (status == eslEOF) return eslOK; /* ok, actually; we'll detect EOF on next sq read */
  else if (status == eslOK)  return eslOK;
  else                       return status;
}
/*------------- end, internal functions (GenBank) ----------------*/




/*****************************************************************
 * 12. Internal functions (daemon)
 *****************************************************************/

/* Special case FASTA format where each sequence is terminated with "//".
 * 
 * The use case is where the sequences are being read from a pipe and a
 * way is needed to signal the end of the sequence so it can be processed.
 * The next sequence might not be in the pipe, so the usual '>' is not
 * present to signal the end of the sequence.  Also, an EOF is not
 * an option, since the daemon might run continuously.
 */

static void
config_daemon(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased = FALSE;
  sqfp->eof_is_ok    = FALSE;
  sqfp->parse_header = &header_fasta;
  sqfp->parse_end    = &end_daemon;
}

static void
inmap_daemon(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
    sqfp->inmap['-']  = eslDSQ_ILLEGAL;                      // don't let the abc_inmap map the gap char; this is an ungapped file format
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  sqfp->inmap['*']  = '*';         /* accept * as a nonresidue/stop codon character */
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;/* DOS eol compatibility */
  sqfp->inmap['\n'] = eslDSQ_EOL;
  sqfp->inmap['/']  = eslDSQ_EOD;
  /* \n is special - fasta reader detects it as an eol */
}


/* end_daemon()
 * 
 * Special case FASTA format where each sequence is terminated with "//".
 * 
 * The use case is were the sequences are being read from a pipe and a
 * way is needed to signal the end of the sequence so it can be processed.
 */
static int 
end_daemon(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  ESL_UNUSED(sq);
  char  c;

  if (sqfp->nc < 3) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Whoops, daemon input stream is corrupted");

  c =  sqfp->buf[sqfp->bpos++];
  if (c != '/') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", sqfp->linenumber);

  c =  sqfp->buf[sqfp->bpos++];
  if (c != '/') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", sqfp->linenumber);

  /* skip to end of line */
  while (c != '\n' && c != '\r' && sqfp->bpos < sqfp->nc) c =  sqfp->buf[sqfp->bpos++];

  /* skip past end of line */
  while ((c == '\n' || c == '\r') && sqfp->bpos < sqfp->nc) c =  sqfp->buf[sqfp->bpos++];

  return eslOK;
}
/*---------------- end, internal functions (daemon) -------------*/


/*****************************************************************
 * 13. Internal functions (hmmpgmd)
 *****************************************************************/

static int
fileheader_hmmpgmd(ESL_SQFILE *sqfp)
{
  char c;
  int  status = eslOK;

  /* We've just loaded first buffer, after an Open. First char should be the # of the hmmpgmd file,
   * but let's tolerate leading whitespace anyway
   */
  c =  sqfp->buf[sqfp->bpos];
  while (status == eslOK && isspace(c)) status = nextchar(sqfp, &c); /* skip space (including \n, \r) */
  if (status == eslEOF) return eslEOF;

  if (c != '#') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "hmmpgmd file expected to start with #");

  /* skip first line; remainder of file is FASTA format */
  while (status == eslOK && (c != '\n' && c != '\r')) status = nextchar(sqfp, &c); 
  if (status == eslEOF) return eslEOF;

  /* next character read should be the '>' of the first FASTA record. We're properly positioned at "start of file". */
  return eslOK;
}
/*------------- end, internal functions (hmmpgmd) ---------------*/




/*****************************************************************
 * 14. Benchmark driver
 *****************************************************************/ 
#ifdef eslSQIO_BENCHMARK
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-d",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use digital sequence input mode",                  0 },
  { "-i",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "benchmark ReadInfo() input",                       0 },
  { "-w",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "benchmark ReadWindow() input",                     0 },
  { "-B",        eslARG_INT,   "4096",  NULL, NULL,  NULL,  NULL, NULL, "buffer size for read, fread tests",                0 },
  { "-C",        eslARG_INT,    "100",  NULL, NULL,  NULL,  NULL, NULL, "context size for ReadWindow()",                    0 },
  { "-W",        eslARG_INT,   "1000",  NULL, NULL,  NULL,  NULL, NULL, "window size for ReadWindow()",                     0 },
  { "-2",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  "-w", NULL, "with ReadWindow(), do both strands",               0 },
  { "--format",  eslARG_STRING,  NULL,  NULL, NULL,  NULL,  NULL, NULL, "assert <seqfile> is in format <s>",                0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet, not DNA",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <DNA FASTA file>";
static char banner[] = "benchmark driver for sqio module";

static int benchmark_read (char *filename, int bufsize, int64_t *ret_magic);
static int benchmark_fread(char *filename, int bufsize, int64_t *ret_magic);
static int benchmark_fgets(char *filename, int bufsize, int64_t *ret_magic);
/*static int benchmark_mmap (char *filename, int bufsize, int64_t *ret_magic);*/

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go       = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w        = esl_stopwatch_Create();
  ESL_ALPHABET  *abc      = NULL;
  ESL_SQ        *sq       = NULL;
  ESL_SQFILE    *sqfp     = NULL;
  char          *filename = esl_opt_GetArg(go, 1);
  int            format   = eslSQFILE_UNKNOWN;
  int            bufsize  = esl_opt_GetInteger(go, "-B");
  int            C        = esl_opt_GetInteger(go, "-C");
  int            W        = esl_opt_GetInteger(go, "-W");
  int            do_crick = esl_opt_GetBoolean(go, "-2");
  int            n        = 0;
  int64_t        magic    = 0;
  int64_t        nr       = 0;

  if (esl_opt_IsOn(go, "--format")) {
    format = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--format"));
    if (format == eslSQFILE_UNKNOWN) esl_fatal("unrecognized database format %s\n", esl_opt_GetString(go, "--format"));
  }

  if (esl_opt_GetBoolean(go, "-d"))
    {
      if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);
      else                                   abc = esl_alphabet_Create(eslDNA);
      sq = esl_sq_CreateDigital(abc);
      if (esl_sqfile_OpenDigital(abc, filename, format, NULL, &sqfp) != eslOK) esl_fatal("failed to open %s", filename);
    } 
  else 
    {
      sq = esl_sq_Create();
      if (esl_sqfile_Open(filename, format, NULL, &sqfp) != eslOK) esl_fatal("failed to open %s", filename);
    }


  /* It's useful to have some baselines of just reading the file without parsing;
   * POSIX read(); C fread(); C fgets(); and POSIX mmap(). 
   * Counting a's (<na>) is just to keep the optimizer from optimizing the 
   * benchmark away.
   */
  esl_stopwatch_Start(w);   benchmark_read (filename, bufsize, &magic);   esl_stopwatch_Stop(w);  printf("magic=%" PRId64 "; ", magic); esl_stopwatch_Display(stdout, w, "read():  "); 
  esl_stopwatch_Start(w);   benchmark_fread(filename, bufsize, &magic);   esl_stopwatch_Stop(w);  printf("magic=%" PRId64 "; ", magic); esl_stopwatch_Display(stdout, w, "fread(): ");
  esl_stopwatch_Start(w);   benchmark_fgets(filename, bufsize, &magic);   esl_stopwatch_Stop(w);  printf("magic=%" PRId64 "; ", magic); esl_stopwatch_Display(stdout, w, "fgets(): ");
  /*  esl_stopwatch_Start(w);   benchmark_mmap (filename, bufsize, &magic);   esl_stopwatch_Stop(w);  printf("magic=%" PRId64 "; ", magic); esl_stopwatch_Display(stdout, w, "mmap():  "); */

  esl_stopwatch_Start(w);
  if (esl_opt_GetBoolean(go, "-i"))
    {
      while (esl_sqio_ReadInfo(sqfp, sq) == eslOK) { n++; nr += sq->L; esl_sq_Reuse(sq); }
    }
  else if (esl_opt_GetBoolean(go, "-w"))
    {
      int wstatus;
      while ((wstatus = esl_sqio_ReadWindow(sqfp, C, W, sq)) != eslEOF)
	{ 
	  if        (wstatus == eslEOD) {
	    if (!do_crick || W < 0) { n++; esl_sq_Reuse(sq); }
	    if (do_crick)           { W = -W; }
	    continue;
	  } else if (wstatus != eslOK) esl_fatal("Error: %s", esl_sqfile_GetErrorBuf(sqfp));
	  nr += sq->W;
	}
    }
  else 
    {
      while (esl_sqio_Read(sqfp, sq) == eslOK)  { n++; nr += sq->L; esl_sq_Reuse(sq); }
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "sqio:  ");
  printf("Read %d sequences; %lld residues.\n", n, (long long int) nr);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

static int
benchmark_read(char *filename, int bufsize, int64_t *ret_magic)
{
  char *buf  = malloc(sizeof(char) * bufsize);
  int   fd   = open(filename, O_RDONLY);
  int   n;
  int64_t magic = 0;
  int   pos;

  if (fd != -1) 
    {
      while ( (n = read(fd, buf, bufsize)) > 0)
	{
	  for (pos = 0; pos < n; pos++) magic += buf[pos];
	}
      close(fd);
    }

  free(buf);
  *ret_magic = magic;
  return eslOK;
}
  
static int
benchmark_fread(char *filename, int bufsize, int64_t *ret_magic)
{
  FILE *fp   = fopen(filename, "r");
  char *buf  = malloc(sizeof(char) * bufsize);
  int64_t magic = 0;
  int   pos;
  int   n;
  
  if (fp != NULL)
    {
      while ( (n = fread(buf, sizeof(char), bufsize, fp)) > 0)
	{
	  for (pos = 0; pos < n; pos++) magic += buf[pos];
	} 
      fclose(fp);
    }

  free(buf);
  *ret_magic = magic;
  return eslOK;
}


static int
benchmark_fgets(char *filename, int bufsize, int64_t *ret_magic)
{
  FILE *fp   = fopen(filename, "r");
  char *buf  = malloc(sizeof(char) * bufsize);
  int64_t magic = 0;
  char *s;
  
  if (fp != NULL)
    {
      while ( fgets(buf, bufsize, fp) != NULL)
	{
	  for (s = buf; *s; s++) magic += *s;
	} 
      fclose(fp);
    }

  free(buf);
  *ret_magic = magic;
  return eslOK;
}


#if 0
/* we can't use mmap anyway; we have to be able to read from streams
 * and pipes.
 */
static int
benchmark_mmap(char *filename, int bufsize, int64_t *ret_magic)
{
  int   fd   = open(filename, O_RDONLY);
  struct stat statbuf;
  char *p;
  uint64_t pos;
  uint64_t magic = 0;
 
  fstat(fd, &statbuf);
  p = mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, fd, 0);
  for (pos = 0; pos < statbuf.st_size; pos++)
    magic += p[pos];

  close(fd);
  *ret_magic = magic;
  return eslOK;
}
#endif  


#endif /*eslSQIO_BENCHMARK*/
/*-------------- end, benchmark driver --------------------------*/



/*****************************************************************
 * 15. Unit tests
 *****************************************************************/ 
#ifdef eslSQIO_TESTDRIVE

#include "esl_keyhash.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

static void
synthesize_testseqs(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int maxL, int N, ESL_SQ ***ret_sqarr)
{
  ESL_SQ **sqarr  = NULL;
  ESL_KEYHASH *kh = esl_keyhash_Create();
  float   *fq     = NULL;
  char    *buf    = NULL;
  int      maxn   = eslSQ_NAMECHUNK*2;
  int      maxa   = eslSQ_ACCCHUNK*2;
  int      maxd   = eslSQ_DESCCHUNK*2;
  int      n      = ESL_MAX( maxn, ESL_MAX(maxa, maxd));
  char     ascii[128];
  float    af[128];
  int      i, pos;
  int      x;

  if ((sqarr = malloc(sizeof(ESL_SQ *) * N))       == NULL) esl_fatal("malloc failed");
  if ((fq    = malloc(sizeof(float)    * abc->Kp)) == NULL) esl_fatal("malloc failed");
  if ((buf   = malloc(sizeof(char)     * (n+1)))   == NULL) esl_fatal("malloc failed");

  /* Set a residue frequency vector that's going to sample degenerate residues too */
  esl_vec_FSet(fq, abc->Kp, 0.0);
  esl_vec_FSet(fq, abc->K,  0.9 / (float) abc->K);
  esl_vec_FSet(fq + abc->K + 1, abc->Kp - abc->K - 2, 0.1 / (float) (abc->Kp - abc->K - 2));

  /* Set an ASCII frequency vector that'll sample all nonspace chars */
  for (x = 0; x < 128; x++) {
    ascii[x] = (char) x;
    if      (isalpha(x))             af[x] = 3.0;
    else if (isdigit(x))             af[x] = 2.0;
    else if (ispunct(x) && x != '%') af[x] = 1.0; /* disallow %; it'll screw up printf()-like Set calls */
    else                             af[x] = 0.0;
  }
  esl_vec_FNorm(af, 128);

  for (i = 0; i < N; i++)
    {
      if ((sqarr[i] = esl_sq_CreateDigital(abc)) == NULL) esl_fatal("failed to allocate seq %d", i);

      do {
	n = esl_rnd_Roll(r, maxn) + 1; /* 1..maxn */
	esl_rsq_fIID(r, ascii, af, 128, n, buf);
	buf[n] = '\0';
      }	while (ispunct(buf[0]) ||                                // #, // are bad things for names to start with, in Stockholm format 
               esl_keyhash_Store(kh, buf, n, NULL) == eslEDUP);  // Make sure names are unique.
      esl_sq_SetName(sqarr[i], buf);

      if (esl_rnd_Roll(r, 2) == 0) { /* 50% chance of an accession */
	n = esl_rnd_Roll(r, maxa) + 1; 
	esl_rsq_fIID(r, ascii, af, 128, n, buf);
	buf[n] = '\0';
	esl_sq_SetAccession(sqarr[i], buf);
      }

      if (esl_rnd_Roll(r, 2) == 0) { /* 50% chance of a description */
	n = esl_rnd_Roll(r, maxd) + 1;
	esl_rsq_fIID(r, ascii, af, 128, n, buf);
	buf[n] = '\0';
	for (pos = 1; pos < n-1; pos++) {                 /* avoid first, last char, and... */
	  if (esl_rnd_Roll(r, 10)  == 0) buf[pos] = ' ';  /* ...sprinkle with spaces... */
	  if (esl_rnd_Roll(r, 100) == 0) buf[pos] = '\t'; /* ...and tabs. */
	}
	esl_sq_SetDesc(sqarr[i], buf);
      }

      n = esl_rnd_Roll(r, (maxL+1)); /* choose seqlen =  0..maxL; 0 length seqs occur in real dbs */
      esl_sq_GrowTo(sqarr[i], n);
      esl_rsq_xfIID(r, fq, abc->Kp, n, sqarr[i]->dsq);

      esl_sq_SetCoordComplete(sqarr[i], n);
    }

  *ret_sqarr = sqarr;
  free(buf);
  free(fq);
  esl_keyhash_Destroy(kh);
  return;
}

/* Write an uglified FASTA file to a stream.
 * Also, remember where the start of the descline and first
 * seq line are, in sq->{roff,doff}. We'll compare against
 * what the input function thinks these locations are.
 */
static void
write_ugly_fasta(ESL_RANDOMNESS *r, FILE *fp, ESL_SQ *sq)
{
  char buf[61];
  int  pos;
  
  sq->roff = ftello(fp);
  fputc('>', fp);
  while (esl_rnd_Roll(r, 10) == 0) fputc(' ', fp);
  fprintf(fp, "%s", sq->name);
  while (esl_rnd_Roll(r, 10) == 0) fputc(' ', fp);
  if (sq->desc[0] != 0) fprintf(fp, " %s", sq->desc);
  fputc('\n', fp);

  sq->doff = ftello(fp);             
  buf[60] = '\0';
  for (pos = 1; pos <= sq->n; pos+=60)
    {
      while (esl_rnd_Roll(r, 10) == 0) fputc(' ', fp);
      esl_dsq_TextizeN(sq->abc, sq->dsq+pos, 60, buf);
      fputs(buf, fp);
      fputc('\n', fp);
    }
  while (esl_rnd_Roll(r, 10) == 0) fputc('\n', fp);

  sq->eoff = ftello(fp) - 1;
  if (sq->n == 0) sq->doff = sq->eoff+1;  // Deals with an edge case, an L=0 seq with multiple newlines.
}

static void
write_spaced_fasta(FILE *fp, ESL_SQ *sq)
{
  char buf[64];
  int  pos;

  sq->roff = ftello(fp);
  fprintf(fp, ">%s", sq->name);
  if (sq->desc[0] != 0) fprintf(fp, " %s", sq->desc);
  fputc('\n', fp);

  sq->doff = ftello(fp);
  buf[10]  = '\0';
  for (pos = 1; pos <= sq->n; pos += 10)
    {
      esl_dsq_TextizeN(sq->abc, sq->dsq+pos, 10, buf);
      fputs(buf, fp);
      if (pos+9 >= sq->n || (pos+9) % 60 == 0) fputc('\n',  fp);
      else                                     fputc(' ', fp);
    }
  sq->eoff = ftello(fp) - 1;
}

/* write_mixed_fasta()
 *
 * Write a FASTA file with deliberately inconsistent line lengths: each
 * record is either 60 residues per line, or all on one line. This tests our
 * detection of "well-formatted" files for bpl/rpl subseq indexing.
 */
static void
write_mixed_fasta(ESL_RANDOMNESS *r, FILE *fp, ESL_SQ *sq)
{
  char *buf = NULL;
  int   w;
  int   pos;

  if ((buf = malloc(sizeof(char) * (sq->n+1))) == NULL) esl_fatal("malloc failed");
  esl_dsq_TextizeN(sq->abc, sq->dsq+1, sq->n, buf);
  buf[sq->n] = '\0';

  w = (esl_rnd_Roll(r, 2) == 0) ? 60 : ESL_MAX(1, sq->n);  // 60 residues per line, or all on one line 

  sq->roff = ftello(fp);
  fprintf(fp, ">%s", sq->name);
  if (sq->desc[0] != 0) fprintf(fp, " %s", sq->desc);
  sq->hoff = ftello(fp);    // unlike the writers above, we have to set hoff: sqarr[] offsets are stale from mode 2 by now 
  fputc('\n', fp);

  sq->doff = ftello(fp);
  for (pos = 0; pos < sq->n; pos += w)
    fprintf(fp, "%.*s\n", (int) ESL_MIN(w, sq->n-pos), buf+pos);
  sq->eoff = ftello(fp) - 1;
  free(buf);
}

static void
check_geometry(int mode, int bpl, int rpl, char *msg)
{
  switch (mode) {
  case 0:  if (bpl != 0)               esl_fatal(msg); break; // uglified: bpl should be invalid (rpl might not be) 
  case 1:  if (rpl != 60 || bpl == 0)  esl_fatal(msg); break; // spaced: bpl, rpl should be valid 
  case 2:  if (rpl != 60 || bpl != 61) esl_fatal(msg); break; // normal: bpl, rpl should be valid, w/ bpl=rpl+1 
  case 3:  break;   // mixed: bpl,rpl could come out either way here; utest_subseq_geometry() pins down the exact cases
  }
}

static void
make_ssi_index(ESL_ALPHABET *abc, const char *tmpfile, int format, char *ssifile, int mode)
{ 
  char       *msg  = "sqio unit testing: failed to make SSI index";
  ESL_NEWSSI *ns   = NULL;
  ESL_SQFILE *sqfp = NULL;
  ESL_SQ     *sq   = esl_sq_CreateDigital(abc);
  uint16_t    fh   = 0;
  int         status;

  int         bpl, rpl;
 
  snprintf(ssifile, 32, "%s.ssi", tmpfile);  // 32 is the allocation size of ssifile, from main()
  if (esl_newssi_Open(ssifile, TRUE, &ns)                       != eslOK) esl_fatal(msg);
  if (esl_newssi_AddFile(ns, tmpfile, format, &fh)              != eslOK) esl_fatal(msg);
  if (esl_sqfile_OpenDigital(abc, tmpfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      if (esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->L)   != eslOK) esl_fatal(msg);
      if (sq->acc[0] != '\0' && esl_newssi_AddAlias(ns, sq->acc, sq->name) != eslOK) esl_fatal(msg);
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal(msg);
  
  bpl = sqfp->bpl;
  rpl = sqfp->rpl;
  if (bpl > 0 && rpl > 0) 
    if ((status = esl_newssi_SetSubseq(ns, fh, bpl, rpl)) != eslOK) esl_fatal(msg);
  
  if (esl_newssi_Write(ns) != eslOK)  esl_fatal(msg);

  check_geometry(mode, sqfp->bpl, sqfp->rpl, msg);

  esl_sqfile_Close(sqfp);
  esl_newssi_Close(ns);
  esl_sq_Destroy(sq);
}

static void
utest_read(ESL_ALPHABET *abc, ESL_SQ **sqarr, int N, char *seqfile, int format, int mode)
{
  char       *msg         = "sqio complete read unit test failed";
  ESL_SQ     *sq          = esl_sq_CreateDigital(abc);
  ESL_SQFILE *sqfp        = NULL;
  int         nseq        = 0;
  int         status;
 
  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      /* FASTA doesn't preserve accessions. Copy it, as a hack, so Compare test succeeds*/
      if (sq->acc[0] == '\0' && esl_sq_SetAccession(sq, sqarr[nseq]->acc) != eslOK) esl_fatal(msg);
      if (esl_sq_Compare(sq, sqarr[nseq])                                 != eslOK) esl_fatal(msg);      
      nseq++;
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal(msg);
  if (nseq   != N)      esl_fatal(msg);

  check_geometry(mode, sqfp->bpl, sqfp->rpl, msg);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
}

static void
utest_read_info(ESL_ALPHABET *abc, ESL_SQ **sqarr, int N, char *seqfile, int format, int mode)
{
  char       *msg         = "sqio info read unit test failed";
  ESL_SQ     *sq          = esl_sq_CreateDigital(abc);
  ESL_SQFILE *sqfp        = NULL;
  int         nseq        = 0;
  int         status;
 
  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      if (strcmp(sq->name,   sqarr[nseq]->name)   != 0) esl_fatal(msg);
      if (format != eslSQFILE_FASTA && 
	  strcmp(sq->acc,    sqarr[nseq]->acc)    != 0) esl_fatal(msg);
      if (strcmp(sq->desc,   sqarr[nseq]->desc)   != 0) esl_fatal(msg);
      if (strcmp(sq->source, sqarr[nseq]->source) != 0) esl_fatal(msg);
      if (sq->n     != 0)  esl_fatal(msg);
      if (sq->start != 0)  esl_fatal(msg);
      if (sq->end   != 0)  esl_fatal(msg);
      if (sq->C     != 0)  esl_fatal(msg);
      if (sq->W     != 0)  esl_fatal(msg);
      if (sq->L     != sqarr[nseq]->L)                  esl_fatal(msg);
      if (sq->roff != -1 && sqarr[nseq]->roff != -1 && sq->roff != sqarr[nseq]->roff) esl_fatal(msg);
      if (sq->doff != -1 && sqarr[nseq]->doff != -1 && sq->doff != sqarr[nseq]->doff) esl_fatal(msg);
  
      nseq++;
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal(msg);
  if (nseq   != N)      esl_fatal(msg);

  check_geometry(mode, sqfp->bpl, sqfp->rpl, msg);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
}

static void
utest_read_window(ESL_ALPHABET *abc, ESL_SQ **sqarr, int N, char *seqfile, int format, int mode)
{
  char       *msg         = "sqio window read unit test failed";
  ESL_SQ     *sq          = esl_sq_CreateDigital(abc);
  ESL_SQ     *rev         = esl_sq_CreateDigital(abc);
  ESL_SQFILE *sqfp        = NULL;
  int         nseq        = 0;
  int         C           = 10;
  int         W           = 50;
  int         nres        = 0;
  int         L;
  int         wstatus;

 
  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  while ((wstatus = esl_sqio_ReadWindow(sqfp, C, W, sq)) == eslOK || wstatus == eslEOD)
    {
      if (wstatus == eslEOD) {
	if (W < 0) {
	  nseq++; 
	  nres = 0;
	  W    = -W;
	  esl_sq_Reuse(sq); 
	  esl_sq_Reuse(rev);
	} else       {
	  /* reverse complement */
	  nres = 0;
	  W    = -W; 	
	  esl_sq_Copy(sqarr[nseq], rev);
	  esl_sq_ReverseComplement(rev);
	}
	continue;
      }

      nres += sq->W;
      if (strcmp(sq->name,   sqarr[nseq]->name)   != 0) esl_fatal(msg);
      if (format != eslSQFILE_FASTA && 
	  strcmp(sq->acc,    sqarr[nseq]->acc)    != 0) esl_fatal(msg);
      if (strcmp(sq->desc,   sqarr[nseq]->desc)   != 0) esl_fatal(msg);

      L = sqfp->L;

      if (W > 0) {
	/* Forward strand coord checks */
	if (L   != nres)                                esl_fatal(msg);
	if (sq->start != nres - sq->n + 1)              esl_fatal(msg);
	if (sq->end   != nres)                          esl_fatal(msg);
	if (sq->C != 0 && sq->C != C)                   esl_fatal(msg);
	if (sq->n != sq->C+sq->W)                       esl_fatal(msg);
	if (sq->start+sq->n-1 > sqarr[nseq]->L)         esl_fatal(msg);
	if (wstatus == eslEOD && sq->L != L)            esl_fatal(msg);
	if (memcmp(sq->dsq + 1, sqarr[nseq]->dsq + sq->start, sq->C+sq->W) != 0) esl_fatal(msg);
      } else {
	/* Reverse strand coord checks */
	if (L    != -1)                                 esl_fatal(msg);
	if (sq->start  != sq->L - nres + sq->W + sq->C) esl_fatal(msg);
	if (sq->end    != sq->L - nres + 1)             esl_fatal(msg);
	if (sq->C != 0 && sq->C != C)                   esl_fatal(msg);
	if (sq->start-sq->n+1 < 1)                      esl_fatal(msg);
	if (wstatus == eslEOD && sq->end != 1)          esl_fatal(msg);
	if (memcmp(sq->dsq + 1, rev->dsq + (sq->L - sq->start + 1), sq->C+sq->W) != 0) esl_fatal(msg);
      }
    }

  check_geometry(mode, sqfp->bpl, sqfp->rpl, msg);

  if (wstatus != eslEOF) esl_fatal(msg);
  if (nseq    != N)      esl_fatal(msg);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(rev);
  esl_sq_Destroy(sq);
}

static void
utest_fetch_subseq(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, ESL_SQ **sqarr, int N, char *seqfile, char *ssifile, int format)
{
  char       *msg         = "sqio subseq read unit test failure";
  ESL_SQ     *sq          = esl_sq_CreateDigital(abc);
  ESL_SQFILE *sqfp        = NULL;
  int         i;
  int         ntest       = 32;
  char       *source;
  int         start;
  int         end;

  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  if (esl_sqfile_OpenSSI(sqfp, ssifile)                         != eslOK) esl_fatal(msg);
  while (ntest--) 
    {
      i      = esl_rnd_Roll(r, N);
      source = sqarr[i]->name; 
      if (sqarr[i]->L == 0) continue;   // Don't try to fetch from empty sequences.
      
      do {
	start = esl_rnd_Roll(r, sqarr[i]->n) + 1;
	end   = esl_rnd_Roll(r, sqarr[i]->n) + 1;
      } while (start > end);

      if (esl_sqio_FetchSubseq(sqfp, source, start, end, sq)        != eslOK) esl_fatal(msg);
      if (memcmp(&(sqarr[i]->dsq[start]), &sq->dsq[1], end-start+1) != 0)     esl_fatal(msg);
      
      esl_sq_Reuse(sq);
    }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
}


/* utest_subseq_geometry()
 *
 * A file only qualifies for bpl/rpl subseq indexing if every sequence
 * line in it is bpl,rpl - except that the last line of a record may be
 * shorter, and (when bpl=rpl+1, single residue resolution) no line may
 * hold anything but residues and its newline. Check that we detect
 * violations wherever they hide: in a record with only one sequence
 * line (which has no nonterminal line to be compared against); in the
 * first record, before bpl,rpl have been set at all; and on a final
 * line with no terminating newline. Fetches from such a file are
 * silently corrupted if we get this wrong.
 */
static void
utest_subseq_geometry(void)
{
  char *msg = "sqio subseq geometry unit test failed";
  struct { char *fasta; int bpl; int rpl; } trial[] = {
    /* uniform 10 residues/line, shorter last line in each record: valid, residue resolution */
    { ">a\nACGTACGTAC\nACGTACGTAC\nAC\n>b\nACGTACGTAC\nACGT\n",               11, 10 },
    /* a record on one line, longer than rpl */
    { ">a\nACGTACGTAC\nACGTACGTAC\nAC\n>b\nACGTACGTACACGTACGTACACGT\n",        0,  0 },
    /* ... and the same, before bpl,rpl have been set by any other record */
    { ">b\nACGTACGTACACGTACGTACACGT\n>a\nACGTACGTAC\nACGTACGTAC\nAC\n",        0,  0 },
    /* a final line longer than bpl, with no terminating newline */
    { ">a\nACGTACGTAC\nACGTACGTAC\nAC\n>b\nACGTACGTAC\nACGTACGTACACGT",        0,  0 },
    /* a space on a last line, which makes bpl=rpl+1 misleading; rpl stays valid, bpl doesn't */
    { ">a\nACGTACGTAC\nACGTACGTAC\nAC\n>b\nACGTACGTAC\nAC GT\n",               0, 10 },
    /* ... and the same, before bpl,rpl have been set */
    { ">b\nAC GT\n>a\nACGTACGTAC\nACGTACGTAC\nAC\n",                           0, 10 },
    /* ... and the same, with no terminating newline */
    { ">a\nACGTACGTAC\nACGTACGTAC\nAC\n>b\nACGTACGTAC\nAC GT",                 0, 10 },
  };
  int         ntrials = sizeof(trial) / sizeof(trial[0]);
  char        tmpfile[32];
  FILE       *fp   = NULL;
  ESL_SQFILE *sqfp = NULL;
  ESL_SQ     *sq   = esl_sq_Create();
  int         i;
  int         status;

  for (i = 0; i < ntrials; i++)
    {
      strcpy(tmpfile, "esltmpXXXXXX");
      if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("failed to make tmpfile");
      fputs(trial[i].fasta, fp);
      fclose(fp);

      if (esl_sqfile_Open(tmpfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal(msg);
      while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK) esl_sq_Reuse(sq);
      if (status != eslEOF)                     esl_fatal("%s: trial %d", msg, i);
      if (sqfp->bpl != trial[i].bpl) esl_fatal("%s: trial %d, bpl %d not %d", msg, i, sqfp->bpl, trial[i].bpl);
      if (sqfp->rpl != trial[i].rpl) esl_fatal("%s: trial %d, rpl %d not %d", msg, i, sqfp->rpl, trial[i].rpl);

      esl_sqfile_Close(sqfp);
      remove(tmpfile);
    }
  esl_sq_Destroy(sq);
}


/* Write the sequences out to a tmpfile in chosen <format>;
 * read them back and make sure they're the same.
 * reposition to beginning, read and check again.
 *
 * The sequences in <sqarr> are in digital mode.
 */
static void
utest_write(ESL_ALPHABET *abc, ESL_SQ **sqarr, int N, int format)
{
  char       *msg         = "sqio write unit test failure";
  char        tmpfile[32] = "esltmpXXXXXX";
  ESL_SQFILE *sqfp        = NULL;
  ESL_SQ     *sq          = esl_sq_CreateDigital(abc);
  FILE       *fp          = NULL;
  int         iterations  = 2;	/* 2: reposition and read again */
  int         i;
  int         require_nonzero_length = FALSE;

  if (esl_sqio_IsAlignment(format)) require_nonzero_length = TRUE;

  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal(msg);
  for (i = 0; i < N; i++)
    if (! require_nonzero_length || sqarr[i]->L > 0)
      esl_sqio_Write(fp, sqarr[i], format, FALSE);
  fclose(fp);

  if (esl_sqfile_OpenDigital(abc, tmpfile, format, NULL, &sqfp)           != eslOK)  esl_fatal(msg);
  while (iterations--)
    {
      for (i = 0; i < N; i++)
	{
          if (require_nonzero_length && sqarr[i]->L == 0) continue;
	  if (esl_sqio_Read(sqfp, sq)                                     != eslOK)  esl_fatal(msg);
	  if (strcmp(sqarr[i]->name,   sq->name)                          != 0)      esl_fatal(msg);
	  if (sqarr[i]->L                                                 !=  sq->L) esl_fatal(msg);
	  if (memcmp(sqarr[i]->dsq, sq->dsq, sizeof(ESL_DSQ) * (sq->L+2)) != 0)      esl_fatal(msg);
	  esl_sq_Reuse(sq);
	}
      esl_sqfile_Position(sqfp, 0); /* rewind and make sure we get same reads again */
    }
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  remove(tmpfile);
}


/* utest_guess_mechanics()
 * SRE H3/70, 8 Apr 17
 *
 * Related to EPN bugfix a61ee23: esl_sqfile_GuessAlphabet() segfaults
 * when file contains 1 seq, file is >4096 bytes, seq is <4000
 * residues, because of a fault in the mechanics of sqio with
 * is_recording TRUE and is_linebased FALSE.
 *
 * This unit test exercises those mechanics, the original bug and
 * more. It is *not* testing GuessAlphabet() itself. The DNA sequences
 * in <sqarr> are so dirty, their alphabet cannot be reliably
 * detected. This unit test is hunting segfaults, not even looking at
 * the return status of _GuessAlphabet().
 */
static void
utest_guess_mechanics(ESL_ALPHABET *abc, ESL_SQ **sqarr, int N)
{
  char       *msg         = "sqio guess_mechanics unit test failure";
  char        tmpfile[32];
  FILE       *fp;          
  ESL_SQFILE *sqfp;
  int         i;
  int         alphatype;

  for (i = 0; i < N; i++)  // for each individual sequence in <sqarr>, one at a time:
    {
      strcpy(tmpfile, "esltmpXXXXXX");
      if (esl_tmpfile_named(tmpfile, &fp)                      != eslOK) esl_fatal(msg);
      if (esl_sqio_Write(fp, sqarr[i], eslSQFILE_FASTA, FALSE) != eslOK) esl_fatal(msg);
      fclose(fp);

      if (esl_sqfile_Open(tmpfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal(msg);

      esl_sqfile_GuessAlphabet(sqfp, &alphatype);
      // generally, the sequences are so dirty, GuessAlphabet() won't make a guess.
      // we're hunting segfaults in this utest.

      esl_sqfile_Close(sqfp);
      remove(tmpfile);
    }
}

/* utest_guess_empty_seq()
 * ML, 9 Oct 21
 *
 * Make sure that esl_sqfile_GuessAlphabet() returns eslNOALPHABET when
 * it tries to guess the alphabet on files containing only empty sequences.
 *
 * Related to bugfix 441a4d3: esl_sqfile_GuessAlphabet() did not handle
 * the case where esl_sqio_ReadWindow() would return eslEOD on empty
 * sequences, and returned an error instead of eslENOALPHABET.
 */
static void
utest_guess_empty_seq()
{
  char       *msg         = "sqio guess_empty_seq unit test failure";
  char        tmpfile[32];
  ESL_SQ*     seqs[2];
  FILE       *fp;
  ESL_SQFILE *sqfp;
  int         alphatype;

  if ((seqs[0] = esl_sq_CreateFrom("seqs0", "", NULL, NULL, NULL)) == NULL) esl_fatal(msg);
  if ((seqs[1] = esl_sq_CreateFrom("seqs1", "", NULL, NULL, NULL)) == NULL) esl_fatal(msg);

  strcpy(tmpfile, "esltmpXXXXXX");
  if (esl_tmpfile_named(tmpfile, &fp)                     != eslOK) esl_fatal(msg);
  if (esl_sqio_Write(fp, seqs[0], eslSQFILE_FASTA, FALSE) != eslOK) esl_fatal(msg);
  if (esl_sqio_Write(fp, seqs[1], eslSQFILE_FASTA, FALSE) != eslOK) esl_fatal(msg);
  fclose(fp);

  if (esl_sqfile_Open(tmpfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal(msg);
  if (esl_sqfile_GuessAlphabet(sqfp, &alphatype) != eslENOALPHABET)    esl_fatal(msg);
  esl_sqfile_Close(sqfp);
  remove(tmpfile);

  esl_sq_Destroy(seqs[0]);
  esl_sq_Destroy(seqs[1]);
}

#endif /*eslSQIO_TESTDRIVE*/
/*------------------ end, unit tests ----------------------------*/



/*****************************************************************
 * 16. Test driver.
 *****************************************************************/

#ifdef eslSQIO_TESTDRIVE
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-L",        eslARG_INT,   "8000",  NULL, NULL,  NULL,  NULL, NULL, "max length of test sequences",                     0 },
  { "-N",        eslARG_INT,    "100",  NULL, NULL,  NULL,  NULL, NULL, "number of test sequences",                         0 },
  { "-s",        eslARG_INT,      "0",  NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sqio module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc      = esl_alphabet_Create(eslDNA); /* DNA because some chars aren't legal in IUPAC DNA */
  ESL_RANDOMNESS *r        = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_SQ        **sqarr    = NULL;
  int             maxL     = esl_opt_GetInteger(go, "-L");
  int             N        = esl_opt_GetInteger(go, "-N");
  int             i;
  int             mode;
  char            tmpfile[32];
  char            ssifile[32];
  FILE           *fp       = NULL;
  char            c;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* Create an array of sequences we'll use for all the tests */
  synthesize_testseqs(r, abc, maxL, N, &sqarr);

  for (mode = 0; mode < 4; mode++) /* 0=ugly 1=spaced 2=normal 3=mixed*/
    {
      /* Write FASTA file to disk, and SSI index it */
      strcpy(tmpfile, "esltmpXXXXXX");
      if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("failed to make tmpfile");
      switch (mode) {
      case 0: for (i = 0; i < N; i++) write_ugly_fasta(r, fp, sqarr[i]); break;
      case 1: for (i = 0; i < N; i++) write_spaced_fasta(fp, sqarr[i]);  break;
      case 2:
	for (i = 0; i < N; i++) {
	  c = sqarr[i]->acc[0];	/* hack: hide the accession, so digital writer doesn't write it. */
	  sqarr[i]->acc[0] = '\0';
	  esl_sqio_Write(fp, sqarr[i], eslSQFILE_FASTA, TRUE); 
	  sqarr[i]->acc[0] = c;
	}
	break;
      case 3: for (i = 0; i < N; i++) write_mixed_fasta(r, fp, sqarr[i]); break;
      }
      fclose(fp);
      make_ssi_index(abc, tmpfile, eslSQFILE_FASTA, ssifile, mode);

      utest_read        (abc, sqarr, N, tmpfile, eslSQFILE_FASTA, mode);
      utest_read_info   (abc, sqarr, N, tmpfile, eslSQFILE_FASTA, mode);
      utest_read_window (abc, sqarr, N, tmpfile, eslSQFILE_FASTA, mode);
      utest_fetch_subseq(r, abc, sqarr, N, tmpfile, ssifile, eslSQFILE_FASTA);

      remove(tmpfile);
      remove(ssifile);
    }  

  utest_guess_mechanics(abc, sqarr, N);
  utest_subseq_geometry();
  utest_write          (abc, sqarr, N, eslMSAFILE_STOCKHOLM);
  utest_guess_empty_seq();

  for (i = 0; i < N; i++) esl_sq_Destroy(sqarr[i]);
  free(sqarr);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);

  fprintf(stderr, "#  status = ok\n");
  return 0;
}
#endif /*eslSQIO_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/



/*****************************************************************
 * 17. Examples
 *****************************************************************/

/* Example 3: using esl_sqio_Parse()
 */
#ifdef eslSQIO_EXAMPLE3
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

int
main(void)
{
  ESL_ALPHABET *abc       = NULL;
  ESL_SQ       *sq        = NULL;
  int           format    = eslSQFILE_FASTA;
  int           alphatype = eslAMINO;
  int           status;

  char *test = ">12345 TEST Test fasta buffer\n"
               "ARVAPVALPSACAPAGTQCLISGWGNTLSNGVNNPDLLQCVDAPVLSQADCEAAYPGEIT\n"
               "SSMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCALPDNPGVYTKVCNFVGWIQ\n"
               "DTIAAN";

  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);

  status = esl_sqio_Parse(test, strlen(test), sq, format);
  if      (status == eslEFORMAT) esl_fatal("Parse failed, invalid format");
  else if (status != eslOK)      esl_fatal("Unexpected error %d", status);

  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  return 0;
}
#endif /*eslSQIO_EXAMPLE3*/

/* Example 2: simplest, text mode
 */
#ifdef eslSQIO_EXAMPLE2
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

int
main(int argc, char **argv)
{
  ESL_SQ     *sq      = esl_sq_Create();
  ESL_SQFILE *sqfp;
  int         format  = eslSQFILE_FASTA;
  char       *seqfile = argv[1];
  int         status;

  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
  {     /* use each sequence for whatever you want */
    printf("%-40s length: %8ld   desclen: %lu\n", sq->name, (long) sq->L, (unsigned long) strlen(sq->desc));
    esl_sq_Reuse(sq);
  }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected read error %d", status);
  
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  return 0;
}
#endif /*eslSQIO_EXAMPLE2*/


/* Example 1: standard idiom for digital seqfile reading
 */
#ifdef eslSQIO_EXAMPLE
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
  /* name          type       default  env  range toggles reqs incomp  help         docgroup */
 { "-h",          eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",    0 },
 { "--dna",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use DNA alphabet",   0 },
 { "--rna",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use RNA alphabet",   0 },
 { "--amino",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use amino alphabet", 0 },
 { "--informat",  eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "set input format",   0 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "example of reading in standard sqio idiom, digital mode";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *seqfile   = esl_opt_GetArg(go, 1);
  ESL_ALPHABET *abc       = NULL;
  ESL_SQ       *sq        = NULL;
  ESL_SQFILE   *sqfp      = NULL;
  int           infmt     = eslSQFILE_UNKNOWN;
  int           alphatype = eslUNKNOWN;
  int           status;

  if (esl_opt_IsOn(go, "--informat")) {
    if ((infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslSQFILE_UNKNOWN)
      esl_fatal("%s is not a valid input sequence file format for --informat", esl_opt_GetString(go, "--informat")); 
  }

  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format couldn't be determined.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
    if      (status == eslENOALPHABET)  esl_fatal("Couldn't guess alphabet");
    else if (status == eslEFORMAT)      esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));     
    else if (status == eslENODATA)      esl_fatal("Sequence file empty?");
    else if (status != eslOK)           esl_fatal("Unexpected error guessing alphabet");
  }
  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);
  esl_sqfile_SetDigital(sqfp, abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {  
      esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, /*update=*/FALSE);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);
  
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslSQIO_EXAMPLE*/


