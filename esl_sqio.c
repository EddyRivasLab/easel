/* Unaligned sequence file i/o.
 * 
 * Contents:
 *    1. An <ESL_SQFILE> object, in text mode.
 *    2. An <ESL_SQFILE> object, in digital mode. [with <alphabet>]
 *    3. Using sequence file format codes
 *    4. Sequence reading (sequential).
 *    5. Sequence/subsequence fetching, random access [with <ssi>]
 *    6. Writing sequences.
 *    7. Internal routines shared by parsers.
 *    8. Internal routines for EMBL format (including Uniprot, TrEMBL)
 *    9. Internal routines for Genbank format
 *   10. Internal routines for FASTA format
 *   11. Internal routines for sq, msa interconversion [with <msa>]
 *   12. Benchmark driver.
 *   13. Unit tests.
 *   14. Test driver.
 *   15. Examples.
 *   16. Copyright and license.
 * 
 * This module shares remote evolutionary homology with Don Gilbert's
 * seminal, public domain ReadSeq package, though the last common
 * ancestor was circa 1991 and no recognizable vestiges are likely to
 * remain. Thanks Don!
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
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"		/* ssi aug adds ability to randomly access sequences/subsequences */
#endif
#include "esl_sqio.h"
#include "esl_sq.h"


/* Internal routines shared by parsers. */
static int  is_blankline(char *s);
static int  loadmem  (ESL_SQFILE *sqfp);
static int  loadbuf  (ESL_SQFILE *sqfp);
static int  nextchar (ESL_SQFILE *sqfp, char *ret_c);
static int  seebuf   (ESL_SQFILE *sqfp, int64_t maxn, int64_t *opt_nres, int64_t *opt_endpos);
static void addbuf   (ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nres);
static void skipbuf  (ESL_SQFILE *sqfp, int64_t nskip);
static int  read_nres(ESL_SQFILE *sqfp, ESL_SQ *sq, int64_t nskip, int64_t nres, int64_t *opt_actual_nres);

/* EMBL format; also Uniprot, TrEMBL */
static void config_embl(ESL_SQFILE *sqfp);
static void inmap_embl (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_embl(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_embl   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* Genbank format; also DDBJ */
static void config_genbank(ESL_SQFILE *sqfp);
static void inmap_genbank (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_genbank   (ESL_SQFILE *sqfp, ESL_SQ *sq);

/* FASTA format */
static void config_fasta(ESL_SQFILE *sqfp);
static void inmap_fasta (ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap);
static int  header_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_fasta   (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  write_fasta (FILE *fp, ESL_SQ *sq, int save_offsets);

/* Optional MSA<->sqio interoperability */
#ifdef eslAUGMENT_MSA
static int convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa);
#endif



/*****************************************************************
 *# 1. An <ESL_SQFILE> object, in text mode.
 *****************************************************************/ 

static int  sqfile_open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);

/* Function:  esl_sqfile_Open()
 * Synopsis:  Open a sequence file for reading.
 * Incept:    SRE, Thu Feb 17 08:22:16 2005 [St. Louis]
 *
 * Purpose:   Open a sequence file <filename> for reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The format of the file is asserted to be <format> (for
 *            example, <eslSQFILE_FASTA>). If <format> is
 *            <eslSQFILE_UNKNOWN> then the routine attempts to
 *            autodetect the file format.
 *            
 *            There are two special cases for <filename>. If
 *            <filename> is "-", the sequence data are read from a
 *            <STDIN> pipe. If <filename> ends in ".gz", the file is
 *            assumed to be compressed with <gzip>, and it is opened
 *            by a pipe from <gzip -dc>. Reading gzip files only works
 *            on POSIX-compliant systems that have pipes
 *            (specifically, the POSIX.2 popen() call). 
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
 *            Returns <eslENOTFOUND> if <filename> can't be found or
 *            opened.  Returns <eslEFORMAT> if the file is empty, or
 *            if autodetection is attempted and the format can't be
 *            determined.  On any error condition, <*ret_sqfp> is
 *            returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqfile_Open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
{
  int status;

  if ((status = sqfile_open(filename, format, env, ret_sqfp)) != eslOK) return status;

  /* text mode inmaps are less thorough than digital mode */
  switch ((*ret_sqfp)->format) {
  case eslSQFILE_EMBL:       inmap_embl(*ret_sqfp,    NULL); break;
  case eslSQFILE_UNIPROT:    inmap_embl(*ret_sqfp,    NULL); break;
  case eslSQFILE_GENBANK:    inmap_genbank(*ret_sqfp, NULL); break;
  case eslSQFILE_DDBJ:       inmap_genbank(*ret_sqfp, NULL); break;
  case eslSQFILE_FASTA:      inmap_fasta(*ret_sqfp,   NULL); break;
    /* stockholm: do nothing; MSAs don't use inmap */
  }
  return eslOK;
}


/* Function:  esl_sqfile_GuessFileFormat()
 * Synopsis:  Guess the format of an open <ESL_SQFILE>.
 * Incept:    SRE, Mon Jun 20 19:07:44 2005 [St. Louis]
 *
 * Purpose:   Try to guess the sequence file format of <sqfp>, and
 *            return the format code in <*ret_fmt>.
 *            
 *            First we attempt to guess based on the <filename>'s
 *            suffix. <*.fa> is assumed to be in FASTA format; <*.gb>
 *            is assumed to be in Genbank format; <*.sto> or <*.stk>
 *            are assumed to be in Stockholm multiple alignment file
 *            format.
 *            
 *            If that fails, we attempt to guess based on peeking at
 *            the first nonblank line of <filename>. If the line
 *            starts with $>$, we assume FASTA format; if the line
 *            starts with <ID>, we assume EMBL format; if the line
 *            starts with <LOCUS> or it contains the string <Genetic
 *            Sequence Data Bank> we assume Genbank format; if it
 *            starts with \verb+# STOCKHOLM+ we assume Stockholm format.
 *            
 *            If that fails too, return an <eslEFORMAT> error, and
 *            <*ret_fmt> is set to <eslSQFILE_UNKNOWN>.
 *            
 * Returns:   <eslOK> on success, and <*ret_fmt> contains
 *            a valid sequence file format code, such as 
 *            <eslSQFILE_FASTA>.
 *            
 *            Returns <eslEFORMAT> if we opened <filename> but it
 *            contains no nonblank lines, or if we peeked at the first
 *            nonblank line and still couldn't guess the format;
 *            <*ret_fmt> is then <eslSQFILE_UNKNOWN>.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 */
int
esl_sqfile_GuessFileFormat(ESL_SQFILE *sqfp, int *ret_fmt)
{
  int   n         = strlen(sqfp->filename);
  const char *sfx = NULL;
  int   is_gzip   = FALSE;
  int   nsfx;
  int   status;

  /* On any premature exit, *ret_fmt is eslSQFILE_UNKNOWN */
  *ret_fmt = eslSQFILE_UNKNOWN;

  /* Is <filename> gzip'ed? Look at suffix. */
  if (n > 3 && strcmp(sqfp->filename+n-3, ".gz") == 0) is_gzip = TRUE;

  /* Locate the suffix that might indicate format (ignore .gz) */
  for (nsfx = 1, sfx = sqfp->filename + n - 1 - (is_gzip ? 3 : 0);
       sfx != sqfp->filename && *sfx != '.'; 
       sfx--) 
    nsfx++;

  /* now sfx points either to filename (we didn't find a suffix) or to the . of the suffix,
   * and nsfx is the suffix length inclusive of the . 
   */
  
  /* Attempt to guess file format based on file name suffix. */
  if      (strncmp(sfx, ".fa",  3) == 0)  { *ret_fmt = eslSQFILE_FASTA;      return eslOK; }
  else if (strncmp(sfx, ".gb",  3) == 0)  { *ret_fmt = eslSQFILE_GENBANK;    return eslOK; }
  else if (strncmp(sfx, ".sto", 4) == 0)  { *ret_fmt = eslMSAFILE_STOCKHOLM; return eslOK; }
  else if (strncmp(sfx, ".stk", 4) == 0)  { *ret_fmt = eslMSAFILE_STOCKHOLM; return eslOK; }
    
  /* If that didn't work, we'll have a peek at the stream; 
   * turn recording on, and set for line based input.
   */

  if (sqfp->is_recording == -1) ESL_EXCEPTION(eslEINVAL, "sq file already too advanced");
  sqfp->is_recording = TRUE;
  sqfp->is_linebased = TRUE;
  loadbuf(sqfp);		/* now sqfp->buf is a line of the file */

  /* get first nonblank line */
  while (is_blankline(sqfp->buf)) {
    status = loadbuf(sqfp);
    if      (status == eslEOF) ESL_XFAIL(eslEFORMAT, sqfp->errbuf, "No data found in file");
    else if (status != eslOK)  goto ERROR;
  } 

  /* formats that can be determined from the first line: */
  if      (*(sqfp->buf) == '>')                                     *ret_fmt = eslSQFILE_FASTA;
  else if (strncmp(sqfp->buf, "ID   ", 5)    == 0)                  *ret_fmt = eslSQFILE_EMBL;
  else if (strncmp(sqfp->buf, "LOCUS   ", 8) == 0)                  *ret_fmt = eslSQFILE_GENBANK;
  else if (strstr(sqfp->buf, "Genetic Sequence Data Bank") != NULL) *ret_fmt = eslSQFILE_GENBANK;
#ifdef eslAUGMENT_MSA
  else if (strncmp(sqfp->buf, "# STOCKHOLM", 11) == 0)              *ret_fmt = eslMSAFILE_STOCKHOLM;
#endif

  /* reset the sqfp */
  sqfp->mpos         = 0;
  sqfp->is_recording = FALSE;
  sqfp->is_linebased = FALSE;
  free(sqfp->buf);
  sqfp->buf    = NULL;
  sqfp->balloc = 0;
  return (*ret_fmt == eslSQFILE_UNKNOWN) ? eslEFORMAT : eslOK;

 ERROR:
  sqfp->mpos         = 0;
  sqfp->is_recording = FALSE;
  sqfp->is_linebased = FALSE;
  if (sqfp->buf != NULL) free(sqfp->buf);
  return status;
}

/* Function:  esl_sqfile_Position()
 * Synopsis:  Reposition an open sequence file to an offset.
 * Incept:    SRE, Tue Mar 28 13:21:47 2006 [St. Louis]
 *
 * Purpose:   Reposition an open <sqfp> to offset <offset>.
 *            <offset> would usually be the first byte of a
 *            desired sequence record.
 *            
 *            Only normal sequence files can be positioned; not
 *            a standard input stream, gunzip stream, or a multiple
 *            alignment file interface.
 *            
 *            After <esl_sqfile_Position()> is called,
 *            <sqfp->linenumber> and other bookkeeping information is
 *            unknown. If caller knows it, it should set it
 *            explicitly.
 *            
 *            See the SSI module for manipulating offsets and indices.
 *
 * Returns:   <eslOK>     on success;
 *            <eslEOF>    if no data can be read from this position.
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails.
 *            <eslEMEM> on (re-)allocation failure.
 *            <eslEINVAL> if the <sqfp> is not positionable.
 */
int
esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset)
{
  if (sqfp->do_stdin)    ESL_EXCEPTION(eslEINVAL, "can't Position() in standard input");
  if (sqfp->do_gzip)     ESL_EXCEPTION(eslEINVAL, "can't Position() in a gzipped file");
  if (sqfp->afp != NULL) ESL_EXCEPTION(eslEINVAL, "can't use esl_sqfile_Position() in an alignment file");

  if (fseeko(sqfp->fp, offset, SEEK_SET) != 0) ESL_EXCEPTION(eslESYS, "fseeko() failed");

  sqfp->currpl     = -1;
  sqfp->curbpl     = -1;
  sqfp->prvrpl     = -1;
  sqfp->prvbpl     = -1;
  sqfp->linenumber = -1;
  sqfp->L          = -1;
  sqfp->mpos       = sqfp->mn;	/* this forces loadbuf to load new data */
  return loadbuf(sqfp);
}

/* Function:  esl_sqfile_Close()
 * Synopsis:  Close a sequence file.
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
  else 
#endif
  if (! sqfp->do_stdin && sqfp->fp != NULL) fclose(sqfp->fp);
  if (sqfp->filename != NULL) free(sqfp->filename);
  if (sqfp->ssifile  != NULL) free(sqfp->ssifile);
  if (sqfp->mem      != NULL) free(sqfp->mem);
  if (sqfp->balloc   > 0)     free(sqfp->buf);
#ifdef eslAUGMENT_SSI
  if (sqfp->ssi      != NULL) esl_ssi_Close(sqfp->ssi);
#endif

#ifdef eslAUGMENT_MSA
  if (sqfp->afp      != NULL) esl_msafile_Close(sqfp->afp);
  if (sqfp->msa      != NULL) esl_msa_Destroy(sqfp->msa);
#endif /*eslAUGMENT_MSA*/
  free(sqfp);
  return;
}


/* sqfile_open():
 * This is the routine that actually opens an ESL_SQFILE.
 * esl_sqfile_Open() and esl_sqfile_OpenDigital() are
 * small wrappers around it.
 */
static int
sqfile_open(const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp    = NULL;
  char       *envfile = NULL;
  int         status;		/* return status from an ESL call */
  int         n;

  ESL_ALLOC(sqfp, sizeof(ESL_SQFILE));
  *ret_sqfp          = NULL;

  sqfp->fp           = NULL;
  sqfp->filename     = NULL;
  sqfp->do_gzip      = FALSE;
  sqfp->do_stdin     = FALSE;
  sqfp->errbuf[0]    = '\0';

  sqfp->mem          = NULL;
  sqfp->allocm       = 0;
  sqfp->mn           = 0;
  sqfp->mpos         = 0;
  sqfp->moff         = -1;
  sqfp->is_recording = FALSE;

  sqfp->buf          = NULL;
  sqfp->boff         = 0;
  sqfp->balloc       = 0;
  sqfp->nc           = 0;
  sqfp->bpos         = 0;
  sqfp->L            = 0;
  sqfp->linenumber   = 1;

  sqfp->bookmark_offset  = 0;
  sqfp->bookmark_linenum = 0;

  sqfp->do_digital   = FALSE;
  sqfp->abc          = NULL;

  sqfp->format       = format;
  sqfp->is_linebased = FALSE;
  sqfp->eof_is_ok    = FALSE;
  sqfp->parse_header = NULL;
  sqfp->parse_end    = NULL;

  sqfp->afp        = NULL;
  sqfp->msa        = NULL;
  sqfp->idx        = -1;

  sqfp->ssifile    = NULL;
  sqfp->rpl        = -1;	/* -1 = not set yet */
  sqfp->bpl        = -1;	/* (ditto) */
  sqfp->prvrpl     = -1;    	/* (ditto) */
  sqfp->prvbpl     = -1;        /* (ditto) */
  sqfp->currpl     = -1;	
  sqfp->curbpl     = -1;	
  sqfp->ssi        = NULL;

  /* Open the file, either in cwd or in a directory listed in <env>.  */
  if (strcmp(filename, "-") == 0) /* stdin special case */
    {
      if ((status = esl_strdup("[STDIN]", -1, &(sqfp->filename))) != eslOK) goto ERROR;
      sqfp->fp       = stdin;
      sqfp->do_stdin = TRUE;
    }
  else
    { /* Check the current working directory first. */
      if ((sqfp->fp = fopen(filename, "r")) != NULL) {
	if ((status = esl_strdup(filename, -1, &(sqfp->filename))) != eslOK) goto ERROR;
      }
      /* if it's not there, then check in directory list provided by <env>. */
      else if (env != NULL && esl_FileEnvOpen(filename, env, &(sqfp->fp), &envfile) == eslOK) {
	if ((status = esl_strdup(envfile, -1, &(sqfp->filename))) != eslOK) goto ERROR;
      }
      else { status = eslENOTFOUND; goto ERROR;}
    }

  /* Deal with the .gz special case: to popen(), "success" only means
   * it found and executed gzip -dc.  If gzip -dc doesn't find our
   * file, popen() still blithely returns success, so we have to be
   * sure the file exists. That's why we fopen()'ed it above, only to
   * close it and popen() it here.
   */                           
#ifdef HAVE_POPEN
  n = strlen(sqfp->filename);
  if (n > 3 && strcmp(sqfp->filename+n-3, ".gz") == 0) 
    {
      char *cmd;
      fclose(sqfp->fp);
      ESL_ALLOC(cmd, sizeof(char) * (n+1+strlen("gzip -dc ")));
      sprintf(cmd, "gzip -dc %s", sqfp->filename);
      sqfp->fp = popen(cmd, "r");
      if (sqfp->fp == NULL) { status = eslENOTFOUND; goto ERROR; }
      sqfp->do_gzip  = TRUE;
      free(cmd);
    }
#endif /*HAVE_POPEN*/

  /* If we don't know the format yet, autodetect it now. */
  if (sqfp->format == eslSQFILE_UNKNOWN &&
      (status = esl_sqfile_GuessFileFormat(sqfp, &(sqfp->format))) != eslOK)
    goto ERROR;

  /* Configure the <sqfp>'s parser for this format. */
  switch (sqfp->format) {
  case eslSQFILE_EMBL:     config_embl(sqfp);    break;
  case eslSQFILE_UNIPROT:  config_embl(sqfp);    break;
  case eslSQFILE_GENBANK:  config_genbank(sqfp); break;
  case eslSQFILE_DDBJ:     config_genbank(sqfp); break;
  case eslSQFILE_FASTA:    config_fasta(sqfp);   break;

#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM: 
    sqfp->is_linebased = TRUE;
    sqfp->eof_is_ok    = FALSE;	/* no-op for msa's */
    sqfp->parse_header = NULL;	/* no-op for msa's */
    sqfp->parse_end    = NULL;	/* no-op for msa's */
    if ((status = esl_msafile_Open(filename, sqfp->format, env, &(sqfp->afp))) != eslOK) goto ERROR;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  /* Preload the first line or chunk of file. */
  if (! esl_sqio_IsAlignment(sqfp->format))
    {
      status = loadbuf(sqfp);
      if      (status == eslEOF) { status = eslEFORMAT; goto ERROR; }
      else if (status != eslOK)  { goto ERROR; }
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
/*------------------- ESL_SQFILE open/close -----------------------*/


/*****************************************************************
 *# 2. An <ESL_SQFILE> object, in digital mode [with <alphabet>]
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  esl_sqfile_OpenDigital()
 * Synopsis:  Open an <ESL_SQFILE> for digital input.
 * Incept:    SRE, Fri May  9 09:17:48 2008 [Janelia]
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

  if ((status = sqfile_open(filename, format, env, ret_sqfp)) != eslOK) return status;
  return esl_sqfile_SetDigital(*ret_sqfp, abc);
}

/* Function:  esl_sqfile_SetDigital()
 * Synopsis:  Set an open <ESL_SQFILE> to read in digital mode.
 * Incept:    SRE, Fri May  9 09:21:31 2008 [Janelia]
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
  switch (sqfp->format) {
  case eslSQFILE_EMBL:       inmap_embl(sqfp,    abc->inmap); break;
  case eslSQFILE_UNIPROT:    inmap_embl(sqfp,    abc->inmap); break;
  case eslSQFILE_GENBANK:    inmap_genbank(sqfp, abc->inmap); break;
  case eslSQFILE_DDBJ:       inmap_genbank(sqfp, abc->inmap); break;
  case eslSQFILE_FASTA:      inmap_fasta(sqfp,   abc->inmap); break;
    /* stockholm: do nothing (no inmap used for MSAs */
  }

  if (esl_sqio_IsAlignment(sqfp->format))
    esl_msafile_SetDigital(sqfp->afp, abc);

  sqfp->do_digital = TRUE;
  sqfp->abc        = abc;
  return eslOK;
}

/* Function:  esl_sqfile_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open <ESL_SQFILE>.
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

  /* Special case: for MSA files, hand this off to msafile_GuessAlphabet. */
#ifdef eslAUGMENT_MSA
  if (esl_sqio_IsAlignment(sqfp->format)) return esl_msafile_GuessAlphabet(sqfp->afp, ret_type);
#endif

  /* set the sqfp to record; we'll rewind afterwards and use the recording */
  sqfp->is_recording = TRUE;

  if ((sq = esl_sq_Create()) == NULL) { status = eslEMEM; goto ERROR; }

  status = esl_sqio_ReadWindow(sqfp, 0, 4000, sq);
  if      (status == eslEOF) { status = eslENODATA; goto ERROR; }
  else if (status != eslOK)  goto ERROR; 

  if ((status = esl_sq_GuessAlphabet(sq, ret_type)) != eslOK) goto ERROR;

  /* reset the sqfp, so it uses the recording next */
  sqfp->mpos         = 0;
  sqfp->is_recording = FALSE;
  if ((status = loadbuf(sqfp)) != eslOK) ESL_EXCEPTION(status, "buffer load failed, but shouldn't have");
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  *ret_type      = eslUNKNOWN;
  return status;
}
#endif /*eslAUGMENT_ALPHABET*/
/*-------------- end, digital mode ESL_SQFILE -------------------*/




/*****************************************************************
 *# 3. Using sequence file format codes
 *****************************************************************/ 

/* Function:  esl_sqio_FormatCode()
 * Synopsis:  Convert a string to an internal format code.
 * Incept:    SRE, Sun Feb 27 09:18:36 2005 [St. Louis]
 *
 * Purpose:   Given <fmtstring>, return format code.  For example, if
 *            <fmtstring> is "fasta", returns <eslSQFILE_FASTA>. Returns 
 *            <eslSQFILE_UNKNOWN> if <fmtstring> doesn't exactly match a 
 *            known format.
 *            
 *            Matching is case insensitive; fasta, FASTA, and FastA
 *            all return <eslSQFILE_FASTA>, for example.
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
 * Synopsis:  Return TRUE for alignment file format codes.
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
 */
int
esl_sqio_IsAlignment(int fmt)
{
  return (fmt >= 100 ? TRUE : FALSE);
}



/*****************************************************************
 *# 4. Sequence reading (sequential)
 *****************************************************************/ 

/* Function:  esl_sqio_Read()
 * Synopsis:  Read the next sequence from a file.
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
esl_sqio_Read(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;
  int64_t epos;
  int64_t n;

#ifdef eslAUGMENT_MSA
  if (esl_sqio_IsAlignment(sqfp->format))
    {
      ESL_SQ *tmpsq = NULL;
      if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
	{ /* we need to load a new alignment? */
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
#endif

  /* Main case: read next seq from sqfp's stream */
  if (sqfp->nc == 0) return eslEOF;
  if ((status = sqfp->parse_header(sqfp, sq)) != eslOK) return status; /* EOF, EFORMAT */

  do {
    if ((status = seebuf(sqfp, -1, &n, &epos)) == eslEFORMAT) return status;
    if (esl_sq_GrowTo(sq, sq->n + n) != eslOK) return eslEMEM;
    addbuf(sqfp, sq, n);
    sqfp->L   += n;
    sq->eoff   = sqfp->boff + epos - 1;
    if (status == eslEOD)     break;
  } while ((status = loadbuf(sqfp)) == eslOK);
    
  if      (status == eslEOF)
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
 * Incept:    SRE, Fri May 16 09:24:21 2008 [Janelia]
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
 */
int
esl_sqio_ReadInfo(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;
  int64_t epos;
  int64_t n;

#ifdef eslAUGMENT_MSA
  if (esl_sqio_IsAlignment(sqfp->format))
    {
      ESL_SQ *tmpsq = NULL;
      if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
	{ /* we need to load a new alignment? */
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
#endif

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
 * Incept:    SRE, Fri May 16 13:42:51 2008 [Janelia]
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
  int     actual_start;
  int64_t nres;
  int64_t line;
  off_t   offset;
  int     status;
  ESL_SQ *tmpsq = NULL;

#ifdef eslAUGMENT_MSA
  if (esl_sqio_IsAlignment(sqfp->format))
    {
      /* special: if we're initializing a revcomp window read, back sqfp->idx up one */
      if (W < 0 && sq->start == 0) sqfp->idx--;

      if (sqfp->msa == NULL || sqfp->idx >= sqfp->msa->nseq)
	{ /* need new alignment? */
	  esl_msa_Destroy(sqfp->msa);
	  status = esl_msa_Read(sqfp->afp, &(sqfp->msa));
	  if (status == eslEFORMAT)
	    { /* oops, a parse error; upload the error info from afp to sqfp */
	      sqfp->linenumber = sqfp->afp->linenumber;
	      strcpy(sqfp->errbuf, sqfp->afp->errbuf); /* errbufs same size! */ 
	      return eslEFORMAT;
	    }
	  else if (status != eslOK) goto ERROR;
	  sqfp->idx = 0;
	}
      
      /* grab appropriate seq from alignment into tmpsq */
      if ((status = esl_sq_FetchFromMSA(sqfp->msa, sqfp->idx, &tmpsq)) != eslOK) goto ERROR;

      /* Figure out tmpsq coords we'll put in sq */
      if (W > 0)
	{			/* forward strand */
	  sq->C     = ESL_MIN(sq->n, C);
	  sq->start = sq->end - sq->C + 1;
	  sq->end   = ESL_MIN(tmpsq->L, sq->end + W);
	  sq->n     = sq->end - sq->start + 1;
	  sq->W     = sq->n - sq->C;
	}
      else 
	{			/* reverse strand */
	  if (sq->L == -1) ESL_XEXCEPTION(eslESYNTAX, "Can't read reverse complement until you've read forward strand");

	  sq->C     = ESL_MIN(sq->n, sq->end + C - 1);
	  sq->end   = (sq->start == 0 ? sq->L : sq->end + sq->C - 1);
	  sq->start = ESL_MAX(1, sq->end + W - sq->C - 1);
	  sq->n     = sq->end - sq->start + 1;
	  sq->W     = sq->n - sq->C;
	}

      if (sq->W == 0)		/* no new sequence? that's the EOD case */
	{
	  sq->start      = 0;
	  sq->end        = 0;
	  sq->C          = 0;
	  sq->W          = 0;
	  sq->n          = 0;
	  sq->L          = tmpsq->L;
	  if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL;
	  else                 sq->seq[0] = '\0';

	  sqfp->idx++;
	  esl_sq_Destroy(tmpsq);
	  return eslEOD;
	}

      /* Copy the sequence frag.  */
      if (tmpsq->ss != NULL && sq->ss == NULL) ESL_ALLOC(sq->ss, sizeof(char) * (sq->salloc)); /* this *must* be for salloc  */
      esl_sq_GrowTo(sq, sq->n);
      if (tmpsq->seq != NULL) 
	{	/* text mode */
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
      if (W < 0 && (status = esl_sq_ReverseComplement(sq)) != eslOK) 
	ESL_XFAIL(eslEINVAL, sqfp->errbuf, "Can't reverse complement that sequence window");
	  
      /* Copy annotation */
      if ((status = esl_sq_SetName     (sq, tmpsq->name))   != eslOK) goto ERROR;
      if ((status = esl_sq_SetSource   (sq, tmpsq->name))   != eslOK) goto ERROR;
      if ((status = esl_sq_SetAccession(sq, tmpsq->acc))    != eslOK) goto ERROR;
      if ((status = esl_sq_SetDesc     (sq, tmpsq->desc))   != eslOK) goto ERROR;
      sq->roff = -1;
      sq->doff = -1;
      sq->eoff = -1;

      esl_sq_Destroy(tmpsq);
      return eslOK;
    }
#endif /* we've completely handled the alignment file case above. */

  /* Now for the normal case: we're reading a normal unaligned seq file, not an alignment. */

  /* Negative W indicates reverse complement direction */
  if (W < 0)	
    {
      if (sq->L == -1) ESL_EXCEPTION(eslESYNTAX, "Can't read reverse complement until you've read forward strand");

      if (sq->end == 1) 
	{ /* last end == 1 means last window was the final one on reverse strand,
	   * so we're EOD; jump back to last forward position. 
	   */
	  if (sqfp->bookmark_offset > 0) {
	    if (esl_sqfile_Position(sqfp, sqfp->bookmark_offset) != eslOK)
	      ESL_EXCEPTION(eslECORRUPT, "Failed to reposition seq file at last forward bookmark");
	    sqfp->linenumber = sqfp->bookmark_linenum;
	  } else {
	    sqfp->nc = 0;	/* signals EOF */
	  }
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

      /* If s == 0, we haven't read any reverse windows yet; 
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
	  sq->C     = ESL_MIN(C, sq->L - sq->end + 1);  /* based on prev window's end */
	  sq->end   = sq->end + sq->C - 1;                /* also based on prev end     */
	  sq->start = ESL_MAX(1, (sq->end - W - sq->C + 1));
	  sq->W     = sq->end - sq->start + 1 - sq->C;
	}

      /* Now position for a subseq fetch of <start..end> on fwd strand, using SSI offset calc  */
      if (sq->doff == 0) ESL_EXCEPTION(eslECORRUPT, "can't happen: sq didn't store data offset");

      if (sqfp->bpl == 0 || sqfp->rpl == 0) /* no help; brute force resolution. */
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
      else		/* line resolution */
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
	  sq->C         = 0;	/* no context in first window                   */
	  sq->L         = -1;	/* won't be known 'til EOD.                     */
	  sqfp->L       = 0;	/* init to 0, so we can count residues as we go */
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
	  if (sq->seq != NULL) memmove(sq->seq,   sq->seq + sq->n - sq->C,     sq->C);
	  else                 memmove(sq->dsq+1, sq->dsq + sq->n - sq->C + 1, sq->C);
	  sq->start = sqfp->L - sq->C + 1;
	  sq->n = C;
	}      

      if ((status = esl_sq_GrowTo(sq, C+W)) != eslOK)                return status; /* EMEM    */
      status = read_nres(sqfp, sq, 0, W, &nres);
      sqfp->L += nres;

      if (status == eslEOD)	
	{ /* Forward strand is done. 0 residues were read. Return eslEOD and an empty (info) <sq>. */
	  if ((status = sqfp->parse_end(sqfp, sq)) != eslOK) return status;

	  sq->start      = 0;
	  sq->end        = 0;
	  sq->C          = 0;
	  sq->W          = 0;
	  sq->L          = sqfp->L;
	  sq->n          = 0;

	  if (sqfp->nc > 0) {
	    sqfp->bookmark_offset  = sqfp->boff+sqfp->bpos; /* remember where the next seq starts. */
	    sqfp->bookmark_linenum = sqfp->bookmark_linenum;
	  } else { 
            sqfp->bookmark_offset  = 0;	                    /* signals for EOF, no more seqs        */
	    sqfp->bookmark_linenum = 0;
	  }

	  if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL; /* erase the saved context */
	  else                 sq->seq[0] = '\0';
	  return eslEOD;
	}
      else if (status == eslOK)
	{ /* Forward strand is still in progress. <= W residues were read. Return eslOK. */
	  sq->end        = sq->start + sq->C + nres - 1;	  
	  sq->W          = nres;	  
	  return eslOK;
	}
      else return status;	/* EFORMAT,EMEM */
    }
  /*NOTREACHED*/
  return eslOK;

 ERROR:
  if (tmpsq != NULL) esl_sq_Destroy(tmpsq);
  return status;
}

/* Function:  esl_sqio_Echo()
 * Synopsis:  Echo a sequence's record onto output stream.
 * Incept:    SRE, Wed Apr  2 16:32:21 2008 [Janelia]
 *
 * Purpose:   Given a complete <sq> that we have read by some means
 *            from an open <sqfp>; echo that sequence's record
 *            onto the output stream <ofp>. 
 *
 *            This allows records to be regurgitated exactly as they
 *            appear, rather than writing the subset of information
 *            stored in an <ESL_SQ>. <esl-sfetch> in the miniapps uses
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
 *            
 *            
 */
int
esl_sqio_Echo(ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)
{
  int     status;
  int64_t save_linenumber;
  int     save_currpl;
  int     save_curbpl;
  int     save_prvrpl;
  int     save_prvbpl;
  int64_t save_L;
  off_t   save_offset;
  int     n;
  int     nwritten;

  if (sqfp->do_stdin)                     ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from standard input");
  if (sqfp->do_gzip)                      ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from a gzipped file");
  if (esl_sqio_IsAlignment(sqfp->format)) ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from an alignment file");
  if (sq->roff == -1 || sq->eoff == -1)   ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence without disk offset info");

  save_linenumber = sqfp->linenumber;
  save_offset     = sqfp->boff;
  save_currpl     = sqfp->currpl;
  save_curbpl     = sqfp->curbpl;
  save_prvrpl     = sqfp->prvrpl;
  save_prvbpl     = sqfp->prvbpl;
  save_L          = sqfp->L;

  status = esl_sqfile_Position(sqfp, sq->roff);
  if      (status == eslEOF) ESL_EXCEPTION(eslECORRUPT, "repositioning failed; bad offset?");
  else if (status != eslOK)  return status;

  while (sqfp->boff + sqfp->nc < sq->eoff)
    {
      if (fwrite(sqfp->buf, sizeof(char), sqfp->nc, ofp) != sqfp->nc) ESL_EXCEPTION(eslESYS, "fwrite() failed");
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
/*------------------ end, sequential sequence input -------------*/


/*****************************************************************
 *# 5. Sequence/subsequence fetching, random access [with <ssi>]
 *****************************************************************/
#ifdef eslAUGMENT_SSI

/* Function:  esl_sqfile_OpenSSI()
 * Synopsis:  Opens an SSI index associated with a sequence file.
 * Incept:    SRE, Wed Apr  2 10:21:04 2008 [Janelia]
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
 *            the sequence filename, while replacing any existing
 *            terminal dot-suffix with <.ssi>. For example, the SSI
 *            index for <foo> is <foo.ssi>, for <./foo.fa> is
 *            <./foo.ssi>, and for </my/path/to/foo.1.fa> is
 *            </my/path/to/foo.1.ssi>.
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
 * Incept:    SRE, Wed Apr  2 09:51:11 2008 [Janelia]
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

  if (sqfp->ssi == NULL)                          ESL_EXCEPTION(eslEINVAL,"Need an open SSI index to call esl_sqfile_PositionByKey()");
  if ((status = esl_ssi_FindName(sqfp->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  return esl_sqfile_Position(sqfp, offset);
}


/* Function:  esl_sqfile_PositionByNumber()
 * Synopsis:  Use SSI to reposition by sequence number
 * Incept:    SRE, Wed Apr  2 17:24:38 2008 [Janelia]
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

  if (sqfp->ssi == NULL)                          ESL_EXCEPTION(eslEINVAL,"Need open SSI index to call esl_sqfile_PositionByNumber()");
  if ((status = esl_ssi_FindNumber(sqfp->ssi, which, &fh, &offset, NULL, NULL, NULL)) != eslOK) return status;
  return esl_sqfile_Position(sqfp, offset);
}


/* Function:  esl_sqio_Fetch()
 * Synopsis:  Fetch a complete sequence, using SSI indexing.
 * Incept:    SRE, Fri May 16 13:25:00 2008 [Janelia]
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

  if (sqfp->ssi == NULL) ESL_FAIL(eslEINVAL, sqfp->errbuf, "No SSI index for %s; can't fetch subsequences", sqfp->filename);
  if ((status = esl_sqfile_PositionByKey(sqfp, key)) != eslOK) return status;
  if ((status = esl_sqio_Read(sqfp, sq))             != eslOK) return status;
  return eslOK;
}
  
/* Function:  esl_sqio_FetchInfo()
 * Synopsis:  Fetch a sequence's info, using SSI indexing.
 * Incept:    SRE, Fri May 16 13:25:00 2008 [Janelia]
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
 * Incept:    SRE, Tue May 13 11:00:04 2008 [Janelia]
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
 * Returns:   <eslOK> on soccess.
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
  uint16_t fh;			/* SSI file handle */
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
  nres  = end - start + 1;	  /* how many residues do we need to read as subseq                  */

  if ((status = esl_sq_GrowTo(sq, nres)) != eslOK) return status;
  status = read_nres(sqfp, sq, nskip, nres, &n);
  if (status != eslOK || n < nres) ESL_EXCEPTION(eslEINCONCEIVABLE, "Failed to fetch subsequence residues -- corrupt coords?");

  /* Set the coords */
  sq->start = start;
  sq->end   = end;
  sq->C     = 0;
  sq->W     = sq->n;
  sq->L     = (L > 0 ? L : -1);
  esl_sq_SetName  (sq, "%s/%d-%d", source, start, end);
  esl_sq_SetSource(sq, source);
  return eslOK;
}  
#endif /*eslAUGMENT_SSI*/
/*------------- end, random sequence access with SSI -------------------*/


/*****************************************************************
 *# 6. Writing sequences.
 *****************************************************************/

/* Function:  esl_sqio_Write()
 * Synopsis:  Write a sequence to a file.
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
  if (esl_sqio_IsAlignment(format))
    {
      if ((status = convert_sq_to_msa(s, &msa)) != eslOK) return status;
      status = esl_msa_Write(fp, msa, format);
      esl_msa_Destroy(msa);
      return status;
    }
#endif
  switch (format) {
  case eslSQFILE_FASTA:   status = write_fasta  (fp, s, FALSE); break;
  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "can't write that format");
  }
  return status;
}
/*-------------------- writing sequences  -----------------------*/



/*****************************************************************
 * 7. Internal routines shared by parsers
 *****************************************************************/

static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++) 
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}

/* loadmem() 
 *
 * Load the next block of data from stream into mem buffer,
 * either concatenating to previous buffer (if we're recording) or
 * overwriting (if not). 
 * 
 * This block is loaded at sqfp->mem + sqfp->mpos.
 * 
 * Upon return:
 * sqfp->mem     now contains up to eslREADBUFSIZE more chars
 * sqfp->mpos    is position of first byte in newly read block
 * sqfp->allocm  may have increased by eslREADBUFSIZE, if we concatenated
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
  int   n;
  int   status;

  if (sqfp->is_recording == TRUE)
    {
      if (sqfp->mem == NULL) sqfp->moff = ftello(sqfp->fp);        /* first time init of the offset */
      ESL_RALLOC(sqfp->mem, tmp, sizeof(char) * (sqfp->allocm + eslREADBUFSIZE));
      sqfp->allocm += eslREADBUFSIZE;
      n = fread(sqfp->mem + sqfp->mpos, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      sqfp->mn += n;
    }
  else
    {
      if (sqfp->mem == NULL) {
	ESL_ALLOC(sqfp->mem, sizeof(char) * eslREADBUFSIZE);
	sqfp->allocm = eslREADBUFSIZE;
      }
      sqfp->is_recording = -1;	/* no more recording is possible now */
      sqfp->mpos = 0;
      sqfp->moff = ftello(sqfp->fp);
      n          = fread(sqfp->mem, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      sqfp->mn   = n;
    }
  return (n == 0 ? eslEOF : eslOK);

 ERROR:
  return status;
}

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
      sqfp->mpos  += sqfp->mn;
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
	    ESL_RALLOC(sqfp->buf, tmp, sizeof(char) * (sqfp->balloc + eslREADBUFSIZE));
	    sqfp->balloc += eslREADBUFSIZE;
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
	  ESL_RALLOC(sqfp->buf, tmp, sizeof(char) * (sqfp->balloc + eslREADBUFSIZE));
	  sqfp->balloc += eslREADBUFSIZE;
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
  int64_t nres2 = 0;	/* an optimization for determining lastrpl from nres, without incrementing lastrpl on every char */
  int     sym;
  ESL_DSQ x;
  int     lasteol = sqfp->bpos - 1;
  int     status  = eslOK;

  if (maxn == -1) maxn = sqfp->nc; /* makes for a more efficient test. nc is a guaranteed upper bound on nres */

  for (bpos = sqfp->bpos; nres < maxn && bpos < sqfp->nc; bpos++)
    {
      sym = sqfp->buf[bpos];
      if (!isascii(sym)) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Non-ASCII character %c in sequence", sym); 
      x   = sqfp->inmap[sym];

      if      (x <= 127) nres++;
      else if (x == eslDSQ_EOL) 
	{
	  if (sqfp->curbpl != -1) sqfp->curbpl += bpos - lasteol;
	  if (sqfp->currpl != -1) sqfp->currpl += nres - nres2;
	  nres2        += nres - nres2;

	  if (sqfp->rpl != 0 && sqfp->prvrpl != -1) { /* need to ignore counts on last line in record, hence cur/prv */
	    if      (sqfp->rpl    == -1)        sqfp->rpl = sqfp->prvrpl; /* init */
	    else if (sqfp->prvrpl != sqfp->rpl) sqfp->rpl = 0;	          /* inval*/
	  }
	  if (sqfp->bpl != 0 && sqfp->prvbpl != -1) {
	    if      (sqfp->bpl    == -1)        sqfp->bpl = sqfp->prvbpl; /* init  */
	    else if (sqfp->prvbpl != sqfp->bpl) sqfp->bpl = 0;            /* inval */
	  }

	  sqfp->prvbpl  = sqfp->curbpl;
	  sqfp->prvrpl  = sqfp->currpl;
	  sqfp->curbpl  = 0;
	  sqfp->currpl  = 0;
	  lasteol       = bpos;
	  if (sqfp->linenumber != -1) sqfp->linenumber++; 
	}
      else if (x == eslDSQ_ILLEGAL) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Illegal character %c", sym); 
      else if (x == eslDSQ_EOD)     { status = eslEOD; break; }
      else if (x != eslDSQ_IGNORED) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "inmap corruption?");
    }

  if (sqfp->curbpl != -1) sqfp->curbpl += bpos - lasteol - 1;
  if (sqfp->currpl != -1) sqfp->currpl += nres - nres2;
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
    if (x <= 127) nskip--;	/* skip IGNORED, EOL. */
  }
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

  if        (status == eslEOF) { 
    if (! sqfp->eof_is_ok) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Premature EOF before end of seq record");
    n = 0;
  } 

  n = ESL_MIN(nres, n); 
  addbuf(sqfp, sq, n);	  /* bpos now at last residue + 1 if OK/EOD, 0 if EOF  */    
  actual_nres += n;

  if (sq->dsq != NULL) sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else                 sq->seq[sq->n]   = '\0';
  
  if (status == eslEOD) { 
    sqfp->bpos = epos; 
  }

  if (opt_actual_nres != NULL) *opt_actual_nres = actual_nres;
  return (actual_nres == 0 ? eslEOD : eslOK);
}


/*--------------- end, buffer-based parsers --------------------*/


/*****************************************************************
 * 8. Internal routines for EMBL format (including Uniprot, TrEMBL)
 *****************************************************************/ 
/* EMBL and Uniprot protein sequence database format.
 * See: http://us.expasy.org/sprot/userman.html
 */
static void
config_embl(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased      = TRUE;
  sqfp->eof_is_ok         = FALSE;	/* records end with // */
  sqfp->parse_header      = &header_embl;
  sqfp->parse_end         = &end_embl;
}

static void
inmap_embl(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
  sqfp->inmap['/']  = eslDSQ_EOD;
}

/* header_embl()
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
  int   toklen;
  int   status;

  /* Find first line:
   * "Each entry must begin with an identification line (ID)..."
   * "The two-character line-type code that begins each line is always
   *  followed by three blanks..."
   */
  if (sqfp->nc == 0) return eslEOF;
  while (is_blankline(sqfp->buf)) {
    if ((status = loadbuf(sqfp)) == eslEOF) return eslEOF; /* normal */
    else if (status != eslOK) return status; /* abnormal */
  } 

  /* ID line is defined as:
   *     ID   ENTRY_NAME DATA_CLASS; MOLECULE_TYPE; SEQUENCE_LENGTH.
   * We're only after the ENTRY_NAME.
   */
  if (strncmp(sqfp->buf, "ID   ", 5) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find ID line");
  
  s = sqfp->buf+5;
  if ((status = esl_strtok(&s, " ", &tok, &toklen)) != eslOK)
    ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to parse name on ID line");
  if ((status = esl_sq_SetName(sq, tok)) != eslOK) return status;
  sq->roff = sqfp->boff;	/* record the offset of the ID line */
  
  /* Look for SQ line; parsing optional info as we go.
   */
  do {
    if ((status = loadbuf(sqfp)) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find SQ line");
						     
    /* "The format of the AC line is:
     *    AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]
     *  Researchers who wish to cite entries in their publications
     *  should always cite the first accession number. This is
     *  commonly referred to as the 'primary accession
     *  number'."
     */
    if (strncmp(sqfp->buf, "AC   ", 5) == 0)
      {
	s = sqfp->buf+5;
	if ((status = esl_strtok(&s, ";", &tok, &toklen)) != eslOK)
	  ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to parse accession on AC line");
	if ((status = esl_sq_SetAccession(sq, tok)) != eslOK) return status;
      }

    /* "The format of the DE line is:
     *    DE   Description.
     * ...In cases where more than one DE line is required, the text is
     * only divided between words and only the last DE line is
     * terminated by a period."
     */
    if (strncmp(sqfp->buf, "DE   ", 5) == 0)
      {
	s = sqfp->buf+5; 
	esl_strchop(s, sqfp->nc);
	if ((status = esl_sq_AppendDesc(sq, s)) != eslOK) 
	  ESL_FAIL(status, sqfp->errbuf, "Failed to parse description on DE line");
      }

    /* "The format of the SQ line is:
     *  SQ   SEQUENCE XXXX AA; XXXXX MW; XXXXXXXXXXXXXXXX CRC64;"
     */
  } while (strncmp(sqfp->buf, "SQ   ", 5) != 0);

  if (loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find any sequence");
  sq->doff = sqfp->boff;
  return eslOK;
}

static int
end_embl(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int status;

  if (strncmp(sqfp->buf, "//", 2) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Did not find // terminator at end of seq record");
  sq->eoff = sqfp->boff + sqfp->nc - 1;
  status = loadbuf(sqfp);
  if      (status == eslEOF) return eslOK; /* ok, actually. */
  else if (status == eslOK)  return eslOK;
  else                       return status;
}

/*---------------------- EMBL format ---------------------------------*/



/*****************************************************************
 * 10. Internal routines for Genbank format 
 *****************************************************************/ 
/* NCBI Genbank sequence database format.
 * See Genbank release notes; for example,
 * ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
 */

static void
config_genbank(ESL_SQFILE *sqfp)
{
  sqfp->is_linebased      = TRUE;
  sqfp->eof_is_ok         = FALSE;	/* records end with //  */
  sqfp->parse_header      = &header_genbank;
  sqfp->parse_end         = &end_genbank;
}

static void
inmap_genbank(ESL_SQFILE *sqfp, const ESL_DSQ *abc_inmap)
{
  int x;

  if (abc_inmap != NULL) {
    for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  for (x = '0'; x <= '9'; x++)
    sqfp->inmap[x] = eslDSQ_IGNORED;
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\n'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
  sqfp->inmap['/']  = eslDSQ_EOD;
}

/* header_genbank()
 * 
 * sqfp->buf is the first (LOCUS) line of the entry, or a line before
 * it (in which case we'll scan forwards to find the LOCUS line - even
 * skipping non-blank lines, because there are sometimes headers at
 * the start of Genbank files).
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
  int   toklen;
  int   status;

  /* Find LOCUS line, allowing for ignoration of a file header.  */
  if (sqfp->nc == 0) return eslEOF;
  while (strncmp(sqfp->buf, "LOCUS   ", 8) != 0) {
    if ((status = loadbuf(sqfp)) == eslEOF) return eslEOF; /* normal   */
    else if (status != eslOK) return status;                /* abnormal */
  } 
  
  s = sqfp->buf+12;
  if ((status = esl_strtok(&s, " ", &tok, &toklen)) != eslOK)
    ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to parse name on LOCUS line");
  if ((status = esl_sq_SetName(sq, tok)) != eslOK) return status;
  sq->roff = sqfp->boff;	/* record the disk offset to the LOCUS line */
  
  /* Look for ORIGIN line, parsing optional info as we go. */
  do {
    if ((status = loadbuf(sqfp)) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find ORIGIN line");

    /* Optional VERSION line is parsed as "accession". */
    if (strncmp(sqfp->buf, "VERSION   ", 10) == 0)
      {
	s = sqfp->buf+12; 
	if ((status = esl_strtok(&s, " ", &tok, &toklen)) != eslOK)
	  ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to parse VERSION line");
	if ((status = esl_sq_SetAccession(sq, tok)) != eslOK) return status;
      }

    /* Optional DEFINITION Line is parsed as "description". */
    if (strncmp(sqfp->buf, "DEFINITION ", 11) == 0)
      {
	s = sqfp->buf+12; 
	esl_strchop(s, sqfp->nc);
	if ((status = esl_sq_AppendDesc(sq, s)) != eslOK) 
	  ESL_FAIL(status, sqfp->errbuf, "Failed to parse desc on DEFINITION line");
      }
  } while (strncmp(sqfp->buf, "ORIGIN", 6) != 0);

  if (loadbuf(sqfp) != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Failed to find any sequence");
  sq->doff = sqfp->boff;
  return eslOK;
}
  
static int
end_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int status;
  if (strncmp(sqfp->buf, "//", 2) != 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Did not find // terminator at end of seq record");
  sq->eoff = sqfp->boff + sqfp->nc - 1;
  status = loadbuf(sqfp);
  if      (status == eslEOF) return eslOK; /* ok, actually; we'll detect EOF on next sq read */
  else if (status == eslOK)  return eslOK;
  else                       return status;
}
/*----------------- end Genbank format -------------------------------*/



/*****************************************************************
 * 11. Internal routines for FASTA format
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
  } else {
    for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
    for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x;
    for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x;
  }
  sqfp->inmap[' ']  = eslDSQ_IGNORED;
  sqfp->inmap['\t'] = eslDSQ_IGNORED;
  sqfp->inmap['\r'] = eslDSQ_IGNORED;	/* DOS eol compatibility */
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
  char  c      = sqfp->buf[sqfp->bpos];
  int   status = eslOK;
  void *tmp;
  int   pos;

  while (status == eslOK && isspace(c)) status = nextchar(sqfp, &c); /* skip space (including \n) */

  if (status == eslEOF) return eslEOF;

  if (status == eslOK && c == '>') {    /* accept the > */
    sq->roff = sqfp->boff + sqfp->bpos; /* store SSI record offset */
    status = nextchar(sqfp, &c);
  } else if (c != '>') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Unexpected char %c seen, expected next FASTA seq", c);
  
  while (status == eslOK && (c == '\t' || c == ' ')) status = nextchar(sqfp, &c); /* skip space */

  /* Store the name (space delimited) */
  pos = 0;
  while (status == eslOK && ! isspace(c))
    {
      sq->name[pos++] = c;
      if (pos == sq->nalloc-1) { ESL_RALLOC(sq->name, tmp, sq->nalloc*2); sq->nalloc*=2; }
      status = nextchar(sqfp, &c); 
    }
  if (pos == 0) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "No FASTA name found");
  sq->name[pos] = '\0';
  
  while (status == eslOK &&  (c == '\t' || c == ' ')) status = nextchar(sqfp, &c);   /* skip space */

  /* Store the description (end-of-line delimited) */
  pos = 0;
  while (status == eslOK && c != '\n' && c != '\r')
    {
      sq->desc[pos++] = c;
      if (pos == sq->dalloc-1) { ESL_RALLOC(sq->desc, tmp, sq->dalloc*2); sq->dalloc*= 2; }
      status = nextchar(sqfp, &c); 
    }
  sq->desc[pos] = '\0';
  
  while (status == eslOK && (c == '\n' || c == '\r')) status = nextchar(sqfp, &c); /* skip past eol (DOS \r\n, MAC \r, UNIX \n */

  if (status != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Premature EOF in parsing FASTA name/description line");
  sq->doff = sqfp->boff + sqfp->bpos;
  sqfp->prvrpl = sqfp->prvbpl = -1;
  sqfp->currpl = sqfp->curbpl = 0;
  sqfp->linenumber++;
  return eslOK;

 ERROR:
  return status;		/* eslEMEM, from failed realloc */
}
      
/* write_fasta():
 *
 * Write sequence <sq> in FASTA format to the open stream <fp>.
 * 
 * If <save_offsets> is TRUE, then store record, data, and end offsets
 * in <sq>; this ability is used by unit tests.
 *
 * Returns <eslOK> on success.
 */
static int
write_fasta(FILE *fp, ESL_SQ *sq, int save_offsets)
{
  char     buf[61];
  int64_t  pos;

  if (save_offsets) sq->roff = ftello(fp);
  fprintf(fp, ">%s", sq->name);
  if (sq->acc[0]  != 0) fprintf(fp, " %s", sq->acc);
  if (sq->desc[0] != 0) fprintf(fp, " %s", sq->desc);
  fputc('\n', fp);

  buf[60] = '\0';
  if (save_offsets) sq->doff = ftello(fp);
  for (pos = 0; pos < sq->n; pos += 60)
    {
      if (sq->dsq != NULL) esl_abc_TextizeN(sq->abc, sq->dsq+pos+1, 60, buf);
      else                 strncpy(buf, sq->seq+pos, 60);
      fprintf(fp, "%s\n", buf);
    }
  if (save_offsets) sq->eoff = ftello(fp) - 1;
  return eslOK;
}

static int 
end_fasta(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  if (sqfp->bpos < sqfp->nc) {
    if (sqfp->buf[sqfp->bpos] != '>') ESL_FAIL(eslEFORMAT, sqfp->errbuf, "Whoops, FASTA reader is corrupted");
    sq->eoff = sqfp->boff + sqfp->bpos - 1;
  } /* else, EOF, and we don't have to do anything. */
  return eslOK;
}
/*------------------- end of FASTA i/o ---------------------------*/	       




/*****************************************************************
 * 12. Functions specific to sqio <-> msa interoperation [with <msa>] 
 *****************************************************************/

#ifdef eslAUGMENT_MSA
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

#ifdef eslAUGMENT_ALPHABET
  if (sq->dsq != NULL) { 
    if ((msa = esl_msa_CreateDigital(sq->abc, 1, sq->n)) == NULL) { status = eslEMEM; goto ERROR; }
  } else 
#endif
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

#ifdef eslAUGMENT_ALPHABET
  if (sq->dsq != NULL) esl_abc_dsqcpy(sq->dsq, sq->n, msa->ax[0]);
  else
#endif
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
 * 13. Benchmark driver
 *****************************************************************/ 
/* Some current results:
 *
 * ./benchmark /misc/data0/genomes/c.elegans/genome/allWS120
 * CPU Time: 0.90u 0.06s 00:00:00.96 Elapsed: 00:00:01
 *
 * /benchmark -i /misc/data0/genomes/c.elegans/genome/allWS120
 * CPU Time: 0.41u 0.04s 00:00:00.44 Elapsed: 00:00:00
 * 
 * ./benchmark -w /misc/data0/genomes/c.elegans/genome/allWS120
 * CPU Time: 0.83u 0.05s 00:00:00.88 Elapsed: 00:00:01
 *
 * ./benchmark -2w /misc/data0/genomes/c.elegans/genome/allWS120
 * CPU Time: 3.55u 0.26s 00:00:03.80 Elapsed: 00:00:04
 *
 * Digital times are comparable (maybe a titch faster), except
 * for -d2w, which runs much faster, because rev complementation is
 * more efficient:
 *
 * ./benchmark -d2w /misc/data0/genomes/c.elegans/genome/allWS120
 * CPU Time: 2.16u 0.31s 00:00:02.47 Elapsed: 00:00:03
 */
/* gcc -std=gnu99 -O3 -fomit-frame-pointer -fstrict-aliasing -mpentiumpro -msse2 -I. -L. -o benchmark -DeslSQIO_BENCHMARK esl_sqio.c -leasel
 * ./benchmark <seqfile>
 */
#ifdef eslSQIO_BENCHMARK
#include <stdlib.h>
#include <stdio.h>


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
  { "-C",        eslARG_INT,    "100",  NULL, NULL,  NULL,  NULL, NULL, "context size for ReadWindow()",                    0 },
  { "-W",        eslARG_INT,   "1000",  NULL, NULL,  NULL,  NULL, NULL, "window size for ReadWindow()",                     0 },
  { "-2",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  "-w", NULL, "with ReadWindow(), do both strands",               0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet, not DNA",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <DNA FASTA file>";
static char banner[] = "benchmark driver for sqio module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go       = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w        = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQ        *sq       = NULL;
  ESL_SQFILE    *sqfp     = NULL;
  char          *filename = esl_opt_GetArg(go, 1);
  int            format   = eslSQFILE_FASTA;
  int            C        = esl_opt_GetInteger(go, "-C");
  int            W        = esl_opt_GetInteger(go, "-W");
  int            do_crick = esl_opt_GetBoolean(go, "-2");
  int            n        = 0;
  uint64_t       nr       = 0;

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
	  } else if (wstatus != eslOK) esl_fatal("%s", sqfp->errbuf);

	  nr += sq->W;
	}
    }
  else 
    {
      while (esl_sqio_Read(sqfp, sq) == eslOK)  { n++; nr += sq->L; esl_sq_Reuse(sq); }
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, NULL);
  printf("Read %d sequences; %lld residues.\n", n, (long long int) nr);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslSQIO_BENCHMARK*/
/*------------------ end of benchmark ---------------------------*/



/*****************************************************************
 * 14. Unit tests
 *****************************************************************/ 
#ifdef eslSQIO_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

static void
synthesize_testseqs(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int maxL, int N, ESL_SQ ***ret_sqarr)
{
  ESL_SQ **sqarr  = malloc(sizeof(ESL_SQ *) * N);
  float   *fq     = malloc(sizeof(float)   * abc->Kp);
  char    *buf    = NULL;
  int      maxn   = eslSQ_NAMECHUNK*2;
  int      maxa   = eslSQ_ACCCHUNK*2;
  int      maxd   = eslSQ_DESCCHUNK*2;
  char     ascii[128];
  float    af[128];
  int      i, pos;
  int      n;
  int      x;

  n = ESL_MAX( maxn, ESL_MAX(maxa, maxd));
  buf = malloc(sizeof(char) * (n+1));

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

      n = esl_rnd_Roll(r, maxn) + 1; /* 1..maxn */
      esl_rsq_fIID(r, ascii, af, 128, n, buf);
      buf[n] = '\0';
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
      esl_abc_TextizeN(sq->abc, sq->dsq+pos, 60, buf);
      fputs(buf, fp);
      fputc('\n', fp);
    }
  while (esl_rnd_Roll(r, 10) == 0) fputc('\n', fp);
  sq->eoff = ftello(fp) - 1;
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
      esl_abc_TextizeN(sq->abc, sq->dsq+pos, 10, buf);
      fputs(buf, fp);
      if (pos+9 >= sq->n || (pos+9) % 60 == 0) fputc('\n',  fp);
      else                                     fputc(' ', fp);
    }
  sq->eoff = ftello(fp) - 1;
}


static void
make_ssi_index(ESL_ALPHABET *abc, const char *tmpfile, int format, char *ssifile, int mode)
{ 
  char       *msg  = "sqio unit testing: failed to make SSI index";
  ESL_NEWSSI *ns   = esl_newssi_Create();
  ESL_SQFILE *sqfp = NULL;
  ESL_SQ     *sq   = esl_sq_CreateDigital(abc);
  FILE       *fp   = NULL;
  uint16_t    fh   = 0;
  int         nseq = 0;
  int         status;
 
  if (esl_newssi_AddFile(ns, tmpfile, format, &fh)              != eslOK) esl_fatal(msg);
  if (esl_sqfile_OpenDigital(abc, tmpfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      nseq++;
      if (esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->L)   != eslOK) esl_fatal(msg);
      if (sq->acc[0] != '\0' && esl_newssi_AddAlias(ns, sq->acc, sq->name) != eslOK) esl_fatal(msg);
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal(msg);
  
  if (sqfp->bpl > 0 && sqfp->rpl > 0) 
    if ((status = esl_newssi_SetSubseq(ns, fh, sqfp->bpl, sqfp->rpl)) != eslOK) esl_fatal(msg);
  
  sprintf(ssifile, "%s.ssi", tmpfile);
  if ((fp = fopen(ssifile, "wb")) == NULL)  esl_fatal(msg);
  if (esl_newssi_Write(fp, ns)   != eslOK)  esl_fatal(msg);

  switch (mode) {
  case 0:  if (sqfp->bpl != 0)                     esl_fatal(msg); break; /* uglified: bpl should be invalid (rpl might not be) */
  case 1:  if (sqfp->rpl != 60 || sqfp->bpl == 0)  esl_fatal(msg); break; /* spaced: bpl, rpl should be valid */
  case 2:  if (sqfp->rpl != 60 || sqfp->bpl != 61) esl_fatal(msg); break; /* normal: bpl, rpl should be valid, w/ bpl=rpl+1 */
  }

  fclose(fp);
  esl_sqfile_Close(sqfp);
  esl_newssi_Destroy(ns);
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

  switch (mode) {
  case 0:  if (sqfp->bpl != 0)                     esl_fatal(msg); break; /* uglified: bpl should be invalid (rpl might not be) */
  case 1:  if (sqfp->rpl != 60 || sqfp->bpl == 0)  esl_fatal(msg); break; /* spaced: bpl, rpl should be valid */
  case 2:  if (sqfp->rpl != 60 || sqfp->bpl != 61) esl_fatal(msg); break; /* normal: bpl, rpl should be valid, w/ bpl=rpl+1 */
  }

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

  switch (mode) {
  case 0:  if (sqfp->bpl != 0)                     esl_fatal(msg); break; /* uglified: bpl should be invalid (rpl might not be) */
  case 1:  if (sqfp->rpl != 60 || sqfp->bpl == 0)  esl_fatal(msg); break; /* spaced: bpl, rpl should be valid */
  case 2:  if (sqfp->rpl != 60 || sqfp->bpl != 61) esl_fatal(msg); break; /* normal: bpl, rpl should be valid, w/ bpl=rpl+1 */
  }

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

      if (W > 0) {
	/* Forward strand coord checks */
	if (sqfp->L   != nres)                            esl_fatal(msg);
	if (sq->start != nres - sq->n + 1)                esl_fatal(msg);
	if (sq->end   != nres)                            esl_fatal(msg);
	if (sq->C != 0 && sq->C != C)                     esl_fatal(msg);
	if (sq->n != sq->C+sq->W)                         esl_fatal(msg);
	if (sq->start+sq->n-1 > sqarr[nseq]->L)           esl_fatal(msg);
	if (wstatus == eslEOD && sq->L != sqfp->L)        esl_fatal(msg);
	if (memcmp(sq->dsq + 1, sqarr[nseq]->dsq + sq->start, sq->C+sq->W) != 0) esl_fatal(msg);
      } else {
	/* Reverse strand coord checks */
	if (sqfp->L    != -1)                           esl_fatal(msg);
	if (sq->start  != sq->L - nres + sq->W + sq->C) esl_fatal(msg);
	if (sq->end    != sq->L - nres + 1)             esl_fatal(msg);
	if (sq->C != 0 && sq->C != C)                   esl_fatal(msg);
	if (sq->start-sq->n+1 < 1)                      esl_fatal(msg);
	if (wstatus == eslEOD && sq->end != 1)          esl_fatal(msg);
	if (memcmp(sq->dsq + 1, rev->dsq + (sq->L - sq->start + 1), sq->C+sq->W) != 0) esl_fatal(msg);
      }
    }

  switch (mode) {
  case 0:  if (sqfp->bpl != 0)                     esl_fatal(msg); break; /* uglified: bpl should be invalid (rpl might not be) */
  case 1:  if (sqfp->rpl != 60 || sqfp->bpl == 0)  esl_fatal(msg); break; /* spaced: bpl, rpl should be valid */
  case 2:  if (sqfp->rpl != 60 || sqfp->bpl != 61) esl_fatal(msg); break; /* normal: bpl, rpl should be valid, w/ bpl=rpl+1 */
  }

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
      i = esl_rnd_Roll(r, N);
      source = sqarr[i]->name; 
      
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


/* Write the sequences out to a tmpfile in chosen <format>;
 * read them back and make sure they're the same.
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
  int         i;

  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal(msg);
  for (i = 0; i < N; i++)
    esl_sqio_Write(fp, sqarr[i], format);
  fclose(fp);

  if (esl_sqfile_OpenDigital(abc, tmpfile, format, NULL, &sqfp) != eslOK) esl_fatal(msg);
  for (i = 0; i < N; i++)
    {
      if (esl_sqio_Read(sqfp, sq) != eslOK) esl_fatal(msg);
      if (strcmp(sqarr[i]->name,   sq->name)   != 0) esl_fatal(msg);
      if (sqarr[i]->L !=  sq->L)                     esl_fatal(msg);
      if (memcmp(sqarr[i]->dsq, sq->dsq, sizeof(ESL_DSQ) * (sq->L+2)) != 0) esl_fatal(msg);
      esl_sq_Reuse(sq);
    }
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  remove(tmpfile);
}
#endif /*eslSQIO_TESTDRIVE*/
/*------------------ end, unit tests ----------------------------*/



/*****************************************************************
 * 15. Test driver.
 *****************************************************************/

/* gcc -g -Wall -I. -L. -o sqio_utest -DeslSQIO_TESTDRIVE esl_sqio.c -leasel -lm
 * ./sqio_utest
 */
#ifdef eslSQIO_TESTDRIVE
#include "esl_config.h"

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
  { "-L",        eslARG_INT,   "1000",  NULL, NULL,  NULL,  NULL, NULL, "max length of test sequences",                     0 },
  { "-N",        eslARG_INT,    "100",  NULL, NULL,  NULL,  NULL, NULL, "number of test sequences",                         0 },
  { "-r",        eslARG_NONE,   NULL,   NULL, NULL,  NULL,  NULL, NULL, "use arbitrary random number seed",                 0 },
  { "-s",        eslARG_INT,     "42",  NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sqio module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc      = esl_alphabet_Create(eslDNA); /* DNA because some chars aren't legal in IUPAC DNA */
  ESL_RANDOMNESS *r        = NULL;
  ESL_SQ        **sqarr    = NULL;
  int             maxL     = esl_opt_GetInteger(go, "-L");
  int             N        = esl_opt_GetInteger(go, "-N");
  int             i;
  int             mode;
  char            tmpfile[32];
  char            ssifile[32];
  FILE           *fp       = NULL;
  char            c;

  if (esl_opt_GetBoolean(go, "-r")) r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  
  /* Create an array of sequences we'll use for all the tests */
  synthesize_testseqs(r, abc, maxL, N, &sqarr);

  for (mode = 0; mode < 3; mode++) /* 0=ugly 1=spaced 2=normal*/
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
	  write_fasta(fp, sqarr[i], TRUE); 
	  sqarr[i]->acc[0] = c;
	}
	break;
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

  utest_write(abc, sqarr, N, eslMSAFILE_STOCKHOLM);

  for (i = 0; i < N; i++) esl_sq_Destroy(sqarr[i]);
  free(sqarr);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslSQIO_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/



/*****************************************************************
 * 16. Examples
 *****************************************************************/
#ifdef eslSQIO_EXAMPLE
/*::cexcerpt::sqio_example::begin::*/
/* compile: gcc -g -Wall -I. -o sqio_example -DeslSQIO_EXAMPLE esl_sqio.c esl_sq.c easel.c
 * run:     ./example <FASTA file>
 */
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

int
main(int argc, char **argv)
{
  ESL_SQ     *sq      = esl_sq_Create();
  ESL_SQFILE *sqfp;
  int         format  = eslSQFILE_UNKNOWN;
  char       *seqfile = argv[1];
  int         status;

  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
  {     /* use each sequence for whatever you want */
    printf("Read %s: length %ld\n", sq->name, (long) sq->L);
    esl_sq_Reuse(sq);
  }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %ld, file %s:\n%s", 
	      (long) sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  return 0;
}
/*::cexcerpt::sqio_example::end::*/
#endif /*eslSQIO_EXAMPLE*/

/* Example 2 shows how to open a file digitally, while guessing its
 * file format and its alphabet. (This is a standard idiom.)
 */
#ifdef eslSQIO_EXAMPLE2
/*::cexcerpt::sqio_example2::begin::*/
/* compile: gcc -g -Wall -I. -L. -o sqio_example2 -DeslSQIO_EXAMPLE2 esl_sqio.c -leasel -lm
 * run:     ./example <sequence file>
 */
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  { "--dna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                        0 },
  { "--rna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                        0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "example for the sqio module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET *abc       = NULL;
  ESL_SQ       *sq        = NULL;
  ESL_SQFILE   *sqfp      = NULL;
  int           format    = eslSQFILE_UNKNOWN;
  int           alphatype = eslUNKNOWN;
  char         *seqfile   = esl_opt_GetArg(go, 1);
  int           status;

  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
    if      (status == eslEAMBIGUOUS) esl_fatal("Couldn't guess alphabet from first sequence in %s", seqfile);
    else if (status == eslEFORMAT)    esl_fatal("Sequence file parse error, line %ld of file %s:\n%s\n",
						(long) sqfp->linenumber, seqfile, sqfp->errbuf);
    else if (status == eslENODATA)    esl_fatal("Sequence file %s contains no data?", seqfile);
    else if (status != eslOK)         esl_fatal("Failed to guess alphabet (error code %d)\n", status);
  }
  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);
  esl_sqfile_SetDigital(sqfp, abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
  {     /* use each sequence for whatever you want */
    printf("Read %s: length %ld\n", sq->name, (long) sq->L);
    esl_sq_Reuse(sq);
  }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %ld, file %s:\n%s", 
	      (long) sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
/*::cexcerpt::sqio_example2::end::*/
#endif /*eslSQIO_EXAMPLE2*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
