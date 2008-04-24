/* Unaligned sequence file i/o.
 * 
 * Contents:
 *    1. The <ESL_SQFILE> object.
 *    2. Sequence input/output.
 *    3. Random sequence file access [augmentation: ssi]         
 *    4. Internal routines for all line-oriented parsers.
 *    5. Internal routines for EMBL format (including Uniprot, TrEMBL)
 *    6. Internal routines for Genbank format
 *    7. Internal routines for FASTA format
 *    8. Internal routines for sq, msa interconversion
 *    9. Benchmark driver.
 *   10 . Unit tests.
 *   11. Test driver.
 *   12. Example.
 *   13. Copyright and license.
 * 
 * Still shares vestiges of remote evolutionary homology with Don
 * Gilbert's seminal, public domain ReadSeq package, though the last
 * common ancestor was 1991 or so.  Thanks Don!
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
#include "esl_sq.h"


/* Forward declarations of stuff that's defined in "internal functions" sections. */

/* Generic functions for line-based parsers. */
static int is_blankline(char *s);
static int loadline(ESL_SQFILE *sqfp);
static int addseq(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int generic_readseq(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int set_name(ESL_SQ *sq, char *s, char *delim);
static int set_accession(ESL_SQ *sq, char *s, char *delim);
static int append_description(ESL_SQ *sq, char *s, char *delim);

/* EMBL format; also Uniprot, TrEMBL */
static void config_embl(ESL_SQFILE *sqfp);
static int  read_embl(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_embl(char *buf);

/* Genbank format; also DDBJ */
static void config_genbank(ESL_SQFILE *sqfp);
static int  read_genbank(ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  end_genbank(char *buf);

/* FASTA format: uses faster character-based fread() i/o. */
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
 *# 1. The <ESL_SQFILE> object.
 *****************************************************************/ 

/* Function:  esl_sqfile_Open()
 * Synopsis:  Open a sequence file for reading.
 * Incept:    SRE, Thu Feb 17 08:22:16 2005 [St. Louis]
 *
 * Purpose:   Open a sequence file <filename> for reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The format of the file is asserted to be <format> (for
 *            example, <eslSQFILE_FASTA>). If <format> is
 *            <eslSQFILE_UNKNOWN> then format autodetection is
 *            invoked.
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
  sqfp->sq_cache   = NULL;	/* only used if we GuessAlphabet() later */
  sqfp->afp        = NULL;
  sqfp->msa        = NULL;
  sqfp->ssifile    = NULL;
  sqfp->rpl        = -1;	/* -1=unset */
  sqfp->bpl        = -1;	/* -1=unset */
  sqfp->lastrpl    = -1;	/* -1=unset */
  sqfp->lastbpl    = 0;		/* -1=unset */
  sqfp->ssi        = NULL;

  /* Open the file. It may either be in the current directory, or in a
   * directory indicated by the <env> argument.  For normal operation
   * (no pipes from stdin, gzip), this section opens the sqfp->fp and
   * stores the filename in sqfp->filename.
   * 
   * stdin special case is handled here. fp is stdin pipe, and filename
   * is "[STDIN]".
   */
  if (strcmp(filename, "-") == 0) /* stdin */
    {
      if ((status = esl_strdup("[STDIN]", -1, &(sqfp->filename))) != eslOK) goto ERROR;
      sqfp->fp       = stdin;
      sqfp->do_stdin = TRUE;
    }
  else
    {
      /* Check the current working directory first. */
      if ((sqfp->fp = fopen(filename, "r")) != NULL) 
	{
	  if ((status = esl_strdup(filename, -1, &(sqfp->filename)))           != eslOK) goto ERROR;
	}
      /* then the env variable. */
      else if (env != NULL && esl_FileEnvOpen(filename, env, &(sqfp->fp), &envfile) == eslOK)
	{
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
    sqfp->linenumber   = 0;	/* line-oriented input */
    sqfp->is_linebased = TRUE;
    sqfp->addfirst     = FALSE;	/* no-op for msa's */
    sqfp->addend       = FALSE;	/* no-op for msa's */
    sqfp->eof_is_ok    = FALSE;	/* no-op for msa's */
    sqfp->endTest      = NULL;	/* no-op for msa's */
    if ((status = esl_msafile_Open(filename, sqfp->format, env, &(sqfp->afp))) != eslOK) goto ERROR;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  /* Preload the first line, chunk of file, or alignment into buf.
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
#endif
  if (! sqfp->do_stdin && sqfp->fp != NULL) fclose(sqfp->fp);
  if (sqfp->filename != NULL) free(sqfp->filename);
  if (sqfp->ssifile  != NULL) free(sqfp->ssifile);
  if (sqfp->buf      != NULL) free(sqfp->buf);
  if (sqfp->sq_cache != NULL) esl_sq_Destroy(sqfp->sq_cache);
#ifdef eslAUGMENT_SSI
  if (sqfp->ssi      != NULL) esl_ssi_Close(sqfp->ssi);
#endif

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
 *# 2. Sequence input/output
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
esl_sqio_Read(ESL_SQFILE *sqfp, ESL_SQ *s)
{
  int     status;
#ifdef eslAUGMENT_MSA
  ESL_SQ *tmpsq = NULL;
#endif

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

/* Function:  esl_sqio_Echo()
 * Synopsis:  Echo the next sequence record onto output stream.
 * Incept:    SRE, Wed Apr  2 16:32:21 2008 [Janelia]
 *
 * Purpose:   Echo the next sequence record in input stream <sqfp> 
 *            onto the output stream <ofp>. 
 *
 *            This largely bypasses parsing. It enables records to be
 *            fetched exactly, after positioning <sqfp> with SSI. 
 *
 *            For example, Easel only parses and stores part of EMBL,
 *            Uniprot, and Genbank sequence records, so calling
 *            <esl_sqio_Read()>, <esl_sqio_Write()> will lose a lot of
 *            annotation. <esl_sqio_Echo()> allows full records record
 *            to be retrieved. (See <esl-sfetch> in the miniapps as an
 *            example.)
 *
 * Returns:   <eslOK> on success.

 *            Returns <eslEFORMAT> if the start or end of the record
 *            cannot be identified.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqio_Echo(FILE *ofp, ESL_SQFILE *sqfp)
{
  int status;

  if (sqfp->is_linebased)
    {
      int done;
      do {
	fputs(sqfp->buf, ofp);
	status = loadline(sqfp);
	if      (status == eslEOF && sqfp->eof_is_ok) done = TRUE;
	else if (status != eslOK) ESL_FAIL(eslEFORMAT, sqfp->errbuf, "File ended prematurely before end of record was found");

	done |= sqfp->endTest(sqfp->buf);
      } while (! done);
    }
  else if (sqfp->format == eslSQFILE_FASTA)
    {
      int   state = 0;
      char *buf   = sqfp->buf;
      int   pos   = sqfp->pos;
      int   nc    = sqfp->nc;
      
      while (state != 3) {
	switch (state) {	/* mini FASTA parsing state machine */
	case 0:			/* state 0: before we've seen a >  */
	  if      (isspace(buf[pos])) pos++;
	  else if (buf[pos] == '>')   { fputc(buf[pos++], ofp); state++; }
	  else                        ESL_FAIL(eslEFORMAT, sqfp->errbuf, "no > found");
	  break;
	  
	case 1:			/* state 1: we're on the name/desc line 'til we see a \n */
	  if (buf[pos] == '\n') state++;
	  fputc(buf[pos++], ofp);
	  break;

	case 2:			/* state 2: we're in sequence 'til we see the next > */
	  if (buf[pos] == '>') state++;
	  else                 fputc(buf[pos++], ofp);
	  break;
	}
	
	if (pos == nc) {
	  if ((nc = fread(buf, sizeof(char), eslREADBUFSIZE, sqfp->fp)) == 0) state = 3;
	  pos = 0;
	}
	
	sqfp->pos = pos;
	sqfp->nc  = nc;
      }
    }
  else esl_fatal("Oops, you've got an unimplemented format, apparently");

  return eslOK;
}


/* Function:  esl_sqio_WhatFormat()
 * Synopsis:  Guess the format of an open file.
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
 * Synopsis:  Convert a string to an internal format code.
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
 *
 */
int
esl_sqio_IsAlignment(int fmt)
{
  if (fmt >= 100) return TRUE;
  else            return FALSE;
}


/* Function:  esl_sqio_Position()
 * Synopsis:  Reposition an open sequence file to an offset.
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
      sqfp->pos  = 0;
      sqfp->nc   = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      if (ferror(sqfp->fp)) { return eslESYS; }
    }
  return eslOK;
}

/* Function:  esl_sqio_Rewind()
 * Synopsis:  Rewind an open sequence file to its beginning.
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


#ifdef eslAUGMENT_ALPHABET
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
  if (esl_sqio_IsAlignment(sqfp->format)) return esl_msafile_GuessAlphabet(sqfp->afp, ret_type);

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
#endif /*eslAUGMENT_ALPHABET*/

/*--------------------- end of i/o API ----------------------------*/




/*****************************************************************
 *# 3. Fast random access in a seqfile  [with SSI augmentation]
 *****************************************************************/
#ifdef eslAUGMENT_SSI
static int
reset_repositioned_file(ESL_SQFILE *sqfp)
{
  int status;

  /* If the <sqfp> had a sequence cached (because we just did a
   * GuessAlphabet(), we have to discard it, even if it's exactly the
   * sequence we're looking for, because we might call an Echo() after
   * repositioning.
   */
  if (sqfp->sq_cache != NULL)
    {
      esl_sq_Destroy(sqfp->sq_cache);
      sqfp->sq_cache = NULL;
    }

  /* Preload the next line or chunk into our input buffer.
   * Note that <linenumber> will now become relative to first line of the record,
   * rather than the first line of the file.
   */
  if (! esl_sqio_IsAlignment(sqfp->format))
    {
      if (sqfp->is_linebased)
	{
	  sqfp->linenumber = 0;		
	  if ((status = loadline(sqfp)) != eslOK) return status;
	}
      else
	{
	  sqfp->linenumber = 1;		
	  sqfp->nc         = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
	  sqfp->pos        = 0;
	  if (ferror(sqfp->fp)) return eslEOF;
	}
    }
  return eslOK;
}

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

  /* set <ssifile>, file name */
  if (ssifile_hint == NULL) {
    if ((status = esl_strdup(sqfp->filename, -1, &(sqfp->ssifile)))           != eslOK) return status;
    if ((status = esl_strcat(&(sqfp->ssifile), -1, ".ssi", 4))                != eslOK) return status;
  } else {
    if ((status = esl_strdup(ssifile_hint, -1, &(sqfp->ssifile)))             != eslOK) return status;
  }

  /* Open the SSI file */
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
  if ((status = esl_ssi_FindName(sqfp->ssi, key, &fh, &offset)) != eslOK) return status;
  if (fseeko(sqfp->fp, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS,  "fseek failed");
  sqfp->boff = offset;
  return reset_repositioned_file(sqfp);
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
  if ((status = esl_ssi_FindNumber(sqfp->ssi, which, &fh, &offset)) != eslOK) return status;
  if (fseeko(sqfp->fp, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS,  "fseek failed");
  sqfp->boff = offset;
  return reset_repositioned_file(sqfp);
}

#endif /*eslAUGMENT_SSI*/
/*------------------ end, SSI augmentation ----------------------*/





/*****************************************************************
 * 4. Internal routines for all line-oriented parsers
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
	  if (sq->dsq != NULL)
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
      if (sq->dsq != NULL && sq->n == (sq->salloc-1))
	{
	  ESL_RALLOC(sq->dsq, tmp, sizeof(ESL_DSQ) * sq->salloc * 2);
	  sq->salloc *= 2; /* doubling */
	}
      else if (sq->n == sq->salloc)
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

  if (sq->dsq != NULL) sq->dsq[sq->n+1] = eslDSQ_SENTINEL;
  else                 sq->seq[sq->n]   = '\0';
  
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
 * 5. Internal routines for EMBL format (including Uniprot and TrEMBL)
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
 * 6. Internal routines for Genbank format 
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
 * 7. Internal routines for FASTA format
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

  if (s->dsq != NULL)  { dsq = s->dsq; seq = NULL;   } 
  else                 { dsq = NULL;   seq = s->seq; }

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
	      else if (sqfp->lastrpl != sqfp->rpl) sqfp->rpl = 0;	    /* inval*/
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
	      if (s->dsq != NULL)
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

  if (s->dsq != NULL)  dsq[spos+1] = eslDSQ_SENTINEL;
  else                 seq[spos]   = '\0';

  /* Reset the data that we took copies of. */
  if (s->dsq != NULL)  s->dsq = dsq;
  else                 s->seq = seq;
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
 * 8. Functions specific to sqio <-> msa interoperation; 
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
 * 9. Benchmark driver
 *****************************************************************/ 
/* gcc -O2 -I. -L. -o benchmark -DeslSQIO_BENCHMARK esl_sqio.c -leasel
 * ./benchmark <seqfile>
 */
#ifdef eslSQIO_BENCHMARK
#include <stdlib.h>
#include <stdio.h>
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

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
/*------------------ end of benchmark ---------------------------*/



/*****************************************************************
 * 10. Unit tests
 *****************************************************************/ 

/*------------------ end, unit tests ----------------------------*/



/*****************************************************************
 * 11. Test driver.
 *****************************************************************/
/* gcc -g -Wall -I. -L. -o testdrive -DeslSQIO_TESTDRIVE esl_sqio.c -leasel -lm
 * ./testdrive
 */
#ifdef eslSQIO_TESTDRIVE
/*::cexcerpt::sqio_test::begin::*/
#include <stdlib.h>
#include <stdio.h>
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"


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
/*------------------ end, test driver ---------------------------*/



/*****************************************************************
 * 12. Example
 *****************************************************************/
#ifdef eslSQIO_EXAMPLE
/*::cexcerpt::sqio_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslSQIO_EXAMPLE esl_sqio.c esl_sq.c easel.c
 * run:     ./example <FASTA file>
 */
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"


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
