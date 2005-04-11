/* sqio.c
 * Sequence file i/o.
 * 
 * SRE, Thu Feb 17 17:45:51 2005
 * SVN $Id$
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <easel.h>
#ifdef eslAUGMENT_ALPHABET
#include <esl_alphabet.h>	/* alphabet aug adds digital sequences */
#endif 
#ifdef eslAUGMENT_MSA
#include <esl_msa.h>		/* msa aug adds ability to read MSAs   */
#endif
#include <esl_sqio.h>


static int read_fasta(ESL_SQFILE *sqfp, ESL_SQ *s);
static int write_fasta(FILE *fp, ESL_SQ *s);
static int check_buffers(FILE *fp, char *buf, int *nc, int *pos, 
			 char **s, int i, int *slen);

#ifdef eslAUGMENT_MSA /* msa module augmentation provides msa<->sqio interop */
static int extract_sq_from_msa(ESL_MSA *msa, int idx, ESL_SQ *s);
static int convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa);
#endif

/*****************************************************************
 * Object manipulation for an ESL_SQ.
 *****************************************************************/ 

/* Function:  esl_sq_Create()
 * Incept:    SRE, Thu Dec 23 11:57:00 2004 [Zaragoza]
 *
 * Purpose:   Creates a new <ESL_SQ> sequence object. 
 *            Implemented as a wrapper around <esl_sq_Inflate()>,
 *            which does everything but allocate the shell
 *            itself. Caller frees with <esl_sq_Destroy()>.
 *            
 * Args:      (void)
 *
 * Returns:   pointer to the new <ESL_SQ>.
 *
 * Throws:    NULL if allocation fails.
 */
ESL_SQ *
esl_sq_Create(void)
{
  ESL_SQ *sq;

  if ((sq = malloc(sizeof(ESL_SQ))) == NULL)    
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  /* Set up the initial allocation sizes.
   * Someday, we may want to allow an application to tune these
   * at runtime, rather than using compiletime defaults.
   */
  sq->nalloc   = eslSQ_NAMECHUNK;	
  sq->aalloc   = eslSQ_ACCCHUNK;
  sq->dalloc   = eslSQ_DESCCHUNK;
  sq->salloc   = eslSQ_SEQCHUNK; 

  sq->name     = NULL;
  sq->seq      = NULL;
  sq->dsq      = NULL;	/* digital seq input currently unimplemented         */
  sq->ss       = NULL;	/* secondary structure input currently unimplemented */
  sq->optmem   = NULL;	/* this stays NULL unless we Squeeze() the structure */

  if ((sq->name = malloc(sizeof(char) * sq->nalloc)) == NULL) 
    { esl_sq_Destroy(sq); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((sq->acc  = malloc(sizeof(char) * sq->aalloc)) == NULL) 
    { esl_sq_Destroy(sq); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((sq->desc = malloc(sizeof(char) * sq->dalloc)) == NULL) 
    { esl_sq_Destroy(sq); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((sq->seq  = malloc(sizeof(char) * sq->salloc)) == NULL) 
    { esl_sq_Destroy(sq); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }

  esl_sq_Reuse(sq);		/* this does the initialization */
  return sq;
}

/* Function:  esl_sq_Reuse()
 * Incept:    SRE, Thu Dec 23 12:23:51 2004 [Zaragoza]
 *
 * Purpose:   Given a sequence object <sq> already in use (Create()'d or
 *            Inflate()'d); reinitialize all its data, so a new seq
 *            may be read into it. This allows sequential sequence
 *            input without a lot of wasted malloc()/free() cycling.
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
  sq->n = 0;
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
  int   nlen, alen, dlen, len;
  char *name, *acc, *desc, *seq, *ss, *dsq;

  nlen = strlen(sq->name);
  alen = strlen(sq->acc);
  dlen = strlen(sq->desc);

  len = nlen + alen + dlen + sq->n + 4; 
  if (sq->ss  != NULL) len += sq->n+1;
  if (sq->dsq != NULL) len += sq->n+2;

  if ((sq->optmem = malloc(sizeof(char) * len)) == NULL) 
    ESL_ERROR(eslEMEM, "allocation failed");
  
  len  = 0;
  name = sq->optmem+len; memcpy(name, sq->name, nlen+1);  len+=nlen+1;
  acc  = sq->optmem+len; memcpy(acc,  sq->acc,  alen+1);  len+=alen+1;
  desc = sq->optmem+len; memcpy(desc, sq->desc, dlen+1);  len+=dlen+1;
  seq  = sq->optmem+len; memcpy(seq,  sq->seq,  sq->n+1); len+=sq->n+1;

  if (sq->ss != NULL)
    { ss  = sq->optmem+len; memcpy(ss,  sq->ss,  sq->n+1); len+=sq->n+1; }
  if (sq->dsq != NULL)
    { dsq = sq->optmem+len; memcpy(dsq, sq->dsq, sq->n+2); len+=sq->n+2; }

  free(sq->name); sq->nalloc = 0; sq->name = name;
  free(sq->acc);  sq->aalloc = 0; sq->acc  = acc;
  free(sq->desc); sq->dalloc = 0; sq->desc = desc;
  free(sq->seq);  sq->salloc = 0; sq->seq  = seq;
  if (sq->ss  != NULL) { free(sq->ss);  sq->ss  = ss; }
  if (sq->dsq != NULL) { free(sq->dsq); sq->dsq = dsq;}
  
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


/*****************************************************************
 * ESL_SQFILE object functions
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
 */
int
esl_sqfile_Open(char *filename, int format, char *env, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp;
  int         status;		/* return status from an ESL call */
  int         x;
  int         n;

  /* Allocate and initialize the structure to base values;
   * only format is set correctly, though.
   */
  ESL_MALLOC(sqfp, sizeof(ESL_SQFILE));
  *ret_sqfp        = NULL;
  sqfp->fp         = NULL;
  sqfp->filename   = NULL;
  sqfp->ssifile    = NULL;
  sqfp->format     = format;
  sqfp->do_gzip    = FALSE;
  sqfp->do_stdin   = FALSE;
  sqfp->errbuf[0]  = '\0';
  sqfp->inmap      = NULL;
  sqfp->buf[0]     = '\0';
  sqfp->nc         = 0;
  sqfp->pos        = 0;
  sqfp->linenumber = 1;
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
      status = esl_strdup("[STDIN]", -1, &(sqfp->filename));
      if (status != eslOK) { esl_sqfile_Close(sqfp); return status; }

      sqfp->fp       = stdin;
      sqfp->do_stdin = TRUE;
    }
  else
    {
      char *envfile;
      n = strlen(filename);  

      /* Check the current working directory first.
       */
      if ((sqfp->fp = fopen(filename, "r")) != NULL)
	{
	  status = esl_FileNewSuffix(filename, "ssi", &(sqfp->ssifile));
	  if (status != eslOK) { esl_sqfile_Close(sqfp); return status;}

	  status = esl_strdup(filename, n, &(sqfp->filename));
	  if (status != eslOK) { esl_sqfile_Close(sqfp); return status;}
	}
      /* then the env variable.
       */
      else if (env != NULL &&
	       esl_FileEnvOpen(filename, env, &(sqfp->fp), &envfile)== eslOK)
	{
	  status = esl_FileNewSuffix(envfile, "ssi", &(sqfp->ssifile));
	  if (status != eslOK) { esl_sqfile_Close(sqfp); return status;}

	  status = esl_strdup(envfile, -1, &(sqfp->filename));
	  if (status != eslOK) { esl_sqfile_Close(sqfp); return status;}

	  free(envfile);
	}
      else
	{ esl_sqfile_Close(sqfp); return eslENOTFOUND; }
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

      if ((cmd = malloc(sizeof(char) * (n+1+strlen("gzip -dc ")))) == NULL)
	{ esl_sqfile_Close(sqfp); ESL_ERROR(eslEMEM, "cmd malloc failed"); }
      sprintf(cmd, "gzip -dc %s", sqfp->filename);

      if ((sqfp->fp = popen(cmd, "r")) == NULL)
	{ esl_sqfile_Close(sqfp); return eslENOTFOUND; }

      status = esl_strdup(sqfp->filename, n, &(sqfp->filename));
      if (status != eslOK)
	{ esl_sqfile_Close(sqfp); return eslEMEM; }

      sqfp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN*/


  /* Create an appropriate default input map for sequence files.
   *   - accept anything alphabetic, case-insensitive;
   *   - ignore whitespace;
   *   - anything else is illegal.
   *
   * (Eventually, we should allow this to be set by the caller if
   *  the caller already knows the correct bio alphabet.)
   * (And eventually, we might want to abstract this away to an API
   *  somewhere, since we're partially duplicating work we're doing
   *  in alphabet.c's inmap too.)
   */
  if ((sqfp->inmap = malloc(sizeof(int) * 256)) == NULL) 
    { esl_sqfile_Close(sqfp); ESL_ERROR(eslEMEM, "inmap malloc failed"); }

  for (x = 0;   x <  256; x++) sqfp->inmap[x] = ESL_ILLEGAL_CHAR;
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x - 'A';
  for (x = 'a'; x <= 'a'; x++) sqfp->inmap[x] = x - 'a';
  sqfp->inmap[' ']  = ESL_IGNORED_CHAR;
  sqfp->inmap['\t'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\n'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\r'] = ESL_IGNORED_CHAR;	/* DOS eol compatibility */


  /* If we don't know the format yet, autodetect it now.
   */
  if (sqfp->format == eslSQFILE_UNKNOWN)
    {
      if (sqfp->do_stdin || sqfp->do_gzip) 
	{ esl_sqfile_Close(sqfp); return eslEINVAL; }

      /* UNFINISHED!! resolve how we're going to do this, w/ msa.  */

      if (sqfp->format == eslSQFILE_UNKNOWN)
	{ esl_sqfile_Close(sqfp); return eslEFORMAT; }
    }

  /* Up 'til now, everything we've done is independent of the format
   * of the sequence file - we've only opened an appropriate stream. 
   *               
   * Now we preload the first data from the file, and how we do that
   * depends on format. There are three possibilities:
   *    - character-based parsers load into a fixed-size buf using fread().
   *      (the FASTA parser, for example.)
   *    
   *    - line-based parsers would need to load into a dynamic buffer using
   *      esl_fgets(). None are implemented yet, though.
   *      
   *    - Multiple alignment files handled specially, through the msa
   *      interface; we create an MSAFILE object and copy our info into 
   *      it; sqfp->msa stays NULL until we try to Read(), so we can
   *      detect and appropriately report parsing problems.
   *      
   * Note on linenumber: character based parsers will bump the linenumber
   * after they see a \n, so they init linenumber to 1. Line-based
   * parsers bump the linenumber as they read a new line, so they should
   * init linenumber to 0. It's 1 now; set to 0 if needed.
   */              
  switch (sqfp->format) {
    /* Character based parsers that use fread();
     * load first block of data.
     */
  case eslSQFILE_FASTA:
    sqfp->nc = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
    if (ferror(sqfp->fp)) {  esl_sqfile_Close(sqfp); return eslENOTFOUND; }
    break;

#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM:
    sqfp->linenumber = 0;	/* line-oriented input */
    sqfp->afp = malloc(sizeof(ESL_MSAFILE));
    if (sqfp->afp == NULL) { esl_sqfile_Close(sqfp); return eslEMEM; }
    sqfp->afp->f          = sqfp->fp;
    sqfp->afp->fname      = sqfp->filename;
    sqfp->afp->linenumber = sqfp->linenumber;
    sqfp->afp->errbuf[0]  = '\0';
    sqfp->afp->buf        = NULL;
    sqfp->afp->buflen     = 0;
    sqfp->afp->do_gzip    = sqfp->do_gzip;
    sqfp->afp->do_stdin   = sqfp->do_stdin;
    sqfp->afp->format     = sqfp->format;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  *ret_sqfp = sqfp;
  return eslOK;
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
  if (sqfp->inmap    != NULL) free(sqfp->inmap);

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





/* Function:  esl_sq_Read()
 * Incept:    SRE, Thu Feb 17 14:24:21 2005 [St. Louis]
 *
 * Purpose:   Reads the next sequence from open sequence file <sqfp> into 
 *            <sq>. Caller provides an allocated <s>, which will be
 *            internally reallocated if its space is insufficient.
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
esl_sq_Read(ESL_SQFILE *sqfp, ESL_SQ *s)
{
  int status;

  switch (sqfp->format) {
  case eslSQFILE_FASTA: status = read_fasta(sqfp, s); break;
    
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
    extract_sq_from_msa(sqfp->msa, sqfp->idx, s);
    sqfp->idx++;
    status = eslOK;
    break;
#endif /*eslAUGMENT_MSA*/
  }

  return status;
}


/* Function:  esl_sq_Write()
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
esl_sq_Write(FILE *fp, ESL_SQ *s, int format)
{
  int status;

#ifdef eslAUGMENT_MSA
  ESL_MSA *msa;
#endif

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
    ESL_ERROR(eslEINCONCEIVABLE, "no such format");
  }

  return status;
}




/*****************************************************************
 * FASTA format i/o
 *****************************************************************/

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
 *            <esl_sio_ReadFASTA()> on an open FASTA file, but
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
  int   npos = 0;	/* position in stored name        */
  int   dpos = 0;	/* position in stored description */
  int   spos = 0;	/* position in stored sequence    */
  int   c;		/* character we're processing     */
  char *buf;		/* ptr to the input buffer        */
  int   nc;		/* number of chars in the buffer  */
  int   pos;		/* position in the buffer         */
  int  *inmap;		/* ptr to the input map           */
  char *seq;            /* ptr to the growing sequence    */
  int   state;		/* state of our FSA parser        */
  int   nsafe;          /* #bytes we can safely move in both input/storage */

  /* If we have no more data, return EOF; we're done.
   */
  if (sqfp->nc == 0 && feof(sqfp->fp)) return eslEOF;

  buf   = sqfp->buf;	/* These are some optimizations, avoid dereferences */
  nc    = sqfp->nc;
  pos   = sqfp->pos;
  inmap = sqfp->inmap;
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
    case 0:     /* START. Accept >, move to state 1. */
      if (buf[pos] == '>') { pos++; state = 1; } else return eslEFORMAT;
      break;

    case 1:     /* NAMESPACE. Switch to NAME state when we see a non-space. */
      c = buf[pos];
      if   (c != ' ' && c != '\t') state = 2; else pos++; 
      break;

    case 2:     /* NAME. Accept/store non-whitespace. Else, move to state 3. */
      while ((nsafe = check_buffers(sqfp->fp, buf, &nc, &pos, 
				    &(s->name), npos, &(s->nalloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if   (!isspace(c)) { s->name[npos++] = c;  pos++;     }
	  else { state = 3; goto NAMEDONE; }
	}
      }
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
      return eslEFORMAT; /* we ran out of data while still in a name. */
    NAMEDONE:
      break;

    case 3:   /* DESCSPACE.  Accepts non-newline space and stays; else 4. */
      c = buf[pos]; 
      if (c != ' ' && c != '\t') state = 4; else pos++;
      break;
	  
    case 4:   /* DESC. Accepts and stores up to \n; accepts \n & moves to 5. */
      while ((nsafe = check_buffers(sqfp->fp, buf, &nc, &pos,
				    &(s->desc), dpos, &(s->dalloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if      (c != '\n') 
	    { s->desc[dpos++] = c;  pos++; sqfp->linenumber++; }
	  else if (c == '\r') 
	    { pos++; } /* ignore \r part of DOS \r\n EOL */
	  else                
	    { state = 5; pos++; goto DESCDONE; }
	}
      }
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
      else return eslEFORMAT;	/* ran out of data while still in desc */
    DESCDONE:
      break;

    case 5:   /* SEQ. Accept/process according to inmap; on '>',  move to 6. */
      while ((nsafe = check_buffers(sqfp->fp, buf, &nc, &pos, 
				    &seq, spos, &(s->salloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if      (inmap[c] >= 0)                
	    { seq[spos++] = c; pos++; }
	  else if (c == '>')               
	    goto FINISH; 
	  else if (c == '\n') 
	    { pos++; sqfp->linenumber++; }
	  else if (inmap[c] == ESL_ILLEGAL_CHAR) 
	    return eslEFORMAT;
	  else
	    { pos++; } /* IGNORED_CHARs, inc. \r */
	}
      }
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
      state = 6;
      break;
    } /* end of switch over FSA states */

    if (pos == nc) {		/* reload the buffer when it empties */
      nc  = fread(buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      pos = 0;
    }
    
  } /* end of while loop waiting to reach state 6 */

 FINISH:
  /* Reset the data that we took copies of.
   */
  s->seq        = seq;
  s->n          = spos;
  sqfp->pos     = pos;
  sqfp->nc      = nc;

  /* check_buffers() was careful to leave at least one
   * free byte on the storage strings, for the NUL, so we don't
   * need to check/reallocate them again.
   */
  s->name[npos] = '\0';
  s->desc[dpos] = '\0';
  s->seq[spos]  = '\0';
  
  return eslOK;
}

/* write_fasta():
 * SRE, Fri Feb 25 16:18:45 2005
 *
 * Write a sequence <s> in FASTA format to the open stream <fp>.
 *
 * Returns <eslOK> on success.
 */
int
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



/* check_buffers()
 * 
 * Given the input fread() buffer, and the storage location for where
 * we're putting the name, acc, desc, seq. If we've run out of input buffer, 
 * fread() a new block. If we're out of storage space, reallocate
 * it by doubling. Return the minimum number of bytes we can safely
 * process before we either run out of input buffer or we run out of storage
 * space, or -1 if a realloc fails.
 *
 * This supports an idiom of
 *     while (nsafe = check_buffers()) {
 *       while (nsafe--) {
 *         process buf[pos];
 *         pos++;
 *       }        
 *     }
 *     if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
 * 
 * which avoids having to check our buffers every character.
 * 
 */
static int
check_buffers(FILE *fp, char *buf, int *nc, int *pos, 
	      char **s, int i, int *slen)
{
  int inlen, savelen;

  inlen = *nc - *pos;  	/* can read this many bytes before reloading buffer */
  if (inlen == 0)	/* if we're at the end, reload now. */
    {
      *nc   = fread(buf, sizeof(char), eslREADBUFSIZE, fp);
      *pos  = 0;
      inlen = *nc - *pos;	/* (if this is still 0, we're at EOF.) */
    }

  /* Note the -1 on savelen, which leaves room for a NUL terminator.
   */
  savelen = *slen - i - 1;	/* can save this many before realloc'ing */
  if (savelen == 0)		/* if we need to reallocate now...       */
    {				/* then double our space. */
      savelen = *slen;		
      *slen  += *slen;		
      *s = realloc(*s, sizeof(char) * *slen);
      if (*s == NULL) return -1;
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

/* Function:  esl_sq_Dealign()
 * Incept:    SRE, Thu Feb 17 15:12:26 2005 [St. Louis]
 *
 * Purpose:   Dealign string <s> in place,  by removing any characters 
 *            aligned to gaps in <aseq>. Gap characters are defined in the 
 *            string <gapstring>; for example, <-_.>. Return the
 *            unaligned length of <s> in characters. 
 *            
 *            <s> can be the same as <aseq> to dealign an aligned
 *            sequence; or <s> may be an aligned annotation string
 *            (secondary structure, surface accessibility codes).
 *           
 *            It is safe to pass a NULL <s> (an unset annotation), in 
 *            which case the function no-ops and returns 0.
 *
 *            Only available when <sqio> is augmented by <msa> module.
 *            
 * Note:      To dealign one or more annotation strings as well as the
 *            sequence itself, dealign the sequence last:
 *                n1 = esl_sq_Dealign(ss,   aseq, gapstring, alen);
 *                n2 = esl_sq_Dealign(sa,   aseq, gapstring, alen);
 *                n3 = esl_sq_Dealign(aseq, aseq, gapstring, alen);
 *            Bonus paranoia if you verify that n1 == n2 == n3, but
 *            this has to be true unless <s> is NULL.
 */
int
esl_sq_Dealign(char *s, char *aseq, char *gapstring, int alen)
{
  int apos, n;
  if (s == NULL) return 0;
  
  for (apos = 0, n = 0; apos < alen; apos++)
    if (strchr(gapstring, aseq[apos]) == NULL)
      s[n++] = s[apos];
  return n;
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
  esl_sq_Dealign(s->ss,  s->seq, "-_.", msa->alen);
  s->n = esl_sq_Dealign(s->seq, s->seq, "-_.", msa->alen);

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

  *ret_msa = NULL;

  msa = esl_msa_Create(1, sq->n);
  if ((status = esl_strdup(sq->name, -1, &(msa->sqname[0]))) != eslOK)
    { esl_msa_Destroy(msa); return status; }
  
  if (*sq->acc != '\0')
    {
      msa->sqacc = malloc(sizeof(char *) * 1);
      if (msa->sqacc == NULL) 
	{ esl_msa_Destroy(msa); ESL_ERROR(eslEMEM, "malloc failed"); }
      if ((status = esl_strdup(sq->acc, -1, &(msa->sqacc[0]))) != eslOK)
	{ esl_msa_Destroy(msa); return status; }
    }

  if (*sq->desc != '\0')
    {
      msa->sqdesc = malloc(sizeof(char *) * 1);
      if (msa->sqdesc == NULL) 
	{ esl_msa_Destroy(msa); ESL_ERROR(eslEMEM, "malloc failed"); }
      if ((status = esl_strdup(sq->desc, -1, &(msa->sqdesc[0]))) != eslOK)
	{ esl_msa_Destroy(msa); return status; }
    }

  strcpy(msa->aseq[0], sq->seq);
  
  if (sq->ss != NULL)
    {
      msa->ss = malloc(sizeof(char *) * 1);
      if (msa->ss == NULL) 
	{ esl_msa_Destroy(msa); ESL_ERROR(eslEMEM, "malloc failed"); }
      if ((status = esl_strdup(sq->ss, -1, &(msa->ss[0]))) != eslOK)
	{ esl_msa_Destroy(msa); return status; }
    }
  
  msa->alen = sq->n;
  msa->nseq = 1;

  *ret_msa = msa;
  return eslOK;
}

#endif /*eslAUGMENT_MSA*/
/*---------- end of msa <-> sqio module interop -----------------*/


/*****************************************************************
 * Miscellaneous API functions
 *****************************************************************/

/* Function:  esl_sqfile_FormatCode()
 * Incept:    SRE, Sun Feb 27 09:18:36 2005 [St. Louis]
 *
 * Purpose:   Given <fmtstring>, return format code.  For example, if
 *            <fmtstring> is "fasta", returns <eslSQFILE_FASTA>. Returns 
 *            <eslSQFILE_UNKNOWN> if <fmtstring> doesn't exactly match a 
 *            known format case-insensitively.
 *            
 *            When augmented by msa, then alignment file formats
 *            are recognized in addition to unaligned file formats.
 */
int
esl_sqfile_FormatCode(char *fmtstring)
{
  if (strcasecmp(fmtstring, "fasta")     == 0) return eslSQFILE_FASTA;
#ifdef eslAUGMENT_MSA
  if (strcasecmp(fmtstring, "stockholm") == 0) return eslMSAFILE_STOCKHOLM;
  if (strcasecmp(fmtstring, "pfam")      == 0) return eslMSAFILE_PFAM;
#endif
  return eslSQFILE_UNKNOWN;
}


/* Function:  esl_sqfile_FormatString()
 * Incept:    SRE, Sun Feb 27 09:24:04 2005 [St. Louis]
 *
 * Purpose:   Given a format code <fmt>, returns a string label for
 *            that format. For example, if <fmt> is <eslSQFILE_FASTA>,
 *            returns "FASTA". 
 *            
 *            When augmented by msa, then alignment file format codes
 *            are recognized in addition to unaligned file format codes.
 *
 * Throws:    NULL if <fmt> is unrecognized.
 *
 * Xref:      
 */
char *
esl_sqfile_FormatString(int fmt)
{
  switch (fmt) {
  case eslSQFILE_UNKNOWN:    return "unknown";
  case eslSQFILE_FASTA:      return "FASTA";
#ifdef eslAUGMENT_MSA
  case eslMSAFILE_STOCKHOLM: return "Stockholm";
  case eslMSAFILE_PFAM:      return "Pfam";
#endif
  default:
    ESL_ERROR_NULL(eslEINVAL, "No such format code");
  }
  /*NOTREACHED*/
  return NULL;
}

/* Function:  esl_sqfile_IsAlignment()
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
esl_sqfile_IsAlignment(int fmt)
{
  if (fmt >= 100) return TRUE;
  else            return FALSE;
}


/*---------- end of miscellaneous API functions -----------------*/












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
  char       *seqfile = argv[1];
  int         status;

  status = esl_sqfile_Open(seqfile, eslSQFILE_FASTA, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz");
  else if (status != eslEOK)       esl_fatal("Open failed, code %d", status);

  sq = esl_sq_Create();
  while ((status = esl_sq_Read(sqfp, sq)) == eslOK)
  {
    /* use the sequence for whatever you want */
    esl_sq_Reuse(sq);
  }
  
  if (status == eslEFORMAT)
    esl_fatal("Parse failed, line %s, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  else if (status != eslEOF)
    esl_fatal("Sequence file read failed with code %d\n", status);
  
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  return 0;
}
/*::cexcerpt::sqio_example::end::*/
#endif /*eslSQIO_EXAMPLE*/


/*****************************************************************
 * Test driver:
 * gcc -g -Wall -I. -o test -DeslSQIO_TESTDRIVE esl_sqio.c easel.c
 * ./test
 *****************************************************************/
#ifdef eslSQIO_TESTDRIVE
/*::cexcerpt::sqio_test::begin::*/
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
  char        filename[] = "tmpxxx.seq";
  int         n;

  /* Create a FASTA file containing two sequences.
   */
  if ((fp = fopen(filename, "w")) == NULL) abort();
  fprintf(fp, ">seq1 seq1's description goes here\n");
  fprintf(fp, "%s\n", seq1);
  fprintf(fp, ">seq2 seq2's description goes here\n");
  fprintf(fp, "%s\n", seq2);	  
  fclose(fp);

  /* Example of the API for opening and reading 
   * seqs from a FASTA file.
   */
  if (esl_sqfile_Open(filename, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) abort();
  sq = esl_sq_Create();

  n=0;
  while (esl_sq_Read(sqfp, sq) == eslOK)
    {
      if (n==0 && strcmp(sq->seq, seq1) != 0) abort();
      if (n==1 && strcmp(sq->seq, seq2) != 0) abort();

      n++;
      esl_sq_Reuse(sq);
    }
  esl_sqfile_Close(sqfp);
  return 0;
}
/*::cexcerpt::sqio_test::end::*/
#endif /*eslSQIO_TESTDRIVE*/
