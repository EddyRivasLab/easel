/* esl_sqio.c
 * Sequence file i/o.
 * 
 * Sections:
 *    1. The ESL_SQ object API.
 *    2. The ESL_SQFILE object API.
 *    3. The sequence input/output API.
 *    4. Internal functions.
 *    5. Test and example code.
 * 
 * Shares remote evolutionary homology with Don Gilbert's seminal,
 * public domain ReadSeq package. Last common ancestor was 1991 or so.
 * Vestiges of that history still remain in the design. Thanks Don!
 * 
 * SRE, Thu Feb 17 17:45:51 2005
 * SVN $Id$
 */

#include <esl_config.h>

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
			 char **s, int i, int *slen);

/* Optional MSA<->sqio interoperability */
#ifdef eslAUGMENT_MSA
static int extract_sq_from_msa(ESL_MSA *msa, int idx, ESL_SQ *s);
static int convert_sq_to_msa(ESL_SQ *sq, ESL_MSA **ret_msa);
#endif




/*****************************************************************
 * 1. Routines for dealing with the ESL_SQ object.
 *****************************************************************/ 

/* Function:  esl_sq_Create()
 * Incept:    SRE, Thu Dec 23 11:57:00 2004 [Zaragoza]
 *
 * Purpose:   Creates an empty <ESL_SQ> sequence object. 
 *            Caller frees with <esl_sq_Destroy()>.
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
  int status;
  ESL_SQ *sq = NULL;

  ESL_ALLOC(sq, sizeof(ESL_SQ));

  sq->name     = NULL;
  sq->acc      = NULL;
  sq->desc     = NULL;
  sq->seq      = NULL;
  sq->ss       = NULL;	/* secondary structure input currently unimplemented */
  sq->dsq      = NULL;	/* digital seq input currently unimplemented         */
  sq->optmem   = NULL;	/* this stays NULL unless we Squeeze() the structure */
  sq->nalloc   = eslSQ_NAMECHUNK;	
  sq->aalloc   = eslSQ_ACCCHUNK;
  sq->dalloc   = eslSQ_DESCCHUNK;
  sq->salloc   = eslSQ_SEQCHUNK; 

  ESL_ALLOC(sq->name, sizeof(char) * sq->nalloc);
  ESL_ALLOC(sq->acc,  sizeof(char) * sq->aalloc);
  ESL_ALLOC(sq->desc, sizeof(char) * sq->dalloc);
  ESL_ALLOC(sq->seq,  sizeof(char) * sq->salloc);
  esl_sq_Reuse(sq);	/* initialization of sq->n, offsets, and strings */
  return sq;

 CLEANEXIT:
  esl_sq_Destroy(sq);
  return NULL;
}

/* Function:  esl_sq_CreateFrom()
 * Incept:    SRE, Wed Mar 22 09:17:04 2006 [St. Louis]
 *
 * Purpose:   Create a new <ESL_SQ> object from elemental data;
 *            this provides an interface between non-Easel code
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
 * Args:      name    -  name of the sequence
 *            seq     -  the sequence (alphabetic)
 *            desc    -  optional: description line [or NULL]
 *            acc     -  optional: accession [or NULL]
 *            ss      -  optional: secondary structure annotation [or NULL]
 *
 * Returns:   a pointer to the new object. Free with
 *            <esl_sq_Destroy()>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_SQ *
esl_sq_CreateFrom(char *name, char *seq, char *desc, char *acc, char *ss)
{
  int status;
  ESL_SQ *sq = NULL;
  int  n;

  if (name == NULL) ESL_DIE(eslEINVAL, "must provide seq name");
  if (seq  == NULL) ESL_DIE(eslEINVAL, "must provide seq");

  ESL_ALLOC(sq, sizeof(ESL_SQ));
  sq->name = sq->acc = sq->desc = sq->seq = sq->ss = sq->dsq = sq->optmem = NULL;
  
  n = strlen(name)+1;
  ESL_ALLOC(sq->name, sizeof(char) * n);
  strcpy(sq->name, name);
  sq->nalloc = n;
  
  n = strlen(seq)+1;
  ESL_ALLOC(sq->seq, sizeof(char) * n);
  strcpy(sq->seq, seq);
  sq->salloc = n;

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
      if (n != sq->salloc) ESL_DIE(eslEINVAL, "ss, seq lengths mismatch");
      ESL_ALLOC(sq->ss, sizeof(char) * n);
      strcpy(sq->ss, ss);
    } 

  sq->doff = -1;
  sq->roff = -1;
  return sq;

 CLEANEXIT: /* on failure: */		  
  esl_sq_Destroy(sq);
  return NULL;
}


/* Function:  esl_sq_Reuse()
 * Incept:    SRE, Thu Dec 23 12:23:51 2004 [Zaragoza]
 *
 * Purpose:   Given a sequence object <sq> already in use;
 *            reinitialize all its data, so a new seq
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
  int   status;
  int   nlen, alen, dlen, len;
  char *name, *acc, *desc, *seq, *ss, *dsq;

  nlen = strlen(sq->name);
  alen = strlen(sq->acc);
  dlen = strlen(sq->desc);

  len = nlen + alen + dlen + sq->n + 4; 
  if (sq->ss  != NULL) len += sq->n+1;
  if (sq->dsq != NULL) len += sq->n+2;

  ESL_ALLOC(sq->optmem, sizeof(char) * len);
  
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

 CLEANEXIT:
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
 */
int
esl_sqfile_Open(char *filename, int format, char *env, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp    = NULL;
  char       *envfile = NULL;
  int         status;		/* return status from an ESL call */
  int         n;

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
  sqfp->inmap      = NULL;
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
      if ((status = esl_strdup("[STDIN]", -1, &(sqfp->filename))) != eslOK) goto CLEANEXIT;
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
	  if ((status = esl_FileNewSuffix(filename, "ssi", &(sqfp->ssifile))) != eslOK) goto CLEANEXIT;
	  if ((status = esl_strdup(filename, n, &(sqfp->filename)))           != eslOK) goto CLEANEXIT;
	}
      /* then the env variable.
       */
      else if (env != NULL &&
	       esl_FileEnvOpen(filename, env, &(sqfp->fp), &envfile)== eslOK)
	{
	  if ((status = esl_FileNewSuffix(envfile, "ssi", &(sqfp->ssifile))) != eslOK) goto CLEANEXIT;
	  if ((status = esl_strdup(envfile, -1, &(sqfp->filename)))          != eslOK) goto CLEANEXIT;
	  free(envfile); envfile = NULL;
	}
      else
	{ status = eslENOTFOUND; goto CLEANEXIT;}
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
      if ((sqfp->fp = popen(cmd, "r")) == NULL)	{ status = eslENOTFOUND; goto CLEANEXIT; }
      if ((status = esl_strdup(sqfp->filename, n, &(sqfp->filename))) != eslOK) goto CLEANEXIT;
      sqfp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN*/


  /* Allocate the input map. config_* functions set this up later,
   * depending on format.
   */
  ESL_ALLOC(sqfp->inmap, sizeof(int) * 256);

  /* If we don't know the format yet, autodetect it now.
   */
  if (sqfp->format == eslSQFILE_UNKNOWN)
    {
      if (sqfp->do_stdin || sqfp->do_gzip)   { status = eslEINVAL;  goto CLEANEXIT; }

      sqfp->format = esl_sqio_WhatFormat(sqfp->fp);

      if (sqfp->format == eslSQFILE_UNKNOWN) { status = eslEFORMAT; goto CLEANEXIT; }
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
	  if (status == eslEOF)     { status = eslEFORMAT; goto CLEANEXIT; }
	  else if (status != eslOK) { goto CLEANEXIT; }
	}
      else
	{
	  sqfp->linenumber = 1;
	  sqfp->balloc = eslREADBUFSIZE;
	  ESL_ALLOC(sqfp->buf, sizeof(char) * sqfp->balloc);
	  sqfp->nc   = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
	  sqfp->boff = 0;
	  if (ferror(sqfp->fp)) { status = eslEFORMAT; goto CLEANEXIT; }
	}
    }

  *ret_sqfp = sqfp;
  return eslOK;

 CLEANEXIT:
  if (envfile != NULL) free(envfile);
  esl_sqfile_Close(sqfp); 
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




/*****************************************************************
 * Section 3. Sequence i/o API
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
    extract_sq_from_msa(sqfp->msa, sqfp->idx, s);
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
  if      (*buf == '>')                      fmt = eslSQFILE_FASTA;
  else if (strncmp(buf, "ID   ", 5)    == 0) fmt = eslSQFILE_EMBL;
  else if (strncmp(buf, "LOCUS   ", 8) == 0) fmt = eslSQFILE_GENBANK;
  else if (strstr(buf, "Genetic Sequence Data Bank") != NULL) fmt = eslSQFILE_GENBANK;

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


/* Function:  esl_sqio_FormatString()
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
esl_sqio_FormatString(int fmt)
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
  default:
    ESL_ERROR_NULL(eslEINVAL, "No such format code");
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
    ESL_ERROR(eslESYS, "fseeko() failed");

  if (sqfp->is_linebased)
    {
      sqfp->linenumber = 0;
      if ((status = loadline(sqfp)) != eslOK) return status;
    }
  else
    {
      sqfp->linenumber = 1;
      sqfp->nc   = fread(sqfp->buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      sqfp->boff = r_off;
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
 * Section 4. Internal routines for line-oriented parsers
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
 * (ESL_IGNORED_CHAR), report an error (ESL_ILLEGAL_CHAR), or store
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
  int n0;

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
  while (sqfp->buf[sqfp->pos] != '\0')
    {
      if (sqfp->inmap[(int) sqfp->buf[sqfp->pos]] >= 0)
	sq->seq[sq->n++] = sqfp->buf[sqfp->pos++];
      else if (sqfp->inmap[(int) sqfp->buf[sqfp->pos]] == ESL_IGNORED_CHAR)
	sqfp->pos++;
      else if (sqfp->inmap[(int) sqfp->buf[sqfp->pos]] == ESL_ILLEGAL_CHAR)
	{
	  sprintf(sqfp->errbuf, "Illegal %c in sequence", sqfp->buf[sqfp->pos]);
	  return eslEFORMAT;
	}

      /* Realloc seq as needed */
      if (sq->n == sq->salloc)
	{
	  sq->salloc += sq->salloc; /* doubling */
	  sq->seq = realloc(sq->seq, sizeof(char) * sq->salloc);
	  if (sq->seq == NULL) ESL_ERROR(eslEMEM, "realloc failed");
	}
    }

  sqfp->lastrpl = sq->n - n0;	/* remember # of residues on this line. */
  sqfp->lastbpl = sqfp->nc;     /* remember # of bytes on this line.    */
  return eslOK;			/* eslOK; eslEFORMAT, eslEMEM           */
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
  char *tok;
  int   toklen;
  int   status;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  if (toklen >= sq->nalloc) {
    sq->nalloc = toklen + eslSQ_NAMECHUNK;
    sq->name   = realloc(sq->name, sizeof(char) * sq->nalloc);
    if (sq->name == NULL) return eslEMEM;
  }
  strcpy(sq->name, tok);
  return eslOK;
}
static int
set_accession(ESL_SQ *sq, char *s, char *delim) 
{
  char *tok;
  int   toklen;
  int   status;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  if (toklen >= sq->aalloc) {
    sq->aalloc = toklen + eslSQ_ACCCHUNK;
    sq->acc    = realloc(sq->acc, sizeof(char) * sq->aalloc);
    if (sq->acc == NULL) return eslEMEM;
  }
  strcpy(sq->acc, tok);
  return eslOK;
}
static int
append_description(ESL_SQ *sq, char *s, char *delim) 
{
  char *tok;
  int   toklen;
  int   status;
  int   dlen;

  status = esl_strtok(&s, delim, &tok, &toklen);
  if (status != eslOK) return eslEFORMAT;

  dlen = strlen(sq->desc);

  if (dlen + toklen + 1 >= sq->dalloc) { /* +1 for \n */
    sq->dalloc = dlen + toklen + eslSQ_DESCCHUNK;
    sq->desc   = realloc(sq->desc, sizeof(char) * sq->dalloc);
    if (sq->desc == NULL) return eslEMEM;
  }

  if (dlen > 0) sq->desc[dlen] = '\n';
  strcpy(sq->desc + dlen + 1, tok);
  return eslOK;
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
  for (x = 0;   x <  256; x++) sqfp->inmap[x] = ESL_ILLEGAL_CHAR;
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x - 'A';
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x - 'a';
  sqfp->inmap[' ']  = ESL_IGNORED_CHAR;
  sqfp->inmap['\t'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\n'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\r'] = ESL_IGNORED_CHAR;	/* DOS eol compatibility */
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


  for (x = 0;   x <  256; x++) sqfp->inmap[x] = ESL_ILLEGAL_CHAR;
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x - 'A';
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x - 'a';
  for (x = '0'; x <= '9'; x++) sqfp->inmap[x] = ESL_IGNORED_CHAR;
  sqfp->inmap[' ']  = ESL_IGNORED_CHAR;
  sqfp->inmap['\t'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\n'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\r'] = ESL_IGNORED_CHAR;	/* DOS eol compatibility */
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
  
  for (x = 0;   x <  256; x++) sqfp->inmap[x] = ESL_ILLEGAL_CHAR;
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x - 'A';
  for (x = 'a'; x <= 'z'; x++) sqfp->inmap[x] = x - 'a';
  sqfp->inmap[' ']  = ESL_IGNORED_CHAR;
  sqfp->inmap['\t'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\n'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\r'] = ESL_IGNORED_CHAR;	/* DOS eol compatibility */
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
  int   at_linestart;	/* TRUE when we're at first char of a data line    */

  /* If we have no more data, return EOF; we're done. (a normal return)
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
				    &(s->name), npos, &(s->nalloc))) > 0) {
	while (nsafe--) {
	  c = buf[pos];
	  if   (!isspace(c)) { s->name[npos++] = c;  pos++;     }
	  else { state = 3; goto NAMEDONE; }
	}
      }
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
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
				    &(s->desc), dpos, &(s->dalloc))) > 0) {
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
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
      sprintf(sqfp->errbuf, "File ended within a description line."); 
      return eslEFORMAT;	/* ran out of data while still in desc */
    DESCDONE:
      break;

    case 5:   /* SEQ. Accept/process according to inmap; on '>',  move to 6. */
      s->doff = sqfp->boff + pos;
      sqfp->lastrpl = sqfp->lastbpl = -1;
      at_linestart  = TRUE;

      while ((nsafe = check_buffers(sqfp->fp, &(sqfp->boff), buf, &nc, &pos, 
				    &seq, spos, &(s->salloc))) > 0) {
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
	  if      (inmap[c] >= 0)                
	    { seq[spos++] = c; pos++; sqfp->lastrpl++; }
	  else if (c == '>')               
	    goto FINISH; 
	  else if (c == '\n') 	/* end of a seq line */
	    { 
	      pos++; 
	      sqfp->linenumber++; 
	      at_linestart = TRUE;
	    }
	  else if (inmap[c] == ESL_ILLEGAL_CHAR) 
	    {
	      sprintf(sqfp->errbuf, "Illegal char %c found in sequence.", c); 
	      return eslEFORMAT;
	    }
	  else
	    { pos++; } /* IGNORED_CHARs, inc. \r */

	  sqfp->lastbpl++;
	}
      }
      if (nsafe == -1) ESL_ERROR(eslEMEM, "realloc failed");
      state = 6;
      break;
    } /* end of switch over FSA states */

    if (pos == nc) {		/* reload the buffer when it empties */
      nc  = fread(buf, sizeof(char), eslREADBUFSIZE, sqfp->fp);
      sqfp->boff = ftello(sqfp->fp);
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
check_buffers(FILE *fp, off_t *boff, char *buf, int *nc, int *pos, 
	      char **s, int i, int *slen)
{
  int inlen, savelen;

  inlen = *nc - *pos;  	/* can read this many bytes before reloading buffer */
  if (inlen == 0)	/* if we're at the end, reload now. */
    {
      *nc   = fread(buf, sizeof(char), eslREADBUFSIZE, fp);
      *boff = ftello(fp);
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

  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  sq = esl_sq_Create();
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
  {
    /* use the sequence for whatever you want */
    printf("Read %12s: length %d\n", sq->name, sq->n);
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
 * Test driver:
 * gcc -g -Wall -I. -o test -DeslSQIO_TESTDRIVE esl_sqio.c easel.c
 * ./test
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
  while (esl_sqio_Read(sqfp, sq) == eslOK)
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
  return 0;
}
#endif /*eslSQIO_BENCHMARK*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
