/* sqio.c
 * Sequence file i/o.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/alphabet.h>
#include <easel/parse.h>
#include <easel/sqio.h>

static int check_buffers(FILE *fp, char *buf, int *nc, int *pos, 
			 char **s, int i, int *slen);

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
    ESL_ERROR_NULL(ESL_EMEM, "malloc failed");
  if (esl_sq_Inflate(sq) != ESL_OK) 
    { free(sq); return NULL; }

  return sq;
}


/* Function:  esl_sq_Inflate()
 * Incept:    SRE, Thu Dec 23 12:20:13 2004 [Zaragoza]
 *
 * Purpose:   Given an already allocated sequence object shell <sq> (usually
 *            on the stack), allocate and initialize its internals.
 *            Caller will free it with <esl_sq_Deflate()>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EMEM> on allocation failure.
 */
int
esl_sq_Inflate(ESL_SQ *sq)
{
  /* Set up the initial allocation sizes.
   * Someday, we may want to allow an application to tune these
   * at runtime, rather than using compiletime defaults.
   */
  sq->nalloc   = ESL_SQ_NAMECHUNK;	
  sq->aalloc   = ESL_SQ_ACCCHUNK;
  sq->dalloc   = ESL_SQ_DESCCHUNK;
  sq->salloc   = ESL_SQ_SEQCHUNK; 

  sq->name     = NULL;
  sq->seq      = NULL;
  sq->dsq      = NULL;	/* digital seq input currently unimplemented         */
  sq->ss       = NULL;	/* secondary structure input currently unimplemented */
  sq->optmem   = NULL;	/* this stays NULL unless we Squeeze() the structure */

  if ((sq->name = malloc(sizeof(char) * sq->nalloc)) == NULL) 
    { esl_sq_Deflate(sq); ESL_ERROR(ESL_EMEM, "malloc failed"); }
  if ((sq->acc  = malloc(sizeof(char) * sq->aalloc)) == NULL) 
    { esl_sq_Deflate(sq); ESL_ERROR(ESL_EMEM, "malloc failed"); }
  if ((sq->desc = malloc(sizeof(char) * sq->dalloc)) == NULL) 
    { esl_sq_Deflate(sq); ESL_ERROR(ESL_EMEM, "malloc failed"); }
  if ((sq->seq  = malloc(sizeof(char) * sq->salloc)) == NULL) 
    { esl_sq_Deflate(sq); ESL_ERROR(ESL_EMEM, "malloc failed"); }

  esl_sq_Reuse(sq);		/* this does the initialization */
  return ESL_OK;
}

/* Function:  esl_sq_Reuse()
 * Incept:    SRE, Thu Dec 23 12:23:51 2004 [Zaragoza]
 *
 * Purpose:   Given a sequence object <sq> already in use (Create()'d or
 *            Inflate()'d); reinitialize all its data, so a new seq
 *            may be read into it. This allows sequential sequence
 *            input without a lot of wasted malloc()/free() cycling.
 *
 * Returns:   <ESL_OK> on success.
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
  return ESL_OK;
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
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EMEM> if reallocation fails.
 */
int
esl_sq_Squeeze(ESL_SQ *sq)
{
  int   nlen, alen, dlen, len;
  char *name, *acc, *desc, *seq;

  nlen = strlen(sq->name);
  alen = strlen(sq->acc);
  dlen = strlen(sq->desc);
  
  len = nlen + alen + dlen + sq->n + 4;

  if ((sq->optmem = malloc(sizeof(char) * len)) == NULL) 
    ESL_ERROR(ESL_EMEM, "allocation failed");
  
  name = sq->optmem;                    memcpy(name, sq->name, nlen+1);
  acc  = sq->optmem+(nlen+1);           memcpy(acc,  sq->acc,  alen+1);
  desc = sq->optmem+(nlen+alen+2);      memcpy(desc, sq->desc, dlen+1);
  seq  = sq->optmem+(nlen+alen+dlen+3); memcpy(seq,  sq->seq,  sq->n+1);

  free(sq->name); sq->nalloc = 0; sq->name = NULL;
  free(sq->acc);  sq->aalloc = 0; sq->acc  = NULL;
  free(sq->desc); sq->dalloc = 0; sq->desc = NULL;
  free(sq->seq);  sq->salloc = 0; sq->seq  = NULL;

  sq->name = name;
  sq->acc  = acc;
  sq->desc = desc;
  sq->seq  = seq;
  
  return ESL_OK;
}




/* Function:  esl_sq_Deflate()
 * Incept:    SRE, Thu Dec 23 12:26:44 2004 [Zaragoza]
 *
 * Purpose:   Free internal memory in an Inflate()'d <sq>, 
 *            leaving its shell.
 */
void
esl_sq_Deflate(ESL_SQ *sq)
{
  if (sq->name   != NULL) { free(sq->name);   sq->name   = NULL; }
  if (sq->acc    != NULL) { free(sq->acc);    sq->acc    = NULL; }
  if (sq->desc   != NULL) { free(sq->desc);   sq->desc   = NULL; }
  if (sq->seq    != NULL) { free(sq->seq);    sq->seq    = NULL; }
  if (sq->dsq    != NULL) { free(sq->dsq);    sq->dsq    = NULL; }
  if (sq->ss     != NULL) { free(sq->ss);     sq->ss     = NULL; }
  if (sq->optmem != NULL) { free(sq->optmem); sq->optmem = NULL; }
  sq->nalloc = 0;
  sq->aalloc = 0;
  sq->dalloc = 0;
  sq->salloc = 0;
  return;
}

/* Function:  esl_sq_Destroy()
 * Incept:    SRE, Thu Dec 23 12:28:07 2004 [Zaragoza]
 *
 * Purpose:   Free a Create()'d <sq>.
 */
void
esl_sq_Destroy(ESL_SQ *sq)
{
  esl_sq_Deflate(sq);
  free(sq);
  return;
}


/* Function:  esl_sqfile_OpenFASTA()
 * Incept:    SRE, Thu Dec 23 13:14:34 2004 [Zaragoza]
 *
 * Purpose:   Opens a FASTA sequence file <seqfile>, in
 *            preparation for <esl_sqfile_ReadFASTA()>; returns
 *            ptr to the <ESL_SQFILE> object via <ret_sqfp>.
 *
 * Args:      seqfile  - name of the file to open for reading
 *            ret_sqfp - RETURN: the opened <ESL_SQFILE>
 *
 * Returns:   <ESL_OK> on success; caller is responsible for
 *            closing the <sqfp> with <esl_sqfile_Close()>.
 *            
 *            Returns <ESL_ENOTFOUND> if the file <seqfile> does not
 *            exist, or cannot be opened for reading (incorrect
 *            permissions, for example), or fread() results in
 *            a read error.
 *
 * Throws:    <ESL_EMEM> if an allocation fails.
 */
int
esl_sqfile_OpenFASTA(char *seqfile, ESL_SQFILE **ret_sqfp)
{
  ESL_SQFILE *sqfp;
  int         x;

  if ((sqfp = malloc(sizeof(ESL_SQFILE))) == NULL) 
    ESL_ERROR(ESL_EMEM, "allocation failed");

  sqfp->fp         = NULL;
  sqfp->filename   = NULL;
  sqfp->format     = ESL_SQFORMAT_FASTA;
  sqfp->do_gzip    = FALSE;
  sqfp->do_stdin   = FALSE;
  sqfp->inmap      = NULL;
  sqfp->buf[0]     = '\0';
  sqfp->nc         = 0;
  sqfp->pos        = 0;
  sqfp->linenumber = 1;

  if ((sqfp->fp = fopen(seqfile,"r")) == NULL) 
    { free(sqfp); return ESL_ENOTFOUND; }
  if ((sqfp->filename = esl_strdup(seqfile, -1)) == NULL) 
    { esl_sqfile_Close(sqfp); ESL_ERROR(ESL_EMEM, "allocation failed"); }

  /* Create an appropriate default input map for FASTA files.
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
    { esl_sqfile_Close(sqfp); ESL_ERROR(ESL_EMEM, "allocation failed"); }

  for (x = 0;   x <  256; x++) sqfp->inmap[x] = ESL_ILLEGAL_CHAR;
  for (x = 'A'; x <= 'Z'; x++) sqfp->inmap[x] = x - 'A';
  for (x = 'a'; x <= 'a'; x++) sqfp->inmap[x] = x - 'a';
  sqfp->inmap[' ']  = ESL_IGNORED_CHAR;
  sqfp->inmap['\t'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\n'] = ESL_IGNORED_CHAR;
  sqfp->inmap['\r'] = ESL_IGNORED_CHAR;	/* DOS eol compatibility */

  /* load the first block of data from the file into memory. 
   */
  sqfp->nc = fread(sqfp->buf, sizeof(char), ESL_READBUFSIZE, sqfp->fp);
  if (ferror(sqfp->fp)) {  esl_sqfile_Close(sqfp); return ESL_ENOTFOUND; }

  *ret_sqfp = sqfp;
  return ESL_OK;
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
  if (sqfp->fp       != NULL) fclose(sqfp->fp);
  if (sqfp->filename != NULL) free(sqfp->filename);
  if (sqfp->inmap    != NULL) free(sqfp->inmap);
  free(sqfp);
  return;
}


/* Function:  esl_sio_ReadFASTA()
 * Incept:    SRE, Thu Dec 23 13:57:59 2004 [Zaragoza]
 *
 * Purpose:   Given an open <sqfp> for a FASTA file; read the next 
 *            sequence into <s>. Caller is responsible for creating
 *            <s> initially; but it will be reallocated here if its space is 
 *            insufficient.
 *            
 *            \verb+sqfp->pos+ is at the first byte in the file (which
 *            must be a $>$ if it's FASTA format); or at a '$>$' 
 *            for a subsequent sequence 2..N; or at EOF, byte B+1
 *            in a B-byte file, in which case we return <ESL_EOF>.
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
 * Returns:   <ESL_OK> on success; the newly read sequence info
 *               is stored in <s>.
 *            <ESL_EOF> when there is no sequence left in the file;
 *               (including first attempt to read an empty file).
 *            <ESL_EFORMAT> if there's a problem with the format,
 *               such as an illegal character. The linenumber that 
 *               the error occurred at is in <sqfp->linenumber>, 
 *               and the line itself is <sqfp->buf>, which an 
 *               application can use to format a useful error
 *               message.
 *
 * Throws:    <ESL_EMEM> on an allocation failure, either in 
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
int
esl_sio_ReadFASTA(ESL_SQFILE *sqfp, ESL_SQ *s)
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
  if (sqfp->nc == 0 && feof(sqfp->fp)) return ESL_EOF;

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
      if (buf[pos] == '>') { pos++; state = 1; } else return ESL_EFORMAT;
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
      if (nsafe == -1) ESL_ERROR(ESL_EMEM, "realloc failed");
      return ESL_EFORMAT; /* we ran out of data while still in a name. */
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
      if (nsafe == -1) ESL_ERROR(ESL_EMEM, "realloc failed");
      else return ESL_EFORMAT;	/* ran out of data while still in desc */
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
	    return ESL_EFORMAT;
	  else
	    { pos++; } /* IGNORED_CHARs, inc. \r */
	}
      }
      if (nsafe == -1) ESL_ERROR(ESL_EMEM, "realloc failed");
      state = 6;
      break;
    } /* end of switch over FSA states */

    if (pos == nc) {		/* reload the buffer when it empties */
      nc  = fread(buf, sizeof(char), ESL_READBUFSIZE, sqfp->fp);
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
  
  return ESL_OK;
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
 *     if (nsafe == -1) ESL_ERROR(ESL_EMEM, "realloc failed");
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
      *nc   = fread(buf, sizeof(char), ESL_READBUFSIZE, fp);
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
	       

/*****************************************************************
 * Test driver:
 * gcc -g -Wall -I. -o sqio_test -DESL_SQIO_TESTDRIVE sqio.c alphabet.c easel.c
 *****************************************************************/
#ifdef ESL_SQIO_TESTDRIVE

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

  if (esl_sqfile_OpenFASTA(filename, &sqfp) != ESL_OK) abort();
  sq = esl_sq_Create();

  n=0;
  while (esl_sio_ReadFASTA(sqfp, sq) == ESL_OK)
    {
      if (n==0 && strcmp(sq->seq, seq1) != 0) abort();
      if (n==1 && strcmp(sq->seq, seq2) != 0) abort();

      n++;
      esl_sq_Reuse(sq);
    }
  esl_sqfile_Close(sqfp);
  return 0;
}
#endif /*ESL_SQIO_TESTDRIVE*/
