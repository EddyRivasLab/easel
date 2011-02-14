/* An input parsing abstraction.
 *
 * Table of contents:
 *   1. ESL_BUFFER object: opening/closing.
 *   2. Manipulating an ESL_BUFFER.
 *   3. Raw access to the buffer.
 *   4. Line-based parsing.
 *   5. Token-based parsing.
 *   6. Binary (fread-like) parsing.
 *   x. Copyright and license.
 *
 * SRE, Fri Dec 31 23:28:58 2010 [Zaragoza]
 * SVN $Id$
 * SVN $URL$
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _POSIX_VERSION
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif /* _POSIX_VERSION */

static int buffer_create           (ESL_BUFFER **ret_bf);
static int buffer_init_file_mmap   (ESL_BUFFER *bf, esl_pos_t filesize);
static int buffer_init_file_slurped(ESL_BUFFER *bf, esl_pos_t filesize);
static int buffer_init_file_basic  (ESL_BUFFER *bf);

static int buffer_refill(ESL_BUFFER *bf, esl_pos_t nmin);

static int buffer_getline(ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n);

/*****************************************************************
 *# 1. ESL_BUFFER object: opening/closing.
 *****************************************************************/


/* Function:  esl_buffer_Open()
 * Synopsis:  Standard Easel idiom for opening a stream by filename.
 * Incept:    SRE, Tue Jan 25 18:10:21 2011 [Janelia]
 *
 * Purpose:   Open <filename> for parsing. Return an open
 *            <ESL_BUFFER> for it.
 *            
 *            If <filename> is '-' (a single dash character), 
 *            capture the standard input stream rather than 
 *            opening a file.
 *
 *            Else, try to find <filename> relative to the current
 *            working directory: <./filename> (that is, <filename>
 *            may be a path relative to <.>). Directory <d> is now
 *            <.>, for the rest of this explanation.
 *            
 *            If not found in current directory, then check the
 *            environment variable <envvar> for a colon-delimited list
 *            of directories. For each directory <d> in that list, try
 *            to find <d/filename>.  (That is, <filename> may be a
 *            path relative to <d>.) Use the first <d> that succeeds.
 *            
 *            If <filename> ends in <.gz>, 'open' it by running <gzip
 *            -dc d/filename 2>/dev/null>, capturing the standard
 *            output from gunzip decompression in the <ESL_BUFFER>.
 *            
 *            Otherwise, open <d/filename> as a file.
 *
 * Args:      filename  - file to open for reading; or '-' for STDIN
 *            envvar    - name of an environment variable in which to
 *                        find a colon-delimited list of directories;
 *                        or <NULL> if none.
 *            ret_bf    - RETURN: new ESL_BUFFER            
 *
 * Returns:   <eslOK> on success; <*ret_bf> is the new <ESL_BUFFER>.
 * 
 *            <eslENOTFOUND> if file isn't found or isn't readable.
 *            <eslFAIL> if gunzip fails on a .gz file, probably 
 *            because gzip executable isn't found in PATH. 
 *            
 *            On any normal error, <*ret_bf> is still returned,
 *            in an unset state, with a user-directed error message
 *            in <*ret_bf->errmsg>.
 *            
 * Throws:    <eslESYS> on system call failures (such as fread()).
 *            <eslEMEM> on allocation failure.
 */
int
esl_buffer_Open(char *filename, char *envvar, ESL_BUFFER **ret_bf)
{
  char *path = NULL;
  int   n;
  int   status;

  /* "-" => stdin */
  if (strcmp(filename, "-") == 0) 
    return esl_buffer_OpenStream(STDIN, ret_bf);

  /* else, a file. find its fully qualified path  */
  if (esl_FileExists(filename))   /* look in current working directory */
    path = filename;
  else {   	       	          /* then search directory list in envvar, if any */
    status = esl_FileEnvOpen(filename, envvar, NULL, &path); 
    if      (status == eslENOTFOUND)  { esl_buffer_OpenFile(filename, ret_bf); return status; } 
    else if (status != eslOK)         { *ret_bf = NULL;                        return status; }
    /* yeah, the esl_buffer_OpenFile() looks weird - we know the file's not there! -
     * but it's a clean way to set our error return status properly, 
     * recording the correct error message in a live ESL_BUFFER's bf->errmsg.
     * note that esl_FileEnvOpen() correctly handles envvar==NULL,
     * returning eslENOTFOUND.
     */
  }

  n = strlen(path);
  if (n > 3 && strcmp(filename+n-3, ".gz") == 0)   /* if .gz => gzip -dc */
    status = esl_buffer_OpenPipe(path, "gunzip -dc %s 2>/dev/null", ret_bf);
  else
    status = esl_buffer_OpenFile(path, ret_bf);

  free(path);
  return status;
}

/* Function:  esl_buffer_OpenFile()
 * Synopsis:  Open a file.
 * Incept:    SRE, Fri Dec 31 23:41:33 2010 [Zaragoza]
 *
 * Purpose:   Open <filename> for reading. Return an open <ESL_BUFFER> in
 *            <*ret_bf>.
 *            
 *            <filename> may be a relative path such as <subdir/foo>
 *            or a full path such as </my/dir/foo>.
 *            
 *            On a POSIX-compliant system, large files are memory 
 *            mapped, and small files are just slurped into memory.
 *
 *            On non-POSIX systems, the file is opened as a stream.
 *            On a short initial read (if the file size is smaller than
 *            the buffer page size), the file is considered to be
 *            completely slurped.
 *
 * Args:      filename  - name of (or path to) file to open
 *           *ret_bf    - RETURN: new ESL_BUFFER
 *                        
 * Returns:   <eslOK> on success; <*ret_bf> is new <ESL_BUFFER>.
 *           
 *            <eslENOTFOUND> if <filename> isn't found or isn't readable.
 *           
 *            On normal errors, a new <*ret_bf> is still returned, in
 *            an unset state, with a user-directed error message in
 *            <*ret_bf->errmsg>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int 
esl_buffer_OpenFile(char *filename, ESL_BUFFER **ret_bf)
{
  ESL_BUFFER *bf = NULL;
#ifdef _POSIX_VERSION
  struct stat fileinfo;
#endif
  esl_pos_t   filesize = -1;
  int         status;

  if ((status = buffer_create(&bf)) != eslOK) goto ERROR;
  
  if ((bf->fp = fopen(filename, "rb")) == NULL) 
    ESL_XFAIL(eslENOTFOUND, bf->errbuf, "couldn't open %s for reading", filename);  

  if ((status = esl_strdup(filename, -1, &(bf->filename))) != eslOK) goto ERROR;

  /* Try to use POSIX fstat() to get file size and optimal read size.
   * Use filesize to decide whether to slurp, mmap, or read normally. 
   * If we don't have fstat(), we'll just read normally, and pagesize
   * will be the Easel default 4096 (set in buffer_create().)
   */
#ifdef _POSIX_VERSION
  status       = fstat(fileno(bf->fp, &fileinfo));
  filesize     = fileinfo.st_size;
  bf->pagesize = fileinfo.st_blksize;
  if (bf->pagesize < 512)     bf->pagesize = 512;      /* I feel paranoid about st_blksize range not being guaranteed to be sensible */
  if (bf->pagesize > 4194304) bf->pagesize = 4194304;
#endif  

  if      (filesize != -1 && filesize <= eslBUFFER_SLURPSIZE)  
    { if ((status = buffer_init_file_slurped(bf, filesize)) != eslOK) goto ERROR; }
#ifdef _POSIX_VERSION
  else if (filesize != -1 && filesize > eslBUFFER_SLURPSIZE) 
    { if ((status = buffer_init_file_mmap(bf, filesize))    != eslOK) goto ERROR; }
#endif
  else
    { if ((status = buffer_init_file_basic(bf))             != eslOK) goto ERROR; }

  *ret_bf = bf;
  return status;

 ERROR:
  if (bf) {
    if (bf->fp)       { fclose(bf->fp);     bf->fp       = NULL; }
    if (bf->filename) { free(bf->filename); bf->filename = NULL; }
    bf->pagesize = eslBUFFER_PAGESIZE;
  }
  *ret_bf = bf;
  return status;
}


/* Function:  esl_buffer_OpenPipe()
 * Synopsis:  Open a file through a command's stdout pipe (e.g. gunzip).
 * Incept:    SRE, Fri Dec 31 23:52:40 2010 [Zaragoza]
 *
 * Purpose:   Run the command <cmdfmt> on <filename> and capture its <stdout>
 *            stream for parsing. Return the open <ESL_BUFFER> in
 *            <*ret_bf>.
 *            
 *            <cmdfmt> has a special format; it is a <printf()>-style
 *            format string with a single <%s>, where <filename> is to
 *            be substituted. An example <cmdfmt> is "gunzip -dc %s
 *            2>/dev/null".
 *            
 *            <filename> is checked for existence and read permission
 *            before a command line is constructed.

 *            <filename> may be <NULL>. In this case, <cmdfmt> is
 *            assumed to be be the complete command, and (obviously)
 *            the diagnostic check for <filename>
 *            existence/readability is skipped. This gives you some
 *            ability to skip the restricted single-argument format of
 *            <cmdfmt>.  If you need to do something fancier with a
 *            pipe, you can always open and manage it yourself and use
 *            <esl_buffer_OpenStream()>.
 *            
 *            <popen()> executes the command under </bin/sh>.
 *            
 *            The <stderr> stream of the command should almost
 *            certainly be redirected (else it will appear on output
 *            of your program) and in general it should be discarded
 *            to </dev/null>. One of the only signs of a command
 *            failure is that the command produces a "short read", of
 *            less than <bf->pagesize> (and often 0, on a complete
 *            failure, if <stderr> has been discarded).  If <stderr>
 *            is longer than the buffer's <pagesize>, we may not
 *            accurately detect error conditions. If you must capture
 *            <stderr> (for example with a <cmdfmt> like
 *            "gunzip -dc %s 2>&1") be aware that the parser may
 *            see that output as "successful" execution.
 *            
 *            The reason to pass <cmdfmt> and <filename> separately is
 *            to enable better error diagnostics. <popen()> itself
 *            tends to "succeed" whether the command or the file exist
 *            or not.  By having <filename>, we can check for its
 *            existence/readability first.
 *            
 *            The reason that error checking <popen()> isn't entirely
 *            straightforward is that we don't see the exit status of
 *            the command until we <pclose()>. We can only <pclose()>
 *            when we're done loading data from the file, and that
 *            only happens here on a short initial read. If we do get
 *            a short read, we <pclose()>, get and check the command's
 *            exit status, and return the <ESL_BUFFER> in an
 *            <eslBUFFER_ALLFILE> state with <bf->cmdline> set.
 *
 * Args:      filename - file name (or path) to plug into <cmdfmt>; or NULL 
 *                       if <cmdfmt> is complete command already
 *            cmdfmt   - command to execute (with /bin/sh) and capture
 *                       stdout from.
 *           *ret_bf   - RETURN: new ESL_BUFFER
 *
 * Returns:   <eslOK> on success, and <*ret_bf> is the new <ESL_BUFFER>.
 * 
 *            <eslENOTFOUND> if <filename> isn't found or isn't readable.
 *
 *            <eslFAIL> if the constructed command fails - which
 *            usually means that the program isn't found or isn't
 *            executable, or that the command returned nonzero
 *            (quickly, i.e. with zero or little output and a 'short
 *            read').
 *            
 *            On any normal error, the <*ret_bf> is returned (in an
 *            <eslBUFFER_UNSET> state) and <bf->errmsg> contains a
 *            user-directed error message.
 *
 * Throws:    <eslESYS> on <*sprintf()> or <fread()> failure.
 *            <eslEMEM> on allocation failure.
 *            
 *            On any exception, <*ret_bf> is NULL.
 */
int
esl_buffer_OpenPipe(char *filename, char *cmdfmt, ESL_BUFFER **ret_bf)
{
  ESL_BUFFER *bf  = NULL;
  char       *cmd = NULL;
  int         status;

  if ((status = buffer_create(&bf)) != eslOK) goto ERROR;

  if (filename && ! esl_FileExists(filename)) 
    ESL_XFAIL(eslENOTFOUND, bf->errmsg, "couldn't read file %s", filename);

  if (filename) { if (status = esl_sprintf(&cmd, cmdfmt, filename) != eslOK) goto ERROR; }
  else          { if (status = esl_strdup(cmdfmt, -1, &cmd)        != eslOK) goto ERROR; }

  if ((bf->fp = popen(cmd, "r")) == NULL) 
    ESL_XFAIL(eslNOTFOUND, bf->errmsg, "couldn't popen() the command: %s\n", cmd);

  if (            (status = esl_strdup(cmd,      -1, &(bf->cmdline)))  != eslOK) goto ERROR;
  if (filename && (status = esl_strdup(filename, -1, &(bf->filename))) != eslOK) goto ERROR;

  ESL_ALLOC(bf->mem, sizeof(char) * bf->pagesize);
  bf->balloc  = bf->pagesize;

  bf->n = fread(bf->mem, sizeof(char), bf->pagesize, bf->fp);
  /* Now check for various errors on a short read. A short read can mean:
   *    - a small file; success, and we have the whole file in one page 
   *    - popen() "succeeded" but the command failed
   *    - an fread() failure
   * Sort it out. The trick here is that we don't get the exit status
   * of <cmd> until we pclose(). So (assuming it isn't fread() itself
   * that failed) we take advantage of the fact that we can set the
   * ESL_BUFFER to a eslBUFFER_ALLFILE state on a short initial read;
   * pclose() and check command exit status. 
   * This pretty much relies on what <stderr> from <cmd> looks like;
   * it needs to either be redirected, or short enough to be a short read.
   */   
  if (bf->n < bf->pagesize) 
    {
      if (ferror(bf->fp))      ESL_XEXCEPTION(eslESYS, "fread() failed");
      if (pclose(bf->fp) != 0) ESL_XFAIL(eslFAIL, "pipe command '%s' did not succeed");

      bf->fp      = NULL;
      bf->balloc  = 0;
      bf->mode_is = eslBUFFER_ALLFILE;
    }
  else
    bf->mode_is = eslBUFFER_PIPE;

  free(cmd);
  *ret_bf = bf;
  return eslOK;

 ERROR:
  if (bf) {	/* restore state to UNSET; w/ error message in errmsg */
    if (bf->mem)      { free(bf->mem);      bf->mem      = NULL; }
    if (bf->fp)       { pclose(bf->fp);     bf->fp       = NULL; }
    if (bf->filename) { free(bf->filename); bf->filename = NULL; }
    if (bf->cmdline)  { free(bf->cmdline);  bf->cmdline  = NULL; }
    bf->n      = 0;
    bf->balloc = 0;
  }
  if (cmd) free(cmd);
  *ret_bf = bf;
  return status;
  
}

/* Function:  esl_buffer_OpenMem()
 * Synopsis:  "Open" an existing string for parsing.
 * Incept:    SRE, Mon Jan 24 08:55:20 2011 [Janelia]
 *
 * Purpose:   Given a buffer or string <p> of length <n>, turn it into
 *            an <ESL_BUFFER>. Return the new buffer in <*ret_bf>.
 *            
 *            The memory for <p> is still managed by the caller. 
 *            Caller should free it, if necessary, only after the 
 *            ESL_BUFFER has been closed. 
 *            
 *            As a special case, if <n> is -1, <p> is assumed to be a
 *            <\0>-terminated string and its length is calculated with
 *            <strlen()>.
 *
 * Args:      p      - ptr to buffer or string
 *            n      - length of buffer/string <p>, in bytes
 *            ret_bf - RETURN: new ESL_BUFFER
 *
 * Returns:   <eslOK> on success, and <*ret_bf> points to new buffer.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            On any exception, <*ret_bf> is <NULL>.
 */
int
esl_buffer_OpenMem(char *p, esl_pos_t n, ESL_BUFFER **ret_bf)
{
  ESL_BUFFER *bf = NULL;
  int         status;

  if ((status = buffer_create(&bf)) != eslOK) goto ERROR;

  bf->mem     = p;
  bf->n       = (n == -1) ? strlen(p) : n;
  bf->mode_is = eslBUFFER_STRING;

  *ret_bf = bf;
  return eslOK;

 ERROR:
  if (bf) { /* on error, restore to UNSET state */				
    bf->mem     = NULL;
    bf->n       = 0;
    bf->mode_is = eslBUFFER_UNSET;
  }
  *ret_bf = bf;
  return status;
}


/* Function:  esl_buffer_OpenStream()
 * Synopsis:  "Open" an existing stream for parsing.
 * Incept:    SRE, Mon Jan 24 11:19:23 2011 [UCSF]
 *
 * Purpose:   Given an open stream <fp> for reading, create an
 *            <ESL_BUFFER> around it.
 *            
 *            <fp> is often <stdin>, for example.
 *            
 *            The caller remains responsible for closing <fp>, if it
 *            opened it. 
 * 
 * Args:      fp     - stream to associate with new ESL_BUFFER
 *           *ret_bf - RETURN: new ESL_BUFFER.
 *
 * Returns:   <eslOK> on success, and <*ret_bf> points to a new <ESL_BUFFER>.
 *
 * Throws:    <eslESYS> : fread() failed
 *            <eslEMEM> : an allocation failed
 */
int 
esl_buffer_OpenStream(FILE *fp, ESL_BUFFER **ret_bf);
{
  ESL_BUFFER *bf = NULL;
  int         status;

  if ((status = buffer_create(&bf)) != eslOK) goto ERROR;

  bf->fp = fp;			/* a copy of <fp>; caller is still responsible for it  */

  ESL_ALLOC(bf->mem, sizeof(char) * bf->pagesize);
  bf->balloc  = bf->pagesize;

  bf->n       = fread(bf->mem, sizeof(char), bf->pagesize, bf->fp);
  if (bf->n < bf->pagesize && ferror(bf->fp))
    ESL_XEXCEPTION(eslESYS, "failed to read first chunk of stream");

  bf->mode_is = eslBUFFER_STREAM;
  *ret_bf = bf;
  return eslOK;

 ERROR:
  if (bf)
    { /* restore bf to UNSET state */
      if (bf->mem)  { free(bf->mem); bf->mem = NULL; }
      bf->fp      = NULL;
      bf->n       = 0;
      bf->balloc  = 0;
      bf->mode_is = eslBUFFER_UNSET;
    }
  *ret_bf = bf; /* NULL, if creation failed; otherwise bf->errmsg may be useful to caller */
  return status;
}

int
esl_buffer_Close(ESL_BUFFER *bf)
{
  char msg[64];
  int  status = eslOK;

  if (bf->mem) 
    {
      if (bf->mode_is != eslBUFFER_MMAP) free(bf->mem);
      if (bf->mode_is == eslBUFFER_MMAP && munmap(bf->mem, bf->n) == -1) { status = eslESYS; strcpy(msg, "munmap() failed"); }
    }	

  if (bf->fp)
    {
      if (bf->mode_is == eslBUFFER_CMDPIPE && pclose(bf->fp) == -1) { status = eslESYS; strcpy(msg, "pclose() failed"); }
      if (bf->mode_is == eslBUFFER_FILE    && fclose(bf->fp) == -1) { status = eslESYS; strcpy(msg, "fclose() failed"); }
    }

  if (bf->filename) free(bf->filename);
  if (bf->cmdline)  free(bf->cmdline);
  free(bf);

  if (status != eslOK) ESL_EXCEPTION(status, msg);
  return eslOK;

 ERROR:
  return status;
}


/* buffer_create()
 * SRE, Sun Jan 23 19:18:40 2011 [UA975 IAD->SFO]
 * 
 * Allocate a new ESL_BUFFER and return it in <*ret_bf>;
 * with all fields initialized; return <eslOK>.
 * 
 * On allocation failure, returns <eslEMEM> and <*ret_bf>
 * is <NULL>.
 */
static int
buffer_create(ESL_BUFFER **ret_bf)
{
  ESL_BUFFER *bf  = NULL;
  int         status;

  ESL_ALLOC(bf, sizeof(ESL_BUFFER));
  bf->mem        = NULL;
  bf->n          = 0;
  bf->balloc     = 0;
  bf->pos        = 0;
  bf->baseoffset = 0;
  bf->anchor     = -1;
  bf->fp         = NULL;
  bf->filename   = NULL;
  bf->cmdline    = NULL;
  bf->pagesize   = eslBUFFER_PAGESIZE;
  bf->errmsg[0]  = '\0';
  bf->mode_is    = eslBUFFER_UNSET;

  *ret_bf = bf;
  return eslOK;

 ERROR:
  if (bf) free(bf);
  *ret_bf = NULL;
  return status;
}
	    
/* buffer_init_file_mmap()
 * SRE, Sun Jan 23 19:02:34 2011 [UA975 IAD->SFO]
 * 
 * On entry, we've already opened the file;
 *   bf->fp       = open stream for reading
 *   bf->filename = name of the file
 *
 * On success, returns eslOK, and *ret_bf has:
 *  bf->mem       mmap'ed file
 *  bf->n         size of the entire file, in bytes
 *  bf->mode_is   eslBUFFER_MMAP
 * 
 * On failure, returns error code, and *ret_bf has:
 *  bf->errmsg    helpful error message
 *  (all other fields at creation defaults)
 *  
 * On exception, returns error code, and *ret_bf
 * is left with all fields at creation defaults. 
 */
static int
buffer_init_file_mmap(ESL_BUFFER *bf, esl_pos_t filesize)
{
  int          status;
  /*    mmap(addr, len,          prot,      flags,                  fd,             offset */
  bf->mem = mmap(0,    filesize, PROT_READ, MAP_FILE | MAP_PRIVATE, fileno(bf->fp), 0);
  if (bf->mem == MAP_FAILED) ESL_XEXCEPTION(eslESYS, "mmap()");

  bf->n       = filesize;
  bf->mode_is = eslBUFFER_MMAP;

  /* open fp no longer needed - close it. */
  fclose(bf->fp);
  bf->fp = NULL;
  return eslOK;

 ERROR:
  if (bf->mem != MAP_FAILED) { munmap(bf->mem, bf->n); bf->mem      = NULL; }
  bf->n       = 0;
  bf->mode_is = eslBUFFER_UNSET;
  return status;
}

int
buffer_init_file_slurped(ESL_BUFFER *bf, esl_pos_t filesize)
{
  ESL_ALLOC(bf->mem, sizeof(char) * filesize);
  bf->balloc = filesize;

  bf->n = fread(bf->mem, sizeof(char), filesize, bf->fp);
  if (bf->n < filesize)
    ESL_XEXCEPTION(eslESYS, "failed to slurp %s\n", bf->filename);

  bf->mode_is = eslBUFFER_ALLFILE;

  fclose(bf->fp);   /* open fp no longer needed - close it. */
  bf->fp = NULL;
  return eslOK;

 ERROR:
  if (bf->mem) { free(bf->mem); bf->mem = NULL; }
  return status;
}

int
buffer_init_file_basic(ESL_BUFFER *bf)
{
  int status;

  ESL_ALLOC(bf->mem, sizeof(char) * bf->pagesize);
  bf->balloc  = bf->pagesize;

  bf->n       = fread(bf->mem, sizeof(char), bf->pagesize, bf->fp);
  if (bf->n < bf->pagesize && ferror(bf->fp))
    ESL_XEXCEPTION(eslESYS, "failed to read first chunk of %s", bf->filename);

  bf->mode_is = eslBUFFER_FILE;
  return eslOK;

 ERROR:
  if (bf->mem)  { free(bf->mem); bf->mem = NULL; }
  return status;
}
/*--------------- end, ESL_BUFFER open/close --------------------*/





/*****************************************************************
 *# 2. Manipulating an ESL_BUFFER
 *****************************************************************/

/* Function:  esl_buffer_GetOffset()
 * Synopsis:  Get the current position of parser in input buffer.
 * Incept:    SRE, Mon Jan 31 13:11:23 2011 [Janelia]
 *
 * Purpose:   Returns the current offset position of the parser
 *            in the input buffer: <bf->baseoffset + bf->pos>.
 */
esl_pos_t
esl_buffer_GetOffset(ESL_BUFFER *bf)
{
  return bf->baseoffset + bf->pos;
}

/* Function:  esl_buffer_SetOffset()
 * Synopsis:  Reposition the input buffer to a new place.
 * Incept:    SRE, Mon Jan 31 13:14:09 2011 [Janelia]
 *
 * Purpose:   Set the buffer's internal state (<bf->pos>) to position
 *            <offset> in the input. Load new data into the buffer if
 *            necessary.
 * 
 *            In modes where <bf->mem> contains the whole input
 *            (ALLFILE, MMAP, STRING), this always works.
 *             
 *            In modes where we're reading a
 *            nonrewindable/nonpositionable stream (STREAM, CMDPIPE),
 *            <offset> may be at or ahead of the current position, but
 *            rewinding to an offset behind the current position only
 *            works if <offset> is within the current buffer
 *            window. If the caller knows it wants to return to some
 *            <offset> later, it should set an anchor to make sure it
 *            stays in the buffer. New data may need to be read into
 *            <bf->mem> to assure <pagesize> bytes are available. If
 *            an anchor is set, this may require reoffset and/or
 *            reallocation of <bf->mem>.
 * 
 *            FILE mode is handled as above, but additionally, if no
 *            anchor is set and <offset> is not in the current buffer,
 *            <fseeko()> is used to reposition in the open file. If
 *            <fseeko()> is unavailable (non-POSIX compliant systems),
 *            FILE mode is handled like other streams, with limited
 *            rewind ability.
 *
 * Args:      bf     - input buffer being manipulated
 *            offset - new position in the input
 *                 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <offset> is invalid, either because it 
 *               would require rewinding the (nonrewindable) stream, 
 *               or because it's beyond the end.
 *            <eslESYS> if a system call fails, such as fread().
 *            <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> if <bf> internal state is corrupt.
 */
int
esl_buffer_SetOffset(ESL_BUFFER *bf, esl_pos_t offset)
{
  int status;

  /* Case 1. We have the entire file in bf->mem (or an mmap of it);
   *         Then this is trivial; we just set bf->pos.
   */
  if (bf->mode_is == eslBUFFER_ALLFILE ||
      bf->mode_is == eslBUFFER_MMAP    || 
      bf->mode_is == eslBUFFER_STRING)
    {
      bf->baseoffset = 0;  	/* (redundant: just to assure you that state is correctly set) */
      bf->pos        = offset;
    }

  /* Case 2. We have an open stream.
   *    Then:
   *     - if offset is behind us, and in our current buffer window,
   *       rewind is always possible and trivial: set bf->pos; or 
   *     - if we're a FILE, and we're on a POSIX system with fseeko(),
   *       and there's no anchor set -- then we can fseeko() to the
   *       desired offset (no matter where it is) and 
   *       reinitialize the buffer; or
   *     - otherwise rewinding a stream is not possible, generating
   *       an <eslEINVAL> error; or
   *     - finally, the remaining possibility is that the offset is
   *       ahead of (or at) the current parser position; fread()
   *       (respecting any set anchor) until <offset> is in the
   *       current buffer window, put bf->pos on it, and call
   *       buffer_refill() to be sure that we either have at least
   *       <bf->pagesize> bytes to parse (inclusive of current pos)
   *       or the stream reaches EOF.
   */
  else if (bf->mode_is == eslBUFFER_STREAM  ||
	   bf->mode_is == eslBUFFER_CMDPIPE ||
	   bf->mode_is == eslBUFFER_FILE)
    {
      if (offset >= bf->baseoffset && offset < bf->baseoffset + bf->pos) /* offset is in our current window and behind our current pos; rewind is trivial */
	{
	  bf->pos = offset-bf->baseoffset;
	}

#ifdef _POSIX_SOURCE
      else if (bf->mode_is == eslBUFFER_FILE && bf->anchor == -1 && offset > bf->baseoffset + bf->n)
	{
	  if (fseeko(bf->fp, offset, SEEK_SET) != 0) ESL_EXCEPTION(eslEINVAL, "fseeko() failed, probably bad offset");
	  bf->baseoffset = offset;
	  bf->n          = 0;
	  bf->pos        = 0;
	  status = buffer_refill(bf, 0);
	  if      (status == eslEOF) ESL_EXCEPTION(eslEINVAL, "requested offset is beyond end of file");
	  else if (status != eslOK)  return status;
	}
#endif

      else if (offset < bf->baseoffset)                /* we've already streamed past the requested offset. */
	ESL_EXCEPTION(eslEINVAL, "can't rewind stream past base offset"); 

      else  /* offset is ahead of pos (or on it); fast forward, put bf->pos on it, reloading bf->mem as needed, respecting any anchors */
	{
	  while (offset >= bf->baseoffset + bf->n)
	    {
	      bf->pos = bf->n;
	      status  = buffer_refill(bf, 0);
	      if      (status == eslEOF) ESL_EXCEPTION(eslEINVAL, "requested offset is beyond end of stream");
	      else if (status != eslOK)  return status;
	    }
	  bf->pos = offset - bf->baseoffset;
	  status  = buffer_refill(bf, 0);
	  if (status != eslEOF || status != eslOK) return status;
	}
    }
      
  else ESL_EXCEPTION(eslEINCONCEIVABLE, "attempting to manipulate an uninitialized buffer");

  return eslOK;
 ERROR:
  return status;
}



/* Function:  esl_buffer_SetAnchor()
 * Synopsis:  Sets an anchor in an input stream.
 * Incept:    SRE, Thu Feb 10 08:54:12 2011 [Janelia]
 *
 * Purpose:   Set an anchor at byte <offset> (in input coords) in
 *            input <bf>: which means, keep everything from this byte
 *            on in buffer memory, until anchor is raised.
 *            
 *            The presence of an anchor affects new reads from <fp>;
 *            <mem[r..n-1]> are protected from overwrite, and may be
 *            moved to <mem[0..n-r-1]> as new data is read from the
 *            stream.  Anchors are only needed for input streams that
 *            we read chunkwise.  If entire input is already in <bf>,
 *            setting an anchor is a no-op.
 *            
 *            In general, the caller should remember what anchor(s) it
 *            sets, so it can raise them later with
 *            <esl_buffer_RaiseAnchor()>.
 * 
 *            Byte <offset> must be in the current buffer window. If
 *            not, an <eslEINVAL> exception is thrown.
 * 
 *            Only one anchor is active at a time. If an anchor is
 *            already set for <bf>, the most upstream one is used.
 *
 * Args:      bf     - input buffer
 *            offset - absolute position in input, <0..inputlen-1>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <offset> is not in current buffer window.
 */
int
esl_buffer_SetAnchor(ESL_BUFFER *bf, esl_pos_t offset)
{
  if (! bf->fp) return eslOK;	/* without an open stream, no-op */
  if (offset < bf->baseoffset || offset > bf->baseoffset + bf->n)
    ESL_EXCEPTION(eslEINVAL, "can't set an anchor outside current buffer");

  if (bf->anchor == -1) bf->anchor = offset-bf->offset;
  else                  bf->anchor = ESL_MIN(offset-bf->offset, bf->anchor);
  return eslOK;
}

/* Function:  esl_buffer_RaiseAnchor()
 * Synopsis:  Raise an anchor.
 * Incept:    SRE, Wed Feb  2 07:37:31 2011 [Janelia]
 *
 * Purpose:   Declare that an anchor previously set at <offset>
 *            in buffer <bf> may be raised. 
 *            
 *            <offset> is in absolute input coordinates (<0..len-1> for
 *            an input of length <len>). Because it's supposed to be
 *            anchored, this position ought to be in the current
 *            buffer window. If an anchor is in effect in <bf>, 
 *            <offset> should be at or distal to that anchor.
 *
 * Args:      bf      - input buffer
 *            offset  - absolute position in input, <0..len-1>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <offset> is outside current buffer window,
 *            or if it is proximal to the active anchor in <bf>.
 */
int
esl_buffer_RaiseAnchor(ESL_BUFFER *bf, esl_pos_t offset)
{
  if (offset < bf->baseoffset || offset > bf->baseoffset + bf->n)
    ESL_EXCEPTION(eslEINVAL, "anchor is outside current buffer window? can't happen.");
  if (bf->anchor > offset - bf->baseoffset)
    ESL_EXCEPTION(eslEINVAL, "anchor is proximal to current active anchor");

  if (bf->anchor ==  offset - bf->baseoffset) bf->anchor = -1;
  return eslOK;
}

/* buffer_refill()      
 * For current buffer position bf->pos, try to assure that
 * we have at least <bf->pagesize> bytes loaded in <bf->mem> to parse.
 *
 *  If new read won't fit (space remaining: balloc-n; read size: memreadsize)
 *  then make room for it:
 *       if no anchor:   ndel = pos
 *       else w/ anchor: ndel = anchor; anchor = 0
 *       n      -= ndel               bytes are moved
 *       pos    -= ndel
 *       if (n) move <n> bytes from <ndel> to 0.
 *       base_offset += ndel
 *       if (n + memreadsize > balloc) reallocate mem
 *       fread() a block into position mem+n
 *       n += memreadsize
 *           
 *  For example, suppose we've completely parsed the buffer (pos == n) and we have no
 *  anchor: then:
 *       ndel        = n
 *       n           = 0   (no bytes are moved)
 *       pos         = 0
 *       base_offset += n
 *       fread a new block into mem+0
 *       n += memreadsize
 *
 * Returns: <eslOK> on success.
 *          <eslEOF> if no data remain to be read. Now pos == n.
 *          
 * Throws:  <eslEMEM> if an allocation fails. 
 *          <eslESYS> if fread() fails mysteriously.
 *          <eslEINCONCEIVABLE> if internal state of <bf> is corrupt.
 */
static int 
buffer_refill(ESL_BUFFER *bf, esl_pos_t nmin)
{
  esl_pos_t readsize = ESL_MAX(nmin, bf->pagesize);
  esl_pos_t ndel;
  esl_pos_t nread;
  int       status = eslOK;

  if (! bf->fp)                    return eslOK;     /* without an open stream, no-op; buffer is already the whole file */
  if (bf->n - bf->pos >= readsize) return eslOK;     /* if we already have enough data in buffer window, no-op       w  */

  if (bf->pos > bf->n) ESL_EXCEPTION(eslEINCONCEIVABLE, "impossible position for buffer <pos>"); 
  
  if (bf->balloc - bf->n < readsize)
    {
      if (bf->anchor == -1)   ndel = bf->pos;                   
      else                  { ndel = bf->anchor; bf->anchor = 0; }
      bf->n   -= ndel;
      bf->pos -= ndel;
      if (bf->n) memmove(bf->mem, bf->mem+ndel, bf->n);
      bf->baseoffset += ndel;
    }

  if (bf->n + readsize > bf->balloc)
    {
      ESL_RALLOC(bf->mem, tmp, sizeof(char) * (bf->n + readsize));
      bf->balloc = bf->n + readsize;
    }

  nread = fread(bf->mem+bf->n, sizeof(char), readsize, bf->fp);
  if (nread == 0 && ferror(bf->fp)) ESL_EXCEPTION(eslESYS, "fread() failure");

  bf->n += nread;
  return (nread == 0 ? eslEOF : eslOK);

 ERROR:
  return status;
}
/*--------------- end, ESL_BUFFER manipulation ------------------*/




/*****************************************************************
 *# 3. Raw access to the buffer
 *****************************************************************/

/* Function:  esl_buffer_Get()
 * Synopsis:  Get a pointer into the current buffer window.
 * Incept:    SRE, Mon Jan 31 20:45:22 2011 [Casa de Gatos]
 *
 * Purpose:   Given a buffer <bf>, return a pointer to the current
 *            parsing position in <*ret_p>, and the number of valid
 *            bytes from that position in <*ret_n>.
 *            
 *            If buffer is at EOF, returns <NULL> in <*ret_p> and
 *            0 in <*ret_n>. 
 *
 *            The buffer's parsing position <bf->pos> is NOT 
 *            changed. Another <Get()> call will return exactly
 *            the same <p> and <n>. Each <Get()> call is generally
 *            followed by a <Set()> call. It's the <Set()> call
 *            that moves <bf->pos> and refills the buffer.
 *            
 *            Assumes that the buffer <bf> is correctly loaded,
 *            with either at least <pagesize> bytes after the 
 *            parser position, or near/at EOF.
 *
 * Args:      bf    - buffer to get ptr from
 *            ret_p - RETURN: pointer to current parser position, or NULL on EOF
 *            ret_n - RETURN: number of valid bytes at *ret_p, or 0 on EOF
 *
 * Returns:   <eslOK> on success;
 *            <eslEOF> if no valid bytes remain in the input.
 */
int
esl_buffer_Get(ESL_BUFFER *bf, char **ret_p, esl_pos_t *ret_n)
{
  int status;

  if (bf->pos < bf->n)
    {
      *ret_p = bf->mem + bf->pos;
      *ret_n = bf->n - bf->pos;
      status = eslOK;
    }
  else
    {
      *ret_p = NULL;
      *ret_n = 0;
      status = eslEOF;
    }
  return status;
}

/* Function:  esl_buffer_Set()
 * Synopsis:  Set position and correct state of the <ESL_BUFFER>.
 * Incept:    SRE, Sun Jan  2 11:56:00 2011 [Zaragoza]
 *
 * Purpose:   Reset the state of buffer <bf>: we were recently
 *            given a pointer <p> by an <esl_buffer_Get()> call
 *            and we parsed <nused> bytes starting at <p[0]>. 
 *            
 *            <bf->pos> is set to point at <p+nused>, and we
 *            reload the buffer (if necessary) to try to have at
 *            least <bf->pagesize> bytes of input following that
 *            position.
 *            
 *            One use is in raw parsing, where we stop parsing
 *            somewhere in the buffer:
 *               esl_buffer_Get(bf, &p, &n);
 *                 (do some stuff on p[0..n-1], using up <nused> bytes)
 *               esl_buffer_Set(bf, p, nused);
 *            This includes the case of nused=n, where we parse the
 *            whole buffer that Get() gave us, and the Set() call may
 *            be needed to load new input data before the next Get().
 *            
 *            Another use is an idiom for peeking at a token, line, or
 *            a number of bytes without moving the parser position:
 *              esl_buffer_GetLine(bf, &p, &n);
 *                (do we like what we see in p[0..n-1]? no? then put it back)
 *              esl_buffer_Set(bf, p, 0);
 *              
 *            Because it is responsible for loading new input as
 *            needed, Set() may reoffset and reallocate <mem>. If the
 *            caller wants an anchor respected, it must make sure that
 *            anchor is still in effect; i.e., a caller that is
 *            restoring state to an ESL_BUFFER should call Set()
 *            BEFORE calling RaiseAnchor().
 * 
 *            As a special case, if <p> is NULL, then <nused> is
 *            ignored, <bf->pos> is left whereever it was, and the
 *            only thing the <Set()> attempts to do is to fulfill the
 *            pagesize guarantee from the current position. If a
 *            <NULL> <p> has been returned by a Get*() call because we
 *            reached EOF, for example in some parsing loop that the
 *            EOF has broken us out of, it is safe to call
 *            <esl_buffer_Set(bf, NULL, 0)>: this is a no-op on a
 *            buffer that is at EOF.
 *
 * Args:      bf    - buffer to set 
 *            p     - pointer that previous Get() gave us
 *            nused - number of bytes we used, starting at p[0]
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslESYS> if fread() fails mysteriously.
 *            <eslEINCONCEIVABLE> if internal state of <bf> is corrupt.
 */
int
esl_buffer_Set(ESL_BUFFER *bf, char *p, esl_pos_t nused)
{
  int status;
  if (p) bf->pos = (p - bf->mem) + nused;
  status = buffer_refill(bf, 0);
  if (status == eslOK || status == eslEOF) return eslOK; 
  else return status;
}
/*----------------- end, lowest level access --------------------*/




/*****************************************************************
 *# 4. Line-based parsing
 *****************************************************************/

/* Function:  esl_buffer_GetLine()
 * Synopsis:  Get ptr to next line in buffer.
 * Incept:    SRE, Tue Feb  1 15:02:33 2011 [Janelia]
 *
 * Purpose:   Get a pointer <*opt_p> to the next line in buffer <bf>,
 *            and the length of the line in <*opt_n> (in bytes, and
 *            exclusive of newline bytes). Advance buffer position
 *            past (one) newline, putting it on the next valid data
 *            byte. Thus <p[0..n-1]> is one data line. It is not
 *            NUL-terminated.
 *            
 *            <bf>'s buffer may be re(al)located as needed, to get
 *            the whole line into the current window.
 *            
 *            Because the caller only gets a pointer into <bf>'s 
 *            internal state, no other <esl_buffer> functions 
 *            should be called until the caller is done with <p>.
 *            
 *            To peek at next line, use Set to restore <bf>'s state:
 *               esl_buffer_GetLine(bf, &p, &n);
 *               esl_buffer_Set(bf, p, 0);           
 *            
 * Args:      bf    - buffer to get line from
 *           *opt_p - optRETURN: pointer to next line
 *           *opt_n - optRETURN: line length, exclusive of newline.
 *
 * Returns:   <eslOK> on success.  <*opt_p> is a valid pointer into <bf>'s buffer,
 *            and <*opt_n> is >=0. (0 would be an empty line.)
 *
 *            <eslEOF> if there's no line (even blank).
 *            On EOF, <*opt_p> is NULL and <*opt_n> is 0.
 *
 * Throws:    <eslEMEM> if allocation fails.
 *            <eslEINVAL> if an anchoring attempt is invalid
 *            <eslESYS> if a system call such as fread() fails unexpectedly
 *            <eslEINCONCEIVABLE> if <bf> internal state is corrupt.
 */
int
esl_buffer_GetLine(ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n)
{
  int       anch_set = FALSE;
  esl_pos_t anch;
  int       status;

  /* The next line starts at offset <baseoffset + pos> */
  anch = bf->baseoffset + bf->pos;
  if ((status = esl_buffer_SetAnchor(bf, anch)) != eslOK) goto ERROR;
  anch_set = TRUE;

  if ( (status = buffer_getline(bf, opt_p, opt_n)) != eslOK) goto ERROR;  /* includes normal EOF. opt_p, opt_n set here! */

  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) goto ERROR; /* accept EOF here, we've already got our line */
  
  anch_set = FALSE;
  if ( (status = esl_buffer_RaiseAnchor(bf, anch)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (anch_set) esl_buffer_RaiseAnchor(bf, anch);
  if (opt_p)   *opt_p = NULL;
  if (opt_n)   *opt_n = 0;
  return status;
}

/* Function:  esl_buffer_FetchLine()
 * Synopsis:  Fetch next line from a buffer.
 * Incept:    SRE, Tue Feb  1 15:37:34 2011 [Janelia]
 *
 * Purpose:   Get the next line from the buffer <bf>, starting from its
 *            current position.  Return an allocated copy of it in
 *            <*ret_p>, and its length in bytes in <*ret_n>.  Advance
 *            the buffer position past (one) newline, putting it on
 *            the next valid byte. The last line in a file does not
 *            need to be terminated by a newline. The returned memory is not
 *            NUL-terminated.
 *            
 *            Caller is responsible for free'ing <*opt_p>. 
 *            
 *            Because <*ret_p> is a copy of <bf>'s internal buffer,
 *            caller may continue to manipulate <bf>, unlike
 *            <esl_buffer_GetLine()>.
 *            
 * Args:      bf      - buffer to get line from
 *            *opt_p  - optRETURN: pointer to allocated copy of next line
 *            *opt_n  - optRETURN: length of <*opt_p>
 *
 * Returns:   <eslOK> on success.  <*opt_p> is an allocated copy
 *            of next line and <*opt_n> is >=0. (0 would be an empty line
 *            terminated by newline, such as "\n".)
 *
 *            <eslEOF> if there's no line (even blank).
 *            On EOF, <*opt_p> is NULL and <*opt_n> is 0.
 *
 * Throws:    <eslEMEM> if allocation fails.
 *            <eslEINVAL> if an anchoring attempt is invalid
 *            <eslESYS> if a system call such as fread() fails unexpectedly
 *            <eslEINCONCEIVABLE> if <bf> internal state is corrupt.
 */
int 
esl_buffer_FetchLine(ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n)
{
  int       anch_set = FALSE;
  char     *p        = NULL;
  esl_pos_t anch;
  char     *ptr;
  int       n;
  int       status;

  anch = bf->baseoffset + bf->pos;
  if ((status = esl_buffer_SetAnchor(bf, anch)) != eslOK) goto ERROR;
  anch_set = TRUE;

  if ( (status = buffer_getline(bf, &ptr, &n)) != eslOK) goto ERROR; /* inc. normal EOF */
  
  ESL_ALLOC(p, sizeof(char) * n);
  memcpy(p, ptr, n);

  anch_set = FALSE;
  if  ((status = esl_buffer_RaiseAnchor(bf, anch)) != eslOK) goto ERROR;

  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) goto ERROR; /* accept EOF here, we've already got our line */

  if (opt_p) *opt_p = p; else free(p);
  if (opt_n) *opt_n = n;
  return eslOK;

 ERROR:
  if (anch_set) esl_buffer_RaiseAnchor(bf, anch);
  if (p)        free(p);
  if (opt_p)    *opt_p = NULL;
  if (opt_n)    *opt_n = 0;
  return status;
}

/* Function:  esl_buffer_FetchLineAsStr()
 * Synopsis:  Fetch next line from buffer, and NUL-terminate it.
 * Incept:    SRE, Thu Feb 10 09:22:47 2011 [Janelia]
 *
 * Purpose:   Same as <esl_buffer_FetchLine()> except the
 *            returned line is NUL-terminated and can be treated
 *            as a string.
 *
 * Args:      bf    - input buffer
 *           *opt_s - optRETURN: pointer to allocated copy of next line
 *           *opt_n - optRETURN: strlen() of <*opt_s>
 *
 * Returns:   <eslOK> on success.  <*opt_p> is an allocated copy
 *            of next line and <*opt_n> is >=0. (0 would be an empty line
 *            terminated by newline, such as "\n".)
 *
 *            <eslEOF> if there's no line (even blank).
 *            On EOF, <*opt_p> is NULL and <*opt_n> is 0.
 *
 * Throws:    <eslEMEM> if allocation fails.
 *            <eslEINVAL> if an anchoring attempt is invalid
 *            <eslESYS> if a system call such as fread() fails unexpectedly
 *            <eslEINCONCEIVABLE> if <bf> internal state is corrupt.
 */
int 
esl_buffer_FetchLineAsStr(ESL_BUFFER *bf, char **opt_s, esl_pos_t *opt_n)
{
  int       anch_set = FALSE;
  char     *s        = NULL;
  esl_pos_t anch;
  char     *ptr;
  int       n;
  int       status;
	  
  anch = bf->baseoffset + bf->pos;
  if ( (status = esl_buffer_SetAnchor(bf, anch)) != eslOK) goto ERROR;
  anch_set = TRUE;

  if ( (status = buffer_getline(bf, &ptr, &n)) != eslOK) goto ERROR; /* inc normal EOF */

  ESL_ALLOC(s, sizeof(char) * (n+1));
  memcpy(s, ptr, n);
  s[n] = '\0';

  anch_set = FALSE;
  if ( (status = esl_buffer_RaiseAnchor(bf, anch)) != eslOK) goto ERROR; 

  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) goto ERROR; /* accept EOF here, we've already got our line */

  if (opt_s) *opt_s = s; else free(s);
  if (opt_n) *opt_n = n;
  return eslOK;

 ERROR:
  if (anch_set) esl_buffer_RaiseAnchor(bf, anch);
  if (s)        free(s);
  if (opt_p)    *opt_p = NULL;
  if (opt_n)    *opt_n = 0;
  return status;

}

/* buffer_getline()
 * The guts of esl_buffer_GetLine().
 * Get a pointer to next line in <bf>'s buffer, and length
 * of the line.
 * This is pulled out on its own so it can also be used in
 * FetchLine and FetchLineAsStr(). They differ in the order
 * of calling RaiseAnchor() and buffer_refill(): if we're
 * fetching, we don't need to keep the memory of the current
 * line in the buffer's window, so we can raise the anchor
 * before refilling the buffer, which can save some memory.
 * 
 * Returns: <eslOK> on success.
 *          <*opt_p> is a valid pointer into <bf>'s buffer,
 *          and <*opt_n> is >=0. (0 would be an empty line.)
 *
 *          <eslEOF> if there's no line (even blank).
 *          On EOF, <*opt_p> is NULL and <*opt_n> is 0.
 *
 * Throws:  <eslEMEM> if allocation fails.
 *          <eslESYS> if a system call such as fread() fails unexpectedly.
 *          <eslEINCONCEIVABLE> if <bf> internal state is corrupt.
 */
static int
buffer_getline(ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n)
{
  esl_pos_t nline;
  int       nterm;
  esl_pos_t nchecked;
  int       status;

  /* Finding the end of line may require buffer to be expanded */
  nchecked = 0;
  status = esl_memnewline(bf->mem + bf->pos, bf->n - bf->pos, &nline, &nterm);
  if (status != eslOK || status != eslEOD) goto ERROR;
  while (nterm == 0)
    {
      nchecked = bf->n - bf->pos;

      status = buffer_refill(bf, (bf->n - bf->pos) + bf->pagesize);
      if (status == eslEOF) break;
      if (status != eslOK)  goto ERROR;

      status = esl_memnewline(bf->mem + bf->pos + nchecked, bf->n - bf->pos - nchecked, &nline, &nterm);
      if (status != eslOK || status != eslEOD) goto ERROR;
    }
  nline += nchecked;

  /* EOF check. If we get here with status == eslEOF, nline = nterm = 0,
   * then esl_memnewline found no data for a line and buffer_refill returned EOF.
   */
  if (status == eslEOF && nline == 0 && nterm == 0) goto ERROR; 

  /* Beware a slippery little danger case: For a Windows newline
   * "\r\n", we may have only seen a '\r' at the end of our current
   * buffer window; esl_memnewline() thinks a bare \r is an old MacOS
   * newline; then when we refill the buffer, we'll think the \n is a
   * UNIX newline. So we'd read two newlines where there's one.  The
   * code below should handle the general case, of any newline string
   * that's a prefix of another, for up to two-char newline strings.
   */
  if (bf->pos + nline + nterm == bf->n && nterm == 1)
    {
      status = buffer_refill(bf, (bf->n - bf->pos) + bf->pagesize);
      if (status != eslOK && status != eslEOF) goto ERROR;
      
      status = esl_memnewline(bf->mem+bf->pos+nline, bf->n-bf->pos-nline, &nchecked, &nterm);
      if (status != eslOK || status != eslEOD) goto ERROR;
    }

   /* Now the line extends from bf->pos..pos+nline-1 in the buffer window 
    * The newline itself is position pos+nline..pos+nline+nterm-1 (one or two chars) 
    * The next line starts at pos+nline+nterm
    * We know that everything up to pos+nline+nterm-1 is in the buffer window.
    * It's safe to set the current buffer position to pos+nline+nterm, which
    * may be bf->n; the next buffer_refill() will do the right thing with it.
    */
  if (opt_p) *opt_p = bf->mem + bf->pos;
  if (opt_n) *opt_n = nline;
  bf->pos += nline + nterm;
  return eslOK;

 ERROR: /* including normal EOF */
  if (opt_p) *opt_p = NULL;
  if (opt_n) *opt_n = 0;
  return status;
}
/*------------------ end, line-based parsing --------------------*/

/*****************************************************************
 *# 5. Token-based parsing
 *****************************************************************/

/* Token-based parsing using an ESL_BUFFER:
 *   int esl_buffer_GetToken (ESL_BUFFER *bf, char *delim, char **opt_tok, esl_pos_t *opt_n)
 *   int esl_buffer_PeekToken(ESL_BUFFER *bf, char *delim, char **opt_tok, esl_pos_t *opt_n)
 *       skip chars in *delim
 *       set anchor
 *          if anchor is on a newline \n or \r, skip \n\r and return EOL
 *       skip chars not in *delim
 *       set *ret_tok, *opt_token
 *       skip chars in *delim
 *       set b->pos to next non-*delim character.
 * (GetLine is the same, with first delim=="" and second delim == "\n\r")
 * 
 * We want to be able to handle newlines in two ways.
 * 
 * 1) We stop at newlines - we know how many tokens we expect on the line.
 * 2) We skip over newlines, looking for next token. The number of
 *    tokens per line isn't important, but the number of tokens per
 *    record or file is.
 *
 * In the first case, we pass a *delim that does not include newline chars.
 * The first *delim skip leaves mem[pos] at the first non-*delim char, 
 * which could be a newline. If we reach EOF, return EOF. If mem[pos]
 * is a newline char, skip forward to start of next line, no token
 * is found, and return EOL. Else set the token. Then the second *delim
 * skips forward. mem[pos] is left at the start of the next non-*delim
 * char, which could be a newline on the current line if no tokens remain.
 * 
 * In the second case, we pass a *delim that includes newline chars: "\n\r".
 * Now, if no token is found on the current line, the first *delim skips
 * newlines looking for the next non-*delim char, starting a token.
 * If it reaches EOF first, return EOF. Else, set the token. Then the
 * second *delim skips forward, including newlines. mem[pos] is set to the start 
 * of the next non-*delim char, which might be on another line. 
 * 
 * To peek at the next token:
 *    esl_buffer_GetToken(bf, delim, &p, &n);
 *    esl_buffer_Set(bf, p, 0);
 */


/* Function:  esl_buffer_GetToken()
 * Synopsis:  
 * Incept:    SRE, Sat Jan  1 19:42:44 2011 [Janelia]
 *
 * Purpose:   A 'token' consists of one or more characters that are
 *            neither in <*delim> nor a newline ('\r' or '\n').
 *
 * Args:      
 *
 * Returns:   <eslOK> if a token is found; <*ret_p> points to it,
 *            and <*opt_n> is its length in chars (> 0). 
 *            
 *            <eslEOF> if the input ends before a token is found.
 *            Now <*ret_p> is <NULL> and <*opt_n> is 0.
 *            
 *            <eslEOL> if a line ends before a token is found.  (This
 *            case only arises if *delim does not contain newline.)
 *            Now <*ret_p> is <NULL> and <*opt_n> is 0.
 *            
 *            <bf->mem[bf->pos]> is left at the next
 *            non-delim/non-newline character; or EOF if no such
 *            character remains in the input.
 *            
 *            <bf->mem> may be modified and/or reallocated, if new
 *            input reads are required to find the entire token.
 *
 * Throws:    <eslEMEM> if an allocation fails.
 *            Now <*ret_p> is <NULL> and <*opt_n> is 0.
 */
int
esl_buffer_GetToken(ESL_BUFFER *bf, char *delim, char **opt_p, esl_pos_t *opt_n)
{
  esl_pos_t anch;
  esl_pos_t nc, nskip;
  int       anch_set = FALSE;
  int       status;

  if ( (status = buffer_skipdelim(bf, delim)) != eslOK) goto ERROR; /* includes EOF */
  /* Now bf->pos is sitting on first char after delims */

  if ( ( status = buffer_newline(bf)) != eslOK) goto ERROR; /* includes EOL */

  anch = bf->baseoffset + bf->pos;
  if ( (status = esl_buffer_SetAnchor(bf, anchor)) != eslOK) goto ERROR;
  anch_set = TRUE;

  if ( (status = buffer_counttok(bf, delim, &nc, &nskip)) != eslOK) goto ERROR;
  /* now we know that pos..pos+nc-1 is a token; pos+nskip is the next parser position */

  /* refill the buffer */
  status = buffer_refill(bf, (bf->n-bf->pos) + bf->pagesize);
  if (status != eslEOF && status != eslOK) goto ERROR; /* accept EOF here, we've already got our line */

  if ( ( status = buffer_newline(bf)) != eslOK) goto ERROR; /* includes EOL */

  if (opt_tok) *opt_tok = bf->mem + bf->pos;
  if (opt_n)   *opt_n   = nc;
  bf->pos     += nskip;
  anch_set     = FALSE;
  if  ((status = esl_buffer_RaiseAnchor(bf, anch)) != eslOK) goto ERROR;
  return eslOK;
  
 ERROR:
  if (anch_set) esl_buffer_RaiseAnchor(bf, anch);
  if (opt_tok) *opt_tok = NULL;
  if (opt_n)   *opt_n   = 0;
  return status;
}

int
esl_buffer_FetchToken(ESL_BUFFER *bf, char *delim, char **ret_p, esl_pos_t *ret_n)
{

}

/* Function:  esl_buffer_FetchTokenAsStr()
 * Synopsis:  
 * Incept:    SRE, Sat Jan  1 19:31:57 2011 [Zaragoza]
 *
 * Purpose:   A 'token' consists of one or more characters that are
 *            neither in <*delim> nor a newline ('\r' or '\n').
 *
 *
 * Args:      
 *
 * Returns:   <eslOK> if a token is found; <*ret_p> points to it,
 *            and <*opt_n> is its length in chars (> 0). 
 *            
 *            <eslEOF> if the input ends before a token is found.
 *            Now <*ret_p> is an empty string "\0" and <*opt_n> is
 *            0.
 *            
 *            <eslEOL> if a line ends before a token is found.  (This
 *            case only arises if *delim does not contain newline.)
 *            Now <*ret_p> is an empty string "\0" and <*opt_n> is 0.
 *            
 *            <bf->mem[bf->pos]> is left at the next
 *            non-delim/non-newline character; or EOF if no such
 *            character remains in the input.
 *            
 *            <bf->mem> may be modified and/or reallocated, if new
 *            input reads are required to find the entire token.
 *
 * Throws:    <eslEMEM> if an allocation fails.
 *            Now <*ret_p> is <NULL> and <*opt_n> is 0.
 */
int
esl_buffer_FetchTokenAsStr(ESL_BUFFER *bf, char *delim, char **ret_p, esl_pos_t *opt_n)
{

}

/* First chunk of token parsing, shared amongst GetToken(), FetchToken(), FetchTokenAsStr()
 * Skip the parser past chars in *delim; return eslEOF if no non-delim char is found. 
 *   Now parser is bf->pos == bf->n, buffer is in EOF state.
 * If the parser stops at a newline, advance past it and return eslEOL.
 *   Now bf->pos is on first byte past the newline.
 * eslOK on success, and bf->pos is on the first nondelim char.
 *
 * Returns:  eslOK, eslEOF, or eslEOL
 * Throws:   eslEMEM, eslESYS, eslEINCONCEIVABLE
 */
static int
buffer_skipdelim(ESL_BUFFER *bf, char *delim)
{
  esl_pos_t nc = 0;

  /* skip characters in delim[], or hit EOF. */
  do {
    if ( (status = buffer_refill(bf, 0)) != eslOK) return status; /* eslEOF, or exceptions. */      
    for ( ; nc < bf->n-bf->pos; nc++)
      if (strchr(delim, bf->mem[bf->pos+nc]) == NULL) break;
  } while (nc == bf->n - bf->pos);
  bf->pos += nc;  /* bf->mem[bf->pos] now sitting at first char not in delim[] */
  
  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) return status;

  /* bf->pos now sitting on first char of a token */
  return eslOK;
}

/* buffer_skipnewline()
 * if bf->pos is on a newline (1 or 2 chars);
 *   advance bf->pos to next byte after newline and return eslEOL
 * else do nothing, and return eslOK.
 */
static int
buffer_newline(ESL_BUFFER *bf)
{
  esl_pos_t nc = bf->n - bf->pos;
  int       status;

  if (nc == 0) 
    return eslEOL;	/* no newline, but EOF is as good as */
  if (nc == 1 && bf->mem[bf->pos] == '\n')  
    { bf->pos += 1; return eslEOL; }
  if (nc >= 2 && memcmp(bf->mem + bf->pos, "\r\n", 2) == 0)
    { bf->pos += 2; return eslEOL; }
  
  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) return status;

  return eslOK;
}

/* bf->pos is sitting on a non-delim, non-newline character, starting
 * a token (i.e., the way buffer_skipdelim() left us on
 * success). Caller has set anchor to be sure this position stays in
 * buffer. Count how many nondelim, nonnewline characters there are,
 * starting here. Expand bf as needed.
 */
static int
buffer_counttok(ESL_BUFFER *bf, char *delim, esl_pos_t *ret_nc, esl_pos_t *ret_nskip)
{
  esl_pos_t nc = 1;

  /* skip chars NOT in delim[]. */
  while (1)
    {
      for ( ; nc < bf->n-bf->pos; nc++)
	if (strchr(delim, bf->mem[bf->pos+nc]) == NULL) continue;
      if (nc < bf->n-bf->pos) break;

      status = buffer_refill(bf, nc + bf->pagesize);
      if (status == eslEOF) break;	/* no more data - token extends to end of input */
      if (status != eslOK)  goto ERROR; /* exceptions */
    }
  /* bf->mem[bf->pos +nc] now sitting on the first char that's a delim */
  *ret_nskip = nc;

  /* check for whether the token ends in newline: back off if it does */
  if      (nc >= 2 && strncmp(bf->mem + bf->pos + nc - 2, "\r\n", 2) == 0) { nc -= 2; }
  else if (bf->mem[bf->pos+nc-1] == '\n') { nc--; }

  status = buffer_refill(bf, 0);
  if (status != eslEOF && status != eslOK) return status;

  *ret_nc = nc;
  return eslOK;
}
  


/*****************************************************************
 *# 6. Binary (fread-like) parsing
 *****************************************************************/

int
esl_buffer_GetBytes(ESL_BUFFER *bf, esl_pos_t nbytes, char **ret_p);
{

}

int 
esl_buffer_CopyBytes(ESL_BUFFER *bf, esl_pos_t nbytes, char *p, esl_pos_t nbytes)
{

}

int
esl_buffer_CopyPeekedBytes(ESL_BUFFER *bf, char *p, esl_pos_t nbytes)
{

}
/*------------ end, binary (fread-like) parsing -----------------*/

/*****************************************************************
 * Footnotes.
 *****************************************************************/
/*
 * [Open() returns int]
 *    Previously, Easel *Open() functions returned a pointer to the
 *    open object, or NULL on error; an exception to the normal Easel
 *    pattern of returning a status code. This is fine so long as the
 *    only error encountered is an eslEMEM, but for complicated
 *    Open()'s, this made it difficult to get back a sensible error
 *    status and error message. A better pattern is to return a status
 *    code and the open object, including the object's error message
 *    buffer. The caller can then use the errmsg and Close() the
 *    incomplete object.
 *    
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

off_t 
get_filesize(char *filename)
{
  struct stat fstats;
  stat(filename, &fstats);
  return fstats.st_size;
}
  
void
mmap_buffer(char *filename, char **ret_buf, off_t *ret_len)
{
  char *buf      = NULL;
  off_t len      = get_filesize(filename);
  int   prot     = PROT_READ;
  int   flags    = MAP_FILE | MAP_PRIVATE;
  int   fd       = -1;
  off_t offset   = 0;

  fd = open(filename, O_RDONLY);
  buf = (char *) mmap(0, len, prot, flags, fd, offset);
  //  posix_madvise(buf, len,  POSIX_MADV_SEQUENTIAL);
  close(fd);
  
  *ret_buf = buf;
  *ret_len = len;
  return;
}

void
read_buffer(char *filename, char **ret_buf, off_t *ret_len)
{
  off_t len      = get_filesize(filename);
  int   fd       = -1;
  char *buf      = malloc(len);

  fd = open(filename, O_RDONLY);
  read(fd, buf, len);
  
  *ret_buf = buf;
  *ret_len = len;
  return;
}


void
fread_buffer(char *filename, char **ret_buf, off_t *ret_len)
{
  FILE *fp       = NULL;
  off_t len      = get_filesize(filename);
  char *buf      = malloc(len);

  fp = fopen(filename, "r");
  fread(buf, 1, len, fp);
  
  *ret_buf = buf;
  *ret_len = len;
  return;
}


int
main(int argc, char **argv)
{
  char *filename = argv[1];
  char *buf;
  off_t len;
  esl_pos_t i;
  int   nseq     = 0;

  mmap_buffer(filename, &buf, &len);
  // fread_buffer(filename, &buf, &len);
  // read_buffer(filename, &buf, &len);

  for (i = 0; i < len; i++)
    if (buf[i] == '>') nseq++;
  printf("nseq = %d\n", nseq);

  return;
}


/*****************************************************************
 * Examples
 *****************************************************************/

int
example_read_fasta(ESL_BUFFER *bf, char **ret_name, char **ret_desc, char **ret_seq, char *ret_seqlen)
{
  char    *seqname = NULL;
  char    *seqdesc = NULL;
  char    *seq     = NULL;
  esl_pos_t  seqlen  = 0;
  esl_pos_t  salloc  = 0;
  char    *p;
  void    *tmp;
  esl_pos_t  n, i;
  int      status;

  ESL_MALLOC(seq, sizeof(char) * 256); salloc = 256;

  status = esl_buffer_GetBytes(bf, sizeof(char), &p);
  if (status == eslEOF) goto ERROR; /* normal EOF */
  if (*p != '>')        ESL_XFAIL(eslEINVAL, bf->errmsg, "Expected FASTA record to start with >");

  status = esl_buffer_GetToken(bf, " \t", &p, &n);
  if (status == eslEOF) ESL_XFAIL(eslEINVAL, bf->errmsg, "Premature eof while trying to parse sequence name");
  if (n == 0)           ESL_XFAIL(eslEINVAL, bf->errmsg, "Failed to parse a sequence name");
  esl_line_CloneAsStr(p, n, &seqname);

  status = esl_buffer_GetLine(bf, &p, &n);
  if (status == eslEOF) goto DONE; /* weird but ok. name, no desc, and a blank sequence. */
  esl_line_CloneAsStr(p, n, &seqdesc);

  while ((status = esl_buffer_GetLine(bf, &p, &n)) == eslOK)
    {
      if (p[0] == '>') { esl_buffer_Set(bf, p, 0); break; }

      if (seqlen+n+1 > salloc) { 
	ESL_RALLOC(seq, tmp, sizeof(char) * (seqlen+n+1));  
	salloc = seqlen+n+1; 
      }

      for (i = 0; i < n; i++) {
	if (! isspace(p[i])) { seq[seqlen++] = p[i]; }
      }
    }
    
 DONE:
  seq[seqlen] = '\0';
  *ret_name   = seqname;
  *ret_desc   = seqdesc;
  *ret_seq    = seq;
  *ret_seqlen = seqlen;
  return eslOK;

 ERROR:
  if (seqname) free(seqname);  *ret_name   = NULL;
  if (seqdesc) free(seqdesc);  *ret_desc   = NULL;
  if (seq)     free(seq);      *ret_seq    = NULL;
  *ret_seqlen = 0;
  return status;
}


int
example_read_lineblock(ESL_BUFFER *bf, char ***ret_lines, esl_pos_t **ret_lens, esl_pos_t nlines)
{
  char      **lines  = NULL;
  esl_pos_t  *lens   = NULL;
  esl_pos_t   nlines = 0;
  char       *p;
  esl_pos_t   n;
  esl_pos_t   start_offset;
  void       *tmp;
  int         status;

  /* skip blank lines */
  do {
    status = esl_buffer_GetPosition(bf, &start_offset);
    status = esl_buffer_GetLine(ESL_BUFFER *bf, &p, &n);
    if (status == eslEOF) return eslEOF; /* normal: no data */
  } while esl_line_IsSpace(p, n);

  /* anchor at current position */
  status = esl_buffer_SetAnchor(bf, start_offset);

  /* set pointers to non-blank lines */
  do {
    ESL_RALLOC(lines, tmp, sizeof(char *) * (nlines+1));
    ESL_RALLOC(lens,  tmp, sizeof(char *) * (nlines+1));
    
    lines[nlines] = p;
    lens[nlines]  = n;
    nlines++;
  } while ( (status = esl_buffer_GetLine(ESL_BUFFER *bf, &p, &n)) == eslOK &&
	    ! esl_line_IsSpace(p, n));
  
  /* restore parser state: clear anchor, and 
   * put it on first char of blank line (or EOF).
   * Set() should be called first, as a matter of idiom, to avoid
   * a reoffset/realloc that doesn't respect the anchor; 
   * though in this case the Set() call is backwards in the buffer,
   * so it will not reoffset/realloc.
   */
  status = esl_buffer_Set(bf, p, 0);
  status = esl_buffer_RaiseAnchor(bf, start_offset);

  *ret_lines  = lines;
  *ret_lens   = lens;
  *ret_nlines = nlines;
  return eslOK;

 ERROR:
  if (lines) free(lines);
  if (lens)  free(lens);
  *ret_lines  = NULL;
  *ret_lens   = NULL;
  *ret_nlines = 0;
  return status;
}
  


  

/*****************************************************************
 * @LICENSE@
 *****************************************************************/




















