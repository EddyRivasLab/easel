/* easel.c
 * SRE, Tue Oct 28 08:29:17 2003 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <easel/easel.h>

static esl_error_handler_f esl_error_handler = NULL;


void
esl_error_SetHandler(void (*handler)(int code, char *file, int line, 
				     char *format, va_list argp))
{
  esl_error_handler = handler;
}

void
esl_error_ResetDefaultHandler(void)
{
  esl_error_handler = NULL;
}

void
esl_error(int code, char *file, int line, char *format, ...)
{
  va_list argp;

  if (esl_error_handler != NULL) {
    va_start(argp, format);
    (*esl_error_handler)(code, file, line, format, argp);
    va_end(argp);
    return;
  } else {
    fprintf(stderr, "Easel fatal error (file %s, line %d):\n", file, line);
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
    abort();
  }
}



/* Function:  esl_Free2D()
 * Incept:    squid's Free2DArray(), 1999.
 *
 * Purpose:   Free a 2D pointer array <p>, where first dimension is
 *            <dim1>. (That is, the array is <p[0..dim1-1][]>.)
 *            Tolerates any of the pointers being NULL, to allow
 *            sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free2D(void **p, int dim1)
{
  int i;

  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) free(p[i]);
    free(p);
  }
  return;
}

/* Function:  esl_Free3D()
 * Incept:    squid's Free3DArray(), 1999.
 *
 * Purpose:   Free a 3D pointer array <p>, where first and second
 *            dimensions are <dim1>,<dim2>. (That is, the array is
 *            <p[0..dim1-1][0..dim2-1][]>.) Tolerates any of the
 *            pointers being NULL, to allow sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free3D(void ***p, int dim1, int dim2)
{
  int i, j;

  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) {
        for (j = 0; j < dim2; j++)
          if (p[i][j] != NULL) free(p[i][j]);
        free(p[i]);
      }
    free(p);
  }
}




/* Function: esl_strdup()
 * Date:     SRE, Wed May 19 17:57:28 1999 [St. Louis]
 *
 * Purpose: Returns a duplicate of string <s>. A version of the common
 *          but non-ANSI strdup() function. Can pass length <n>, if it's known,
 *          to save a strlen() call; else pass -1 to have the string length
 *          determined.
 *
 * Args:     s  - string to duplicate (NUL-terminated)
 *           n  - length of string, if known; -1 if unknown.
 *                
 * Returns:  allocated copy of string; NULL on failure.
 */
char *
esl_strdup(char *s, int n)
{
  char *new;

  if (s == NULL) return NULL;
  if (n < 0) n = strlen(s);
  if ((new = malloc(sizeof(char) * (n+1))) == NULL) return NULL;
  strcpy(new, s);
  return new;
}


/* Function: esl_fgets()
 * Date:     SRE, Thu May 13 10:56:28 1999 [St. Louis]
 *
 * Purpose:  Dynamic allocation version of fgets(),
 *           capable of reading unlimited line lengths.
 *
 * Args:     buf - ptr to a string (may be reallocated)
 *           n   - ptr to current allocated length of buf,
 *                 (may be changed)
 *           fp  - open file ptr for reading
 *           
 *           Before the first call to esl_parse_fgets(), 
 *           initialize buf to NULL and n to 0.
 *           They're a linked pair, so don't muck with the
 *           allocation of buf or the value of n while
 *           you're still doing esl_parse_fgets() calls with them.
 *
 * Returns:  ESL_OK   on success. 
 *           ESL_EOF  on normal end-of-file.
 *
 *           when ESL_OK:
 *           *buf points to a NUL-terminated line from the file.
 *           *n contains the current alloc'ed length for *buf.
 * 
 *           Caller must free *buf eventually. 
 *
 * Throws:   ESL_EMEM on malloc/realloc failure.
 *
 * Example:  char *buf;
 *           int   n;
 *           FILE *fp;
 *           
 *           fp  = fopen("my_file", "r");
 *           buf = NULL;
 *           n   = 0;
 *           while (esl_fgets(&buf, &n, fp) == ESL_OK) 
 *           {
 *             do stuff with buf;
 *           }
 */
int
esl_fgets(char **buf, int *n, FILE *fp)
{
  char *s;
  int   len;
  int   pos;

  if (*n == 0) 
    {
      if ((*buf = malloc(sizeof(char) * 128)) == NULL)
	ESL_ERROR(ESL_EMEM, "malloc failed");
      *n   = 128;
    }

  /* Simple case 1. We're sitting at EOF, or there's an error.
   *                fgets() returns NULL, so we return EOF.
   */
  if (fgets(*buf, *n, fp) == NULL) return ESL_EOF;

  /* Simple case 2. fgets() got a string, and it reached EOF doing it.
   *                return success status, so caller can use
   *                the last line; on the next call we'll
   *                return the 0 for the EOF.
   */
  if (feof(fp)) return ESL_OK;

  /* Simple case 3. We got a complete string, with \n,
   *                and don't need to extend the buffer.
   */
  len = strlen(*buf);
  if ((*buf)[len-1] == '\n') return ESL_OK;

  /* The case we're waiting for. We have an incomplete string,
   * and we have to extend the buffer one or more times. Make
   * sure we overwrite the previous fgets's \0 (hence +(n-1)
   * in first step, rather than 128, and reads of 129, not 128).
   */
  pos = (*n)-1;
  while (1) {
    *n  += 128;
    if ((*buf = realloc(*buf, sizeof(char) * (*n))) == NULL) 
      ESL_ERROR(ESL_EMEM, "realloc failed");
    s = *buf + pos;
    if (fgets(s, 129, fp) == NULL) return ESL_OK;
    len = strlen(s);
    if (s[len-1] == '\n') return ESL_OK;
    pos += 128;
  } 
  /*NOTREACHED*/
  return ESL_OK;
}


/* Function: esl_strtok()
 * Date:     SRE, Wed May 19 16:30:20 1999 [St. Louis]
 *
 * Purpose:  Thread-safe version of strtok().
 *
 *           Returns ptr to next token in a string: skips
 *            until it reaches a character that is not in the delim
 *            string, and sets beginning of token. Skips to
 *            next delim character (or '\0') to set the end; replaces that
 *            character with '\0'.
 *           If there's still more string left, sets s to point to next 
 *            character after the '\0' that was written, so successive 
 *            calls extract tokens in succession. If there was no string
 *            left, s points at the terminal '\0'. 
 *            
 *           If no token is found, returns NULL.
 *            
 *           Also returns the length of the token, which
 *           may save us a strlen() call in some applications.
 *           
 * Limitations:
 *           *s can't be a constant string, since we write to it.
 *                      
 * Example:  
 *           char *tok;
 *           int   len;
 *           char *s;             
 *           char  buf[50] = "This is  a sentence.";
 *           
 *           s = buf;  
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is "This"; s is "is  a sentence."; len is 4.
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is "is"; s is " a sentence."; len is 2.
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is "a"; s is "sentence."; len is 1.
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is "sentence."; s is "\0"; len is 9.
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is NULL; s is "\0", len is undefined.
 *       
 * Args:     s     - a tmp, modifiable ptr to string
 *           delim - characters that delimits tokens
 *           tok   - RETURN: ptr to \0-terminated token string
 *           len   - RETURN: length of token; pass NULL if not wanted
 *
 * Returns:  ESL_OK  on success: token points to next token, toklen is its len.
 *           ESL_EOL on end of line.
 *
 * Throws:   (no abnormal error conditions)
 */
int
esl_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen)
{
  char *begin, *end;
  int   n;

  if (ret_tok    != NULL) *ret_tok    = NULL;
  if (ret_toklen != NULL) *ret_toklen = 0;

  begin = *s;
  begin += strspn(begin, delim);
  if (! *begin) return ESL_EOL;

  n = strcspn(begin, delim);
  end  = begin + n;
  if (*end == '\0') { *s = end;}
  else {
    *end = '\0';
    *s   = end+1;
  }

  if (ret_tok    != NULL) *ret_tok    = begin;
  if (ret_toklen != NULL) *ret_toklen = n;
  return ESL_OK;
}





/******************************************************************************
 * File path/name manipulation functions                                      *
 *                                                                            *
 * Sufficiently widespread in the modules that we make them core.             *
 *****************************************************************************/

/* Function:  esl_FileExists()
 * Incept:    SRE, Sat Jan 22 09:07:24 2005 [St. Louis]
 *
 * Purpose:   Returns TRUE if <filename> exists, else FALSE.
 *     
 * Note:      Testing a read-only fopen() is the only portable ANSI C     
 *            I'm aware of. We could also use a POSIX func here, since
 *            we have a ESL_POSIX_AUGMENTATION flag in the code.
 *            
 * Xref:      squid's FileExists().
 */
int
esl_FileExists(char *filename)
{
  FILE *fp;
  if ((fp = fopen(filename, "r"))) { fclose(fp); return TRUE; }
  return FALSE;
}


/* Function:  esl_FileConcat()
 * Incept:    SRE, Sat Jan 22 07:28:46 2005 [St. Louis]
 *
 * Purpose:   Concatenates directory path prefix <dir> and a filename
 *            <file>, and returns the new full pathname through
 *            <ret_path>. If <dir> does not already end in the
 *            appropriate delimiter (e.g. / for UNIX), one is added.
 *            
 *            If <dir> is NULL, then <ret_path> is just the same as
 *            <file>. Similarly, if <file> already appears to be a
 *            full path (because its first character is a /), then
 *            <dir> is ignored and <ret_path> is the same as
 *            <file>. It wouldn't normally make sense for a caller to
 *            call this function with such arguments.
 *            
 *            <file> may be a relative path. For example, 
 *            if <dir> is "/usr/local" and <file> is "lib/myapp/data",
 *            <ret_path> will be "/usr/local/lib/myapp/data".
 *
 * Returns:   <ESL_OK> on success, and puts the path
 *            in <ret_path>; this string is allocated here, 
 *            and must be free'd by caller with <free()>.
 *
 * Throws:    <ESL_EMEM>   on allocation failure.
 *            <ESL_EINVAL> on bad argument.
 *
 * Xref:      squid's FileConcat().
 */
int
esl_FileConcat(char *dir, char *file, char **ret_path)
{
  char *path;
  int   nd, nf;

  if (ret_path != NULL) *ret_path = NULL;
  if (nf == NULL) ESL_ERROR(ESL_EINVAL, "null file");

  nd   = (dir  != NULL)? strlen(dir)  : 0;
  nf   = strlen(file);
  path = malloc (sizeof(char) * (nd+nf+2));
  if (path == NULL) ESL_ERROR(ESL_EMEM, "malloc failed");

  if (dir == NULL)		     /* 1. silly caller didn't give a path */
    strcpy(path, file);
  else if (*file == eslDIRSLASH)     /* 2. <file> is already a path?   */
    strcpy(path, file); 
  else if (dir[nd-1] == eslDIRSLASH) /* 3. <dir><file> (dir is / terminated) */
    sprintf(path, "%s%s", dir, file);
  else				     /* 4. <dir>/<file> (usual case)   */
    sprintf(path, "%s%c%s", dir, eslDIRSLASH, file);	

  if (ret_path != NULL) *ret_path = path; else free(path);
  return ESL_OK;
}


/* Function:  esl_FileNewSuffix()
 * Incept:    SRE, Sat Jan 22 10:04:08 2005 [St. Louis]
 *
 * Purpose:   Add a file suffix <sfx> to <filename>; or if <filename>
 *            already has a suffix, replace it with <sfx>. A suffix is
 *            usually 2-4 letters following a '.' character. Returns
 *            an allocated string containing the result in <ret_newpath>.
 *            
 *            For example, if <filename> is "foo" and <sfx> is "ssi",
 *            returns "foo.ssi". If <filename> is "foo.db" and <sfx>
 *            is "idx", returns "foo.idx". 
 *
 * Returns:   <ESL_OK> on success, and <ret_newpath> is set
 *            string "<base_filename>.<sfx>". Caller must <free()>
 *            this string.
 *
 * Throws:    <ESL_EMEM> on allocation failure.
 *
 * Xref:      squid's FileAddSuffix().
 */
int 
esl_FileNewSuffix(char *filename, char *sfx, char **ret_newpath)
{
  char *new;
  char *lastdot;
  int   nf;

  if (ret_newpath != NULL) *ret_newpath = NULL;

  lastdot   = strrchr(filename, '.'); /* check for suffix to replace */
  if (strchr(lastdot, eslDIRSLASH) != NULL) lastdot = NULL; /*foo.1/filename */
  nf = (lastdot == NULL)? strlen(filename) : lastdot-filename;
  
  new = malloc(sizeof(char) * (nf+strlen(sfx)+2)); /* '.' too */
  if (new == NULL) ESL_ERROR(ESL_EMEM, "malloc failed");
  strncpy(new, filename, nf);
  *(new+nf) = '.';
  strcpy(new+nf+1, sfx);

  if (ret_newpath != NULL) *ret_newpath = new; else free(new);
  return new;
}



/* Function:  esl_FileEnvOpen()
 * Incept:    SRE, Sat Jan 22 08:41:48 2005 [St. Louis]
 *
 * Purpose:   Looks for a file <fname> in a colon-separated list of
 *            directories that is configured in an environment variable
 *            <env>. The first occurrence of file <fname> in this directory 
 *            list is opened read-only. The open file ptr is returned
 *            through <ret_fp>, and the full path name to the file
 *            that was opened is returned through <ret_path>. 
 *            Caller can pass NULL in place of <ret_fp> or <ret_path>
 *            if it is not interested in one or both of these.
 *            
 *            Does not look in the current directory unless "." is
 *            explicitly in the directory list provided by <env>.
 *            
 * Note:      One reason to pass <ret_path> back to the caller is that
 *            sometimes we're opening the first in a group of files
 *            (for example, a database and its SSI index), and we want
 *            to make sure that after we find the main file, the
 *            caller can look for the auxiliary file(s) in exactly the
 *            same directory.
 *            
 * Examples:  % setenv BLASTDB /nfs/databases/blast-db:/nfs/databases/nr/
 *           
 *            FILE *fp;
 *            char *path;
 *            int   status;
 *            status = esl_FileEnvOpen("swiss42", "BLASTDB", &fp, &path);
 * 
 * Returns:   <ESL_OK> on success, and provides <ret_fp> and <ret_path>;
 *            <ret_fp> is opened here, and must be <fclose()>'d by caller;
 *            <ret_path> is allocated here, and must be <free()>'d by caller.
 *
 *            Returns <ESL_ENOTFOUND> if the file not found in any directory,
 *            or if <env> does not contain any directories to look in.
 *
 * Throws:    <ESL_EMEM> on allocation error.
 *
 * Xref:      squid's EnvFileOpen().
 */
int
esl_FileEnvOpen(char *fname, char *env, FILE **ret_fp, char **ret_path)
{
  FILE *fp;
  char *dirlist;		/* :-separated list of directories */
  char *s, *s2;                 /* ptrs into elems in env list */
  char *path;
  int   np;

  fp = NULL;
  if (ret_fp   != NULL) *ret_fp   = NULL;
  if (ret_path != NULL) *ret_path = NULL;

  if (env == NULL)               return ESL_ENOTFOUND;
  if ((s = getenv(env)) == NULL) return ESL_ENOTFOUND;
  if ((dirlist = esl_strdup(s, -1)) == NULL) /* need a modifiable copy */
    ESL_ERROR(ESL_EMEM, "malloc failed");

  np   = strlen(fname) + strlen(s) + 2; /* upper bound on full path len */
  path = malloc(sizeof(char) * np);
  if (path == NULL) { free(dirlist); ESL_ERROR(ESL_EMEM, "malloc failed");}

  s  = dirlist;
  while (s != NULL) 
    {
      if ((s2 = strchr(s, ':')) != NULL) { *s2 = '\0'; s2++;} /* ~=strtok(). */
      sprintf(path, "%s%c%s", s, eslDIRSLASH, fname); /* // won't hurt */
      if ((fp = fopen(path, "r")) != NULL) break;      
      s = s2;
    }
  if (fp == NULL) { free(path); free(dirlist); return ESL_ENOTFOUND; }

  if (ret_dir != NULL) { *ret_dir = path; } else free(path);
  if (ret_fp  != NULL) { *ret_fp  = fp; }   else fclose(fp);
  free(dirlist);
  return ESL_OK;
}

/*----------------- end of file path/name functions ------------------------*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/  
