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
esl_error_SetHandler(void (*handler)(int code, char *file, int line, char *format, va_list argp))
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

void *
esl_malloc(char *file, int line, size_t size)
{
  void *ptr;

  if ((ptr = malloc (size)) == NULL) {
    esl_error(ESL_EMEM, file, line, "malloc of %ld bytes failed\n");
    return NULL;
  }
  return ptr;
}

void *
esl_realloc(char *file, int line, void *p, size_t size)
{
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL) {
    esl_error(ESL_EMEM, file, line, "realloc of %ld bytes failed\n", size);
    return NULL;
  }
  return ptr;
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
 * Returns:  ESL_OK  on success: token points to next token, toklen is its length.
 *           ESL_EOL on end of line.
 *
 * Throws:   (no abnormal error conditions)
 */
int
esl_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen)
{
  char *begin, *end;
  int   n;

  begin = *s;
  begin += strspn(begin, delim);
  if (! *begin) {*ret_tok = NULL; return ESL_EOL;}

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

/*****************************************************************
 * @LICENSE@
 *****************************************************************/  
