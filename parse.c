/* parse.c
 * 
 * SRE, Tue Jul 13 14:41:52 2004 [St. Louis]
 * SVN $Id$
 */

#include <string.h>
#include <stdlib.h>
#include <easel/easel.h>
#include <easel/parse.h>

/* Function: esl_parse_fgets()
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
 *           while (esl_parse_fgets(&buf, &n, fp) != ESL_OK) 
 *           {
 *             do stuff with buf;
 *           }
 */
int
esl_parse_fgets(char **buf, int *n, FILE *fp)
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


/* Function: esl_parse_strtok()
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
 *           esl_parse_strtok(&s, " ", &tok, &len);
 *                tok is "This"; s is "is  a sentence."; len is 4.
 *           esl_parse_strtok(&s, " ", &tok, &len);
 *                tok is "is"; s is " a sentence."; len is 2.
 *           esl_parse_strtok(&s, " ", &tok, &len);
 *                tok is "a"; s is "sentence."; len is 1.
 *           esl_parse_strtok(&s, " ", &tok, &len);
 *                tok is "sentence."; s is "\0"; len is 9.
 *           esl_parse_strtok(&s, " ", &tok, &len);
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
esl_parse_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen)
{
  char *begin, *end;
  int   n;

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


/* Function:  esl_fileparse_create()
 * Incept:    SRE, Fri Jul  9 12:50:29 2004 [St. Louis]
 *
 * Purpose:   Take an open file (fp), and transform it to
 *            a fileparser object -- preparing to parse it
 *            one field at a time.
 *
 * Args:      fp        - open FILE to parse
 *            ret_efp   - RETURN: open ESL_FILEPARSER object
 *
 * Returns:   ESL_OK on success. ret_efp points to a new
 *                   ESL_FILEPARSER object, which must be 
 *                   destroyed with esl_fileparse_close().
 *
 * Throws:    ESL_EMEM:  malloc/realloc failed.
 *            
 * Xref:      STL8 p.56.
 */
int
esl_fileparse_create(FILE *fp, ESL_FILEPARSER **ret_efp)
{
  ESL_FILEPARSER *efp;

  *ret_efp = NULL;
  if ((efp = malloc(sizeof(ESL_FILEPARSER))) == NULL)
    ESL_ERROR(ESL_EMEM, "malloc failed");
  efp->fp          = fp;
  efp->buf         = NULL;
  efp->buflen      = 0;
  efp->s           = NULL;
  efp->commentchar = '\0';

  *ret_efp = efp;
  return ESL_OK;
}

/* Function:  esl_fileparse_set_commentchar()
 * Incept:    SRE, Sat Jul 10 10:18:35 2004 [St. Louis]
 *
 * Purpose:   Defines a single character for comments; anything
 *            on a line following this character is ignored
 *            when parsing.
 *            
 *            '#' is a common convention.
 *
 * Args:      efp  - open fileparser
 *            c    - comment character ('#')        
 *
 * Returns:   ESL_OK.
 */
int
esl_fileparse_set_commentchar(ESL_FILEPARSER *efp, char c)
{
  efp->commentchar = c;
  return ESL_OK;
}


/* Function:  esl_fileparse_nextline()
 * Incept:    SRE, Fri Jul  9 12:58:32 2004 [St. Louis]
 *
 * Purpose:   Skip the file parser to the next line (for instance,
 *            if an end-of-line comment is found).
 *
 * Args:      efp  - open file parser
 *
 * Returns:   ESL_OK:   success
 *            ESL_EOF:  normal end of file
 *
 * Throws:    ESL_EMEM: malloc/realloc failed in fgets()
 *
 * Xref:      STL8 p.56
 */
int
esl_fileparse_nextline(ESL_FILEPARSER *efp)
{
  int eslcode;

  if ((eslcode = esl_parse_fgets(&(efp->buf), &(efp->buflen), efp->fp)) != ESL_OK) return eslcode;
  efp->s = efp->buf;
  return ESL_OK;
}

/* Function:  esl_fileparse_token()
 * Incept:    SRE, Fri Jul  9 13:03:50 2004 [St. Louis]
 *
 * Purpose:   Returns a pointer to the next field in the 
 *            file we're parsing.
 *
 * Args:      efp        - open fileparser
 *            ret_tok    - RETURN: ptr to next field
 *            ret_toklen - RETURN: length of tok.       
 *
 * Returns:   ESL_OK:  tok, toklen contain valid data.
 *            ESL_EOF: normal end-of-file.
 *            
 * Throws:    ESL_EMEM: malloc/realloc failure somewhere.
 *
 * Xref:      STL8 p.56.
 */
int
esl_fileparse_token(ESL_FILEPARSER *efp, char **ret_tok, int *ret_toklen)
{
  char *tok    = NULL;
  int   toklen = 0;
  int   tokcode, fcode;
  int   goodtok;

  if (ret_tok != NULL)    *ret_tok    = NULL;
  if (ret_toklen != NULL) *ret_toklen = 0;

  if (efp->buf == NULL) {
    fcode = esl_fileparse_nextline(efp);
    if (fcode != ESL_OK) return fcode;
  }

  do {
    goodtok = FALSE;
    tokcode = esl_parse_strtok(&(efp->s), " \t\n", &tok, &toklen);
    if (tokcode == ESL_EOL ||
	(tokcode == ESL_OK && *tok == efp->commentchar)) 
      {
	fcode = esl_fileparse_nextline(efp);
	if (fcode != ESL_OK) return fcode;
      } 
    else if (tokcode == ESL_OK) goodtok = TRUE;
    else return tokcode;
  } while (! goodtok);

  if (ret_tok != NULL)    *ret_tok    = tok;
  if (ret_toklen != NULL) *ret_toklen = toklen;
  return ESL_OK;
}

/* Function:  esl_fileparse_free()
 * Incept:    SRE, Fri Jul  9 13:22:36 2004 [St. Louis]
 *
 * Purpose:   Close an open ESL_FILEPARSER. The original fp is
 *            still open - whoever opened it is still
 *            responsible for closing it.
 *
 * Args:      efp - fileparser to shut down.
 *
 * Returns:   ESL_OK, success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      STL8 p.56.
 */
int
esl_fileparse_free(ESL_FILEPARSER *efp)
{
  if (efp->buf != NULL) free(efp->buf);
  free(efp);
  return ESL_OK;
}
