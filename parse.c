/* esl_fileparser.c
 * A simple token-based file parsing system.
 * 
 * SRE, Tue Jul 13 14:41:52 2004 [St. Louis]
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <easel.h>
#include <esl_parse.h>

static int nextline(ESL_FILEPARSER *efp);

/* Function:  esl_fileparser_Create()
 * Incept:    SRE, Fri Jul  9 12:50:29 2004 [St. Louis]
 *
 * Purpose:   Take an open file <fp>, and transform it to
 *            a fileparser object -- preparing to parse it
 *            one whitespace-delimited field at a time.
 *
 * Args:      fp  - open FILE to parse
 *
 * Returns:   a new <ESL_FILEPARSER> object, which must be 
 *            free'd by the caller with <esl_fileparser_Destroy()>.
 *
 * Throws:    <eslEMEM>:  malloc/realloc failed.
 *            
 * Xref:      STL8 p.56.
 */
ESL_FILEPARSER *
esl_fileparser_Create(FILE *fp)
{
  ESL_FILEPARSER *efp;

  if ((efp = malloc(sizeof(ESL_FILEPARSER))) == NULL)
    ESL_ERROR_NULL(eslEMEM, "malloc failed");
  efp->fp          = fp;
  efp->buf         = NULL;
  efp->buflen      = 0;
  efp->s           = NULL;
  efp->commentchar = '\0';

  *ret_efp = efp;
  return efp;
}

/* Function:  esl_fileparser_SetCommentChar()
 * Incept:    SRE, Sat Jul 10 10:18:35 2004 [St. Louis]
 *
 * Purpose:   Defines a single character <c> for comments. Anything
 *            on a line following this character is ignored
 *            when parsing.
 *            
 *            <c> = '#' is a common convention.
 *
 * Args:      efp - open fileparser
 *            c    - comment character ('#')        
 *
 * Returns:   eslOK.
 */
int
esl_fileparser_SetCommentChar(ESL_FILEPARSER *efp, char c)
{
  efp->commentchar = c;
  return eslOK;
}


/* Function:  esl_fileparser_GetToken()
 * Incept:    SRE, Fri Jul  9 13:03:50 2004 [St. Louis]
 *
 * Purpose:   Sets a pointer to the next field in the 
 *            file we're parsing.
 *
 * Args:      efp        - open fileparser
 *            ret_tok    - RETURN: ptr to next field
 *            ret_toklen - RETURN: length of tok.       
 *
 * Returns:   <eslOK> if <tok>, <toklen> contain valid data.
 *            <eslEOF> on normal end-of-file.
 *            
 * Throws:    <eslEMEM> malloc/realloc failure somewhere.
 *
 * Xref:      STL8 p.56.
 */
int
esl_fileparser_GetToken(ESL_FILEPARSER *efp, char **ret_tok, int *ret_toklen)
{
  char *tok    = NULL;
  int   toklen = 0;
  int   tokcode, fcode;
  int   goodtok;

  if (ret_tok != NULL)    *ret_tok    = NULL;
  if (ret_toklen != NULL) *ret_toklen = 0;

  if (efp->buf == NULL) {
    fcode = nextline(efp);
    if (fcode != eslOK) return fcode;
  }

  do {
    goodtok = FALSE;
    tokcode = esl_strtok(&(efp->s), " \t\n", &tok, &toklen);
    if (tokcode == eslEOL ||
	(tokcode == eslOK && *tok == efp->commentchar)) 
      {
	fcode = nextline(efp);
	if (fcode != eslOK) return fcode;
      } 
    else if (tokcode == eslOK) goodtok = TRUE;
    else return tokcode;
  } while (! goodtok);

  if (ret_tok != NULL)    *ret_tok    = tok;
  if (ret_toklen != NULL) *ret_toklen = toklen;
  return eslOK;
}

/* Function:  esl_fileparser_Destroy()
 * Incept:    SRE, Fri Jul  9 13:22:36 2004 [St. Louis]
 *
 * Purpose:   Close an open <ESL_FILEPARSER>. The original fp is
 *            still open - whoever opened it is still
 *            responsible for closing it.
 *
 * Args:      efp - fileparser to shut down.
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      STL8 p.56.
 */
int
esl_fileparser_Destroy(ESL_FILEPARSER *efp)
{
  if (efp->buf != NULL) free(efp->buf);
  free(efp);
  return eslOK;
}




/* nextline()
 *
 * Purpose:   Skip the file parser to the next line (for instance,
 *            if an end-of-line comment is found).
 *
 * Args:      efp  - open file parser
 *
 * Returns:   eslOK:   success
 *            eslEOF:  normal end of file
 *
 * Throws:    eslEMEM: malloc/realloc failed in fgets()
 *
 * Xref:      STL8 p.56
 */
static int
nextline(ESL_FILEPARSER *efp)
{
  int eslcode;

  if ((eslcode = esl_fgets(&(efp->buf), &(efp->buflen), efp->fp)) != eslOK) return eslcode;
  efp->s = efp->buf;
  return eslOK;
}
