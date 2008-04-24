/* A simple token-based file parsing system.
 * 
 * Contents:
 *    1. The ESL_FILEPARSER object and its API.
 *    2. Private functions.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Examples.
 *    6. Copyright/license information.
 * 
 * SRE, Tue Jul 13 14:41:52 2004 [St. Louis]
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_fileparser.h"

static int nextline(ESL_FILEPARSER *efp);

/*****************************************************************
 * 1. The ESL_FILEPARSER object and its API.
 *****************************************************************/

/* Function:  esl_fileparser_Open()
 * Incept:    SRE, Tue Apr  3 08:09:56 2007 [Janelia]
 *
 * Purpose:   Opens <filename> for reading. 
 * 
 *            As a special case, if <filename> is "-", set up the
 *            fileparser to read and parse <stdin>.
 *
 * Returns:   <eslOK> on success, and <ret_fp> points
 *            to a new <ESL_FILEPARSER> object.
 *            
 *            Returns <eslENOTFOUND> if <filename> can't
 *            be opened for reading, and <ret_fp> is set
 *            to <NULL>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_fileparser_Open(char *filename, ESL_FILEPARSER **ret_efp)
{
  int             status;
  ESL_FILEPARSER *efp = NULL;

  if ((efp = esl_fileparser_Create(NULL)) == NULL) { status = eslEMEM;      goto ERROR; }

  if (strcmp(filename, "-") == 0)
    efp->fp = stdin;
  else {
    if ((efp->fp = fopen(filename, "r")) == NULL)    { status = eslENOTFOUND; goto ERROR; }
  }
  *ret_efp = efp;
  return eslOK;

 ERROR:
  esl_fileparser_Close(efp);
  *ret_efp = NULL;
  return status;
}


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
 * Throws:    <eslEMEM> if an allocation failed.
 *            
 * Xref:      STL8 p.56.
 */
ESL_FILEPARSER *
esl_fileparser_Create(FILE *fp)
{
  int status;
  ESL_FILEPARSER *efp = NULL;

  ESL_ALLOC(efp, sizeof(ESL_FILEPARSER));
  efp->fp          = fp;
  efp->buf         = NULL;
  efp->buflen      = 0;
  efp->s           = NULL;
  efp->commentchar = '\0';
  efp->tok         = NULL;
  efp->toklen      = 0;
  efp->linenumber  = 0;
  efp->errbuf[0]   = '\0';

  return efp;
  
 ERROR:
  esl_fileparser_Destroy(efp);
  return NULL;
}




/* Function:  esl_fileparser_SetCommentChar()
 * Incept:    SRE, Sat Jul 10 10:18:35 2004 [St. Louis]
 *
 * Purpose:   Defines a single character <c> for comments. Anything
 *            on a line following this character is ignored
 *            when parsing.
 *
 * Args:      efp - open fileparser
 *            c    - comment character ('#', for example)        
 *
 * Returns:   <eslOK> on success.
 */
int
esl_fileparser_SetCommentChar(ESL_FILEPARSER *efp, char c)
{
  efp->commentchar = c;
  return eslOK;
}


/* Function:  esl_fileparser_NextLine()
 * Incept:    SRE, Tue Apr  3 08:27:22 2007 [Janelia]
 *
 * Purpose:   Advance the parser to the next non-blank, non-comment
 *            data line that contains at least one token. 
 *
 * Returns:   <eslOK> on success.
 *            <eslEOF> if no more tokens remain in the file.  
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_fileparser_NextLine(ESL_FILEPARSER *efp)
{
  int   status;
  char *tok    = NULL;
  int   toklen = 0;
  int   tokcode;

  while ((status = nextline(efp)) == eslOK) 
    {
      tokcode = esl_strtok(&(efp->s), " \t\r\n", &tok, &toklen);
      if (tokcode == eslEOL || (tokcode == eslOK && *tok == efp->commentchar)) continue; /* no tokens on line */
      if (tokcode != eslOK) ESL_XFAIL(tokcode, efp->errbuf, "esl_strtok() failed");
      break;
    } 
  if (status == eslEOF) return eslEOF;
  if (status != eslOK)  ESL_XFAIL(status, efp->errbuf, "nextline() failed");
  
  /* Remember this token. The next GetToken call will regurgitate it instead
   * of finding its own. efp->s now points to the character after tok ends. 
   */
  efp->tok    = tok;
  efp->toklen = toklen;
  return eslOK;

 ERROR:
  efp->tok    = NULL;
  efp->toklen = 0;
  return status;
}  


/* Function:  esl_fileparser_GetToken()
 * Incept:    SRE, Fri Jul  9 13:03:50 2004 [St. Louis]
 *
 * Purpose:   Sets a pointer to the next field in the 
 *            file we're parsing.
 *            
 *            The <ret_tok> pointer is into an internal line buffer
 *            that may be invalidated upon the next call to a
 *            <fileparser> function. If you want to store it, make a
 *            copy.
 *
 * Args:      efp        - open fileparser
 *            ret_tok    - RETURN: ptr to next field
 *            ret_toklen - RETURN: length of tok.       
 *
 * Returns:   <eslOK> if <tok>, <toklen> contain valid data.
 *            <eslEOF> on normal end-of-file.
 *            
 * Throws:    <eslEMEM> if an allocation fails.
 *
 * Xref:      STL8 p.56.
 */
int
esl_fileparser_GetToken(ESL_FILEPARSER *efp, char **ret_tok, int *ret_toklen)
{
  char *tok    = NULL;
  int   toklen = 0;
  int   tokcode;
  int   fcode;
  int   goodtok;

  if (ret_tok != NULL)    *ret_tok    = NULL;
  if (ret_toklen != NULL) *ret_toklen = 0;

  /* Do we already have a token from a NextLine() call? */
  if (efp->tok != NULL) {
    if (ret_tok    != NULL) *ret_tok    = efp->tok;
    if (ret_toklen != NULL) *ret_toklen = efp->toklen;
    efp->tok    = NULL;
    efp->toklen = 0;
    return eslOK;
  }

  /* If not, then find next token.
   */

  /* First, make sure we have a line loaded. 
   * On the first call to GetToken, we won't.
   */
  if (efp->buf == NULL) {
    fcode = nextline(efp);
    if (fcode != eslOK) return fcode;
  }

  /* Start strtok()'ing this line to try to find token.
   * If we don't find one, keep loading lines until we
   * do, or we run out of data.
   */
  do {
    goodtok = FALSE;
    tokcode = esl_strtok(&(efp->s), " \t\r\n", &tok, &toklen);
    if (tokcode == eslEOL ||
	(tokcode == eslOK && *tok == efp->commentchar)) 
      {
	fcode = nextline(efp);
	if (fcode != eslOK) return fcode;
      } 
    else if (tokcode == eslOK) goodtok = TRUE;
    else 
      { sprintf(efp->errbuf, "esl_strtok() failed"); return tokcode;}
  } while (! goodtok);

  if (ret_tok != NULL)    *ret_tok    = tok;
  if (ret_toklen != NULL) *ret_toklen = toklen;
  return eslOK;
}

/* Function:  esl_fileparser_GetTokenOnLine()
 * Incept:    SRE, Tue Apr  3 08:46:59 2007 [Janelia]
 *
 * Purpose:   Same as <esl_fileparser_GetToken()>, except that it only
 *            retrieves tokens from the line that the parser is
 *            on. When it runs out of tokens on the line, it returns
 *            <eslEOL>. This allows a caller to count the tokens on a
 *            line (whereas <GetToken()> reads through newlines
 *            silently).
 *            
 *            The <ret_tok> pointer is into an internal line buffer
 *            that may be invalidated upon the next call to a
 *            <fileparser> function. If you want to store it, make a
 *            copy.
 *            
 *            Normally, a call to <esl_fileparser_GetTokenOnLine()>
 *            would be preceded by <esl_fileparser_NextLine()> to
 *            position the parser on the next data line with at least
 *            one token on it. However, you could also conceivably
 *            call <esl_fileparser_GetTokenOnLine()> after one or more
 *            calls to <esl_fileparser_GetToken()>, to get remaining
 *            tokens from a given line. What you can't do is to call
 *            <esl_fileparser_GetTokenOnLine()> immediately after 
 *            opening a file; the parser won't have a line loaded yet.
 *            (In this case, it would return <eslEOL>.)
 *
 * Returns:   <eslOK> on success, and the token and its length are
 *            in <ret_tok> and <ret_toklen>.
 *            
 *            Returns <eslEOL> if no more tokens exist on the line;
 *            in this case <ret_tok> is set to <NULL> and <ret_toklen>
 *            to 0.
 */
int
esl_fileparser_GetTokenOnLine(ESL_FILEPARSER *efp, char **ret_tok, int *ret_toklen)
{
  int status;
  char *tok    = NULL;
  int   toklen = 0;

  /* Do we already have a token from a NextLine() call? */
  if (efp->tok != NULL) {
    if (ret_tok    != NULL) *ret_tok    = efp->tok;
    if (ret_toklen != NULL) *ret_toklen = efp->toklen;
    efp->tok    = NULL;
    efp->toklen = 0;
    return eslOK;
  }

  /* No line loaded? Then we can't find any token on it. */
  if (efp->buf == NULL) { status = eslEOL;  goto ERROR; }

  /* Find next token in the line loaded in the parser. */
  status = esl_strtok(&(efp->s), " \t\r\n", &tok, &toklen);
  if (status == eslEOL) goto ERROR;
  if (status != eslOK)  goto ERROR;
  if (status == eslOK && *tok == efp->commentchar) { status = eslEOL; goto ERROR; }

  if (ret_tok    != NULL) *ret_tok    = tok;
  if (ret_toklen != NULL) *ret_toklen = toklen;
  return eslOK;

 ERROR:
  if (ret_tok    != NULL) *ret_tok    = NULL;
  if (ret_toklen != NULL) *ret_toklen = 0;
  return status;
}


/* Function:  esl_fileparser_Destroy()
 * Incept:    SRE, Fri Jul  9 13:22:36 2004 [St. Louis]
 *
 * Purpose:   Frees an open <ESL_FILEPARSER>. The original fp is
 *            still open - whoever opened it is still
 *            responsible for closing it.
 *
 * Xref:      STL8 p.56.
 */
void
esl_fileparser_Destroy(ESL_FILEPARSER *efp)
{
  if (efp->buf != NULL) free(efp->buf);
  free(efp);
}

/* Function:  esl_fileparser_Close()
 * Incept:    SRE, Tue Apr  3 08:18:11 2007 [Janelia]
 *
 * Purpose:   Closes an open <ESL_FILEPARSER>, including the 
 *            file it opened. 
 */
void
esl_fileparser_Close(ESL_FILEPARSER *efp)
{
  if (efp == NULL) return;
  
  if (efp->fp != NULL && efp->fp != stdin) fclose(efp->fp);
  esl_fileparser_Destroy(efp);
}



/*****************************************************************
 * 2. Private functions
 *****************************************************************/

/* nextline()
 *
 * Purpose:   Skip the file parser to the next line (for instance,
 *            if an end-of-line comment is found). The new line might
 *            have no tokens on it.
 *
 * Args:      efp  - open file parser
 *
 * Returns:   eslOK:   success
 *            eslEOF:  normal end of file
 *
 * Throws:    <eslEMEM> if a reallocation failed in fgets()
 *
 * Xref:      STL8 p.56
 */
static int
nextline(ESL_FILEPARSER *efp)
{
  int status;

  if ((status = esl_fgets(&(efp->buf), &(efp->buflen), efp->fp)) != eslOK) 
    { sprintf(efp->errbuf, "esl_fgets() failed"); return status;}
  efp->s = efp->buf;
  efp->linenumber++;
  return eslOK;
}



/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef eslFILEPARSER_TESTDRIVE
/* test the interface for getting all tokens in a file, regardless
 * of newlines. Also, uses the Create/Destroy interface instead of
 * Open/Close.
 */
static void
utest_GetToken(char *filename)
{
  int status;
  ESL_FILEPARSER *efp = NULL;
  FILE           *fp  = NULL;
  char           *tok = NULL;
  int             toklen = 0;
  int             ntok   = 0;

  if ((fp  = fopen(filename, "r"))      == NULL)  esl_fatal("File open failed");
  if ((efp = esl_fileparser_Create(fp)) == NULL)  esl_fatal("Failed to associate stream with fileparser");
  esl_fileparser_SetCommentChar(efp, '#');
  
  while ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK)
    {
      if (toklen != 6)                   esl_fatal("bad token %s", tok);
      if (strncmp(tok, "token", 5) != 0) esl_fatal("bad token %s", tok);
      ntok++;
    }
  if (status != eslEOF)  esl_fatal("Abnormal parse termination");
  if (ntok != 5)         esl_fatal("bad total token number %d\n", ntok);
  
  esl_fileparser_Destroy(efp);
  fclose(fp);
  return;
}

/* test the NextLine and GetTokenOnLine interface, as well as the
 * Open/Close interface.
 */
static void
utest_GetTokenOnLine(char *filename)
{
  int status;
  ESL_FILEPARSER *efp = NULL;
  char           *tok = NULL;
  int             toklen = 0;
  int             ntok   = 0;
  int             nlines = 0;
  char            expect[32];

  if (esl_fileparser_Open(filename, &efp) != eslOK) esl_fatal("File open failed");
  esl_fileparser_SetCommentChar(efp, '#');

  while ((status = esl_fileparser_NextLine(efp)) == eslOK)
    {
      nlines++;
      while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) == eslOK)
	{
	  ntok++;
	  sprintf(expect, "token%d", ntok);
	  if (toklen != 6)               esl_fatal("bad token length for %s", tok);
	  if (strcmp(expect, tok) != 0)  esl_fatal("bad token %s", tok);
	}
      if (status != eslEOL) esl_fatal("Unexpected code in place of end-of-line");
    }
  if (status != eslEOF) esl_fatal("Unexpected code in place of end-of-file.");

  if (nlines != 3) esl_fatal("expected to parse 3 lines; parsed %d", nlines);
  if (ntok   != 5) esl_fatal("expected to parse 5 tokens; parsed %d", ntok);
  
  esl_fileparser_Close(efp);
  return;
}
#endif /*eslFILEPARSER_TESTDRIVE*/

/*****************************************************************
 * 4. Test driver.
 *****************************************************************/

/*
    gcc -g -Wall -I. -o test -DeslFILEPARSER_TESTDRIVE esl_fileparser.c easel.c
    ./test
*/
#ifdef eslFILEPARSER_TESTDRIVE
#include <stdio.h>
#include <string.h>
#include "easel.h"
#include "esl_fileparser.h"

int 
main(int argc, char **argv)
{
  char  tmpfile[32] = "esltmpXXXXXX";
  FILE *fp;

  /* Create a test file to read.
   */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("File open failed");
  fprintf(fp, "# Full line comment\n");
  fprintf(fp, "token1  # Trailing comment\n");
  fprintf(fp, "\n");		/* blank line */
  fprintf(fp, "   \n");		/* whitespace line */
  fprintf(fp, "   # sowing comment/whitespace confusion...\n"); 
  fprintf(fp, "token2\ttoken3  token4\n");
  fprintf(fp, "token5");	/* file ends w/ no \n */
  fclose(fp);

  /* Run unit tests using that file.
   * Unit tests have hardwired knowledge of what's supposed to be in the file.
   */
  utest_GetToken(tmpfile);
  utest_GetTokenOnLine(tmpfile);

  remove(tmpfile);
  return 0;
}
#endif /*eslFILEPARSER_TESTDRIVE*/


/*****************************************************************
 * 5. Examples.
 *****************************************************************/

/* The first example shows the simplest interface: get all tokens
 * in the file, one at a time.
 *
     gcc -g -Wall -I. -o example -DeslFILEPARSER_EXAMPLE esl_fileparser.c easel.c
     ./example <any file>
 */
#ifdef eslFILEPARSER_EXAMPLE
/*::cexcerpt::fileparser_example::begin::*/
#include <stdio.h>
#include "easel.h"
#include "esl_fileparser.h"

int 
main(int argc, char **argv)
{
  char           *filename = argv[1];
  int             ntok     = 1;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;

  if (esl_fileparser_Open(filename, &efp) != eslOK) esl_fatal("File open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_GetToken(efp, &tok, &toklen) == eslOK) { 
    printf("%5d %3d %s\n", ntok, toklen, tok); 
    ntok++;  
  }
  esl_fileparser_Close(efp);
  return 0;
}
/*::cexcerpt::fileparser_example::end::*/
#endif /*eslFILEPARSER_EXAMPLE*/

/* The second example shows the more line-oriented interface
 * of NextLine(), GetTokenOnLine().
     gcc -g -Wall -I. -o example -DeslFILEPARSER_EXAMPLE2 esl_fileparser.c easel.c
     ./example <any file>
 */
#ifdef eslFILEPARSER_EXAMPLE2
/*::cexcerpt::fileparser_example2::begin::*/
#include <stdio.h>
#include "easel.h"
#include "esl_fileparser.h"

int 
main(int argc, char **argv)
{
  char           *filename = argv[1];
  int             nline    = 1;
  int             ntok;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;

  if (esl_fileparser_Open(filename, &efp) != eslOK) esl_fatal("File open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK)
  {
    ntok = 0;
    while (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen) == eslOK)
      ntok++;
    printf("Line %d in the file (%d non-blank, non-comment) contains %d tokens...\n", 
	   efp->linenumber, nline, ntok);
    nline++;
  }
  esl_fileparser_Close(efp);
  return 0;
}
/*::cexcerpt::fileparser_example2::end::*/
#endif /*eslFILEPARSER_EXAMPLE*/






/*****************************************************************
 * @LICENSE@
 *****************************************************************/

