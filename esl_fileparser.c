/* esl_fileparser.c
 * A simple token-based file parsing system.
 * 
 * SRE, Tue Jul 13 14:41:52 2004 [St. Louis]
 * SVN $Id$
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <easel.h>
#include <esl_fileparser.h>

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
  efp->linenumber  = 0;
  efp->errbuf[0]   = '\0';

  return efp;
  
 FAILURE:
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

  if (efp->buf == NULL) {
    fcode = nextline(efp);
    if (fcode != eslOK) return fcode;
  }

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
 * Throws:    <eslEMEM> if a reallocation failed in fgets()
 *
 * Xref:      STL8 p.56
 */
static int
nextline(ESL_FILEPARSER *efp)
{
  int eslcode;

  if ((eslcode = esl_fgets(&(efp->buf), &(efp->buflen), efp->fp)) != eslOK) 
    { sprintf(efp->errbuf, "esl_fgets() failed"); return eslcode;}
  efp->s = efp->buf;
  efp->linenumber++;
  return eslOK;
}



/*****************************************************************
 * Example:
 *    gcc -g -Wall -I. -o example -DeslFILEPARSER_EXAMPLE esl_fileparser.c easel.c
 *    ./example <any file>
 * Reads whitespace-delimited tokens from a file, and prints them
 * out one at a time.
 *****************************************************************/
#ifdef eslFILEPARSER_EXAMPLE
/*::cexcerpt::fileparser_example::begin::*/
#include <stdio.h>
#include <easel.h>
#include <esl_fileparser.h>

int 
main(int argc, char **argv)
{
  ESL_FILEPARSER *efp;
  char *filename;
  FILE *fp;
  char *tok;
  int   toklen;
  int   status;
  int   ntok;

  filename = argv[1];           
  if ((fp = fopen(filename, "r")) == NULL) 
    esl_fatal("File open failed");
  
  if ((efp = esl_fileparser_Create(fp)) == NULL) 
    esl_fatal("Failed to associate stream with fileparser");
  esl_fileparser_SetCommentChar(efp, '#');
  
  ntok = 1;
  while ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK)
    {
      printf("%5d %3d %s\n", ntok, toklen, tok);
      ntok++;
    }
  if (status != eslEOF)
    esl_fatal("Abnormal parse termination at line %d of file %s;\n%s",
	      efp->linenumber, filename, efp->errbuf);
  
  esl_fileparser_Destroy(efp);
  fclose(fp);
  return 0;
}
/*::cexcerpt::fileparser_example::end::*/
#endif /*eslFILEPARSER_EXAMPLE*/




/*****************************************************************
 * Test driver:
 *    gcc -g -Wall -I. -o test -DeslFILEPARSER_TESTDRIVE esl_fileparser.c easel.c
 *    ./test
 * Creates a test file "tmpxxx", then reads it back in.
 *****************************************************************/
#ifdef eslFILEPARSER_TESTDRIVE
#include <stdio.h>
#include <string.h>
#include <easel.h>
#include <esl_fileparser.h>

int 
main(int argc, char **argv)
{
  char *filename = "tmpxxx";
  ESL_FILEPARSER *efp;
  FILE *fp;
  char *tok;
  int   toklen;
  int   status;
  int   ntok;

  /* Create a test file to read.
   */
  if ((fp = fopen(filename, "w")) == NULL)
    esl_fatal("File open failed");
  fprintf(fp, "# Full line comment\n");
  fprintf(fp, "token1  # Trailing comment\n");
  fprintf(fp, "\n");		/* blank line */
  fprintf(fp, "   \n");		/* whitespace line */
  fprintf(fp, "   # sowing comment/whitespace confusion...\n"); 
  fprintf(fp, "token2\ttoken3  token4\n");
  fprintf(fp, "token5");	/* file ends w/ no \n */
  fclose(fp);

  /* Read it back in. Should consist of 5 tokens, all of length 6.
   */
  if ((fp = fopen(filename, "r")) == NULL) 
    esl_fatal("File open failed");
  
  if ((efp = esl_fileparser_Create(fp)) == NULL) 
    esl_fatal("Failed to associate stream with fileparser");
  esl_fileparser_SetCommentChar(efp, '#');
  
  ntok = 0;
  while ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK)
    {
      if (toklen != 6)                   esl_fatal("bad token %s", tok);
      if (strncmp(tok, "token", 5) != 0) esl_fatal("bad token %s", tok);
      ntok++;
    }
  if (status != eslEOF)
    esl_fatal("Abnormal parse termination");
  if (ntok != 5) esl_fatal("bad token number %d\n", ntok);
  
  esl_fileparser_Destroy(efp);
  fclose(fp);
  return 0;
}
#endif /*eslFILEPARSER_TESTDRIVE*/


