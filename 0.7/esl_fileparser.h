/* esl_fileparser.h
 * A simple token-based file parsing system.
 * 
 * SRE, Tue Jul 13 14:40:35 2004 [St. Louis]
 * SVN $Id$
 */

#ifndef ESL_PARSE_INCLUDED
#define ESL_PARSE_INCLUDED

#include <stdio.h>
#include <easel.h>

typedef struct {
  FILE *fp;
  char *buf;
  int   buflen;
  char *s;
  char  commentchar;		/* often '#' */

  int   linenumber;		/* what line is loaded into buf; 1..nlines */
  char  errbuf[eslERRBUFSIZE];  /* for holding error diagnostics           */
} ESL_FILEPARSER;

extern ESL_FILEPARSER *esl_fileparser_Create(FILE *fp);
extern int esl_fileparser_SetCommentChar(ESL_FILEPARSER *efp, char c);
extern int esl_fileparser_GetToken(ESL_FILEPARSER *efp, 
				   char **ret_tok, int *ret_toklen);
extern int esl_fileparser_Destroy(ESL_FILEPARSER *efp);

#endif /* ESL_PARSE_INCLUDED */
