/* A simple token-based file parsing system.
 * 
 * SRE, Tue Jul 13 14:40:35 2004 [St. Louis]
 * SVN $Id$
 */

#ifndef ESL_FILEPARSER_INCLUDED
#define ESL_FILEPARSER_INCLUDED

#include <stdio.h>
#include "easel.h"

typedef struct {
  FILE *fp;			/* open file pointer, for reading                  */
  char *buf;			/* current line; will be modified by esl_strtok(). */
  int   buflen;			/* current allocated length of buf                 */
  char *s;			/* used by esl_strtok(); current position in buf.  */
  char  commentchar;		/* often '#'                                       */

  char *tok;			/* _NextLine() may remember a token...             */
  int   toklen;			/* ... and its length                              */

  int   linenumber;		/* what line is loaded into buf; 1..nlines         */
  char  errbuf[eslERRBUFSIZE];  /* for holding error diagnostics                   */
} ESL_FILEPARSER;

extern int  esl_fileparser_Open(const char *filename, ESL_FILEPARSER **ret_efp);
extern ESL_FILEPARSER *esl_fileparser_Create(FILE *fp);
extern int  esl_fileparser_SetCommentChar(ESL_FILEPARSER *efp, char c);
extern int  esl_fileparser_NextLine(ESL_FILEPARSER *efp);
extern int  esl_fileparser_GetToken(ESL_FILEPARSER *efp, 
				   char **opt_tok, int *opt_toklen);
extern int  esl_fileparser_GetTokenOnLine(ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
extern void esl_fileparser_Destroy(ESL_FILEPARSER *efp);
extern void esl_fileparser_Close(ESL_FILEPARSER *efp);

#endif /*ESL_FILEPARSER_INCLUDED */
