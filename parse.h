/* parse.h
 * 
 * SRE, Tue Jul 13 14:40:35 2004 [St. Louis]
 * SVN $Id$
 */

#ifndef ESL_PARSE_INCLUDED
#define ESL_PARSE_INCLUDED

#include <stdio.h>

struct esl_fileparser_s {
  FILE *fp;
  char *buf;
  int   buflen;
  char *s;
  char  commentchar;		/* often '#' */
};
typedef struct esl_fileparser_s ESL_FILEPARSER;

extern int esl_parse_fgets(char **buf, int *n, FILE *fp);
extern int esl_parse_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen);
extern int esl_fileparse_create(FILE *fp, ESL_FILEPARSER **ret_efp);
extern int esl_fileparse_set_commentchar(ESL_FILEPARSER *efp, char c);
extern int esl_fileparse_nextline(ESL_FILEPARSER *efp);
extern int esl_fileparse_token(ESL_FILEPARSER *efp, char **ret_tok, int *ret_toklen);
extern int esl_fileparse_free(ESL_FILEPARSER *efp);

#endif /* ESL_PARSE_INCLUDED */
