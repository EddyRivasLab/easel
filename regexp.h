/* regexp.h
 * Regular expression matching on strings.
 * 
 * SRE, Sun Jan  2 10:52:34 2005 [Zaragoza]
 * SVN $Id$
 ******************************************************************
 * The regexp module is a wrapper around a modified version of Henry
 * Spencer's regex library. Spencer's copyright notice appears below,
 * after my wrappers, prefacing the section that includes his code. I
 * believe you can obtain the original code from:
 *    ftp://ftp.zoo.toronto.edu/pub/bookregex.tar.Z
 * Thanks, Henry!
 *****************************************************************
 */    

/* ESL_REX_NSUBEXP specifies the maximum number of () expressions
 * in a regexp. The whole regexp counts as one, so 16 allows for 
 * parsing out up to 15 substrings from the match.
 */
#define ESL_REX_NSUBEXP 16



/* The esl__regexp structure is from the original Spencer code.
 * It's wrapped by the ESL_REGEXP structure, below.
 */
struct esl__regexp {
  char *startp[ESL_REX_NSUBEXP]; /* ptrs to starts of submatches on target string */
  char *endp[ESL_REX_NSUBEXP];   /* ptrs to end pts of submatches on target string */
  char regstart;		 /* Internal use only. */
  char reganch;		         /* Internal use only. */
  char *regmust;		 /* Internal use only. */
  int regmlen;		         /* Internal use only. */
  char program[1];	         /* Unwarranted chumminess with compiler. */  
};


struct esl_regexp_s {
  struct esl__regexp *r;	 /* a compiled regexp */

  /* We maintain internal allocated buffers to temporarily 
   * hold matched substrings.
   */
  char *submatch[ESL_REX_NSUBEXP];
  int   alloc[ESL_REX_NSUBEXP];
};
typedef struct esl_regexp_s ESL_REGEXP;



