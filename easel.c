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




/*****************************************************************
 * @LICENSE@
 *****************************************************************/  
