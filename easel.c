/* easel.c
 * SRE, Tue Oct 28 08:29:17 2003 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <easel/easel.h>

static esl_error_handler_f esl_error_handler = NULL;

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
    fprintf(stderr, "Easel fatal error:\n");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\nAborted at file %s, line %d.\n", file, line);
    fflush(stderr);
    abort();
  }
}
  
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

void *
esl_malloc(char *file, int line, size_t size)
{
  void *ptr;

  SQD_DPRINTF3(("MALLOC: %d bytes (file %s line %d)\n", size, file, line));
  if ((ptr = malloc (size)) == NULL)
    Die("malloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/  
