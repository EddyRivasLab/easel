/* easel.c
 * SRE, Tue Oct 28 08:29:17 2003 [St. Louis]
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <easel/easel.h>

esl_error_handler_f esl_error_handler = NULL;

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
    fprintf(stderr, "Easel fatal error: ");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\nExiting at file %s, line %d.\n", file, line);
    fflush(stderr);
    abort();
  }
}
  

/* Function: esl_Die()
 * 
 * Purpose:  Print an error message and die. The arguments
 *           are formatted exactly like arguments to printf().
 *           
 * Return:   None. Exits the program.
 */          
/* VARARGS0 */
void
esl_Die(char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
                                /* exit  */
  exit(ESL_EUNKNOWN);
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/  
