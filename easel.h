/* easel.h
 * SRE, Wed Jul  7 09:43:28 2004 [St. Louis]
 * SVN $Id$
 * 
 * Core functionality of easel: errors and memory allocations.
 */
#ifndef ESL_EASEL_INCLUDED
#define ESL_EASEL_INCLUDED

#include <stdarg.h>

/* Error handling.
 * Originally modeled on GNU Scientific Library (GSL).
 */
#define ESL_ERROR(code, mesg)  do {\
     esl_error(code, __FILE__, __LINE__, mesg);\
     return code; }\
     while (0)

#define ESL_ERROR_NULL(mesg)  do {\
     esl_error(ESL_EMEM, __FILE__, __LINE__, mesg);\
     return NULL; }\
     while (0)

typedef void (*esl_error_handler_f)(int code, char *file, int line, char *format, va_list argp);
extern esl_error_handler_f esl_error_handler;


#define ESL_MALLOC(x)     esl_malloc(__FILE__, __LINE__, (x))
#define ESL_REALLOC(x,y)  esl_realloc(__FILE__, __LINE__, (x), (y))   






/* error codes used in Easel.
 */
#define ESL_OK         0	/* no error                     */
#define ESL_EOL        1	/* end-of-line (often normal)   */
#define ESL_EOF        2	/* end-of-file (often normal)   */
#define ESL_EMEM       3	/* malloc or realloc failed     */
#define ESL_ENOFILE    4	/* file not found               */
#define ESL_EFORMAT    5	/* file format not recognized   */
#define ESL_EPARAM     6	/* bad parameter passed to func */
#define ESL_EDIVZERO   7	/* attempted div by zero        */
#define ESL_EINCOMPAT  8	/* incompatible parameters      */
#define ESL_EINVAL     9	/* invalid argument             */
#define ESL_ETESTFAIL  10	/* calculated test failure      */
#define ESL_EUNKNOWN   127      /* generic error, unidentified  */

/* Making sure TRUE/FALSE are defined, for convenience
 */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Function declarations 
 */
extern void esl_error(int code, char *file, int line, char *format, ...);
extern void esl_error_SetHandler(void (*handler)(int code, char *file, int line, char *format, va_list argp));
extern void esl_error_ResetDefaultHandler(void);


#endif /*ESL_EASEL_INCLUDED*/
