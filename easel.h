/* easel.h
 * SRE, Wed Jul  7 09:43:28 2004 [St. Louis]
 * SVN $Id$
 * 
 * Core functionality of easel: errors and memory allocations.
 */
#ifndef ESL_EASEL_INCLUDED
#define ESL_EASEL_INCLUDED

#include <stdlib.h>
#include <stdarg.h>


/* Error handling.
 * Originally modeled on GNU Scientific Library (GSL).
 *
 * Wrapping in the while(0) loop allows one to write
 *  if (something) ESL_ERROR(code,mesg); 
 * without the trailing semicolon being a null statement
 * after macro expansion.
 */
#define ESL_ERROR(code, mesg)  do {\
     esl_error(code, __FILE__, __LINE__, mesg);\
     return code; }\
     while (0)

#define ESL_ERROR_NULL(code, mesg)  do {\
     esl_error(code, __FILE__, __LINE__, mesg);\
     return NULL; }\
     while (0)



#define ESL_OK              0	/* no error/success/TRUE        */
#define ESL_EOL             1	/* end-of-line (often normal)   */
#define ESL_EOF             2	/* end-of-file (often normal)   */
#define ESL_EOD             3 	/* end-of-data (often normal)   */
#define ESL_EMEM            4	/* malloc or realloc failed     */
#define ESL_ENOTFOUND       5	/* file or key not found        */
#define ESL_EFORMAT         6	/* file format not correct      */
#define ESL_EAMBIGUOUS      7   /* an ambiguity of some sort    */   
#define ESL_EDIVZERO        8	/* attempted div by zero        */
#define ESL_EINCOMPAT       9	/* incompatible parameters      */
#define ESL_EINVAL         10	/* invalid argument/parameter   */
#define ESL_ESYS           11	/* generic system call failure  */
#define ESL_ECORRUPT       12	/* unexpected data corruption   */
#define ESL_EINCONCEIVABLE 13   /* "can't happen" error         */
#define ESL_ESYNTAX        14   /* invalid syntax in input data */
#define ESL_ERANGE         15   /* value out of allowed range   */


typedef void (*esl_error_handler_f)(int code, char *file, int line, char *format, va_list argp);
extern esl_error_handler_f esl_error_handler;

extern void esl_error(int code, char *file, int line, char *format, ...);
extern void esl_error_SetHandler(esl_error_handler_f);
extern void esl_error_ResetDefaultHandler(void);

extern void  esl_Free2D(void  **p, int dim1);
extern void  esl_Free3D(void ***p, int dim1, int dim2);

extern char *esl_strdup(char *s, int n);
extern int   esl_fgets(char **buf, int *n, FILE *fp);
extern int   esl_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen);


/* See esl_msa_Expand() for first use example.
 */
#define ESL_REALLOC(p, tmp, newsize) {\
     (tmp) = realloc((p), (newsize));\
     if ((tmp) != NULL) (p) = (tmp);\
     else {\
       esl_error(ESL_EMEM, __FILE__, __LINE__, "realloc failed");\
       return ESL_EMEM;\
     }}


/* Making sure TRUE/FALSE are defined, for convenience
 */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Some basic constants.
 */
#define ESL_CONSTANT_E    2.71828182845904523536028747135
#define ESL_CONSTANT_PI   3.14159265358979323846264338328

#endif /*ESL_EASEL_INCLUDED*/
