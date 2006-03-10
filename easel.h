/* easel.h
 * SRE, Wed Jul  7 09:43:28 2004 [St. Louis]
 * SVN $Id$
 *
 * Core functionality of easel: errors, memory allocations, constants,
 * and configuration for portability.
 */
#ifndef eslEASEL_INCLUDED
#define eslEASEL_INCLUDED

#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>		/* for FILE */
#include <stdarg.h>		/* for va_list */
#ifdef HAVE_STDINT_H
#include <stdint.h>		/* for uint32_t and the like, on C99 systems */
#elif  HAVE_INTTYPES_H
#include <inttypes.h>		/* some systems supposedly put uints here */
#endif

/*****************************************************************
 * Available augmentations.
 * 
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 * 
 * If you grab additional Easel .c files, you can enable their
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 * 
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure. So, if you steal files
 * from an installed library, just set these all back to #undef
 * (depending on which files you have).
 *****************************************************************/

#undef eslAUGMENT_ALPHABET
#undef eslAUGMENT_FILEPARSER
#undef eslAUGMENT_GEV
#undef eslAUGMENT_GUMBEL
#undef eslAUGMENT_HISTOGRAM
/* #undef eslAUGMENT_KEYHASH  */ /* not yet! no keyhash module yet. */
#undef eslAUGMENT_MINIMIZER
#undef eslAUGMENT_MSA		
#undef eslAUGMENT_RANDOM
#undef eslAUGMENT_STATS

/*------------------------------------------------------------------------*/
/*------------ no user serviceable parts below this line -----------------*/
/*------------------------------------------------------------------------*/



/*****************************************************************
 * Version info, set by the ./configure script.
 *****************************************************************/

#undef EASEL_VERSION    
#undef EASEL_DATE      
#undef EASEL_COPYRIGHT 
#undef EASEL_LICENSE   



/*****************************************************************
 * Error handling and allocation macros,
 * including a garbage collection convention.
 *****************************************************************/

/* Macro: ESL_DIE()
 * 
 * The error-throwing convention. Requires that you have an
 * <int status> variable in scope, and that you have a <CLEANEXIT>
 * target for the goto. Allows an error handler to catch an error,
 * but return control (eventually) to the application, without
 * leaking memory inside Easel.
 * 
 * Wrapping macros in <while(0)> loops allows one to write
 *     <if (something) ESL_DIE(code,mesg);>
 * without the trailing semicolon being a null statement
 * after macro expansion. 
 */
#define ESL_DIE(code, mesg)  do {\
     status = code;\
     esl_error(code, __FILE__, __LINE__, mesg);\
     goto CLEANEXIT; }\
     while (0)

/* Macros: ESL_ERROR(), ESL_ERROR_NULL()
 * 
 * Error-throwing conventions for simpler cases,
 * with no garbage collection to worry about. No
 * <int status> variable or <CLEANEXIT> target need
 * to be in scope.
 */
#define ESL_ERROR(code, mesg)  do {\
     esl_error(code, __FILE__, __LINE__, mesg);\
     return code; }\
     while (0)

#define ESL_ERROR_NULL(code, mesg)  do {\
     esl_error(code, __FILE__, __LINE__, mesg);\
     return NULL; }\
     while (0)


/* Macros: ESL_ALLOC(), ESL_RALLOC()
 * 
 * Allocation and reallocation wrappers, including convention for
 * error handling and error recovery. Like ESL_DIE(), they require
 * that you have an <int status> variable in scope, and a <CLEANEXIT>
 * target for the goto. ESL_REALLOC() additionally requires a
 * <void *> ptr to be provided as <tmp>.
 */
#define ESL_ALLOC(p, size) do {\
     if (((p) = malloc(size)) == NULL) {\
       status = eslEMEM;\
       esl_error(eslEMEM, __FILE__, __LINE__, "malloc of size %d failed", size);\
       goto CLEANEXIT;\
     }} while (0)

#define ESL_RALLOC(p, tmp, newsize) do {\
     (tmp) = realloc((p), (newsize));\
     if ((tmp) != NULL) (p) = (tmp);\
     else {\
       status = eslEMEM;\
       esl_error(eslEMEM, __FILE__, __LINE__, "realloc for size %d failed", newsize);\
       goto CLEANEXIT;\
     }} while (0)

/* Macros: ESL_MALLOC(), ESL_REALLOC()
 * [Deprecated: use ESL_ALLOC(), ESL_RALLOC()]
 * 
 * Deprecated allocation wrappers, with no garbage collection.
 */
#define ESL_MALLOC(p, size) do {\
     (p) = malloc(size);\
     if ((p) == NULL) {\
       esl_error(eslEMEM, __FILE__, __LINE__, "malloc failed");\
       return eslEMEM;\
     }} while (0)

#define ESL_REALLOC(p, tmp, newsize) do {\
     (tmp) = realloc((p), (newsize));\
     if ((tmp) != NULL) (p) = (tmp);\
     else {\
       esl_error(eslEMEM, __FILE__, __LINE__, "realloc failed");\
       return eslEMEM;\
     }} while (0)
     
/* Return codes for error handling
 */
#define eslOK              0	/* no error/success             */
#define eslFAIL            1    /* failure                      */
#define eslEOL             2	/* end-of-line (often normal)   */
#define eslEOF             3	/* end-of-file (often normal)   */
#define eslEOD             4 	/* end-of-data (often normal)   */
#define eslEMEM            5	/* malloc or realloc failed     */
#define eslENOTFOUND       6	/* file or key not found        */
#define eslEFORMAT         7	/* file format not correct      */
#define eslEAMBIGUOUS      8    /* an ambiguity of some sort    */
#define eslEDIVZERO        9	/* attempted div by zero        */
#define eslEINCOMPAT      10	/* incompatible parameters      */
#define eslEINVAL         11	/* invalid argument/parameter   */
#define eslESYS           12	/* generic system call failure  */
#define eslECORRUPT       13	/* unexpected data corruption   */
#define eslEINCONCEIVABLE 14    /* "can't happen" error         */
#define eslESYNTAX        15    /* invalid syntax in input data */
#define eslERANGE         16    /* value out of allowed range   */
#define eslEDUP           17    /* saw a duplicate of something */
#define eslECONVERGENCE   18    /* a failure to converge        */      

/* Debugging hooks, w/ three levels (1-3).
 */
#undef eslDEBUGLEVEL

#if eslDEBUGLEVEL >= 1		/* for ESL_DASSERT() macros */
#include <assert.h>
#endif

#if (eslDEBUGLEVEL >= 1)
#define ESL_DPRINTF1(x)  printf x
#define ESL_DASSERT1(x)  assert x
#else
#define ESL_DPRINTF1(x)
#define ESL_DASSERT1(x)
#endif
#if (eslDEBUGLEVEL >= 2)
#define ESL_DPRINTF2(x)  printf x
#define ESL_DASSERT2(x)  assert x
#else
#define ESL_DPRINTF2(x)
#define ESL_DASSERT2(x)
#endif
#if (eslDEBUGLEVEL >= 3)
#define ESL_DPRINTF3(x)  printf x
#define ESL_DASSERT3(x)  assert x
#else
#define ESL_DPRINTF3(x)
#define ESL_DASSERT3(x)
#endif

/* File parsers all contain a fixed length "errbuf" for failure
 * diagnostics. 
 */
#define eslERRBUFSIZE 128





/* Header includes/configuration
 */
#undef HAVE_UNISTD_H
#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H



/* Exact-size integer types.
 * If you see "FIXME" below, the configure script failed to
 * find appropriate types. This shouldn't happen, but if it
 * does, you can try to replace the FIXME(s) with appropriate 
 * type(s).
 */
#undef ESL_UINT16
#undef ESL_UINT32
#undef ESL_UINT64
typedef ESL_UINT16 esl_uint16;
typedef ESL_UINT32 esl_uint32;
typedef ESL_UINT64 esl_uint64;


/* Function portability/configuration
 */
#undef HAVE_POPEN
#undef HAVE_STRCASECMP
#undef HAVE_TIMES

typedef void (*esl_error_handler_f)(int code, char *file, int line,
				    char *format, va_list argp);

extern void esl_error(int code, char *file, int line, char *format, ...);
extern void esl_error_SetHandler(esl_error_handler_f);
extern void esl_error_ResetDefaultHandler(void);
extern void esl_fatal(char *format, ...);

extern void esl_Free2D(void  **p, int dim1);
extern void esl_Free3D(void ***p, int dim1, int dim2);

extern void esl_banner(FILE *fp, char *banner);

extern int  esl_strdup(char *s, int n, char **ret_dup);
extern int  esl_strcat(char **dest, int ldest, char *src, int lsrc);
extern int  esl_fgets(char **buf, int *n, FILE *fp);
extern int  esl_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen);

#ifndef HAVE_STRCASECMP
#ifdef _MSC_VER
#define strcasecmp stricmp
#else
extern int  esl_strcasecmp(const char *s1, const char *s2);
#define strcasecmp esl_strcasecmp
#endif
#endif

extern int  esl_FileExists(char *filename);
extern int  esl_FileTail(char *path, int nosuffix, char *ret_file);
extern int  esl_FileConcat(char *dir, char *file, char **ret_path);
extern int  esl_FileNewSuffix(char *filename, char *sfx, char **ret_newpath);
extern int  esl_FileEnvOpen(char *fname, char *env,
			    FILE **ret_fp, char **ret_path);



/* Making sure TRUE/FALSE are defined, for convenience
 */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/* Some basic constants.
 */
#define eslCONST_E     2.71828182845904523536028747135
#define eslCONST_PI    3.14159265358979323846264338328
#define eslCONST_EULER 0.57721566490153286060651209008
#define eslCONST_GOLD  1.61803399

/* Define <eslINFINITY> portably. Harder than it looks. 
 * We assume we're in an IEEE 754 environment.
 * We assume that HUGE_VAL in a IEEE754 environment is infinity.
 * If we don't have HUGE_VAL set, we assume we can get infinity
 * by division by zero. (But if we don't have HUGE_VAL, we probably
 * have other problems; HUGE_VAL is required by ANSI spec).
 * We can't portably get infinity by overflow (e.g. 1e9999);
 * some compilers (Microsoft) will complain.
 */
#ifdef HUGE_VAL
#define eslINFINITY    HUGE_VAL	 /* assume IEEE754 HUGE_VAL = infinity. ok? */
#else
#define eslINFINITY    (1.0/0.0) /* portable? */
#endif
#define eslNaN         (eslINFINITY/eslINFINITY) /* portably make a IEEE754 NaN */


/* Define some crossovers for numerical approximations.
 */
/* log(1+x) ~ x and  1-e^x = -x approximation.
 * Same threshold appears to be optimal for float or double x. xref STL9/138.
 */
#define eslSMALLX1    5e-9



/* A placeholder for helping w/ portability of filenames/paths.
 * I think, but have not tested, that:
 *   VMS:    #define DIRSLASH ']'
 *   MacOS:  #define DIRSLASH ':'
 *   DOS:    #define DIRSLASH '\\'
 * Setting DIRSLASH correctly is probably not the only thing
 * that would need to be done to port to other OS's, but it's
 * probably a start.
 *
 * The code assumes that '.' is used for file name extensions,
 * such as "foo.bar".
 *
 * This gets used in easel.c's *_File*() functions.
 */
#define eslDIRSLASH '/'           /* UNIX directory paths have /foo/bar */


/* The simple concept of an "inmap" (input map) is shared between
 * the alphabet, msa, and sqio modules, so we put it here to keep
 * these modules separated.
 */
/* Flags in an <inmap>, input map.
 */ 
#define ESL_ILLEGAL_CHAR -2
#define ESL_IGNORED_CHAR -1




#endif /*eslEASEL_INCLUDED*/
