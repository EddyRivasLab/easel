/* easel.h
 *
 * Core functionality of easel: errors, memory allocations, constants,
 * and configuration for portability.
 *
 * SRE, Wed Jul  7 09:43:28 2004 [St. Louis]
 * SVN $Id$
 */
#ifndef eslEASEL_INCLUDED
#define eslEASEL_INCLUDED

#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>		/* for FILE */
#include <stdarg.h>		/* for va_list */
#ifdef HAVE_STDINT_H
#include <stdint.h>		/* for uint32_t and the like (C99) */
#elif  HAVE_INTTYPES_H
#include <inttypes.h>		/* some systems allegedly put uints here */
#endif

/*****************************************************************
 * Macros implementing Easel's error handling conventions
 *****************************************************************/
/* ESL_FAIL()       - return an error message, without cleanup.
 * ESL_XFAIL()      - return an error message, with cleanup.
 * ESL_EXCEPTION()  - throwing an exception, without cleanup.
 * ESL_XEXCEPTION() - throwing an exception, with cleanup.
 * ESL_FWD()        - percolate an exception outwards, without cleanup.
 * ESL_XFWD()       - percolate an exception outwards, with cleanup.
 * 
 * The X versions (with cleanup) require the caller to have an
 * <int status> variable and a <ERROR:> goto target in scope.
 *
 * Wrapping these macros in <while(0)> loops allows a statement:
 *       if (something) ESL_XEXCEPTION(code,mesg);
 * without the trailing semicolon becoming a null statement after 
 * macro expansion.
 */
/*::cexcerpt::error_macros::begin::*/
#define ESL_FAIL(code, errbuf, ...) do {\
     if (errbuf != NULL) sprintf(errbuf, __VA_ARGS__);\
     return code; }\
     while (0)

#define ESL_XFAIL(code, errbuf, ...) do {\
     status = code;\
     if (errbuf != NULL) sprintf(errbuf, __VA_ARGS__);\
     goto ERROR; }\
     while (0)

#define ESL_EXCEPTION(code, ...)  do {\
     esl_exception(code, __FILE__, __LINE__, __VA_ARGS__);\
     return code; }\
     while (0)

#define ESL_XEXCEPTION(code, ...)  do {\
     status = code;\
     esl_exception(code, __FILE__, __LINE__, __VA_ARGS__);\
     goto ERROR; }\
     while (0)

#define ESL_FWD(code)  do {\
     esl_exception(code, __FILE__, __LINE__, "[percolated]");\
     return code; }\
     while (0);

#define ESL_XFWD(code) do {\
     status = code;\
     esl_exception(code, __FILE__, __LINE__, "[percolated]");\
     goto ERROR; }\
     while (0);
/*::cexcerpt::error_macros::end::*/


/* ESL_ALLOC(), ESL_RALLOC():
 * 
 * Allocation and reallocation wrappers.
 * Both require <int status> in scope, and <ERROR:> goto target.
 * ESL_RALLOC() also requires <void *> ptr to be provided as <tmp>.
 */
/*::cexcerpt::alloc_macros::begin::*/
#define ESL_ALLOC(p, size) do {\
     if (((p) = malloc(size)) == NULL) {\
       status = eslEMEM;\
       esl_exception(eslEMEM, __FILE__, __LINE__, "malloc of size %d failed", size);\
       goto ERROR;\
     }} while (0)

#define ESL_RALLOC(p, tmp, newsize) do {\
     if ((p) == NULL) { (tmp) = malloc(newsize);         }\
     else             { (tmp) = realloc((p), (newsize)); }\
     if ((tmp) != NULL) (p) = (tmp);\
     else {\
       status = eslEMEM;\
       esl_exception(eslEMEM, __FILE__, __LINE__, "realloc for size %d failed", newsize);\
       goto ERROR;\
     }} while (0)
/*::cexcerpt::alloc_macros::end::*/

     
/* Return codes for error handler
 */
/*::cexcerpt::statuscodes::begin::*/
#define eslOK              0    /* no error/success             */
#define eslFAIL            1    /* failure                      */
#define eslEOL             2    /* end-of-line (often normal)   */
#define eslEOF             3    /* end-of-file (often normal)   */
#define eslEOD             4    /* end-of-data (often normal)   */
#define eslEMEM            5    /* malloc or realloc failed     */
#define eslENOTFOUND       6    /* file or key not found        */
#define eslEFORMAT         7    /* file format not correct      */
#define eslEAMBIGUOUS      8    /* an ambiguity of some sort    */
#define eslEDIVZERO        9    /* attempted div by zero        */
#define eslEINCOMPAT      10    /* incompatible parameters      */
#define eslEINVAL         11    /* invalid argument/parameter   */
#define eslESYS           12    /* generic system call failure  */
#define eslECORRUPT       13    /* unexpected data corruption   */
#define eslEINCONCEIVABLE 14    /* "can't happen" error         */
#define eslESYNTAX        15    /* invalid syntax in input data */
#define eslERANGE         16    /* value out of allowed range   */
#define eslEDUP           17    /* saw a duplicate of something */
#define eslECONVERGENCE   18    /* a failure to converge        */      
#define eslECONTRACT      19    /* an API contract violation    */
#define eslENORESULT      20    /* no result was obtained       */
/*::cexcerpt::statuscodes::end::*/

/* File parsers all contain a fixed length "errbuf" for failure
 * diagnostics. 
 */
#define eslERRBUFSIZE 128


/* Debugging hooks, w/ three levels (1-3).
 */
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


typedef void (*esl_exception_handler_f)(int code, char *file, int line,
					char *format, va_list argp);

extern void esl_exception(int code, char *file, int line, char *format, ...);
extern void esl_exception_SetHandler(esl_exception_handler_f);
extern void esl_exception_ResetDefaultHandler(void);
extern void esl_fatal(char *format, ...);
extern void esl_nonfatal_handler(int code, char *file, int line, char *format, va_list argp);

extern void esl_Free2D(void  **p, int dim1);
extern void esl_Free3D(void ***p, int dim1, int dim2);

extern void esl_banner(FILE *fp, char *banner);

extern int  esl_DCompare(double a, double b, double tol);
extern int  esl_FCompare(float  a, float  b, float  tol);

extern int  esl_strdup(char *s, int n, char **ret_dup);
extern int  esl_strcat(char **dest, int ldest, char *src, int lsrc);
extern int  esl_fgets(char **buf, int *n, FILE *fp);
extern int  esl_strtok(char **s, char *delim, char **ret_tok, int *ret_toklen);
extern int  esl_strchop(char *s, int n);
#ifndef HAVE_STRCASECMP
#ifdef _MSC_VER
#define strcasecmp stricmp
#else
extern int  esl_strcasecmp(const char *s1, const char *s2);
#define strcasecmp esl_strcasecmp
#endif
#endif

extern int  esl_FileExists(char *filename);
extern int  esl_FileTail(char *path, int nosuffix, char **ret_file);
extern int  esl_FileConcat(char *dir, char *file, char **ret_path);
extern int  esl_FileNewSuffix(char *filename, char *sfx, char **ret_newpath);
extern int  esl_FileEnvOpen(char *fname, char *env,
			    FILE **ret_fp, char **ret_path);
extern int  esl_tmpfile(char *template, FILE **ret_fp);
extern int  esl_tmpfile_named(char *template, FILE **ret_fp);

/* Making sure TRUE/FALSE are defined, for convenience
 */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Some basic constants. 
 * Assuming IEEE754 math with 64-bit doubles (53-bit mantissas), we 
 * want 17 significant decimal digits in our constants. More is
 * a waste (but we do it anyway).
 */
#define eslCONST_E     2.71828182845904523536028747135
#define eslCONST_PI    3.14159265358979323846264338328
#define eslCONST_EULER 0.57721566490153286060651209008
#define eslCONST_GOLD  1.618033988749894
#define eslCONST_LOG2  0.69314718055994529

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


/* Digitized sequences.
 * Most of this support is in the alphabet module, but we externalize it
 * into the easel foundation because ESL_INMAP is used in unaugmented
 * sqio, msa modules.
 * 
 * A digital sequence residue (ESL_DSQ) is an unsigned 8-bit type
 * (0..255).  A valid digital residue has a value in the range 0..127
 * (Easel can represent alphabets of up to 128 different characters).
 * Values 128..255 are reserved for flags.
 *
 * An "inmap" is ESL_DSQ[128], or *ESL_DSQ allocated for 128 values;
 * it is a many-to-one construct for mapping 7-bit ASCII chars (in
 * range 0..127) either to new ASCII chars (in the case of raw
 * sequence input in sqio, msa) or to digital codes (in the alphabet
 * module).  Valid mapped values are 0..127; any value in range
 * 128..255 is some kind of flag.
 */
typedef uint8_t ESL_DSQ;
#define eslDSQ_SENTINEL 255	/* sentinel bytes 0,L+1 in a dsq */
#define eslDSQ_ILLEGAL  254	/* input symbol is unmapped and unexpected */
#define eslDSQ_IGNORED  253     /* input symbol is unmapped and ignored */

/* If you try to test sym > 0 && sym <= 127 below, instead of isascii(sym),
 * you'll get a compiler warning for an always-successful test regardless
 * of whether a char is signed or unsigned. So we trust that isascii() is
 * doing the Right Thing.
 */
#define esl_inmap_IsValid(inmap, sym)  (isascii(sym) && (inmap)[(int)sym] <= 127)

/* Some generic macros for swapping, min, and max.
 */
#define ESL_SWAP(x, y, type)  do { type tmpxyz = (x); (x) = (y); (y) = tmpxyz; } while (0)
#define ESL_MIN(a,b)          (((a)<(b))?(a):(b))
#define ESL_MAX(a,b)          (((a)>(b))?(a):(b))

#endif /*eslEASEL_INCLUDED*/
