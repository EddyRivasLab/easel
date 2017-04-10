/* esl_config.h.  Generated from esl_config.h.in by configure.  */
/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by autoconf.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 */
#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

/* Version info.
 */
#define EASEL_VERSION "0.44"
#define EASEL_DATE "January 2017"
#define EASEL_COPYRIGHT "Copyright (C) 2017 Howard Hughes Medical Institute"
#define EASEL_LICENSE "Freely distributed under the BSD open source license."

/* Large file support
 * Must precede any header file inclusion.
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* Debugging verbosity (0=none;3=most verbose)
 */
#define eslDEBUGLEVEL 1

/* System headers
 */
#define HAVE_ENDIAN_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_STRINGS_H 1

#define HAVE_SYS_PARAM_H 1
#define HAVE_SYS_SYSCTL_H 1

#define HAVE_EMMINTRIN_H 1
#define HAVE_PMMINTRIN_H 1
#define HAVE_XMMINTRIN_H 1

/* #undef HAVE_ALTIVEC_H */

/* Types
 */
/* #undef WORDS_BIGENDIAN */
/* #undef int8_t */
/* #undef int16_t */
/* #undef int32_t */
/* #undef int64_t */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Optional packages
 */
/* #undef HAVE_LIBGSL */

/* Optional parallel implementation support
 */
#define HAVE_SSE2 1
/* #undef HAVE_VMX */
/* #undef HAVE_MPI */
#define HAVE_PTHREAD 1
/* #undef HAVE_NEON */
/* #undef HAVE_ARM64 */
/* #undef HAVE_SSE2_CAST */

/* Programs */
#define HAVE_GZIP 1

/* Functions */
#define HAVE_MKSTEMP 1
#define HAVE_POPEN 1
#define HAVE_PUTENV 1
#define HAVE_STRCASECMP 1
#define HAVE_STRSEP 1
#define HAVE_TIMES 1
#define HAVE_GETPID 1
#define HAVE_SYSCTL 1
#define HAVE_SYSCONF 1
#define HAVE_GETCWD 1
#define HAVE_CHMOD 1
#define HAVE_STAT 1
#define HAVE_FSTAT 1
#define HAVE_ERFC 1
#define HAVE_FSEEKO 1

#define HAVE_FUNC_ATTRIBUTE_NORETURN 1
#define HAVE_FUNC_ATTRIBUTE_FORMAT 1
 
/* Function behavior */
#define eslSTOPWATCH_HIGHRES

/*****************************************************************
 * Available augmentations.
 * 
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 * 
 * If you grab additional Easel .c files, you can enable any
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 * 
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure, plus the eslLIBRARY flag
 * which means the full library with all augmentations is
 * available. So, if you steal files from an installed library, just
 * set these all back to #undef (depending on which files you have).
 *****************************************************************/
#define eslLIBRARY 1

#ifndef eslLIBRARY
#define eslAUGMENT_ALPHABET 1
#define eslAUGMENT_NCBI 1
#define eslAUGMENT_DMATRIX 1
#define eslAUGMENT_FILEPARSER 1
#define eslAUGMENT_GEV 1
#define eslAUGMENT_GUMBEL 1
#define eslAUGMENT_HISTOGRAM 1
#define eslAUGMENT_KEYHASH 1
#define eslAUGMENT_MINIMIZER 1
#define eslAUGMENT_MSA 1
#define eslAUGMENT_RANDOM 1
/* #undef eslAUGMENT_RANDOMSEQ */
#define eslAUGMENT_SSI 1
#define eslAUGMENT_STATS 1
#endif

#ifdef eslLIBRARY
#define eslAUGMENT_ALPHABET 1
#define eslAUGMENT_NCBI 1
#define eslAUGMENT_DMATRIX 1
#define eslAUGMENT_FILEPARSER 1
#define eslAUGMENT_GEV 1
#define eslAUGMENT_GUMBEL 1
#define eslAUGMENT_HISTOGRAM 1
#define eslAUGMENT_KEYHASH 1
#define eslAUGMENT_MINIMIZER 1
#define eslAUGMENT_MSA 1
#define eslAUGMENT_RANDOM 1
#define eslAUGMENT_RANDOMSEQ
#define eslAUGMENT_SSI 1
#define eslAUGMENT_STATS 1
#endif


#endif /*eslCONFIG_INCLUDED*/

