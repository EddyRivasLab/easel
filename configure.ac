# configure.ac
# Process this file with autoconf to produce the Easel configuration script.
#
# reminders to save re-reading autoconf manual for the n'th time:
#   - output variables:
#        are defined here as normal shell variables, e.g. FOO="my string"
#        are made into output variables by calling AC_SUBST(FOO)
#        any occurrence of @FOO@ in an output file is then substituted
#        I think this only happens in files we assign w/ AC_CONFIG_FILES;
#        that's the Makefile.in's.
#
#   - C preprocessor symbols:
#        are defined here by calling AC_DEFINE(FOO) or AC_DEFINE(FOO, [42])
#        then #undef FOO lines in easel.h.in become #define FOO or #define FOO 42
#        I think this only happens in header files that we assign
#        w/ AC_CONFIG_HEADERS -- which means, easel.h.in
#
# SRE, Mon Feb 14 10:51:07 2005
# SVN $Id$
#
# xref autoconf macro archive at:  http://www.gnu.org/software/ac-archive/
#
# GNU recommends the following order:
#   1. autoconf requirements
#   2. AC_INIT
#   3. info on the package
#   4. checks for programs
#   5. checks for libraries
#   6. checks for header files
#   7. checks for types
#   8. checks for structures
#   9. checks for compiler characteristics
#  10. checks for library functions
#  11. checks for system services
#  12. AC_CONFIG_FILES
#  13. AC_OUTPUT


# AC_INIT information
#
# The four AC_INIT args set the following output variables and preprocessor 
# symbols:
#     PACKAGE_NAME      <package>     e.g. "Easel"
#     PACKAGE_VERSION   <version>     e.g. "0.1"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddy@genetics.wustl.edu"
#     PACKAGE_TARNAME   <tarname>     e.g. "easel"
# From them, AC_INIT automatically derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "Easel 0.1"
# Then we define some of our own:
#     EASEL_RELCODE     unique package/version SVN tag w/o .'s: e.g. "easel0_1"
#     EASEL_DATE        release date: e.g. "February 2005"
#     EASEL_COPYRIGHT   one-line copyright string
#     EASEL_LICENSE     one-line license string
#     EASEL_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#     EASEL_URL         URL home for Easel.
#
# Because Easel is designed to be a library and must coexist with
# packages that include it as a subdirectory (HMMER, Infernal...),
# we don't want to use AC_INIT's PACKAGE_ variables as preprocessor
# symbols to put version info into executables; we'll get name clashes,
# plus we might want to have both Easel version info and main package
# version info. So, we use the PACKAGE_ stuff to make a
# new preprocessor symbol of our own:
#     EASEL_VERSION    e.g. "1.9a"
#
# The preprocessor symbols (version info) goes into easel.h.
# Output variables are used in Makefiles.
#


AC_PREREQ(2.57)
AC_INIT(Easel, 0.1, eddy@genetics.wustl.edu, easel)
AC_MSG_NOTICE([Configuring the Easel library for your system.])

################################################################
# 3. Info on the package
################################################################

EASEL_RELCODE="easel0_1"
EASEL_DATE="February 2005"
EASEL_COPYRIGHT="Copyright (C) 2004-2005 HHMI/Washington University School of Medicine"
EASEL_LICENSE="Licensed and distributed under the Creative Commons Attribution License."
EASEL_LICENSETAG=cc
EASEL_VERSION=$PACKAGE_VERSION
EASEL_URL="http://selab.wustl.edu/easel/"

# Output variables.
# (AC_OUTPUT replaces instances of @var@ in input files.)
AC_SUBST(EASEL_RELCODE)
AC_SUBST(EASEL_DATE)
AC_SUBST(EASEL_COPYRIGHT)
AC_SUBST(EASEL_LICENSE)
AC_SUBST(EASEL_LICENSETAG)
AC_SUBST(EASEL_VERSION)
AC_SUBST(EASEL_URL)

# Preprocessor symbols.
# adds a -Dvar=value to @DEFS@;
# replace #undef's in easel.h.in
AC_DEFINE_UNQUOTED(EASEL_DATE,      "$EASEL_DATE")
AC_DEFINE_UNQUOTED(EASEL_COPYRIGHT, "$EASEL_COPYRIGHT")
AC_DEFINE_UNQUOTED(EASEL_LICENSE,   "$EASEL_LICENSE")
AC_DEFINE_UNQUOTED(EASEL_VERSION,   "$EASEL_VERSION")

################################################################
# 3b. ./configure command line options, for optional config
#     of the package
################################################################

# --enable-modular       - configure for use of standalone modules
#
# Default is to enable all augmentations and compile entire library.
# With --enable-modular, any desired augmentations must be set
# at compile-time by the user.
#
AC_ARG_ENABLE(modular, 
  [  --enable-modular        Compile and prepare modules for standalone use], 
  [ 
    AC_MSG_NOTICE([Compiling for standalone (independent) use of modules.])
  ], 
  [
    AC_MSG_NOTICE([Compiling the full Easel library, with all augmentations.])
    AC_DEFINE(eslAUGMENT_ALPHABET)
    AC_DEFINE(eslAUGMENT_FILEPARSER)
    AC_DEFINE(eslAUGMENT_GEV)
    AC_DEFINE(eslAUGMENT_GUMBEL)
    AC_DEFINE(eslAUGMENT_HISTOGRAM)
    AC_DEFINE(eslAUGMENT_KEYHASH)
    AC_DEFINE(eslAUGMENT_MINIMIZER)
    AC_DEFINE(eslAUGMENT_MSA)
    AC_DEFINE(eslAUGMENT_RANDOM)
    AC_DEFINE(eslAUGMENT_STATS)
  ])

# --enable-debugging=x    - set debugging level to <x> (1-3)
#
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the eslDEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) (3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(eslDEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(eslDEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
         ;;
esac])


################################################################
# 4. Checks for programs
################################################################

# Set output variable "RANLIB"
AC_PROG_RANLIB

# Set output variable "AR"
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])

# GNU make check; format of makefile dependency lines for executables.
# original from John Darrington <j.darrington@elvis.murdoch.edu.au> 
# w/ heavy modifications.
#
# We need this because GNU make and SYSV make use different systems
# specifying variables for dependencies: $$@ in sysv, %: %.o in GNU.
# Would love to hear a better way of doing this.
# 
# I use two different conventions in my Makefiles. Sometimes 
# executable "foo" has a file "foo.c" - this is the HMMER convention.
# Sometimes executable "foo" has a file "foo_main.c" - this is
# the SQUID convention. The configure script sets the
# EXEC_DEPENDENCY appropriately: here, HMMER style.
#
# This creates a function CHECK_GNU_MAKE, which we immediately call.
# It sets an output variable EXEC_DEPENDENCY. 
# This is used in the src/Makefile.in.
#
AC_DEFUN(CHECK_GNU_MAKE,[ 
  AC_MSG_CHECKING(whether your make is GNU make)
  foundGNUmake='nope, assuming sysv make.' ;
  EXEC_DEPENDENCY=[\$\$\@.o] ;
  if ( make --version nothing 2> /dev/null | grep GNU > /dev/null ) ;  then
     foundGNUmake='yes, it is.' ;
     EXEC_DEPENDENCY='%: %.o' ;
  fi
  AC_MSG_RESULT($foundGNUmake)
  AC_SUBST(EXEC_DEPENDENCY)
])
CHECK_GNU_MAKE



# 5. Checks for libraries - none


# 6. Checks for header files.
#    Defines preprocessor symbols like HAVE_UNISTD_H
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([inttypes.h])


# 7. Checks for types.
#    ssi module needs exact-size integer types.
#    Sets output variables ESL_UINT16, ESL_UINT32, and ESL_UINT64.
#    We try to use ISO C uint*_t's, if stdint.h is present; 
#    else, we start guessing different ANSI C unsigned types.
#
AC_CHECK_SIZEOF(unsigned short)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(unsigned long long)
AC_CHECK_TYPE(uint16_t,
   [AC_DEFINE(ESL_UINT16, uint16_t)],
   [
     if test "$ac_cv_sizeof_unsigned_short" -eq 2; then
       AC_DEFINE(ESL_UINT16, [unsigned short])
     elif test "$ac_cv_sizeof_unsigned_int" -eq 2; then
       AC_DEFINE(ESL_UINT16, [unsigned int])
     else	
       AC_MSG_WARN("You will need to manually edit easel.h to typedef esl_uint16.")
     fi
   ]
)
AC_CHECK_TYPE(uint32_t, 
   [AC_DEFINE(ESL_UINT32, uint32_t)],
   [
     if test "$ac_cv_sizeof_unsigned_int" -eq 4; then
       AC_DEFINE(ESL_UINT32, [unsigned int])
     elif test "$ac_cv_sizeof_unsigned_long" -eq 4; then
       AC_DEFINE(ESL_UINT32, [unsigned long])
     elif test "$ac_cv_sizeof_unsigned_short" -eq 4; then
       AC_DEFINE(ESL_UINT32, [unsigned short])
     else	
       AC_MSG_WARN("You will need to manually edit easel.h to typedef esl_uint32.")
     fi
   ]
)
AC_CHECK_TYPE(uint64_t, 
   [AC_DEFINE(ESL_UINT64, uint64_t)],
   [
     if test "$ac_cv_sizeof_unsigned_long" -eq 8; then
       AC_DEFINE(ESL_UINT64, [unsigned long])
     elif test "$ac_cv_sizeof_unsigned_long_long" -eq 8; then
       AC_DEFINE(ESL_UINT64, [unsigned long long])
     elif test "$ac_cv_sizeof_unsigned_int" -eq 8; then
       AC_DEFINE(ESL_UINT64, [unsigned int])
     else	
       AC_MSG_WARN("You will need to manually edit easel.h to typedef esl_uint64.")
     fi
   ]
)




# 8. Checks for structures - none


################################################################
# 9. Checks for compiler characteristics
################################################################

# We can probably get fancier & more comprehensive about setting
# optimization flags here. The IBM xlc settings are at request of IBM.
# 
AC_PROG_CC
if test "$CC" = xlc; then
   CFLAGS=$CFLAGS" -O3 -qarch=auto -qtune=auto"
fi




################################################################
# 10. Checks for library functions: define HAVE_FOO
################################################################

AC_CHECK_FUNCS(popen)
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(times)
AC_FUNC_FSEEKO

# 11. Checks for system services 
AC_SYS_LARGEFILE




# 12. Write out easel.h header and the Makefiles
AC_CONFIG_HEADERS([easel.h])
AC_CONFIG_HEADERS([esl_ssi.h])
AC_CONFIG_FILES([miniapps/Makefile])
AC_CONFIG_FILES([Makefile])


# 13. voila!
AC_OUTPUT
