# configure.ac
# Process this file with autoconf to produce the Easel configuration script.
#
# SRE, Mon Feb 14 10:51:07 2005
# SVN $Id$
#
# xref autoconf macro archive at:  http://www.gnu.org/software/ac-archive/
#


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
AC_MSG_NOTICE([configuring the Easel library for your system.])

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


# Checks for programs.
#
AC_PROG_CC
if test "$CC" = xlc; then
   CFLAGS=$CFLAGS" -O3 -qarch=auto -qtune=auto"
fi
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])


# msa.c and sqio.c can read gzip'ed datafiles via popen(),
# if popen() is available: set preprocessor variable HAVE_POPEN.
# 
AC_CHECK_FUNCS(popen);




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


# Write out easel.h header
AC_CONFIG_HEADERS(easel.h)
AC_CONFIG_FILES(miniapps/Makefile)
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
