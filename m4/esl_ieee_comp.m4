# SYNOPSIS
#
#   ESL_IEEE_COMP
#
# DESCRIPTION
#
#   Let <x> be a floating point variable. IEEE 754-1985 establishes that "x == x" is
#   always false when x represents a NaN (not-a-number). This function will define
#   HAVE_IEEE_COMPARISONS when this is indeed the case.
#
# ACKNOWLEDGMENT
#
#   Based on the source code of the GSL - GNU Scientific Library.

AC_DEFUN([ESL_IEEE_COMP], [
  AC_MSG_CHECKING([whether it complies with IEEE comparisons])

  AC_RUN_IFELSE([AC_LANG_SOURCE([[
    #include <math.h>
    #ifndef NAN
    #define NAN (0./0.)
    #endif
    int main (void)
    {
      int status;
      float nanf = NAN;
      double nand = NAN;
      status = (nanf == nanf) | (nand == nand);
      exit (status);
    }
  ]])],
    [c_ieee_comparisons="yes"; AC_MSG_RESULT(yes)],
    [c_ieee_comparisons="no"; AC_MSG_RESULT(no)])

  if test "$c_ieee_comparisons" != no ; then
    AC_DEFINE(HAVE_IEEE_COMPARISONS,1,[Define this if IEEE comparisons work correctly (e.g. NaN != NaN)])
  fi
])
