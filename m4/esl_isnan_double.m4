# SYNOPSIS
#
#   ESL_ISNAN_DOUBLE
#
# DESCRIPTION
#
#   Check whether isnan works for double NaNs. This function will define
#   ISNAN_DOUBLE_WORKS when this is indeed the case.

AC_DEFUN([ESL_ISNAN_DOUBLE], [
  AC_MSG_CHECKING([whether isnan works for double NaNs])

  AC_RUN_IFELSE([AC_LANG_SOURCE([[
    #include <math.h>
    #ifndef NAN
    #define NAN (0./0.)
    #endif
    int main (void)
    {
      int status;
      double nan = NAN;
      status = !isnan(nan);
      exit (status);
    }
  ]])],
    [c_isnan_double_works="yes"; AC_MSG_RESULT(yes)],
    [c_isnan_double_works="no"; AC_MSG_RESULT(no)])

  if test "$c_isnan_double_works" != no ; then
    AC_DEFINE(ISNAN_DOUBLE_WORKS,1,[Define this if isnan works for double NaNs])
  fi
])
