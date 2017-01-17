# ESL_NEON([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check whether compiler supports features we need in our ARM NEON
# implementations.
#
# Tries to compile and link a test program using the current CC,
# CFLAGS, and (optionally) any NEON_CFLAGS passed by the user. If no
# NEON_CFLAGS are provided, then we try to determine them, by trying
# nothing (i.e. the compiler deals with ARM NEON intrinsics by
# default), then trying -mfpu=neon.
# 
# Sets $esl_have_neon = yes | no
# Sets $esl_neon_cflags to any needed CFLAGS, such as -mfpu-neon
# 
# A typical ACTION-IF-FOUND might be:
#      AC_DEFINE(HAVE_NEON)
#      NEON_CFLAGS=$esl_neon_cflags
#      AC_SUBST(NEON_CFLAGS)
#
# (Note that this is correct, even if user provided NEON_CFLAGS;
#  because if they did, all the macro does is make sure that they
#  work, and sets esl_neon_cflags=$NEON_CFLAGS.)
#
AC_DEFUN([ESL_NEON], [
  AC_MSG_CHECKING([whether $CC can compile our NEON code])
  esl_have_neon=no

  if test "x$NEON_CFLAGS" != x; then 
    esl_neon_try_flags=$NEON_CFLAGS
  else
    esl_neon_try_flags="none -mfpu=neon"
  fi

  save_CFLAGS=$CFLAGS
  for esl_neon_cflags in $esl_neon_try_flags; do 
    case $esl_neon_cflags in 
      none) CFLAGS=$save_CFLAGS;;
      *)    CFLAGS="$save_CFLAGS $esl_neon_cflags";;
    esac

    AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[
#include <stdint.h>
#include <arm_neon.h>
int stub_neon(void) {
  int32x4_t v1 = vdupq_n_s32(42);
  int32x4_t v2 = vdupq_n_s32(86);
  union { int32x4_t v; int32_t x[4]; } v3;
  v3.v = vaddq_s32(v1, v2);
  return (int) v3.x[0];
}
int main(void) { if (stub_neon() != 128) return 1; else return 0; }
  ]])], [ esl_have_neon=yes; break; ], [])
  done

  CFLAGS=$save_CFLAGS
  if test "$esl_have_neon" = "no"; then
    AC_MSG_RESULT([no])
    esl_neon_cflags=
  else
    if test "$esl_neon_cflags" = "none"; then
      AC_MSG_RESULT([yes])
      esl_neon_cflags=
    else
      AC_MSG_RESULT([yes, with $esl_neon_cflags])
    fi
  fi

  AS_VAR_IF([esl_have_neon],[yes],[$1],[$2])
])



