# ESL_ALTIVEC([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check whether compiler supports features we need in our Altivec/VMX
# implementations.
#
# Tries to compile and link a test program using the current CC,
# CFLAGS, and (optionally) any ALTIVEC_CFLAGS passed by the user. If no
# ALTIVEC_CFLAGS are provided, then we try to determine them, by trying
# nothing, then trying ...
# 
# Sets $esl_have_altivec = yes | no
# Sets $esl_altivec_cflags to any needed CFLAGS, such as ...
# 
# A typical ACTION-IF-FOUND might be:
#      AC_DEFINE(HAVE_ALTIVEC)
#      ALTIVEC_CFLAGS=$esl_altivec_cflags
#      AC_SUBST(ALTIVEC_CFLAGS)
#
AC_DEFUN([ESL_ALTIVEC], [
  AC_MSG_CHECKING([whether $CC can compile our Altivec/VMX code])
  esl_have_altivec=no

  if test "x$ALTIVEC_CFLAGS" != x; then 
    esl_altivec_try_flags=$ALTIVEC_CFLAGS
  else
    esl_altivec_try_flags="none -maltivec -faltivec"
  fi

  save_CFLAGS=$CFLAGS
  for esl_altivec_cflags in $esl_altivec_try_flags; do 
    case $esl_altivec_cflags in 
      none) CFLAGS=$save_CFLAGS;;
      *)    CFLAGS="$save_CFLAGS $esl_altivec_cflags";;
    esac

    AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[
#include <altivec.h>
#include <stdint.h>
int stub_altivec(void) {
vector signed int v1,v2;   
union { vector signed int v; int32_t x[4]; } v3;
v1   = vec_splat_s32(6);
v2   = vec_splat_s32(12);
v3.v = vec_add(v1,v2);
return (int) v3.x[0];}
int main(void) { if (stub_altivec() != 18) return 1; else return 0; }
  ]])], [ esl_have_altivec=yes; break; ], [])
  done
  CFLAGS=$save_CFLAGS

  if test "$esl_have_altivec" = "no"; then
    AC_MSG_RESULT([no])
    esl_altivec_cflags=
  else
    if test "$esl_altivec_cflags" = "none"; then
      AC_MSG_RESULT([yes])
      esl_altivec_cflags=
    else
      AC_MSG_RESULT([yes, with $esl_altivec_cflags])
    fi
  fi

  AS_VAR_IF([esl_have_altivec],[yes],[$1],[$2])
])



