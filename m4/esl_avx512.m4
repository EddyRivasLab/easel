# ESL_AVX512([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Checks whether compiler supports features we need in our AVX-512
# implementations.
#
# We call this "AVX-512" in the generic sense of the lab's vector
# implementations in Easel, HMMER, Infernal, etc; more precisely,
# we're checking for AVX-512 Foundation, Exponential/Reciprocal, and
# Byte/Word instruction subsets. Xeon Phi200 "Knights Landing" does
# not support the AVX512BW instructions. We await Skylake EX/EP
# processors for that.
#
# Tries to compile and link a test program using the current CC,
# CFLAGS, and (optionally) any AVX512_CFLAGS passed by the user. If
# AVX512_CFLAGS are not provided, then we try to determine them, by
# trying nothing (i.e. the compiler deals with AVX-512 intrinsics by
# default), then trying each of:
#    "-mavx512f -mavx512er -mavx512bw"
#    "-march=skylake-avx512"
#
# Sets $esl_have_avx512 = yes | no
# Sets $esl_avx512_cflags to any extra needed CFLAGS, such as "-mavx512f -mavx512er -mavx512bw"
# 
# A typical ACTION-IF-FOUND might be:
#      AC_DEFINE(HAVE_AVX512)
#      AVX512_CFLAGS=$esl_avx512_cflags
#      AC_SUBST(AVX512_CFLAGS)
#
# (Note that this is correct, even if user provided AVX512_CFLAGS;
#  because if they did, all the macro does is make sure that they
#  work, and sets esl_avx512_cflags=$AVX512_CFLAGS.)
#
AC_DEFUN([ESL_AVX512],[
  AC_MSG_CHECKING([whether $CC can compile our AVX-512 code])
  esl_have_avx512=no
 
  # need to loop over array of whitespace-containing strings;
  # here's one way to do it, that I think is portable, by 
  # temporarily munging IFS to use comma-delimited lists.
  save_IFS=$IFS
  IFS=,
  if test "x$AVX512_CFLAGS" != x; then 
    esl_avx512_try_flags=$AVX_CFLAGS
  else
    esl_avx512_try_flags="none,-mavx512f -mavx512er -mavx512bw,-march=skylake-avx512"
  fi

  for esl_avx512_cflags in $esl_avx512_try_flags; do 
    case $esl_avx512_cflags in 
      none) CFLAGS=$save_CFLAGS;;
      *)    CFLAGS="$save_CFLAGS $esl_avx512_cflags";;
    esac

    # The test program could be improved. It's only using Foundation instructions.
    IFS=$save_IFS
    AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <x86intrin.h>
#include <stdint.h>
int stub_avx512(void) {
__m512i v1 = _mm512_set1_epi32(42);
__m512i v2 = _mm512_set1_epi32(470);
union { __m512i v; int32_t x[16]; } v3;
v3.v = _mm512_add_epi32(v1, v2);
return (int) v3.x[0];
}
int main(void) { if (stub_avx512() != 512) return 1; else return 0;}
    ]])], [ esl_have_avx512=yes; break; ], [])
    IFS=,
  done

  IFS=$save_IFS
  CFLAGS=$save_CFLAGS

  if test "$esl_have_avx512" = "no"; then
    AC_MSG_RESULT([no])
    esl_avx512_cflags=
  else
    if test "$esl_avx512_cflags" = "none"; then
      AC_MSG_RESULT([yes])
      esl_avx512_cflags=
    else
      AC_MSG_RESULT([yes, with $esl_avx512_cflags])
    fi    
  fi

  AS_VAR_IF([esl_have_avx512],[yes],[$1],[$2])
])



