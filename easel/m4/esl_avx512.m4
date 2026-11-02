# ESL_AVX512([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check whether the compiler supports the AVX-512F/BW intrinsics used by
# the HMMER3 impl_avx512 implementation (SSV filter).
#
# Tries to compile and link a test program using the current CC, CFLAGS,
# and (optionally) any AVX512_CFLAGS passed by the user. If none are
# provided, tries nothing first (e.g. -march=cascadelake already enables
# AVX-512), then -mavx512f -mavx512bw.
#
# Sets $esl_have_avx512 = yes | no
# Sets $esl_avx512_cflags to any needed CFLAGS.
#
AC_DEFUN([ESL_AVX512], [
  AC_MSG_CHECKING([whether $CC can compile our AVX-512 code])
  esl_have_avx512=no

  if test "x$AVX512_CFLAGS" != x; then
    esl_avx512_try_flags=$AVX512_CFLAGS
  else
    esl_avx512_try_flags="none -mavx512f -mavx512bw"
  fi

  save_CFLAGS=$CFLAGS
  for esl_avx512_cflags in $esl_avx512_try_flags; do
    case $esl_avx512_cflags in
      none) CFLAGS=$save_CFLAGS;;
      *)    CFLAGS="$save_CFLAGS $esl_avx512_cflags";;
    esac

    AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[
#include <immintrin.h>
#include <stdint.h>
int stub_avx512(void) {
__m512i v1 = _mm512_set1_epi8((char) 3);
__m512i v2 = _mm512_set1_epi8((char) 86);
__m512i v3 = _mm512_subs_epi8(v1, v2);
__m512i v4 = _mm512_max_epu8(v1, v3);
__m512i v5 = _mm512_alignr_epi8(v4, _mm512_mask_shuffle_i32x4(_mm512_setzero_si512(), 0xE, v4, v4, 0x90), 15);
__m128i r = _mm_max_epu8(_mm512_castsi512_si128(v5), _mm512_extracti32x4_epi32(v5, 3));
return (int) _mm_cvtsi128_si32(r);
}
int main(void) { if (stub_avx512() < 0) return 1; else return 0; }
  ]])], [ esl_have_avx512=yes; break; ], [])
  done

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
