# ESL_AVX2([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check whether the compiler supports the AVX2 intrinsics used by the
# HMMER3 impl_avx2 implementation (SSV filter).
#
# Tries nothing first (e.g. -march=cascadelake already enables AVX2),
# then -mavx2. Sets $esl_have_avx512 = no by symmetry is NOT used here;
# sets $esl_have_avx2 = yes | no and $esl_avx2_cflags.
#
AC_DEFUN([ESL_AVX2], [
  AC_MSG_CHECKING([whether $CC can compile our AVX2 code])
  esl_have_avx2=no

  if test "x$AVX2_CFLAGS" != x; then
    esl_avx2_try_flags=$AVX2_CFLAGS
  else
    esl_avx2_try_flags="none -mavx2"
  fi

  save_CFLAGS=$CFLAGS
  for esl_avx2_cflags in $esl_avx2_try_flags; do
    case $esl_avx2_cflags in
      none) CFLAGS=$save_CFLAGS;;
      *)    CFLAGS="$save_CFLAGS $esl_avx2_cflags";;
    esac

    AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[
#include <immintrin.h>
#include <stdint.h>
int stub_avx2(void) {
__m256i v1 = _mm256_set1_epi8((char) 3);
__m256i v2 = _mm256_set1_epi8((char) 86);
__m256i v3 = _mm256_subs_epi8(v1, v2);
__m256i v4 = _mm256_max_epu8(v1, v3);
__m256i v5 = _mm256_alignr_epi8(v4, _mm256_permute2x128_si256(v4, v4, 0x08), 15);
__m128i r = _mm_max_epu8(_mm256_castsi256_si128(v5), _mm256_extracti128_si256(v5, 1));
return (int) _mm_cvtsi128_si32(r);
}
int main(void) { if (stub_avx2() < 0) return 1; else return 0; }
  ]])], [ esl_have_avx2=yes; break; ], [])
  done

  CFLAGS=$save_CFLAGS
  if test "$esl_have_avx2" = "no"; then
    AC_MSG_RESULT([no])
    esl_avx2_cflags=
  else
    if test "$esl_avx2_cflags" = "none"; then
      AC_MSG_RESULT([yes])
      esl_avx2_cflags=
    else
      AC_MSG_RESULT([yes, with $esl_avx2_cflags])
    fi
  fi

  AS_VAR_IF([esl_have_avx2],[yes],[$1],[$2])
])
