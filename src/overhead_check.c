/* overhead_check: bitwise validation of length-dependent fixed-cost caches.
 *
 * Build (in src/):
 *   gcc -O3 -march=cascadelake -I. -I../easel -I./impl_sse \
 *     -o overhead_check overhead_check.c libhmmer.a ../easel/libeasel.a -lpthread -lm
 * Run: ./overhead_check
 */
#include <p7_config.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "hmmer.h"

int
main(void)
{
  ESL_ALPHABET *abc = esl_alphabet_Create(eslAMINO);
  P7_BG        *bg  = p7_bg_Create(abc);
  int           L;
  int           checks = 0;
  int           fails  = 0;

  for (L = 0; L <= 20010; L++)
    {
      float actual;
      float expected;
      float p1;

      p7_bg_SetLength(bg, L);
      p1 = (float) L / (float) (L+1);
      expected = (float) L * log(p1) + log(1.-p1);
      p7_bg_NullOne(bg, NULL, L, &actual);
      checks++;
      if (memcmp(&actual, &expected, sizeof(float)) != 0)
        {
          if (fails < 10) fprintf(stderr, "standard mismatch: L=%d actual=%a expected=%a\n", L, actual, expected);
          fails++;
        }
    }

  for (L = 1; L <= 20001; L *= 10)
    {
      float actual;
      float expected;

      bg->p1 = 0.75f;
      expected = (float) L * log(bg->p1) + log(1.-bg->p1);
      p7_bg_NullOne(bg, NULL, L, &actual);
      checks++;
      if (memcmp(&actual, &expected, sizeof(float)) != 0)
        {
          fprintf(stderr, "nonstandard mismatch: L=%d actual=%a expected=%a\n", L, actual, expected);
          fails++;
        }
    }

  printf("NullOne bitwise checks=%d mismatches=%d\n", checks, fails);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  return fails ? 1 : 0;
}
