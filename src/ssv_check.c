/* ssv_check: cross-validate AVX-512 p7_SSVFilter against SSE p7_MSVFilter.
 *
 * Invariant: when SSV returns eslOK (J state unused), its score must equal
 * the MSV score exactly. Both live in the same impl_avx512 build, so this
 * isolates the AVX-512 SSV port (layout, shift, arithmetic) per model size.
 *
 * Build (in src/): gcc -O3 -march=cascadelake -I. -I../easel -I./impl_avx512 \
 *   -o ssv_check ssv_check.c libhmmer.a ../easel/libeasel.a -lpthread -lm
 */
#include <p7_config.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"
#include "impl_avx512.h"

int
main(int argc, char **argv)
{
  int          seed = (argc > 1 ? atoi(argv[1]) : 42);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(seed);
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg  = p7_bg_Create(abc);
  int  Ms[] = { 1, 2, 3, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 116, 127, 128, 129,
                150, 200, 255, 256, 257, 300, 400, 500, 700, 1000, 1500 };
  int  Ls[] = { 1, 2, 50, 200, 1000 };
  int  nseq = 40;
  int  fails = 0, compared = 0, fallback = 0;
  int  mi, li, n;

  p7_bg_SetLength(bg, 200);

  for (mi = 0; mi < (int)(sizeof(Ms)/sizeof(int)); mi++)
    {
      int M = Ms[mi];
      P7_HMM      *hmm = NULL;
      P7_PROFILE  *gm  = NULL;
      P7_OPROFILE *om  = NULL;
      P7_OMX      *ox  = NULL;
      int          mfails = 0;

      if (p7_oprofile_Sample(r, abc, bg, M, 200, &hmm, &gm, &om) != eslOK)
	esl_fatal("sample failed M=%d", M);
      ox = p7_omx_Create(M, 0, 0);

      for (li = 0; li < (int)(sizeof(Ls)/sizeof(int)); li++)
	{
	  int L = Ls[li];
	  ESL_DSQ *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
	  p7_oprofile_ReconfigLength(om, L);

	  for (n = 0; n < nseq; n++)
	    {
	      float sc_ssv = 0.0f, sc_msv = 0.0f;
	      int   st1, st2;

	      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
	      st1 = p7_SSVFilter(dsq, L, om, &sc_ssv);
	      st2 = p7_MSVFilter(dsq, L, om, ox, &sc_msv);

	      if (st1 == eslENORESULT) { fallback++; continue; }
	      if (st1 == eslERANGE || st2 == eslERANGE) continue;

	      compared++;
	      if (fabsf(sc_ssv - sc_msv) > 1e-4)
		{
		  if (mfails < 3)
		    fprintf(stderr, "MISMATCH M=%d L=%d n=%d ssv=%.6f msv=%.6f\n",
			    M, L, n, sc_ssv, sc_msv);
		  mfails++;
		  fails++;
		}
	    }
	  free(dsq);
	}
      if (mfails) fprintf(stderr, "  --> M=%d: %d mismatches\n", M, mfails);

      p7_omx_Destroy(ox);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
    }

  printf("compared=%d mismatches=%d (ssv fallback to MSV: %d)\n", compared, fails, fallback);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return (fails ? 1 : 0);
}
