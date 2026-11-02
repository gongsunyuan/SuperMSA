/* ssv_check3: validate AVX-512 p7_SSVFilter against GENERIC scalar DP.
 *
 * check2 was invalid: p7_MSVFilter() itself calls p7_SSVFilter first, so
 * it compared AVX-SSV against itself. Here the reference is p7_GViterbi()
 * on a generic profile adjusted by p7_profile_SameAsMF() (same rounding,
 * no gaps, same NCJ correction) -- the approach of msvfilter.c's utest.
 *
 * When SSV returns eslOK, sc1 must equal sc2 (within 1e-3).
 *
 * Build (in src/): gcc -O3 -march=cascadelake -I. -I../easel -I./impl_avx512 \
 *   -o ssv_check3 ssv_check3.c libhmmer.a ../easel/libeasel.a -lpthread -lm
 * Run: ./ssv_check3 <query.fasta> <db.fasta> [max_seqs]
 */
#include <p7_config.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

int
main(int argc, char **argv)
{
  const char *qfile, *dbfile;
  ESL_ALPHABET *abc = esl_alphabet_Create(eslAMINO);
  ESL_SQFILE   *qfp = NULL, *dbfp = NULL;
  ESL_SQ       *qsq = NULL, *sq = NULL;
  P7_BG        *bg  = NULL;
  P7_BUILDER   *bld = NULL;
  P7_HMM       *hmm = NULL;
  P7_PROFILE   *gm  = NULL;
  P7_OPROFILE  *om  = NULL;
  P7_GMX       *gx  = NULL;
  P7_TRACE     *qtr = NULL;
  int           compared = 0, fails = 0, fallback = 0, nseq = 0;
  int           max_seqs = 0;
  int           status;

  if (argc < 3) esl_fatal("usage: ssv_check3 <query.fasta> <db.fasta> [max_seqs]");
  qfile = argv[1]; dbfile = argv[2];
  if (argc > 3) max_seqs = atoi(argv[3]);

  bg  = p7_bg_Create(abc);
  bld = p7_builder_Create(NULL, abc);
  if (p7_builder_LoadScoreSystem(bld, "BLOSUM62", 0.45f, 0.02f, bg) != eslOK)
    esl_fatal("load score system failed");

  if (esl_sqfile_OpenDigital(abc, qfile, eslSQFILE_UNKNOWN, NULL, &qfp) != eslOK)
    esl_fatal("cannot open query %s", qfile);
  qsq = esl_sq_CreateDigital(abc);
  if (esl_sqio_Read(qfp, qsq) != eslOK) esl_fatal("cannot read query");

  /* build the HMM, then derive gm/om exactly like msvfilter_example */
  if (p7_SingleBuilder(bld, qsq, bg, &hmm, &qtr, NULL, NULL) != eslOK)
    esl_fatal("single builder failed");
  fprintf(stderr, "query %s L=%ld\n", qsq->name, (long) qsq->n);

  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, qsq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_profile_SameAsMF(om, gm);
  fprintf(stderr, "profile M=%d bias_b=%d base_b=%d\n", om->M, om->bias_b, om->base_b);

  if (esl_sqfile_OpenDigital(abc, dbfile, eslSQFILE_UNKNOWN, NULL, &dbfp) != eslOK)
    esl_fatal("cannot open db %s", dbfile);
  sq = esl_sq_CreateDigital(abc);
  gx = p7_gmx_Create(gm->M, 200);

  while ((status = esl_sqio_Read(dbfp, sq)) == eslOK)
    {
      float sc1 = 0.0f, sc2 = 0.0f;
      int   st;

      nseq++;
      if (max_seqs && nseq > max_seqs) break;

      /* No per-sequence ReconfigLength on om/gm: keep the initial L0
       * configuration (same pattern as msvfilter.c's own utest).
       * SameAsMF is not idempotent and must run only once. */
      p7_gmx_GrowTo(gx, gm->M, sq->n);

      st = p7_SSVFilter(sq->dsq, sq->n, om, &sc1);
      if (st == eslENORESULT) { fallback++; esl_sq_Reuse(sq); continue; }
      if (st == eslERANGE)    { esl_sq_Reuse(sq); continue; }

      p7_GViterbi(sq->dsq, sq->n, gm, gx, &sc2);
      sc2 = sc2 / om->scale_b - 3.0f;

      compared++;
      if (fabsf(sc1 - sc2) > 1e-3)
	{
	  if (fails < 10)
	    fprintf(stderr, "MISMATCH %s L=%ld ssv_avx=%.6f generic=%.6f\n",
		    sq->name, (long) sq->n, sc1, sc2);
	  fails++;
	}
      esl_sq_Reuse(sq);
    }

  printf("query=%s M=%d nseq=%d compared=%d mismatches=%d fallback=%d\n",
	 qsq->name, om->M, nseq, compared, fails, fallback);

  esl_sq_Destroy(qsq);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(qfp);
  esl_sqfile_Close(dbfp);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_trace_Destroy(qtr);
  p7_builder_Destroy(bld);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  return (fails ? 1 : 0);
}
