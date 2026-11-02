/* ssv_check2: cross-validate AVX-512 SSV vs SSE MSV on a REAL query profile.
 *
 * Builds the round-1 profile exactly as jackhmmer does (p7_SingleBuilder),
 * then compares scores for every sequence of a real database (demodb).
 *
 * Build (in src/): gcc -O3 -march=cascadelake -I. -I../easel -I./impl_avx512 \
 *   -o ssv_check2 ssv_check2.c libhmmer.a ../easel/libeasel.a -lpthread -lm
 * Run: ./ssv_check2 ../data/query.fasta ../data/demodb.fasta
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
#include "impl_avx512.h"

int
main(int argc, char **argv)
{
  const char *qfile, *dbfile;
  ESL_ALPHABET *abc = esl_alphabet_Create(eslAMINO);
  ESL_SQFILE   *qfp = NULL, *dbfp = NULL;
  ESL_SQ       *qsq = NULL, *sq = NULL;
  P7_BG        *bg  = NULL;
  P7_BUILDER   *bld = NULL;
  P7_OPROFILE  *om  = NULL;
  P7_OMX       *ox  = NULL;
  P7_TRACE     *qtr = NULL;
  int           compared = 0, fails = 0, fallback = 0, nseq = 0;
  int           status;

  if (argc != 3) esl_fatal("usage: ssv_check2 <query.fasta> <db.fasta>");
  qfile = argv[1]; dbfile = argv[2];

  bg = p7_bg_Create(abc);
  bld = p7_builder_Create(NULL, abc);
  if (bld == NULL) esl_fatal("builder create failed");
  if (p7_builder_LoadScoreSystem(bld, "BLOSUM62", 0.45f, 0.02f, bg) != eslOK)
    esl_fatal("load score system failed");

  if (esl_sqfile_OpenDigital(abc, qfile, eslSQFILE_UNKNOWN, NULL, &qfp) != eslOK)
    esl_fatal("cannot open query %s", qfile);
  qsq = esl_sq_CreateDigital(abc);
  if (esl_sqio_Read(qfp, qsq) != eslOK) esl_fatal("cannot read query");
  fprintf(stderr, "query %s L=%ld\n", qsq->name, (long) qsq->n);

  if (p7_SingleBuilder(bld, qsq, bg, NULL, &qtr, NULL, &om) != eslOK)
    esl_fatal("single builder failed");
  fprintf(stderr, "profile M=%d bias_b=%d base_b=%d tjb=%d tbm=%d tec=%d\n",
	  om->M, om->bias_b, om->base_b, om->tjb_b, om->tbm_b, om->tec_b);
  ox = p7_omx_Create(om->M, 0, 0);

  if (esl_sqfile_OpenDigital(abc, dbfile, eslSQFILE_UNKNOWN, NULL, &dbfp) != eslOK)
    esl_fatal("cannot open db %s", dbfile);
  sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(dbfp, sq)) == eslOK)
    {
      float sc_ssv = 0.0f, sc_msv = 0.0f;
      int   st1, st2;

      p7_oprofile_ReconfigLength(om, sq->n);
      st1 = p7_SSVFilter(sq->dsq, sq->n, om, &sc_ssv);
      st2 = p7_MSVFilter(sq->dsq, sq->n, om, ox, &sc_msv);
      nseq++;

      if (st1 == eslENORESULT) { fallback++; esl_sq_Reuse(sq); continue; }
      if (st1 == eslERANGE || st2 == eslERANGE) { esl_sq_Reuse(sq); continue; }

      compared++;
      if (fabsf(sc_ssv - sc_msv) > 1e-4)
	{
	  if (fails < 10)
	    fprintf(stderr, "MISMATCH %s L=%ld ssv=%.6f msv=%.6f\n",
		    sq->name, (long) sq->n, sc_ssv, sc_msv);
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
  p7_omx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_trace_Destroy(qtr);
  p7_builder_Destroy(bld);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  return (fails ? 1 : 0);
}
