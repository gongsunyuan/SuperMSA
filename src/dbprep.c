/* dbprep: convert a (FASTA) sequence database to PMSABIN1 binary format.
 *
 * Usage: dbprep <in.fasta> <out.bin>
 *
 * Digitization is done by easel's own reader with the eslAMINO alphabet,
 * so the dsq bytes stored here are exactly what jackhmmer would get from
 * parsing the FASTA file itself. Layout is documented in
 * easel/esl_sqio_pmsa.h: 64-byte header (backfilled at the end), then
 * variable-length records, then the 16-byte-per-sequence index at the
 * file tail. One pass over the input; the index is held in memory
 * (16 bytes/sequence) and written last.
 */
#include <esl_config.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#define HDR_SIZE 64
#define IDX_SIZE 16

static void put_u64(uint8_t *p, uint64_t v) { memcpy(p, &v, sizeof(uint64_t)); }
static void put_u32(uint8_t *p, uint32_t v) { memcpy(p, &v, sizeof(uint32_t)); }

static void
write_bytes(FILE *out, const void *buf, size_t n, const char *outfn)
{
  if (n > 0 && fwrite(buf, 1, n, out) != n)
    esl_fatal("write failed on %s (disk full?)", outfn);
}

int
main(int argc, char **argv)
{
  const char      *infn;
  const char      *outfn;
  ESL_ALPHABET    *abc;
  ESL_SQFILE      *sqfp;
  ESL_SQ          *sq;
  FILE            *out;
  uint8_t         *idx       = NULL;
  uint64_t         idx_cap   = 1 << 20;
  uint64_t         nseq      = 0;
  uint64_t         total_res = 0;
  uint64_t         max_len   = 0;
  uint64_t         index_off;
  uint64_t         file_size;
  uint8_t          hdr[HDR_SIZE];
  time_t           t0        = time(NULL);
  int              status;

  if (argc != 3)
    esl_fatal("usage: dbprep <in.fasta> <out.bin>");
  infn  = argv[1];
  outfn = argv[2];

  if (sizeof(ESL_DSQ) != 1) esl_fatal("PMSABIN1 requires 1-byte ESL_DSQ");

  abc = esl_alphabet_Create(eslAMINO);
  if (abc == NULL) esl_fatal("failed to create amino alphabet");

  status = esl_sqfile_OpenDigital(abc, infn, eslSQFILE_UNKNOWN, NULL, &sqfp);
  if (status != eslOK) esl_fatal("failed to open %s (status %d)", infn, status);

  sq = esl_sq_CreateDigital(abc);
  if (sq == NULL) esl_fatal("failed to allocate sequence");

  if ((out = fopen(outfn, "wb")) == NULL) esl_fatal("failed to create %s", outfn);

  if ((idx = malloc((size_t) idx_cap * IDX_SIZE)) == NULL) esl_fatal("out of memory");

  /* reserve the header; it is backfilled once counts are known */
  memset(hdr, 0, HDR_SIZE);
  write_bytes(out, hdr, HDR_SIZE, outfn);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      off_t rec_off = ftello(out);
      int   nl      = (int) strlen(sq->name);
      int   al      = (int) strlen(sq->acc);
      int   dl      = (int) strlen(sq->desc);

      if (sq->n > UINT32_MAX) esl_fatal("sequence %s too long for PMSABIN1", sq->name);

      write_bytes(out, sq->dsq,  (size_t) sq->n + 2, outfn);
      write_bytes(out, sq->name, (size_t) nl + 1,    outfn);
      write_bytes(out, sq->acc,  (size_t) al + 1,    outfn);
      write_bytes(out, sq->desc, (size_t) dl + 1,    outfn);

      if (nseq == idx_cap)
	{
	  uint8_t *tmp;
	  idx_cap *= 2;
	  if ((tmp = realloc(idx, (size_t) idx_cap * IDX_SIZE)) == NULL) esl_fatal("out of memory");
	  idx = tmp;
	}
      put_u64(idx + nseq * IDX_SIZE,     (uint64_t) rec_off);
      put_u32(idx + nseq * IDX_SIZE + 8, (uint32_t) sq->n);
      put_u32(idx + nseq * IDX_SIZE + 12, 0);

      nseq++;
      total_res += (uint64_t) sq->n;
      if ((uint64_t) sq->n > max_len) max_len = (uint64_t) sq->n;

      if (nseq % 1000000 == 0)
	fprintf(stderr, "dbprep: %llu sequences, %llu residues so far...\n",
		(unsigned long long) nseq, (unsigned long long) total_res);

      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal("parse error reading %s: %s", infn, esl_sqfile_GetErrorBuf(sqfp));

  /* index goes at the tail, then backfill the header */
  index_off = (uint64_t) ftello(out);
  write_bytes(out, idx, (size_t)(nseq * IDX_SIZE), outfn);
  file_size = (uint64_t) ftello(out);

  memset(hdr, 0, HDR_SIZE);
  memcpy(hdr, "PMSABIN1", 8);
  put_u64(hdr + 8,  nseq);
  put_u64(hdr + 16, index_off);
  put_u64(hdr + 24, HDR_SIZE);      /* records_off */
  put_u64(hdr + 32, file_size);
  put_u64(hdr + 40, total_res);
  put_u64(hdr + 48, max_len);
  if (fseeko(out, 0, SEEK_SET) != 0) esl_fatal("failed to rewind %s", outfn);
  write_bytes(out, hdr, HDR_SIZE, outfn);

  if (fclose(out) != 0) esl_fatal("close failed on %s", outfn);

  printf("dbprep: %s -> %s\n", infn, outfn);
  printf("  sequences: %llu\n", (unsigned long long) nseq);
  printf("  residues:  %llu\n", (unsigned long long) total_res);
  printf("  max_len:   %llu\n", (unsigned long long) max_len);
  printf("  file size: %.2f GB\n", (double) file_size / (1024.0 * 1024.0 * 1024.0));
  printf("  elapsed:   %ld s\n", (long)(time(NULL) - t0));

  free(idx);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  return 0;
}
