/* Unaligned sequence i/o: PMSA binary format (mmap'd, pre-digitized).
 *
 * PMSABIN1 is a SuperMSA competition format: a FASTA database converted
 * once by src/dbprep.c into digital (ESL_DSQ) records, so that repeated
 * database scans skip ASCII parsing entirely. The whole file is mmap'd
 * read-only; concurrent jackhmmer processes share one page cache.
 *
 * Layout (all integers native little-endian, x86-64):
 *   Header, 64 bytes:
 *     [0..7]   magic "PMSABIN1"
 *     [8]      u64 nseq
 *     [16]     u64 index_off   (byte offset of index block, at file tail)
 *     [24]     u64 records_off (= 64, informational)
 *     [32]     u64 file_size
 *     [40]     u64 total_res   (informational)
 *     [48]     u64 max_len     (informational)
 *     [56]     u64 reserved (0)
 *   Records, from offset 64, variable length:
 *     len+2 bytes of dsq (including both sentinel bytes), then
 *     name\0 acc\0 desc\0 (acc/desc may be empty strings)
 *   Index, nseq entries of 16 bytes each at index_off:
 *     +0  u64 rec_off (absolute file offset of record)
 *     +8  u32 len     (number of residues)
 *     +12 u32 reserved (0)
 *
 * Digital-mode only; the amino acid alphabet mapping is baked in by
 * dbprep using easel's own FASTA reader, so sequences are byte-identical
 * to what esl_sqascii would produce.
 */
#ifndef eslSQIO_PMSA_INCLUDED
#define eslSQIO_PMSA_INCLUDED
#include <esl_config.h>

#include <stdio.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <stdint.h>

#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

/* forward declaration */
struct esl_sqio_s;

typedef struct esl_sqpmsa_s {
  int      fd;         /* underlying fd for the mmap                    */
  uint8_t *mem;        /* mmap base                                     */
  uint64_t memsize;    /* file size in bytes                            */
  uint64_t nseq;       /* total number of sequence records              */
  uint64_t next;       /* index of next record to read, 0..nseq         */
  uint64_t index_off;  /* byte offset of the index block                */

  char     errbuf[eslERRBUFSIZE]; /* error message; size matches ascii   */
} ESL_SQPMSA_DATA;

extern int esl_sqpmsa_Open(const char *seqfile, int format, struct esl_sqio_s *sqfp);
extern int esl_sqpmsa_ReadLazy(struct esl_sqio_s *sqfp, ESL_SQ *sq);
extern int esl_sqpmsa_ReadBlockLazy(struct esl_sqio_s *sqfp, ESL_SQ_BLOCK *sqBlock);

#endif /*eslSQIO_PMSA_INCLUDED*/
