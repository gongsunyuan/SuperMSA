/* PMSABIN1 sequence file backend: mmap'd, pre-digitized database.
 *
 * See esl_sqio_pmsa.h for the file layout. Only the interfaces used by
 * jackhmmer are implemented: read, read_info, read_block, position(0),
 * set_digital, is_rewindable, get_error, close. Everything else stays
 * NULL in the ESL_SQFILE method table.
 */
#include <esl_config.h>

#include <fcntl.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_sqio_pmsa.h"

#define PMSA_MAGIC      "PMSABIN1"
#define PMSA_HDR_SIZE   64
#define PMSA_IDX_SIZE   16

static int  pmsa_Position    (ESL_SQFILE *sqfp, off_t offset);
static void pmsa_Close       (ESL_SQFILE *sqfp);
static int  pmsa_SetDigital  (ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
static int  pmsa_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);
static int  pmsa_Read        (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  pmsa_ReadInfo    (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  pmsa_ReadBlock   (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target);
static int  pmsa_IsRewindable(const ESL_SQFILE *sqfp);
static const char *pmsa_GetError(const ESL_SQFILE *sqfp);

/* pmsa_load():
 * Shared by Read/ReadInfo. On success advances <next> and fills name,
 * acc, desc, L, idx and record offsets; with <with_seq> also digitizes
 * (memcpy's) the residues into <sq>.
 */
static int
pmsa_load(ESL_SQFILE *sqfp, ESL_SQ *sq, int with_seq)
{
  ESL_SQPMSA_DATA *p   = &sqfp->data.pmsa;
  const uint8_t   *ient;
  const uint8_t   *rec;
  const char      *name;
  const char      *acc;
  const char      *desc;
  uint64_t         rec_off;
  uint32_t         len;
  int              status;

  if (p->next >= p->nseq) return eslEOF;
  if (! sqfp->do_digital || ! esl_sq_IsDigital(sq)) return eslEINVAL;

  ient = p->mem + p->index_off + p->next * PMSA_IDX_SIZE;
  memcpy(&rec_off, ient,     sizeof(uint64_t));
  memcpy(&len,     ient + 8, sizeof(uint32_t));

  rec  = p->mem + rec_off;
  name = (const char *) rec + len + 2;
  acc  = name + strlen(name) + 1;
  desc = acc  + strlen(acc)  + 1;

  if ((status = esl_sq_SetName(sq, name))      != eslOK) return status;
  if ((status = esl_sq_SetAccession(sq, acc))  != eslOK) return status;
  if ((status = esl_sq_SetDesc(sq, desc))      != eslOK) return status;

  if (with_seq)
    {
      if ((status = esl_sq_GrowTo(sq, (int64_t) len)) != eslOK) return status;
      memcpy(sq->dsq, rec, ((size_t) len + 2) * sizeof(ESL_DSQ));
      sq->n = (int64_t) len;
    }
  else
    sq->n = 0;

  sq->L     = (int64_t) len;
  sq->start = 1;
  sq->end   = (int64_t) len;
  sq->C     = 0;
  sq->W     = (with_seq ? (int64_t) len : 0);
  sq->idx   = (int64_t) p->next;
  sq->roff  = (off_t) rec_off;
  /* eoff is only consumed by the MPI block planner, which this format
   * does not support; skip the three strlen() calls its computation
   * would otherwise need on every sequence. */
  sq->eoff  = -1;
  sq->hoff  = sq->doff = -1;

  p->next++;
  return eslOK;
}

static int
pmsa_Read(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  return pmsa_load(sqfp, sq, TRUE);
}

static int
pmsa_ReadInfo(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  return pmsa_load(sqfp, sq, FALSE);
}

/* pmsa_ReadBlock():
 * Mirrors the !long_target branch of sqascii_ReadBlock(): whole
 * sequences, up to <max_sequences> and MAX_RESIDUE_COUNT residues.
 * Long-target (windowed DNA) reading is not supported for this format.
 */
static int
pmsa_ReadBlock(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int max_init_window, int long_target)
{
  int i      = 0;
  int size   = 0;
  int status = eslOK;

  if (long_target) return eslEINVAL;

  sqBlock->count = 0;
  if (max_sequences < 1 || max_sequences > sqBlock->listSize)
    max_sequences = sqBlock->listSize;

  for (i = 0; i < max_sequences && size < MAX_RESIDUE_COUNT; ++i)
    {
      status = pmsa_Read(sqfp, sqBlock->list + i);
      if (status != eslOK) break;
      size += sqBlock->list[i].n;
      ++sqBlock->count;
    }

  /* EOF will be returned only in the case where no sequences were read */
  if (status == eslEOF && i > 0) status = eslOK;

  sqBlock->complete = TRUE;
  return status;
}

static int
pmsa_Position(ESL_SQFILE *sqfp, off_t offset)
{
  ESL_SQPMSA_DATA *p = &sqfp->data.pmsa;

  if (offset != 0) return eslEINVAL;  /* only rewind is supported */
  p->next = 0;
  return eslOK;
}

static void
pmsa_Close(ESL_SQFILE *sqfp)
{
  ESL_SQPMSA_DATA *p = &sqfp->data.pmsa;

  if (p->mem != NULL) munmap(p->mem, (size_t) p->memsize);
  if (p->fd  >= 0)    close(p->fd);
  p->mem     = NULL;
  p->memsize = 0;
  p->fd      = -1;
}

static int
pmsa_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
{
  /* Sequences are already digital; esl_sqfile_SetDigital() sets
   * do_digital and abc itself. No inmap is needed. */
  return (abc == NULL ? eslEINVAL : eslOK);
}

static int
pmsa_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)
{
  *ret_type = eslAMINO;   /* dbprep digitizes with eslAMINO */
  return eslOK;
}

static int
pmsa_IsRewindable(const ESL_SQFILE *sqfp)
{
  return TRUE;
}

static const char *
pmsa_GetError(const ESL_SQFILE *sqfp)
{
  return sqfp->data.pmsa.errbuf;
}

/* Function:  esl_sqpmsa_Open()
 *
 * Called from sqfile_open() before the ascii/ncbi backends, whenever
 * format is eslSQFILE_UNKNOWN or eslSQFILE_PMSA. If the file does not
 * carry the PMSABIN1 magic, returns eslENOTFOUND so that the remaining
 * backends get their turn (unless the caller explicitly asserted
 * eslSQFILE_PMSA, in which case eslEFORMAT).
 */
int
esl_sqpmsa_Open(const char *filename, int format, ESL_SQFILE *sqfp)
{
  ESL_SQPMSA_DATA *p = &sqfp->data.pmsa;
  uint8_t          hdr[PMSA_HDR_SIZE];
  uint64_t         nseq;
  uint64_t         index_off;
  struct stat      st;
  ssize_t          nread;
  int              fd;

  if (format != eslSQFILE_UNKNOWN && format != eslSQFILE_PMSA) return eslENOTFOUND;

  p->fd        = -1;
  p->mem       = NULL;
  p->memsize   = 0;
  p->nseq      = 0;
  p->next      = 0;
  p->index_off = 0;
  p->errbuf[0] = '\0';

  if ((fd = open(filename, O_RDONLY)) < 0) return eslENOTFOUND;

  nread = read(fd, hdr, PMSA_HDR_SIZE);
  if (nread < PMSA_HDR_SIZE || memcmp(hdr, PMSA_MAGIC, 8) != 0)
    {
      close(fd);
      return (format == eslSQFILE_PMSA ? eslEFORMAT : eslENOTFOUND);
    }

  if (fstat(fd, &st) != 0) { close(fd); return eslENOTFOUND; }

  memcpy(&nseq,      hdr + 8,  sizeof(uint64_t));
  memcpy(&index_off, hdr + 16, sizeof(uint64_t));

  if (index_off > (uint64_t) st.st_size ||
      nseq > ((uint64_t) st.st_size - index_off) / PMSA_IDX_SIZE)
    {
      close(fd);
      snprintf(p->errbuf, sizeof(p->errbuf), "corrupt PMSA header in %s", filename);
      return eslEFORMAT;
    }

  p->mem = mmap(NULL, (size_t) st.st_size, PROT_READ, MAP_SHARED, fd, 0);
  if (p->mem == MAP_FAILED)
    {
      close(fd);
      p->mem = NULL;
      snprintf(p->errbuf, sizeof(p->errbuf), "mmap failed for %s", filename);
      return eslFAIL;
    }

  p->fd        = fd;
  p->memsize   = (uint64_t) st.st_size;
  p->nseq      = nseq;
  p->index_off = index_off;

  sqfp->format          = eslSQFILE_PMSA;
  sqfp->position        = &pmsa_Position;
  sqfp->close           = &pmsa_Close;
  sqfp->set_digital     = &pmsa_SetDigital;
  sqfp->guess_alphabet  = &pmsa_GuessAlphabet;
  sqfp->is_rewindable   = &pmsa_IsRewindable;
  sqfp->read            = &pmsa_Read;
  sqfp->read_info       = &pmsa_ReadInfo;
  sqfp->read_block      = &pmsa_ReadBlock;
  sqfp->get_error       = &pmsa_GetError;

  return eslOK;
}
