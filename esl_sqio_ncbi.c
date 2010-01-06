/* Unaligned ncbi sequence file i/o.
 * 
 * Contents:
 *    1. An <ESL_SQFILE> object, in text mode.
 *    2. An <ESL_SQFILE> object, in digital mode. [with <alphabet>]
 *    3. Miscellaneous routines.
 *    4. Sequence reading (sequential).
 *    5. Parsing routines
 *    6. Copyright and license.
 * 
 * MSF, Thu Feb 17 17:45:51 2005
 * SVN $Id: esl_sqio.c 463 2009-11-30 19:55:01Z eddys $
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"	/* alphabet aug adds digital sequences */
#endif 
#include "esl_sqio.h"
#include "esl_sq.h"

#ifndef HAVE_ENDIAN_H
#ifndef WORDS_BIGENDIAN
#define htobe32(x) (x)
#else
#define htobe32(x) \
     ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) |		      \
      (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24))
#endif
#endif

/* format specific routines */
static int   sqncbi_Position       (ESL_SQFILE *sqfp, off_t offset);
static void  sqncbi_Close          (ESL_SQFILE *sqfp);
static int   sqncbi_SetDigital     (ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
static int   sqncbi_GuessAlphabet  (ESL_SQFILE *sqfp, int *ret_type);
static int   sqncbi_Read           (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadInfo       (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadSequence   (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadWindow     (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
static int   sqncbi_ReadBlock      (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock);
static int   sqncbi_Echo           (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);

static int   sqncbi_IsRewindable   (const ESL_SQFILE *sqfp);
static const char *sqncbi_GetError (const ESL_SQFILE *sqfp);

/* common routines for processing ncbi database */
static int  get_offsets         (ESL_SQNCBI_DATA *ncbi, int inx, off_t *hdr, off_t *seq);
static int  inmap_ncbi          (ESL_SQFILE *sqfp);

/* parsing routines */
static int  parse_header              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_def_line            (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_seq_id              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_textseq_id          (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_object_id           (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_dbtag               (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_patent_seq_id       (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_id_pat              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_pdb_seq_id          (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_date_std            (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_string              (ESL_SQNCBI_DATA *ncbi, int max, char **str);
static int  parse_integer             (ESL_SQNCBI_DATA *ncbi, int *value);
static int  ignore_sequence_of_integer(ESL_SQNCBI_DATA *ncbi);

#define INDEX_TABLE_SIZE      1024
#define INIT_HDR_BUFFER_SIZE  2048

/* set the max residue count to 1 meg when reading a block */
#define MAX_RESIDUE_COUNT (1024 * 1024)

/*****************************************************************
 *# 1. An <ESL_SQFILE> object, in text mode.
 *****************************************************************/ 

/* Function:  esl_sqfile_Open()
 * Synopsis:  Open a sequence file for reading.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Open a sequence file <filename> for reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The .pin, .phr and .psq files are required for the
 *            open function to succeed.  Only protien version 4
 *            databases are currently supported.
 *            
 * Returns:   <eslOK> on success, and <*ret_sqfp> points to a new
 *            open <ESL_SQFILE>. Caller deallocates this object with
 *            <esl_sqfile_Close()>. 
 *            
 *            Returns <eslENOTFOUND> if <filename> can't be found or
 *            opened.  Returns <eslEFORMAT> if the file is empty, or
 *            if autodetection is attempted and the format can't be
 *            determined.  On any error condition, <*ret_sqfp> is
 *            returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqncbi_Open(char *filename, int format, ESL_SQFILE *sqfp)
{
  int         status = eslOK;	/* return status from an ESL call */
  int         len;

  uint32_t    info[4];
  char       *name = NULL;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  /* before we go any further, make sure we can handle the format */
  if (format != eslSQFILE_NCBI && format != eslSQFILE_UNKNOWN) return eslENOTFOUND;

  ncbi->fppin        = NULL;
  ncbi->fpphr        = NULL;
  ncbi->fppsq        = NULL;

  ncbi->title        = NULL;
  ncbi->timestamp    = NULL;

  ncbi->index        = 0;

  ncbi->cur_indexes  = -1;
  ncbi->hdr_indexes  = NULL;
  ncbi->seq_indexes  = NULL;

  ncbi->hdr_buf      = NULL;

  ncbi->alphasym     = NULL;

  len = strlen(filename);
  ESL_ALLOC(name, sizeof(char) * (len+5));
  strcpy(name, filename);

  /* Check the current working directory first. */
  strcpy(name+len, ".pin");
  if ((ncbi->fppin = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".phr");
  if ((ncbi->fpphr = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".psq");
  if ((ncbi->fppsq = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }

  /* make sure we are looking at a version 4 protien db.
   * the values are stored in big endian, so we will just
   * against the values in big endian format
   */

  if (fread(&info[0], sizeof(uint32_t), 3, ncbi->fppin) != 3) status = eslFAIL;
  if (info[0] != 0x4000000)                                   status = eslEFORMAT;
  if (info[1] != 0x1000000)                                   status = eslEUNIMPLEMENTED;

  if (status != eslOK) goto ERROR;
  ncbi->version = htobe32(info[0]);

  /* read the database title */
  len = htobe32(info[2]);
  ESL_ALLOC(ncbi->title, sizeof(char) * (len + 1));
  if (fread(ncbi->title, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->title[len] = 0;

  /* read the database time stamp */
  if (fread(&info[0], sizeof(uint32_t), 1, ncbi->fppin) != 1) { status = eslFAIL; goto ERROR; }
  len = htobe32(info[0]);
  ESL_ALLOC(ncbi->timestamp, sizeof(char) * (len + 1));
  if (fread(ncbi->timestamp, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->timestamp[len] = 0;

  /* read in database stats */
  if (fread(&info[0], sizeof(uint32_t), 4, ncbi->fppin) != 4) { status = eslFAIL; goto ERROR; }
  ncbi->num_seq   = htobe32(info[0]);
  ncbi->total_res = *(uint64_t *)(info+1);
  ncbi->max_seq   = htobe32(info[3]);

  /* save the offsets to the index tables */
  ncbi->hdr_off = ftell(ncbi->fppin);
  ncbi->seq_off = ncbi->hdr_off + sizeof(uint32_t) * (ncbi->num_seq + 1);

  /* allocate buffers used in parsing the database files */
  ESL_ALLOC(ncbi->hdr_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);
  ESL_ALLOC(ncbi->seq_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);

  ncbi->hdr_alloced = INIT_HDR_BUFFER_SIZE;
  ESL_ALLOC(ncbi->hdr_buf, sizeof(char) * INIT_HDR_BUFFER_SIZE);

  /* skip the first sentinal byte in the .psq file */
  fgetc(ncbi->fppsq);

  sqfp->format = eslSQFILE_NCBI;
  if ((status = inmap_ncbi(sqfp)) != eslOK) return status;

  /* initialize the function pointers for the ncbi routines */
  sqfp->position          = &sqncbi_Position;
  sqfp->close             = &sqncbi_Close;

  sqfp->set_digital       = &sqncbi_SetDigital;
  sqfp->guess_alphabet    = &sqncbi_GuessAlphabet;

  sqfp->is_rewindable     = &sqncbi_IsRewindable;

  sqfp->read              = &sqncbi_Read;
  sqfp->read_info         = &sqncbi_ReadInfo;
  sqfp->read_seq          = &sqncbi_ReadSequence;
  sqfp->read_window       = &sqncbi_ReadWindow;
  sqfp->echo              = &sqncbi_Echo;

  sqfp->read_block        = &sqncbi_ReadBlock;

  sqfp->get_error         = &sqncbi_GetError;

  if (name != NULL) free(name);

  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  sqncbi_Close(sqfp); 

  return status;
}


/* Function:  esl_sqfile_Position()
 * Synopsis:  Reposition an open sequence file to an offset.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reposition an open <sqfp> to offset <offset>.
 *            <offset> for the ncbi db format specified the sequence
 *            index, not file offset.  Both the sequence and header
 *            files are repositioned.
 *            
 * Returns:   <eslOK>     on success;
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails.
 *            On errors, the state of <sqfp> is indeterminate, and
 *            it should not be used again.
 */
static int
sqncbi_Position(ESL_SQFILE *sqfp, off_t offset)
{
  off_t    hdr_start;
  off_t    seq_start;

  int      status;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if ((status = get_offsets(ncbi, offset, &hdr_start, &seq_start)) != eslOK) return status;

  if (fseek(ncbi->fpphr, hdr_start, SEEK_SET) != 0) return eslESYS;
  if (fseek(ncbi->fppsq, seq_start, SEEK_SET) != 0) return eslESYS;

  ncbi->index = offset;

  return eslOK;
}


/* Function:  esl_sqfile_Close()
 * Synopsis:  Close a sequence file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Closes an open <sqfp>.
 *
 * Returns:   (void).
 */
static void
sqncbi_Close(ESL_SQFILE *sqfp)
{
  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->title != NULL)       free(ncbi->title);
  if (ncbi->timestamp != NULL)   free(ncbi->timestamp);

  if (ncbi->hdr_buf != NULL)     free(ncbi->hdr_buf);

  if (ncbi->hdr_indexes != NULL) free(ncbi->hdr_indexes);
  if (ncbi->seq_indexes != NULL) free(ncbi->seq_indexes);

  if (ncbi->alphasym != NULL)    free(ncbi->alphasym);

  if (ncbi->fppin != NULL) fclose(ncbi->fppin);
  if (ncbi->fpphr != NULL) fclose(ncbi->fpphr);
  if (ncbi->fppsq != NULL) fclose(ncbi->fppsq);

  ncbi->fppin        = NULL;
  ncbi->fpphr        = NULL;
  ncbi->fppsq        = NULL;

  ncbi->title        = NULL;
  ncbi->timestamp    = NULL;

  ncbi->index        = -1;

  ncbi->cur_indexes  = -1;
  ncbi->hdr_indexes  = NULL;
  ncbi->seq_indexes  = NULL;

  ncbi->hdr_buf      = NULL;

  ncbi->alphasym     = NULL;

  return;
}
/*------------------- ESL_SQFILE open/close -----------------------*/


/*****************************************************************
 *# 2. An <ESL_SQFILE> object, in digital mode [with <alphabet>]
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  esl_sqfile_SetDigital()
 * Synopsis:  Set an open <ESL_SQFILE> to read in digital mode.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Given an <ESL_SQFILE> that's already been opened,
 *            configure it to expect subsequent input to conform
 *            to the digital alphabet <abc>.
 *            
 *            Calling <esl_sqfile_Open(); esl_sqfile_SetDigital()> is
 *            equivalent to <esl_sqfile_OpenDigital()>. The two-step
 *            version is useful when you need a
 *            <esl_sqfile_GuessAlphabet()> call in between, guessing
 *            the file's alphabet in text mode before you set it to
 *            digital mode.
 *
 * Returns:   <eslOK> on success.
 */
static int
sqncbi_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
{
  return eslOK;
}

/* Function:  esl_sqfile_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open <ESL_SQFILE>.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   The only ncbi db format supported is protein.
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set to <eslAMINO>.
 */
static int
sqncbi_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)
{
  *ret_type = eslAMINO;
  return eslOK;
}
#endif /*eslAUGMENT_ALPHABET*/
/*-------------- end, digital mode ESL_SQFILE -------------------*/




/*****************************************************************
 *# 3. Miscellaneous routines 
 *****************************************************************/ 

/* Function:  esl_sqfile_IsRewindable()
 * Synopsis:  Return <TRUE> if <sqfp> can be rewound.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Returns <TRUE> if <sqfp> can be rewound (positioned 
 *            to an offset of zero), in order to read it a second
 *            time.
 */
static int
sqncbi_IsRewindable(const ESL_SQFILE *sqfp)
{
  return TRUE;
}

/* Function:  sqncbi_GetError()
 * Synopsis:  Returns error buffer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Return a pointer to the error buffer.
 */
static const char *
sqncbi_GetError(const ESL_SQFILE *sqfp)
{
  return sqfp->data.ncbi.errbuf;
}




/*****************************************************************
 *# 4. Sequence reading (sequential)
 *****************************************************************/ 

/* Function:  esl_sqio_Read()
 * Synopsis:  Read the next sequence from a file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads the next sequence from open sequence file <sqfp> into 
 *            <sq>. Caller provides an allocated and initialized <s>, which
 *            will be internally reallocated if its space is insufficient.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character; the line number that the parse
 *            error occurs on is in <sqfp->linenumber>, and an informative
 *            error message is placed in <sqfp->errbuf>. 
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
static int
sqncbi_Read(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     inx;
  int     status;

  void   *tmp;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  int     size;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end)) != eslOK) return status;

  /* figure out the sequence length */
  size = seq_end - seq_start;
  if (esl_sq_GrowTo(sq, size) != eslOK) return eslEMEM;

  /* figure out if the sequence is in digital mode or not */
  if (sq->dsq != NULL) {
    ESL_DSQ *ptr = sq->dsq + 1;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      ++ptr;
    }
    *ptr = eslDSQ_SENTINEL;
  } else {
    char *ptr = sq->seq;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      *ptr = ncbi->alphasym[(int) *ptr];
      ++ptr;
    }
    *ptr = '\0';
  }

  sq->start = 1;
  sq->end   = size - 1;
  sq->C     = 0;
  sq->W     = size - 1;
  sq->L     = size - 1;
  sq->n     = size - 1;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->eoff = -1;

  /* read in the header data */
  size = hdr_end - hdr_start;
  if (ncbi->hdr_alloced < size) {
    while (ncbi->hdr_alloced < size) ncbi->hdr_alloced += ncbi->hdr_alloced;
    ESL_RALLOC(ncbi->hdr_buf, tmp, sizeof(char) * ncbi->hdr_alloced);
  }
  if (fread(ncbi->hdr_buf, sizeof(char), size, ncbi->fpphr) != size) return eslEFORMAT;
  ncbi->hdr_ptr  = ncbi->hdr_buf;
  ncbi->hdr_fpos = hdr_start;
  ncbi->hdr_size = size;

  /* parse the ncbi header */
  if ((status = parse_header(ncbi, sq)) != eslOK) return status;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;

 ERROR:
  return eslEMEM;
}


/* Function:  esl_sqio_ReadInfo()
 * Synopsis:  Read sequence info, but not the sequence itself.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Read the next sequence from open sequence file <sqfp>,
 *            but don't store the sequence (or secondary structure).
 *            Upon successful return, <s> holds all the available 
 *            information about the sequence -- its name, accession,
 *            description, and overall length <sq->L>. 
 *            
 *            This is useful for indexing sequence files, where
 *            individual sequences might be ginormous, and we'd rather
 *            avoid reading complete seqs into memory.
 *
 * Returns:   <eslOK> on success.
 */
static int
sqncbi_ReadInfo(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;

  void   *tmp;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  int     size;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end)) != eslOK) return status;

  /* advance the sequence file pointer to point to the next sequence in case
   * the user calls a Read after a ReadInfo.  this will garuntee that the
   * header and sequence match up for the Read.
   */
  if (fseek(ncbi->fppsq, seq_end, SEEK_SET) != 0) return eslEFORMAT;

  /* figure out the sequence length */
  sq->L = seq_end - seq_start - 1;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->eoff = -1;

  /* read in the header data */
  size = hdr_end - hdr_start;
  if (ncbi->hdr_alloced < size) {
    while (ncbi->hdr_alloced < size) ncbi->hdr_alloced += ncbi->hdr_alloced;
    ESL_RALLOC(ncbi->hdr_buf, tmp, sizeof(char) * ncbi->hdr_alloced);
  }
  if (fread(ncbi->hdr_buf, sizeof(char), size, ncbi->fpphr) != size) return eslEFORMAT;
  ncbi->hdr_ptr  = ncbi->hdr_buf;
  ncbi->hdr_fpos = hdr_start;
  ncbi->hdr_size = size;

  /* parse the ncbi header */
  if ((status = parse_header(ncbi, sq)) != eslOK) return status;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;

 ERROR:
  return eslEMEM;
}


/* Function:  esl_sqio_ReadSequence()
 * Synopsis:  Read the sequence, not the sequence header.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Read the next sequence from open sequence file <sqfp>,
 *            but not the header information.  Upon successful return,
 *            <s> holds all the sequence.
 *            
 *            This is useful reading binary formats and delaying the
 *            over heads of reading the sequence name until needed by
 *            the report generator.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 *
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
sqncbi_ReadSequence(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     inx;
  int     status;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  int     size;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end)) != eslOK) return status;

  /* figure out the sequence length */
  size = seq_end - seq_start;
  if (esl_sq_GrowTo(sq, size) != eslOK) return eslEMEM;

  /* figure out if the sequence is in digital mode or not */
  if (sq->dsq != NULL) {
    ESL_DSQ *ptr = sq->dsq + 1;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      ++ptr;
    }
    *ptr = eslDSQ_SENTINEL;
  } else {
    char *ptr = sq->seq;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      *ptr = ncbi->alphasym[(int) *ptr];
      ++ptr;
    }
    *ptr = '\0';
  }

  sq->start = 1;
  sq->end   = size - 1;
  sq->C     = 0;
  sq->W     = size - 1;
  sq->L     = size - 1;
  sq->n     = size - 1;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->eoff = -1;

  /* advance the header file pointer to point to the next header in case
   * the user calls a Read after a ReadSequence.  this will garuntee that
   * the header and sequence match up for the Read.
   */
  if (fseek(ncbi->fpphr, hdr_end, SEEK_SET) != 0) return eslEFORMAT;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;
}


/* Function:  esl_sqio_ReadWindow()
 * Synopsis:  Read next window of sequence.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Returns:   <eslEUNIMPLEMENTED>.
 */
static int
sqncbi_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq)
{
  return eslEUNIMPLEMENTED;
}

/* Function:  esl_sqio_ReadBlock()
 * Synopsis:  Read the next block of sequences from a file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads a block of sequences from open sequence file <sqfp> into 
 *            <sqBlock>.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sqBlock>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character;
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
static int
sqncbi_ReadBlock(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock)
{
  int     i;
  int     size = 0;
  int     status = eslOK;

  sqBlock->count = 0;
  for (i = 0; i < sqBlock->listSize && size < MAX_RESIDUE_COUNT; ++i)
    {
      status = sqncbi_Read(sqfp, sqBlock->list + i);
      if (status != eslOK) break;
      size += sqBlock->list[i].n;
      ++sqBlock->count;
    }

  /* EOF will be returned only in the case were no sequences were read */
  if (status == eslEOF && i > 0) status = eslOK;

  return status;
}

/* Function:  esl_sqio_Echo()
 * Synopsis:  Echo a sequence's record onto output stream.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Returns:   <eslEUNIMPLEMENTED>.
 */
static int
sqncbi_Echo(ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)
{
  return eslEUNIMPLEMENTED;
}
/*------------------ end, sequential sequence input -------------*/



/* Function:  get_offsets()
 * Synopsis:  Return the header and sequence offsets
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   For sequence <inx> reads the offsets in the sequence
 *            and header files.  If <hdr> or <seq> are non-null, the
 *            offsets for the .phr and .psq files will be set.
 *
 *            Before reading the offsets from the file, check if the
 *            current offset is cached.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there is an error reading the index file.
 */
static int
get_offsets(ESL_SQNCBI_DATA *ncbi, int inx, off_t *hdr, off_t *seq)
{
  int        cnt;
  off_t      offset;

  uint32_t   start;
  uint32_t   end;

  if (inx < 0 || inx > ncbi->num_seq) return eslEINVAL;

  start = ncbi->cur_indexes;
  end   = start + INDEX_TABLE_SIZE - 1;

  if (ncbi->cur_indexes == -1 || inx < start || inx > end) {

    /* when calculating the count be sure to take into account the fact that the
     * index tables contain one index more that the number of sequences and this
     * last index is used to point to the end of the last header and sequences.
     */
    cnt = ncbi->num_seq - inx + 1;
    cnt = (cnt > INDEX_TABLE_SIZE) ? INDEX_TABLE_SIZE : cnt;

    offset = ncbi->hdr_off + (sizeof(uint32_t) * inx);
    if (fseek(ncbi->fppin, offset, SEEK_SET) != 0) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error seeking header index %ld\n", offset);
    }
    if (fread(ncbi->hdr_indexes, sizeof(uint32_t), cnt, ncbi->fppin) != cnt) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error reading header index %ld(%d)\n", offset, cnt);
    }

    offset = ncbi->seq_off + (sizeof(uint32_t) * inx);
    if (fseek(ncbi->fppin, offset, SEEK_SET) != 0) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error seeking sequence index %ld\n", offset);
    }
    if (fread(ncbi->seq_indexes, sizeof(uint32_t), cnt, ncbi->fppin) != cnt) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error reading sequence index %ld(%d)\n", offset, cnt);
    }

    ncbi->cur_indexes = inx;
  }

  inx -= ncbi->cur_indexes;
  if (hdr != NULL) *hdr = htobe32(ncbi->hdr_indexes[inx]);
  if (seq != NULL) *seq = htobe32(ncbi->seq_indexes[inx]);

  return eslOK;
}

/* Function:  inmap_ncbi()
 * Synopsis:  Set up a translation map
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Initialize the translation map used to translate a ncbi
 *            protein sequence to the internal representation used in
 *            hmmer.
 *
 * Returns:   <eslOK> on success;
 * 
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
inmap_ncbi(ESL_SQFILE *sqfp)
{
  int x, y;
  const char *ncbisym = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

  ESL_ALPHABET    *abc  = NULL;
  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) return eslEMEM;

  for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;

  /* for each letter in the ncbi alphabet, find that letter in the
   * hmmer alphabet and map the translation.
   */
  for (x = 0; x < strlen(ncbisym); ++x) {
    for (y = 0; y < strlen(abc->sym); ++y) {
      if (ncbisym[x] == abc->sym[y]) {
	sqfp->inmap[x] = y;
	break;
      }
    }

    /* there is a problem if a translation does not exist */
    if (y >= strlen(abc->sym)) return eslEFORMAT;
  }

  if (ncbi->alphasym == NULL) esl_strdup(abc->sym, -1, &ncbi->alphasym);

  esl_alphabet_Destroy(abc);

  return eslOK;
}



/*****************************************************************
 *# 5. Parsing routines
 *****************************************************************/ 

/* Function:  parse_expect()
 * Synopsis:  Expect the next bytes to parse match
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Match if the next <len> bytes to parse match the bytes
 *            in <str>.  If the bytes do not match, throw <eslEFORMAT>
 *            error.  Advance the parsers pointer.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header or if the data to parse does not match
 *            what is expected.
 */
static int
parse_expect(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %ld : 0x%X(%d)\n",
	       len, ncbi->hdr_ptr - ncbi->hdr_buf, ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  /* check the buffer matches the token string */
  c = (unsigned char *) str;
  while (len--) {
    if (*ncbi->hdr_ptr != *c) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting 0x%X found 0x%X at %ld : 0x%X(%d)\n",
	       *ncbi->hdr_ptr, *c, ncbi->hdr_ptr - ncbi->hdr_buf, ncbi->hdr_fpos, ncbi->hdr_size); 
    }
    ncbi->hdr_ptr++;
    c++;
  }

  return eslOK;
}

/* Function:  parse_accept()
 * Synopsis:  Check if the next bytes to parse match
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Check if the next <len> bytes to parse match the bytes
 *            in <str>.  If the bytes match, they are consumed and the
 *            parsers pointer is advanced.
 *
 * Returns:   <eslOK> on success
 *            <eslEFORMAT> if the bytes to not match.
 */
static int
parse_accept(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  int i;
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* check the buffer matches the token string */
  if (ncbi->hdr_ptr + len > limit)  return eslEFORMAT;

  /* verify the buffer matches the token string without advancing
   * the buffer pointers until we have a complete match.
   */
  c = (unsigned char *) str;
  for (i = 0; i < len; ++i) {
    if (ncbi->hdr_ptr[i] != c[i])   return eslEFORMAT;
  }
  ncbi->hdr_ptr += len;

  return eslOK;
}

/* Function:  parse_peek()
 * Synopsis:  Peek at the next byte to parse
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Return the next characer to be parsed without advancing the
 *            parsers pointer.
 *
 * Returns:   <eslOK> on success
 *            <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_peek(ESL_SQNCBI_DATA *ncbi, unsigned char *c)
{
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + 1 > limit)    return eslEFORMAT;

  *c = *ncbi->hdr_ptr;

  return eslOK;
}

/* Function:  parse_consume()
 * Synopsis:  Copies bytes from the header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Copies <len> bytes from the header to the buffer supplied by
 *            <str> if non-null.  Adcance the parser pointer.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_consume(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  int i;
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %ld : 0x%X(%d)\n",
	       len, ncbi->hdr_ptr - ncbi->hdr_buf, ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  /* copy the characters in the buffer to <str> */
  c = (unsigned char *) str;
  for (i = 0; i < len; ++i) {
    if (c != NULL) *c++ = *ncbi->hdr_ptr++;
  }

  return eslOK;
}

/* Function:  parse_advance()
 * Synopsis:  Advance the parser pointer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Advance the parser pointer <len> bytes.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_advance(ESL_SQNCBI_DATA *ncbi, int len)
{
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %ld : 0x%X(%d)\n",
	       len, ncbi->hdr_ptr - ncbi->hdr_buf, ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  ncbi->hdr_ptr += len;

  return eslOK;
}

/* Function:  parse_header()
 * Synopsis:  Parse the ncbi db header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Parse a ncbi database header.  This routine implements
 *            a recursive descent parser for the ASN.1 definition of
 *            a blast database header filling in <sq>.
 *
 *            The blast db header can have multiple definitions defined
 *            within it.  Only the information from the first usable
 *            defition will be used.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character.
 */
static int
parse_header(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int status;

  unsigned char c;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  if (parse_peek(ncbi, &c) != eslOK)                          return eslEFORMAT;

  /* parse the different seq id structures */
  while (c != 0x00) {
    if ((status = parse_def_line(ncbi, sq)) != eslOK)         return status;

    if (parse_peek(ncbi, &c) != eslOK)                        return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}

/* Function:  parse_def_line()
 * Synopsis:  Parse the Blast-def-line definition
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Blast-def-line ::= SEQUENCE {
 * 	title       VisibleString       OPTIONAL,  -- simple title
 * 	seqid       SEQUENCE OF Seq-id,            -- Regular NCBI Seq-Id
 * 	taxid       INTEGER             OPTIONAL,  -- taxonomy id
 * 	memberships SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
 * 	links       SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
 * 	other-info  SEQUENCE OF INTEGER OPTIONAL   -- for future use (probably genomic sequences)
 * }
 */
static int
parse_def_line(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  char *buf;
  int   taxid;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional title */
  sq->desc[0] = 0;
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    free(sq->desc);
    sq->dalloc = strlen(buf) + 1;
    sq->desc   = buf;
  }

  /* look for sequence id structure */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_seq_id(ncbi, sq)) != eslOK)             return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional taxonomy id */
  sq->tax_id = -1;
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, &taxid)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    sq->tax_id = taxid;
  }

  /* look for an optional memberships */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional links */
  if (parse_accept(ncbi, "\xa4\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional other info */
  if (parse_accept(ncbi, "\xa5\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_seq_id()
 * Synopsis:  Parse the Blast-def-line definition
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Seq-id ::= CHOICE {
 *     local             Object-id ,       -- local use
 *     gibbsq            INTEGER ,         -- Geninfo backbone seqid
 *     gibbmt            INTEGER ,         -- Geninfo backbone moltype
 *     giim              Giimport-id ,     -- Geninfo import id
 *     genbank           Textseq-id ,
 *     embl              Textseq-id ,
 *     pir               Textseq-id ,
 *     swissprot         Textseq-id ,
 *     patent            Patent-seq-id ,
 *     other             Textseq-id ,      -- for historical reasons, 'other' = 'refseq'
 *     general           Dbtag ,           -- for other databases
 *     gi                INTEGER ,         -- GenInfo Integrated Database
 *     ddbj              Textseq-id ,      -- DDBJ
 *     prf               Textseq-id ,      -- PRF SEQDB
 *     pdb               PDB-seq-id ,      -- PDB sequence
 *     tpg               Textseq-id ,      -- Third Party Annot/Seq Genbank
 *     tpe               Textseq-id ,      -- Third Party Annot/Seq EMBL
 *     tpd               Textseq-id ,      -- Third Party Annot/Seq DDBJ
 *     gpipe             Textseq-id ,      -- Internal NCBI genome pipeline processing ID
 *     named-annot-track Textseq-id        -- Internal named annotation tracking ID
 * }
 */
static int
parse_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  unsigned char c;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)           return eslEFORMAT;

  if (parse_consume(ncbi, &c, 1) != eslOK)                  return eslEFORMAT;

  /* parse the different seq id structures */
  while (c != 0x00) {
    if (parse_expect(ncbi, "\x80", 1) != eslOK)             return eslEFORMAT;
    switch (c) {
    case 0xa0: /* LOCAL */
      status = parse_object_id(ncbi, sq);
      break;
    case 0xa1: /* GIBBSQ */
    case 0xa2: /* GIBBMT */
      status = parse_integer(ncbi, NULL);
      break;
    case 0xa3: /* GIIM */
      return eslEFORMAT;
      break;
    case 0xa4: /* GENBANK */
    case 0xa5: /* EMBL */
    case 0xa6: /* PIR */
    case 0xa7: /* SWISSPROT */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xa8: /* PATENT */
      status = parse_patent_seq_id(ncbi, sq);
      break;
    case 0xa9: /* OTHER */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xaa: /* GENERAL */
      status = parse_dbtag(ncbi, sq);
      break;
    case 0xab: /* GI */
      status = parse_integer(ncbi, NULL);
      break;
    case 0xac: /* DDBJ */
    case 0xad: /* PRF */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xae: /* PDB */
      status = parse_pdb_seq_id(ncbi, sq);
      break;
    case 0xaf: /* TPG */
    case 0xb0: /* TPE */
    case 0xb1: /* TPD */
    case 0xb2: /* GPIPE */
    case 0xb3: /* NAMED ANNOT TRACK */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    default:
      status = eslEFORMAT;
    }

    if (status != eslOK)                                    return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)         return eslEFORMAT;
    if (parse_consume(ncbi, &c, 1)        != eslOK)         return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (c != 0x00 || parse_expect(ncbi, "\x00", 1) != eslOK)  return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_textseq_id()
 * Synopsis:  Parse the general text header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Textseq-id ::= SEQUENCE {
 *     name      VisibleString OPTIONAL ,
 *     accession VisibleString OPTIONAL ,
 *     release   VisibleString OPTIONAL ,
 *     version   INTEGER       OPTIONAL
 * }
 */
static int
parse_textseq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  char *buf;
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional name */
  if (sq != NULL) sq->name[0] = 0;
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
    if (sq != NULL) {
      free(sq->name);
      sq->nalloc = strlen(buf) + 1;
      sq->name   = buf;
    }
  }

  /* look for an optional accession */
  if (sq != NULL) sq->acc[0] = 0;
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    if (sq != NULL) {
      free(sq->acc);
      sq->aalloc = strlen(buf) + 1;
      sq->acc    = buf;
    }
  }

  /* look for an optional release */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional version */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_dbtag()
 * Synopsis:  Parse the a general db tag
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Dbtag ::= SEQUENCE {
 *     db  VisibleString ,     -- name of database or system
 *     tag Object-id           -- appropriate tag
 * }
 */
static int
parse_dbtag(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an db name */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;

  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for a tag object */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_object_id(ncbi, sq)) != eslOK)          return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_patent_seq_id()
 * Synopsis:  Parse the patent header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Patent-seq-id ::= SEQUENCE {
 *     seqid INTEGER ,          -- number of sequence in patent
 *     cit   Id-pat             -- patent citation
 * }
 */
static int
parse_patent_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a seqid */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_integer(ncbi, NULL)) != eslOK)          return status;

  /* look for a patent citation object */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_id_pat(ncbi, sq)) != eslOK)             return status;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_id_pat()
 * Synopsis:  Parse the patent citation
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Id-pat ::= SEQUENCE {                         -- just to identify a patent
 *     country  VisibleString ,                  -- Patent Document Country
 *     id       CHOICE {
 *         number     VisibleString ,            -- Patent Document Number
 *         app-number VisibleString              -- Patent Doc Appl Number
 *     } ,
 *     doc-type VisibleString         OPTIONAL   -- Patent Doc Type
 * }
 */
static int
parse_id_pat(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a country */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;

  /* look for an id */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;

  /* the id is a choice of two strings */

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional taxonomy id */
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    status = parse_string(ncbi, 0, NULL);
  } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    status = parse_string(ncbi, 0, NULL);
  } else {
    status = eslEFORMAT;
  }
  if (status != eslOK)                                        return status;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for a doc type */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_object_id()
 * Synopsis:  Parse a generic sequence id
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Object-id ::= CHOICE {
 *     id  INTEGER ,
 *     str VisibleString
 * }
 */
static int
parse_object_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* look for an optional taxonomy id */
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    status = parse_integer(ncbi, NULL);
  } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    status = parse_string(ncbi, 0, NULL);
  } else {
    status = eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (status == eslOK) {
    status = parse_expect(ncbi, "\x00\x00", 2);
  }

  return status;
}


/* Function:  parse_pdb_seq_id()
 * Synopsis:  Parse a PDB sequence
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * PDB-seq-id ::= SEQUENCE {
 *     mol   PDB-mol-id ,              -- the molecule name
 *     chain INTEGER ,                 -- a single ASCII character, chain id
 *     rel   Date         OPTIONAL }   -- release date, month and year
 *
 * Date ::= CHOICE {
 *     str   VisibleString ,           -- for those unparsed dates
 *     std   Date-std                  -- use this if you can
 * }
 */
static int
parse_pdb_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an pdb mol id */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for chain */
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional date */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
      status = parse_string(ncbi, 0, NULL);
    } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
      status = parse_date_std(ncbi, sq);
    } else {
      status = eslEFORMAT;
    }
    if (status != eslOK)                                      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_date_std()
 * Synopsis:  Parse the data structure
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Date-std ::= SEQUENCE {              -- NOTE: this is NOT a unix tm struct
 *     year   INTEGER ,                 -- full year (including 1900)
 *     month  INTEGER       OPTIONAL ,  -- month (1-12)
 *     day    INTEGER       OPTIONAL ,  -- day of month (1-31)
 *     season VisibleString OPTIONAL ,  -- for "spring", "may-june", etc
 *     hour   INTEGER       OPTIONAL ,  -- hour of day (0-23)
 *     minute INTEGER       OPTIONAL ,  -- minute of hour (0-59)
 *     second INTEGER       OPTIONAL    -- second of minute (0-59)
 * }
 */
static int
parse_date_std(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a year */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_integer(ncbi, NULL)) != eslOK)          return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional month */
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional day */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional season */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional hour */
  if (parse_accept(ncbi, "\xa4\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional minute */
  if (parse_accept(ncbi, "\xa5\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional second */
  if (parse_accept(ncbi, "\xa6\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_string()
 * Synopsis:  Parse a visible string
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads an string from the header stream.  The arguement <max>
 *            specified the maximum number of characters to save.  If <max>
 *            is -1, the entire string will be saved.
 *
 *            The string will always be zero terminated.
 *
 *            If <str> is non null, the parsed integer will be placed
 *            in the pointer.  The calling routine is responsible for
 *            freeing the allocated memory.
 *
 * Returns:   <eslOK> on success.
 *            <eslEMEM> if there's a memory allocation error.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
parse_string(ESL_SQNCBI_DATA *ncbi, int max, char **str)
{
  int n;
  int len;
  int status;

  char *v  = NULL;

  unsigned char  x;
  unsigned char  c;
  unsigned char *ptr;

  if (parse_expect(ncbi, "\x1a", 1) != eslOK)  return eslEFORMAT;

  /* the next byte is the length of the string.  if the length is
   * less than 128, then this is the true length; otherwise this
   * length describes the number of bytes necessary to hold the
   * true length of the string in the lower 7 bits.
   */
  if (parse_consume(ncbi, &c, 1) != eslOK)     return eslEFORMAT;
  if (c < 128) {
    n = c;
  } else {
    c = c & 0x7f;
    if (c > sizeof(n))                                 return eslEFORMAT;

    n = 0;
    while (c > 0) {
      if (parse_consume(ncbi, &x, 1) != eslOK) return eslEFORMAT;
      n = (n << 8) + (unsigned int) x;
      --c;
    }
  }

  /* validate the length of the string */
  ptr = ncbi->hdr_ptr;
  if (parse_advance(ncbi, n) != eslOK)         return eslEFORMAT;

  /* now that we have the length of the string, check how much
   * of it (if any) we need to save.
   */
  if (str != NULL && max != 0) {
    if (max == -1 || max > n)  len = n;
    else                       len = max - 1;

    ESL_ALLOC(v, sizeof(char) * (len + 1));
    memcpy(v, ptr, len);
    v[len] = 0;

    *str = v;
  }

  return eslOK;

 ERROR:
  if (v != NULL) free(v);
  return status;
}


/* Function:  parse_integer()
 * Synopsis:  Parse an integer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads an integer from the header stream.  If the integer is
 *            more bytes than the native int format, the most significant
 *            bytes will be lost.
 *
 *            If <value> is non null, the parsed integer will be placed
 *            in the pointer.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
parse_integer(ESL_SQNCBI_DATA *ncbi, int *value)
{
  int n;

  unsigned char  c;
  unsigned char *ptr;

  if (parse_expect(ncbi, "\x02", 1) != eslOK) return eslEFORMAT;

  /* get the length of the integer */
  if (parse_peek(ncbi, &c) != eslOK)          return eslEFORMAT;
  ptr = ncbi->hdr_ptr + 1;

  /* advance past the integer to make sure the buffer holds all
   * of the integer.  the pointer <ptr> points the the start of
   * the integer.
   */
  if (parse_advance(ncbi, c + 1) != eslOK)    return eslEFORMAT;

  n = 0;
  while (c > 0) {
    n = (n << 8) + (unsigned int) *ptr++;
    --c;
  }

  if (value != NULL) *value = n;

  return eslOK;
}


/* Function:  ignore_sequence_of_integer()
 * Synopsis:  Skip over the sequence of integers
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Skip over a sequence of integers.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
ignore_sequence_of_integer(ESL_SQNCBI_DATA *ncbi)
{
  int status;
  unsigned char c;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)      return eslEFORMAT;

  if (parse_peek(ncbi, &c) != eslOK)                   return eslEFORMAT;
  while (c == 0x02) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK) return status;
    if (parse_peek(ncbi, &c) != eslOK)                 return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)      return eslEFORMAT;

  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
