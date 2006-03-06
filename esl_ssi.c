/* esl_ssi.c
 * simple sequence indices: fast lookup in large sequence files by keyword.
 * 
 * SVN $Id$
 * adapted from squid's ssi.c
 * SRE, Thu Mar  2 18:46:02 2006 [St. Louis]
 */

#include <esl_config.h>

#include <easel.h>
#include <esl_ssi.h>

static uint32_t v20magic = 0xf3f3e9b1; /* SSI 1.0: "ssi1" + 0x80808080 */
static uint32_t v20swap  = 0xb1e9f3f3; /* byteswapped */

static int  load_indexfile(ESL_SSI *ssi);
static void free_ssi(ESL_SSI *ssi);
static int  binary_search(ESL_SSI *ssi, char *key, uint32_t klen, off_t base, 
			  uint32_t recsize, uint32_t maxidx);

/*****************************************************************
 * Using an existing SSI index
 *****************************************************************/ 

/* Function:  esl_ssi_Open()
 * Incept:    SRE, Mon Mar  6 10:52:42 2006 [St. Louis]
 *
 * Purpose:   Open the SSI file <filename>; return ptr to a new
 *            <ESL_SSI> object in <ret_ssi>.
 *
 * Returns:   <eslOK> on success, and <ret_ssi> contains a new object;
 *            caller is responsible for free'ing it with <esl_ssi_Close()>.
 *
 *            <eslENOTFOUND> if <filename> is not found, or can't be opened
 *            for reading.
 *            <eslEFORMAT> if it's not an SSI file, according to the magic
 *            number at the start of the file.
 *            <eslERANGE> if it has 64-bit file offsets, but we're on a system
 *            that doesn't support 64-bit file offsets.
 *            <eslECORRUPT> if some other problem is found in initialization.
 *
 *            On any of these "normal" errors, <ret_ssi> is returned NULL.
 *            
 * Throws:    <eslEMEM> on allocation errors.
 */
int
esl_ssi_Open(char *filename, ESL_SSI **ret_ssi)
{
  SSI  *ssi;
  int   status;

  *ret_ssi = NULL;

  ESL_MALLOC(ssi, sizeof(SSI));
  if ((ssi->fp = fopen(filename, "rb")) == NULL) {
    free(ssi);
    return eslENOTFOUND;
  }

  if ((status = load_indexfile(ssi)) == eslOK)
    { *ret_ssi = ssi; }
  return status;
}

/* load_indexfile():
 *    given an SSI with an open and positioned 
 *    stream <ssi->fp> -- but no other data loaded -- read the next index
 *    in from disk. We use this routine without its <_Open()> wrapper
 *    as part of the external mergesort when creating large indices.
 */
static int
load_indexfile(ESL_SSI *ssi)
{
  uint32_t     magic;		/* magic number that starts the SSI file */
  uint16_t     i;		/* counter over files */
  int          status;		/* overall return status if an error is thrown */

  status = eslECORRUPT;    /* pessimistic default: almost every kind of error is a broken file error */

  ssi->filename   = NULL;
  ssi->fileformat = NULL;
  ssi->fileflags  = NULL;
  ssi->bpl        = NULL;
  ssi->rpl        = NULL;
  ssi->nfiles     = 0;          

  /* Read the magic number: make sure it's an SSI file, and determine
   * whether it's byteswapped.
   */
  if (esl_fread_i32(ssi->fp, &magic) != eslOK)   {status = eslEFORMAT;  goto FAILURE; }
  if (magic != v20magic && magic != v20swap)     {status = eslEFORMAT;  goto FAILURE; }

  /* Determine what kind of offsets (32 vs. 64 bit) are stored in the file.
   * If we can't deal with 64-bit file offsets, get out now. 
   */
  if (esl_fread_i32(ssi->fp, &(ssi->flags)) != eslOK) goto FAILURE; 
  ssi->imode = (sfp->flags & eslSSI_USE64_INDEX) ? 64 : 32;
  ssi->smode = (sfp->flags & eslSSI_USE64) ?       64 : 32;

  if (sizeof(off_t) != 8 && (ssi->imode == 64 || ssi->smode == 64))
    { status = eslERANGE; goto FAILURE; }

  /* The header data.
   */
  if (esl_fread_i16(ssi->fp, &(ssi->nfiles))     != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->nprimary))   != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->nsecondary)) != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->flen))       != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->plen))       != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->slen))       != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->frecsize))   != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->precsize))   != eslOK) goto FAILURE;
  if (esl_fread_i32(ssi->fp, &(ssi->srecsize))   != eslOK) goto FAILURE;
  
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->foffset)) != eslOK) goto FAILURE;
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->poffset)) != eslOK) goto FAILURE;
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->soffset)) != eslOK) goto FAILURE;

  /* The file information.
   *
   * We expect the number of files to be small, so reading it once
   * should be advantageous overall. If SSI ever had to deal with
   * large numbers of files, you'd probably want to read file
   * information on demand.
   * 
   * Failures of malloc's are internal errors (thrown), not "normal"
   * return codes; throw them immediately. That requires free'ing the
   * ssi structure.
   */
  if (ssi->nfiles == 0) goto FAILURE;

  ssi->filename = malloc(sizeof(char *) * ssi->nfiles);
  if (ssi->filename == NULL) 
    { esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }
  for (i = 0; i < ssi->nfiles; i++)
    ssi->filename[i] = NULL; 

  ssi->fileformat = malloc(sizeof(uint32_t) * ssi->nfiles);
  if (ssi->fileformat == NULL) 
    { esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }

  ssi->fileflags = malloc(sizeof(uint32_t) * ssi->nfiles);
  if (ssi->fileflags == NULL)
    { esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }

  ssi->bpl = malloc(sizeof(uint32_t) * ssi->nfiles);
  if (ssi->bpl == NULL)
    { esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }
  
  ssi->rpl = malloc(sizeof(uint32_t) * ssi->nfiles);
  if (ssi->rpl== NULL)
    { esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }

  /* (most) mallocs done, now we read.
   */
  for (i = 0; i < ssi->nfiles; i++) 
    {
      /* We have to explicitly position, because header and file 
       * records may expand in the future; frecsize and foffset 
       * give us forwards compatibility. 
       */ 
      if (fseeko(ssi->fp, ssi->foffset + (n * ssi->frecsize), SEEK_SET) != 0) goto FAILURE;

      ssi->filename[i] = malloc(sizeof(char)* ssi->flen);
      if (ssi->filename[i] == NULL) 
	{ esl_ssi_Close(ssi); ESL_ERROR(eslEMEM, "malloc failed"); }

      if (fread(ssi->filename[i],sizeof(char),ssi->flen, ssi->fp)!=ssi->flen) goto FAILURE;
      if (esl_fread_i32(ssi->fp, &(ssi->fileformat[i])))                      goto FAILURE;
      if (esl_fread_i32(ssi->fp, &(ssi->fileflags[i])))                       goto FAILURE;
      if (esl_fread_i32(ssi->fp, &(ssi->bpl[i])))                             goto FAILURE;
      if (esl_fread_i32(ssi->fp, &(ssi->rpl[i])))                             goto FAILURE;
    }
  
  /* Success.
   */
  return eslOK;			

  /* "Normal" failure (something wrong with the file - in userland)
   */
 FAILURE:
  esl_ssi_Close(ssi);
  return status;
}

/* Function:  esl_ssi_Close()
 * Incept:    SRE, Mon Mar  6 13:40:17 2006 [St. Louis]
 *
 * Purpose:   Close an open SSI index <ssi>.
 */
void
esl_ssi_Close(ESL_SSI *ssi)
{
  if (ssi != NULL) { 
    free_ssi(ssi);
    if (ssi->fp != NULL) fclose(ssi->fp);
    free(ssi);
  }
}  
/* free_ssi():
 * free the innards of <ssi>, without 
 * destroying the structure or closing the stream.
 */
static void
free_ssi(ESL_SSI *ssi)
{
  int i;

  if (ssi->filename != NULL) {
    for (i = 0; i < ssi->nfiles; i++) 
      if (ssi->filename[i] != NULL) free(ssi->filename[i]);
    free(ssi->filename);
  }
  if (ssi->fileformat != NULL) free(ssi->fileformat);
  if (ssi->fileflags  != NULL) free(ssi->fileflags);
  if (ssi->bpl        != NULL) free(ssi->>bpl);
  if (ssi->rpl        != NULL) free(ssi->rpl);
}


/* Function: esl_ssi_GetOffsetByName()
 * Date:     SRE, Sun Dec 31 13:55:31 2000 [St. Louis]
 *
 * Purpose:  Looks up the string <key> in index <ssi>.
 *           <key> can be either a primary or secondary key. If <key>
 *           is found, <ret_fh> contains a unique handle on
 *           the file that contains <key> (suitable for an <esl_ssi_FileInfo()>
 *           call, or for comparison to the handle of the last file
 *           that was opened for retrieval), and <ret_offset> contains
 *           the offset of the sequence record in that file.
 *           
 * Args:     ssi         - open index file
 *           key         - name to search for
 *           ret_fh      - RETURN: handle on file that key is in
 *           ret_offset  - RETURN: offset of the start of that key's record
 *
 * Returns:  <eslOK> on success.
 *           <eslECORRUPT> if an fread() or fseeko() fails, which almost
 *           certainly reflects some kind of corruption of the file.
 */
int
esl_ssi_GetOffsetByName(ESL_SSI *ssi, char *key, int *ret_fh, off_t *ret_offset)
{
  int         status;
  sqd_uint16  fnum;

  /* Look in the primary keys.
   */
  status = binary_search(ssi, key, ssi->plen, ssi->poffset, ssi->precsize,
			 ssi->nprimary);
  if (status == eslOK) {		
    /* We found it as a primary key; get our data & return.
     */
    if (esl_fread_i16(ssi->fp, &fnum) != eslOK) return eslECORRUPT;
    *ret_fh = (int) fnum;
    if (esl_fread_offset(ssi->fp, ssi->smode, ret_offset) != eslOK) return eslECORRUPT;
    return 0;	/* success! (we don't need the other key data) */
  } else if (status == eslENOTFOUND) {
    /* Not in the primary keys? OK, try the secondary keys.
     */
    if (ssi->nsecondary > 0) {
      char *pkey;
      status = binary_search(ssi, key, ssi->slen, ssi->soffset, ssi->srecsize,
			     ssi->nsecondary);
      if (status != eslOK) return status;

      /* We have the secondary key; flip to its primary key, then look that up.
       */
      ESL_MALLOC(pkey, sizeof(char) * ssi->plen);
      if (fread(pkey, sizeof(char), ssi->plen, ssi->fp) != ssi->plen) 
	{ free(pkey); return eslECORRUPT; }
      status = SSIGetOffsetByName(ssi, pkey, ret_fh, ret_offset);
      free(pkey);
    }
    return status;
  } else 
    return status;		
  /*NOTREACHED*/
}

int
esl_ssi_GetOffsetByNumber(ESL_SSI *ssi, int nkey, int *ret_fh, off_t *ret_offset)
{

}

int
esl_ssi_GetSubseqOffset(ESL_SSI *ssi, char key, long requested_start,
			int *ret_fh, off_t *record_offset, off_t *data_offset, 
			long *ret_actual_start)
{

}

int
esl_ssi_SetFilePosition(FILE *fp, off_t *offset)
{

}

int
esl_ssi_FileInfo(ESL_SSI *ssi, int fh, char **ret_filename, int *ret_format)
{

}



/* binary_search()
 * Date:     SRE, Sun Dec 31 16:05:03 2000 [St. Louis]
 *
 * Purpose:  Find <key> in an SSI index, by a binary search
 *           in an alphabetically sorted list of keys. If successful,
 *           return <eslOK>, and the index file is positioned to read
 *           the rest of the data for that key. If unsuccessful, 
 *           return <eslFAIL>.
 *
 * Args:     ssi    - an open ESL_SSI
 *           key    - key to find
 *           klen   - key length to allocate (plen or slen from ssi)
 *           base   - base offset (poffset or soffset)
 *           recsize - size of each key record in bytes (precsize or srecsize)
 *           maxidx  - # of keys (nprimary or nsecondary)
 *
 * Returns:  <eslOK> on success, and leaves file positioned for reading remaining
 *           data for the key. 
 *           
 *           <eslENOTFOUND> if <key> is not found.
 *           <eslECORRUPT>  if an fread() or fseeko() fails, probably indicating
 *           some kind of corruption of the index file.
 *
 * Throws:   <eslEMEM> on allocation failure.
 *           
 */
static int
binary_search(ESL_SSI *ssi, char *key, uint32_t klen, off_t base, 
	      uint32_t recsize, uint32_t maxidx)
{
  char        *name;
  uint32_t     left, right, mid;
  int          cmp;
  int          status;
  
  if (maxidx == 0) return eslENOTFOUND; /* special case: empty index */

  ESL_MALLOC(name, (sizeof(char)*klen));

  left  = 0;
  right = maxidx-1;
  while (1) {			/* A binary search: */
    mid   = (left+right) / 2;	/* careful here. left+right potentially overflows if
				   we didn't limit unsigned vars to signed ranges. */
    if (fseeko(ssi->fp, base + recsize*mid, SEEK_SET) != 0)
      { free(name); return eslECORRUPT; }
    if (fread(name, sizeof(char), klen, ssi->fp) != klen) 
      { free(name); return eslECORRUPT; }

    cmp = strcmp(name, key);
    if      (cmp == 0) break;	          /* found it!              */
    else if (left >= right)	          /* oops, missed it; fail  */
      { free(name); return eslENOTFOUND; }
    else if (cmp < 0)       left  = mid+1; /* it's right of mid     */
    else if (cmp > 0) {
      if (mid == 0) { free(name); return eslENOTFOUND; } /* special case, beware */
      else right = mid-1;                  /* it's left of mid      */
    }
  }

  free(name);
  return eslOK; /* and ssi->fp is positioned... */
}



/*****************************************************************
 * Functions for platform-independent binary files
 *****************************************************************/ 

/* Function:  esl_byteswap()
 *
 * Purpose:   Swap between big-endian and little-endian, in place.
 */
void
esl_byteswap(char *swap, int nbytes)
{
  int  x;
  char byte;
  
  for (x = 0; x < nbytes / 2; x++)
    {
      byte = swap[nbytes - x - 1];
      swap[nbytes - x - 1] = swap[x];
      swap[x] = byte;
    }
}

/* Function:  esl_ntoh16()
 *
 * Purpose:   Convert a 2-byte integer from network-order to host-order,
 *            and return it.
 *            
 *            <esl_ntoh32()> and <esl_ntoh64()> do the same, but for 4-byte
 *            and 8-byte integers, respectively.
 */
uint16_t
esl_ntoh16(uint16_t netshort)
{
#ifdef WORDS_BIGENDIAN
  return netshort;
#else
  esl_byteswap((char *) &netshort, 2);
  return netshort;
#endif
}
uint32_t
esl_ntoh32(uint32_t netlong)
{
#ifdef WORDS_BIGENDIAN
  return netlong;
#else
  esl_byteswap((char *) &netlong, 4);
  return netlong;
#endif
}
uint64_t
esl_ntoh64(uint64_t net_int64)
{
#ifdef WORDS_BIGENDIAN
  return net_int64;
#else
  esl_byteswap((char *) &net_int64, 8);
  return net_int64;
#endif
}

/* Function:  esl_hton16()
 *
 * Purpose:   Convert a 2-byte integer from host-order to network-order, and
 *            return it.
 * 
 *            <esl_hton32()> and <esl_hton64()> do the same, but for 4-byte
 *            and 8-byte integers, respectively.
 */
uint16_t
esl_hton16(uint16_t hostshort)
{
#ifdef WORDS_BIGENDIAN
  return hostshort;
#else
  esl_byteswap((char *) &hostshort, 2);
  return hostshort;
#endif
}
uint32_t
esl_hton32(uint32_t hostlong)
{
#ifdef WORDS_BIGENDIAN
  return hostlong;
#else
  esl_byteswap((char *) &hostlong, 4);
  return hostlong;
#endif
}
uint64_t
esl_hton64(uint64_t host_int64)
{
#ifdef WORDS_BIGENDIAN
  return host_int64;
#else
  esl_byteswap((char *) &host_int64, 8);
  return host_int64;
#endif
}


/* Function:  esl_fread_i16()
 *
 * Purpose:   Read a 2-byte network-order integer from <fp>, convert to
 *            host order, leave it in <ret_result>.
 *            
 *            <esl_fread_i32()> and <esl_fread_i64()> do the same, but
 *            for 4-byte and 8-byte integers, respectively.
 *
 * Returns:   <eslOK> on success, and <eslFAIL> on <fread()> failure.
 */
int
esl_fread_i16(FILE *fp, uint16_t *ret_result)
{
  uint16_t result;
  if (fread(&result, sizeof(uint16_t), 1, fp) != 1) return eslFAIL;
  *ret_result = esl_ntoh16(result);
  return eslOK;
}
int
esl_fread_i32(FILE *fp, uint32_t *ret_result)
{
  uint32_t result;
  if (fread(&result, sizeof(uint32_t), 1, fp) != 1) return eslFAIL;
  *ret_result = esl_ntoh32(result);
  return eslOK;
}
int
esl_fread_i64(FILE *fp, uint64_t *ret_result)
{
  uint64_t result;
  if (fread(&result, sizeof(uint64_t), 1, fp) != 1) return eslFAIL;
  *ret_result = esl_ntoh64(result);
  return eslOK;
}


/* Function:  esl_fwrite_i16()
 *
 * Purpose:   Write a 2-byte host-order integer <n> to stream <fp>
 *            in network order.
 *            
 *            <esl_fwrite_i32()> and <esl_fwrite_i64()> do the same, but
 *            for 4-byte and 8-byte integers, respectively.
 *
 * Returns:   <eslOK> on success, and <eslFAIL> on <fwrite()> failure.
 */
int
esl_fwrite_i16(FILE *fp, uint16_t n)
{
  n = esl_hton16(n);
  if (fwrite(&n, sizeof(uint16_t), 1, fp) != 1) return eslFAIL;
  return eslOK;
}
int
esl_fwrite_i32(FILE *fp, uint32_t n)
{
  n = esl_hton32(n);
  if (fwrite(&n, sizeof(uint32_t), 1, fp) != 1) return eslFAIL;
  return eslOK;
}
int
esl_fwrite_i64(FILE *fp, uint64_t n)
{
  n = esl_hton64(n);
  if (fwrite(&n, sizeof(uint64_t), 1, fp) != 1) return eslFAIL;
  return eslOK;
}

/* Function:  esl_fread_offset()
 * Incept:    SRE, Fri Mar  3 13:19:41 2006 [St. Louis]
 *
 * Purpose:   Read a file offset from the stream <fp> (which would usually
 *            be a save file), and store it in <ret_offset>.
 *            
 *            Offsets may have been saved by a different machine
 *            than the machine that reads them. The writer and the reader
 *            may differ in byte order and in width (<sizeof(off_t)>). 
 *            
 *            Byte order is dealt with by saving offsets in 
 *            network byte order, and converting them to host byte order
 *            when they are read (if necessary). 
 *            
 *            Width is dealt with by the <mode> argument, which must
 *            be either 32 or 64, specifying that the saved offset is a
 *            32-bit versus 64-bit <off_t>. If the reading host <off_t> width 
 *            matches the <mode> of the writer, no problem. If <mode> is
 *            32 but the reading host has 64-bit <off_t>, this is also
 *            no problem; the conversion is handled. If <mode> is 64
 *            but the reading host has only 32-bit <off_t>, we cannot
 *            guarantee that we have sufficient dynamic range to represent
 *            the offset, so we throw a fatal error.
 *
 * Returns:   <eslOK> on success; <eslFAIL> on a read failure.
 *
 * Throws:    <eslEINVAL> if mode is something other than 32 or 64;
 *            <eslEINCOMPAT> if mode is 64 but host <off_t> is only 32.
 */
int			
esl_fread_offset(FILE *fp, int mode, off_t *ret_offset)
{

  if      (mode == 64 && sizeof(off_t) == 8)
    return esl_fread_i64(fp, ret_offset);
  else if (mode == 32 && sizeof(off_t) == 4)
    return esl_fread_i32(fp, ret_offset);
  else if (mode == 32 && sizeof(off_t) == 8)
    {
      esl_uint32 x;
      if (esl_fread_i32(fp, &x) != eslOK) return eslFAIL;
      *ret_offset = (uint64_t) x;
      return eslOK;
    }

  if (mode != 32 && mode != 64)
    ESL_ERROR(eslEINVAL, "mode must be 32 or 64");
  else
    ESL_ERROR(eslEINCOMPAT, "can't read 64-bit off_t on this 32-bit host");
  /*UNREACHED*/
  return 1;
}

/* Function:  esl_fwrite_offset()
 * Incept:    SRE, Fri Mar  3 13:35:04 2006 [St. Louis]
 *
 * Purpose:   Portably write (save) <offset> to the stream <fp>, in network
 *            byte order. 
 *
 * Returns:   <eslOK> on success; <eslFAIL> on write failure.
 *
 * Throws:    <eslESYS> if <off_t> is something other than a 32-bit or
 *            64-bit integer on this machine, in which case we don't know
 *            how to deal with it portably.
 */
int
esl_fwrite_offset(FILE *fp, off_t *offset)
{
  if      (sizeof(off_t) == 4) return esl_fwrite_i32(fp, offset);
  else if (sizeof(off_t) == 8) return esl_fwrite_i64(fp, offset);
  else ESL_ERROR(eslESYS, "off_t is neither 32-bit nor 64-bit");
  /*UNREACHED*/
  return 1;
}









/*****************************************************************
 * @LICENSE@
 *****************************************************************/
