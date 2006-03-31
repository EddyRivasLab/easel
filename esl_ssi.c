/* esl_ssi.c
 * simple sequence indices: fast lookup in large sequence files by keyword.
 *
 *  1. Using an existing SSI index
 *  2. Creating new SSI files
 *  3. Functions for platform-independent binary i/o
 *  4. Example code
 *  5. Test driver
 * 
 * SVN $Id$
 * adapted from squid's ssi.c
 * SRE, Thu Mar  2 18:46:02 2006 [St. Louis]
 */

#include <esl_config.h>

#include <stdio.h>
#include <string.h>

#include <easel.h>
#include <esl_ssi.h>

static uint32_t v20magic = 0xf3f3e9b1; /* SSI 1.0: "ssi1" + 0x80808080 */
static uint32_t v20swap  = 0xb1e9f3f3; /* byteswapped */


/*****************************************************************
 * 1. Using an existing SSI index
 *****************************************************************/ 

static int  binary_search(ESL_SSI *ssi, char *key, uint32_t klen, off_t base, 
			  uint32_t recsize, uint32_t maxidx);

/* Function:  esl_ssi_Open()
 * Incept:    SRE, Mon Mar  6 10:52:42 2006 [St. Louis]
 *
 * Purpose:   Open the SSI index file <filename>, and returns a pointer
 *            to the new <ESL_SSI> object in <ret_ssi>.
 *            
 *            Caller is responsible for closing the SSI file with
 *            <esl_ssi_Close()>.
 *
 * Args:      <filename>   - name of SSI index file to open.       
 *            <ret_ssi>    - RETURN: the new <ESL_SSI>.
 *                        
 * Returns:   <eslOK>        on success;
 *            <eslENOTFOUND> if <filename> cannot be opened for reading;
 *            <eslEFORMAT>   if it's not in correct SSI file format;
 *            <eslERANGE>    if it uses 64-bit file offsets, and we're on a system
 *                           that doesn't support 64-bit file offsets.
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_ssi_Open(char *filename, ESL_SSI **ret_ssi)
{
  ESL_SSI *ssi = NULL;
  int      status;
  uint32_t magic;	/* magic number that starts the SSI file */
  uint16_t i;		/* counter over files */

  /* Initialize the SSI structure, null'ing so we can autocleanup.
   */
  ESL_ALLOC(ssi, sizeof(ESL_SSI));
  ssi->fp         = NULL;
  ssi->filename   = NULL;
  ssi->fileformat = NULL;
  ssi->fileflags  = NULL;
  ssi->bpl        = NULL;
  ssi->rpl        = NULL;
  ssi->nfiles     = 0;          

  /* Open the file.
   */
  status = eslENOTFOUND; 
  if ((ssi->fp = fopen(filename, "rb")) == NULL) goto CLEANEXIT; 

  /* Read the magic number: make sure it's an SSI file, and determine
   * whether it's byteswapped.
   */
  status = eslEFORMAT;
  if (esl_fread_i32(ssi->fp, &magic) != eslOK) goto CLEANEXIT;
  if (magic != v20magic && magic != v20swap)   goto CLEANEXIT;

  /* Determine what kind of offsets (32 vs. 64 bit) are stored in the file.
   * If we can't deal with 64-bit file offsets, get out now. 
   */
  status = eslEFORMAT;
  if (esl_fread_i32(ssi->fp, &(ssi->flags)) != eslOK) goto CLEANEXIT;
  ssi->imode = (ssi->flags & eslSSI_USE64_INDEX) ? 64 : 32;
  ssi->smode = (ssi->flags & eslSSI_USE64) ?       64 : 32;

  status = eslERANGE;
  if (sizeof(off_t) != 8 && (ssi->imode == 64 || ssi->smode == 64)) goto CLEANEXIT;

  /* The header data.
   */
  status = eslEFORMAT;
  if (esl_fread_i16(ssi->fp, &(ssi->nfiles))     != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->nprimary))   != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->nsecondary)) != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->flen))       != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->plen))       != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->slen))       != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->frecsize))   != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->precsize))   != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &(ssi->srecsize))   != eslOK) goto CLEANEXIT;
  
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->foffset)) != eslOK) goto CLEANEXIT;
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->poffset)) != eslOK) goto CLEANEXIT;
  if (esl_fread_offset(ssi->fp, ssi->imode, &(ssi->soffset)) != eslOK) goto CLEANEXIT;

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
  status = eslEFORMAT;
  if (ssi->nfiles == 0) goto CLEANEXIT;

  ESL_ALLOC(ssi->filename,   sizeof(char *) * ssi->nfiles);
  for (i = 0; i < ssi->nfiles; i++)  ssi->filename[i] = NULL; 
  ESL_ALLOC(ssi->fileformat, sizeof(uint32_t) * ssi->nfiles);
  ESL_ALLOC(ssi->fileflags,  sizeof(uint32_t) * ssi->nfiles);
  ESL_ALLOC(ssi->bpl,        sizeof(uint32_t) * ssi->nfiles);
  ESL_ALLOC(ssi->rpl,        sizeof(uint32_t) * ssi->nfiles);

  /* (most) mallocs done, now we read.
   */
  for (i = 0; i < ssi->nfiles; i++) 
    {
      ESL_ALLOC(ssi->filename[i], sizeof(char)* ssi->flen);

      /* We have to explicitly position, because header and file 
       * records may expand in the future; frecsize and foffset 
       * give us forwards compatibility. 
       */ 
      status = eslEFORMAT;
      if (fseeko(ssi->fp, ssi->foffset + (i * ssi->frecsize), SEEK_SET) != 0) goto CLEANEXIT;
      if (fread(ssi->filename[i],sizeof(char),ssi->flen, ssi->fp)!=ssi->flen) goto CLEANEXIT;
      if (esl_fread_i32(ssi->fp, &(ssi->fileformat[i])))                      goto CLEANEXIT;
      if (esl_fread_i32(ssi->fp, &(ssi->fileflags[i])))                       goto CLEANEXIT;
      if (esl_fread_i32(ssi->fp, &(ssi->bpl[i])))                             goto CLEANEXIT;
      if (esl_fread_i32(ssi->fp, &(ssi->rpl[i])))                             goto CLEANEXIT;
    }
  
  /* Success.
   */
  status = eslOK;			

 CLEANEXIT:
  if (status != eslOK && ssi != NULL) { esl_ssi_Close(ssi); ssi = NULL; } 
  *ret_ssi = ssi;
  return status;
}

/* Function:  esl_ssi_Close()
 * Incept:    SRE, Mon Mar  6 13:40:17 2006 [St. Louis]
 *
 * Purpose:   Close an open SSI index <ssi>.
 * 
 * Args:      <ssi>   - an open SSI index file.
 */
void
esl_ssi_Close(ESL_SSI *ssi)
{
  int i;

  if (ssi == NULL) return;

  if (ssi->fp != NULL) fclose(ssi->fp);
  if (ssi->filename != NULL) {
    for (i = 0; i < ssi->nfiles; i++) 
      if (ssi->filename[i] != NULL) free(ssi->filename[i]);
    free(ssi->filename);
  }
  if (ssi->fileformat != NULL) free(ssi->fileformat);
  if (ssi->fileflags  != NULL) free(ssi->fileflags);
  if (ssi->bpl        != NULL) free(ssi->bpl);
  if (ssi->rpl        != NULL) free(ssi->rpl);
  free(ssi);
}  


/* Function: esl_ssi_FindName()
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
 * Args:     <ssi>         - open index file
 *           <key>         - name to search for
 *           <ret_fh>      - RETURN: handle on file that key is in
 *           <ret_offset>  - RETURN: offset of the start of that key's record
 *
 * Returns:  <eslOK>        on success;
 *           <eslENOTFOUND> if no such key is in the index;
 *           <eslEFORMAT>   if an fread() or fseeko() fails, which almost
 *                          certainly reflects some kind of misformatting of
 *                          the index.
 *
 * Throws:   <eslEMEM>      on allocation error.
 */
int
esl_ssi_FindName(ESL_SSI *ssi, char *key, uint16_t *ret_fh, off_t *ret_offset)
{
  int       status;
  char     *pkey   = NULL;
  uint16_t  fh     = 0;
  off_t     offset = 0;

  /* Look in the primary keys.
   */
  status = binary_search(ssi, key, ssi->plen, ssi->poffset, ssi->precsize,
			 ssi->nprimary);

  if (status == eslOK) 
    {		
      /* We found it as a primary key; get our data & return.
       */
      status = eslEFORMAT;
      if (esl_fread_i16(ssi->fp, &fh) != eslOK)                    goto CLEANEXIT;
      if (esl_fread_offset(ssi->fp, ssi->smode, &offset) != eslOK) goto CLEANEXIT;
    } 
  else if (status == eslENOTFOUND) 
    {
      /* Not in the primary keys? OK, try the secondary keys.
       */
      if (ssi->nsecondary > 0) {
	status = binary_search(ssi, key, ssi->slen, ssi->soffset, ssi->srecsize,
			       ssi->nsecondary);
	if (status != eslOK) goto CLEANEXIT;

      /* We have the secondary key; flip to its primary key, then look that up.
       */
	ESL_ALLOC(pkey, sizeof(char) * ssi->plen);
	status = eslEFORMAT;
	if (fread(pkey, sizeof(char), ssi->plen, ssi->fp) != ssi->plen) goto CLEANEXIT;
	status = esl_ssi_FindName(ssi, pkey, &fh, &offset);
	if (status != eslOK) goto CLEANEXIT;
      }
    }
  else goto CLEANEXIT;	/* status from binary search was an error code. */

  status = eslOK;

 CLEANEXIT:
  if (pkey != NULL) free(pkey);

  if (status == eslOK) { 
    *ret_fh     = fh;
    *ret_offset = offset;
  } else {
    *ret_fh     = 0;
    *ret_offset = 0;
  }
  return status;
}



/* Function:  esl_ssi_FindNumber()
 * Incept:    SRE, Mon Jan  1 19:42:42 2001 [St. Louis]
 *
 * Purpose:   Looks up primary key number <nkey> in the open index
 *            <ssi>.  <nkey> ranges from <0..ssi->nprimary-1>. When
 *            key <nkey> is found, <ret_fh> contains a unique handle
 *            on the file that contains that key (suitable for an
 *            <esl_ssi_FileInfo()> call, or for comparison to the handle of
 *            the last file that was opened for retrieval), and
 *            <ret_offset> is filled in with the offset in that file.
 *           
 * Args:      <ssi>        - open index file
 *            <n>          - primary key number to retrieve.
 *            <ret_fh>     - RETURN: handle on file that key is in
 *            <ret_offset> - RETURN: offset of the start of that key's record
 *
 * Returns:   <eslOK>        on success;
 *            <eslENOTFOUND> if there is no sequence record <nkey>;
 *            <eslEFORMAT>   if a read or a seek fails, probably indicating
 *                           some kind of file misformatting.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_ssi_FindNumber(ESL_SSI *ssi, int nkey, uint16_t *ret_fh, off_t *ret_offset)
{
  int      status;
  uint16_t fh;
  off_t    offset;
  char    *pkey = NULL;


  if (nkey >= ssi->nprimary) { status = eslENOTFOUND; goto CLEANEXIT; }
  ESL_ALLOC(pkey, sizeof(char) * ssi->plen);

  status = eslEFORMAT;
  if (fseeko(ssi->fp, ssi->poffset+ssi->precsize*nkey, SEEK_SET)!= 0)         goto CLEANEXIT;
  if (fread(pkey, sizeof(char), ssi->plen, ssi->fp)             != ssi->plen) goto CLEANEXIT;
  if (esl_fread_i16(ssi->fp, &fh)                               != eslOK)     goto CLEANEXIT;
  if (esl_fread_offset(ssi->fp, ssi->smode, ret_offset)         != eslOK)     goto CLEANEXIT;

  status = eslOK;

 CLEANEXIT:
  if (pkey != NULL) free(pkey);

  if (status == eslOK) {
    *ret_fh     = fh;
    *ret_offset = offset;
  } else {
    *ret_fh     = 0;
    *ret_offset = 0;
  }
  return status;
}


/* Function: esl_ssi_FindSubseq()
 * Date:     SRE, Mon Jan  1 19:49:31 2001 [St. Louis]
 *
 * Purpose:  Fast subsequence retrieval:
 *           look up a primary or secondary <key> in the open
 *           index <ssi>. Ask for the nearest data offset to a
 *           subsequence starting at residue <requested_start>
 *           in the sequence (numbering the sequence <1..L>). 
 *           If <key> is found, on return, <ret_fh>
 *           contains a unique handle on the file that contains 
 *           <key>; <record_offset> contains the
 *           disk offset to the start of the sequence record; <data_offset>
 *           contains the disk offset either exactly at the requested
 *           residue, or at the start of the line containing the
 *           requested residue; <ret_actual_start> contains the 
 *           coordinate (1..L) of the first valid residue at or
 *           after <data_offset>. <ret_actual_start> is $\leq$ 
 *           <requested_start>. 
 *
 * Args:     <ssi>             - open index file
 *           <key>             - primary or secondary key to find
 *           <requested_start> - residue we'd like to start at (1..L)
 *           <ret_fh>          - RETURN: handle for file the key is in
 *           <record_offset>   - RETURN: offset to start of sequence record
 *           <data_offset>     - RETURN: offset to start of subseq (see above)
 *           <ret_actual_start>- RETURN: coord (1..L) of residue at data_offset
 *
 * Returns:  <eslOK>      on success;
 *           <eslEINVAL>  if we can't do fast subsequence lookup on the file;
 *           <eslEFORMAT> on a read or seek failure, presumably meaning that
 *                        the file is misformatted somehow;
 *           <eslERANGE>  if <requested_start> isn't somewhere in the range
 *                        <1..len> for the target sequence.
 */
int
esl_ssi_FindSubseq(ESL_SSI *ssi, char *key, long requested_start,
		   uint16_t *ret_fh, off_t *record_offset, off_t *data_offset, 
		   long *ret_actual_start)
{
  int      status;
  uint16_t fh;
  off_t    r_off, d_off;
  long     actual_start;
  uint32_t len;
  int      r, b, i, l;	/* tmp variables for "clarity", to match docs */
  
  /* Look up the key. Rely on the fact that esl_ssi_FindName()
   * leaves the index file positioned at the rest of the data for this key.
   */
  status = esl_ssi_FindName(ssi, key, &fh, &r_off);
  if (status != eslOK) goto CLEANEXIT;

  /* Check that we're allowed to do subseq lookup on that file.
   */
  if (! (ssi->fileflags[fh] & eslSSI_FASTSUBSEQ))
    { status = eslEINVAL; goto CLEANEXIT; }

  /* Read the rest of the index record for this primary key:
   * the data offset, and seq length.
   */
  status = eslEFORMAT;
  if (esl_fread_offset(ssi->fp, ssi->smode, &d_off) != eslOK) goto CLEANEXIT;
  if (esl_fread_i32(ssi->fp, &len)                  != eslOK) goto CLEANEXIT;

  /* Set up tmp variables for clarity of equations below,
   * and to make them match tex documentation 
   */
  r = ssi->rpl[fh];         /* residues per line */
  b = ssi->bpl[fh];         /* bytes per line    */
  i = requested_start;	    /* start position 1..L */
  l = (i-1)/r;		    /* data line # (0..) that the residue is on */
  if (r == 0 || b == 0) { status = eslEINVAL; goto CLEANEXIT; }
  if (i < 0 || i > len) { status = eslERANGE; goto CLEANEXIT; }
  
  /* When b = r+1, there's nothing but sequence on each data line (and the \0).
   * In this case, we know we can find each residue precisely.
   */
  if (b == r+1) 
    {
      d_off       += l*b + (i-1)%r;
      actual_start = requested_start;
    } 
  /* else, there's other stuff on seq lines - probably spaces - so the best
   * we can do (without figuring out the spacing pattern and checking that
   * it's consistent everywhere) is to position at start of relevant line.
   */
  else
    { 
      d_off       += l*b;
      actual_start = 1 + l*r;
    }

  status = eslOK;

 CLEANEXIT:
  if (status == eslOK) {
    *ret_fh           = fh;
    *record_offset    = r_off;
    *data_offset      = d_off;
    *ret_actual_start = actual_start;
  } else {
    *ret_fh           = 0;
    *record_offset    = 0;
    *data_offset      = 0;
    *ret_actual_start = 0;
  }
  return status;
}


/* Function: esl_ssi_FileInfo()
 * Date:     SRE, Tue Jan  2 10:31:01 2001 [St. Louis]
 *
 * Purpose:  Given a file number <fh> in an open index file
 *           <ssi>, retrieve file name <ret_filename> and
 *           the file format <ret_format>. 
 *           
 *           <ret_filename> is a pointer to a string maintained
 *           internally by <ssi>. It should not be free'd; 
 *           <esl_ssi_Close(ssi)> will take care of it.
 *
 * Args:     <ssi>          - open index file
 *           <fh>           - handle on file to look up
 *           <ret_filename> - RETURN: name of file n
 *           <ret_format>   - RETURN: format code for file n
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEINVAL> if there is no such file number <fh>.
 */
int
esl_ssi_FileInfo(ESL_SSI *ssi, uint16_t fh, char **ret_filename, int *ret_format)
{
  int status;

  if (fh >= ssi->nfiles) ESL_DIE(eslEINVAL, "no such file number");
  *ret_filename = ssi->filename[fh];
  *ret_format   = ssi->fileformat[fh];
  status = eslOK;

 CLEANEXIT:
  return status;
}


/* binary_search()
 * Date:     SRE, Sun Dec 31 16:05:03 2000 [St. Louis]
 *
 * Purpose:  Find <key> in an SSI index, by a binary search
 *           in an alphabetically sorted list of keys. If successful,
 *           return <eslOK>, and the index file is positioned to read
 *           the rest of the data for that key. If unsuccessful, 
 *           return <eslFAIL>, and the positioning of the index file
 *           is left in an undefined state.
 *
 * Args:     <ssi>     - an open ESL_SSI
 *           <key>     - key to find
 *           <klen>    - key length to allocate (plen or slen from ssi)
 *           <base>    - base offset (poffset or soffset)
 *           <recsize> - size of each key record in bytes (precsize or srecsize)
 *           <maxidx>  - # of keys (nprimary or nsecondary)
 *
 * Returns:  <eslOK> on success, and leaves file positioned for reading remaining
 *           data for the key. 
 *           
 *           <eslENOTFOUND> if <key> is not found.
 *           <eslEFORMAT>   if an fread() or fseeko() fails, probably indicating
 *                          some kind of misformatting of the index file.
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

  ESL_ALLOC(name, (sizeof(char)*klen));

  left  = 0;
  right = maxidx-1;
  while (1) {			/* A binary search: */
    mid   = (left+right) / 2;	/* careful here. left+right potentially overflows if
				   we didn't limit unsigned vars to signed ranges. */
    status = eslEFORMAT;
    if (fseeko(ssi->fp, base + recsize*mid, SEEK_SET) != 0)    goto CLEANEXIT;
    if (fread(name, sizeof(char), klen, ssi->fp)      != klen) goto CLEANEXIT;

    status = eslENOTFOUND;
    cmp = strcmp(name, key);
    if      (cmp == 0) break;	             /* found it!               */
    else if (left >= right) goto CLEANEXIT;  /* no such key             */
    else if (cmp < 0)       left  = mid+1;   /* it's still right of mid */
    else if (cmp > 0) {
      if (mid == 0) goto CLEANEXIT;          /* beware left edge case   */
      else right = mid-1;                    /* it's left of mid        */
    }
  }
  status = eslOK;

 CLEANEXIT:
  if (name != NULL) free(name);
  return status; /* and if eslOK, ssi->fp is positioned to read the record. */
}


/*****************************************************************
 * 2. Creating new SSI files
 *****************************************************************/ 
static int current_newssi_size(ESL_NEWSSI *ns);
static int activate_external_sort(ESL_NEWSSI *ns);
static int parse_pkey(char *buf, ESL_PKEY *pkey);
static int parse_skey(char *buf, ESL_SKEY *skey);
static int pkeysort(const void *k1, const void *k2);
static int skeysort(const void *k1, const void *k2);

/* Function:  esl_newssi_Create()
 * Incept:    SRE, Tue Jan  2 11:23:25 2001 [St. Louis]
 *
 * Purpose:   Creates and returns a <ESL_NEWSSI>, in order to create a 
 *            new SSI index file.
 *
 * Returns:   a pointer to the <ESL_NEWSSI> structure.
 *            
 * Throws:    <NULL> on allocation error.
 */
ESL_NEWSSI *
esl_newssi_Create(void)
{
  int status;
  ESL_NEWSSI *ns = NULL;

  ESL_ALLOC(ns, sizeof(ESL_NEWSSI));

  ns->external   = FALSE;	    /* we'll switch to external sort if...       */
  ns->max_ram    = eslSSI_MAXRAM;   /* ... if we exceed this memory limit in MB. */
  ns->filenames  = NULL;
  ns->fileformat = NULL;
  ns->bpl        = NULL;
  ns->rpl        = NULL;
  ns->flen       = 0;
  ns->nfiles     = 0;
  ns->pkeys      = NULL;
  ns->plen       = 0;
  ns->nprimary   = 0;
  ns->ptmpfile   = ".ssi.tmp.1"; /* hardcoded, for now */
  ns->ptmp       = NULL;
  ns->skeys      = NULL;
  ns->slen       = 0;
  ns->nsecondary = 0;
  ns->stmpfile   = ".ssi.tmp.2"; /* hardcoded, for now */
  ns->stmp       = NULL;

  /* All mallocs follow NULL initializations, because of the cleanup strategy:
   * we free anything non-NULL if any malloc fails.
   * filenames[] array doesn't have to be NULL'ed because we're tracking nfiles.
   */
  ESL_ALLOC(ns->filenames,  sizeof(char *)   * eslSSI_FCHUNK);
  ESL_ALLOC(ns->fileformat, sizeof(uint32_t) * eslSSI_FCHUNK);
  ESL_ALLOC(ns->bpl,        sizeof(uint32_t) * eslSSI_FCHUNK);
  ESL_ALLOC(ns->rpl,        sizeof(uint32_t) * eslSSI_FCHUNK);
  ESL_ALLOC(ns->pkeys,      sizeof(ESL_PKEY) * eslSSI_KCHUNK);
  ESL_ALLOC(ns->skeys,      sizeof(ESL_SKEY) * eslSSI_KCHUNK);
  status = eslOK;

 CLEANEXIT:
  if (status == eslOK)  return ns;
  else {
    esl_newssi_Destroy(ns);	/* free the damaged structure */
    return NULL;
  } 
}


/* Function:  esl_newssi_AddFile()
 * Incept:    SRE, Tue Mar  7 08:57:39 2006 [St. Louis]
 *
 * Purpose:   Registers the file <filename> into the new index <ns>,
 *            along with its format code <fmt>. The index assigns it
 *            a unique handle, which it returns in <ret_fh>. This
 *            handle is needed when registering primary keys.
 *
 *            Caller should make sure that the same file isn't registered
 *            twice; this function doesn't check.
 *            
 * Args:      <ns>         - new ssi index under construction.
 *            <filename>   - filename to add to the index.
 *            <fmt>        - format code to associate with <filename> (or 0)
 *            <ret_fh>     - RETURN: filehandle associated with <filename>        
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if registering this file would exceed the
 *            maximum number of indexed files.
 *
 * Throws:    <eslEMEM> on allocation or reallocation error.
 */
int
esl_newssi_AddFile(ESL_NEWSSI *ns, char *filename, int fmt, uint16_t *ret_fh)
{
  int      status;
  uint16_t fh;
  int      n;

  if (ns->nfiles >= eslSSI_MAXFILES) { status = eslERANGE; goto CLEANEXIT; }

  n = strlen(filename);
  if ((n+1) > ns->flen) ns->flen = n+1;

  status = esl_FileTail(filename, FALSE, &(ns->filenames[ns->nfiles]));
  if (status != eslOK) ESL_DIE(status, "esl_FileTail() failed");
  
  ns->fileformat[ns->nfiles] = fmt;
  ns->bpl[ns->nfiles]        = 0;
  ns->rpl[ns->nfiles]        = 0;
  fh                         = ns->nfiles;   /* handle is simply = file number */
  ns->nfiles++;

  if (ns->nfiles % eslSSI_FCHUNK == 0) {
    void  *tmp;
    ESL_RALLOC(ns->filenames,  tmp, sizeof(char *)   * (ns->nfiles+eslSSI_FCHUNK));
    ESL_RALLOC(ns->fileformat, tmp, sizeof(uint32_t) * (ns->nfiles+eslSSI_FCHUNK));
    ESL_RALLOC(ns->bpl,        tmp, sizeof(uint32_t) * (ns->nfiles+eslSSI_FCHUNK));
    ESL_RALLOC(ns->rpl,        tmp, sizeof(uint32_t) * (ns->nfiles+eslSSI_FCHUNK));
  }
  status = eslOK;

 CLEANEXIT:
  if (status == eslOK) *ret_fh = fh;
  else                 *ret_fh = 0;
  return status;
}



/* Function:  esl_newssi_SetSubseq()
 * Incept:    SRE, Tue Mar  7 09:03:59 2006 [St. Louis]
 *
 * Purpose:   Declare that the file associated with handle <fh> is
 *            suitable for fast subsequence lookup, because it has
 *            a constant number of residues and bytes per (nonterminal)
 *            data line, <rpl> and <bpl>, respectively.
 *            
 *            Caller is responsible for this being true: <rpl> and
 *            <bpl> must be constant for every nonterminal line of 
 *            every sequence in this file.
 *            
 * Args:      <ns>   - ssi index under construction
 *            <fh>   - handle on file to set fast subseq lookup on
 *            <bpl>  - constant bytes per nonterminal line in <fh>                   
 *            <rpl>  - constant residues per nonterminal line in <fh>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> on invalid argument(s).
 */
int
esl_newssi_SetSubseq(ESL_NEWSSI *ns, uint16_t fh, int bpl, int rpl)
{
  int status;

  if (fh >= ns->nfiles)      ESL_DIE(eslEINVAL, "invalid file number");
  if (bpl <= 0 || rpl <= 0)  ESL_DIE(eslEINVAL, "invalid bpl or rpl");
  ns->bpl[fh] = bpl;
  ns->rpl[fh] = rpl;
  status = eslOK;

 CLEANEXIT:
  return status;
}


/* Function: esl_newssi_AddKey()
 * Date:     SRE, Tue Jan  2 11:50:54 2001 [St. Louis]
 *
 * Purpose:  Register primary key <key> in new index <ns>, while telling
 *           the index that this primary key is in the file associated
 *           with filehandle <fh> (the handle returned by a previous call
 *           to <esl_newssi_AddFile()>), and that its record starts at 
 *           offset <r_off> in the file.
 *           
 *           <d_off> and <L> are optional; they may be left unset by
 *           passing <NULL> and <0>, respectively. (If one is
 *           provided, both must be provided.) If they are provided,
 *           <d_off> gives the position of the first line of sequence
 *           data in the record, and <L> gives the length of the
 *           sequence in residues. These are necessary when
 *           <eslSSI_FASTSUBSEQ> is set for this file, for fast
 *           subsequence lookup. If <eslSSI_FASTSUBSEQ> is not set for
 *           the file, <d_off> and <L> will be ignored by the index
 *           reading API, so it doesn't hurt to provide them;
 *           typically an indexing program won't know whether it's safe to set
 *           <eslSSI_FASTSUBSEQ> for the whole file until the whole file
 *           has been read and every key has already been added to the
 *           index.
 *           
 * Args:     <ns>     - active index
 *           <key>    - primary key to add
 *           <fh>     - handle on file that this key's in 
 *           <r_off>  - offset to start of record
 *           <d_off>  - offset to start of sequence data, or 0
 *           <L>      - length of sequence, or 0
 *
 * Returns:  <eslOK>        on success;
 *           <eslERANGE>    if registering this key would exceed the maximum
 *                          number of primary keys;
 *           <eslEDUP>      if we needed to open external tmp files for a 
 *                          large index, but they already existed;
 *           <eslENOTFOUND> if we needed to open external tmp files, but
 *                          the attempt to open them failed.
 *           
 * Throws:   <eslEINVAL> on an invalid argument;
 *           <eslEMEM>   on allocation failure.       
 */
int
esl_newssi_AddKey(ESL_NEWSSI *ns, char *key, uint16_t fh, 
		  off_t r_off, off_t d_off, uint32_t L)
{
  int status;
  int n;			/* a string length */
  
  if (fh >= eslSSI_MAXFILES)           ESL_DIE(eslEINVAL, "invalid fh");
  if (ns->nprimary >= eslSSI_MAXKEYS)  return eslERANGE;

  /* Before adding the key: check how big our index is.
   * If it's getting too large, switch to external mode.
   */
  if (!ns->external && current_newssi_size(ns) >= ns->max_ram) 
    if ((status = activate_external_sort(ns)) != eslOK) goto CLEANEXIT;

  /* Update maximum pkey length, if needed. (Inclusive of '\0').
   */
  n = strlen(key)+1;
  if (n > ns->plen) ns->plen = n;

  /* External mode? Simply append to disk... 
   */
  if (ns->external) 
    {
      if (sizeof(off_t) == 4) {
	fprintf(ns->ptmp, "%s\t%d\t%lu\t%lu\t%lu\n", 
		key, 
		fh,
		(unsigned long) r_off, 
		(unsigned long) d_off,
		(unsigned long) L);
      } else {
	fprintf(ns->ptmp, "%s\t%d\t%llu\t%llu\t%lu\n", 
		key, 
		fh, 
		(unsigned long long) r_off,
		(unsigned long long) d_off,
		(unsigned long) L);
      }
      ns->nprimary++;
    }
  else
    {
      /* Else: internal mode, keep keys in memory...
       */
      if (esl_strdup(key, n, &(ns->pkeys[ns->nprimary].key)) != eslOK)
	ESL_DIE(eslEMEM, "esl_strdup failed");
      ns->pkeys[ns->nprimary].fnum  = fh;
      ns->pkeys[ns->nprimary].r_off = r_off;
      ns->pkeys[ns->nprimary].d_off = d_off;
      ns->pkeys[ns->nprimary].len   = L;
      ns->nprimary++;

      /* Reallocate as needed.
       */
      if (ns->nprimary % eslSSI_KCHUNK == 0) {
	void *tmp;
	ESL_RALLOC(ns->pkeys, tmp, sizeof(ESL_PKEY) * (ns->nprimary+eslSSI_KCHUNK));
      }
    }
  status = eslOK;

 CLEANEXIT:
  return status;
}

/* Function:  esl_newssi_AddAlias()
 * Incept:    SRE, Tue Mar  7 15:49:43 2006 [St. Louis]
 *
 * Purpose:   Registers secondary key <alias> in index <ns>, and 
 *            map it to the primary key <key>. <key> must already
 *            have been registered. That is, when someone looks up <alias>,
 *            we'll retrieve record <key>. 
 *            
 * Args:      <ns>    - ssi index being constructed
 *            <alias> - secondary key to register
 *            <key>   - primary key to associate with <skey>.                  
 *
 * Returns:   <eslOK>        on success;
 *            <eslERANGE>    if registering this key would exceed the maximum
 *                           number of secondary keys that can be stored;
 *            <eslEDUP>      if we needed to open external tmp files for a 
 *                           large index, but they already existed;
 *            <eslENOTFOUND> if we needed to open external tmp files, but
 *                           the attempt to open them failed.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_newssi_AddAlias(ESL_NEWSSI *ns, char *alias, char *key)
{
  int status;
  int n;			/* a string length */
  
  if (ns->nsecondary >= eslSSI_MAXKEYS) return eslERANGE;

  /* Before adding the key: check how big our index is.
   * If it's getting too large, switch to external mode.
   */
  if (!ns->external && current_newssi_size(ns) >= ns->max_ram) 
    if ((status = activate_external_sort(ns)) != eslOK) goto CLEANEXIT;

  /* Update maximum secondary key length, if necessary.
   */
  n = strlen(alias)+1;
  if (n > ns->slen) ns->slen = n;

  /* if external mode: write info to disk.
   */
  if (ns->external) 
    {
      fprintf(ns->stmp, "%s\t%s\n", alias, key);
      ns->nsecondary++;
    }
  else
    {
      /* else, internal mode... store info in memory.
       */
      if ((status = esl_strdup(alias, n, &(ns->skeys[ns->nsecondary].key))) != eslOK) goto CLEANEXIT;
      if ((status = esl_strdup(key, -1, &(ns->skeys[ns->nsecondary].pkey))) != eslOK) goto CLEANEXIT;
      ns->nsecondary++;

      if (ns->nsecondary % eslSSI_KCHUNK == 0) {
	void *tmp;
	ESL_RALLOC(ns->skeys, tmp, sizeof(ESL_SKEY) * (ns->nsecondary+eslSSI_KCHUNK));
      }
    }

  status = eslOK;
 CLEANEXIT:
  return status;
}


/* Function:  esl_newssi_Write()
 * Incept:    SRE, Tue Mar  7 16:06:09 2006 [St. Louis]
 *
 * Purpose:   Writes the complete index <ns> in SSI format to a binary
 *            stream <fp>, which the caller has already opened.
 *            
 *            Handles all necessary overhead of sorting the primary and
 *            secondary keys, including any externally sorted tmpfiles that
 *            may have been needed for large indices.
 *            
 * Args:      <fp>  - open file stream to write the index to
 *            <ns>  - new SSI index to write                   
 *            
 * Returns:   <eslOK>    on success;
 *            <eslEFAIL> if any write fails, or if index 
 *                       size exceeds system's maximum file size;
 *            <eslESYS>  if any of the steps of an external sort fail.
 *
 * Throws:    <eslEINVAL> on invalid argument, including too-long tmpfile names;
 *            <eslEMEM>   on buffer allocation failure.
 */
int
esl_newssi_Write(FILE *fp, ESL_NEWSSI *ns)
{
  int      status, 		/* convention                               */
           i;			/* counter over files, keys                 */
  uint32_t header_flags,	/* bitflags in the header                   */
           file_flags,		/* bitflags for a file record               */
           frecsize, 		/* size of a file record (bytes)            */
           precsize, 		/* size of a primary key record (bytes)     */
           srecsize;		/* size of a secondary key record (bytes)   */
  off_t    foffset, 		/* offset to file section                   */
           poffset, 		/* offset to primary key section            */
           soffset;		/* offset to secondary key section          */
  char    *fk       = NULL,     /* fixed-width (flen) file name             */
          *pk       = NULL, 	/* fixed-width (plen) primary key string    */
          *sk       = NULL,	/* fixed-width (slen) secondary key string  */
          *buf      = NULL;	/* esl_fgets() growable buffer              */
  int      n        = 0;	/* esl_fgets() buffer size                  */
  ESL_PKEY pkey;		/* primary key info from external tmpfile   */
  ESL_SKEY skey;		/* secondary key info from external tmpfile */

  /* We need fixed-width buffers to get our keys fwrite()'ten in their
   * full binary lengths; pkey->key (for instance) is not guaranteed
   * to be allocated for the final maximum plen.
   */
  ESL_ALLOC(fk, sizeof(char) * ns->flen);
  ESL_ALLOC(pk, sizeof(char) * ns->plen);
  ESL_ALLOC(sk, sizeof(char) * ns->slen);

  /* How big is the index? If it's going to be > 2GB, we better have
   * 64-bit offsets. (2047 (instead of 2048) gives us
   * some slop room.) If not, abort here.
   *
   * aborting here is pretty brutal - we've processed hundreds of
   * millions of keys for nothing. Ah well.
   */
  if (current_newssi_size(ns) >= 2047 && sizeof(off_t) != 8)
    { status = eslFAIL; goto CLEANEXIT; }

  /* Magic-looking numbers come from adding up sizes 
   * of things in bytes: they match current_newssi_size().
   */
  frecsize = 16 + ns->flen;
  precsize = 2*sizeof(off_t) + 6 + ns->plen;
  srecsize = ns->slen + ns->plen;

  header_flags = 0;
  if (sizeof(off_t) == 8) 
    {
      header_flags |= eslSSI_USE64;
      header_flags |= eslSSI_USE64_INDEX;
    }

  /* Magic-looking numbers again come from adding up sizes 
   * of things in bytes: matches current_newssi_size()
   */
  foffset = 3*sizeof(off_t) + 42; /* the answer: of course */
  poffset = foffset + frecsize*ns->nfiles;
  soffset = poffset + precsize*ns->nprimary;
  
  /* Sort the keys.
   * If external mode, make system calls to UNIX/POSIX "sort" in place, then
   * open new sorted files for reading thru ptmp and stmp handles.
   * If internal mode, call qsort. 
   * 
   * Note that you'd better force a POSIX locale for the sort; else,
   * some silly distro (e.g. Mandrake Linux >=8.1) may have specified
   * LC_COLLATE=en_US, and this'll give a sort "bug" in which it doesn't
   * sort by byte order.
   */
  if (ns->external) 
    {
      char cmd[1024];

      /* A last minute security check: make sure we won't overflow
       * sprintf() with the tmpfile names. They're hardcoded now, so
       * we know they don't overflow, but they might be configurable 
       * in the future, and we wouldn't want a security hole to open
       * up.
       */
      if (strlen(ns->ptmpfile) > 256 || strlen(ns->ptmpfile) > 256) 
	ESL_DIE(eslEINVAL, "tmpfile name too long"); 

      status = eslESYS;	/* any premature return now is ESYS error */
      fclose(ns->ptmp);
      ns->ptmp = NULL;	
      sprintf(cmd, "env LC_ALL=POSIX sort -o %s %s\n", ns->ptmpfile, ns->ptmpfile);
      if (system(cmd) != 0)                              goto CLEANEXIT;
      if ((ns->ptmp = fopen(ns->ptmpfile, "r")) == NULL) goto CLEANEXIT;

      fclose(ns->stmp);
      ns->stmp = NULL;
      sprintf(cmd, "env LC_ALL=POSIX sort -o %s %s\n", ns->stmpfile, ns->stmpfile);
      if (system(cmd) != 0)                              goto CLEANEXIT;
      if ((ns->stmp = fopen(ns->stmpfile, "r")) == NULL) goto CLEANEXIT;
    }
  else 
    {
      qsort((void *) ns->pkeys, ns->nprimary,   sizeof(ESL_PKEY), pkeysort); 
      qsort((void *) ns->skeys, ns->nsecondary, sizeof(ESL_SKEY), skeysort); 
    }

  /* Write the header
   */
  status = eslFAIL;		/* any write error is a FAIL */
  if (esl_fwrite_i32(fp, v20magic)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, header_flags)  != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i16(fp, ns->nfiles)    != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, ns->nprimary)  != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, ns->nsecondary)!= eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, ns->flen)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, ns->plen)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, ns->slen)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, frecsize)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, precsize)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_i32(fp, srecsize)      != eslOK) goto CLEANEXIT;
  if (esl_fwrite_offset(fp, foffset)    != eslOK) goto CLEANEXIT;
  if (esl_fwrite_offset(fp, poffset)    != eslOK) goto CLEANEXIT;
  if (esl_fwrite_offset(fp, soffset)    != eslOK) goto CLEANEXIT;

  /* Write the file section
   */
  for (i = 0; i < ns->nfiles; i++)
    {
      file_flags = 0;
      if (ns->bpl[i] > 0 && ns->rpl[i] > 0) file_flags |= eslSSI_FASTSUBSEQ;
      strcpy(fk, ns->filenames[i]);

      status     = eslFAIL;
      if (fwrite(fk, sizeof(char), ns->flen, fp) != ns->flen) goto CLEANEXIT;
      if (esl_fwrite_i32(fp, ns->fileformat[i])  != eslOK)    goto CLEANEXIT;              
      if (esl_fwrite_i32(fp, file_flags)         != eslOK)    goto CLEANEXIT;             
      if (esl_fwrite_i32(fp, ns->bpl[i])         != eslOK)    goto CLEANEXIT;
      if (esl_fwrite_i32(fp, ns->rpl[i])         != eslOK)    goto CLEANEXIT;
    }

  /* Write the primary key section
   */
  if (ns->external) 
    {
      for (i = 0; i < ns->nprimary; i++) 
	{
	  status = eslESYS;		/* any external read error is an ESYS */
	  if (esl_fgets(&buf, &n, ns->ptmp)  != eslOK)    goto CLEANEXIT;
	  if (parse_pkey(buf, &pkey)         != eslOK)    goto CLEANEXIT;
	  strcpy(pk, pkey.key);

	  status = eslFAIL;		/* any write error is an EFAIL */
	  if (fwrite(pk,sizeof(char),ns->plen,fp) != ns->plen) goto CLEANEXIT;
	  if (esl_fwrite_i16(   fp, pkey.fnum)    != eslOK)    goto CLEANEXIT;   
	  if (esl_fwrite_offset(fp, pkey.r_off)   != eslOK)    goto CLEANEXIT;
	  if (esl_fwrite_offset(fp, pkey.d_off)   != eslOK)    goto CLEANEXIT;
	  if (esl_fwrite_i32(   fp, pkey.len)     != eslOK)    goto CLEANEXIT;
	}
    } 
  else 
    {
      for (i = 0; i < ns->nprimary; i++)
	{
	  strcpy(pk, ns->pkeys[i].key);

	  status = eslFAIL;
	  if (fwrite(pk,sizeof(char),ns->plen,fp)      != ns->plen) goto CLEANEXIT;
	  if (esl_fwrite_i16(   fp, ns->pkeys[i].fnum)  != eslOK)    goto CLEANEXIT;
	  if (esl_fwrite_offset(fp, ns->pkeys[i].r_off) != eslOK)    goto CLEANEXIT;
	  if (esl_fwrite_offset(fp, ns->pkeys[i].d_off) != eslOK)    goto CLEANEXIT;
	  if (esl_fwrite_i32(   fp, ns->pkeys[i].len)   != eslOK)    goto CLEANEXIT;
	}
    }

  /* Write the secondary key section
   */
  if (ns->external) 
    {
      for (i = 0; i < ns->nsecondary; i++)
	{
	  status = eslESYS;
	  if (esl_fgets(&buf, &n, ns->stmp) != eslOK) goto CLEANEXIT;
	  if (parse_skey(buf, &skey)        != eslOK) goto CLEANEXIT;
	  strcpy(sk, skey.key);
	  strcpy(pk, skey.pkey);

	  status = eslFAIL;
	  if (fwrite(sk, sizeof(char), ns->slen, fp) != ns->slen) goto CLEANEXIT;
	  if (fwrite(pk, sizeof(char), ns->plen, fp) != ns->plen) goto CLEANEXIT;
	}
    } 
  else 
    {
      for (i = 0; i < ns->nsecondary; i++)
	{
	  strcpy(sk, ns->skeys[i].key);
	  strcpy(pk, ns->skeys[i].pkey);

	  status = eslFAIL;
	  if (fwrite(sk, sizeof(char), ns->slen, fp) != ns->slen) goto CLEANEXIT;
	  if (fwrite(pk, sizeof(char), ns->plen, fp) != ns->plen) goto CLEANEXIT;
	} 
    }

  status = eslOK;

 CLEANEXIT: 
  if (fk  != NULL)       free(fk);
  if (pk  != NULL)       free(pk);
  if (sk  != NULL)       free(sk);
  if (buf != NULL)       free(buf);
  if (ns->ptmp != NULL)  { fclose(ns->ptmp); ns->ptmp = NULL; }
  if (ns->stmp != NULL)  { fclose(ns->stmp); ns->stmp = NULL; }
  return status;
}

/* Function:  esl_newssi_Destroy()
 * Incept:    SRE, Tue Mar  7 08:13:27 2006 [St. Louis]
 *
 * Purpose:   Frees a <ESL_NEWSSI>.
 */
void
esl_newssi_Destroy(ESL_NEWSSI *ns)
{
  int i;
  if (ns == NULL) return;

  if (ns->external == FALSE) 
    {
      if (ns->pkeys != NULL) 
	{
	  for (i = 0; i < ns->nprimary; i++) 
	    if (ns->pkeys[i].key != NULL) free(ns->pkeys[i].key);
	  free(ns->pkeys);       	
	}
      if (ns->skeys != NULL) 
	{
	  for (i = 0; i < ns->nsecondary; i++) 
	    {
	      if (ns->skeys[i].key  != NULL) free(ns->skeys[i].key);
	      if (ns->skeys[i].pkey != NULL) free(ns->skeys[i].pkey);
	    }
	  free(ns->skeys);       
	}
    }
  else 
    {
      if (ns->ptmp != NULL) fclose(ns->ptmp);
      if (ns->stmp != NULL) fclose(ns->stmp);       
    }

  if (ns->filenames   != NULL)  
    {
      for (i = 0; i < ns->nfiles; i++) 
	if (ns->filenames[i] != NULL) free(ns->filenames[i]);
      free(ns->filenames);
    }
  if (ns->fileformat  != NULL)     free(ns->fileformat);
  if (ns->bpl         != NULL)     free(ns->bpl);       
  if (ns->rpl         != NULL)     free(ns->rpl);       
  free(ns);
}




/* current_newssi_size()
 *
 * Calculates the size of the current index, in megabytes, in
 * its disk version (which is essentially the same as the
 * RAM it takes, modulo some small overhead for the structures
 * and ptrs).
 *  
 * The header costs 10 uint32, 1 uint16, and 3 off_t's: 42 + (12 | 24).
 * Each file record costs 4 uint32 and flen chars;
 * each primary key costs us 2 off_t, 1 uint16, 1 uint32, and plen chars;
 * each sec key costs us  plen+slen chars.
 */
static int
current_newssi_size(ESL_NEWSSI *ns) 
{
  uint64_t frecsize, precsize, srecsize;
  uint64_t total;

  /* Magic-looking numbers come from adding up sizes 
   * of things in bytes
   */
  frecsize = 16 + ns->flen;
  precsize = 2*sizeof(off_t) + 6 + ns->plen;
  srecsize = ns->plen + ns->slen;
  total = (42 +		               /* header size, if 64bit index offsets */
	   3 * sizeof(off_t) + 
	   frecsize * ns->nfiles +      /* file section size                   */
	   precsize * ns->nprimary +    /* primary key section size            */
	   srecsize * ns->nsecondary) / /* secondary key section size          */
          1048576L;
  return (int) total;
}

/* activate_external_sort()
 * 
 * Switch to external sort mode.
 * Open file handles for external index files (ptmp, stmp).
 * Flush current index information to these files.
 * Free current memory, turn over control to the tmpfiles.
 *           
 * Return <eslOK>        on success; 
 *        <eslEDUP>      if one of the external tmpfiles already exists;
 *        <eslENOTFOUND> if we can't open a tmpfile for writing.
 */
static int
activate_external_sort(ESL_NEWSSI *ns)
{
  int status;
  int i;

  if (ns->external)                   return eslOK; /* we already are external, fool */

  status = eslEDUP;
  if (esl_FileExists(ns->ptmpfile)) goto CLEANEXIT;
  if (esl_FileExists(ns->stmpfile)) goto CLEANEXIT;
  
  status = eslENOTFOUND;
  if ((ns->ptmp = fopen(ns->ptmpfile, "w")) == NULL) goto CLEANEXIT;
  if ((ns->stmp = fopen(ns->stmpfile, "w")) == NULL) goto CLEANEXIT;

  /* Flush the current indices.
   */
  ESL_DPRINTF1(("Switching to external sort - flushing new ssi to disk...\n"));
  for (i = 0; i < ns->nprimary; i++) {
    if (sizeof(off_t) == 4) {
      fprintf(ns->ptmp, "%s\t%u\t%lu\t%lu\t%lu\n", 
	      ns->pkeys[i].key, 
	      ns->pkeys[i].fnum,
	      (unsigned long) ns->pkeys[i].r_off, 
	      (unsigned long) ns->pkeys[i].d_off, 
	      (unsigned long) ns->pkeys[i].len);
    } else {
      fprintf(ns->ptmp, "%s\t%u\t%llu\t%llu\t%lu\n", 
	      ns->pkeys[i].key, 
	      ns->pkeys[i].fnum,
	      (unsigned long long) ns->pkeys[i].r_off, 
	      (unsigned long long) ns->pkeys[i].d_off, 
	      (unsigned long) ns->pkeys[i].len);
    }
  }
  for (i = 0; i < ns->nsecondary; i++)
    fprintf(ns->stmp, "%s\t%s\n", ns->skeys[i].key, ns->skeys[i].pkey);
  
  /* Free the memory now that we've flushed our lists to disk
   */
  for (i = 0; i < ns->nprimary;   i++) free(ns->pkeys[i].key);
  for (i = 0; i < ns->nsecondary; i++) free(ns->skeys[i].key);
  for (i = 0; i < ns->nsecondary; i++) free(ns->skeys[i].pkey);
  if (ns->pkeys != NULL) free(ns->pkeys);       	
  if (ns->skeys != NULL) free(ns->skeys);       
  ns->pkeys = NULL;
  ns->skeys = NULL;

  status  = eslOK;

  /* Turn control over to external accumulation mode.
   */
 CLEANEXIT:
  if (status == eslOK) {
    ns->external = TRUE;
    return status;
  } else {
    if (ns->ptmp != NULL) { fclose(ns->ptmp); ns->ptmp = NULL; }
    if (ns->stmp != NULL) { fclose(ns->stmp); ns->stmp = NULL; }
    return status;
  }
}

/* parse_pkey(), parse_skey()
 * 
 * Given a <buf> containing a line read from the external
 * primary-key or secondary-key tmpfile; parse it, and fill in the fields of
 * <pkey> or <skey>
 * 
 * <?key> is a ptr to a structure on the stack. It is assumed
 * to be in use only transiently.
 * <?key>'s strings become ptrs into <buf>'s space, so we don't have to
 * allocate new space for them. This means that the transient <?key> structure
 * is only usable until <buf> is modified or free'd.
 * 
 * Returns <eslOK> on success.
 * 
 * Throws  <eslEFORMAT>        on parse error (shouldn't happen; we created it!)
 *         <eslEINCONCEIVABLE> if we can't deal with off_t's size.     
 */
static int
parse_pkey(char *buf, ESL_PKEY *pkey)
{
  int   status;
  char *s, *tok;
  int   n;
  
  s = buf;
  if (esl_strtok(&s, "\t\n", &(pkey->key), &n) != eslOK)
    ESL_DIE(eslEFORMAT, "parse failed");
  if (esl_strtok(&s, "\t\n", &tok,         &n) != eslOK) 
    ESL_DIE(eslEFORMAT, "parse failed");

  pkey->fnum = (uint16_t) atoi(tok);

  if (esl_strtok(&s, "\t\n", &tok, &n) != eslOK) ESL_DIE(eslEFORMAT, "parse failed");
  if      (sizeof(off_t) == 4) pkey->r_off  = (off_t) strtoul (tok, NULL, 10);
  else if (sizeof(off_t) == 8) pkey->r_off  = (off_t) strtoull(tok, NULL, 10);
  else                         ESL_DIE(eslEINCONCEIVABLE, "whoa - weird off_t");

  if (esl_strtok(&s, "\t\n", &tok, &n) != eslOK) ESL_DIE(eslEFORMAT, "parse failed");
  pkey->len = (uint32_t) strtoul(tok, NULL, 10);
  status = eslOK;

 CLEANEXIT:
  return status;
}
static int
parse_skey(char *buf, ESL_SKEY *skey)
{
  int   status;
  char *s;
  int   n;
  
  s = buf;
  if (esl_strtok(&s, "\t\n", &(skey->key),  &n) != eslOK) ESL_DIE(eslEFORMAT, "parse failed");
  if (esl_strtok(&s, "\t\n", &(skey->pkey), &n) != eslOK) ESL_DIE(eslEFORMAT, "parse failed");
  status = eslOK;

CLEANEXIT:
  return status;
}

/* ordering functions needed for qsort() */
static int 
pkeysort(const void *k1, const void *k2)
{
  ESL_PKEY *key1;
  ESL_PKEY *key2;
  key1 = (ESL_PKEY *) k1;
  key2 = (ESL_PKEY *) k2;
  return strcmp(key1->key, key2->key);
}
static int 
skeysort(const void *k1, const void *k2)
{
  ESL_SKEY *key1;
  ESL_SKEY *key2;
  key1 = (ESL_SKEY *) k1;
  key2 = (ESL_SKEY *) k2;
  return strcmp(key1->key, key2->key);
}


/*****************************************************************
 * 3. Functions for platform-independent binary i/o
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
  uint32_t  x32;
  uint64_t  x64;

  if      (mode == 64 && sizeof(off_t) == 8) 
    {
      if (esl_fread_i64(fp, &x64) != eslOK) { *ret_offset = 0; return eslFAIL; }
      *ret_offset = (off_t) x64;
      return eslOK;
    }
  else if (mode == 32)
    {
      if (esl_fread_i32(fp, &x32) != eslOK) { *ret_offset = 0; return eslFAIL; }
      *ret_offset = (off_t) x32;
      return eslOK;
    }
  else if (mode != 32 && mode != 64) 
    ESL_ERROR(eslEINVAL, "mode must be 32 or 64");
  else 
    ESL_ERROR(eslEINCOMPAT, "can't read 64-bit off_t on this 32-bit host");
  /*NOTREACHED*/
  return eslEINCONCEIVABLE;
}

/* Function:  esl_fwrite_offset()
 * Incept:    SRE, Fri Mar  3 13:35:04 2006 [St. Louis]
 *
 * Purpose:   Portably write (save) <offset> to the stream <fp>, in network
 *            byte order. 
 *
 * Returns:   <eslOK> on success; <eslFAIL> on write failure.
 *
 * Throws:    <eslEINVAL> if <off_t> is something other than a 32-bit or
 *            64-bit integer on this machine, in which case we don't know
 *            how to deal with it portably.
 */
int
esl_fwrite_offset(FILE *fp, off_t offset)
{
  if      (sizeof(off_t) == 4) return esl_fwrite_i32(fp, offset);
  else if (sizeof(off_t) == 8) return esl_fwrite_i64(fp, offset);
  else ESL_ERROR(eslEINVAL, "off_t is neither 32-bit nor 64-bit");
  /*UNREACHED*/
  return eslEINCONCEIVABLE;
}



/*****************************************************************
 * 4. Example code.
 ****************************************************************/
#ifdef eslSSI_EXAMPLE
/*::cexcerpt::ssi_example::begin::*/
#include <stdio.h>
#include <easel.h>
#include <esl_ssi.h>

int main(int argc, char **argv)
{
  ESL_NEWSSI *ns;
  char    *fafile;              /* name of FASTA file                   */
  FILE    *fp;                  /* opened FASTA file for reading        */
  char    *ssifile;             /* name of SSI file                     */
  FILE    *sfp;                 /* opened SSI file for writing          */
  uint16_t fh;                  /* file handle SSI associates w/ fafile */
  char    *buf = NULL;          /* growable buffer for esl_fgets()      */
  int      n   = 0;             /* length of buf                        */
  char    *s, *seqname;		
  off_t    seq_offset;

  /* Collect the sequence names from a FASTA file into an index */
  fafile = argv[1];
  ns = esl_newssi_Create();
  if ((fp = fopen(fafile, "r"))                  == NULL)  esl_fatal("failed to open %s", fafile);
  if (esl_newssi_AddFile(ns, fafile, 1, &fh)     != eslOK) esl_fatal("failed to add %s to index", fafile);

  seq_offset = ftello(fp);
  while (esl_fgets(&buf, &n, fp) == eslOK)
    {
      if (*buf == '>') {
	s = buf+1;                              /* skip past >                */
	esl_strtok(s, " \t\n", &seqname, NULL); /* name = 1st token on > line */
	esl_newssi_AddKey(ns, seqname, fh, seq_offset, 0, 0);
      }
      seq_offset = ftello(fp);				 
    }
  free(buf);
  fclose(fp);

  /* Save the index to disk */
  if (esl_FileNewSuffix(fafile, "ssi", &ssifile) != eslOK) esl_fatal("failed to name ssi file for %s", fafile);
  if ((sfp = fopen(ssifile, "wb"))               == NULL)  esl_fatal("failed to open SSI file %s", ssifile);
  if (esl_newssi_Write(sfp, ns)                  != eslOK) esl_fatal("failed to write ssi file");
  fclose(sfp);
  free(ssifile);
  esl_newssi_Destroy(ns);  

  return 0;
}
/*::cexcerpt::ssi_example::end::*/
#endif /*eslSSI_EXAMPLE*/


#ifdef eslSSI_EXAMPLE2
/*::cexcerpt::ssi_example2::begin::*/
#include <stdio.h>
#include <easel.h>
#include <esl_ssi.h>

int main(int argc, char **argv)
{
  ESL_SSI *ssi;
  char    *seqname;             /* name of sequence to retrieve         */
  char    *ssifile;             /* name of SSI file                     */
  uint16_t fh;                  /* file handle SSI associates w/ fafile */
  char    *fafile;              /* name of FASTA file                   */
  int      fmt;                 /* format code (1, in this example)     */
  off_t    offset;              /* disk offset of seqname in fafile     */
  FILE    *fp;                  /* opened FASTA file for reading        */
  char    *buf = NULL;          /* growable buffer for esl_fgets()      */
  int      n = 0;               /* size of buffer                       */

  seqname = argv[1];
  ssifile = argv[2];

  if (esl_ssi_Open(ssifile, &ssi)                  != eslOK)  esl_fatal("failed to open %s", ssifile);
  if (esl_ssi_FindName(ssi, seqname, &fh, &offset) != eslOK)  esl_fatal("failed to find key %s", seqname);
  if (esl_ssi_FileInfo(ssi, fh, &fafile, &fmt)     != eslOK)  esl_fatal("failed to get filename %d\n", fh);
  esl_ssi_Close(ssi);  

  if ((fp = fopen(fafile, "r"))                    == NULL)   esl_fatal("failed to open file %s", fafile);
  if (fseeko(fp, offset, SEEK_SET)                 != 0)      esl_fatal("failed to position file %s", fafile);
  if (esl_fgets(&buf, &n, fp)                      != eslOK)  esl_fatal("failed to get name/desc line");
  do {
    puts(buf); 
  } while (esl_fgets(&buf, &n, fp) != eslOK && *buf != '>');
  
  fclose(fp);  
  free(buf);
  return 0;
}
/*::cexcerpt::ssi_example2::end::*/
#endif /*eslSSI_EXAMPLE2*/

/*****************************************************************
 * 5. Test driver
 *****************************************************************/ 

/* gcc -g -Wall -o testdrive -L. -I. -DeslSSI_TESTDRIVE esl_ssi.c -leasel -lm */
#ifdef eslSSI_TESTDRIVE
#include <stdio.h>
#include <string.h>
#include <easel.h>
#include <esl_sqio.h>
#include <esl_ssi.h>
#include <esl_random.h>

int main(int argc, char **argv)
{
  ESL_RANDOMNESS *r    = NULL;
  ESL_NEWSSI     *ns   = NULL;
  ESL_SSI        *ssi  = NULL;
  ESL_SQ         *sq   = NULL;
  ESL_SQFILE     *sqfp = NULL;
  FILE  *fp;
  FILE  *ssifp;
  uint16_t fh;
  int    nseq, nfiles, i, j, maxL;
  char **sqfile;
  char **seqname;
  char **seq;
  int   *seqlen;
  int    nq;
  char   query[32];
  char  *qfile;
  int    qfmt;
  off_t  roff;
  double p[4] = { 0.25, 0.25, 0.25, 0.25 };
  char  *ssifile = NULL;
  int    status;
  int    be_verbose;

  nfiles = 3;
  nseq   = 10;
  maxL   = 1000;
  nq     = 10;
  be_verbose = FALSE;
  
  /* Create <nfiles> sequence file names.
   */
  ESL_MALLOC(sqfile, sizeof(char *) * nfiles);
  for (j = 0; j < nfiles; j++)
    {
      ESL_MALLOC(sqfile[j], sizeof(char) * 32);
      sprintf(sqfile[j], "xxxssi-test%d.fa", j);
    } 

  /* Create <nfiles*nseq> sequences with random 
   * lengths up to 1000.
   */
  /* r = esl_randomness_CreateTimeseeded(); */
  r = esl_randomness_Create(39); 
  ESL_MALLOC(seq,    sizeof(char *) * nseq * nfiles);
  ESL_MALLOC(seqname,sizeof(char *) * nseq * nfiles);
  ESL_MALLOC(seqlen, sizeof(int)    * nseq * nfiles);
  for (i = 0; i < nseq*nfiles; i++)
    {
      seqlen[i] = 1 + esl_rnd_Choose(r, maxL); /* 1..maxL */
      seq[i]    = esl_rnd_IID(r, "ACGT", p, 4, seqlen[i]);
      ESL_MALLOC(seqname[i],sizeof(char) * 64);
      sprintf(seqname[i], "seq%d-file%d", i, i/nseq);
    }

  /* Save them to FASTA files.
   */
  for (j = 0; j < nfiles; j++)
    {
      if ((fp = fopen(sqfile[j], "w")) == NULL) 
	esl_fatal("failed to open %s", sqfile[j]);
      for (i = j*nseq; i < (j+1)*nseq; i++)
	{
	  sq = esl_sq_CreateFrom(seqname[i], seq[i], NULL, NULL, NULL);
	  esl_sqio_Write(fp, sq, eslSQFILE_FASTA);
	  esl_sq_Destroy(sq);
	}
      fclose(fp);
    }

  /* Create an ssi index of all the FASTA files.
   */
  ns = esl_newssi_Create();
  sq = esl_sq_Create();
  for (j = 0; j < nfiles; j++)
    {
      if (esl_sqfile_Open(sqfile[j], eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK)
	esl_fatal("failed to open fasta file %s", sqfile[j]);
      esl_newssi_AddFile(ns, sqfile[j], sqfp->format, &fh);
      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  if (be_verbose) printf("%16s  %ld  %ld  %d\n", sq->name, (long) sq->roff, (long) sq->doff, sq->n);
	  esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->n);
	  esl_sq_Reuse(sq);
	}
      if (status != eslEOF) esl_fatal("sequence read failure");
      esl_sqfile_Close(sqfp);
    }
  esl_sq_Destroy(sq);

  /* Save the SSI index to a file.
   */
  esl_strdup("xxxssi-test.ssi", -1, &ssifile);
  if ((ssifp = fopen(ssifile, "w")) == NULL) esl_fatal("failed to open ssi file");
  esl_newssi_Write(ssifp, ns);
  fclose(ssifp);
  esl_newssi_Destroy(ns);
  
  /* Open the SSI index - now we'll use it to retrieve
   * <nq> random sequences.
   */
  if (esl_ssi_Open(ssifile, &ssi) != eslOK) esl_fatal("failed to open ssi index");
  sq = esl_sq_Create();
  while (nq--)
    {
      /* Choose a seq and file */
      i = esl_rnd_Choose(r, nseq*nfiles);
      j = i/nseq;
      sprintf(query, "seq%d-file%d", i, j);

      /* Retrieve it */
      status = esl_ssi_FindName(ssi, query, &fh, &roff);
      if (status != eslOK) esl_fatal("didn't find %s in index", query);
      esl_ssi_FileInfo(ssi, fh, &qfile, &qfmt);      
      if (esl_sqfile_Open(qfile, qfmt, NULL, &sqfp) != eslOK)
	esl_fatal("failed to open fasta file %s", qfile);
      esl_sqio_Position(sqfp, roff);
      if (esl_sqio_Read(sqfp, sq) != eslOK) esl_fatal("failed to read seq %s", query);

      /* Check that it's the right one */
      if (strcmp(sq->name, query) != 0)  esl_fatal("sought %s, retrieved %s", query, sq->name);
      if (sq->n != seqlen[i])            esl_fatal("wrong sequence length retrieved");
      if (strcmp(sq->seq,  seq[i]) != 0) esl_fatal("unexpected sequence retrieved");
      if (strcmp(qfile, sqfile[j]) != 0) esl_fatal("file names %s and %s differ", qfile, sqfile[j]);

      esl_sq_Reuse(sq);
      esl_sqfile_Close(sqfp);
    }
  esl_sq_Destroy(sq);
  esl_ssi_Close(ssi);

  free(ssifile);
  free(seqlen);
  esl_Free2D((void **) seqname, nseq*nfiles);
  esl_Free2D((void **) seq,     nseq*nfiles);
  esl_randomness_Destroy(r);
  esl_Free2D((void **) sqfile, nfiles);
  return 0;
}
#endif /*eslSSI_TESTDRIVE*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
