/* esl_ssi.c
 * simple sequence indices: fast lookup in large sequence files by keyword.
 * 
 * SVN $Id$
 * adapted from squid's ssi.c
 * SRE, Thu Mar  2 18:46:02 2006 [St. Louis]
 */

#include <esl_config.h>

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
