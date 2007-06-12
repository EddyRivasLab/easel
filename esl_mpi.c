/* Support for MPI parallelization.
 * 
 * SRE, Sat Jun  2 08:16:14 2007 [Janelia] [Hitchhiker's Guide, Tertiary Phase]
 * SVN $Id$
 */

#include "esl_config.h"		
#ifdef HAVE_MPI
#include <string.h>
#include "mpi.h"

#include "easel.h"
#include "esl_mpi.h"


/* Function:  esl_mpi_PackOpt()
 * Synopsis:  Pack an optional, variable-sized array (or string).
 * Incept:    SRE, Sat Jun  2 08:40:39 2007 [Janelia]
 *
 * Purpose:   Pack data array <inbuf> of <incount> elements of type <type> into
 *            an MPI packed buffer <pack_buf> of total size <pack_buf_size> destined
 *            for MPI communicator <comm> that is currently filled to position <*position>.
 *            
 *            <inbuf> may be <NULL>, in which case <incount> is
 *            assumed to be 0, and a `null array' is packed that
 *            <esl_mpi_UnpackOpt()> knows how to decode as a <NULL>
 *            pointer.
 *            
 *            As a special case for strings, if <type> is <MPI_CHAR>,
 *            <incount> may be passed as <-1> to indicate `unknown';
 *            the routine will use <strlen(inbuf)+1> to determine the
 *            size of the string including its <NUL> terminator.
 *
 * Returns:   <eslOK> on success, the array is packed into <pack_buf>, 
 *            and the <*position> counter is updated to point to the next byte
 *            in <pack_buf> after the packed array.
 * 
 * Throws:    <eslESYS> if an MPI call fails.
 */
int
esl_mpi_PackOpt(void *inbuf, int incount, MPI_Datatype type, void *pack_buf, int pack_buf_size, int *position, MPI_Comm comm)
{
  if (inbuf == NULL) {
    incount = 0;
    if (MPI_Pack(&incount,    1, MPI_INT, pack_buf, pack_buf_size, position, comm) != 0)  ESL_EXCEPTION(eslESYS, "MPI pack failed");
  } else {
    if (incount == -1 && type == MPI_CHAR) incount = strlen(inbuf) + 1;
    if (MPI_Pack(&incount,    1, MPI_INT, pack_buf, pack_buf_size, position, comm) != 0)  ESL_EXCEPTION(eslESYS, "MPI pack failed");
    if (MPI_Pack(inbuf, incount,    type, pack_buf, pack_buf_size, position, comm) != 0)  ESL_EXCEPTION(eslESYS, "MPI pack failed");
  }
  return eslOK;
}

/* Function:  esl_mpi_PackOptSize()
 * Synopsis:  Determine the size of a packed optional, variable-sized array.
 * Incept:    SRE, Sat Jun  2 10:09:16 2007 [Janelia]
 *
 * Purpose:   Determine an upper bound on the size (in bytes) required
 *            to pack an array <inbuf> of <incount> elements of type
 *            <type> destined for MPI communicator <comm> using
 *            <esl_mpi_PackOpt()>, and return it in <*ret_n>.
 *            
 *            If <inbuf> is non-<NULL>, the packed message consists
 *            of 1 integer (the length, <incount>) followed by the array.
 *            If <inbuf> is <NULL>, the packed message consists of one
 *            integer (0). 
 *            
 *            As a special case for strings, if <type> is <MPI_CHAR>,
 *            <incount> may be passed as <-1> to indicate `unknown';
 *            in this case, the routine uses <strlen(inbuf)+1> to determine the
 *            size of the string including its <NUL> terminator.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the upper limit size in
 *            bytes.
 *            
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.           
 */
int
esl_mpi_PackOptSize(void *inbuf, int incount, MPI_Datatype type, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;

  *ret_n = 0;
  if (inbuf == NULL) {
    status = MPI_Pack_size(1,    MPI_INT, MPI_COMM_WORLD, &sz);  *ret_n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
  } else {
    if (incount == -1 && type == MPI_CHAR) incount = strlen(inbuf) + 1;
    status = MPI_Pack_size(1,    MPI_INT, MPI_COMM_WORLD, &sz);  *ret_n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
    status = MPI_Pack_size(incount, type, MPI_COMM_WORLD, &sz);  *ret_n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
  }
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}



/* Function:  esl_mpi_UnpackOpt()
 * Synopsis:  Unpack an optional, variable-sized array (or string).
 * Incept:    SRE, Sat Jun  2 08:39:39 2007 [Janelia]
 *
 * Purpose:   Unpack a packed MPI message in buffer <pack_buf>, of total size
 *            <pack_buf_size>, at current position <*pos> in <pack_buf>,
 *            for MPI communicator <comm>, where the next packed element is an optional
 *            array of type <type>, consisting of a <(n,data)> pair, with <n=0>
 *            indicating no data. 
 *            
 *            If array data is present (<n>0>), allocate <*outbuf>,
 *            put the array in it, and optionally return <n> in
 *            <*opt_n>. The caller is responsible for free'ing this
 *            <*outbuf>.
 *
 *            If data are not present (<n=0>), no allocation is done,
 *            <*outbuf> is set to <NULL>, and the optional <*opt_n> is
 *            0.
 *            
 *            <*pos> is updated to point at the next element in <pack_buf>
 *            that needs to be unpacked.
 *
 *            This routine is designed for an optional-array idiom in
 *            which <array==NULL> means the array isn't available, and
 *            otherwise the array contains valid data. For instance,
 *            this is used for optional annotation on multiple
 *            alignments. 
 *            
 * Returns:   <eslOK> on success; <*pos> is updated; <*outbuf> is either a newly allocated 
 *            array (that caller is responsible for freeing) and optional <*opt_n>
 *            is its length, or <*outbuf> is <NULL> and optional <*opt_n> is 0.
 *
 * Throws:    <eslESYS> on an MPI call failure; <eslEINVAL> if something's wrong
 *            with the arguments; <eslEMEM> on allocation failure. 
 *            In either case, <*outbuf> is <NULL> and optional <*opt_n> is 0.
 */
int
esl_mpi_UnpackOpt(void *pack_buf, int pack_buf_size, int *pos, void **outbuf, int *opt_n, MPI_Datatype type, MPI_Comm comm)
{
  int sz;
  int status;
  
  if (MPI_Unpack(pack_buf, pack_buf_size, pos, &sz, 1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if (sz == 0) {
    *outbuf = NULL;
  } else {
    if      (type == MPI_CHAR)           ESL_ALLOC(*outbuf, sizeof(char)           * sz); 
    else if (type == MPI_SHORT)          ESL_ALLOC(*outbuf, sizeof(short)          * sz); 
    else if (type == MPI_INT)            ESL_ALLOC(*outbuf, sizeof(int)            * sz);  
    else if (type == MPI_LONG)           ESL_ALLOC(*outbuf, sizeof(long)           * sz);  
    else if (type == MPI_UNSIGNED_CHAR)  ESL_ALLOC(*outbuf, sizeof(unsigned char)  * sz);
    else if (type == MPI_UNSIGNED_SHORT) ESL_ALLOC(*outbuf, sizeof(unsigned short) * sz);  
    else if (type == MPI_UNSIGNED)       ESL_ALLOC(*outbuf, sizeof(unsigned int)   * sz);
    else if (type == MPI_UNSIGNED_LONG)  ESL_ALLOC(*outbuf, sizeof(unsigned long)  * sz); 
    else if (type == MPI_FLOAT)          ESL_ALLOC(*outbuf, sizeof(float)          * sz); 
    else if (type == MPI_DOUBLE)         ESL_ALLOC(*outbuf, sizeof(double)         * sz);  
    else if (type == MPI_LONG_DOUBLE)    ESL_ALLOC(*outbuf, sizeof(long double)    * sz); 
    else if (type == MPI_BYTE)           ESL_ALLOC(*outbuf, sz);                       
    else if (type == MPI_PACKED)         ESL_ALLOC(*outbuf, sz);
    else ESL_XEXCEPTION(eslEINVAL, "no such MPI datatype");

    if (MPI_Unpack(pack_buf, pack_buf_size, pos, *outbuf, sz, type, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  }
  if (opt_n != NULL) *opt_n = sz;
  return eslOK;

 ERROR:
  if (*outbuf != NULL) free(*outbuf);
  *outbuf = NULL;
  if (opt_n != NULL) *opt_n = 0;
  return status;
}

#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
