/* Support for MPI parallelization.
 * 
 * SRE, Sat Jun  2 09:07:25 2007 [Janelia]
 * SVN $Id$
 */

#ifdef HAVE_MPI
#ifndef eslMPI_INCLUDED
#define eslMPI_INCLUDED

extern int esl_mpi_PackOpt(void *inbuf, int incount, MPI_Datatype type, void *pack_buf, 
			   int pack_buf_size, int *position, MPI_Comm comm);
extern int esl_mpi_PackOptSize(void *inbuf, int incount, MPI_Datatype type, MPI_Comm comm, int *ret_n);
extern int esl_mpi_UnpackOpt(void *pack_buf, int pack_buf_size, int *pos, void **outbuf, 
			     int *opt_n, MPI_Datatype type, MPI_Comm comm);


#endif /*eslMPI_INCLUDED*/
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
