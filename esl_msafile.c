/* Multiple sequence alignment file i/o 
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_mem.h"
#include "esl_msafile.h"

/* Function:  eslx_msafile_Open()
 * Synopsis:  Open a multiple sequence alignment file for input.
 *
 * Purpose:   Open a multiple sequence alignment file <msafile> for input.
 *
 *            Caller asserts that <msafile> is in format code
 *            <format>, such as <eslMSAFILE_STOCKHOLM>,
 *            <eslMSAFILE_AFA>, <eslMSAFILE_CLUSTAL>.  If <format> is
 *            <eslMSAFILE_UNKNOWN>, format autodetection is performed.
 *            
 *            Optionally, caller can provide in <env> the name of an
 *            environment variable ("PFAMDB", perhaps), in which the
 *            routine can find a colon-delimited list of directories.
 *            Then, if <msafile> is not found in the current working
 *            directory, we look for it in these directories, in the
 *            order they're listed. Otherwise caller sets <env> to
 *            <NULL>.
 *
 *            <msafile> is usually the name of a file. Alignments may
 *            also be read from standard input, or from
 *            gzip-compressed files.  If <msafile> is ``-'', alignment
 *            input is taken from the standard input stream. If
 *            <msafile> ends in ``.gz'', alignment input is read
 *            through a pipe from <gzip -dc>.
 *
 * Args:      msafile   - name of alignment file, or ``-''.
 *            format    - format code, such as <eslMSAFILE_STOCKHOLM>;
 *                        or <eslMSAFILE_UNKNOWN> to autodetect format.
 *            env       - <NULL>, or the name of an environment variable
 *                        containing colon-delimited list of directories
 *                        in which to search for <msafile>.
 *            *ret_afp  - RETURN: open MSA input stream.
 *
 * Returns:   <eslOK> on success, and <*ret_afp> is the newly opened msa file.
 *
 *            <eslENOTFOUND> if <msafile> doesn't exist or can't be opened for
 *            reading; or (in the case of a <.gz> file) if a <gzip> executable
 *            doesn't exist in user's <PATH> or can't be executed. 
 *
 *            <eslFAIL> in the case of a <.gz> file and the <gzip -dc> command
 *            fails on it.
 *            
 *            <eslEFORMAT> if we tried to autodetect the file format (caller
 *            provided <format=eslMSAFILE_UNKNOWN>), and failed.
 *
 *            On any of these normal errors, <*ret_afp> is returned in an error
 *            state, containing a user-directed error message in <afp->errmsg>
 *            and (if relevant) the full path to <msafile> that we attempted to
 *            open in <afp->bf->filename>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> on a system call failure, such as <fread()>.
 *            <eslEINVAL> if we tried to use <stdin> but the <stdin> stream was
 *            invalid (in an error state, <NULL>, or at <EOF>).
 *            On thrown exceptions, <*ret_afp> is <NULL>.
 */
int
eslx_msafile_Open(const char *msafile, int format, const char *env, ESLX_MSAFILE **ret_afp)
{
  ESLX_MSAFILE *afp = NULL;
  int           status;
  int           x;

  ESL_ALLOC(afp, sizeof(ESLX_MSAFILE));
  afp->bf         = NULL;
  afp->linenumber = 0;
  afp->abc        = NULL;
  afp->ssi        = NULL;

  for (x = 0; x < 128; x++)
    afp->inmap[x] = eslDSQ_ILLEGAL;

  if ((status = esl_buffer_Open(msafile, env, &(afp->bf))) != eslOK)
    ESL_XFAIL(status, afp->errmsg, "%s", afp->bf->errmsg); /* ENOTFOUND; FAIL are normal here */

  if (format == eslMSAFILE_UNKNOWN && (status = eslx_msafile_GuessFileFormat(afp->bf, &format)) != eslOK)
    ESL_XFAIL(status, afp->errmsg, "couldn't determine alignment input format"); /* EFORMAT is normal failure */
  afp->format = format;

  *ret_afp = afp;
  return eslOK;

 ERROR:
  if (status == eslENOTFOUND || status == eslFAIL || status == eslEFORMAT) {
    *ret_afp = afp;
    return status;
  } else {
    eslx_msafile_Close(afp);
    *ret_afp = NULL;
    return status;
  }
}

/* Function:  eslx_msafile_OpenFailure()
 * Synopsis:  Report diagnostics of normal error in opening MSA file, and exit.
 *
 * Purpose:   Report user-directed diagnostics of a normal error in opening
 *            an MSA input. Print a message to <stderr>, then exit. 
 */
void
eslx_msafile_OpenFailure(ESLX_MSAFILE *afp, int status)
{
  fprintf(stderr, "alignment input open failed:\n   %s\n", afp->errmsg); 
  
  eslx_msafile_Close(afp);
  exit(status);
}


/* Function:  eslx_msafile_ReadFailure()
 * Synopsis:  Report diagnostics of a normal error in parsing MSA file, and exit.
 *
 * Purpose:   Report user-directed diagnostics of a normal error from
 *            parsing an MSA file.  Output the error message to
 *            <stderr>, along with information about what we were
 *            parsing (filename, if it was a file) and where we were
 *            in the input (linenumber, if we know it). This
 *            information is all available in <afp>. Then close <afp>
 *            and exit with the <status> provided by the caller.
 *
 * Args:      afp    - open ESLX_MSAFILE, containing information about
 *                     the error and the input source.
 *            status - exit status; generally eslEFORMAT. 
 *
 * Returns:   no return. Exits here with <status>.
 */
void
eslx_msafile_ReadFailure(ESLX_MSAFILE *afp, int status)
{
  switch (status) {
  case eslEFORMAT:  fprintf(stderr, "alignment input parse error: %s\n", afp->errmsg);           break;
  case eslEOF:      fprintf(stderr, "alignment input appears empty?\n");                         break;
  default:          fprintf(stderr, "alignment input read error; unexpected code %d\n", status); break;
  }
  
  switch (afp->bf->mode_is) {
  case eslBUFFER_STREAM:   fprintf(stderr, "   while reading from an input stream (not a file)\n");   break;
  case eslBUFFER_CMDPIPE:  fprintf(stderr, "   while reading through a pipe (not a file)\n");         break;
  case eslBUFFER_FILE:     
  case eslBUFFER_ALLFILE:
  case eslBUFFER_MMAP:     fprintf(stderr, "   while reading file %s\n", afp->bf->filename);          break;
  case eslBUFFER_STRING:   fprintf(stderr, "   while reading from a provided string (not a file)\n"); break;
  default:                 break; 
  }

  if (afp->linenumber > 0) fprintf(stderr, "   at or near line %" PRIu64 "\n", afp->linenumber);
  else                     fprintf(stderr, "   at or near byte %" PRIu64 "\n", esl_buffer_GetOffset(afp->bf));

  eslx_msafile_Close(afp);
  exit(status);
}


/* Function:  eslx_msafile_GuessFileFormat()
 * Synopsis:  Guess the MSA file format of an open buffer.
 *
 * Purpose:   Peek into an open buffer, and try to determine what
 *            alignment file format (if any) its input is in. If a
 *            format can be determined, return <eslOK> and set
 *            <*ret_fmtcode> to the format code.  If not, return
 *            <eslEFORMAT> and set <*ret_fmtcode> to
 *            <eslMSAFILE_UNKNOWN.  In either case, the buffer <bf> is
 *            restored to its original position upon return.
 *
 *            If the <bf> corresponds to an open file with a file
 *            name, we attempt to use the suffix as a clue. Suffix
 *            rules for alignment files are as follows:
 *                 | Stockholm     |  .sto .sth .stk |
 *                 | Aligned FASTA |  .afa .afasta   |
 *                 | CLUSTAL       |  .aln           |
 *                 | Pfam          |  .pfam          |
 *                 | MSF           |  .msf           |
 *                 | A2M           |  .a2m           | 
 *                 | SELEX         |  .slx .selex    |   
 *
 * Args:      bf          - the open buffer to read input from
 *            ret_fmtcode - RETURN: format code that's determined
 *
 * Returns:   <eslOK> on success, and <*ret_fmtcode> contains the format code.
 *            <eslEFORMAT> if format can't be guessed, and <*ret_fmtcode> contains
 *            <eslMSAFILE_UNKNOWN>.
 *
 * Throws:    (no abnormal error conditions)
 */
int
eslx_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode)
{
  esl_pos_t  initial_offset;
  char      *p;
  esl_pos_t  n;
  int        fmt_bysuffix    = eslMSAFILE_UNKNOWN;
  int        fmt_byfirstline = eslMSAFILE_UNKNOWN;
  int        status;

  /* As we start, save parser status:
   *   remember the offset where we started (usually 0, but not necessarily)
   *   set an anchor to be sure that this offset stays in the buffer's memory
   */
  initial_offset = esl_buffer_GetOffset(bf);
  esl_buffer_SetAnchor(bf, initial_offset);

  /* First we try to guess based on the filename suffix.
   * (if there's a filename, and if it has a suffix, anyway.)
   */
  if (bf->filename)
    {
      esl_file_Extension(bf->filename, 0, &p, &n);
      if (esl_memstrcmp(p, n, ".gz"))
	esl_file_Extension(bf->filename, 3, &p, &n);

      if (p)
	{
	  if      (esl_memstrcmp(p, n, ".sto"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".stk"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".sth"))    fmt_bysuffix = eslMSAFILE_STOCKHOLM;
	  else if (esl_memstrcmp(p, n, ".pfam"))   fmt_bysuffix = eslMSAFILE_PFAM;
	  else if (esl_memstrcmp(p, n, ".a2m"))    fmt_bysuffix = eslMSAFILE_A2M;
	  else if (esl_memstrcmp(p, n, ".psi"))    fmt_bysuffix = eslMSAFILE_PSIBLAST;
	  else if (esl_memstrcmp(p, n, ".slx"))    fmt_bysuffix = eslMSAFILE_SELEX;
	  else if (esl_memstrcmp(p, n, ".selex"))  fmt_bysuffix = eslMSAFILE_SELEX;
	  else if (esl_memstrcmp(p, n, ".afa"))    fmt_bysuffix = eslMSAFILE_AFA;
	  else if (esl_memstrcmp(p, n, ".afasta")) fmt_bysuffix = eslMSAFILE_AFA;
	  else if (esl_memstrcmp(p, n, ".aln"))    fmt_bysuffix = eslMSAFILE_CLUSTAL;
	}
    }

  /* Second, we peek at the first non-blank line of the file.
   * Multiple sequence alignment files are generally identifiable by 
   * a token on this line.
   */
  /* Skip blank lines, get first non-blank one */
  do   { 
    status = esl_buffer_GetLine(bf, &p, &n);
  } while (status == eslOK && esl_memspn(p, n, " \t\r\n") == n);

  if      (esl_memstrpfx(p, n, "# STOCKHOLM")) fmt_byfirstline = eslMSAFILE_STOCKHOLM;
  else if (esl_memstrpfx(p, n, ">"))           fmt_byfirstline = eslMSAFILE_AFA;
  else if (esl_memstrpfx(p, n, "CLUSTAL"))     fmt_byfirstline = eslMSAFILE_CLUSTAL;
  else if (esl_memstrpfx(p, n, "MUSCLE"))      fmt_byfirstline = eslMSAFILE_MUSCLE;  

  if      (fmt_byfirstline != eslMSAFILE_UNKNOWN && fmt_bysuffix == eslMSAFILE_UNKNOWN)  *ret_fmtcode = fmt_byfirstline;
  else if (fmt_byfirstline == eslMSAFILE_UNKNOWN && fmt_bysuffix != eslMSAFILE_UNKNOWN)  *ret_fmtcode = fmt_bysuffix;     
  else if (fmt_byfirstline == fmt_bysuffix)                                              *ret_fmtcode = fmt_byfirstline;   /* including MSAFILE_UNKNOWN */
  else if (fmt_byfirstline == eslMSAFILE_MUSCLE)                                         *ret_fmtcode = fmt_byfirstline;   /* MUSCLE files can have any suffix */
  else    *ret_fmtcode = fmt_bysuffix;

  /* As we return, restore parser status:
   *  put it back where it was when we started;
   *  clear the anchor that made sure that position stayed in memory.
   */
  esl_buffer_SetOffset  (bf, initial_offset);
  esl_buffer_RaiseAnchor(bf, initial_offset);
  return ((*ret_fmtcode == eslMSAFILE_UNKNOWN) ? eslEFORMAT: eslOK);
}

void
eslx_msafile_Close(ESLX_MSAFILE *afp)
{
  if (! afp) return;
  if (afp->bf)  esl_buffer_Close(afp->bf);
  if (afp->ssi) esl_ssi_Close(afp->ssi);
  free(afp);
}

      

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

