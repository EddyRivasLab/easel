/* i/o of multiple sequence alignment files in PHYLIP format
 * 
 * 
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#include "esl_mem.h"
#include "esl_msa.h"
#include "esl_msafile.h"
//#include "esl_msafile_phylip.h"


/* Function:  esl_msafile_phylip_SetInmap()
 * Synopsis:  Configure input map for PHYLIP formats.
 *
 * Purpose:   Set the <afp->inmap> for PHYLIP formats.
 * 
 *            Phylip documentation states that DNA programs accept
 *            'ABCDGHKMNORSTUVWXY?-', that 'a period was previously
 *            allowed' and that O means a deletion. Protein programs
 *            accept 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-', and while JOU
 *            are accepted, they are unused. 
 *
 *            So: in text mode, we accept any alphabetic character
 *            plus '-*?.', verbatim. '~_', which Easel would normally
 *            accept, are illegal. Whitespace and numbers are ignored.
 *            
 *            In digital mode, we modify the digital alphabet by
 *            demapping '~_' and making them illegal; '?' is mapped to
 *            missing data; whitespace and numbers are ignored;
 *            and ONLY in <eslDNA> or <eslRNA> alphabets, 'O' is
 *            mapped to a gap.
 *            
 *            The inconsistent mapping of 'O' poses potential
 *            problems. In text mode (where we don't know the
 *            alphabet, and thus don't know what to do with O), we
 *            input the O verbatim. In digital mode, in a DNA or RNA
 *            alphabet, we map O to a gap; in other digital alphabets,
 *            we use the default digital alphabet mapping of O.
 *            
 * Xref:      http://evolution.genetics.washington.edu/phylip/doc/sequence.html
 */
int 
esl_msafile_phylip_SetInmap(ESLX_MSAFILE *afp)
{
  int sym;

#ifdef eslAUGMENT_ALPHABET
  if (afp->abc)
    {
      for (sym = 1;   sym < 128; sym++) afp->inmap[sym] = afp->abc->inmap[sym];
      for (sym = '0'; sym < '9'; sym++) afp->inmap[sym] = eslDSQ_IGNORED;
      afp->inmap['?']  = esl_abc_XGetMissing(afp->abc);
      afp->inmap['~']  = eslDSQ_ILLEGAL;
      afp->inmap['_']  = eslDSQ_ILLEGAL;
      afp->inmap[' ']  = eslDSQ_IGNORED;
      afp->inmap['\t'] = eslDSQ_IGNORED;
      afp->inmap[0]    = esl_abc_XGetUnknown(afp->abc);

      if (afp->abc->type == eslDNA || afp->abc->type == eslRNA) 
	afp->inmap['O'] = esl_abc_XGetGap(afp->abc);
    }
#endif
  if (! afp->abc)
    {
      for (sym = 1; sym < 128; sym++)   afp->inmap[sym] = eslDSQ_ILLEGAL;
      for (sym = 'a'; sym < 'z'; sym++) afp->inmap[sym] = sym;
      for (sym = 'A'; sym < 'Z'; sym++) afp->inmap[sym] = sym;
      for (sym = '0'; sym < '9'; sym++) afp->inmap[sym] = eslDSQ_IGNORED;
      afp->inmap['-']  = '-';
      afp->inmap['*']  = '*';
      afp->inmap['?']  = '?';
      afp->inmap['.']  = '.';
      afp->inmap[' ']  = eslDSQ_IGNORED;
      afp->inmap['\t'] = eslDSQ_IGNORED;
      afp->inmap[0]    = '?';
    } 
  return eslOK;
}
  
  
/* Function:  esl_msafile_phylip_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open PHYLIP MSA input.
 *
 * Purpose:   Guess the alphabet of the sequences in open 
 *            PHYLIP format MSA file <afp>.
 *            
 *            On normal return, <*ret_type> is set to <eslDNA>,
 *            <eslRNA>, or <eslAMINO>, and <afp> is reset to its
 *            original point.
 *
 * Args:      afp      - open PHYLIP format MSA file
 *           *ret_type - RETURN: <eslDNA>, <eslRNA>, <eslAMINO>; or <eslUNKNOWN>.
 *
 * Returns:   <eslOK> on success.
 *            <eslENOALPHABET> if autodetection fails.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslESYS> on failures of fread() or other system calls
 */
int
esl_msafile_phylip_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type)
{
  int alphatype = eslUNKNOWN;
  int x;
  esl_pos_t anchor;
  int64_t ct[26];
  char *p;
  esl_pos_t n;
  int status;

  for (x = 0; x < 26; x++) ct[x] = 0;

  anchor = esl_buffer_GetOffset(afp->bf);
  if ((status = esl_buffer_SetAnchor(afp->bf, anchor)) != eslOK) { status = eslEINCONCEIVABLE; goto ERROR; } /* [eslINVAL] can't happen here */

  /* Find the first nonblank line, which says " <nseq> <alen>" and may also have options */
  while ( (status = eslx_msafile_GetLine(afp, &p, &n)) == eslOK  && esl_memspn(afp->line, afp->n, " \t") == afp->n) ;
  if      (status == eslEOF) ESL_XFAIL(eslENOALPHABET, afp->errmsg, "can't determine alphabet: no alignment data found");
  else if (status != eslOK)  goto ERROR;

 DONE:
  esl_buffer_SetOffset(afp->bf, anchor);   /* Rewind to where we were. */
  esl_buffer_RaiseAnchor(afp->bf, anchor);
  *ret_type = alphatype;
  return status;


 ERROR:
   if (anchor != -1) {
    esl_buffer_SetOffset(afp->bf, anchor);
    esl_buffer_RaiseAnchor(afp->bf, anchor);
  }
  *ret_type = eslUNKNOWN;
  return status;
}




static int
phylip_autodetect_details(ESL_BUFFER *bf, int *ret_isinterleaved, int *ret_namewidth)
{
  esl_pos_t anchor        = -1;
  char     *p, *tok;
  esl_pos_t n,  toklen, pos;
  int32_t   nseq, alen;		/* int32_t is because we're using esl_mem_strtoi32() to parse them out */
  int      *nci, *nci0;		/* number of chars per sequence if the format is interleaved: nci[0..nseq-1] */
  int      *ncs, *ncs0;		/* number of chars per sequence if the format is sequential:  ncs[0..nseq-1] */
  int       nc;
  int       nblock, nline;
  int       ncpb;               /* number of chars per interleaved block > 0  */
  int       idxi;		/* sequence index if the format is interleaved */
  int       idxs;		/* sequence index if the format is sequential */
  int       nillegal;
  int       is_interleaved, is_sequential;
  int       status;

  is_interleaved = is_sequential = TRUE;

  anchor = esl_buffer_GetOffset(bf);
  if ((status = esl_buffer_SetAnchor(bf, anchor)) != eslOK) { status = eslEINCONCEIVABLE; goto ERROR; } /* [eslINVAL] can't happen here */

  /* Find the first nonblank line, which says " <nseq> <alen>" and may also have options */
  while ( (status = esl_buffer_GetLine(bf, &p, &n)) == eslOK  && esl_memspn(p, n, " \t") == n) ;
  if      (status == eslEOF) ESL_XFAIL(eslENOALPHABET, bf->errmsg, "can't determine alphabet: no alignment data found");
  else if (status != eslOK)  goto ERROR;
  
  esl_memtok(&p, &n, " \t", &tok, &toklen);
  if (esl_mem_strtoi32(tok, toklen, 0, NULL, &nseq)  != eslOK) ESL_XFAIL(eslENOALPHABET, bf->errmsg, "can't determine alphabet: input is not in PHYLIP format");
  if (esl_memtok(&p, &n, " \t", &tok, &toklen)       != eslOK) ESL_XFAIL(eslENOALPHABET, bf->errmsg, "can't determine alphabet: input is not in PHYLIP format");
  if (esl_mem_strtoi32(tok, toklen, 0, NULL, &alen)  != eslOK) ESL_XFAIL(eslENOALPHABET, bf->errmsg, "can't determine alphabet: input is not in PHYLIP format");

  ESL_ALLOC(nci,  sizeof(int) * nseq);
  ESL_ALLOC(nci0, sizeof(int) * nseq);
  ESL_ALLOC(ncs,  sizeof(int) * nseq);
  ESL_ALLOC(ncs0, sizeof(int) * nseq);
  for (idxi = 0; idxi < nseq; idxi++) { nci0[idxi] = nci[idxi] = 0; }
  for (idxs = 0; idxs < nseq; idxs++) { ncs0[idxs] = ncs[idxs] = 0; }

  idxi = idxs = 0;
  nblock = nline = 0;
  while ( (status = esl_buffer_GetLine(bf, &p, &n)) == eslOK)
    {
      /* number of characters on this line */
      for (nillegal = 0, nc = 0, pos = 0; pos < n; pos++) {
	if (isspace(p[pos]) || isdigit(p[pos])) continue;
	if (! isalpha(p[pos]) && strchr("-*?.", p[pos]) == NULL) { nillegal++; continue; }
	nc++;
      }

      if (!nc) 		/* blank line? */
	{
	  if (idxi)   is_interleaved = FALSE; 
	  if (nline)  is_sequential  = FALSE;
	  continue;
	}

      if (nblock == 0) 
	nci0[idxi++] = nc;
      else 
	{
	  if      (idxi == 0) ncpb = nc;
	  else if (nc != ncpb) is_interleaved = FALSE; 
	  if (nillegal)        is_interleaved = FALSE; 
          nci[idxi++] += nc;
	}
      if (idxi == nseq) { idxi = 0; nblock++; ncpb = 0; }        /* advance to next putative block in interleaved format */

      if   (nline == 0) ncs0[idxs] = nc;
      else {
	if (nillegal) is_sequential = FALSE; 
	ncs[idxs] += nc; 
      }
      nline++;
      if (ncs0[idxs] + ncs[idxs] > alen) { idxs++; nline = 0; } /* advance to next sequence in sequential format */
    }

  if (is_interleaved == FALSE) printf("rejected interleaved\n");
  if (is_sequential  == FALSE) printf("rejected sequential\n");

  for (idxi = 0; idxi < nseq; idxi++)
    {
      printf("%3d %5d %5d %d\n", idxi, nci0[idxi], nci[idxi], nci[idxi] + nci0[idxi] - alen);
      printf("    %5d %5d %d\n",       ncs0[idxi], ncs[idxi], ncs[idxi] + ncs0[idxi] - alen);
    }


 ERROR:
  esl_fatal("piss off");
  return status;
}


int 
main(int argc, char **argv)
{
  char *filename = argv[1];
  ESL_BUFFER *bf = NULL;
  int   status;

  esl_buffer_Open(filename, NULL, &bf);
  phylip_autodetect_details(bf, NULL, NULL);
}
  
