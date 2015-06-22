/* Translate DNA sequence into six frames, into individual ORFs.
 * 
 */
#include "esl_trans.h"

/*****************************************************************
 * 1. Functions for creating and using a stateful structure,
 * ESL_TRANS_WORKSTATE, to support both ReadSeq and ReadWindow()
 *****************************************************************/

void
esl_trans_WorkstateDestroy(ESL_TRANS_WORKSTATE *wrk)
{
  int f;
  if (wrk)
    {
      for (f = 0; f < 3; f++) esl_sq_Destroy(wrk->psq[f]);
      free(wrk);
    }
}

ESL_TRANS_WORKSTATE *
esl_trans_WorkstateCreate(ESL_GETOPTS *go, ESL_GENCODE *gcode)
{
  ESL_TRANS_WORKSTATE *wrk = NULL;
  int    f;
  int    status;

  ESL_ALLOC(wrk, sizeof(ESL_TRANS_WORKSTATE));
  for (f = 0; f < 3; f++) wrk->psq[f] = NULL;

  for (f = 0; f < 3; f++)
    {
      wrk->psq[f]         = esl_sq_CreateDigital(gcode->aa_abc);
      wrk->psq[f]->dsq[0] = eslDSQ_SENTINEL;
      wrk->in_orf[f]      = FALSE;
    }

  wrk->apos             = 1;
  wrk->frame            = 0;
  wrk->codon            = 0;
  wrk->inval            = 0;
  wrk->is_revcomp       = FALSE;
  wrk->orfcount         = 0;

  wrk->do_watson        = (esl_opt_GetBoolean(go, "--crick")  ? FALSE : TRUE);
  wrk->do_crick         = (esl_opt_GetBoolean(go, "--watson") ? FALSE : TRUE);
  wrk->using_initiators = ((esl_opt_GetBoolean(go, "-m") || esl_opt_GetBoolean(go, "-M")) ? TRUE : FALSE);
  wrk->minlen           = esl_opt_GetInteger(go, "-l");
  wrk->outfp            = stdout;
  wrk->outformat        = eslSQFILE_FASTA;

  return wrk;

 ERROR:
  esl_trans_WorkstateDestroy(wrk);
  return NULL;
}




/*****************************************************************
 *  2. Functions shared by the esl-translate's full or windowed reads style,
 *  as well as nhmmscant
 *****************************************************************/

int
esl_trans_ProcessOrf(ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq)
{
  ESL_SQ *psq = wrk->psq[wrk->frame];

  psq->end = (wrk->is_revcomp ? wrk->apos+1 : wrk->apos-1);

  if (wrk->in_orf[wrk->frame] && psq->n >= wrk->minlen)
    {
      wrk->orfcount++;
      if (psq->n+2 > psq->salloc) 
	esl_sq_Grow(psq, /*opt_nsafe=*/NULL);
      psq->dsq[1+psq->n] = eslDSQ_SENTINEL;
      
      esl_sq_FormatName(psq, "orf%d", wrk->orfcount);
      esl_sq_FormatDesc(psq, "source=%s coords=%d..%d length=%d frame=%d  %s", psq->source, psq->start, psq->end, psq->n, wrk->frame + 1 + (wrk->is_revcomp ? 3 : 0), sq->desc);
      esl_sqio_Write(wrk->outfp, psq, wrk->outformat, /*sq ssi offset update=*/FALSE);
    }

  esl_sq_Reuse(psq);
  esl_sq_SetSource(psq, sq->name);
  wrk->in_orf[wrk->frame] = FALSE;
  return eslOK;
}

void
esl_trans_ProcessStart(ESL_GENCODE *gcode, ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq)
{
  int f;

  ESL_DASSERT1(( sq->n >= 3 ));     

  for (f = 0; f < 3; f++)
    {
      esl_sq_SetSource(wrk->psq[f], sq->name);
      wrk->in_orf[f] = FALSE;
    }
  wrk->frame      = 0;
  wrk->codon      = 0;
  wrk->inval      = 0;
  wrk->is_revcomp = (sq->end > sq->start ? FALSE : TRUE  );   // this test fails for seqs of length 1, but we know that L>=3
  wrk->apos       = (wrk->is_revcomp ?     sq->L : 1     );

  if (esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[1])) wrk->codon += 4 * sq->dsq[1]; else wrk->inval = 1;
  if (esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[2])) wrk->codon +=     sq->dsq[2]; else wrk->inval = 2;
}


int
esl_trans_ProcessPiece(ESL_GENCODE *gcode, ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq)
{
  ESL_DSQ aa;
  int     rpos;

  for (rpos = 1; rpos <= sq->n-2; rpos++)
    {
      wrk->codon = (wrk->codon * 4) % 64;
      if   ( esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[rpos+2])) wrk->codon += sq->dsq[rpos+2]; 
      else wrk->inval = 3;

      /* Translate the current codon starting at <pos>;
       * see if it's an acceptable initiator 
       */
      if (wrk->inval > 0) // degenerate codon: needs special, tedious handling
      {
        aa =  esl_gencode_GetTranslation(gcode, sq->dsq+rpos);                         // This function can deal with any degeneracy
        if (! wrk->in_orf[wrk->frame] && esl_gencode_IsInitiator(gcode, sq->dsq+rpos)) //   ...as can IsInitiator.
          {
            if (wrk->using_initiators)  // If we're using initiation codons, initial codon translates to M even if it's something like UUG or CUG
              aa = esl_abc_DigitizeSymbol(gcode->aa_abc, 'M');
            wrk->in_orf[wrk->frame]     = TRUE;
            wrk->psq[wrk->frame]->start = wrk->apos;
          }
        wrk->inval--;
      }
      else          
      {
        aa = gcode->basic[wrk->codon];                             // If we know the digitized codon has no degeneracy, translation is a simple lookup
        if (gcode->is_initiator[wrk->codon] && ! wrk->in_orf[wrk->frame])
          {
            if (wrk->using_initiators)  // If we're using initiation codons, initial codon translates to M even if it's something like UUG or CUG
              aa = esl_abc_DigitizeSymbol(gcode->aa_abc, 'M');
            wrk->psq[wrk->frame]->start = wrk->apos;
            wrk->in_orf[wrk->frame]     = TRUE;
          }
      }

      /* Stop codon: deal with this ORF sequence and reinitiate */
      if ( esl_abc_XIsNonresidue(gcode->aa_abc, aa))
        esl_trans_ProcessOrf(wrk, sq);
      
      /* Otherwise: we have a residue. If we're in an orf (if we've
       * seen a suitable initiator), add this residue, reallocating as needed. 
       */
      if (wrk->in_orf[wrk->frame])
      {
        if (wrk->psq[wrk->frame]->n + 2 > wrk->psq[wrk->frame]->salloc)
          esl_sq_Grow(wrk->psq[wrk->frame], /*opt_nsafe=*/NULL);
        wrk->psq[wrk->frame]->dsq[1+ wrk->psq[wrk->frame]->n] = aa;
        wrk->psq[wrk->frame]->n++;
      }

      /* Advance +1 */
      if (wrk->is_revcomp) wrk->apos--; else wrk->apos++;
      wrk->frame = (wrk->frame + 1) % 3;
    }
  return eslOK;
}


int
esl_trans_ProcessEnd(ESL_TRANS_WORKSTATE *wrk, ESL_SQ *sq)
{
  int f;

  /* Done with the sequence. Now terminate all the orfs we were working on.
   * <apos> is sitting at L-1 (or 2, if revcomp) and we're in some <frame> 
   * there.
   */
  ESL_DASSERT1(( (wrk->is_revcomp && wrk->apos == 2) || (! wrk->is_revcomp && wrk->apos == sq->L-1) ));
  for (f = 0; f < 3; f++) // f counts 0..2, but it is *not* the <frame> index; <frame> is stateful
    {
      esl_trans_ProcessOrf(wrk, sq);
      if (wrk->is_revcomp) wrk->apos--; else wrk->apos++;
      wrk->frame = (wrk->frame + 1) % 3;
    }  
  return eslOK;
}

