/* Translate DNA sequence into six frames, into individual ORFs.
 * 
 */
#include "esl_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"


static ESL_OPTIONS options[] = {
  /* name           type        default  env  range toggles reqs incomp  help                                          docgroup*/
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "show brief help on version and usage",             0 },
  { "-l",         eslARG_INT,      "20", NULL, NULL, NULL,  NULL, NULL,  "minimum ORF length",                               0 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL,  NULL, NULL,  "specify that input file is in format <s>",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <nucleicseqfile>";
static char banner[] = "six-frame translation of nucleic acid seq to ORFs";

static int output_orf(ESL_GETOPTS *go, FILE *outfp, int outformat, ESL_SQ *psq, int *upd_orfcounter);

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go          = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *nt_abc      = esl_alphabet_Create(eslDNA);
  ESL_ALPHABET  *aa_abc      = esl_alphabet_Create(eslAMINO);
  ESL_GENCODE   *gcode       = NULL;
  ESL_SQFILE    *sqfp        = NULL;
  FILE          *outfp       = stdout;
  char          *seqfile     = esl_opt_GetArg(go, 1);
  int            informat    = eslSQFILE_UNKNOWN;
  int            outformat   = eslSQFILE_FASTA;
  int            windowsize  = 4092;                // can be any value, but a multiple of 3 makes most sense. windowsize can be +/-; + means reading top strand; - means bottom strand.
  int            contextsize = 2;                   // contextsize (adjacent window overlap) must be 2, or translation won't work properly.
  ESL_SQ        *sq;                                // DNA sequence window, digital
  ESL_SQ        *psq[3];                            // Current ORF sequences, for frames 0..2
  int8_t         in_orf[3];
  int            frame, f;        // 0..2, counter for which frame we're in.
  int32_t        codon;           // 0..63, digitally encoded canonical codon
  int            inval;           // counter for how many nt's we need to get past before we can have a fully canonical triplet again (no degeneracy syms)
  int            apos;            // absolute position in the DNA sequence, 1..L 
  int            rpos;            // "relative" position, 1..W in current window
  int            orfcounter = 0;  // >=1; number of ORFs output so far.
  ESL_DSQ        aa;              // amino acid residue, translated from current codon
  int            wstatus;         // return status specifically for DNA sequence window read
  int            status;
  
  ESL_DASSERT1(( windowsize  % 3 == 0 ));  

  if (esl_opt_GetString(go, "--informat") != NULL) 
    {
      informat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
      if (informat == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat", esl_opt_GetString(go, "--informat")); 
    }

  status = esl_sqfile_OpenDigital(nt_abc, seqfile, informat, /*env=*/NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to find (or open) sequence file %s", seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Failed to recognize format of sequence file %s", seqfile);
  else if (status != eslOK)        esl_fatal("Failure in opening sequence file %s; code %d", seqfile, status);


  /*  A default standard genetic code */
  gcode = esl_gencode_Create(nt_abc, aa_abc);

  sq    = esl_sq_CreateDigital(nt_abc);

  for (frame = 0; frame < 3; frame++)
    {
      psq[frame]         = esl_sq_CreateDigital(aa_abc);
      psq[frame]->dsq[0] = eslDSQ_SENTINEL;
      in_orf[frame]      = 0;
    }
  
  while (( wstatus = esl_sqio_ReadWindow(sqfp, contextsize, windowsize, sq)) != eslEOF)
    {
      if (wstatus == eslEOD)
	{
	  /* Terminate all orfs we were working on.
	   * <apos> is sitting at L-1 (or 2, if revcomp), and we're in some <frame> there
	   */
	  ESL_DASSERT1((  (windowsize > 0 && apos == sq->L-1) || (windowsize < 0 && apos == 2)  ));
	  for (f = 0; f < 3; f++)
	    {
	      if (in_orf[frame])
		{ 
		  psq[frame]->end = (windowsize > 0 ? apos-1 : apos+1);
		  output_orf(go, outfp, outformat, psq[frame], &orfcounter); 
		}
	      frame = (frame+1)%3;
	      apos  = (windowsize > 0 ? apos+1 : apos-1);
	    }
	  if (windowsize < 0) esl_sq_Reuse(sq); // Finished revcomp, so we're done with this seq.
	  windowsize = -windowsize;             // Switch either to next seq, or to revcomp of this one.
	  continue;
	}
      else if (wstatus == eslEFORMAT) esl_fatal("Parsing failed in sequence file %s:\n%s",          sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (wstatus == eslEINVAL)  esl_fatal("Invalid residue(s) found in sequence file %s\n%s", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (wstatus != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s", wstatus, sqfp->filename);

      /* If we're the first window in this input DNA sequence 
       * (or the first window in its revcomp), then initialize.
       * sq->C is the actual context overlap; 0=1st window; 2 (i.e. C)= subsequent.
       */
      if (sq->C == 0) 
	{
	  if (sq->n < 3) continue; // DNA sequence too short; skip it, don't even bother to revcomp, go to next sequence.
	  for (frame = 0; frame < 3; frame++)
	    {
	      esl_sq_SetSource(psq[frame], sq->name);
	      in_orf[frame] = FALSE;
	    }
	  frame = 0;
	  codon = 0;
	  inval = 0;
	  if (esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[1])) codon += 4 * sq->dsq[1]; else inval = 1;
	  if (esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[2])) codon +=     sq->dsq[2]; else inval = 2;
	  apos  = (windowsize > 0 ? 1 : sq->L);  // apos is "absolute position": i.e. position in sequence
	}

      /* Process the window, starting at position 1.
       * That's either the 1st nucleotide of the seq (or its revcomp) (if this is the first window);
       * or it's the 1st nucleotide of the 2nt context overlap.
       */
      for (rpos = 1; rpos <= sq->n-2; rpos++)  // rpos is "relative position": i.e. position in current window
	{
	  codon = (codon * 4) % 64;
	  if   ( esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[rpos+2])) codon += sq->dsq[rpos+2]; else inval = 3;
	  
	  /* Translate the codon. 
	   * See if it's an acceptable initiator.
	   */
	  if (inval > 0) 
	    {
	      aa = esl_gencode_TranslateCodon(gcode, sq->dsq+rpos);  // This function can deal with any degeneracy
	      inval--; 
	    }
	  else          
	    {
	      aa = gcode->basic[codon];                             // If we know the digitized codon has no degeneracy, translation is a simple lookup
	      if (gcode->is_initiator[codon] && ! in_orf[frame]) 
		{
		  aa = esl_abc_DigitizeSymbol(aa_abc, 'M');
		  psq[frame]->start = apos;
		  in_orf[frame] = TRUE;
		}
	    }

	  /* Stop codon: deal with this ORF sequence and reinitiate */
	  if ( esl_abc_XIsNonresidue(gcode->aa_abc, aa))
	    {
	      psq[frame]->end = (windowsize > 0 ? apos-4 : apos+4);     // stop codon itself doesn't count in coords.
	      if (in_orf[frame])
		output_orf(go, outfp, outformat, psq[frame], &orfcounter);  // orfcounter is bumped +1 by the call, and psq[frame] is Reuse'd.
	      in_orf[frame] = FALSE;
	    }

	  /* Otherwise: we have a residue. If we're in an orf (if we've
	   * seen a suitable initiator), add this residue, and reallocate
	   * as needed.
	   */
	  if (in_orf[frame]) 
	    {
	      esl_sq_Grow(psq[frame], /*opt_nsafe=*/NULL);
	      psq[frame]->dsq[1+psq[frame]->n] = aa;
	      psq[frame]->n++;
	    }

	  if (windowsize > 0) apos++; else apos--;
	  frame = (frame + 1) % 3;
	}
      
    }

  esl_sq_Destroy(sq);
  for (f = 0; f < 3; f++) esl_sq_Destroy(psq[f]);
  esl_sqfile_Close(sqfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_getopts_Destroy(go);
  return 0;
}


static int
output_orf(ESL_GETOPTS *go, FILE *outfp, int outformat, ESL_SQ *psq, int *upd_orfcounter)
{
  int minlen = esl_opt_GetInteger(go, "-l");

  if (psq->n >= minlen)
    {
      *upd_orfcounter = (*upd_orfcounter)+1;
      esl_sq_Grow(psq, /*opt_nsafe=*/NULL);   // Make sure we have room for the sentinel
      psq->dsq[1+psq->n] = eslDSQ_SENTINEL;   
      esl_sq_FormatName(psq, "orf%d", *upd_orfcounter);
      esl_sq_FormatDesc(psq, "source=%s coords=%d..%d length=%d", psq->source, psq->start, psq->end, psq->n);
      esl_sqio_Write(outfp, psq, outformat, /*sq ssi offset update=*/FALSE);
    }
  esl_sq_Reuse(psq);
  return eslOK;
}
