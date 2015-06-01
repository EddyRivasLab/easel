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
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "show brief help on version and usage",          0 },
  { "-c",         eslARG_INT,       "1", NULL, NULL, NULL,  NULL, NULL,  "use alt genetic code of NCBI transl table <n>", 0 },
  { "-l",         eslARG_INT,      "20", NULL, NULL, NULL,  NULL, NULL,  "minimum ORF length",                            0 },
  { "-m",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-M",  "ORFs must initiate with AUG only",              0 },
  { "-M",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-m",  "ORFs must start with allowed initiation codon", 0 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL,  NULL, NULL,  "specify that input file is in format <s>",      0 },
  { "--watson",   eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate top strand",                     0 },
  { "--crick",    eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate bottom strand",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "six-frame translation of nucleic acid seq to ORFs";

static int output_orf(ESL_GETOPTS *go, int on_topstrand, int frame, FILE *outfp, int outformat, ESL_SQ *sq, ESL_SQ *psq, int *upd_orfcounter);

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go          = NULL;
  ESL_ALPHABET  *nt_abc      = esl_alphabet_Create(eslDNA);
  ESL_ALPHABET  *aa_abc      = esl_alphabet_Create(eslAMINO);
  ESL_GENCODE   *gcode       = NULL;
  ESL_SQFILE    *sqfp        = NULL;
  FILE          *outfp       = stdout;
  char          *seqfile     = NULL;
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

  /***************************************************************** 
   * command line parsing
   *****************************************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  
  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, /*docgroup=*/0, /*indent=*/2, /*textwidth=*/80);

      puts("\nAvailable NCBI genetic code tables (for -c <id>):");
      esl_gencode_DumpCodeOptions(stdout);

      puts("\n<seqfile> must generally be a file, not a stdin pipe, because reverse");
      puts("complement requires file repositioning. If <seqfile> is - (i.e. signifying");
      puts("input from stdin), must also use --watson to restrict translation to top strand.");
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
 
  seqfile = esl_opt_GetArg(go, 1);

  if (esl_opt_IsOn(go, "--informat") &&
      (informat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid input sequence file format for --informat", esl_opt_GetString(go, "--informat")); 

  status = esl_sqfile_OpenDigital(nt_abc, seqfile, informat, /*env=*/NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to find (or open) sequence file %s", seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Failed to recognize format of sequence file %s", seqfile);
  else if (status != eslOK)        esl_fatal("Failure in opening sequence file %s; code %d", seqfile, status);

  /* A limitation. The esl_sqio_ReadWindow() interface needs to use SSI positioning
   * to read reverse complement, and that doesn't work on nonrewindable streams.
   */
  if (! esl_sqfile_IsRewindable(sqfp) && 
      ! esl_opt_GetBoolean(go, "--watson"))
    esl_fatal("esl-translate can't read reverse complement from a nonrewindable stream (stdin pipe, .gz file, etc).");

  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any aa
   */
  gcode = esl_gencode_Create(nt_abc, aa_abc);
  esl_gencode_Set(gcode, esl_opt_GetInteger(go, "-c"));  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m or -M are set

  /* Set up the sq (DNA) and psq (protein) objects
   */
  sq    = esl_sq_CreateDigital(nt_abc);
  for (frame = 0; frame < 3; frame++)
    {
      psq[frame]         = esl_sq_CreateDigital(aa_abc);
      psq[frame]->dsq[0] = eslDSQ_SENTINEL;
      in_orf[frame]      = 0;
    }
  
  /*****************************************************************
   * Main loop.
   * Read DNA sequence file, in windows (thus memory efficiently),
   * and translate as we go.
   *****************************************************************/
  while (( wstatus = esl_sqio_ReadWindow(sqfp, contextsize, windowsize, sq)) != eslEOF)
    {
      if (wstatus == eslEOD)
	{
	  /*  An edge case: it's possible to have L=0 for an empty sequence. */
	  if (sq->L == 0) { esl_sq_Reuse(sq); continue; }

	  /* Terminate all orfs we were working on.
	   * <apos> is sitting at L-1 (or 2, if revcomp), and we're in some <frame> there.
	   */
	  ESL_DASSERT1(( (windowsize > 0 && apos == sq->L-1) || (windowsize < 0 && apos == 2)  ));
	  for (f = 0; f < 3; f++)  // f counts 0..2; it is *not* the <frame> index, it's just going to make <frame> bump to frame, frame+1, frame+2
	    {
	      if (in_orf[frame])    
		{ 
		  psq[frame]->end = (windowsize > 0 ? apos-1 : apos+1);
		  output_orf(go, (windowsize > 0), frame, outfp, outformat, sq, psq[frame], &orfcounter); 
		  esl_sq_Reuse(psq[frame]);
		}
	      frame = (frame+1)%3;  //   ... and the <frame> index advances, wrapping around when needed.
	      apos  = (windowsize > 0 ? apos+1 : apos-1);
	    }
	  if (windowsize < 0) esl_sq_Reuse(sq); // Finished revcomp, so we're done with this DNA seq.

	  if (windowsize > 0 && esl_opt_GetBoolean(go, "--watson")) { esl_sq_Reuse(sq); continue; } // if we're only doing top strand, never use ReadWindow() on revcomp
	  if (windowsize < 0) esl_sq_Reuse(sq); // If we finished the crick strand, we're done with <sq> and we'll read a new one in next ReadWindow() call
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
		  if (esl_opt_GetBoolean(go, "-m") || esl_opt_GetBoolean(go, "-M"))  // If we're using initiation codons, initial codon translates to M even if it's something like UUG or CUG
		    aa = esl_abc_DigitizeSymbol(aa_abc, 'M');
		  psq[frame]->start = apos;
		  in_orf[frame] = TRUE;
		}
	    }

	  /* Stop codon: deal with this ORF sequence and reinitiate */
	  if ( esl_abc_XIsNonresidue(gcode->aa_abc, aa))
	    {
	      psq[frame]->end = (windowsize > 0 ? apos-1 : apos+1);           // stop codon itself doesn't count in coords. apos is on first nt of the stop codon, so back it up one.
	      if (in_orf[frame])
		{
		  output_orf(go, (windowsize > 0), frame, outfp, outformat, sq, psq[frame], &orfcounter);  // orfcounter is bumped +1 by the call if an orf is output
		  esl_sq_Reuse(psq[frame]);		  
		  esl_sq_SetSource(psq[frame], sq->name);
		}
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

#if 0
static int
do_by_sequences(ESL_SQFILE *sqfp, ESL_GENCODE *gcode)
{
  ESL_SQ *sq = esl_sq_CreateDigital(gcode->nt_abc);

  while (( status = esl_sqio_Read(sqfp, sq )) == eslOK)
    {
      


      esl_sq_ReverseComplement(sq);


      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n",
					   sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					   status, sqfp->filename);
  
  esl_sq_Destroy(sq);
}

static int
process_sequence(ESL_GENCODE *gcode, ESL_SQ *sq)
{
  int     codon = 0;
  int     inval = 0;
  int     frame = 0;
  int8_t  in_orf[3];
  ESL_DSQ aa;
  int     pos;
  int     f;

  for (f = 0; f < 3; f++)
    in_orf[f] = FALSE;

  for (pos = 1; pos <= sq->n-2; pos++)
    {
      codon = (codon * 4) % 64;
      if   ( esl_abc_XIsCanonical(gcode->nt_abc, sq->dsq[pos+2])) codon += sq->dsq[pos+2]; 
      else inval = 3;

      /* Translate the current codon starting at <pos>;
       * see if it's an acceptable initiator 
       */
      if (inval > 0)
	{
	  aa =  esl_gencode_TranslateCodon(gcode, sq->dsq+rpos);  // This function can deal with any degeneracy
	  inval--; 
	}
      else          
	{
	  aa = gcode->basic[codon];                             // If we know the digitized codon has no degeneracy, translation is a simple lookup
	  if (gcode->is_initiator[codon] && ! in_orf[frame]) 
	    {
	      if (esl_opt_GetBoolean(go, "-m") || esl_opt_GetBoolean(go, "-M"))  // If we're using initiation codons, initial codon translates to M even if it's something like UUG or CUG
		aa = esl_abc_DigitizeSymbol(aa_abc, 'M');
	      psq[frame]->start = apos;
	      in_orf[frame] = TRUE;
	    }
	}
	  

    }

}
#endif /*0*/


static int
output_orf(ESL_GETOPTS *go, int on_topstrand, int frame, FILE *outfp, int outformat, ESL_SQ *sq, ESL_SQ *psq, int *upd_orfcounter)
{
  int minlen    = esl_opt_GetInteger(go, "-l");
  int strand_ok = ((on_topstrand && ! esl_opt_GetBoolean(go, "--crick")) || (! on_topstrand && ! esl_opt_GetBoolean(go, "--watson")));

  if ( strand_ok && psq->n >= minlen)
    {
      *upd_orfcounter = (*upd_orfcounter)+1;
      esl_sq_Grow(psq, /*opt_nsafe=*/NULL);   // Make sure we have room for the sentinel
      psq->dsq[1+psq->n] = eslDSQ_SENTINEL;   
      esl_sq_FormatName(psq, "orf%d", *upd_orfcounter);
      esl_sq_FormatDesc(psq, "source=%s coords=%d..%d length=%d frame=%d  %s", psq->source, psq->start, psq->end, psq->n, frame + 1 + (on_topstrand ? 0 : 3), sq->desc);
      esl_sqio_Write(outfp, psq, outformat, /*sq ssi offset update=*/FALSE);
    }
  return eslOK;
}
