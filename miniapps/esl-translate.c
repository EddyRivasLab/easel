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

/*****************************************************************
 * 1. A stateful structure, workstate_s, to support both ReadSeq and ReadWindow()
 *****************************************************************/

/* struct workstate_s 
 *   keeps state in DNA sequence <sq>, allowing us to process a sequence
 *   either in a single gulp (using ReadSeq) or in overlapping windows 
 *   (using ReadWindow).
 *
 *   also contains one-time configuration information
 */
struct workstate_s {
  /* stateful info (which may get updated with each new seq, strand, and/or window): */
  ESL_SQ *psq[3];     // Growing ORFs in each frame
  int8_t  in_orf[3];  // TRUE|FALSE: TRUE if we're growing an ORF in this frame
  int     apos;       // 1..L:  current nucleotide we're on (starting a codon) in <sq>
  int     frame;      // 0..2:  which frame <apos> is in
  int     codon;      // 0..63: Digitized codon for apos,apos+1,apos+2
  int     inval;      // 0..3:  how many apos increments we need to get past an ambiguous nucleotide
  int     is_revcomp; // TRUE|FALSE: TRUE if we're doing reverse complement strand
  int     orfcount;   // >=0:   How many ORFs we've processed so far

  /* one-time configuration information (from options) */
  int     do_watson;         // TRUE|FALSE:  TRUE if we translate the top strand
  int     do_crick;          // TRUE|FALSE:  TRUE if we translate the reverse complement strand
  int     using_initiators;  // TRUE|FALSE : TRUE if -m or -M, only valid initiators can start an ORF, and initiator codon always translates to Met
  int     minlen;            // >=0: minimum orf length that process_orf will deal with
  FILE   *outfp;             // default stdout: where to write output ORF data
  int     outformat;         // default eslSQFILE_FASTA: sqfile format to write ORFs in
};

static void
workstate_destroy(struct workstate_s *wrk)
{
  int f;
  if (wrk)
    {
      for (f = 0; f < 3; f++) esl_sq_Destroy(wrk->psq[f]);
      free(wrk);
    }
}

static struct workstate_s *
workstate_create(ESL_GETOPTS *go, ESL_GENCODE *gcode)
{
  struct workstate_s *wrk = NULL;
  int    f;
  int    status;

  ESL_ALLOC(wrk, sizeof(struct workstate_s));
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
  workstate_destroy(wrk);
  return NULL;
}




/*****************************************************************
 * 2. Components shared by the two styles, full or windowed reads
 *****************************************************************/

static int
process_orf(struct workstate_s *wrk, ESL_SQ *sq)
{
  ESL_SQ *psq = wrk->psq[wrk->frame];

  psq->end = (wrk->is_revcomp ? wrk->apos+1 : wrk->apos-1);

  if (wrk->in_orf[wrk->frame] && psq->n >= wrk->minlen)
    {
      wrk->orfcount++;
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

static void
process_start(ESL_GENCODE *gcode, struct workstate_s *wrk, ESL_SQ *sq)
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


static int
process_piece(ESL_GENCODE *gcode, struct workstate_s *wrk, ESL_SQ *sq)
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
      if (wrk->inval > 0)
	{
	  aa =  esl_gencode_TranslateCodon(gcode, sq->dsq+rpos);  // This function can deal with any degeneracy
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
	process_orf(wrk, sq);  
      
      /* Otherwise: we have a residue. If we're in an orf (if we've
       * seen a suitable initiator), add this residue, reallocating as needed. 
       */
      if (wrk->in_orf[wrk->frame])
	{
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


static int
process_end(struct workstate_s *wrk, ESL_SQ *sq)
{
  int f;

  /* Done with the sequence. Now terminate all the orfs we were working on.
   * <apos> is sitting at L-1 (or 2, if revcomp) and we're in some <frame> 
   * there.
   */
  ESL_DASSERT1(( (wrk->is_revcomp && wrk->apos == 2) || (! wrk->is_revcomp && wrk->apos == sq->L-1) ));
  for (f = 0; f < 3; f++) // f counts 0..2, but it is *not* the <frame> index; <frame> is stateful
    {
      process_orf(wrk, sq);
      if (wrk->is_revcomp) wrk->apos--; else wrk->apos++;
      wrk->frame = (wrk->frame + 1) % 3;
    }  
  return eslOK;
}


/*****************************************************************
 * 3. Main loop for reading complete sequences with ReadSeq()
 *****************************************************************/

static int
do_by_sequences(ESL_GENCODE *gcode, struct workstate_s *wrk, ESL_SQFILE *sqfp)
{
  ESL_SQ *sq = esl_sq_CreateDigital(gcode->nt_abc);
  int     status;

  while (( status = esl_sqio_Read(sqfp, sq )) == eslOK)
    {
      if (sq->n < 3) continue;

      if (wrk->do_watson) {
	process_start(gcode, wrk, sq);
	process_piece(gcode, wrk, sq);
	process_end(wrk, sq);
      }

      if (wrk->do_crick) {
	esl_sq_ReverseComplement(sq);
	process_start(gcode, wrk, sq);
	process_piece(gcode, wrk, sq);
	process_end(wrk, sq);
      }

      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n",
					   sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					   status, sqfp->filename);
  
  esl_sq_Destroy(sq);
  return eslOK;
}


static int 
do_by_windows(ESL_GENCODE *gcode, struct workstate_s *wrk, ESL_SQFILE *sqfp)
{
  ESL_SQ *sq = esl_sq_CreateDigital(gcode->nt_abc);
  int     windowsize  = 4092;                // can be any value, but a multiple of 3 makes most sense. windowsize can be +/-; + means reading top strand; - means bottom strand.
  int     contextsize = 2;                   // contextsize (adjacent window overlap) must be 2, or translation won't work properly.
  int     wstatus;

  ESL_DASSERT1(( windowsize  % 3 == 0 ));  

  while (( wstatus = esl_sqio_ReadWindow(sqfp, contextsize, windowsize, sq)) != eslEOF)
    {
      if (wstatus == eslEOD)
	{
	  if ( (windowsize > 0 && wrk->do_watson) || (windowsize < 0 && wrk->do_crick))
	    process_end(wrk, sq);

	  if (windowsize > 0 && ! wrk->do_crick) { esl_sq_Reuse(sq); continue; } // Don't switch to revcomp if we don't need do. Allows -W --watson to work on nonrewindable streams
	  if (windowsize < 0) esl_sq_Reuse(sq);             // Do not Reuse the sq on the switch from watson to crick; ReadWindow needs sq->L
	  windowsize = -windowsize;                         // switch to other strand.
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
	  if ( (windowsize > 0 && wrk->do_watson) || (windowsize < 0 && wrk->do_crick))
	    process_start(gcode, wrk, sq);
	}

      if ( (windowsize > 0 && wrk->do_watson) || (windowsize < 0 && wrk->do_crick))      
	process_piece(gcode, wrk, sq);
    }
  esl_sq_Destroy(sq);
  return eslOK;
}


/*****************************************************************
 * x. main() for the esl-translate program
 *****************************************************************/

static ESL_OPTIONS options[] = {
  /* name           type        default  env  range toggles reqs incomp  help                                          docgroup*/
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "show brief help on version and usage",          0 },
  { "-c",         eslARG_INT,       "1", NULL, NULL, NULL,  NULL, NULL,  "use alt genetic code of NCBI transl table <n>", 0 },
  { "-l",         eslARG_INT,      "20", NULL, NULL, NULL,  NULL, NULL,  "minimum ORF length",                            0 },
  { "-m",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-M",  "ORFs must initiate with AUG only",              0 },
  { "-M",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-m",  "ORFs must start with allowed initiation codon", 0 },
  { "-W",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "use windowed, memory-efficient seq reading",    0 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL,  NULL, NULL,  "specify that input file is in format <s>",      0 },
  { "--watson",   eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate top strand",                     0 },
  { "--crick",    eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate bottom strand",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "six-frame translation of nucleic acid seq to ORFs";


int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go          = NULL;
  ESL_ALPHABET  *nt_abc      = esl_alphabet_Create(eslDNA);
  ESL_ALPHABET  *aa_abc      = esl_alphabet_Create(eslAMINO);
  ESL_GENCODE   *gcode       = NULL;
  struct workstate_s *wrk    = NULL;
  char          *seqfile     = NULL;
  int            informat    = eslSQFILE_UNKNOWN;
  ESL_SQFILE    *sqfp        = NULL;
  int            status;
  
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
  if ( esl_opt_GetBoolean(go, "-W") && ! esl_sqfile_IsRewindable(sqfp) && ! esl_opt_GetBoolean(go, "--watson"))
    esl_fatal("esl-translate can't read reverse complement from a nonrewindable stream (stdin pipe, .gz file, etc).");

  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any aa
   */
  gcode = esl_gencode_Create(nt_abc, aa_abc);
  esl_gencode_Set(gcode, esl_opt_GetInteger(go, "-c"));  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m or -M are set


  /* Set up the workstate structure, which contains both stateful 
   * info about our position in <sqfp> and the DNA <sq>, as well as
   * one-time config info from options
   */
  wrk = workstate_create(go, gcode);


  /* The two styles of main processing loop:
   */
  if (esl_opt_GetBoolean(go, "-W"))  do_by_windows(gcode, wrk, sqfp);
  else                               do_by_sequences(gcode, wrk, sqfp);


  workstate_destroy(wrk);
  esl_sqfile_Close(sqfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_getopts_Destroy(go);
  return 0;
}




