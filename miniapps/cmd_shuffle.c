/* `easel shuffle` : shuffling/randomizing sequences
 *
 * Usage:   easel shuffle <seqfile>
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"

#define SHUFOPTS "-d,-k,-0,-1,-r,-w"      // toggle group: alternative seq shuffling methods

static ESL_OPTIONS cmd_options[] = {
  /* name         type           default   env range      togs  reqs  incomp      help                                      docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL,     NULL, NULL, NULL, "help; show brief info on version and usage",          1 },
  { "-o",         eslARG_OUTFILE,  NULL, NULL, NULL,     NULL, NULL, NULL, "direct output data to file <f>",                      1 },
  { "--seed",     eslARG_INT,       "0", NULL,"n>=0",    NULL, NULL, NULL, "set random number generator seed to <n>",             1 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL,     NULL, NULL, NULL, "assert that input file is in format <s>",             1 },

  /* Choice of shuffling method */
  { "-d",         eslARG_NONE,    FALSE, NULL, NULL, SHUFOPTS, NULL, NULL, "shuffle preserving mono- and di-residue composition", 2 },
  { "-k",         eslARG_INT,     FALSE, NULL,"n>0", SHUFOPTS, NULL, NULL, "shuffle nonoverlapping <n>-mers",                     2 },
  { "-0",         eslARG_NONE,    FALSE, NULL, NULL, SHUFOPTS, NULL, NULL, "generate with 0th order Markov properties per input", 2 },
  { "-1",         eslARG_NONE,    FALSE, NULL, NULL, SHUFOPTS, NULL, NULL, "generate with 1st order Markov properties per input", 2 },
  { "-r",         eslARG_NONE,    FALSE, NULL, NULL, SHUFOPTS, NULL, NULL, "reverse each input",                                  2 },
  { "-w",         eslARG_INT,     FALSE, NULL,"n>0", SHUFOPTS, NULL, NULL, "regionally shuffle inputs in window size <n>",        2 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


/* There's multiple sections of options, so we provide a customized function
 * to esl_subcmd_CreateDefaultApp() for showing option help 
 */
static int
show_opthelp(const ESL_GETOPTS *go)
{
  if ( esl_printf("\nwhere general options are:\n")                                             != eslOK) return eslFAIL;
  if ( esl_opt_DisplayHelp(stdout, go, /*docgroup=*/1, /*indent=*/2, /*textwidth=*/80)          != eslOK) return eslFAIL;

  if ( esl_printf("\noptions for alternative shuffling methods (default is monoshuffling):\n")  != eslOK) return eslFAIL;
  if ( esl_opt_DisplayHelp(stdout, go, /*docgroup=*/2, /*indent=*/2, /*textwidth=*/80)          != eslOK) return eslFAIL;
  return eslOK;
}


/* esl_cmd_shuffle()
 *   
 *   <topcmd> : argv[0] for the main call to `easel`; e.g. `easel` or `./miniapps/easel`
 *   <sub>    : ptr to ESL_SUBCMD struct for esl_cmd_shuffle, including .func|.subcmd="shuffle"|.nargs|.usage|.description
 *   <argc>   : # of args passed to subcommand; original argc minus whatever was skipped to get to the subcmd
 *   <argv>   : ptr to the start of the subcmd `shuffle` in cmdline args
 */
int
esl_cmd_shuffle(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, &show_opthelp);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *outfile = esl_opt_GetString(go, "-o");
  FILE           *ofp     = stdout;
  int             outfmt  = eslSQFILE_FASTA;
  char           *seqfile = esl_opt_GetArg(go, 1);
  int             infmt   = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = esl_sq_Create();
  int             status;

  if (outfile) {
    if ((ofp = fopen(outfile, "w")) == NULL)
      esl_fatal("Failed to open output file %s\n", outfile);
  }

  if (esl_opt_GetString(go, "--informat") != NULL) {
    infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }

  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such seqfile");
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile couldn't be determined");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      if      (esl_opt_GetBoolean(go, "-d"))  esl_rsq_CShuffleDP     (rng, sq->seq, sq->seq);                                // diresidue shuffling 
      else if (esl_opt_GetBoolean(go, "-0"))  esl_rsq_CMarkov0       (rng, sq->seq, sq->seq);                                // 0th order Markov 
      else if (esl_opt_GetBoolean(go, "-1"))  esl_rsq_CMarkov1       (rng, sq->seq, sq->seq);                                // 1st order Markov
      else if (esl_opt_GetBoolean(go, "-r"))  esl_rsq_CReverse       (     sq->seq, sq->seq);                                // reverse 
      else if (esl_opt_IsOn      (go, "-w"))  esl_rsq_CShuffleWindows(rng, sq->seq, esl_opt_GetInteger(go, "-w"), sq->seq);  // monoshuffling in nonoverlapping windows
      else if (esl_opt_IsOn      (go, "-k"))  esl_rsq_CShuffleKmers  (rng, sq->seq, esl_opt_GetInteger(go, "-k"), sq->seq);  // shuffle nonoverlapping kmers
      else                                    esl_rsq_CShuffle       (rng, sq->seq, sq->seq);                                // default: monoresidue shuffling 

      esl_sq_FormatName(sq, "%s-shuffled", sq->name);
      esl_sqio_Write(ofp, sq, outfmt, FALSE);

      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					   sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    status, sqfp->filename);

  if (outfile) fclose(ofp);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
} 


