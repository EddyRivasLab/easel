/* easel synth :: generating random sequences
 *
 * Usage:    easel synth <alphabet> <N> <L>
 * Example:  easel synth dna 100 10000       # sample 100 DNA sequences of length 10K
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"
#include "esl_vectorops.h"

static ESL_OPTIONS cmd_options[] = {
  /* name         type           default   env range      togs  reqs  incomp      help                                      docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL,     NULL, NULL, NULL, "help; show brief info on version and usage",          0 },
  { "-o",         eslARG_OUTFILE,  NULL, NULL, NULL,     NULL, NULL, NULL, "direct output data to file <f>",                      0 },
  { "--seed",     eslARG_INT,       "0", NULL,"n>=0",    NULL, NULL, NULL, "set random number generator seed to <n>",             0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
  
/* esl_cmd_synth()
 *   
 *   <topcmd> : argv[0] for the main call to `easel`; e.g. `easel` or `./miniapps/easel`
 *   <sub>    : ptr to ESL_SUBCMD struct for esl_cmd_synth, including .func|.subcmd="synth"|.nargs|.usage|.description
 *   <argc>   : # of args passed to subcommand; original argc minus whatever was skipped to get to the subcmd
 *   <argv>   : ptr to the start of the subcmd `synth` in cmdline args
 */
int
esl_cmd_synth(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, NULL);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQ         *sq      = NULL;
  int             atype   = esl_abc_EncodeType(esl_opt_GetArg(go, 1));
  int             N       = atoi(esl_opt_GetArg(go, 2));
  int             L       = atoi(esl_opt_GetArg(go, 3));
  char           *outfile = esl_opt_GetString(go, "-o");
  FILE           *ofp     = stdout;
  double         *fq      = NULL;
  int             outfmt  = eslSQFILE_FASTA;
  int             i;
  int             status;

  if (atype == eslUNKNOWN) esl_fatal("<alphatype> argument needs to be e.g. rna|dna|amino; not %s", esl_opt_GetArg(go, 1));
  if (N <= 0)              esl_fatal("<N> argument (number of seqs) is an integer > 0; not %s",     esl_opt_GetArg(go, 2));
  if (L <= 0)              esl_fatal("<L> argument (seq length) is an integer > 0; not %s",         esl_opt_GetArg(go, 3));

  abc = esl_alphabet_Create(atype);
  sq  = esl_sq_CreateDigital(abc);
  esl_sq_GrowTo(sq, L);

  if (outfile) 
    { if ((ofp = fopen(outfile, "w")) == NULL) esl_fatal("Failed to open output file %s\n", outfile); }

  /* Pick the iid frequency distribution to use */
  ESL_ALLOC(fq, sizeof(double) * abc->K);
  switch (atype) {
  case eslRNA:    esl_vec_DSet(fq, 4, 0.25); break;
  case eslDNA:    esl_vec_DSet(fq, 4, 0.25); break;
  case eslAMINO:  esl_composition_SW34(fq);  break;
  default:        esl_vec_DSet(fq, abc->K, 1.0 / (double) abc->K); break;
  }

  /* generate */
  for (i = 0; i < N; i++)
    {
      esl_rsq_xIID(rng, fq, abc->K, L, sq->dsq);
      if (N > 1) esl_sq_FormatName(sq, "random%d", i);
      else       esl_sq_SetName(sq, "random");
      sq->n = L;
      esl_sqio_Write(ofp, sq, outfmt, FALSE);
    }

  if (outfile) fclose(ofp);
  free(fq);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  esl_fatal("allocation failed");
}
  
  
