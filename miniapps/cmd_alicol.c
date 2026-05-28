/* `easel alicol` : remove/trim columns from an MSA
 *
 * Currently requires one of three options:
 *   easel alicol --mingap  : remove columns containing only gaps
 *   easel alicol --nogap   : remove columns containing *any* gaps (usually means deleting residues!)
 *   easel alicol --rfonly  : remove non-consensus columns, according to #=GC RF annotation
 *
 * To do:
 *     - when there's more operations than --mingap|--nogap, organize
 *       options into operations (exclusive), modifiers of operations,
 *       alphabet/format. For now, one of --mingap|--nogap is required.
 */
#include <esl_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_subcmd.h"

#define OPERATIONS "--mingap,--nogap,--rfonly"
#define ALPHOPTS   "--amino,--dna,--rna"

static ESL_OPTIONS cmd_options[] = {
  /* name         type         default  env   range togs  reqs        incomp       help                                                        docgroup */
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       NULL,        "help; show brief info on version and usage",                      0 },
  { "--mingap",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       OPERATIONS,  "remove columns that contain only gaps",                           0 },
  { "--nogap",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       OPERATIONS,  "remove columns containing any gaps",                              0 },
  { "--rfonly",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       OPERATIONS,  "remove nonconsensus columns, according to #=GC RF annotation",    0 },
  { "-o",         eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL,       NULL,        "send output MSA to file <f>, not stdout",                         0 },
  { "--keeprf",   eslARG_NONE,   FALSE, NULL, NULL, NULL, "--mingap", NULL,        "with --mingap, keep all consensus columns (nongaps in #=GC RF)",  0 },
  { "--amino",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       ALPHOPTS,    "assert <msafile> is protein (don't autodetect)",                  0 },
  { "--dna",      eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       ALPHOPTS,    "   ... <msafile> is DNA ...",                                     0 },
  { "--rna",      eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       ALPHOPTS,    "   ... <msafile> is RNA ...",                                     0 },
  { "--informat", eslARG_STRING, FALSE, NULL, NULL, NULL, NULL,       NULL,        "assert <msafile> is in format <s> (no autodetection)",            0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
esl_cmd_alicol(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go           = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, /*custom opthelp_f=*/NULL);
  ESL_ALPHABET   *abc          = NULL;
  char           *msafile      = esl_opt_GetArg(go, 1);
  ESL_MSAFILE    *afp          = NULL;
  ESL_MSA        *msa          = NULL;
  int             fmt          = eslMSAFILE_UNKNOWN;
  int             do_mingap    = esl_opt_GetBoolean(go, "--mingap");
  int             do_nogap     = esl_opt_GetBoolean(go, "--nogap");
  int             do_rfonly    = esl_opt_GetBoolean(go, "--rfonly");
  int             do_keeprf    = esl_opt_GetBoolean(go, "--keeprf");
  FILE           *ofp          = NULL;
  int             i;
  int             status;

  if (esl_opt_IsOn(go, "--informat") &&
      (fmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  /* Open output file, if we were given one */
  ofp = (esl_opt_GetString (go, "-o") == NULL ? stdout : fopen(esl_opt_GetString(go, "-o"), "w"));
  if (! ofp)  esl_fatal("Failed to open output file %s\n", esl_opt_GetString(go, "-o"));

  if (( status = esl_msafile_Open(&abc, msafile, /*env=*/NULL, fmt, /*fmtd=*/NULL, &afp)) != eslOK)
    esl_msafile_OpenFailure(afp, status);

  while ((status = esl_msafile_Read(afp, &msa)) != eslEOF)
    {
      if (status != eslOK) esl_msafile_ReadFailure(afp, status);

      if      (do_mingap) { if ( esl_msa_MinimGaps(msa, do_keeprf)             != eslOK) esl_fatal("Failed to remove all-gap columns from MSA"); }
      else if (do_nogap)  { if ( esl_msa_NoGaps   (msa, /*consider_rf=*/FALSE) != eslOK) esl_fatal("Failed to remove any-gap columns from MSA"); }
      else if (do_rfonly)
        {
          int *useme = NULL;
          if (! msa->rf) esl_fatal("--rfonly requires MSAs to have #=GC RF annotation");
          ESL_ALLOC(useme, sizeof(int) * msa->alen);

          for (i = 0; i < msa->alen; i++)
            useme[i] = (esl_abc_CIsGap(msa->abc, msa->rf[i]) || esl_abc_CIsMissing(msa->abc, msa->rf[i])) ? FALSE : TRUE;

          esl_msa_ColumnSubset(msa, useme);

          free(useme);
        }
      
      esl_msafile_Write(ofp, msa, afp->format);
      esl_msa_Destroy(msa);
    }

  if (ofp != stdout) fclose(ofp);
  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR: // UNREACHED
  esl_fatal("allocation failures were already fatal here");
}
