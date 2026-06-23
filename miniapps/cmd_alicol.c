/* `easel alicol` : select a subset of columns from an MSA
 *
 * 1. esl_cmd_alicol and command line parsing
 * 2. two code paths: default and --small
 * 3. selection options: creating masks
 * 4. selection options specialized for --small
 * 5. other internal functions
 */
#include <esl_config.h>

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_arr2.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile2.h"
#include "esl_regexp.h"
#include "esl_subcmd.h"
#include "esl_vectorops.h"

/*****************************************************************
 * 1. esl_cmd_alicol and command line parsing 
 *****************************************************************/

#define ALPHOPTS   "--amino,--dna,--rna"

static ESL_OPTIONS cmd_options[] = {
  /* name                 type     default  env  range      togs     reqs   incomp       help                                              docgroup */
  { "-h",             eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "show brief help information and exit",                     1 },
  { "-o",             eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,    NULL, "send output MSA(s) to file <f>, not stdout",               1 },
  { "-q",             eslARG_NONE,   FALSE, NULL, NULL,      NULL,    "-o",    NULL, "quiet: with -o, suppress verbose diagnostic output",       1 },
  { "--amino",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,ALPHOPTS, "assert <msafile> is protein (bypass autodetection)",       1 },
  { "--dna",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,ALPHOPTS, "   ... <msafile> is DNA ...",                              1 },
  { "--rna",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,ALPHOPTS, "   ... <msafile> is RNA ...",                              1 },
  { "--informat",     eslARG_STRING,  NULL, NULL, NULL,      NULL,    NULL,    NULL, "assert <msafile> is in format <s> (bypass autodetection)", 1 },
  { "--outformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,    NULL,"--small","write output MSAs in format <s> (default: same as input)", 1 },

  { "--gapfrac",      eslARG_REAL,    NULL, NULL, "0<x<=1",  NULL,    NULL,    NULL, "keep columns where gap char fraction is < x",              2 },
  { "--mask",         eslARG_INFILE,  NULL, NULL, NULL,      NULL,    NULL,    NULL, "specify columns to keep using mask string in file <f>",    2 },
  { "--mingap",       eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "remove columns that contain only gaps",                    2 },
  { "--nogap",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "remove columns containing any gaps",                       2 },
  { "--ppavg",        eslARG_REAL,    NULL, NULL, "0<=x<=1", NULL,    NULL,    NULL, "keep columns with mean individual post probs >= x",        2 },
  { "--ppcons",       eslARG_REAL,    NULL, NULL, "0<=x<=1", NULL,    NULL,    NULL, "keep columns with consensus post prob >= x",               2 },
  { "--ppfrac",       eslARG_REAL,    NULL, NULL, "0<=x<=1", NULL,    NULL,    NULL, "keep cols with >= x fraction of postprobs >= ppfracthresh",2 },
  { "--rfonly",       eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "keep consensus columns, remove nonconsensus columns",      2 },
  { "--span",         eslARG_STRING,  NULL, NULL, NULL,      NULL,    NULL,    NULL, "keep columns from '<pos1>:<pos2>', in 1..alen coords",     2 },
  { "--symfrac",      eslARG_REAL,    NULL, NULL, "0<x<=1",  NULL,    NULL,    NULL, "keep columns where residue char fraction is >= x",         2 },

  { "-c",             eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "use consensus coords 1..C (not alen) for --mask|--span",   3 },
  { "-w",             eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,"--small","use MSA weights when calculating fractions/averages",      3 },
  { "--keepins",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--ppcons",  NULL, "with --ppcons, keep all nonconsensus columns too",         3 },
  { "--keeprf",       eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    NULL, "keep all consensus columns (nongaps in #=GC RF)",          3 },
  { "--outmask",      eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,    NULL, "write the final all-column mask to file <f>",              3 },
  { "--outrfmask",    eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,    NULL, "write the final consensus column mask to file <f>",        3 },
  { "--ppfracthresh", eslARG_REAL,  "0.95", NULL, "0<=x<=1", NULL,"--ppfrac",  NULL, "sets postprob threshold for --ppfrac",                     3 },
  { "--small",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    "-w", "small memory mode",                                        3 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


typedef struct {
  int  alphatype;        // eslUNKNOWN to init. Need to know eslRNA|eslDNA MSAs.
  int  do_small;         // TRUE to use small-memory protocol
  
  /* Selection options. At least one is set. Selections are AND'ed. */
  int  do_gapfrac;
  int  do_maskfile;
  int  do_ppavg;
  int  do_ppcons;
  int  do_ppfrac;
  int  do_rfonly;
  int  do_span;
  int  do_symfrac;

  /* Options controlling selection options */
  float   gapfrac;         
  char   *maskfile;        
  float   ppavg;           
  float   ppcons;
  float   ppfrac;
  float   ppfracthresh;
  int64_t span_start;      // 1..alen, from user. Off by one from 0..alen-1 internal MSA and mask coords.
  int64_t span_end;        //  (ditto)
  float   symfrac;
  int     use_conscoords;
  int     use_weights;
  int     do_keepins;
  int     do_keeprf;

  /* Options controlling input/output */
  int   infmt;             // asserted input MSA file format, or eslMSAFILE_UNKNOWN. Actual file format after autodetection is in afp->format.
  int   outfmt;            // requested MSA output format from --outformat; or eslMSAFILE_UNKNOWN to default to same as input format.
  char *outfile;           // NULL, or output MSA file from -o. If NULL, a "verbose output" table is written to stdout instead.
  int   be_quiet;          // with output MSA going to a non-NULL outfile: don't write the "verbose output" table, write nothing to stdout.
  char *outmask_file;
  char *outrfmask_file;

} ESL_ALICOL_CFG;


/* 1. esl_cmd_alicol and command line parsing */
static int  show_opthelp       (const ESL_GETOPTS *go);
static void process_commandline(const ESL_GETOPTS *go, ESL_ALICOL_CFG *cfg);

/* 2. two code paths: default and --small */
static void default_protocol(const char *msafile, ESL_ALICOL_CFG *cfg, FILE *ofp);
static void small_protocol  (const char *msafile, ESL_ALICOL_CFG *cfg, FILE *ofp);

/* 3. selection options: creating masks */
static void mask_by_gapfrac (const ESL_MSA *msa, double gapfrac, int use_wgts, const int *rfmask, int *mask);
static void mask_by_symfrac (const ESL_MSA *msa, double symfrac, int use_wgts, const int *rfmask, int *mask);
static void mask_by_maskfile(int64_t alen, const char *maskfile, const int64_t *rfmap, int64_t C, int *mask);
static void mask_by_span    (int64_t alen, int64_t pos1, int64_t pos2, const int64_t *rfmap, int *mask);
static void mask_by_ppavg   (const ESL_MSA *msa, double ppavg,  int use_wgts,   const int *rfmask, int *mask);
static void mask_by_ppcons  (const ESL_MSA *msa, double ppcons, int do_keepins, const int *rfmask, int *mask);
static void mask_by_ppfrac  (const ESL_MSA *msa, double ppfrac, double ppfracthresh, int use_wgts, const int *rfmask, int *mask);

/* 4. selection options specialized for --small */
static void small_mask_by_gapfrac(double **abc_ct, int64_t alen, int K, double gapfrac,                     const int *rfmask, int *mask);
static void small_mask_by_symfrac(double **abc_ct, int64_t alen, int K, double symfrac,                     const int *rfmask, int *mask);
static void small_mask_by_ppavg  (double **pp_ct,  int64_t alen,        double ppavg,                       const int *rfmask, int *mask);
static void small_mask_by_ppfrac (double **pp_ct,  int64_t alen,        double ppfrac, double ppfracthresh, const int *rfmask, int *mask);

/* 5. other internal functions */
static void   create_rfmask(const char *rf, int64_t alen, int **ret_rfmask, int64_t *ret_C);
static void   create_rfmap(const int *rfmask, int64_t alen, int64_t C, int64_t **ret_rfmap);
static void   add_mask(const int *new_mask, int64_t alen, int *overall_mask);
static void   read_maskfile   (const char *filename, int **ret_mask, int64_t *ret_mlen);
static void   write_maskfile  (const char *filename, const int *mask, int64_t alen);
static void   write_rfmaskfile(const char *filename, const int *mask, const int64_t *rfmap, int64_t C);
static void   write_verbose_header(FILE *ofp);
static void   write_verbose_line  (FILE *ofp, const char *selection_option, const int *mask, int64_t alen, const int *rfmask, int64_t C);
static double pp_to_mean(char pp);
static double pp_to_min (char pp);

/*****************************************************************
 * 1. esl_cmd_alicol and command line parsing
 *****************************************************************/

/* esl_cmd_alicol()
 *
 * Follows standard interface for an Easel miniapp.
 */
int
esl_cmd_alicol(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, /*custom opthelp_f=*/&show_opthelp);
  char           *msafile = esl_opt_GetArg(go, 1);
  FILE           *ofp     = stdout;
  ESL_ALICOL_CFG  cfg;

  process_commandline(go, &cfg);

  if (cfg.outfile) {
    if ( (ofp = fopen(cfg.outfile, "w")) == NULL)
      esl_fatal("Failed to open output file %s", cfg.outfile);
  }

  if (cfg.do_small) small_protocol  (msafile, &cfg, ofp);
  else              default_protocol(msafile, &cfg, ofp);
  
  if (cfg.outfile) fclose(ofp);
  esl_getopts_Destroy(go);
  return eslOK;
}


/* show_opthelp()
 *
 * Formats the help display in sections.
 * 
 * Thrown errors are fatal, exiting the program with abnormal status 1.
 */
static int
show_opthelp(const ESL_GETOPTS *go)
{
  esl_printf("\nwhere general options are:\n");
  esl_opt_DisplayHelp(stdout, go, /*docgroup=*/1, /*indent=*/2, /*textwidth=*/80);

  esl_printf("\ncolumn selection options (at least one must be used):\n");
  esl_opt_DisplayHelp(stdout, go, /*docgroup=*/2, /*indent=*/2, /*textwidth=*/80);

  esl_printf("\nother options:\n");
  esl_opt_DisplayHelp(stdout, go, /*docgroup=*/3, /*indent=*/2, /*textwidth=*/80);
  return eslOK;
}


/* process_commandline()
 *
 * Parse the command line options into configuration structure <cfg>,
 * for passing into either default_protocol() or small_protocol().
 */
static void
process_commandline(const ESL_GETOPTS *go, ESL_ALICOL_CFG *cfg)
{
  int status;

  /* -h was already handled in esl_subcmd_CreateDefaultApp()
   */
  cfg->outfile  = esl_opt_GetString (go, "-o");
  cfg->be_quiet = esl_opt_GetBoolean(go, "-q");

  /* if cfg->alphatype stays eslUNKNOWN, caller will attempt
   * autodetection on input MSA file
   */
  if      (esl_opt_GetBoolean(go, "--rna"))   cfg->alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   cfg->alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) cfg->alphatype = eslAMINO;
  else                                        cfg->alphatype = eslUNKNOWN;

  /* if cfg->infmt stays eslMSAFILE_UNKNOWN, caller will attempt
   * autodetection
   */
  cfg->infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg->infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  cfg->outfmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--outformat") &&
      (cfg->outfmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--outformat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --outformat", esl_opt_GetString(go, "--outformat"));

  /* Selection options.
   */
  cfg->do_gapfrac = cfg->do_maskfile  = cfg->do_ppavg   = cfg->do_ppcons  = FALSE;
  cfg->do_ppfrac  = cfg->do_span      = cfg->do_symfrac = FALSE;
  cfg->gapfrac    = cfg->ppavg        = cfg->ppcons     = -1.0;
  cfg->ppfrac     = cfg->ppfracthresh = cfg->symfrac    = -1.0;
  cfg->maskfile   = NULL;

  cfg->do_rfonly = esl_opt_GetBoolean(go, "--rfonly");

  if ( esl_opt_IsOn(go, "--gapfrac")) { cfg->do_gapfrac  = TRUE; cfg->gapfrac  = esl_opt_GetReal  (go, "--gapfrac"); }
  if ( esl_opt_IsOn(go, "--mask"))    { cfg->do_maskfile = TRUE; cfg->maskfile = esl_opt_GetString(go, "--mask");    }
  if ( esl_opt_IsOn(go, "--mingap"))  { cfg->do_gapfrac  = TRUE; cfg->gapfrac  = 1.0; }
  if ( esl_opt_IsOn(go, "--nogap"))   { cfg->do_symfrac  = TRUE; cfg->symfrac  = 1.0; }
  if ( esl_opt_IsOn(go, "--symfrac")) { cfg->do_symfrac  = TRUE; cfg->symfrac  = esl_opt_GetReal  (go, "--symfrac"); }

  if ( esl_opt_IsOn(go, "--ppavg"))   { cfg->do_ppavg    = TRUE; cfg->ppavg    = esl_opt_GetReal  (go, "--ppavg");   }
  if ( esl_opt_IsOn(go, "--ppcons"))  { cfg->do_ppcons   = TRUE; cfg->ppcons   = esl_opt_GetReal  (go, "--ppcons");  }
  if ( esl_opt_IsOn(go, "--ppfrac"))  { cfg->do_ppfrac   = TRUE; cfg->ppfrac   = esl_opt_GetReal  (go, "--ppfrac");  }

  if ( esl_opt_IsOn(go, "--span")) {
    char *coordstring = esl_opt_GetString(go, "--span");
    status = esl_regexp_ParseCoordString(coordstring, &(cfg->span_start), &(cfg->span_end));
    if (status == eslESYNTAX)
      esl_fatal("For --span <s>, coords arg is '<pos1><sep><pos2>', like '42:100'; couldn't parse %s", coordstring);
    if ( (cfg->span_end != 0 && (cfg->span_start > cfg->span_end)) || cfg->span_start == 0 ) // ParseCoordString() allows revcomp coords, including "suffix"; we don't
      esl_fatal("For --span <s>, coords must be pos1 <= pos2 (if pos2 given); couldn't parse %s", coordstring);
    cfg->do_span = TRUE;
    // caller needs to validate <cfg->span_end> against alen once alen is known
  } else {
    cfg->span_start = -1;
    cfg->span_end   = -1;
  }

  if (! cfg->do_gapfrac && ! cfg->do_maskfile && ! cfg->do_ppavg && ! cfg->do_ppcons  &&
      ! cfg->do_ppfrac  && ! cfg->do_rfonly   && ! cfg->do_span  && ! cfg->do_symfrac)
    esl_fatal("At least one column selection option must be used.");

  cfg->use_conscoords = esl_opt_GetBoolean(go, "-c");              // requires --span || --mask; we'll check below
  cfg->use_weights    = esl_opt_GetBoolean(go, "-w");              // requires --gapfrac || --mingap || --nogap || --ppavg || --ppfrac || --symfrac
  cfg->ppfracthresh   = esl_opt_GetReal   (go, "--ppfracthresh");  // requires --ppfrac, getopts already checked
  cfg->do_keepins     = esl_opt_GetBoolean(go, "--keepins");       // requires --ppcons, getopts already checked
  cfg->do_keeprf      = esl_opt_GetBoolean(go, "--keeprf");        // requires a thresholded selection option; we'll check below
  cfg->outmask_file   = esl_opt_GetString (go, "--outmask");
  cfg->outrfmask_file = esl_opt_GetString (go, "--outrfmask");
  cfg->do_small       = esl_opt_GetBoolean(go, "--small");

  /* Easel's getopts doesn't let us do OR's for "one of these options
   * is required"; it only does AND's
   */
  if (cfg->use_conscoords && (! cfg->do_span && ! cfg->do_maskfile))
    esl_fatal("-c requires >=1 selection option that uses consensus coords: --span|--mask");
  if (cfg->use_weights    && (! esl_opt_IsOn(go, "--gapfrac") && ! cfg->do_ppavg && ! cfg->do_ppfrac && ! esl_opt_IsOn(go, "--symfrac")))
    esl_fatal("-w requires >=1 selection option that calculates fractions or averages");     // test directly for --gapfrac|--symfrac because --mingap|--nogap aren't accepted with -w
  if (cfg->do_keeprf && (! cfg->do_gapfrac && ! cfg->do_symfrac && ! cfg->do_ppfrac && ! cfg->do_ppavg))
    esl_fatal("--keeprf requires >=1 selection option that might delete consensus columns");

  /* Some pathological combinations that would remove nothing or
   * everything
   */
  if (cfg->do_keeprf && cfg->do_rfonly)
    esl_fatal("Using --rfonly with a thresholded selection + --keeprf is a no-op, would mask nothing");

  /* Additional --small requirements.
   *   - getopts already enforced no --outformat (must be eslMSAFILE_PFAM, one-block Stockholm)
   *   - getops also already enforced no -w. legacy esl_msafile2_RegurgitatePfam doesn't parse weights
   *   - input must also be eslMSAFILE_PFAM; if user asserted --informat, it must be that
   *   - count arrays require a known alphabet, so an alphabet must be asserted to use --small
   *   - we make two passes over the file, so it must be a real file, not stdin stream (checked later)
   */
  if (cfg->do_small)
    {
      if (cfg->infmt != eslMSAFILE_UNKNOWN && cfg->infmt != eslMSAFILE_PFAM)
        esl_fatal("--small requires Pfam (one-block Stockholm) format input");
      cfg->infmt = eslMSAFILE_PFAM;

      if (cfg->alphatype == eslUNKNOWN)
        esl_fatal("--small requires alphabet to be asserted with --amino, --dna, or --rna");
    }

  // Some options require particular input formats e.g. Stockholm. These are checked after we open the MSA file and see its format.
  // Several options require consensus column annotation. These are checked after we read the MSA.
}



/*****************************************************************
 * 2. two code paths: default and --small
 *****************************************************************/

static void
default_protocol(const char *msafile, ESL_ALICOL_CFG *cfg, FILE *ofp)
{
  ESL_MSAFILE *afp          = NULL;
  ESL_MSA     *msa          = NULL;
  int          nmsa         = 0;
  int         *onemask      = NULL;
  int         *overall_mask = NULL;
  int         *rfmask       = NULL;   // [0..alen-1] of 1:0 flags for consensus columns, or NULL
  int64_t     *rfmap        = NULL;   // [0..C-1] of coords 0..alen-1, mapping consensus cols to overall 0..alen coords
  int64_t      C;                     // consensus length (if msa->rf, else 0)
  int64_t      i;
  int          status;

  /* Open <msafile> in text mode.
   * cfg->alphatype is only used for knowing when to fix RNA basepair annotation
   */
  if (( status = esl_msafile_Open(/*abc=*/NULL, msafile, /*env=*/NULL, cfg->infmt, /*fmtd=*/NULL, &afp)) != eslOK)
    esl_msafile_OpenFailure(afp, status);

  /* Autodetect the alphabet, even though we're going to read in text
   * model. If this is an RNA MSA with secondary structure annotation,
   * we'll need to worry about fixing any base pairs broken by columns
   * downselection.
   */
  if (cfg->alphatype == eslUNKNOWN) {
    status = esl_msafile_GuessAlphabet(afp, &(cfg->alphatype));
    if (status != eslOK && status != eslENOALPHABET)
      esl_fatal("alphabet autodetection failed for unexpected reasons");
    // if we fail to guess alphabet, we leave alphatype == eslUNKNOWN; we're in text mode so that's fine
  }

  /* Check for format requirements
   */
  if ( (cfg->do_ppavg || cfg->do_ppcons || cfg->do_ppfrac) &&
       (afp->format != eslMSAFILE_STOCKHOLM && afp->format != eslMSAFILE_PFAM))
    esl_fatal("To use --ppavg|--ppcons|--ppfrac, input MSA file must be in Stockholm format"); // because only Stockholm has PP markups

  /* Default output format is the same as the input
   */
  cfg->outfmt = (cfg->outfmt == eslMSAFILE_UNKNOWN ? afp->format : cfg->outfmt);

  /* Main loop: all MSAs in file (though usually there would be just one)
   */
  while ((status = esl_msafile_Read(afp, &msa)) != eslEOF)
    {
      if (status != eslOK) esl_msafile_ReadFailure(afp, status);
      nmsa++;

      /* Multi-MSA checks: some selection options can't apply to >1 MSA, and only some 
       * MSA file formats allow >1 MSA.
       */
      if (nmsa > 1)
        {
          if (cfg->do_span)        esl_fatal("--span can only be applied to a single MSA; this msafile has >1 MSAs");
          if (cfg->do_maskfile)    esl_fatal("--mask can only be applied to a single MSA; this msafile has >1 MSAs");
          if (cfg->outmask_file)   esl_fatal("--outmask can only be applied to a single MSA; msafile has >1 MSAs");
          if (cfg->outrfmask_file) esl_fatal("--outrfmask can only be applied to a single MSA; msafile has >1 MSAs");
          if (! esl_msafile_IsMultiRecord(cfg->outfmt))
            esl_fatal("your output MSA file format (--outformat) doesn't allow >1 MSA per file");
        }

      /* Check for consensus column annotation, for options that need it
       */
      if (! msa->rf && cfg->do_rfonly)      esl_fatal("--rfonly selection requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_ppcons)      esl_fatal("--ppcons selection requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->use_conscoords) esl_fatal("-c requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_keeprf)      esl_fatal("--keeprf requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_keepins)     esl_fatal("--keepins requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->outrfmask_file) esl_fatal("--outrfmask requires consensus column annotation on the MSA");

      /* Check for posterior probability annotation, for options that need it
       */
      if (cfg->do_ppavg || cfg->do_ppfrac)
        {
          if (! msa->pp)                    esl_fatal("--ppavg|--ppfrac require posterior probability annotation (#=GR PP) for the MSA");
          for (i = 0; i < msa->nseq; i++)
            if (! msa->pp[i])               esl_fatal("--ppavg|--ppfrac require post prob annotation (#=GR PP) for every sequence");
        }
      if (cfg->do_ppcons && ! msa->pp_cons) esl_fatal("--ppcons requires consensus post prob annotation (#=GC PP_cons) for MSA");

      /* If using consensus column information, create <rfmask> and <rfmap>
       */
      if (msa->rf)
        {
          create_rfmask(msa->rf, msa->alen,    &rfmask, &C);
          create_rfmap (rfmask,  msa->alen, C, &rfmap);
        }
      else C = 0;

      /* Check --span against <alen> or <C>, now that both are known
       */
      if (cfg->do_span)
        {
          if (cfg->use_conscoords) {
            if      (cfg->span_end == 0) cfg->span_end = C;
            else if (cfg->span_end > C)  esl_fatal("--span end consensus coord %" PRId64 " is too large (C=%" PRId64 ")", cfg->span_end, C);  
          } else {
            if      (cfg->span_end == 0)        cfg->span_end = msa->alen;   // suffix rule for coordstring like "42:"
            else if (cfg->span_end > msa->alen) esl_fatal("--span end coord %" PRId64 " is off end of MSA (alen=%" PRId64 ")", cfg->span_end, msa->alen);
          }
        }

      /* Allocation and initialization of masks */
      ESL_ALLOC(overall_mask, sizeof(int) * msa->alen);
      ESL_ALLOC(onemask,      sizeof(int) * msa->alen);
      esl_vec_ISet(overall_mask, msa->alen, TRUE);

      /* Verbose output table header for this MSA, if MSA is going to an outfile
       */
      if (cfg->outfile && ! cfg->be_quiet)
        {
          if (msa->name) esl_printf("MSA #%d : %s\n", nmsa, msa->name);
          else           esl_printf("MSA %d\n",       nmsa);
          write_verbose_header(stdout);
        }

      /* Selection options: all masks are AND'ed together
       */
      if (cfg->do_gapfrac) {
        mask_by_gapfrac(msa, cfg->gapfrac, cfg->use_weights, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "gapfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_symfrac) {
        mask_by_symfrac(msa, cfg->symfrac, cfg->use_weights, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "symfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_maskfile) {
        mask_by_maskfile(msa->alen, cfg->maskfile, (cfg->use_conscoords ? rfmap : NULL), C, onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "maskfile", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_span) {
        mask_by_span(msa->alen, cfg->span_start, cfg->span_end, (cfg->use_conscoords ? rfmap : NULL), onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "span", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppavg) {
        mask_by_ppavg(msa, cfg->ppavg, cfg->use_weights, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "ppavg", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppcons) {
        mask_by_ppcons(msa, cfg->ppcons, cfg->do_keepins, rfmask, onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "ppcons", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppfrac) {
        mask_by_ppfrac(msa, cfg->ppfrac, cfg->ppfracthresh, cfg->use_weights, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "ppfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_rfonly) {
        add_mask(rfmask, msa->alen, overall_mask);
        if (cfg->outfile && ! cfg->be_quiet)
          write_verbose_line(stdout, "rfonly", rfmask, msa->alen, rfmask, C);
      }

      /* We have our final <overall_mask>.
       * Take the subset of columns and write output(s).
       */
      if (cfg->alphatype == eslRNA || cfg->alphatype == eslDNA)
        esl_msa_RemoveBrokenBasepairs(msa, overall_mask);

      /* Write the verbose output and the masks *before* we do column
       * selection, while we still have original MSA and its <alen>
       */
      if (cfg->outfile && ! cfg->be_quiet)
        write_verbose_line(stdout, NULL, overall_mask, msa->alen, rfmask, C);  // NULL makes it write the final summary

      if (cfg->outmask_file)
        write_maskfile(cfg->outmask_file, overall_mask, msa->alen);

      if (cfg->outrfmask_file)
        write_rfmaskfile(cfg->outrfmask_file, overall_mask, rfmap, C);


      /* Finally, downselect and write.
       */
      esl_msa_ColumnSubset(msa, overall_mask);
      esl_msafile_Write(ofp, msa, cfg->outfmt);

      free(rfmask);         rfmask       = NULL;
      free(rfmap);          rfmap        = NULL;
      free(onemask);        onemask      = NULL;
      free(overall_mask);   overall_mask = NULL;
      esl_msa_Destroy(msa); msa          = NULL;
    }
  esl_msafile_Close(afp);
  return;

 ERROR:
  esl_fatal("allocation failed"); //NOTREACHED
}



/* small_protocol()
 *
 * Small-memory flow: never read a whole MSA into memory. Instead make
 * two passes over the file with the legacy ESL_MSAFILE2
 * RegurgitatePfam() parser.
 *
 *   Pass 1: esl_msafile2_ReadInfoPfam() reads each MSA's non-sequence
 *           info (comments, #=GF, #=GC including RF and PP_cons) into an
 *           "info-only" MSA, plus optional per-column count arrays
 *           (abc_ct for residues/gaps, pp_ct for posterior probs). From
 *           these we compute and store each MSA's final column mask.
 *
 *   Pass 2: esl_msafile2_RegurgitatePfam() streams each MSA back out,
 *           deleting masked columns from aseqs/#=GR/#=GC as it goes.
 *           For RNA|DNA MSAs with RNA secondary structure annotation,
 *           we fix any broken base pairs in SS_cons/#=GR SS annotation.
 */
static void
small_protocol(const char *msafile, ESL_ALICOL_CFG *cfg, FILE *ofp)
{
  ESL_ALPHABET *abc        = esl_alphabet_Create(cfg->alphatype);  // asserted; used for count arrays and for SS basepair fixing in pass 2
  ESL_MSAFILE2 *afp2       = NULL;
  ESL_MSA      *msa        = NULL;                  // info-only MSA (GF|GC|comments), one at a time
  int         **masks      = NULL;                  // masks[0..nmsa-1][0..alen-1] : final column masks, stored across passes
  int64_t      *alens      = NULL;                  // alens[0..nmsa-1] : alen of each MSA (masks are ragged)
  int           msa_nalloc = 1;                     // <masks>,<alens> grown by doubling from here
  int          *onemask    = NULL;
  int          *rfmask     = NULL;                  // [0..alen-1] consensus column flags, or NULL
  int64_t      *rfmap      = NULL;                  // [0..C-1] consensus coords in 0..alen-1
  int64_t       C;                                  // consensus length, or 0
  double      **abc_ct     = NULL;                  // [0..alen-1][0..K] residue/gap counts, or NULL
  double      **pp_ct      = NULL;                  // [0..alen-1][0..11] posterior prob counts, or NULL
  int           need_abc_ct = (cfg->do_gapfrac || cfg->do_symfrac);
  int           need_pp_ct  = (cfg->do_ppavg   || cfg->do_ppfrac);
  int           nseq;
  int64_t       alen;
  int           nmsa       = 0;
  int           ai;
  int64_t       i;
  int           status;

  /* The two passes can't rewind a stdin stream, so the input must be a real file. */
  if (strcmp(msafile, "-") == 0)
    esl_fatal("--small can't read from a stdin pipe (-); it needs a file it can read twice");

  ESL_ALLOC(masks, sizeof(int *)   * msa_nalloc);
  ESL_ALLOC(alens, sizeof(int64_t) * msa_nalloc);

  /* Pass 1: read info + counts for each MSA, compute and store its mask.
   */
  status = esl_msafile2_Open(msafile, /*env=*/NULL, &afp2);
  if      (status == eslENOTFOUND) esl_fatal("MSA file %s doesn't exist or is not readable", msafile);
  else if (status == eslEFORMAT)   esl_fatal("Failed to parse MSA file %s : needs to be Pfam (one-block Stockholm) format", msafile);
  else if (status != eslOK)        esl_fatal("MSA file %s open failed with error code %d", msafile, status);

  while ((status = esl_msafile2_ReadInfoPfam(afp2, /*listfp=*/NULL, abc, /*known_alen=*/-1, /*known_rf=*/NULL, /*known_ss_cons=*/NULL,
                                             &msa, &nseq, &alen, /*opt_ngs=*/NULL, /*opt_maxname=*/NULL, /*opt_maxgf=*/NULL, /*opt_maxgc=*/NULL, /*opt_maxgr=*/NULL,
                                             (need_abc_ct ? &abc_ct : NULL),
                                             (need_pp_ct  ? &pp_ct  : NULL),
                                             /*opt_bp_ct=*/NULL, /*opt_spos_ct=*/NULL, /*opt_epos_ct=*/NULL)) != eslEOF)
    {
      if      (status == eslEFORMAT) esl_fatal("MSA file %s parse failed:\n%s", msafile, afp2->errbuf);
      else if (status != eslOK)      esl_fatal("MSA file %s read failed with error code %d", msafile, status);

      msa->alen = alen;   // ReadInfoPfam doesn't set alen on the info-only MSA; we need it for mask_by_span/maskfile

      /* Multi-MSA checks: some selection options can't apply to >1 MSA */
      if (nmsa >= 1)   // i.e. this is the 2nd or later MSA (nmsa is the count stored so far)
        {
          if (cfg->do_span)        esl_fatal("--span can only be applied to a single MSA; this msafile has >1 MSAs");
          if (cfg->do_maskfile)    esl_fatal("--mask can only be applied to a single MSA; this msafile has >1 MSAs");
          if (cfg->outmask_file)   esl_fatal("--outmask can only be applied to a single MSA; msafile has >1 MSAs");
          if (cfg->outrfmask_file) esl_fatal("--outrfmask can only be applied to a single MSA; msafile has >1 MSAs");
        }

      /* Annotation requirements, same as the default protocol */
      if (! msa->rf && cfg->do_rfonly)      esl_fatal("--rfonly selection requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_ppcons)      esl_fatal("--ppcons selection requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->use_conscoords) esl_fatal("-c requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_keeprf)      esl_fatal("--keeprf requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->do_keepins)     esl_fatal("--keepins requires consensus column annotation on the MSA");
      if (! msa->rf && cfg->outrfmask_file) esl_fatal("--outrfmask requires consensus column annotation on the MSA");
      if (cfg->do_ppcons && ! msa->pp_cons) esl_fatal("--ppcons requires consensus post prob annotation (#=GC PP_cons) for MSA");
      if (need_pp_ct && nseq > 0 && esl_vec_DSum(pp_ct[0], 12) == 0.0)  // pp_ct sums to nseq per column if PP present; 0 means no #=GR PP at all
        esl_fatal("--ppavg|--ppfrac require posterior probability annotation (#=GR PP) for the MSA");

      /* Consensus column structures */
      if (msa->rf)
        {
          create_rfmask(msa->rf, msa->alen,    &rfmask, &C);
          create_rfmap (rfmask,  msa->alen, C, &rfmap);
        }
      else C = 0;

      /* Validate --span against alen or C, now known */
      if (cfg->do_span)
        {
          if (cfg->use_conscoords) {
            if      (cfg->span_end == 0) cfg->span_end = C;
            else if (cfg->span_end > C)  esl_fatal("--span end consensus coord %" PRId64 " is too large (C=%" PRId64 ")", cfg->span_end, C);
          } else {
            if      (cfg->span_end == 0)        cfg->span_end = msa->alen;
            else if (cfg->span_end > msa->alen) esl_fatal("--span end coord %" PRId64 " is off end of MSA (alen=%" PRId64 ")", cfg->span_end, msa->alen);
          }
        }

      /* Allocate this MSA's masks; initialize overall mask to keep-all */
      ESL_ALLOC(onemask,     sizeof(int) * msa->alen);
      ESL_ALLOC(masks[nmsa], sizeof(int) * msa->alen);
      esl_vec_ISet(masks[nmsa], msa->alen, TRUE);

      /* Verbose output table header for this MSA, if MSA is going to an outfile */
      if (cfg->outfile && ! cfg->be_quiet)
        {
          if (msa->name) esl_printf("MSA #%d : %s\n", nmsa+1, msa->name);
          else           esl_printf("MSA %d\n",       nmsa+1);
          write_verbose_header(stdout);
        }

      /* Selection options: AND each mask into masks[nmsa] */
      if (cfg->do_gapfrac) {
        small_mask_by_gapfrac(abc_ct, msa->alen, abc->K, cfg->gapfrac, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "gapfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_symfrac) {
        small_mask_by_symfrac(abc_ct, msa->alen, abc->K, cfg->symfrac, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "symfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_maskfile) {
        mask_by_maskfile(msa->alen, cfg->maskfile, (cfg->use_conscoords ? rfmap : NULL), C, onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "maskfile", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_span) {
        mask_by_span(msa->alen, cfg->span_start, cfg->span_end, (cfg->use_conscoords ? rfmap : NULL), onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "span", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppavg) {
        small_mask_by_ppavg(pp_ct, msa->alen, cfg->ppavg, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "ppavg", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppcons) {
        mask_by_ppcons(msa, cfg->ppcons, cfg->do_keepins, rfmask, onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "ppcons", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_ppfrac) {
        small_mask_by_ppfrac(pp_ct, msa->alen, cfg->ppfrac, cfg->ppfracthresh, (cfg->do_keeprf ? rfmask : NULL), onemask);
        add_mask(onemask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "ppfrac", onemask, msa->alen, rfmask, C);
      }
      if (cfg->do_rfonly) {
        add_mask(rfmask, msa->alen, masks[nmsa]);
        if (cfg->outfile && ! cfg->be_quiet) write_verbose_line(stdout, "rfonly", rfmask, msa->alen, rfmask, C);
      }

      /* Final summary line + mask file output, before we lose this MSA's alen */
      if (cfg->outfile && ! cfg->be_quiet)
        write_verbose_line(stdout, NULL, masks[nmsa], msa->alen, rfmask, C);
      if (cfg->outmask_file)
        write_maskfile(cfg->outmask_file, masks[nmsa], msa->alen);
      if (cfg->outrfmask_file)
        write_rfmaskfile(cfg->outrfmask_file, masks[nmsa], rfmap, C);

      alens[nmsa] = msa->alen;
      nmsa++;
      if (nmsa == msa_nalloc) {
        msa_nalloc *= 2;
        ESL_REALLOC(masks, sizeof(int *)   * msa_nalloc);
        ESL_REALLOC(alens, sizeof(int64_t) * msa_nalloc);
      }

      esl_arr2_Destroy((void **) abc_ct, msa->alen); abc_ct = NULL;
      esl_arr2_Destroy((void **) pp_ct,  msa->alen); pp_ct  = NULL;
      free(onemask);        onemask = NULL;
      free(rfmask);         rfmask  = NULL;
      free(rfmap);          rfmap   = NULL;
      esl_msa_Destroy(msa); msa     = NULL;
    }
  esl_msafile2_Close(afp2); afp2 = NULL;

  if (nmsa == 0) esl_fatal("Failed to read any MSAs from file %s", msafile);

  /* Pass 2: reopen the file, regurgitate each MSA
   * with its stored mask, deleting columns as it streams.
   *
   * Legacy esl_msafile2_RegurgitatePfam() uses afp2->abc, but only to
   * check for whether to fix broken basepairs in RNA secondary
   * structure annotation. It actually does its reading in text mode.
   * To make this look slightly less confusing, we open in text mode
   * and set afp2->abc without setting afp2->do_digital.
   */
  if ((status = esl_msafile2_Open(msafile, /*env=*/NULL, &afp2)) != eslOK)
    esl_fatal("MSA file %s reopen failed on second pass with error code %d", msafile, status);
  afp2->abc = abc;

  for (ai = 0; ai < nmsa; ai++)
    {
      status = esl_msafile2_RegurgitatePfam(afp2, ofp,
                                            /*maxname=*/-1, /*maxgf=*/-1, /*maxgc=*/-1, /*maxgr=*/-1,  // use margins from the file
                                            TRUE, TRUE, TRUE, TRUE, TRUE,  // header, trailer, blanks, comments, GF
                                            TRUE, TRUE, TRUE, TRUE,        // GS, GC, GR, aseq
                                            /*seqs2regurg=*/NULL, /*seqs2skip=*/NULL,
                                            masks[ai],                     // useme: which columns to keep
                                            /*add2me=*/NULL,
                                            alens[ai],                     // expected alen
                                            /*gapchar2add=*/'.',
                                            /*opt_nseq_read=*/NULL, /*opt_nseq_written=*/NULL);
      if      (status == eslEOF)     esl_fatal("MSA file %s: ran out of MSAs on second pass (expected %d)", msafile, nmsa);
      else if (status == eslEFORMAT) esl_fatal("MSA file %s regurgitation failed at MSA %d:\n%s", msafile, ai+1, afp2->errbuf);
      else if (status != eslOK)      esl_fatal("MSA file %s regurgitation failed at MSA %d with error code %d", msafile, ai+1, status);
    }
  esl_msafile2_Close(afp2);

  for (i = 0; i < nmsa; i++) free(masks[i]);
  free(masks);
  free(alens);
  esl_alphabet_Destroy(abc);
  return;

 ERROR:
  esl_fatal("allocation failed"); //NOTREACHED
}



/*****************************************************************
 * 3. selection options: creating masks
 *****************************************************************/

/* mask_by_gapfrac()
 *
 * msa      : the multiple sequence alignment
 * gapfrac  : keep columns with < this fraction of gaps
 * use_wgts : TRUE to calculate sequence-weighted fractions (MSA has weight annotations)
 * rfmask   : if non-NULL, implement do_keeprf: [0..alen-1] 0|1 flags marking consensus columns
 * mask     : RESULT, caller provided space : mask[0..alen-1] 1:0 flags marking cols to keep
 *
 */
static void
mask_by_gapfrac(const ESL_MSA *msa, double gapfrac, int use_wgts, const int *rfmask, int *mask)
{
  double   ntot = (use_wgts ? esl_vec_DSum(msa->wgt, msa->nseq) : msa->nseq);
  double  *ngap = NULL;
  int64_t  s,i;
  int      status;

  ESL_ALLOC(ngap, sizeof(double) * msa->alen);
  esl_vec_DSet(ngap, msa->alen, 0.);

  /* collect <ngap> for all columns efficiently, in stride with aseq[s][i] */
  for (s = 0; s < msa->nseq; s++)
    for (i = 0; i < msa->alen; i++)
      if (strchr("-._~", msa->aseq[s][i]))
        ngap[i] += (use_wgts ? msa->wgt[s] : 1.0);

  for (i = 0; i < msa->alen; i++)
    {
      if (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)    mask[i] = FALSE;  // pathological; no sequences in MSA
      else                     mask[i] = (ngap[i] / ntot < gapfrac ? TRUE : FALSE);
    }
  free(ngap);
  return;

 ERROR:
  esl_fatal("allocation error"); //NOTREACHED
}

/* mask_by_symfrac()
 *
 * msa      : the multiple sequence alignment
 * symfrac  : keep columns with >= this fraction of residues
 * use_wgts : TRUE to calculate sequence-weighted fractions (MSA has weight annotations)
 * rfmask   : if non-NULL, implement do_keeprf: [0..alen-1] 0|1 flags marking consensus columns
 * mask     : RESULT, caller provided space : mask[0..alen-1] 1:0 flags marking cols to keep
 */
static void
mask_by_symfrac(const ESL_MSA *msa, double symfrac, int use_wgts, const int *rfmask, int *mask)
{
  double  ntot  = (use_wgts ? esl_vec_DSum(msa->wgt, msa->nseq) : msa->nseq);
  int64_t s,i;
  double *nsym = NULL;
  int     status;

  ESL_ALLOC(nsym, sizeof(double) * msa->alen);
  esl_vec_DSet(nsym, msa->alen, 0.);

  /* collect <nsym> for all columns in stride with aseq[s][i]*/
  for (s = 0; s < msa->nseq; s++)
    for (i = 0; i < msa->alen; i++)
      if (! strchr("-._~", msa->aseq[s][i]))
        nsym[i] += (use_wgts ? msa->wgt[s] : 1.0);

  for (i = 0; i < msa->alen; i++)
    {
      if (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)    mask[i] = FALSE;     // pathological; no sequences in MSA
      else                     mask[i] = (nsym[i] / ntot >= symfrac ? TRUE : FALSE);
    }
  free(nsym);
  return;

 ERROR:
  esl_fatal("allocation error"); //NOTREACHED
}

/* mask_by_maskfile()
 *
 * alen     : number of columns in MSA
 * maskfile : file to read user's mask from
 * rfmap    : optional : if non-NULL, use_conscoords: rfmap[0..C-1] gives coord for each cons column in MSA coords [0..alen-1] 
 * C        : if rfmap is non-NULL, number of consensus columns in MSA ; else 0
 * mask     : RESULT, caller provided space : mask[0..alen-1] 1:0 flags marking cols to keep
 */
static void
mask_by_maskfile(int64_t alen, const char *maskfile, const int64_t *rfmap, int64_t C, int *mask)
{
  int     *usermask;
  int64_t  mlen;
  int64_t  c;

  read_maskfile(maskfile, &usermask, &mlen);

  if (rfmap)
    {
      if (mlen != C) esl_fatal("expected consensus mask of len C=%" PRId64 " in mask file, but it's %" PRId64 " long", C, mlen);
      esl_vec_ISet(mask, alen, FALSE);    // insert columns are FALSE
      for (c = 0; c < C; c++)
        mask[rfmap[c]] = usermask[c];
    }
  else
    {
      if (mlen != alen) esl_fatal("expected mask of alen=%" PRId64 " in mask file, but it's %" PRId64 " long", alen, mlen);
      esl_vec_ICopy(usermask, alen, mask);
    }
  free(usermask);
}

/* mask_by_span()
 *
 * User provided <pos1>, <pos2> in 1..alen coords, or 1..C if
 * use_conscoords option is set (i.e. if <rfmap> is
 * provided). Internally, we work on 0..alen and 0..C-1 coords in our
 * arrays.
 *
 * Caller has already verified that 1 <= pos1 <= pos2 <= alen (or C).
 *
 * We don't need to know C, even though rfmap (if provided) is [0..C-1].
 *
 * alen  : number of columns in MSA
 * pos1  : start coord; either 0..alen-1 or (if rfmap non-NULL) 0..C-1
 * pos2  : end coord; ditto ""
 * rfmap : optional : if non-NULL, use_conscoords: rfmap[0..C-1] gives coord for each cons column in MSA coords [0..alen-1]
 * mask  : RESULT, caller provided space : mask[0..alen-1] 1:0 flags marking cols to keep
 */
static void
mask_by_span(int64_t alen, int64_t pos1, int64_t pos2, const int64_t *rfmap, int *mask)
{
  int64_t i;

  esl_vec_ISet(mask, alen, FALSE);

  if (rfmap) { pos1 = rfmap[pos1-1]; pos2 = rfmap[pos2-1]; }  // -1's here because of the shift from user 1.. coords to internal 0..
  else       { pos1 = pos1-1;        pos2 = pos2-1;        }

  for (i = pos1; i <= pos2; i++)
    mask[i] = TRUE;
}

/* mask_by_ppavg()
 *
 * msa       : the MSA; already checked that all seqs have individual PP annotation.
 * ppavg     : keep columns with >= this threshold ppavg
 * use_wgts  : TRUE to calculate sequence-weighted fractions (MSA has weight annotations)
 * rfmask    : optional : if non-NULL, [0..alen-1] 1|0 flags marking columns that are always kept (implements do_keeprf, for keeping consensus cols)
 * mask      : RESULT, caller provided space : [0..alen-1] 1:0 flags marking which cols to keep
 *
 * discretized postprob symbols [0-9*] are interpreted as the mean of their bin range.
 */
static void
mask_by_ppavg(const ESL_MSA *msa, double ppavg, int use_wgts, const int *rfmask, int *mask)
{
  double  *ppsum = NULL;
  double  *ntot  = NULL;
  int64_t  s,i;
  int      status;

  ESL_DASSERT1(( msa->pp ));

  ESL_ALLOC(ppsum, sizeof(double) * msa->alen);
  ESL_ALLOC(ntot,  sizeof(double) * msa->alen);
  esl_vec_DSet(ppsum, msa->alen, 0.0);
  esl_vec_DSet(ntot,  msa->alen, 0.0);

  for (s = 0; s < msa->nseq; s++)
    {
      ESL_DASSERT1(( msa->pp[s] ));
      for (i = 0; i < msa->alen; i++)
        {
          if (! strchr("-._~", msa->aseq[s][i]))
            {
              if ( ! strchr("0123456789*", msa->pp[s][i]))
                esl_fatal("invalid PP residue annotation char %c: sequence %s, apos %" PRId64, msa->pp[s][i], msa->sqname[s], i+1);

              ntot[i]  += (use_wgts ? msa->wgt[s] : 1.0);
              ppsum[i] += (use_wgts ? msa->wgt[s] * pp_to_mean(msa->pp[s][i]) : pp_to_mean(msa->pp[s][i]));
            }
          else if (! strchr("-._~", msa->pp[s][i]))
            esl_fatal("invalid PP gap annotation char %c: sequence %s, apos %" PRId64, msa->pp[s][i], msa->sqname[s], i+1);
        }
    }

  for (i = 0; i < msa->alen; i++)
    {
      if (rfmask && rfmask[i]) mask[i] = TRUE;    // do_keeprf option to keep consensus columns
      else if (ntot[i] == 0.0) mask[i] = FALSE;   // no residues in column; column is removed
      else                     mask[i] = (ppsum[i] / ntot[i] >= ppavg ? TRUE : FALSE);
    }
  free(ppsum);
  free(ntot);
  return;

 ERROR:
  esl_fatal("allocation error"); //NOTREACHED
}


/* mask_by_ppcons()
 *
 * msa        : the MSA. Already checked that it has PP_cons annotation line (and RF too)
 * ppcons     : keep columns with consensus postprob >= ppcons, 0 <= ppcons <= 1
 * do_keepins : TRUE to keep insert columns and only remove consensus columns. 
 * rfmask     : [0..alen-1] 1|0 flags for consensus columns
 * mask       : RESULT, caller provided space : [0..alen-1] 1|0 flags marking which cols to keep
 *
 * discretized postprob symbols [0-9*] are interpreted as the min of their bin range.
 */
static void
mask_by_ppcons(const ESL_MSA *msa, double ppcons, int do_keepins, const int *rfmask, int *mask)
{
  int64_t i;

  for (i = 0; i < msa->alen; i++)
    if (rfmask[i])  mask[i] = (pp_to_min(msa->pp_cons[i]) >= ppcons ? TRUE : FALSE);
    else            mask[i] = (do_keepins ? TRUE : FALSE);
}

/* mask_by_ppfrac()
 *
 * Keep columns that have at least a fraction >= ppfrac of
 * well-aligned residues, where well-aligned is defined as a postprob
 * of >= ppfracthresh.
 *
 * Discretized postprob symbols [0-9*] are interpreted as the minimum of
 * their bin range.
 *
 * msa          : the MSA; already checked that all seqs have individual PP annotation
 * ppfrac       : keep columns with >= this fraction of well-aligned residues
 * ppfracthresh : well-aligned residues are defined as postprob >= this threshold
 * use_wgts     : TRUE to calculate sequence-weighted fractions (MSA has weight annotations)
 * rfmask       : optional : if non-NULL, [0..alen-1] 1|0 flags marking columns that are always kept (implements do_keeprf, for keeping consensus cols)
 * mask         : RESULT, caller provided space: [0..alen-1] 1|0 flags marking which cols to keep
 */
static void
mask_by_ppfrac(const ESL_MSA *msa, double ppfrac, double ppfracthresh, int use_wgts, const int *rfmask, int *mask)
{
  double *nsat = NULL;
  double *ntot = NULL;
  int64_t s,i;
  int     status;

  ESL_DASSERT1(( msa->pp ));

  ESL_ALLOC(nsat, sizeof(double) * msa->alen);
  ESL_ALLOC(ntot, sizeof(double) * msa->alen);
  esl_vec_DSet(nsat, msa->alen, 0.0);
  esl_vec_DSet(ntot, msa->alen, 0.0);

  /* collect nsat[i] and ntot[i] for all columns: count residues that satisfy threshold, and total residue count */
  for (s = 0; s < msa->nseq; s++)
    {
      ESL_DASSERT1(( msa->pp[s] ));
      for (i = 0; i < msa->alen; i++)
        {
          if (! strchr("-._~", msa->aseq[s][i]))
            {
              if ( ! strchr("0123456789*", msa->pp[s][i]))
                esl_fatal("invalid PP residue annotation char %c: sequence %s, apos %" PRId64, msa->pp[s][i], msa->sqname[s], i+1);

              ntot[i] +=  (use_wgts ? msa->wgt[s] : 1.0);

              if ( pp_to_min(msa->pp[s][i]) >= ppfracthresh )
                nsat[i] += (use_wgts ? msa->wgt[s] : 1.0);
            }
          else if (! strchr("-._~", msa->pp[s][i]))
            esl_fatal("invalid PP gap annotation char %c: sequence %s, apos %" PRId64, msa->pp[s][i], msa->sqname[s], i+1);
        }
    }

  /* now check against threshold <ppfrac> */
  for (i = 0; i < msa->alen; i++)
    {
      if (rfmask && rfmask[i]) mask[i] = TRUE;   // do_keeprf option to keep consensus columns
      else if (ntot[i] == 0.0) mask[i] = FALSE;  // no residues in column; column is removed
      else                     mask[i] = (nsat[i] / ntot[i] >= ppfrac ? TRUE : FALSE);
    }

  free(nsat);
  free(ntot);
  return;

 ERROR:
  esl_fatal("allocation error"); //NOTREACHED

}

/***************************************************************** 
 * 4. selection options specialized for --small
 *****************************************************************/

/* small_mask_by_gapfrac()
 * small_mask_by_symfrac()
 *
 * Small-memory variants of mask_by_gapfrac()/mask_by_symfrac() that
 * use the histogram <abc_ct> produced by a first pass of
 * esl_msafile2_ReadInfoPfam(): abc_ct[i][K] is the gap count at
 * column i, and abc_ct[i][0..K-1] are the (possibly fractional, for
 * degenerate residues) per-symbol residue counts. 
 *
 * These are unweighted counts; -w is disallowed with --small.  The
 * legacy ESL_MSAFILE2 parser doesn't parse weights.
 *
 * Note one subtle difference from the default path: esl_abc_DCount() doesn't
 * count missing-data (~) or nonresidue (*) symbols at all, whereas the
 * default path's strchr("-._~",...) lumps ~ in with gaps. For ordinary Pfam
 * MSAs (no ~) the two paths agree exactly.
 *
 * abc_ct  : [0..alen-1][0..K] per-column residue/gap counts
 * alen    : number of columns
 * K       : alphabet size; abc_ct[i][K] is the gap count
 * gapfrac : (gapfrac variant) keep columns with gap fraction < this
 * symfrac : (symfrac variant) keep columns with residue fraction >= this
 * rfmask  : if non-NULL, implement do_keeprf: [0..alen-1] flags for always-kept consensus cols
 * mask    : RESULT, caller-provided : mask[0..alen-1] 1:0 flags marking cols to keep
 */
static void
small_mask_by_gapfrac(double **abc_ct, int64_t alen, int K, double gapfrac, const int *rfmask, int *mask)
{
  int64_t i;
  double  ngap, ntot;

  for (i = 0; i < alen; i++)
    {
      ngap = abc_ct[i][K];
      ntot = esl_vec_DSum(abc_ct[i], K+1);   // residues (0..K-1) + gaps (K)
      if      (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)         mask[i] = FALSE;  // pathological; no sequences
      else                          mask[i] = (ngap / ntot < gapfrac ? TRUE : FALSE);
    }
}

static void
small_mask_by_symfrac(double **abc_ct, int64_t alen, int K, double symfrac, const int *rfmask, int *mask)
{
  int64_t i;
  double  nsym, ntot;

  for (i = 0; i < alen; i++)
    {
      nsym = esl_vec_DSum(abc_ct[i], K);     // residues only, 0..K-1
      ntot = nsym + abc_ct[i][K];            // + gaps
      if      (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)         mask[i] = FALSE;  // pathological; no sequences
      else                          mask[i] = (nsym / ntot >= symfrac ? TRUE : FALSE);
    }
}

/* small_mask_by_ppavg()
 * small_mask_by_ppfrac()
 *
 * Small-memory variants of mask_by_ppavg()/mask_by_ppfrac() that use
 * the histogram <pp_ct> from esl_msafile2_ReadInfoPfam():
 * pp_ct[i][0..9] count symbols '0'-'9', pp_ct[i][10] counts '*', and
 * pp_ct[i][11] counts gaps.  Bins 0..10 are the residues with PP
 * annotation; their sum is the per-column residue count.
 *
 * The discretized bins are interpreted via the same pp_to_mean()/pp_to_min()
 * used by the default path: ppavg uses the bin mean, ppfrac the bin min.
 * There's no <use_wgts>: counts are unweighted, -w is disallowed with --small.
 *
 * pp_ct        : [0..alen-1][0..11] per-column PP symbol counts
 * alen         : number of columns
 * ppavg        : (ppavg variant)  keep columns with mean PP >= this
 * ppfrac       : (ppfrac variant) keep columns with >= this fraction of residues PP >= ppfracthresh
 * ppfracthresh : (ppfrac variant) threshold defining a confidently-aligned residue
 * rfmask       : if non-NULL, implement do_keeprf: always-kept consensus cols
 * mask         : RESULT, caller-provided : mask[0..alen-1] 1:0 flags
 */
static const char pp_bin_char[11] = {'0','1','2','3','4','5','6','7','8','9','*'};

static void
small_mask_by_ppavg(double **pp_ct, int64_t alen, double ppavg, const int *rfmask, int *mask)
{
  int64_t i;
  int     b;
  double  ppsum, ntot;

  for (i = 0; i < alen; i++)
    {
      ppsum = ntot = 0.0;
      for (b = 0; b < 11; b++)
        {
          ntot  += pp_ct[i][b];
          ppsum += pp_ct[i][b] * pp_to_mean(pp_bin_char[b]);
        }
      if      (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)         mask[i] = FALSE;  // no residues in column; removed
      else                          mask[i] = (ppsum / ntot >= ppavg ? TRUE : FALSE);
    }
}

static void
small_mask_by_ppfrac(double **pp_ct, int64_t alen, double ppfrac, double ppfracthresh, const int *rfmask, int *mask)
{
  int64_t i;
  int     b;
  double  nsat, ntot;

  for (i = 0; i < alen; i++)
    {
      nsat = ntot = 0.0;
      for (b = 0; b < 11; b++)
        {
          ntot += pp_ct[i][b];
          if (pp_to_min(pp_bin_char[b]) >= ppfracthresh) nsat += pp_ct[i][b];
        }
      if      (rfmask && rfmask[i]) mask[i] = TRUE;
      else if (ntot == 0.0)         mask[i] = FALSE;  // no residues in column; removed
      else                          mask[i] = (nsat / ntot >= ppfrac ? TRUE : FALSE);
    }
}



/*****************************************************************
 * 4. other internal functions
 *****************************************************************/

/* create_rfmask()
 *
 * rfmask[0..alen-1] is an array of 1|0 flags, 1's marking consensus
 * columns, 0's marking insert columns, based on nongap chars vs. gap
 * chars in the input <rf> annotation line.
 *
 * rf         : reference annotation string for the MSA
 * alen       : # of columns in MSA
 * ret_rfmask : RETURN, allocated here : rfmask[0..alen-1]
 * ret_C      : RETURN : # of annotated consensus columns
 */
static void
create_rfmask(const char *rf, int64_t alen, int **ret_rfmask, int64_t *ret_C)
{
  int    *rfmask = NULL;
  int64_t C      = 0;
  int64_t i;
  int     status;

  ESL_ALLOC(rfmask, sizeof(int) * alen);

  for (i = 0; i < alen; i++)
    if (strchr("-._~", rf[i]))  rfmask[i] = FALSE;
    else                      { rfmask[i] = TRUE;  C++; }
  *ret_rfmask = rfmask;
  *ret_C      = C;
  return;

 ERROR:
  esl_fatal("allocation failed"); //NOTREACHED
}

/* create_rfmap()
 *
 * rfmap[0..C-1] is an array of coordinates 0..alen-1, mapping
 * consensus-only column coords to full MSA coords; i.e.  rfmap[c] is
 * the column index (0..alen-1) of consensus position c.
 *
 * It's possible to have C=0, in which case we still allocate a tiny
 * rfmap of size 1, just so it's non-NULL. Some code checks `if
 * (rfmap)` as a switch for consensus coords.
 *
 * rfmask    : [0..alen-1] array of 1|0 flags marking consensus|insert cols in MSA
 * alen      : total # of cols in MSA
 * C         : # of consensus cols; known from # of 1's in <rfmask>
 * ret_rfmap : RETURN : rfmap[0..C-1] array of (0..alen-1) coords for each consensus col
 */
static void
create_rfmap(const int *rfmask, int64_t alen, int64_t C, int64_t **ret_rfmap)
{
  int64_t *rfmap = NULL;
  int64_t  c     = 0;
  int64_t  i;
  int      status;

  if (C) ESL_ALLOC(rfmap, sizeof(int64_t) * C);
  else   ESL_ALLOC(rfmap, sizeof(int64_t) * 1);  // the edge case

  for (i = 0; i < alen; i++)
    if (rfmask[i]) rfmap[c++] = i;
  ESL_DASSERT1(( c == C ));

  *ret_rfmap = rfmap;
  return;

 ERROR:
  esl_fatal("allocation failed"); //NOTREACHED
}

/* add_mask()
 *
 * <overall_mask> = <new_mask> & <overall_mask>
 */
static void
add_mask(const int *new_mask, int64_t alen, int *overall_mask)
{
  int i;
  for (i = 0; i < alen; i++)
    overall_mask[i] = new_mask[i] & overall_mask[i];
}


/* read_maskfile()
 *
 * Read a mask consisting of 1|0 characters from <filename>, ignoring
 * whitespace and comments ('#' and anything that follows a '#' on a
 * line). Return the corresponding mask of 1|0 flags in <*ret_mask>
 * (allocated here), and its length in <*ret_mlen>. Caller can then
 * check that it's the expected length (alen or C).
 *
 * filename : file to read mask from
 * ret_mask : RETURN, allocated here : [0..mlen-1] array of 1|0 flags
 * ret_mlen : RETURN: length of mask
 */
static void
read_maskfile(const char *filename, int **ret_mask, int64_t *ret_mlen)
{
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  int            *mask = NULL;
  int64_t         mlen = 0;
  int             i;
  int             status;

  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK)
    esl_fatal("Failed to open mask file %s", filename);
  esl_fileparser_SetCommentChar(efp, '#');

  while ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK)
    {
      if (mlen == 0) ESL_ALLOC  (mask, sizeof(int) * toklen);
      else           ESL_REALLOC(mask, sizeof(int) * (mlen+toklen));

      for (i = 0; i < toklen; i++)
        {
          if      (tok[i] == '0') mask[mlen++] = FALSE;
          else if (tok[i] == '1') mask[mlen++] = TRUE;
          else esl_fatal("Input mask must consist only of 0|1's. Character %" PRId64 " is a %c.", mlen+1, tok[i]);
        }
    }

  *ret_mask = mask;
  *ret_mlen = mlen;
  esl_fileparser_Close(efp);
  return;

 ERROR:
  esl_fatal("Allocation failure"); //NOTREACHED
}

/* write_maskfile()
 *
 * Write <mask> of length <alen> to <filename> as a string of 1's and 0's.
 */
static void
write_maskfile(const char *filename, const int *mask, int64_t alen)
{
  FILE    *ofp = NULL;
  int64_t  i;

  if ((ofp = fopen(filename, "w")) == NULL)
    esl_fatal("failed to open mask file %s for writing", filename);

  for (i = 0; i < alen; i++)
    fputc( (mask[i] ? '1' : '0'), ofp);
  fputc('\n', ofp);
  fclose(ofp);
}

/* write_rfmaskfile()
 *
 * Write a mask string (a string of 1 and 0 characters) for consensus columns
 * only to <filename>, using mask[0..alen-1] and an rfmap[0..C-1].
 *
 * filename : file to write to
 * mask     : full mask for the MSA, [0..alen-1], 1:0 flags
 * rfmap    : [0..C-1] map of consensus columns to coords (0,alen-1) in full MSA
 * C        : # of consensus cols; length of <rfmap>
 */
static void
write_rfmaskfile(const char *filename, const int *mask, const int64_t *rfmap, int64_t C)
{
  FILE    *ofp = NULL;
  int64_t  c;

  if ((ofp = fopen(filename, "w")) == NULL)
    esl_fatal("failed to open consensus mask file %s for writing", filename);

  for (c = 0; c < C; c++)
    fputc( (mask[rfmap[c]] ? '1' : '0'), ofp);
  fputc('\n', ofp);
  fclose(ofp);
}

/* write_verbose_header()
 *
 * Write the table header for the verbose output.
 */
static void
write_verbose_header(FILE *ofp)
{
  esl_fprintf(ofp, "%-19s  %7s  %7s  %16s  %16s\n",          "",                    "",        "",        "all columns",        "consensus cols");
  esl_fprintf(ofp, "%-19s  %7s  %7s  %16s  %16s\n",          "",                    "",        "",        "----------------",   "----------------");
  esl_fprintf(ofp, "%-19s  %7s  %7s  %7s  %7s  %7s  %7s\n",  "",                    "",        "consens", "num",     "num",     "num",     "num");
  esl_fprintf(ofp, "%-19s  %7s  %7s  %7s  %7s  %7s  %7s\n",  "selection mode",      "aln len", "len",     "kept",    "removed", "kept",    "removed");
  esl_fprintf(ofp, "%-19s  %7s  %7s  %7s  %7s  %7s  %7s\n",  "-------------------", "-------", "-------", "-------", "-------", "-------", "-------");
}

/* write_verbose_line()
 *
 * ofp              : output stream to write to (stdout)
 * selection_option : name of selection option, or NULL for final overall_mask
 * mask             : mask [0..alen-1]
 * alen             : length of mask in columns
 * rfmask           : [0..alen-1] 1|0 flags for consensus cols; or NULL if no consensus col annotation
 * C                : number of consensus columns marked in <rfmask>, or 0
 */
static void
write_verbose_line(FILE *ofp, const char *selection_option, const int *mask, int64_t alen, const int *rfmask, int64_t C)
{
  int64_t nkept  = esl_vec_ISum(mask, alen);
  int64_t rfkept = 0;
  int64_t i;

  if (!selection_option)  // for overall_mask, insert a summing up line
    esl_fprintf(ofp, "%-19s  %7s  %7s  %7s  %7s  %7s  %7s\n",  "-------------------", "-------", "-------", "-------", "-------", "-------", "-------");

  if (rfmask)
    {
      for (i = 0; i < alen; i++)
        if (rfmask[i] && mask[i]) rfkept++;

      esl_fprintf(ofp, "%-19s  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "\n",
                  selection_option ? selection_option : "overall", alen, C, nkept, alen-nkept, rfkept, C-rfkept);
    }
  else
    esl_fprintf(ofp, "%-19s  %7" PRId64 "  %7s  %7" PRId64 "  %7" PRId64 "  %7s  %7s\n",
                selection_option ? selection_option : "overall", alen, "-", nkept, alen-nkept, "-", "-");
}

/* pp_to_mean()
 *
 * Convert a discretized PP annotation symbol "0-9*" to the mean of
 * its bin.
 */
static double
pp_to_mean(char pp)
{
  switch (pp) {
  case '0': return 0.025;
  case '1': return 0.1;
  case '2': return 0.2;
  case '3': return 0.3;
  case '4': return 0.4;
  case '5': return 0.5;
  case '6': return 0.6;
  case '7': return 0.7;
  case '8': return 0.8;
  case '9': return 0.9;
  case '*': return 0.975;
  default: esl_fatal("invalid PP annotation");
  }
}

/* pp_to_min()
 *
 * Convert a discretized PP annotation symbol "0-9*" to the minimum of
 * its bin.
 */
static double
pp_to_min(char pp)
{
  switch (pp) {
  case '0': return 0.0;
  case '1': return 0.05;
  case '2': return 0.15;
  case '3': return 0.25;
  case '4': return 0.35;
  case '5': return 0.45;
  case '6': return 0.55;
  case '7': return 0.65;
  case '8': return 0.75;
  case '9': return 0.85;
  case '*': return 0.95;
  default: esl_fatal("invalid PP annotation");
  }
}
