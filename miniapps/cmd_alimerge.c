/* `easel alimerge`: merge MSAs into one, via reference (RF) annotation
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile2.h"
#include "esl_subcmd.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#define ALPHOPTS "--amino,--dna,--rna"

static ESL_OPTIONS cmd_options[] = {
  /* name         type          default  env   range togs reqs  incomp     help                                                     docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "help; show brief info on version and usage",                  1 },
  { "-o",         eslARG_OUTFILE,  NULL, NULL, NULL, NULL,NULL, NULL,      "output merged MSA to file <f>, not stdout",                   1 },
  { "-v",         eslARG_NONE,    FALSE, NULL, NULL, NULL,"-o", NULL,      "print info about merge to stdout; requires -o",               1 },
  { "--rfonly",   eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "remove all columns that are gaps in GC RF annotation",        1 },
  { "--small",    eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "use minimal RAM (make RAM usage independent of MSA sizes)",   1 },
  { "--outformat",eslARG_STRING,   NULL, NULL, NULL, NULL,NULL,"--small",  "write merged output MSA in format <s>",                       1 },
  { "--amino",    eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, ALPHOPTS,  "assert that merged MSAs are protein sequences",               1 },
  { "--dna",      eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, ALPHOPTS,  "assert that merged MSAs are DNA sequences",                   1 },
  { "--rna",      eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, ALPHOPTS,  "assert that merged MSAs are RNA sequences",                   1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


static void    show_verbose_header (char **msafile, int nfile, int *ret_namewidth);
static void    default_protocol    (char **msafile, int nfile, int infmt, int namewidth, int outfmt, FILE *ofp, int be_verbose);
static int     update_maxgap_maxmis(const ESL_MSA *msa, int *maxgap, int *maxmis, char *errbuf);
static void    transfer_annotation (ESL_MSA *mmsa,       ESL_MSA **msa, int nmsa, const int *maxgap, const int *maxmis, int C, int be_verbose);
static void    merge_msa           (ESL_MSA *mmsa, const ESL_MSA  *msa,           const int *maxgap, const int *maxmis, int C);
static void    inflate_string      (const char *s, const char *rf, int L, const int *maxgap, const int *maxmis, int C, int do_local_rna_tildes, int do_frag_tildes, int new_alen, char **ret_as);

static int     get_msa_consensus_len(const ESL_MSA *msa);
static int *   create_msa_consensus_mask(const ESL_MSA *msa);
static int     first_gc_ok(const char *gc,  const char *rf);
static int     other_gc_ok(const char *gc1, const char *rf1, const char *gc2, const char *rf2);
static int     other_rf_ok(const char *rf1, const char *rf2);

int
esl_cmd_alimerge(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS   *go         = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, /*custom opthelp_f=*/NULL);
  int            nfile      = esl_opt_ArgNumber(go);
  char         **msafile    = NULL;                               // array of MSA filenames; msafile[0..nfile-1]
  int            be_verbose = esl_opt_GetBoolean(go, "-v");       // with -v, we output a summary table to stdout, and merged MSA to -o <f> file
  int            namewidth  = 0;                                  // stays 0 if be_verbose==FALSE. `if (namewidth)` and `if (be_verbose)` do the same
  int            do_small   = esl_opt_GetBoolean(go, "--small");
  int            do_rfonly  = esl_opt_GetBoolean(go, "--rfonly");
  char          *outfile    = esl_opt_GetString (go, "-o");
  FILE          *ofp        = stdout;
  int            infmt      = (do_small? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM);
  int            outfmt     = (do_small? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM);
  ESL_STOPWATCH *w          = NULL;  
  int            i;
  int            status;
  
  if (nfile == 0)
    esl_fatal("Incorrect number of cmdline arguments.\nProvide at least one <msafile> containing MSAs to merge.");

  if ( esl_opt_IsOn(go, "--outformat")) {  // --small and --outformat are incompatible; outfmt == eslMSAFILE_PFAM if --small. getopts already checked this.
    if ((outfmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--outformat"))) == eslMSAFILE_UNKNOWN)
      esl_fatal("%s is not a valid output MSA file format for --outformat\n", esl_opt_GetString(go, "--outformat"));
  }

  if (outfile) {
    if ((ofp = fopen(outfile, "w")) == NULL) 
      esl_fatal("Failed to open -o output MSA file %s\n", outfile);
  } 

  ESL_ALLOC(msafile, sizeof(char *) * nfile);   // each msafile[i] is just a ptr into <go> cmdline parser.  
  for (i = 1; i <= nfile; i++)                  // msafiles[i] names are not allocated.
    msafile[i-1] = esl_opt_GetArg(go, i);       // watch off-by-one: cmdline args are 1..nfile, but msafile[0..nfile-1]

  if (be_verbose) {
    show_verbose_header(msafile, nfile, &namewidth);
    w = esl_stopwatch_Create();
    esl_stopwatch_Start(w);
  }

  if (! do_small)
    default_protocol(msafile, nfile, infmt, namewidth, outfmt, ofp, be_verbose);
  else
    esl_fatal("--small isn't implemented yet");

  free(msafile);
  if (ofp != stdout) fclose(ofp);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  return status;
}

static void
show_verbose_header(char **msafile, int nfile, int *ret_namewidth)
{
  int   namewidth  = 9;
  char *namedashes = NULL;
  char *basename   = NULL;
  int   fi, i;
  int   status;

  /* Determine longest filename in the <msafile> list */
  for (fi = 0; fi < nfile; fi++)
    {
      esl_FileTail(msafile[fi], /*nosuffix=*/FALSE, &basename);
      namewidth = ESL_MAX(namewidth, strlen(basename));
      free(basename);
    }

  /* Variable-length table underscoring for the filename field */
  ESL_ALLOC(namedashes, sizeof(char) * (namewidth+1));
  for (i = 0; i < namewidth; i++) namedashes[i] = '-';
  namedashes[namewidth] = '\0';

  /* Format the header for the verbose output table */
  esl_printf("# Reading %d alignment files\n", nfile);
  esl_printf("#\n");
  esl_printf("# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "",        namewidth, "",          "",        "",          "",            "cumulative",    "ncols");
  esl_printf("# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "file #",  namewidth, "file name", "ali #",   "#seq/ali",  "ncols/ali",   "# seq total",   "required");
  esl_printf("# %7s  %*s  %7s  %9s  %9s  %13s  %8s\n",  "-------", namewidth, namedashes,  "-------", "---------", "---------",   "-------------", "--------");

  free(namedashes);
  *ret_namewidth = namewidth;
  return;

 ERROR: // NOTREACHED
  esl_fatal("allocation error");
}


/* default_protocol()
 *
 * The standard flow: read all input MSAs into memory, then merge them
 * one at a time into the output MSA. (As opposed to the small memory
 * protocol with --small.)
 */
static void
default_protocol(char **msafile, int nfile, int infmt, int namewidth, int outfmt, FILE *ofp, int be_verbose)
{
  ESL_MSAFILE *afp        = NULL;
  ESL_MSA    **msa        = NULL;    // array of all input MSAs in all files
  ESL_MSA     *mmsa       = NULL;    // output merged MSA
  int          msa_nalloc = 2;       // <msa> is dynamically resized, starting from this
  int          C          = 0;       // consensus length, according to all #=GC RF annotations on all MSAs
  int         *maxgap     = NULL;    // maxgap[0..c..C] = max # of gaps seen before consensus column <c>. maxgap[0] = leading, maxgap[C] = trailing
  int         *maxmis     = NULL;    // maxmis[0..c..C] = max # of missing data '~' seen before consensus column <c>
  int          nseq_tot   = 0;       // total # of seqs over all MSAs
  int          nmsa       = 0;       // total # of MSAs in all files
  int          ai;                   // counter over all MSAs in all files
  int          ai2;                  // counter over MSAs in current file (only needed for error reporting)
  int          fi;                   // counter over MSA files
  int          i;                    // generic loop counter
  char         errbuf[eslERRBUFSIZE];
  int          status;

  /* Initial allocation for MSAs */
  ESL_ALLOC(msa, sizeof(ESL_MSA *) * msa_nalloc);
  for (ai = 0; ai < msa_nalloc; ai++) msa[ai] = NULL;

  /* Pass #1 over the MSA files */
  for (ai = 0, fi = 0; fi < nfile; fi++)
    {
      status = esl_msafile_Open(/*abc=*/NULL, msafile[fi], /*env=*/NULL, infmt, NULL, &afp);
      if (status != eslOK) esl_msafile_OpenFailure(afp, status);

      ai2 = 0;  // counter over MSAs in current file; only needed for error reporting
      while ((status = esl_msafile_Read(afp, &(msa[ai]))) == eslOK)
        {
          if (msa[ai]->rf == NULL)
            esl_fatal("All MSAs to merge must have #=GC RF consensus annotation.\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                      ai2, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai);
          
          if (ai == 0) {
            C = get_msa_consensus_len(msa[0]);
            ESL_ALLOC(maxgap, sizeof(int) * (C+1));               
            ESL_ALLOC(maxmis, sizeof(int) * (C+1));
            esl_vec_ISet(maxgap, C+1, 0);
            esl_vec_ISet(maxmis, C+1, 0);
          } else {
            if (! other_rf_ok(msa[0]->rf, msa[ai]->rf) || get_msa_consensus_len(msa[ai]) != C)  // checking C is redundant, since it already passed other_rf_ok, but doesn't hurt
              esl_fatal("All MSAs must have compatible consensus annotation (#=GC RF) lines.\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                        ai2, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai);
          }

          if ((status = update_maxgap_maxmis(msa[ai], maxgap, maxmis, errbuf)) != eslOK)
            esl_fatal("Problem with consensus annotation line (#=GC RF)\n  in MSA #%d (%s) in file %s (#%d MSA overall)\n  %s", errbuf);

          nseq_tot += msa[ai]->nseq;

          /* Output to summary table, if we're printing one */
          if (be_verbose)
            {
              char *basename;
              esl_FileTail(msafile[fi], /*nosuffix=*/FALSE, &basename);
              esl_printf("  %7d  %-*s  %7d  %9d  %9" PRId64 "  %13d  %8d\n",
                         (fi+1),  namewidth, basename,
                         (ai+1),  msa[ai]->nseq,  msa[ai]->alen, nseq_tot,
                         C + esl_vec_ISum(maxgap, C+1) + esl_vec_ISum(maxmis, C+1));
              free(basename);
            }

          /* Reallocation for more MSAs as needed */
          ai++;
          ai2++;            
          if (ai == msa_nalloc)
            {
              msa_nalloc *= 2;
              ESL_REALLOC(msa,  sizeof(ESL_MSA *) * msa_nalloc);
              for (i = ai; i < msa_nalloc; i++) msa[i]  = NULL;
            }
        }
      if (status != eslEOF) esl_msafile_ReadFailure(afp, status);
      esl_msafile_Close(afp);
    }
  nmsa = ai;

  /* First read pass finished.
   * Now we have all MSAs stored in msa[], and maxgap[],maxmis[] to tell us how to expand each alignment for the merge.
   * <nali>     is total number of input MSAs
   * <nseq_tot> is total number of seqs in merged MSA
   */
  mmsa       = esl_msa_Create(nseq_tot, -1);
  mmsa->alen =  C + esl_vec_ISum(maxgap, C+1) + esl_vec_ISum(maxmis, C+1);
  /* by setting <alen> *after* creating the <mmsa>, we're responsible for all <alen+1> aligned string allocations */


  /* Per-file and per-column annotations need to be checked for
   * consistency across all individual MSAs.
   */
  transfer_annotation(mmsa, msa, nmsa, maxgap, maxmis, C, be_verbose);
  
  /* Merge input MSAs into output MSA one at a time, free'ing as we go. */
  for (ai = 0; ai < nmsa; ai++)
    {
      merge_msa(mmsa, msa[ai], maxgap, maxmis, C);
      esl_msa_Destroy(msa[ai]);
    }

  /* Output the merged MSA */
  esl_msafile_Write(ofp, mmsa, outfmt);
  
  esl_msa_Destroy(mmsa);
  free(maxgap);
  free(maxmis);
  free(msa);
  return;

 ERROR: // UNREACHED
  return;
}


/* update_maxgap_maxmis()
 * 
 * The maxgap[] and maxmis[] arrays keep track of maximum number of
 * gap columns and missing data columns annotated in #=GC RF reference
 * annotation over all MSAs.
 *
 * For each <msa> in turn, with <msa->rf> annotation, given the
 * current <maxgap> and <maxmis>, update them as needed.
 *
 * There are <C> consensus columns, indexed 0..c..C-1. maxgap[0..c..C]/maxmis[0..c..C]
 * are the max number of gaps/~ *before* column c, where maxgap[C]/maxmis[C]
 * are for the final columns trailing the last consensus column. 
 *
 * Missing data ~ columns must precede gap columns between any two consensus columns
 * in the RF annotation (including leading and trailing nonconsensus columns). 
 * If this is not true, an <eslEINVAL> error is returned; if caller provides
 * an <errbuf> allocated for at least <eslERRBUFSIZE> chars, an informative error
 * message is put there if this failure happens.
 */
static int
update_maxgap_maxmis(const ESL_MSA *msa, int *maxgap, int *maxmis, char *errbuf)
{
  int apos, cpos;
  int ngap = 0;
  int nmis = 0;

  for (apos = 0, cpos = 0; apos <= msa->alen; apos++)  // on last iteration, apos = alen and rf[alen] is on its '\0' terminator, and we'll set maxgap and maxmis for trailing flank.
    {
      if (strchr("-_.", msa->rf[apos]))
        ngap++;
      else if (msa->rf[apos] == '~') {
        nmis++;
        if (ngap) ESL_FAIL(eslEINVAL, errbuf, "gaps precede missing data ~ in RF annotation at col=%d", apos+1);
      } else {
        maxgap[cpos] = ESL_MAX(maxgap[cpos], ngap);
        maxmis[cpos] = ESL_MAX(maxmis[cpos], nmis);
        cpos++;
        ngap = nmis = 0;
      }
    } 
  return eslOK;
}
  

static int
get_msa_consensus_len(const ESL_MSA *msa)
{
  int C = 0;
  int apos;

  for (apos = 0; apos < msa->alen; apos++)
    if (! strchr("-_.~", msa->rf[apos]))
      C++;
  return C;
}

static int *
create_msa_consensus_mask(const ESL_MSA *msa) 
{
  int *useme = malloc(sizeof(int) * msa->alen);
  int  apos;

  if (useme == NULL) esl_fatal("malloc failed");

  for (apos = 0; apos < msa->alen; apos++)
    useme[apos] = (strchr("-_.~", msa->rf[apos]) ? FALSE : TRUE);
  return useme;
}
  

/* transfer_annotation()
 *
 * Only per-file and per-column annotations that are identical across
 * all MSAs is transferred to the output merged MSA.
 */
static void
transfer_annotation(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, const int *maxgap, const int *maxmis, int C, int be_verbose)
{
  char *tmpstr = NULL;
  int   ai;
  int   i, i2;
  int   status;

  /* name annotation */
  if (msa[0]->name) {
    for (ai = 1; ai < nmsa; ai++)
      if (esl_strcmp(msa[0]->name, msa[ai]->name) != 0) break;
    if (ai == nmsa) esl_msa_SetName(mmsa, msa[0]->name, -1);
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# msa name:    absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# msa name:    not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# msa name:    identical in all MSAs, transferred to merged MSA\n");
  }

  /* description */
  if (msa[0]->desc) {
    for (ai = 1; ai < nmsa; ai++)
      if (esl_strcmp(msa[0]->desc, msa[ai]->desc) != 0) break;
    if (ai == nmsa) esl_msa_SetDesc(mmsa, msa[0]->desc, -1);
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# description: absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# description: not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# description: identical in all MSAs, transferred to merged MSA\n");
  }

  /* accession */
  if (msa[0]->acc) {
    for (ai = 1; ai < nmsa; ai++)
      if (esl_strcmp(msa[0]->acc, msa[ai]->acc) != 0) break;
    if (ai == nmsa) esl_msa_SetAccession(mmsa, msa[0]->acc, -1);
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# accession:   absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# accession:   not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# accession:   identical in all MSAs, transferred to merged MSA\n");
  }

  /* author */
  if (msa[0]->au) {
    for (ai = 1; ai < nmsa; ai++)
      if (esl_strcmp(msa[0]->au, msa[ai]->au) != 0) break;
    if (ai == nmsa) esl_msa_SetAuthor(mmsa, msa[0]->au, -1);
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# author:      absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# author:      not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# author:      identical in all MSAs, transferred to merged MSA\n");
  }

  /* other unparsed per-file #=GF annotations
   * Must be identically in order and content in all msa's to merge.
   *
   * Must do it this way, as opposed to checking individual tags,
   * because Xfam doesn't obey the format, and has multiple GF lines
   * with the same tag.
   */
  if (msa[0]->ngf)
    {
      for (ai = 1; ai < nmsa; ai++) {
        if (msa[ai]->ngf != msa[0]->ngf) break;
        for (i = 0; i < msa[0]->ngf; i++) {
          if (esl_strcmp(msa[0]->gf_tag[i], msa[ai]->gf_tag[i]) != 0) break;
          if (esl_strcmp(msa[0]->gf[i],     msa[ai]->gf[i])     != 0) break;
        }
        if (i < msa[0]->ngf) break;
      }
      if (ai == nmsa)
        for (i = 0; i < msa[0]->ngf; i++)
          esl_msa_AddGF(mmsa, msa[0]->gf_tag[i], -1, msa[0]->gf[i], -1);

      if (be_verbose) {
        if (ai < nmsa) esl_printf("# other GF:    not identical in all MSAs, not transferred to merged MSA\n");
        else           esl_printf("# other GF:    identical in all MSAs, transferred to merged MSA\n");
      }
    }

  /* unparsed comments: must be identically ordered and identical in all msa's to merge */
  if (msa[0]->ncomment)
    {
      for (ai = 1; ai < nmsa; ai++) {
        if (msa[ai]->ncomment != msa[0]->ncomment) break;
        for (i = 0; i < msa[0]->ncomment; i++)
          if (esl_strcmp(msa[0]->comment[i], msa[ai]->comment[i]) != 0) break;
        if (i < msa[0]->ncomment) break;
      }
      if (ai == nmsa)
        for (i = 0; i < msa[0]->ncomment; i++)
          esl_msa_AddComment(mmsa, msa[0]->comment[i], -1);

      if (be_verbose) {
        if (ai < nmsa) esl_printf("# comments:      not identical in all MSAs, not transferred to merged MSA\n");
        else           esl_printf("# comments:      identical in all MSAs, transferred to merged MSA\n");
      }
    }

  /* SS_cons
   */
  if (msa[0]->ss_cons && first_gc_ok(msa[0]->ss_cons, msa[0]->rf))
    {
      for (ai = 1; ai < nmsa; ai++)
        if (! other_gc_ok(msa[0]->ss_cons,  msa[0]->rf, msa[ai]->ss_cons, msa[ai]->rf)) break;
      if (ai == nmsa) {
        ESL_ALLOC(mmsa->ss_cons, sizeof(char) * (mmsa->alen+1));  
        inflate_string(msa[0]->ss_cons, msa[0]->rf, msa[0]->alen, maxgap, maxmis, C, /*do_local_rna_tildes=*/TRUE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->ss_cons));
      }
    } else ai=0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# SS_cons:     absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# SS_cons:     not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# SS_cons:     identical in all MSAs, transferred to merged MSA\n");
  }

  /* SA_cons
   */
  if (msa[0]->sa_cons && first_gc_ok(msa[0]->sa_cons, msa[0]->rf))
    {
      for (ai = 1; ai < nmsa; ai++)
        if (! other_gc_ok(msa[0]->sa_cons,  msa[0]->rf,  msa[ai]->sa_cons, msa[ai]->rf)) break;
      if (ai == nmsa) {
        ESL_ALLOC(mmsa->sa_cons, sizeof(char) * (mmsa->alen+1)); 
        inflate_string(msa[0]->sa_cons, msa[0]->rf, msa[0]->alen, maxgap, maxmis, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->sa_cons));
      }
    } else ai=0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# SA_cons:     absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# SA_cons:     not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# SA_cons:     identical in all MSAs, transferred to merged MSA\n");
  }

  /* PP_cons
   */
  if (msa[0]->pp_cons && first_gc_ok(msa[0]->pp_cons, msa[0]->rf))
    {
      for (ai = 1; ai < nmsa; ai++)
        if (! other_gc_ok(msa[0]->pp_cons,  msa[0]->rf,  msa[ai]->pp_cons, msa[ai]->rf)) break;
      if (ai == nmsa) {
        ESL_ALLOC(mmsa->pp_cons, sizeof(char) * (mmsa->alen+1)); 
        inflate_string(msa[0]->pp_cons, msa[0]->rf, msa[0]->alen, maxgap, maxmis, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->pp_cons));
      }
    } else ai=0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# PP_cons:     absent (at least from first MSA); not included in merge\n");
    else if (ai < nmsa) esl_printf("# PP_cons:     not identical in all MSAs, not transferred to merged MSA\n");
    else                esl_printf("# PP_cons:     identical in all MSAs, transferred to merged MSA\n");
  }

  /* Other unparsed GC annotations */
  if (msa[0]->ngc > 0)
    for (i = 0; i < msa[0]->ngc; i++) 
      {
        if ( first_gc_ok(msa[0]->gc[i], msa[0]->rf )) {
          for (ai = 1; ai < nmsa; ai++) {
            for (i2 = 0; i2 < msa[ai]->ngc; i2++)
              if (esl_strcmp(msa[0]->gc_tag[i], msa[ai]->gc_tag[i2]) == 0) break;
            if (i2 == msa[ai]->ngf || ! other_gc_ok(msa[0]->gc[i],  msa[0]->rf, msa[ai]->gc[i], msa[ai]->rf))
              break;
          }
        } else ai=0;

        if (ai == nmsa) {
          ESL_ALLOC(tmpstr, sizeof(char) * (mmsa->alen+1));
          inflate_string(msa[0]->gc[i], msa[0]->rf, msa[0]->alen, maxgap, maxmis, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &tmpstr);
          esl_msa_AppendGC(mmsa, msa[0]->gc_tag[i], tmpstr);   // would be nice to have a esl_msa_SetGC() here to save the tmpstr allocation
          free(tmpstr);
          tmpstr = NULL;
        }
        if (be_verbose) {
          if (ai < nmsa) esl_printf("# GC %s:     not identical in all MSAs, not transferred to merged MSA\n", msa[0]->gc_tag[i]);
          else           esl_printf("# GC %s:     identical in all MSAs, transferred to merged MSA\n", msa[0]->gc_tag[i]);
        }
      }

  /* Finally, RF itself.
   * Caller already validated that all RF lines exist and are compatible.
   */
  inflate_string(msa[0]->rf, msa[0]->rf, msa[0]->alen, maxgap, maxmis, C, /*do_local_rna_tildes=*/TRUE, FALSE, mmsa->alen, &(mmsa->rf));
  if (be_verbose)
    esl_printf("# Identical RF annotation from all alignments transferred to merged alignment.\n");
  return;

 ERROR: //NOTREACHED
  return;
}

/* first_gc_ok()
 * other_gc_ok()
 * other_rf_ok()
 * 
 * Helper functions used by transfer_annotation to check whether a
 * #=GC annotation can be transferred to the merged MSA.
 *
 * first_gc_ok() verifies that GC annotation has no nongap
 * character in a nonconsensus column. Only consensus columns can have
 * annotations.
 *
 * other_gc_ok() then verifies that for all other GC annotations, plus
 * verifying that the GC annotation on consensus positions is
 * identical between the two GC's.
 *
 * other_rf_ok() is specifically for RF annotation. It only checks
 * that two RF annotations are identical at consensus positions.
 */
static int
first_gc_ok(const char *gc, const char *rf)
{
  int i;

  for (i = 0; rf[i] != '\0'; i++)
    if (strchr("-_.~", rf[i]) && ! strchr("-_.~", gc[i])) return FALSE;
  return TRUE;
}

static int
other_gc_ok(const char *gc1, const char *rf1, const char *gc2, const char *rf2)
{
  int i = 0;
  int j = 0;

  while (rf2[j] != '\0')  // rf1[i] is also '\0' when we end
    {
      while (strchr("-_.~", rf2[j])) { if (! strchr("-_.~", gc2[j])) return FALSE; j++; }  // advance j to next consensus position, while checking for nongap GC annotation
      while (strchr("-_.~", rf1[i]))   i++;                                                // advance i to next consensus position 
      if (gc1[i] != gc2[j]) return FALSE;         // check for gc identity at consensus positions. Includes a final comparison of \0:\0 at ends of strings.
      i++;
      j++;
    }
  return TRUE;
}

static int
other_rf_ok(const char *rf1, const char *rf2)
{
  int i = 0;
  int j = 0;

  while (rf2[j] != '\0')  // rf1[i] is also '\0' when we end
    {
      while (strchr("-_.~", rf2[j])) j++;  
      while (strchr("-_.~", rf1[i])) i++;  
      if (rf1[i] != rf2[j]) return FALSE;   
      i++;
      j++;
    }
  return TRUE;
}




/* merge_msa()
 *
 * Merge one input MSA into the new merged output MSA, given precalculated
 * arrays <maxgap> and <maxmis> that tell us how many padding gaps we need
 * to add at each insert region.
 *
 * Sequence-specific annotations are transferred here (#=GS and #=GR). These
 * are always transferred; since they're seq specific, they don't conflict
 * between input MSAs. Per-file and per-column annotations (#=GF and #=GC)
 * that might conflict are handled by transfer_annotation(). 
 *
 * <mmsa>    : new, merged MSA. Allocated by caller. 
 *             <mmsa->sqalloc> is already the total # seqs in all individual MSAs.
 *             <mmsa->alen>    is already the new length of the merged MSA.
 *             <mmsa->nseq>    is the running current # of seqs; we update it here.
 * <msa>     : MSA to merge into <mmsa>
 * <maxgap>  : array of lengths of each insert region in the merged MSA, [0..C].  
 * <maxmis>  : array of lengths of each missing data ~ region in the merged MSA, [0..C]
 * <C>       : number of consensus positions marked in <rf>; dictates the C+1 length of <maxgap|maxmis> arrays
 */
static void
merge_msa(ESL_MSA *mmsa, const ESL_MSA *msa, const int *maxgap, const int *maxmis, int C)
{
  int   prv_nseq = mmsa->nseq;
  int   i,a;
  int   mi;
  char *tmpstr = NULL;
  int   status;

  /* When we first encounter a #=GR PP, SS, or SA aligned residue annotation,
   * allocate in <mmsa>. This only happens once. 
   */
  if (msa->pp && mmsa->pp == NULL)
    {
      ESL_ALLOC(mmsa->pp, sizeof(char *) * mmsa->sqalloc);
      for (mi = 0; mi < mmsa->sqalloc; mi++) mmsa->pp[mi] = NULL;
    }
  if (msa->ss && mmsa->ss == NULL) 
    {
      ESL_ALLOC(mmsa->ss, sizeof(char *) * mmsa->sqalloc);
      for (mi = 0; mi < mmsa->sqalloc; mi++) mmsa->ss[mi] = NULL;
    }
  if (msa->sa && mmsa->sa == NULL) 
    {
      ESL_ALLOC(mmsa->sa, sizeof(char *) * mmsa->sqalloc);
      for (mi = 0; mi < mmsa->sqalloc; mi++) mmsa->sa[mi] = NULL;
    }

  for (i = 0, mi = prv_nseq; i < msa->nseq; i++, mi++)
    {
      esl_msa_SetSeqName(mmsa, mi, msa->sqname[i], -1);
      if (msa->sqacc  && msa->sqacc[i])  esl_msa_SetSeqAccession  (mmsa, mi, msa->sqacc[i],  -1);
      if (msa->sqdesc && msa->sqdesc[i]) esl_msa_SetSeqDescription(mmsa, mi, msa->sqdesc[i], -1);

      if (msa->ngs > 0)
        for (a = 0; a < msa->ngs; a++)
          if (msa->gs[a][i]) esl_msa_AddGS(mmsa, msa->gs_tag[a], -1, mi, msa->gs[a][i], -1);

      inflate_string(msa->aseq[i], msa->rf, msa->alen, maxgap, maxmis, C,
                     /*do_local_rna_tildes=*/ FALSE,
                     /*do_frag_tildes=*/      TRUE,
                     mmsa->alen, &(mmsa->aseq[mi]));

      if (msa->pp && msa->pp[i])  inflate_string(msa->pp[i], msa->rf, msa->alen, maxgap, maxmis, C, FALSE, FALSE, mmsa->alen, &(mmsa->pp[mi]));
      if (msa->ss && msa->ss[i])  inflate_string(msa->ss[i], msa->rf, msa->alen, maxgap, maxmis, C, FALSE, FALSE, mmsa->alen, &(mmsa->ss[mi]));
      if (msa->sa && msa->sa[i])  inflate_string(msa->sa[i], msa->rf, msa->alen, maxgap, maxmis, C, FALSE, FALSE, mmsa->alen, &(mmsa->sa[mi]));

      if (msa->gr != NULL)
        for (a = 0; a < msa->ngr; a++)
          if (msa->gr[a][i] != NULL)
            {
              inflate_string(msa->gr[a][i], msa->rf, msa->alen, maxgap, maxmis, C, FALSE, FALSE, mmsa->alen, &(tmpstr));
              esl_msa_AppendGR(mmsa, msa->gr_tag[a], mi, tmpstr);
              free(tmpstr); 
            }
    }

  mmsa->nseq += msa->nseq;
  return;

 ERROR: //NOTREACHED
  return;  
}



/* inflate_string()
 *
 * Inflate an aseq or GC annotation line, adding extra padding as
 * required for the overall merged MSA.
 *
 * This padding is usually gap characters (`.`), but missing data
 * (`~`) characters may also be needed, in two situations:
 *
 *   1. HMMER and Infernal annotate local alignment fragments in MSAs by
 *      using ~ for terminal gaps. (I don't like this convention any more,
 *      but can't do anything about it at the moment.) This only applies
 *      to aligned sequences, not annotation lines.
 *
 *   2. Infernal annotates RNA local structure alignment (residues treated
 *      as inserted/unaligned to the structure model because of local structure
 *      alignment) using ~ in #=GC RF and #=GC SS_cons annotation lines.
 *      This only applies to those two types of annotation lines, not aseqs
 *      or other annotation lines (not even #=GR SS lines).
 *
 * How this function treats `~` missing data characters is determined
 * by the <do_local_rna_tildes> and <do_frag_tildes> flags. For an
 * aseq that may be using the HMMER|Infernal local sequence fragment
 * convention, set <do_frag_tildes> TRUE. For a #=GC RF or SS_cons
 * annotation that may be using the Infernal local RNA structure
 * alignment convention, set <do_local_rna_tildes> TRUE.
 *
 *   s       - aligned string to inflate, [0..L-1], \0-terminated.
 *   rf      - consensus annotation aligned to s, also [0..L-1], \0-terminated. non-gap/non-missing positions are consensus: there are C of them.
 *   L       - length of both <s>, <rf>
 *   maxgap  - array of lengths of each insert region in the merged MSA, [0..C].  
 *   maxmis  - array of lengths of each missing data ~ region in the merged MSA, [0..C]
 *   C       - number of consensus positions marked in <rf>; dictates the C+1 length of <maxgap|maxmis> arrays
 *   do_local_rna_tildes - TRUE if this is a #=GC RF | SS_cons line, and we're to expand Infernal RNA local alignment ~ regions with ~'s.
 *   do_frag_tildes      - TRUE if this is an aseq that could be a fragment, and fragment ends need to be ~'s.
 *   new_alen - length of new alignment after inflation: caller has already calculated clen + sum_{0..clen} maxgap[] + sum_{0..clen} maxmis[]
 *   ret_as   - RETURN: fully aligned string (allocated here)
 *
 * If necessary, this could be optimized by preprocessing <rf> into
 * 1|0's, which would save a ton of strchr() calls.
 */
static void
inflate_string(const char *s, const char *rf, int L, const int *maxgap, const int *maxmis, int C, int do_local_rna_tildes, int do_frag_tildes, int new_alen, char **ret_as)
{
  char *as      = NULL;
  int   i       = 0;  // position in <s> and <rf>, 0..L-1. len(s) = len(rf) = L <= len(as)
  int   j       = 0;  // position in <as>, 0.. C+sum[maxgap]+sum[maxmis]-1, though we don't need to know that len
  int   c       = 0;  // position in consensus (nongap RF) positions, 0..C-1,C.  c=C is special case: right-trailing nonconsensus positions
  int   n       = 0;  // number of characters copied from <s> in current nonconsensus region <c>
  char  mischar = (do_local_rna_tildes ? '~' : '.');
  char  gapchar = '.';
  int   lb      = 0;
  int   rb      = L;   // init to "this is not a fragment": i > rb is always FALSE because i=0..L-1, i=L on the last iteration of the for (c..) loop
  int   status;
  
  ESL_DASSERT1(( strlen(s) == L && strlen(rf) == L ));

  ESL_ALLOC(as, sizeof(char) * (new_alen+1));

  if (do_frag_tildes && (s[0] == '~' || s[L-1] == '~'))
    {
      for (lb  = 0;   lb  <  L; lb++) if (s[lb] != '~') break;
      for (rb  = L-1; rb  >= 0; rb--) if (s[rb] != '~') break;
    }

  for (c = 0; c <= C; c++)
    {  // on last iteration of this loop, i=L, rf[i] and s[i] may already be '\0'; we'll only add the right terminal gaps, if any
      n = 0;
      while (rf[i] == '~')          { as[j++] = (i < lb || i > rb ? '~' : s[i++]);  n++; }
      while (n < maxmis[c])         { as[j++] = (i < lb || i > rb ? '~' : mischar); n++; }
      n = 0;
      while (rf[i] == '.' || rf[i] == '-' || rf[i] == '_')                                   // can't use strchr, because strchr(".-_", c) for c='\0' is not NULL
                                    { as[j++] = (i < lb || i > rb ? '~' : s[i++]);  n++; }   // copy insertion in <s> that's before consensus pos <c>. Now s[i] is consensus pos <c>.
      while (n < maxgap[c])         { as[j++] = (i < lb || i > rb ? '~' : gapchar); n++; }

      as[j++] = (i < lb || i > rb ? '~' : s[i++]);  // on the final c=C iteration, i=L. On a non-frag, we'll copy the terminal \0 from s[L]. 
    }                                               // but on a fragment, we'll set as[j] to '~' instead of '\0', so we need to fix that.
  as[j-1] = '\0';                                   // <== ... which is why we're doing this.
  ESL_DASSERT1(( i == L+1 ));
  ESL_DASSERT1(( j == new_alen+1 ));
  *ret_as = as;
  return;

 ERROR: //NOTREACHED
  return;
}
