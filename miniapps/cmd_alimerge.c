/* `easel alimerge`: merge MSAs into one, via reference (RF) annotation
 *
 * Contents:
 *    1. esl_cmd_alimerge
 *    2. two code paths: default and --small
 *    3. internal functions for default path (or both)
 *    4. internal functions specific to --small
 *    5. exegesis: notes on how MSAs are merged.
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
#include "esl_vectorops.h"

static ESL_OPTIONS cmd_options[] = {
  /* name         type          default  env   range togs reqs  incomp     help                                                     docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "help; show brief info on version and usage",                  1 },
  { "-o",         eslARG_OUTFILE,  NULL, NULL, NULL, NULL,NULL, NULL,      "output merged MSA to file <f>, not stdout",                   1 },
  { "-v",         eslARG_NONE,    FALSE, NULL, NULL, NULL,"-o", NULL,      "print info about merge to stdout; requires -o",               1 },
  { "--rfonly",   eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "remove all columns that are gaps in GC RF annotation",        1 },
  { "--small",    eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,      "use minimal RAM (make RAM usage independent of MSA sizes)",   1 },
  { "--outformat",eslARG_STRING,   NULL, NULL, NULL, NULL,NULL,"--small",  "write merged output MSA in format <s>",                       1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


static void    default_protocol(char **msafile, int nfile, int outfmt, FILE *ofp, int do_rfonly, int be_verbose, int max_filename_width);
static void    small_protocol  (char **msafile, int nfile,             FILE *ofp, int do_rfonly, int be_verbose, int max_filename_width);

static void    show_verbose_header (char **msafile, int nfile, int *ret_max_filename_width);
static int     update_maxmis_maxgap(const ESL_MSA *msa, int *maxmis, int *maxgap, char *errbuf);
static int     get_msa_consensus_len(const ESL_MSA *msa);
static int *   create_msa_consensus_mask(const ESL_MSA *msa);
static void    transfer_annotation (ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, const int *maxmis, const int *maxgap, int C, int be_verbose);
static void    transfer_threshold_info_double(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, int which, const char *tag, int be_verbose);
static void    transfer_threshold_info_single(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, int which, const char *tag, int be_verbose);
static int     first_gc_ok(const char *gc,  const char *rf);
static int     other_gc_ok(const char *gc1, const char *rf1, const char *gc2, const char *rf2);
static int     other_rf_ok(const char *rf1, const char *rf2);
static void    merge_msa(ESL_MSA *mmsa, const ESL_MSA *msa, const int *maxmis, const int *maxgap, int C);
static void    inflate_string(const char *s, const char *rf, int L, const int *maxmis, const int *maxgap, int C, int do_local_rna_tildes, int do_frag_tildes, int new_alen, char **ret_as);

static void    small_write_msa_top(FILE *ofp, int max_gf_taglen, ESL_MSA *mmsa);
static void    small_determine_gap_inflation(const ESL_MSA *msa, const int *maxmis, const int *maxgap, int C, int *extra_gaps);
static void    small_write_msa_bottom(FILE *ofp, ESL_MSA *mmsa, int max_namelen, int max_gc_taglen, int max_gr_taglen);


/*****************************************************************
 * 1. esl_cmd_alimerge()
 *****************************************************************/

/* esl_cmd_alimerge()
 *
 * Follows standard interface for an Easel miniapp.
 */
int
esl_cmd_alimerge(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS   *go                  = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, /*custom opthelp_f=*/NULL);
  int            nfile               = esl_opt_ArgNumber(go);
  char         **msafile             = NULL;                               // array of MSA filenames; msafile[0..nfile-1]
  int            be_verbose          = esl_opt_GetBoolean(go, "-v");       // with -v, we output a summary table to stdout, and merged MSA to -o <f> file
  int            max_filename_width  = 0;                                  // stays 0 if be_verbose==FALSE. `if (max_filename_width)` and `if (be_verbose)` do the same
  int            do_small            = esl_opt_GetBoolean(go, "--small");
  int            do_rfonly           = esl_opt_GetBoolean(go, "--rfonly");
  char          *outfile             = esl_opt_GetString (go, "-o");
  FILE          *ofp                 = stdout;
  int            outfmt              = (do_small? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM);
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

  if (be_verbose) show_verbose_header(msafile, nfile, &max_filename_width);

  if (! do_small) default_protocol(msafile, nfile, outfmt, ofp, do_rfonly, be_verbose, max_filename_width);
  else            small_protocol  (msafile, nfile,         ofp, do_rfonly, be_verbose, max_filename_width);

  free(msafile);
  if (ofp != stdout) fclose(ofp);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  return status;
}


/*****************************************************************
 * 2. two code paths: default and --small
 *****************************************************************/

/* default_protocol()
 *
 * The standard flow: read all input MSAs into memory, then merge them
 * one at a time into the output MSA. (As opposed to the small memory
 * protocol with --small.)
 *
 * MSA files must be in Stockholm format, because we need the #=GC RF
 * consensus column annotation.  They can be one-block ("Pfam") or
 * multiblock.
 *
 *    msafile             : array of MSA file names to merge. (These are just ptrs into cmdline **argv.)
 *    nfile               : # of MSA files in <msafile>
 *    outfmt              : format of output merged MSA
 *    ofp                 : stream to write the merged MSA to; often stdout, unless <be_verbose>
 *    do_rfonly           : TRUE to remove insertions, only use consensus columns from each MSA
 *    be_verbose          : TRUE to print verbose progress information; in this case ofp != stdout, MSA is going somewhere else
 *    max_filename_width  : width of longest msa name, for formatting <be_verbose> output; or 0 if !be_verbose
 */
static void
default_protocol(char **msafile, int nfile, int outfmt, FILE *ofp, int do_rfonly, int be_verbose, int max_filename_width)
{
  ESL_MSAFILE *afp        = NULL;
  ESL_MSA    **msa        = NULL;    // array of all input MSAs in all files
  ESL_MSA     *mmsa       = NULL;    // output merged MSA
  int          msa_nalloc = 2;       // <msa> is dynamically resized, starting from this
  int          C          = 0;       // consensus length, according to all #=GC RF annotations on all MSAs
  int         *maxmis     = NULL;    // maxmis[0..c..C] = max # of missing data '~' seen before consensus column <c>
  int         *maxgap     = NULL;    // maxgap[0..c..C] = max # of gaps seen before consensus column <c>. maxgap[0] = leading, maxgap[C] = trailing
  int         *useme      = NULL;    // when do_rfonly: useme[0..alen-1] = TRUE | FALSE flags for consensus | insert columns
  int          nseq_tot   = 0;       // total # of seqs over all MSAs
  int          nmsa       = 0;       // total # of MSAs in all files
  int          ai;                   // counter over all MSAs in all files
  int          ai2;                  // counter over MSAs in current file (only needed for error reporting)
  int          fi;                   // counter over MSA files
  int          i;                    // generic loop counter
  char         errbuf[eslERRBUFSIZE];
  int          status;

  /* Initial allocation for array of input MSAs */
  ESL_ALLOC(msa, sizeof(ESL_MSA *) * msa_nalloc);
  for (ai = 0; ai < msa_nalloc; ai++) msa[ai] = NULL;

  /* Pass #1 over the MSA files */
  for (ai = 0, fi = 0; fi < nfile; fi++)
    {
      status = esl_msafile_Open(/*abc=*/NULL, msafile[fi], /*env=*/NULL, eslMSAFILE_STOCKHOLM, NULL, &afp);
      if (status != eslOK) esl_msafile_OpenFailure(afp, status);

      ai2 = 0;  // counter over MSAs in current file; only needed for error reporting
      while ((status = esl_msafile_Read(afp, &(msa[ai]))) == eslOK)
        {
          if (msa[ai]->rf == NULL)
            esl_fatal("All MSAs to merge must have #=GC RF consensus annotation.\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                      ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1);
          
          if (ai == 0) {
            C = get_msa_consensus_len(msa[0]);
            ESL_ALLOC(maxmis, sizeof(int) * (C+1));
            ESL_ALLOC(maxgap, sizeof(int) * (C+1));               
            esl_vec_ISet(maxmis, C+1, 0);
            esl_vec_ISet(maxgap, C+1, 0);
          } else {
            if (! other_rf_ok(msa[0]->rf, msa[ai]->rf) || get_msa_consensus_len(msa[ai]) != C)  // checking C is redundant, since it already passed other_rf_ok, but doesn't hurt
              esl_fatal("All MSAs must have compatible consensus annotation (#=GC RF) lines.\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                        ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1);
          }

          if (do_rfonly)
            {
              useme = create_msa_consensus_mask(msa[ai]);
              esl_msa_ColumnSubset(msa[ai], useme);
              free(useme);
            }
          else
            {
              if ((status = update_maxmis_maxgap(msa[ai], maxmis, maxgap, errbuf)) != eslOK)
                esl_fatal("Problem with consensus annotation line (#=GC RF)\n  in MSA #%d (%s) in file %s (#%d MSA overall)\n  %s",
                          ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1, errbuf);
            }

          nseq_tot += msa[ai]->nseq;

          /* Output to summary table, if we're printing one */
          if (be_verbose)
            {
              char *base_filename;        // tmp var for filename without path prefix

              esl_FileTail(msafile[fi], /*nosuffix=*/FALSE, &base_filename);
              esl_printf("  %7d  %-*s  %7d  %9d  %9" PRId64 "  %13d  %8d\n",
                         (fi+1),  max_filename_width, base_filename,
                         (ai+1),  msa[ai]->nseq,  msa[ai]->alen, nseq_tot,
                         C + esl_vec_ISum(maxmis, C+1) + esl_vec_ISum(maxgap, C+1));
              free(base_filename);
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
   * Now we have all MSAs stored in msa[], and maxmis[],maxgap[] to tell us how to expand each alignment for the merge.
   * <nmsa>     is total number of input MSAs
   * <nseq_tot> is total number of seqs in merged MSA
   */
  mmsa       = esl_msa_Create(nseq_tot, -1);
  mmsa->alen = C + esl_vec_ISum(maxmis, C+1) + esl_vec_ISum(maxgap, C+1);
  /* by setting <alen> *after* creating the <mmsa>, we're responsible for all <alen+1> aligned string allocations */


  /* Per-file and per-column annotations need to be checked for
   * consistency across all individual MSAs.
   */
  transfer_annotation(mmsa, msa, nmsa, maxmis, maxgap, C, be_verbose);
  
  /* Merge input MSAs into output MSA one at a time, free'ing as we go. */
  if (be_verbose) { esl_printf("#\n# Merging alignments ... "); fflush(stdout); }
  for (ai = 0; ai < nmsa; ai++)
    {
      merge_msa(mmsa, msa[ai], maxmis, maxgap, C);
      esl_msa_Destroy(msa[ai]);
    }
  if (be_verbose) { esl_printf("done.\n"); fflush(stdout); }

  /* Output the merged MSA */
  if (be_verbose) { esl_printf("# Writing merged alignment to file ... "); fflush(stdout); }
  esl_msafile_Write(ofp, mmsa, outfmt);
  if (be_verbose) { esl_printf("done.\n"); }
  
  esl_msa_Destroy(mmsa);
  free(maxmis);
  free(maxgap);
  free(msa);
  return;

 ERROR: // UNREACHED
  esl_fatal("allocation failure (but it already failed)");
}


/* small_protocol()
 *
 * Small-memory flow: avoid ever reading an entire MSA into memory,
 * and instead do multiple passes through the MSA files, using the
 * legacy ESL_MSAFILE2 small-memory parser.
 *
 * MSA files must be in "Pfam" format (one-block Stockholm). They all
 * must be files, no stdin streams, because we need to make multiple
 * passes and we can't rewind a stream. It's ok if they're gzip'ed
 * files though.
 *
 * All MSAs must have #=GC RF annotation and the number of annotated
 * consensus columns must be identical in all of them. We merge using
 * this annotation.
 *
 *    msafile            : array of MSA file names to merge. These are just pointers into cmdline **argv. 
 *    nfile              : number of MSA files in <msafile>
 *    outfmt             : format of output merged MSA
 *    ofp                : stream to write the merged MSA to; often stdout, unless <be_verbose>
 *    do_rfonly          : TRUE to remove insertions, only use consensus columns from each MSA
 *    be_verbose         : TRUE to print verbose progress information; in this case ofp != stdout, MSA is going somewhere else
 *    max_filename_width : width of longest msa name, for formatting <be_verbose> output; or 0 if !be_verbose
 */
static void
small_protocol(char **msafile, int nfile, FILE *ofp, int do_rfonly, int be_verbose, int max_filename_width)
{
  ESL_MSAFILE2 *afp2       = NULL;	       // open alignment file, small-mem interface 
  ESL_MSA     **msa        = NULL;             // array of info-only MSAs in all files (GF|GC|comments)
  ESL_MSA      *mmsa       = NULL;             // merged MSA, info-only (GF|GC|comments)
  int           msa_nalloc = 2;                // <msa> is dynamically resized, starting from this
  int           C          = 0;                // consensus length, according to all #=GC RF annotations
  int          *maxmis     = NULL;             // maxmis[0..c..C] = max # of missing data '~' seen before consensus column <c>
  int          *maxgap     = NULL;             // maxgap[0..c..C] = max # of gaps seen before consensus column <c>. maxgap[0] = leading, maxgap[C] = trailing
  int          *alen       = NULL;             // for do_rfonly: need to remember original MSA alens before ColumnSubset(), for Regurgitate()
  int         **useme      = NULL;             // for do_rfonly: useme[0..nmsa-1][0..alen-1]: TRUE|FALSE flags marking consensus|insert columns in each MSA
  int          *extra_gaps = NULL;             // [0..c..C] = # of all-gap cols added after each consensus col of a particular MSA
  int           nseq_tot   = 0;                // total # of seqs over all MSAs
  int           nmsa       = 0;                // total # of MSAs in all files
  int           ai;                            // counter over all MSAs in all files
  int           ai2;                           // counter over MSAs in current file (only needed for error reporting)
  int           fi;                            // counter over MSA files
  int           i;                             // generic loop counter
  int          *nali_per_file;                 // number of MSAs in each file, [0..nfile-1]

  /* sizing info for current MSA: */
  int           nseq_cur;                      // number of aseqs in current MSA
  int64_t       alen_cur;                      // number of columns ""
  int           max_namelen_cur;               // max length of aseq name ""
  int           max_gf_taglen_cur;             // max length of #=GF tag length ""
  int           max_gc_taglen_cur;             // max length of #=GC tag length ""
  int           max_gr_taglen_cur;             // max length of #=GR tag length ""
  int           ngslines_cur;                  // number of GS lines ""
  /* overall sizing info for merged MSA: */
  int           max_namelen     = 0;
  int           max_gf_taglen   = 0;
  int           max_gc_taglen   = 0;
  int           max_gr_taglen   = 0;
  int           has_gs          = FALSE;       // TRUE if at least one MSA has #=GS annotation, necessitating an extra pass
  char          errbuf[eslERRBUFSIZE];
  int           status;
  
  /* Initial allocation for MSAs (info-only, no aseqs or GR annotations) */
  ESL_ALLOC(msa,   sizeof(ESL_MSA *) * msa_nalloc);
  ESL_ALLOC(useme, sizeof(int *)     * msa_nalloc);
  ESL_ALLOC(alen,  sizeof(int)       * msa_nalloc); 
  for (ai = 0; ai < msa_nalloc; ai++) { msa[ai] = NULL; useme[ai] = NULL; alen[ai] = 0; }

  ESL_ALLOC(nali_per_file, sizeof(int) * nfile);
  esl_vec_ISet(nali_per_file, nfile, 0);

  /* Pass #1. 
   * Read and store "info" for each MSA: GF, GC, comments; but not aseqs, GS
   * Determine consensus length C.
   * Allocate maxmis, maxgap arrays for C+1.
   */
  ai = 0;
  for (fi = 0; fi < nfile; fi++)
    {
      status = esl_msafile2_Open(msafile[fi], /*env=*/NULL, &afp2);   // open in text mode; afp2->abc stays NULL. msafile2_Open() requires Pfam format
      if      (status == eslENOTFOUND) esl_fatal("MSA file %s doesn't exist or is not readable\n", msafile[fi]);
      else if (status != eslOK)        esl_fatal("MSA file %s open failed with error code %d\n",   msafile[fi], status);

      /* Read all MSAs in the file - there might be more than one */
      ai2 = 0;   // counter over all MSAs in file; only needed for error reporting
      while (( status = esl_msafile2_ReadInfoPfam(afp2,       /*listfp=*/NULL, /*abc=*/NULL, /*known_alen=*/-1, /*known_rf=*/NULL, /*known_ss_cons=*/NULL,
                                                  &(msa[ai]), 
                                                  &nseq_cur, &alen_cur, &ngslines_cur,
                                                  &max_namelen_cur, &max_gf_taglen_cur, &max_gc_taglen_cur, &max_gr_taglen_cur,
                                                  /*five optional counts EPN uses somewhere:*/ NULL, NULL, NULL, NULL, NULL)) != eslEOF)
        {
          if      (status == eslEFORMAT) esl_fatal("MSA file %s parsing failed:\n%s\n",            msafile[fi], afp2->errbuf);
          else if (status != eslOK)      esl_fatal("MSA file %s read failed with error code %d\n", msafile[fi], status);

          if (msa[ai]->rf == NULL)
            esl_fatal("All MSAs to merge must have #=GC RF consensus annotation.\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                      ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1);
   
          max_namelen   = ESL_MAX(max_namelen,   max_namelen_cur);
          max_gf_taglen = ESL_MAX(max_gf_taglen, max_gf_taglen_cur);
          max_gc_taglen = ESL_MAX(max_gc_taglen, max_gc_taglen_cur);
          max_gr_taglen = ESL_MAX(max_gr_taglen, max_gr_taglen_cur);
          if (ngslines_cur > 0) has_gs = TRUE;
          nseq_tot     += nseq_cur;
          msa[ai]->alen = alen_cur;    // normally we wouldn't keep alen in an info-only MSA, but this is safe here

          /* First MSA: store consensus length C; allocate for <maxmis>, <maxgap> arrays.
           * Subsequent MSAs: their C must match
           */
          if (ai == 0) {
            C = get_msa_consensus_len(msa[0]);
            ESL_ALLOC(maxmis, sizeof(int) * (C+1));
            ESL_ALLOC(maxgap, sizeof(int) * (C+1));               
            esl_vec_ISet(maxmis, C+1, 0);
            esl_vec_ISet(maxgap, C+1, 0);
          } else {
            if (! other_rf_ok(msa[0]->rf, msa[ai]->rf) || get_msa_consensus_len(msa[ai]) != C)  // checking C is redundant, since it already passed other_rf_ok, but doesn't hurt
              esl_fatal("All MSAs to merge must have compatible consensus annotation (#=GC RF).\nMSA #%d (%s) in file %s (#%d MSA overall) does not.",
                        ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1);
          }

          /* With do_rfonly: leave maxmis, maxgap all 0's, inserts will be removed;
           *                 instead, construct and store useme[] for each MSA, for Regurgitate to use on aseqs|GR annotation later;
           *                 and ColumnSubset() each info-only MSA now, which does the subsetting of the GC annotations.
           *
           * Normally:       update maxmis|maxgap counters, which we can do on an info-only MSA because we 
           *                 only need the RF annotation.
           */
          if (do_rfonly)
            {
              useme[ai] = create_msa_consensus_mask(msa[ai]);
              alen[ai]  = msa[ai]->alen;   // remember this before subsetting. Regurgitate needs it.
              esl_msa_ColumnSubset(msa[ai], useme[ai]);
            }
          else
            { 
              if ((status = update_maxmis_maxgap(msa[ai], maxmis, maxgap, errbuf)) != eslOK)
                esl_fatal("Problem with consensus annotation line (#=GC RF)\n  in MSA #%d (%s) in file %s (#%d MSA overall)\n  %s",
                          ai2+1, msa[ai]->name ? msa[ai]->name : "unnamed MSA", msafile[fi], ai+1, errbuf);
            }

          /* Output to summary table, if we're printing one */
          if (be_verbose)
            {
              char *base_filename;  // tmp var for filename without path prefix

              esl_FileTail(msafile[fi], /*nosuffix=*/FALSE, &base_filename);
              esl_printf("  %7d  %-*s  %7d  %9d  %9" PRId64 "  %13d  %8d\n",
                         (fi+1),  max_filename_width, base_filename,
                         (ai+1),  msa[ai]->nseq,  msa[ai]->alen, nseq_tot,
                         C + esl_vec_ISum(maxmis, C+1) + esl_vec_ISum(maxgap, C+1));
              free(base_filename);
            }

          /* Reallocate for more MSAs as needed */
          ai++;
          ai2++;
          if (ai == msa_nalloc)
            {
              msa_nalloc *= 2;
              ESL_REALLOC(msa,   sizeof(ESL_MSA *) * msa_nalloc);
              ESL_REALLOC(useme, sizeof(int *)     * msa_nalloc);
              ESL_REALLOC(alen,  sizeof(int)       * msa_nalloc);
              for (i = ai;  i < msa_nalloc; i++) {
                msa[i]   = NULL;
                useme[i] = NULL;
                alen[i]  = 0;
              }
            }

        } // end loop over reading all MSAs from open <afp2> for msafile[fi] 

      nali_per_file[fi] = ai2;
      if (ai2 == 0) esl_fatal("Failed to read any MSAs from file %s", msafile[fi]);
      esl_msafile2_Close(afp2);
    } // end loop over all msa files
  nmsa = ai;

  /* End of pass #1.
   * Now we have MSA info (GF, GC, comments) stored in <msa>;
   *     maxmis[], maxgap[] tell us how to expand each MSA for the merge;
   *     <nmsa> is the total # of input MSAs to merge;
   * and <nseq_tot> is the total # of seqs in the merged MSA.
   */

  /* Allocate for merged MSA, info-only */
  mmsa       = esl_msa_Create(nseq_tot, -1);
  mmsa->alen = C + esl_vec_ISum(maxmis, C+1) + esl_vec_ISum(maxgap, C+1);
  /* by setting <alen> *after* creating the <mmsa>, we're responsible for all <alen+1> aligned string allocations */

  /* Per-file and per-column (GF, GC, comment) annotations need to be checked for
   * consistency across all individual MSAs.
   */
  transfer_annotation(mmsa, msa, nmsa, maxmis, maxgap, C, be_verbose);
 
  /* Write the header, merged comments, and merged #=GF information */
  small_write_msa_top(ofp, max_gf_taglen, mmsa);

  if (be_verbose) { esl_printf("# Writing merged alignment to file ... "); fflush(stdout); }

  /* Pass #2: an extra pass is needed if any of the MSA files contain #=GS annotations.
   */
  if (has_gs)
    {
      for (fi = 0; fi < nfile; fi++)
        {
          status = esl_msafile2_Open(msafile[fi], /*env=*/NULL, &afp2);   // open in text mode; afp2->abc stays NULL
          if      (status == eslENOTFOUND) esl_fatal("MSA file %s doesn't exist or is not readable (at extra pass)\n", msafile[fi]);
          else if (status != eslOK)        esl_fatal("MSA file %s open failed with error code %d (at extra pass)\n",   msafile[fi], status);

          while ( ( status = esl_msafile2_RegurgitatePfam(afp2, ofp, 
                                                          max_namelen, max_gf_taglen, max_gc_taglen, max_gr_taglen, /* max width of a seq name, gf tag, gc tag, gr tag (irrelevant here) */
                                                          FALSE, FALSE, FALSE, FALSE, FALSE, // no header, trailer, blank lines, comments, GF yet
                                                          TRUE,                              // only #=GS lines now
                                                          FALSE, FALSE, FALSE,               // no GC, GR, aseq
                                                          NULL, NULL,                        // no seqs2regurg or seqs2skip keyhash downselection
                                                          NULL, NULL,                        // no useme or add2me column downselect or addition
                                                          -1,                                // no need to provide <alen>
                                                          '.',                               // gap char is irrelevant, but we have to give one
                                                          NULL, NULL)) != eslEOF)            // not asking for nseq_read or nseq_regurg on return
            {
              if      (status == eslEFORMAT) esl_fatal("MSA file %s parsing failed (at extra pass):\n%s\n",            msafile[fi], afp2->errbuf);
              else if (status != eslOK)      esl_fatal("MSA file %s read failed (at extra pass) with error code %d\n", msafile[fi], status);
            }  // end loop over MSAs in msafile <ai>
          esl_msafile2_Close(afp2);
        } // end loop over all MSA files
      fprintf(ofp, "\n"); /* a single blank line to separate GS annotation from aligned data */
    }


  /* Pass #3: finally, the aseqs - realigned w/ legacy msafile2 regurgitation interface.
   */
  ai = 0;
  ESL_ALLOC(extra_gaps, sizeof(int) * (mmsa->alen+1));
  for (fi = 0; fi < nfile; fi++)
    {
      status = esl_msafile2_Open(msafile[fi], /*env=*/NULL, &afp2);   // open in text mode; afp2->abc stays NULL
      if      (status == eslENOTFOUND) esl_fatal("MSA file %s doesn't exist or is not readable (at final pass)\n", msafile[fi]);
      else if (status != eslOK)        esl_fatal("MSA file %s open failed with error code %d (at final pass)\n",   msafile[fi], status);

      for (ai2 = 0; ai2 < nali_per_file[fi]; ai2++)
        {
          if (! do_rfonly)
            small_determine_gap_inflation(msa[ai], maxmis, maxgap, C, extra_gaps);

          status = esl_msafile2_RegurgitatePfam(afp2, ofp,
                                                max_namelen, max_gf_taglen, max_gc_taglen, max_gr_taglen,  // overall field widths
                                                FALSE,         // no stockholm header 
                                                FALSE,         // no // trailer
                                                FALSE,         // no blank lines 
                                                FALSE,         // no comments
                                                FALSE,         // no GF annotation
                                                FALSE,         // no GS annotation
                                                FALSE,         // no GC annotation
                                                TRUE,          // regurgitate and align GR annotations
                                                TRUE,          // regurgitate and align aseqs
                                                NULL,          // regurgitate all seqs, not a subset 
                                                NULL,          // regurgitate all seqs, not a subset 
                                                (do_rfonly ? useme[ai] : NULL),          // do_rfonly: remove inserts altogether
                                                (do_rfonly ? NULL      : extra_gaps),    // overall realignment uses this information
                                                (do_rfonly ? alen[ai]  : msa[ai]->alen), // alignment length, as we read it in first pass (do_rfonly needs this to be the original alen, not C)
                                                '.',           // gap character to use
                                                NULL,          // don't return num seqs read 
                                                NULL);         // don't return num seqs regurgitated 
          if      (status == eslEFORMAT) esl_fatal("MSA file %s regurgitation failed at msa %d (at final pass):\n%s\n",               msafile[fi], ai2+1, afp2->errbuf);
          else if (status != eslOK)      esl_fatal("MSA file %s regurgitation failed at msa %d (at final pass) with error code %d\n", msafile[fi], ai2+1, status);
        }
      ai++;
      esl_msafile2_Close(afp2);
    }

  /* Finally, the GC annotation; we already aligned it in the info-only <mmsa> in transfer_annotation()
   */
  small_write_msa_bottom(ofp, mmsa, max_namelen, max_gc_taglen, max_gr_taglen);
  if (be_verbose) { esl_printf("done.\n"); }

  /* Cleanup 
   */
  for (ai = 0; ai < nmsa; ai++) { esl_msa_Destroy(msa[ai]); free(useme[ai]); }
  free(msa);
  free(useme);
  free(alen);
  esl_msa_Destroy(mmsa);
  free(nali_per_file);
  free(maxmis);
  free(maxgap);
  free(extra_gaps);
  return;
  
 ERROR: // UNREACHED
  esl_fatal("allocation failure (but it already failed)");
}


/***************************************************************** 
 * 3. internal functions for default path (or both)
 *****************************************************************/


/* show_verbose_header()
 *
 * When <be_verbose> is used, we're writing progress information to
 * <stdout> while writing the merged MSA to an output file. Here we
 * start that progress information by writing the top of a tabular
 * output for each input MSA file.
 *
 * Args:      msafile                : array of <nfile> input MSA file names, [0..nfile-1]
 *            nfile                  : length of <msafile> array
 *            ret_max_filename_width : RETURN: strlen of longest file name in <msafile>.
 *                                     (Used for table formatting later.)
 */
static void
show_verbose_header(char **msafile, int nfile, int *ret_max_filename_width)
{
  int   max_filename_width  = 9;
  char *namedashes          = NULL;
  char *base_filename       = NULL;
  int   fi, i;
  int   status;

  /* Determine longest filename in the <msafile> list */
  for (fi = 0; fi < nfile; fi++)
    {
      esl_FileTail(msafile[fi], /*nosuffix=*/FALSE, &base_filename);
      max_filename_width = ESL_MAX(max_filename_width, strlen(base_filename));
      free(base_filename);
    }

  /* Variable-length table underscoring for the filename field */
  ESL_ALLOC(namedashes, sizeof(char) * (max_filename_width+1));
  for (i = 0; i < max_filename_width; i++) namedashes[i] = '-';
  namedashes[max_filename_width] = '\0';

  /* Format the header for the verbose output table */
  esl_printf("# Reading %d alignment files\n", nfile);
  esl_printf("# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "",        max_filename_width, "",          "",        "",          "",            "cumulative",    "ncols");
  esl_printf("# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "file #",  max_filename_width, "file name", "ali #",   "#seq/ali",  "ncols/ali",   "# seq total",   "required");
  esl_printf("# %7s  %*s  %7s  %9s  %9s  %13s  %8s\n",  "-------", max_filename_width, namedashes,  "-------", "---------", "---------",   "-------------", "--------");

  free(namedashes);
  *ret_max_filename_width = max_filename_width;
  return;

 ERROR: // NOTREACHED
  esl_fatal("allocation error");
}


/* update_maxmis_maxgap()
 * 
 * The maxmis[] and maxgap[] arrays keep track of maximum number of
 * gap columns and missing data columns annotated in #=GC RF reference
 * annotation over all MSAs.
 *
 * For each <msa> in turn, with <msa->rf> annotation, given the
 * current <maxgap> and <maxmis>, update them as needed.
 *
 * There are <C> consensus columns, indexed 1..C. maxmis[0..C]/maxgap[0..C]
 * are the max number of ~/gaps after column c; maxmis[0]/maxgap[0]
 * are for the leading left side.
 *
 * Missing data ~ columns must precede gap columns between any two consensus columns
 * in the RF annotation (including leading and trailing nonconsensus columns). 
 * If this is not true, an <eslEINVAL> error is returned; if caller provides
 * an <errbuf> allocated for at least <eslERRBUFSIZE> chars, an informative error
 * message is put there if this failure happens.
 */
static int
update_maxmis_maxgap(const ESL_MSA *msa, int *maxmis, int *maxgap, char *errbuf)
{
  int apos, cpos;
  int ngap = 0;
  int nmis = 0;

  for (apos = 0, cpos = 0; apos <= msa->alen; apos++)  // on last iteration, apos = alen and rf[alen] is on its '\0' terminator;
    {                                                  //  we'll set maxgap and maxmis for trailing flank.
      if (msa->rf[apos] == '~')
        {
          nmis++;
          if (ngap) ESL_FAIL(eslEINVAL, errbuf, "gaps precede missing data ~ in RF annotation at col=%d", apos+1);
        }
      else if (memchr("-_.", msa->rf[apos], 3))
        ngap++;
      else // includes rf[alen] = `\0` to terminate
        {
          maxmis[cpos] = ESL_MAX(maxmis[cpos], nmis);
          maxgap[cpos] = ESL_MAX(maxgap[cpos], ngap);
          cpos++;
          ngap = nmis = 0;
        }
    } 
  return eslOK;
}
  
/* get_msa_consensus_len()
 *
 * Given an <msa> with <msa->rf> annotation, return the
 * number of consensus positions marked in that annotation.
 */
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

/* create_msa_consensus_mask()
 *
 * Create and return a <useme[0..alen-1]> array of TRUE|FALSE flags
 * marking consensus columns, based on nongap/gap characters in
 * <msa->rf> annotation.
 */
static int *
create_msa_consensus_mask(const ESL_MSA *msa) 
{
  int *useme = NULL; 
  int  apos;
  int  status;

  ESL_ALLOC(useme, sizeof(int) * msa->alen);

  for (apos = 0; apos < msa->alen; apos++)
    useme[apos] = (memchr("-_.~", msa->rf[apos], 4) ? FALSE : TRUE);
  return useme;

 ERROR: // NOTREACHED
  esl_fatal("allocation failed");
}
  

/* transfer_annotation()
 *
 * Only per-file and per-column annotations that are identical across
 * all MSAs is transferred to the output merged MSA.
 */
static void
transfer_annotation(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, const int *maxmis, const int *maxgap, int C, int be_verbose)
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

  /* GA|TC|NC score thresholds. 
   * Pfam has one. Rfam has two.
   */
  if ( (msa[0]->cutset[eslMSA_TC1] && msa[0]->cutset[eslMSA_TC2]) ||
       (msa[0]->cutset[eslMSA_NC1] && msa[0]->cutset[eslMSA_NC2]) ||
       (msa[0]->cutset[eslMSA_GA1] && msa[0]->cutset[eslMSA_GA2]))
    {
      transfer_threshold_info_double(mmsa, msa, nmsa, eslMSA_TC1, "TC", be_verbose);
      transfer_threshold_info_double(mmsa, msa, nmsa, eslMSA_NC1, "NC", be_verbose);
      transfer_threshold_info_double(mmsa, msa, nmsa, eslMSA_GA1, "GA", be_verbose);
    }
  else if ( msa[0]->cutset[eslMSA_TC1] ||
            msa[0]->cutset[eslMSA_NC1] ||
            msa[0]->cutset[eslMSA_GA1]) 
    {
      transfer_threshold_info_single(mmsa, msa, nmsa, eslMSA_TC1, "TC", be_verbose);
      transfer_threshold_info_single(mmsa, msa, nmsa, eslMSA_NC1, "NC", be_verbose);
      transfer_threshold_info_single(mmsa, msa, nmsa, eslMSA_GA1, "GA", be_verbose);
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
        inflate_string(msa[0]->ss_cons, msa[0]->rf, msa[0]->alen, maxmis, maxgap, C, /*do_local_rna_tildes=*/TRUE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->ss_cons));
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
        inflate_string(msa[0]->sa_cons, msa[0]->rf, msa[0]->alen, maxmis, maxgap, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->sa_cons));
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
        inflate_string(msa[0]->pp_cons, msa[0]->rf, msa[0]->alen, maxmis, maxgap, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &(mmsa->pp_cons));
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
            if (i2 == msa[ai]->ngc || ! other_gc_ok(msa[0]->gc[i],  msa[0]->rf, msa[ai]->gc[i], msa[ai]->rf))
              break;
          }
        } else ai=0;

        if (ai == nmsa) {
          ESL_ALLOC(tmpstr, sizeof(char) * (mmsa->alen+1));
          inflate_string(msa[0]->gc[i], msa[0]->rf, msa[0]->alen, maxmis, maxgap, C, /*do_local_rna_tildes=*/FALSE, /*do_frag_tildes=*/FALSE, mmsa->alen, &tmpstr);
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
  inflate_string(msa[0]->rf, msa[0]->rf, msa[0]->alen, maxmis, maxgap, C, /*do_local_rna_tildes=*/TRUE, FALSE, mmsa->alen, &(mmsa->rf));
  if (be_verbose)
    esl_printf("# Identical RF annotation from all alignments transferred to merged alignment.\n");
  return;

 ERROR: //NOTREACHED
  return;
}


static void
transfer_threshold_info_double(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, int which, const char *tag, int be_verbose)
{
  int which2 = which + 1;   // this assumes e.g. eslMSA_GA2 = eslMSA_GA1+1, which is documented in esa_msa.h
  int ai;
  
  if (msa[0]->cutset[which] && msa[0]->cutset[which2]) {
    for (ai = 1; ai < nmsa; ai++)
      if (! msa[ai]->cutset[which] || ! msa[ai]->cutset[which2] ||
          esl_FCompare(msa[ai]->cutoff[which],  msa[0]->cutoff[which],  /*r_tol*/1e-5, /*a_tol*/1e-5) != eslOK ||
          esl_FCompare(msa[ai]->cutoff[which2], msa[0]->cutoff[which2], /*r_tol*/1e-5, /*a_tol*/1e-5) != eslOK)
        break;
    if (ai == nmsa) {
      mmsa->cutoff[which]  = msa[0]->cutoff[which];   mmsa->cutset[which]  = TRUE;
      mmsa->cutoff[which2] = msa[0]->cutoff[which2];  mmsa->cutset[which2] = TRUE;
    }
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# %s1|%s2:     absent (at least from first MSA); not included in merge\n",  tag, tag);
    else if (ai < nmsa) esl_printf("# %s1|%s2:     not identical in all MSAs, not transferred to merged MSA\n", tag, tag);
    else                esl_printf("# %s1|%s2:     identical in all MSAs, transferred to merged MSA\n",         tag, tag);
  }
}

static void
transfer_threshold_info_single(ESL_MSA *mmsa, ESL_MSA **msa, int nmsa, int which, const char *tag, int be_verbose)
{
  int ai;
  
  if (msa[0]->cutset[which]) {
    for (ai = 1; ai < nmsa; ai++)
      if (! msa[ai]->cutset[which] || esl_FCompare(msa[ai]->cutoff[which],  msa[0]->cutoff[which],  /*r_tol*/1e-5, /*a_tol*/1e-5) != eslOK)
        break;
    if (ai == nmsa) {
      mmsa->cutoff[which]  = msa[0]->cutoff[which];
      mmsa->cutset[which]  = TRUE;
    }
  } else ai = 0;
  if (be_verbose) {
    if      (ai == 0)   esl_printf("# %s1:         absent (at least from first MSA); not included in merge\n",  tag, tag);
    else if (ai < nmsa) esl_printf("# %s1:         not identical in all MSAs, not transferred to merged MSA\n", tag, tag);
    else                esl_printf("# %s1:         identical in all MSAs, transferred to merged MSA\n",         tag, tag);
  }
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

  while (1)
    {
      while (memchr("-_.~", rf2[j], 4)) { if (! memchr("-_.~", gc2[j], 4)) return FALSE; j++; }  // advance j to next consensus position, while checking for nongap GC annotation
      while (memchr("-_.~", rf1[i], 4)) i++;                                                     // advance i to next consensus position 
      if      (gc1[i] != gc2[j]) return FALSE;   // check for gc identity at consensus positions. Includes a final comparison of \0:\0 at ends of strings.
      else if (gc2[j] == '\0')   return TRUE;    // gc1[i] must also be \0 here, and so are rf1[i],rf2[j]
      i++;
      j++;
    }
}

static int
other_rf_ok(const char *rf1, const char *rf2)
{
  int i = 0;
  int j = 0;

  while (1) 
    {
      while (memchr("-_.~", rf2[j], 4)) j++;    // memchr not strchr because we don't want terminator \0 to match when we hit it in rf2[j]
      while (memchr("-_.~", rf1[i], 4)) i++;  
      if      (rf1[i] != rf2[j]) return FALSE;   
      else if (rf2[j] == '\0')   return TRUE;   // rf1[i] also must be \0 here
      i++;
      j++;
    }
}




/* merge_msa()
 *
 * Merge one input MSA into the new merged output MSA, given precalculated
 * arrays <maxmis> and <maxgap> that tell us how many padding gaps we need
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
 * <maxmis>  : array of lengths of each missing data ~ region in the merged MSA, [0..C]
 * <maxgap>  : array of lengths of each insert region in the merged MSA, [0..C].  
 * <C>       : number of consensus positions marked in <rf>; dictates the C+1 length of <maxgap|maxmis> arrays
 */
static void
merge_msa(ESL_MSA *mmsa, const ESL_MSA *msa, const int *maxmis, const int *maxgap, int C)
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

      inflate_string(msa->aseq[i], msa->rf, msa->alen, maxmis, maxgap, C,
                     /*do_local_rna_tildes=*/ FALSE,
                     /*do_frag_tildes=*/      TRUE,
                     mmsa->alen, &(mmsa->aseq[mi]));

      if (msa->pp && msa->pp[i])  inflate_string(msa->pp[i], msa->rf, msa->alen, maxmis, maxgap, C, FALSE, FALSE, mmsa->alen, &(mmsa->pp[mi]));
      if (msa->ss && msa->ss[i])  inflate_string(msa->ss[i], msa->rf, msa->alen, maxmis, maxgap, C, FALSE, FALSE, mmsa->alen, &(mmsa->ss[mi]));
      if (msa->sa && msa->sa[i])  inflate_string(msa->sa[i], msa->rf, msa->alen, maxmis, maxgap, C, FALSE, FALSE, mmsa->alen, &(mmsa->sa[mi]));

      if (msa->gr != NULL)
        for (a = 0; a < msa->ngr; a++)
          if (msa->gr[a][i] != NULL)
            {
              inflate_string(msa->gr[a][i], msa->rf, msa->alen, maxmis, maxgap, C, FALSE, FALSE, mmsa->alen, &(tmpstr));
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
 *   maxgap  - array of lengths of each missing data ~ region in the merged MSA, [0..C]
 *   maxmis  - array of lengths of each insert region in the merged MSA, [0..C].  
 *   C       - number of consensus positions marked in <rf>; dictates the C+1 length of <maxmis|maxgap> arrays
 *   do_local_rna_tildes - TRUE if this is a #=GC RF | SS_cons line, and we're to expand Infernal RNA local alignment ~ regions with ~'s.
 *   do_frag_tildes      - TRUE if this is an aseq that could be a fragment, and fragment ends need to be ~'s.
 *   new_alen - length of new alignment after inflation: caller has already calculated clen + sum_{0..clen} maxmis[] + sum_{0..clen} maxgap[]
 *   ret_as   - RETURN: fully aligned string (allocated here)
 *
 * If necessary, this could be optimized by preprocessing <rf> into
 * 1|0's, which would save a ton of strchr() calls.
 */
static void
inflate_string(const char *s, const char *rf, int L, const int *maxmis, const int *maxgap, int C, int do_local_rna_tildes, int do_frag_tildes, int new_alen, char **ret_as)
{
  char *as       = NULL;           // newly allocated aligned string, as[0..new_alen-1, \0]
  int    j       = 0;              // position in as[0..new_alen-1]
  int    i       = 0;              // position of current consensus column in s[i], rf[i] (or 0 at start); then advances past it when we copy a consensus position
  int    c       = 0;              // index of last consensus column advanced to: 0,1..C,C+1
  int    ip      = -1;             // i coord of the previous consensus column we saw, such that insert block = [ip+1, i-1]
  char   mischar = (do_local_rna_tildes ? '~' : '.');  // ~ chars for Infernal local structure ali are only used in #=GC SS_cons and RF lines
  char   gapchar = '.';            // all other new gaps we need to open in as[] are '.'
  int    lb      = 0;              // position of first non-'~' in s[], for fragment marking; 0 = no left ~'s
  int    rb      = L;              // position of last non-'~' in s[]; L = no right ~'s
  int    nmis;                     // # of ~ chars in local structure alignment part of current insert block, according to rf[]
  int    ngap;                     // # of gap chars in insert part of current insert block, according to rf[]
  int    x,i2;
  int    status;

  ESL_DASSERT1(( strlen(s) == L && strlen(rf) == L ));

  ESL_ALLOC(as, sizeof(char) * (new_alen+1));

  /* find first and last non-'~' chars in s[], if as[] is an aligned seq that may need to be fragment-marked */
  if (do_frag_tildes && (s[0] == '~' || s[L-1] == '~'))
    {
      for (lb  = 0;   lb  <  L; lb++) if (s[lb] != '~') break;
      for (rb  = L-1; rb  >= 0; rb--) if (s[rb] != '~') break;
    }

  while (c <= C)
    {
      /* Advance i on s and rf to next consensus position rf[i] (or rf[L] = '\0' at end)
       * The previous insertion block is then s[ip+1..i-1]
       *     with ~ part :     [ip+1,      ip+nmis]
       *     and  . part :     [ip+nmis+1, i-1]
       *
       * <nmis> is Infernal-specific. Only Infernal MSAs have nmis > 0, maxmis[]> 0.
       */
      nmis = ngap = 0;
      while (rf[i] == '~')            { nmis++; i++; }
      while (memchr(".-_", rf[i], 3)) { ngap++; i++; }   // strchr() would match on `\0` and we don't want that
      c++;

      if (c == 1)  // Leading insertion block is right-flush
        {          // Infernal MSAs can't have maxmis[0] > 0: local structure ali can only have unaligned insert between two consensus cols
          ESL_DASSERT1( (maxmis[0] == 0) );
          for (x  = 0; x < maxgap[0] - ngap; x++)  { as[j++] = (0  < lb ? '~' : gapchar); }
          for (i2 = 0; i2 <= i-1;            i2++) { as[j++] = (i2 < lb ? '~' : s[i2]);   }
        }
      else if (c == C+1)  // Trailing insert block is left-flush
        {                 // Infernal MSAs can't have maxmis[C] > 0, see above
          ESL_DASSERT1( (maxmis[c-1] == 0) );
          for (i2 = ip+1; i2 <= ip + nmis + ngap;  i2++) { as[j++] = (i2 > rb ? '~' : s[i2]);   }
          for (x  = 0;    x  < maxgap[c-1] - ngap; x++)  { as[j++] = (i2 > rb ? '~' : gapchar); }
        }
      else
        {               // All other insert blocks are split, left- and right-justified; expansion happens in middles
          for (i2 = ip+1; i2 <= ip + nmis/2;        i2++)  { as[j++] = (i2 < lb || i2 > rb ? '~' : s[i2]);   }
          for (x  = 0;    x  <  maxmis[c-1] - nmis; x++)   { as[j++] = (i2 < lb || i2 > rb ? '~' : mischar); }
          for (      ;    i2 <= ip + nmis;          i2++)  { as[j++] = (i2 < lb || i2 > rb ? '~' : s[i2]);   }

          for (      ; i2 <= ip + nmis + ngap/2; i2++)     { as[j++] = (i2 < lb || i2 > rb ? '~' : s[i2]);   }
          for (x  = 0; x  <  maxgap[c-1] - ngap; x++)      { as[j++] = (i2 < lb || i2 > rb ? '~' : gapchar); }
          for (      ; i2 <= i-1;                i2++)     { as[j++] = (i2 < lb || i2 > rb ? '~' : s[i2]);   }
        }

      ip = i;
      as[j++] = s[i++];     // copy this consensus position ... including the terminal '\0' 
    }   
  ESL_DASSERT1(( i == L+1 ));
  ESL_DASSERT1(( j == new_alen+1 ));
  *ret_as = as;
  return;

 ERROR: //NOTREACHED   (because a failed alloc was already immediately fatal above)
  return;
}




/*****************************************************************
 * 4. internal functions specific to --small 
 *****************************************************************/

/* small_write_msa_top()
 *
 * Write the top of the output merged MSA <mmsa>:
 * Stockholm header, comments, and #=GF per-file annotation.  
 *
 * <mmsa> is an info-only <mmsa> containing #=GF, #=GC, and comments.
 *
 * The code is cribbed from esl_msa.c::actually_write_stockholm().
 */
static void
small_write_msa_top(FILE *ofp, int max_gf_taglen, ESL_MSA *mmsa)
{
  int i;
  if (max_gf_taglen < 2) max_gf_taglen = 2;  // Always need to accommodate standard tags like ID|DE|AC|AU

  /* Magic Stockholm header */
  esl_fprintf(ofp, "# STOCKHOLM 1.0\n");

  /* Free text comments */
  for (i = 0;  i < mmsa->ncomment; i++)
    esl_fprintf(ofp, "# %s\n", mmsa->comment[i]);
  if (mmsa->ncomment > 0) esl_fprintf(ofp, "\n");

  /* GF section: per-file annotation */
  if (mmsa->name) esl_fprintf(ofp, "#=GF %-*s %s\n", max_gf_taglen, "ID", mmsa->name);
  if (mmsa->acc)  esl_fprintf(ofp, "#=GF %-*s %s\n", max_gf_taglen, "AC", mmsa->acc);
  if (mmsa->desc) esl_fprintf(ofp, "#=GF %-*s %s\n", max_gf_taglen, "DE", mmsa->desc);
  if (mmsa->au)   esl_fprintf(ofp, "#=GF %-*s %s\n", max_gf_taglen, "AU", mmsa->au);
  
  /* Thresholds are hacky. Pfam has two. Rfam has one. */
  if      (mmsa->cutset[eslMSA_GA1] && mmsa->cutset[eslMSA_GA2])
    esl_fprintf(ofp, "#=GF %-*s %.1f %.1f\n", 
                max_gf_taglen, "GA", mmsa->cutoff[eslMSA_GA1], mmsa->cutoff[eslMSA_GA2]);
  else if (mmsa->cutset[eslMSA_GA1])
    esl_fprintf(ofp, "#=GF %-*s %.1f\n", 
                max_gf_taglen, "GA", mmsa->cutoff[eslMSA_GA1]);

  if      (mmsa->cutset[eslMSA_NC1] && mmsa->cutset[eslMSA_NC2])
    esl_fprintf(ofp, "#=GF %-*s %.1f %.1f\n",
                max_gf_taglen, "NC", mmsa->cutoff[eslMSA_NC1], mmsa->cutoff[eslMSA_NC2]);
  else if (mmsa->cutset[eslMSA_NC1])
    esl_fprintf(ofp, "#=GF %-*s %.1f\n",
                max_gf_taglen, "NC", mmsa->cutoff[eslMSA_NC1]);

  if      (mmsa->cutset[eslMSA_TC1] && mmsa->cutset[eslMSA_TC2])
    esl_fprintf(ofp, "#=GF %-*s %.1f %.1f\n",
                max_gf_taglen, "TC", mmsa->cutoff[eslMSA_TC1], mmsa->cutoff[eslMSA_TC2]);
  else if (mmsa->cutset[eslMSA_TC1])
    esl_fprintf(ofp, "#=GF %-*s %.1f\n", 
                max_gf_taglen, "TC", mmsa->cutoff[eslMSA_TC1]);

  for (i = 0; i < mmsa->ngf; i++)
    esl_fprintf(ofp, "#=GF %-*s %s\n", max_gf_taglen, mmsa->gf_tag[i], mmsa->gf[i]); 

  esl_fprintf(ofp, "\n");
  return;
}


/* small_determine_gap_inflation()
 *
 * For a given <msa>, determine the number of extra all-gap columns
 * that need to be opened up on this MSA to align it to all the
 * others. Store this in <extra_gaps[0..msa->alen]>. Caller provides
 * the <maxmis[0..C]> and <maxgap[0..C]> arrays, which give the number
 * of gap/~ columns after each consensus column 1..C in the final
 * merged MSA (with 0 being prepended ones).
 *
 * <extra_gaps[i]> is the number of extra all-gap columns to add after
 * column <i> of the <msa>, where the columns in a digital MSA are
 * numbered 1..L. <extra_gaps[0]> is how many are prepended.
 *
 * The extra gap columns are added to the right of existing insertion
 * blocks (to the immediate left of consensus columns), keeping
 * insertion blocks left-flushed. There are C+1 possible places where
 * extra gap columns may be added, for C consensus columns.
 *
 * This number of extra gap columns is maxmis[c] + maxgap[c] - # of
 * residues after consensus position c, and it can be determined
 * entirely from the RF annotation line on the msa, without needing to
 * look at any aseqs.
 *
 * This information will then be provided by the caller to the legacy
 * esl_msafile2 regurgitation interface.
 *
 * <extra_gaps[0..alen]> is allocated by the caller. It is sufficient
 * for the caller to allocate once for the final mmsa->alen+1 alignment
 * length; this is guaranteed to be >= the largest alen of all input
 * MSAs.
 *
 * (The coordinate systems are confusing; notes above have been
 * written carefully. In the digital MSA and in the extra_gaps[]
 * result, columns are indexed i=1..alen, but the MSA's RF annotation
 * string is indexed 0..alen-1, unfortunately.  In maxmis[] and
 * maxgap[], consensus columns are indexed c=1..C.)
 *
 * Returns: <eslOK> on success; extra_gaps[0..alen] is the result.
 *
 * Throws:  <eslEMEM> on allocation failure. <extra_gaps> is indeterminate.
 */
static void
small_determine_gap_inflation(const ESL_MSA *msa, const int *maxmis, const int *maxgap, int C, int *extra_gaps)
{
  int  i = 0;   // index on MSA columns: msa->rf[0..alen-1,alen] and extra_gaps[0..alen]
  int  c = 0;   // index on consensus columns: maxmis[0..C] and maxmat[0..C]
  int  nmis, ngap;

  esl_vec_ISet(extra_gaps, msa->alen+1, 0);

  while (c <= C)
    {
      nmis = ngap = 0;
      while (msa->rf[i] == '~')            { nmis++; i++; }
      while (memchr(".-_", msa->rf[i], 3)) { ngap++; i++; }
      c++;

      if      (c == 1)   extra_gaps[0]         = maxgap[0] - ngap;
      else if (c == C+1) extra_gaps[msa->alen] = maxgap[C] - ngap;
      else
        {
          extra_gaps[i - ngap - (nmis+1)/2]  = maxmis[c-1] - nmis;
          extra_gaps[i -        (ngap+1)/2] += maxgap[c-1] - ngap;
        }
      i++; // advance past RF consensus position
    }
}

/* small_write_msa_bottom()
 *
 * Write the aligned GC annotation, and the final "//".
 */
static void
small_write_msa_bottom(FILE *ofp, ESL_MSA *mmsa, int max_namelen, int max_gc_taglen, int max_gr_taglen)
{
  int margin;
  int i;

  margin = max_namelen + 1;
  if (max_gc_taglen) margin = ESL_MAX(margin, max_gc_taglen + 6);                // "#=GC <longest_tag> "
  if (max_gr_taglen) margin = ESL_MAX(margin, max_namelen + max_gr_taglen + 7);  // "#=GR <longest_name> <longest_tag> "

  if (mmsa->ss_cons) esl_fprintf(ofp, "#=GC %-*s %s\n", margin-6, "SS_cons", mmsa->ss_cons); 
  if (mmsa->sa_cons) esl_fprintf(ofp, "#=GC %-*s %s\n", margin-6, "SA_cons", mmsa->sa_cons); 
  if (mmsa->pp_cons) esl_fprintf(ofp, "#=GC %-*s %s\n", margin-6, "PP_cons", mmsa->pp_cons); 
  if (mmsa->rf)      esl_fprintf(ofp, "#=GC %-*s %s\n", margin-6, "RF",      mmsa->rf); 
  for (i = 0; i < mmsa->ngc; i++)
    esl_fprintf(ofp, "#=GC %-*s %s\n", margin-6, mmsa->gc_tag[i], mmsa->gc[i]); 

  esl_fprintf(ofp, "//\n");
}


/***************************************************************** 
 * 5. exegesis: notes on how MSAs are merged.
 *****************************************************************
 *
 * `easel alimerge` accomodates two idiosyncrasies of Stockholm format
 * MSAs generated by HMMER and Infernal.
 *
 * One is the use of missing data `~` characters as special gap
 * symbols, in two circumstances:
 *
 *   1. Sequence fragments may be annotated using ~ for terminal
 *      gaps. (I don't like this convention any more, but can't do
 *      anything about it at the moment.) This only applies to aligned
 *      sequences, not annotation lines.
 *
 *   2. Infernal annotates RNA local structure alignment (residues
 *      treated as inserted/unaligned to the structure model) using ~
 *      in RF and SS_cons annotation lines. This only applies to those
 *      two annotation lines, not aseqs or other annotations (not even
 *      #=GR SS lines).
 *
 * A second idiosyncrasy is that insertions are considered unaligned,
 * and are split-justified. For m inserted residues in an insertion
 * subblock, m/2 are flushed left and (m+1)/2 are flushed right. For
 * example, for m=5, 2 residue goes left, 3 residues go right.
 *
 * Local structure insertions (marked with ~ in the RF and SS_cons
 * annotation) are considered to be a distinct kind of insertion, as
 * opposed to normal insertions. This means that insert blocks in the
 * MSA are composed of two subblocks: first the local structure
 * insertions, then the normal insertions. Local structure insertions
 * can only occur between consensus columns, not in the leading or
 * trailing left and right flanks.
 *
 * MSAs are handled in `easel alimerge` in text mode, so aseqs are
 * 0..alen-1, and so is the (required) msa->rf consensus annotation.
 * We consider the consensus columns to be numbered 1..C. There are
 * C+1 insert blocks, numbered 0..C: block 0 is the leading left
 * flank, and 1..C follow each consensus position.
 *
 * An illustrative example of an Infernal alignment of phage tRNAs,
 * where I've changed insert gaps from - to . and inserted residues to
 * lower case for clarity (which is also what HMMER does):
 *
 * seq1         ~CCGUAGUAGUGUAGU..GGU..A...............ACACGACACC.UUGCCAA.............GGUGUAA...................U.UGCG.AGUUCGAUU.........CUCGUCUACGG~~
 * seq2         ~CCGUAGUGGUGUAGU..GGU..A...............ACAUACUUCC.UUGCCAA.............GGAAGCG...................U.UGCG.AGUUCGAUU.........CUCGUCUACGG~~
 * seq3         GCGUGGAUAGCUCAGU..GGCcgA...............GAGCACCCGA.C------gc........aa.UCGGGGG...................A.CGCG.GGUUCAAAU.........CCCGCUCCA~~~~
 * seq4         GGACCUCUAGCUCAGA..AGC..A...............GAGCACCCGC.CUCUUAA.............GCGGGGG...................GuCGGG.A--------uggcagaau-UCCCGGGGUCCA
 * seq5         GCCCCUAUAGCUCAGU.uGGU..A...............GAGCUAUCGC.CUUUUAA.............GCGACAG...................GuCGCA.GGUUCGAGU.........C~~~~~~~~~~~~
 * seq6         ~CGGUUAUCGUUUAAC..GGA..A...............AAACAAUUCA.CUCAUAA.............UGAAUCA...................C.UCCG.GGUUCGAGU.........CCCGGU~~~~~~~
 * seq7         ~CGGUUAUCGUUUAAC..GGA..A...............AAACAAUUCA.CUCAUAA.............UGAAUCA...................C.UCCG.GGUUCGAGU.........CCCGGU~~~~~~~
 * seq8         GGCGCGGUAGGCAAU-..GGA..A...............GACCACUGAUGCUCAGAA............cAUCAGAU...................GcUGCG.GGUUCGAAU.........CCCGUCUGCGCUA
 * seq9         GCCACGUUGGAGUAGUcuGGU.uA...............UCUCGUCAGA.UUCUCAA.............UCUGAAG...................AcCACG.GGUUCGAAU.........CCC~~~~~~~~~~
 * seq10        AGUCGCAUAGUUUAAC.uGGU.gA...............GAACAGC---.-------uuccc.uguggg.---GCGA...................G.UGUC.GGUUCGAAU.........CCGGCUGCGACGC
 * seq11        UCCCCCGUGG------..---..-uacgcugguuugcca---C-CCUGA.CUUUGAA.............UCAGGAA...................CaCGUA.GGUUCGAUC.........CCUGCCGGGGGAG
 * seq12        UGCUACGUAGCUUAAUauGGC..A...............AAGCAA----.-------cucuuggccagg.----GUG...................U.CACGuGGUUCGAAC.........CCCGUCGUAGCAG
 * seq13        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CGCCGGA.CUUUGAC.............UCCGGAG...................G.UGCA.GGUUCGAGU.........CCUGCUACCCC~~
 * seq14        ~~~~~~~~~~~GUAAA.uGGC..A...............ACACACUUCC.CUCGUAA.............GGAAGAA...................U.UCUC.GGUUCGAAU.........CCGAG~~~~~~~~
 * seq15        ~~~~~~~~~~~GUAAA.uGGC..A...............ACACACUUCC.CUCGUAA.............GGAAGAA...................U.UCUC.GGUUCGAAU.........CCGAG~~~~~~~~
 * seq16        GCCCUAGUA-------..---..-uc...........gc----CACUGG.CUUCUAA.............CCAGUC-guuaacagcguaauggaguA.UGCA.GGUUCGAGU.........CCUGUCUAGGGUA
 * seq17        UACGGAGUAGGGGAGUucGGA.gU...............CCCCGCUGGU.CUCAUAC.............GCCAGA-ggcuuaca..aacgcgcauA.AGUG.G-UUCAAAU.........CCCGCCUCCGUAA
 * seq18        GCAUGUAAGC------..---..-caggga..uaagcua---GGCGGGC.CUGUAAA.............GCCCUAC...................CuCGAU.GGUUCGAAU.........CCGUCUACAUGCA
 * #=GC SS_cons (((((((,,<<<<___..___.._~~~~~~~~~~~~~~~>>>>,<<<<<._______~~~~~~~~~~~~.>>>>>,,...................,.,<<<.<<_______~~~~~~~~~>>>>>))))))):
 * #=GC RF      ggggaugUAGCuuAaU..GGU..A~~~~~~~~~~~~~~~aaGCauugGa.CUuauAA~~~~~~~~~~~~.uCcaaAg...................g.cgcg.GGUUCgAaU~~~~~~~~~CCcgccauccccA
 */


