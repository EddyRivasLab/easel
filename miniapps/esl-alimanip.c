/* Manipulate a multiple sequence alignment in some useful ways.
 *
 * Derived from easel's esl-alistat which was from squid's alistat (1995)
 * EPN, Fri Aug 10 08:52:30 2007
 * SVN $Id$
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "easel.h"
#include "esl_distance.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_msa.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_wuss.h"

static char banner[] = "manipulate a multiple sequence alignment file";
static char usage[]  = "[options] <msafile>\n\
The <msafile> must be in Stockholm format.";

#define OTHERMSAOPTS  "--merge,--morph,--map,--submap"          /* Exclusive choice for options involving an additional msa */
#define CLUSTOPTS     "--cn-id,--cs-id,--cx-id,--cn-ins,--cs-ins,--cx-ins" /* Exclusive choice for clustering */
#define CHOOSESEQOPTS "--seq-k,--seq-r,--seq-ins,--seq-del" /* Exclusive choice for choosing which seqs to keep/remove */

static int  keep_or_remove_rf_gaps(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int keep_flag, int remove_flag);
static int  write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, float gapthresh);
static int  write_rf_given_alen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *rfmask);
static int  write_rf_given_rflen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *rfmask);
static int  write_rf_given_useme(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int *useme);
static int  morph_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **newmsa1);
static int  merge_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_merged_msa);
static int  add_gap_columns_to_msa(char *errbuf, ESL_MSA *msa, int *toadd, ESL_MSA **ret_msa, int do_treat_as_rf_gap);
static int  cp_and_add_gaps_to_aseq(char *new_aseq, char *orig_aseq, int alen, int *toadd, int nnew, char gapchar);
static int  is_flush_left(int *ngaps, int astart, int aend);
static int  is_flush_right(int *ngaps, int astart, int aend);
static int  pick_gappiest_columns(int *ngaps, int astart, int aend, int nkeep, int **ret_cols_to_keep);
static int  get_gaps_per_column(ESL_MSA *msa, int **ret_ngaps);
static int  map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen);
static int  individualize_consensus(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  read_sqfile(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc, int nseq, ESL_SQ ***ret_sq);
static int  trim_msa(ESL_MSA *msa, ESL_SQ **sq, char *errbuf);
static int  dump_insert_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  dump_residue_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  dump_cres_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  dump_delete_info(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  plot_inserts(FILE *fp, ESL_MSA *msa, int do_log, char *errbuf);
static int  dump_infocontent(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  plot_gaps(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  get_tree_order(ESL_TREE *T, char *errbuf, int **ret_order);
static int  reorder_msa(ESL_MSA *msa, int *order, char *errbuf);
static int  dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask);
static int  map_msas(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, char **ret_msa1_to_msa2_map);
static int  map_sub_msas(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, char **ret_msa1_to_msa2_mask);
static int  handle_post_opts(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  output_rf_as_mask(FILE *fp, char *errbuf, ESL_MSA *msa);
static int  expand_msa2mask(char *errbuf, ESL_MSA *msa1, char *xmask, ESL_MSA **newmsa1);
static int  compare_ints(const void *el1, const void *el2);
static int  msa_median_length(ESL_MSA *msa);
static int  msa_remove_seqs_below_minlen(ESL_MSA *msa, float minlen, ESL_MSA **ret_new_msa);
static int  msa_remove_truncated_seqs(ESL_MSA *msa, char *errbuf, int ntrunc, ESL_MSA **ret_new_msa);
static int  number_columns(ESL_MSA *msa, int do_all, char *errbuf);
static char digit_to_char(int digit);
static int  int_ndigits(int i);
static char get_char_digit_x_from_int(int i, int place);
static int  keep_contiguous_column_block(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
static int  read_seq_name_file(char *filename, char *errbuf, char ***ret_seqlist, int *ret_seqlist_n);
static int  msa_keep_or_remove_seqs(ESL_MSA *msa, char *errbuf, char **seqlist, int seqlist_n, int do_keep, ESL_MSA **ret_new_msa);
static int  insert_x_diffmx(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int do_length_weight, int do_only_internal_inserts, ESL_DMATRIX **ret_D);
static int  insert_x_pair_shared(ESL_MSA *msa, int i, int j, int cfirst, int clast, double *opt_pshared, int *opt_nshared, int *opt_nins);
static int  insert_x_pair_shared_length(ESL_MSA *msa, int i, int j, int cfirst, int clast, double *opt_pshared, double *opt_nshared, int *opt_nins);
static int  MSADivide(ESL_MSA *mmsa, ESL_DMATRIX *D, int do_mindiff, int do_nc, int do_nsize, float mindiff, int target_nc, int target_nsize, int *ret_num_msa, ESL_MSA ***ret_cmsa, int *ret_xsize, char *errbuf);
static int  select_node(ESL_TREE *T, double *diff, double mindiff, int **ret_clust, int *ret_nc, int *ret_xsize, int *ret_best, char *errbuf);
static float find_mindiff(ESL_TREE *T, double *diff, int do_nsize, int target, int **ret_clust, int *ret_nc, int *ret_xsize, int *ret_best, float *ret_mindiff, char *errbuf);
static int  determine_first_last_consensus_columns(ESL_MSA *msa, char *errbuf, int **ret_fA, int **ret_lA, int *ret_clen);
static int  dst_nongap_XPairId(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2, double *opt_distance, int *opt_nid, int *opt_n);
static int  dst_nongap_XDiffMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D);
static int  find_seqs_with_given_insert(ESL_MSA *msa, char *errbuf, int target, int min, int max, int **ret_useme);
static int  minorize_msa(const ESL_GETOPTS *go, ESL_MSA *msa, char *errbuf, FILE *fp, char *tag);
static int  remove_gc_markup(ESL_MSA *msa, char *errbuf, char *tag);

static ESL_OPTIONS options[] = {
  /* name          type        default  env   range      togs reqs  incomp                      help                                                       docgroup */
  { "-h",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "help; show brief info on version and usage",                     1 },
  { "-o",          eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output the alignment to file <f>, not stdout",                   1 },
  { "-1",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "output alignment in Pfam (non-interleaved, 1 line/seq) format",  1 },
  { "--list",      eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output list of sequence names in alignment to file <f>",         1 },
  { "--devhelp",   eslARG_NONE,  NULL,  NULL, NULL,      NULL,NULL, NULL,                       "show list of undocumented developer options",                    1 },
  { "-g",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "add/rewrite #=GC RF markup based on gap frequency in each col",  2 },
  { "--gapthresh", eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,NULL, NULL,                       "with -g, fraction of gaps allowed in non-gap RF columns [0.5]",  2 },
  { "--mask-all",  eslARG_INFILE,FALSE,NULL, NULL,      NULL,NULL, NULL,                        "set #=GC RF as x=1, gap=0 from 1/0s in 1-line <f> (len=alen)",   2 },
  { "--mask-rf",   eslARG_INFILE, FALSE,NULL, NULL,      NULL,NULL, NULL,                       "set #=GC RF as x=1, gap=0 from 1/0s in 1-line <f> (len=rf len)", 2 },
  { "--pfract",    eslARG_REAL,  NULL,  NULL, "0<=x<=1", NULL,NULL, NULL,                       "set #=GC RF as cols w/<x> fraction of seqs w/POST >= --pthresh", 2 },
  { "--pthresh",   eslARG_REAL,  "0.9", NULL, "0<=x<=1", NULL,"--pfract", NULL,                 "set #=GR POST threshold for --pfract as <x> [default=0.9]",      2 },
  { "--p-rf",      eslARG_NONE,  NULL,  NULL, NULL,      NULL,"--pfract", NULL,                 "with --pfract options, ignore gap #=GC RF columns",              2 },
  { "-k",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "keep  only columns w/(possibly post -g) non-gap #=GC RF markup", 3 },
  { "-r",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "remove all columns w/(possibly post -g) non-gap #=GC RF markup", 3 },
  { "--start-all", eslARG_INT,   NULL,  NULL, NULL,      NULL,"--end-all",  "--start-rf",       "keep columns starting at column <n>", 3 },
  { "--end-all",   eslARG_INT,   NULL,  NULL, NULL,      NULL,"--start-all","--start-rf",       "keep columns ending   at column <n>", 3 },
  { "--start-rf",  eslARG_INT,   NULL,  NULL, NULL,      NULL,"--end-rf",   "--start-all",      "keep columns starting at non-gap RF column <n>", 3 },
  { "--end-rf",    eslARG_INT,   NULL,  NULL, NULL,      NULL,"--start-rf", "--start-all",      "keep columns ending   at non-gap RF column <n>", 3 },
  { "--rm-gc",     eslARG_STRING,NULL,  NULL, NULL,      NULL,NULL, NULL,                       "remove GC <s> markup, <s> must be RF|SS_cons|SA_cons|PP_cons", 3},
  { "--tree",      eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL,OTHERMSAOPTS,                "reorder MSA to tree order following SLC, save Newick tree to <f>", 4 },
  { "--lfract",    eslARG_REAL,  NULL,  NULL, "0<=x<=1", NULL,NULL, NULL,                       "remove sequences w/length < <x> fraction of median length",      4 },
  { "--lmin",      eslARG_INT,   NULL,  NULL, "n>0",     NULL,NULL, NULL,                       "remove sequences w/length < <n> residues",                       4 },
  { "--detrunc",   eslARG_INT,   NULL,  NULL, "n>0",     NULL,NULL, NULL,                       "remove seqs w/gaps in >= <n> 5' or 3'-most non-gap #=GC RF cols",4 },
  { "--seq-r",     eslARG_INFILE,NULL,  NULL, NULL,      NULL,NULL, CHOOSESEQOPTS,              "remove sequences with names listed in file <f>",                 4 },
  { "--seq-k",     eslARG_INFILE,NULL,  NULL, NULL,      NULL,NULL, CHOOSESEQOPTS,              "remove all seqs *except* those listed in <f>, reorder seqs also",4 },
  { "--seq-ins",   eslARG_INT,   NULL,  NULL, NULL,      NULL,NULL, CHOOSESEQOPTS,              "keep only seqs w/an insert after non-gap RF col <n>",            4 },
  { "--seq-del",   eslARG_INT,   NULL,  NULL, NULL,      NULL,NULL, CHOOSESEQOPTS,              "keep only seqs w/a  delete in non-gap RF col <n>",               4 },
  { "--seq-ni",    eslARG_INT,    "1",  NULL, "n>0",     NULL,"--seq-ins", NULL,                "w/--seq-ins require at least <n> residue insertions",            4 },
  { "--seq-xi",    eslARG_INT,"1000000",NULL, "n>0",     NULL,"--seq-ins", NULL,                "w/--seq-ins require at most  <n> residue insertions",            4 },
  { "--trim",      eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "trim aligned seqs in <msafile> to subseqs in <f>",               4 },
  { "--iinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,                "print info on # of insertions b/t all non-gap RF cols to <f>",  5 },
  { "--icinfo",    eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on information content of each non-gap RF column",    5 },
  { "--rinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on # of residues in each col of alignment to <f>",    5 },
  { "--cresinfo",  eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on # of columns with 1 residue due to each seq",      5 },
  { "--dinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on # of deletes in non-gap RF cols of aln to <f>",    5 },
  { "--pinfo",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "print info on posterior probabilities in <msafile> to <f>",      5 },
  { "--sindi",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, "-g,-k,-r,--morph",         "annotate individual secondary structures by imposing consensus", 7 },
  { "--num-all",   eslARG_NONE,   NULL, NULL, NULL,      NULL,NULL, NULL,                       "add annotation numbering all columns",                          11 },
  { "--num-rf",    eslARG_NONE,   NULL, NULL, NULL,      NULL,NULL, NULL,                       "add annotation numbering the non-gap RF columns",               11 },
  { "--omask",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output RF annotation as 1/0 mask to file <f>",                   9 },
  { "--amino",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--dna,--rna",               "<msafile> contains protein alignments",                         10 },
  { "--dna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--rna",             "<msafile> contains DNA alignments",                             10 },
  { "--rna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--dna",             "<msafile> contains RNA alignments",                             10 },

  /* All options below are developer options, only shown if --devhelp invoked */
  { "--iplot",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL,OTHERMSAOPTS,                "plot heatmap of # of insertions b/t all non-gap RF cols to <f>", 101 },
  { "--ilog",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,"--iplot", NULL,                  "w/--iplot, use log scale for heatmap of insert counts",          101 },
  { "--gplot",     eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "plot checkerboard grid of # of gaps in non-gap RF cols to <f>",  101 },
  { "--morph",     eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "morph msa in <msafile> to msa in <f>'s gap structure",           101 },
  { "--merge",     eslARG_INFILE,FALSE, NULL, NULL,      NULL,NULL, "--morph,-g,-k,-r",         "merge msa in <msafile> with msa in <f>",                         101 },

  { "--map",       eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "map msa in <msafile> to msa in <f>, output mask (1s and 0s)",    102 },
  { "--submap",    eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, OTHERMSAOPTS,               "map msa in <msafile> to msa in <f> (<f> is subaln of <msafile>", 102 },
  { "--omap",      eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "with --map/--submap, output file for 1/0 mask map is <f>",       102 },
  { "--xmask",     eslARG_INFILE, NULL, NULL, NULL,      NULL,NULL, NULL,                       "for each 0 column in <f>, add a 100% gap column to <msafile>",   102 },
  { "--verbose",   eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "be verbose (usually with --morph, --merge or --map)",            102 },

  { "--cn-id",      eslARG_INT,   NULL,   NULL, "n>0",    NULL,NULL, CLUSTOPTS,                 "split MSA into <n> clusters based on sequence identity",          103 },
  { "--cs-id",      eslARG_INT,   NULL,   NULL, "n>0",    NULL,NULL, CLUSTOPTS,                 "split MSA into clusters on id s.t max cluster has <n> seqs",      103 },
  { "--cx-id",      eslARG_REAL,  NULL,   NULL, "0.<x<1.",NULL,NULL, CLUSTOPTS,                 "split MSA into clusters s.t. no seq b/t 2 clusters > <x> seq id", 103},
  { "--cn-ins",     eslARG_INT,   NULL,   NULL, "n>0",    NULL,NULL, CLUSTOPTS,                 "split MSA into <n> clusters based on insert similarity",          103 },
  { "--cs-ins",     eslARG_INT,   NULL,   NULL, "n>0",    NULL,NULL, CLUSTOPTS,                 "split MSA into clusters on inserts s.t. max cluster has <n> seqs",103 },
  { "--cx-ins",     eslARG_REAL,  NULL,   NULL, "0.<x<1.",NULL,NULL, CLUSTOPTS,                 "split MSA into clusters s.t. no seq b/t 2 clusters > <x> ins id", 103},
  { "--c-nmin",     eslARG_INT,   NULL,   NULL, "n>0",    NULL,NULL, NULL,                      "only keep the cluster(s) with number of seqs > <n>",              103},
  { "--c-mx",       eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                      "output identity matrix to file <f>",                              103 },
  { "-M",           eslARG_STRING,NULL,  NULL, NULL,      NULL,NULL, "--seq-r,--seq-k",         "use #=GS tag <s> to define minor alignments, and output them",    103 },
  { "--M-rf",       eslARG_NONE,  NULL,  NULL, NULL,      NULL,"-M", NULL,                      "w/-M, impose major #=GC RF onto all minor alns",                  103 },

  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  ESL_ALPHABET *abc     = NULL;	/* biological alphabet             */
  char         *alifile = NULL;	/* alignment file name             */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *msa     = NULL;	/* one multiple sequence alignment */
  int           status;		/* easel return code               */
  int           nali;		/* number of alignments read       */
  FILE         *ofp;		/* output file (default is stdout) */
  char          errbuf[eslERRBUFSIZE*4];
  int           write_ali = FALSE; /* set to TRUE if we should print a new MSA */
  /* --merge, --morph, --map, --submap related vars */
  ESL_MSAFILE  *otherafp = NULL;	/* other input alignment file (with --morph) */
  ESL_MSA      *othermsa = NULL;	/* other input alignment      (with --morph) */
  /* --trim related vars */
  ESL_SQFILE   *trimfp = NULL;  /* sequence file with subsequences for --trim */
  /* --iinfo, --iplot, --gplot --rinfo, --dinfo related vars */
  FILE *treefp  = NULL;  /* output file for --tree */
  FILE *iinfofp = NULL;  /* output file for --iinfo */
  FILE *iplotfp = NULL;  /* output file for --iplot */
  FILE *gplotfp = NULL;  /* output file for --gplot */
  FILE *rinfofp = NULL;  /* output file for --rinfo */
  FILE *cresinfofp = NULL; /* output file for --cresinfo */
  FILE *dinfofp = NULL;  /* output file for --dinfo */
  FILE *icinfofp = NULL; /* output file for --icinfo */
  FILE *listfp = NULL;   /* output file for --list */
  /* --mask-all */
  char *amask = NULL;
  /* --mask-all */
  char *rfmask = NULL;
  /* --xmask */
  char *xmask = NULL;
  /* --map, --omap */
  FILE *omapfp;            /* output file for --omap */
  char *msa1_to_msa2_mask; /* the map from <msafile> to <f> from --map, a 1/0 mask */
  /* --omask */
  FILE *omaskfp;

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "--devhelp") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for adding/rewriting #=GC RF annotation:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for removing columns:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions for numbering columns:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 
      puts("\noptions for reordering/removing/trimming sequences:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for displaying info on inserts/gaps/posterior probabilities:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noptions for manipulating secondary structure annotation:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      puts("\noptions for outputting a lanemask file:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for specifying input alphabet:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\nundocumented, experimental developer options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\noptions for comparison/modification based on another MSA file:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80); 
      puts("\noptions for partitioning MSA into clusters:");
      esl_opt_DisplayHelp(stdout, go, 103, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for adding/rewriting #=GC RF annotation:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for removing columns:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions for numbering columns:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 
      puts("\noptions for reordering/removing/trimming sequences:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for displaying info on inserts/gaps/posterior probabilities:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noptions for manipulating secondary structure annotation:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      puts("\noptions for outputting a lanemask file:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for specifying input alphabet:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  alifile = esl_opt_GetArg(go, 1);

  fmt             = eslMSAFILE_STOCKHOLM;

  /***********************************************
   * Open the MSA file; determine alphabet; set for digital input
   ***********************************************/

  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if (status == eslENOTFOUND) 
    esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT) 
    esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK) 
    esl_fatal("Alignment file open failed with error %d\n", status);

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	esl_fatal("Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else ofp = stdout;

  if      (esl_opt_GetBoolean(go, "--amino"))   abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     abc = esl_alphabet_Create(eslRNA);
  else {
    int type;
    status = esl_msafile_GuessAlphabet(afp, &type);
    if (status == eslEAMBIGUOUS)    esl_fatal("Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", alifile);
    else if (status == eslEFORMAT)  esl_fatal("Alignment file parse failed: %s\n", afp->errbuf);
    else if (status == eslENODATA)  esl_fatal("Alignment file %s is empty\n", alifile);
    else if (status != eslOK)       esl_fatal("Failed to read alignment file %s\n", alifile);
    abc = esl_alphabet_Create(type);
  }
  esl_msafile_SetDigital(afp, abc);
  if((esl_opt_GetBoolean(go, "--sindi")) && (abc->type != eslRNA && abc->type != eslDNA))
    esl_fatal("--sindi option pertains to base pairs and only makes sense with DNA or RNA alphabets.");

  /* optionally, open --morph, --merge or --map msa file for reading, --merge, --morph and --map are all incompatible
   * with each other, so we'll never try to do open othermsafile more than once.
   */
  if(esl_opt_GetString(go, "--morph") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--morph"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    esl_fatal("--morph alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of --morph alignment %s\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status != eslOK)      esl_fatal("Alignment file open failed with error %d\n", status);
    }
  if(esl_opt_GetString(go, "--merge") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--merge"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    esl_fatal("--merge alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of --merge alignment %s\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status != eslOK)      esl_fatal("Alignment file open failed with error %d\n", status);
    }
  if(esl_opt_GetString(go, "--map") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--map"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    esl_fatal("--map alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--map"));
      else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of --map alignment %s\n", 
					      esl_opt_GetString(go, "--map"));
      else if (status != eslOK)      esl_fatal("Alignment file open failed with error %d\n", status);
    }
  if(esl_opt_GetString(go, "--submap") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--submap"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    esl_fatal("--submap alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--submap"));
      else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of --submap alignment %s\n", 
					      esl_opt_GetString(go, "--submap"));
      else if (status != eslOK)      esl_fatal("Alignment file open failed with error %d\n", status);
    }

  /* read --mask-all file, if nec */
  if(esl_opt_GetString(go, "--mask-all") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-all"), errbuf, &amask)) != eslOK)
      esl_fatal("--mask-all input file: %s open failed.\n", esl_opt_GetString(go, "--mask-all"));
  }
  /* read --mask-rf file, if nec */
  if(esl_opt_GetString(go, "--mask-rf") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-rf"), errbuf, &rfmask)) != eslOK)
      esl_fatal("--mask-rf input file: %s open failed.\n", esl_opt_GetString(go, "--mask-rf"));
  }
  /* read --xmask file, if nec */
  if(esl_opt_GetString(go, "--xmask") != NULL) {
    if((status = read_mask_file(esl_opt_GetString(go, "--xmask"), errbuf, &xmask)) != eslOK)
      esl_fatal("--xmask input file: %s open failed.\n", esl_opt_GetString(go, "--xmask"));
  }
  /***********************************************
   * Read MSAs one at a time.
   ***********************************************/

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      /* first handle the --lfract option if enabled, all subsequent manipulations will omit any short seqs removed here */
      if (esl_opt_IsOn(go, "--lfract")) {
	int median   = msa_median_length(msa);
	float minlen = esl_opt_GetReal(go, "--lfract") * (float) median;
	ESL_MSA *new_msa;
	msa_remove_seqs_below_minlen(msa, minlen, &new_msa);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }

      /* handle the --lmin option if enabled, all subsequent manipulations will omit any short seqs removed here */
      if (esl_opt_IsOn(go, "--lmin")) {
	float minlen = esl_opt_GetInteger(go, "--lmin");
	ESL_MSA *new_msa;
	msa_remove_seqs_below_minlen(msa, minlen, &new_msa);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }

      /* handle the --detrunc option if enabled, all subsequent manipulations will omit any seqs removed here */
      if( esl_opt_IsOn(go, "--detrunc")) {
	ESL_MSA *new_msa;
	if((status =msa_remove_truncated_seqs(msa, errbuf, esl_opt_GetInteger(go, "--detrunc"), &new_msa)) != eslOK) esl_fatal(errbuf);
	/* new_msa is msa without seqs below minlen, swap ptrs */
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }

      /* handle the --seq-k and --seq-r options if enabled, all subsequent manipulations will omit any seqs removed here */
      if ( esl_opt_IsOn(go, "--seq-k") || esl_opt_IsOn(go, "--seq-r")) {
	ESL_MSA *new_msa;
	char   **seqlist;
	int      seqlist_n, n;
	if( esl_opt_IsOn(go, "--seq-k")) { 
	  if((status = read_seq_name_file(esl_opt_GetString(go, "--seq-k"), errbuf, &seqlist, &seqlist_n)) != eslOK) esl_fatal(errbuf);	  
	  if((status = msa_keep_or_remove_seqs(msa, errbuf, seqlist, seqlist_n, TRUE, &new_msa)) != eslOK)        esl_fatal(errbuf);	  
	  /* new_msa is msa but only with seqs listed in --seq-k <f> file */
	}
	else { /* --seq-r enabled */
	  if((status = read_seq_name_file(esl_opt_GetString(go, "--seq-r"), errbuf, &seqlist, &seqlist_n)) != eslOK) esl_fatal(errbuf);	  
	  if((status = msa_keep_or_remove_seqs(msa, errbuf, seqlist, seqlist_n, FALSE, &new_msa)) != eslOK)        esl_fatal(errbuf);	  
	  /* new_msa is msa but without seqs listed in --seq-r <f> file */
	}
	esl_msa_Destroy(msa);
	msa = new_msa;
	for(n = 0; n < seqlist_n; n++) free(seqlist[n]); 
	free(seqlist);
	write_ali = TRUE;
      }
      
      /* handle the --seq-ins, --seq-del options if enabled, all subsequent manipulations will omit any seqs removed here */
      if(  esl_opt_IsOn(go, "--seq-ins") || esl_opt_IsOn(go, "--seq-del")) {
	ESL_MSA *new_msa;
	int     *useme;
	if( esl_opt_IsOn(go, "--seq-ins")) { 
	  if((status = find_seqs_with_given_insert(msa, errbuf, esl_opt_GetInteger(go, "--seq-ins"), esl_opt_GetInteger(go, "--seq-ni"), esl_opt_GetInteger(go, "--seq-xi"), &useme)) != eslOK) esl_fatal(errbuf);	  
	  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK)  esl_fatal(errbuf);	  
	  /* new_msa is msa but without seqs that do not have an insert of length <a>..<b> (from --seq-ni <a> and --seq-xi <b>) after consensus column <n> from --seq-ins <n> file */
	}
	esl_msa_Destroy(msa);
	msa = new_msa;
	write_ali = TRUE;
      }      

      /* read other msa if --morph, --merge, or --map (which are incompatible with each other) is enabled */
      if(((esl_opt_GetString(go, "--morph") != NULL) || (esl_opt_GetString(go, "--merge") != NULL))
	 || ((esl_opt_GetString(go, "--map") != NULL) || (esl_opt_GetString(go, "--submap") != NULL)))
	{
	  if ((status = esl_msa_Read(otherafp, &othermsa)) != eslOK) {
	    if(status == eslEFORMAT) 
	      esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
			otherafp->linenumber, otherafp->fname, otherafp->errbuf, otherafp->buf);	
	    else if (status == eslEOF)
	      esl_fatal("No alignments read in %s.", esl_opt_GetString(go, "--morph"));
	  }
	}

      /* if nec, handle --trim option */
      if(esl_opt_GetString(go, "--trim") != NULL) { 
	/* open seq file for --trim */
	status = esl_sqfile_Open(esl_opt_GetString(go, "--trim"), eslSQFILE_UNKNOWN, NULL, &(trimfp));
	if (status == eslENOTFOUND)    esl_fatal("File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEINVAL)  esl_fatal("Can’t autodetect stdin or .gz."); 
	else if (status != eslOK)      esl_fatal("Sequence file open failed with error %d\n", status);
	/* read the sequences */
	ESL_SQ **sq;
	read_sqfile(trimfp, msa->abc, msa->nseq, &sq); /* dies on failure */
	/* trim the msa */
	if((status = trim_msa(msa, sq, errbuf)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }

      /* if nec, morph <msafile> into gap structure in <f> (from --morph <f>)*/
      if(esl_opt_GetString(go, "--morph") != NULL)
	{
	  ESL_MSA *newmsa;
	  if((status = morph_msa(go, errbuf, msa, othermsa, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  /*status = esl_msa_Write(stdout, othermsa, eslMSAFILE_STOCKHOLM);
	    status = esl_msa_Write(stdout, newmsa, eslMSAFILE_STOCKHOLM);*/
	  msa = newmsa;
	}

      /* if nec, merge <msafile> and <f> (from --merge <f>) */
      if(esl_opt_GetString(go, "--merge") != NULL)
	{
	  ESL_MSA *newmsa;
	  if((status = merge_msa(go, errbuf, msa, othermsa, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  /*status = esl_msa_Write(stdout, othermsa, eslMSAFILE_STOCKHOLM);
	    status = esl_msa_Write(stdout, newmsa, eslMSAFILE_STOCKHOLM);*/
	  msa = newmsa;
	}

      /* rewrite RF annotation, if nec */
      if(esl_opt_GetBoolean(go, "-g")) {
	if((status = write_rf_gapthresh(go, errbuf, msa, esl_opt_GetReal(go, "--gapthresh"))) != eslOK) goto ERROR;
	write_ali = TRUE;
      }
      if(amask != NULL) { /* --mask-all enabled */
	if((status = write_rf_given_alen(go, errbuf, msa, amask)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }
      if(rfmask != NULL) { /* --mask-rf enabled */
	if((status = write_rf_given_rflen(go, errbuf, msa, rfmask)) != eslOK) goto ERROR;
	write_ali = TRUE;
      }

      /* handle posterior (--p*) options, if nec */
      if( esl_opt_IsOn(go, "--pfract") && ! esl_opt_IsOn(go, "--pinfo")) {
	if((status = handle_post_opts(go, errbuf, msa) != eslOK)) goto ERROR;
	if( esl_opt_IsOn(go, "--pfract")) write_ali = TRUE;
      }

      /* Remove columns based on --start-all --end-all, --start-rf --end-rf, if nec */
      if( esl_opt_IsOn(go, "--start-all") || esl_opt_IsOn(go, "--start-rf"))
	{
	  if((status = keep_contiguous_column_block(go, errbuf, msa) != eslOK)) goto ERROR;
	  write_ali = TRUE;
	}

      /* keep or remove columns based on RF annotation, if nec */
      if(esl_opt_GetBoolean(go, "-k") || esl_opt_GetBoolean(go, "-r"))
	{
	  if((status = keep_or_remove_rf_gaps(go, errbuf, msa, 
					      esl_opt_GetBoolean(go, "-k"),
					      esl_opt_GetBoolean(go, "-r"))) != eslOK) goto ERROR;
	  write_ali = TRUE;
	}

      /* if nec, map <msafile> to <f> (from --map <f>)
       * this is purposefully done after RF annotation is potentially rewritten and
       * some columns are potentially removed. 
       */
      if(esl_opt_GetString(go, "--map") != NULL)
	{
	  if((status = map_msas(go, errbuf, msa, othermsa, &msa1_to_msa2_mask)) != eslOK) goto ERROR;
	  if(esl_opt_GetString(go, "--omap") != NULL) { 
	    if ((omapfp = fopen(esl_opt_GetString(go, "--omap"), "w")) == NULL) 
	      esl_fatal("Failed to open --omap output file %s\n", esl_opt_GetString(go, "--omap"));
	    fprintf(omapfp, "%s\n", msa1_to_msa2_mask);
	    fclose(omapfp);
	    /*printf("# Mask of 1/0s with 1 indicating aln column in %s maps to non-gap RF column in %s saved to file %s.\n", alifile, esl_opt_GetString(go, "--map"), esl_opt_GetString(go, "--omap")); */
	  }
	  else printf("%s\n", msa1_to_msa2_mask);
	  free(msa1_to_msa2_mask);
	}

      /* --submap: if nec, map <msafile> to a subset of it's own columns in <f> (from --map <f>) */
      if(esl_opt_GetString(go, "--submap") != NULL)
	{
	  if((status = map_sub_msas(go, errbuf, msa, othermsa, &msa1_to_msa2_mask)) != eslOK) goto ERROR;
	  if(esl_opt_GetString(go, "--omap") != NULL) { 
	    if ((omapfp = fopen(esl_opt_GetString(go, "--omap"), "w")) == NULL) 
	      esl_fatal("Failed to open --omap output file %s\n", esl_opt_GetString(go, "--omap"));
	    fprintf(omapfp, "%s\n", msa1_to_msa2_mask);
	    fclose(omapfp);
	    /*printf("# Mask of 1/0s with 1 indicating aln column in %s maps to non-gap RF column in %s saved to file %s.\n", alifile, esl_opt_GetString(go, "--map"), esl_opt_GetString(go, "--omap")); */
	  }
	  else printf("%s\n", msa1_to_msa2_mask);
	  free(msa1_to_msa2_mask);
	}

      /* impose consensus structure to get individual secondary structures, if nec */
      if(esl_opt_GetBoolean(go, "--sindi"))
	{
	  if((status = individualize_consensus(go, errbuf, msa) != eslOK)) goto ERROR;
	  write_ali = TRUE;
	}

      /* handle the --tree option, if enabled */
      if( esl_opt_IsOn(go, "--tree"))
	{
	  if ((treefp = fopen(esl_opt_GetString(go, "--tree"), "w")) == NULL) 
	    esl_fatal("Failed to open --tree output file %s\n", esl_opt_GetString(go, "--tree"));

	  ESL_TREE    *T = NULL;/* the tree, created by Single-Linkage Clustering */
	  ESL_DMATRIX *D = NULL;/* the distance matrix */
	  
	  /* Create distance matrix and infer tree by single linkage clustering */
	  esl_dst_XDiffMx(msa->abc, msa->ax, msa->nseq, &D);
	  esl_tree_SingleLinkage(D, &T);
	  esl_tree_SetTaxaParents(T);
	  esl_tree_SetTaxonlabels(T, msa->sqname);

	  esl_tree_WriteNewick(treefp, T); 
	  fclose(treefp);
	  /*printf("# Tree saved in Newick format to file %s.\n", esl_opt_GetString(go, "--tree")); */

	  esl_tree_Validate(T, NULL);

	  /* Get new order for seqs in the MSA based on the tree */
	  int *order;
	  if((status = get_tree_order(T, errbuf, &order)) != eslOK) goto ERROR;
	  /*for(i = 0; i < msa->nseq; i++) {
	    printf("new MSA idx: %3d | orig MSA idx: %3d\n", i, order[i]);
	    }*/
	  esl_tree_Destroy(T);
	  esl_dmatrix_Destroy(D);
	  if((status = reorder_msa(msa, order, errbuf)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  free(order);
	}	  

      /* --xmask option: expand the alignment to fit lanemask in xmask <f>, number of TOTAL msa
       * columns must equal number of 1s in <f>.
       */
      if(xmask != NULL) { 
	ESL_MSA *newmsa;
	if((status = expand_msa2mask(errbuf, msa, xmask, &newmsa)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  msa = newmsa;
      }

      /* handle the --iinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--iinfo")) {
	if ((iinfofp = fopen(esl_opt_GetString(go, "--iinfo"), "w")) == NULL) 
	  esl_fatal("Failed to open --iinfo output file %s\n", esl_opt_GetString(go, "--iinfo"));
	if((status = dump_insert_info(iinfofp, msa, errbuf) != eslOK)) goto ERROR;
	/*printf("# Insert information saved to file %s.\n", esl_opt_GetString(go, "--iinfo")); */
	fclose(iinfofp);
      }

      /* handle the --iplot option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--iplot")) {
	if ((iplotfp = fopen(esl_opt_GetString(go, "--iplot"), "w")) == NULL) 
	  esl_fatal("Failed to open --iplot output file %s\n", esl_opt_GetString(go, "--iplot"));
	if((status = plot_inserts(iplotfp, msa, esl_opt_GetBoolean(go, "--ilog"), errbuf) != eslOK)) goto ERROR;
	fclose(iplotfp);
      }

      /* handle the --icinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--icinfo")) {
	if ((icinfofp = fopen(esl_opt_GetString(go, "--icinfo"), "w")) == NULL) 
	  esl_fatal("Failed to open --icinfo output file %s\n", esl_opt_GetString(go, "--iplot"));
	if((status = dump_infocontent(icinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(icinfofp);
      }

      /* handle the --gplot option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--gplot")) {
	if ((gplotfp = fopen(esl_opt_GetString(go, "--gplot"), "w")) == NULL) 
	  esl_fatal("Failed to open --gplot output file %s\n", esl_opt_GetString(go, "--gplot"));
	if((status = plot_gaps(gplotfp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(gplotfp);
      }

      /* handle the --rinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--rinfo")) {
	if ((rinfofp = fopen(esl_opt_GetString(go, "--rinfo"), "w")) == NULL) 
	  esl_fatal("Failed to open --rinfo output file %s\n", esl_opt_GetString(go, "--rinfo"));
	if((status = dump_residue_info(rinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(rinfofp);
      }

      /* handle the --cresinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--cresinfo")) {
	if ((cresinfofp = fopen(esl_opt_GetString(go, "--cresinfo"), "w")) == NULL) 
	  esl_fatal("Failed to open --cresinfo output file %s\n", esl_opt_GetString(go, "--cresinfo"));
	if((status = dump_cres_info(cresinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(cresinfofp);
      }

      /* handle the --dinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--dinfo")) {
	if ((dinfofp = fopen(esl_opt_GetString(go, "--dinfo"), "w")) == NULL) 
	  esl_fatal("Failed to open --dinfo output file %s\n", esl_opt_GetString(go, "--dinfo"));
	if((status = dump_delete_info(dinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(dinfofp);
      }

      /* handle the --num-rf and --num-all options, if enabled, do this after all MSA has been manipulated due to other options */
      if( esl_opt_IsOn(go, "--num-rf")) { 
	if((status = number_columns(msa, FALSE, errbuf) != eslOK)) goto ERROR;
	write_ali = TRUE;
      }
      if( esl_opt_IsOn(go, "--num-all")) { 
	if((status = number_columns(msa, TRUE, errbuf) != eslOK)) goto ERROR;
	write_ali = TRUE;
      }

      /* handle the -M option, if enabled */
      if( esl_opt_IsOn(go, "-M")) { 
	if((status = minorize_msa(go, msa, errbuf, ofp, esl_opt_GetString(go, "-M")) != eslOK)) goto ERROR;
	esl_msa_Destroy(msa);
	goto END;
      }

      /* handle the --c* options, if enabled */
      int do_id_cluster     = ((esl_opt_IsOn(go, "--cn-id"))  || (esl_opt_IsOn(go, "--cs-id"))  || (esl_opt_IsOn(go, "--cx-id"))) ? TRUE : FALSE;
      int do_insert_cluster = ((esl_opt_IsOn(go, "--cn-ins")) || (esl_opt_IsOn(go, "--cs-ins")) || (esl_opt_IsOn(go, "--cx-ins")))? TRUE : FALSE;
      int do_ctarget_nc, do_ctarget_nsize, do_cmindiff, nmsa, m, nc, nsize, xsize, nmin;
      float mindiff;

      if(do_id_cluster || do_insert_cluster) { 
	if(msa->rf == NULL) esl_fatal("--c* options require #=GC RF annotation marking consensus columns.");
	ESL_DMATRIX *D = NULL;/* the distance matrix */
	ESL_MSA     **cmsa;
	if(do_id_cluster) { 
	  /* create distance matrix and infer tree by single linkage clustering */
	  /* first, remove all non-consensus columns */
	  ESL_MSA *rfmsa;
	  rfmsa = esl_msa_Clone(msa);
	  if((status = keep_or_remove_rf_gaps(go, errbuf, rfmsa, TRUE, FALSE)) != eslOK) goto ERROR;
	  dst_nongap_XDiffMx(rfmsa->abc, rfmsa->ax, rfmsa->nseq, &D);
	  esl_msa_Destroy(rfmsa);
	  rfmsa = NULL;
	  do_ctarget_nc    = esl_opt_IsOn(go, "--cn-id");
	  do_ctarget_nsize = esl_opt_IsOn(go, "--cs-id");
	  do_cmindiff      = esl_opt_IsOn(go, "--cx-id");
	  nc               = esl_opt_IsOn(go, "--cn-id") ? esl_opt_GetInteger(go, "--cn-id")   : 0;
	  nsize            = esl_opt_IsOn(go, "--cs-id") ? esl_opt_GetInteger(go, "--cs-id")   : 0;
	  mindiff          = esl_opt_IsOn(go, "--cx-id") ? 1. - esl_opt_GetReal(go, "--cx-id") : 0; 
	}
	else { /* do_insert_cluster, create insert distance matrix and infer tree by SLC */ 
	  if((status = insert_x_diffmx(go, errbuf, msa, TRUE, TRUE, &D)) != eslOK) goto ERROR;
	  do_ctarget_nc    = esl_opt_IsOn(go, "--cn-ins");
	  do_ctarget_nsize = esl_opt_IsOn(go, "--cs-ins");
	  do_cmindiff      = esl_opt_IsOn(go, "--cx-ins");
	  nc               = esl_opt_IsOn(go, "--cn-ins") ? esl_opt_GetInteger(go, "--cn-ins")   : 0;
	  nsize            = esl_opt_IsOn(go, "--cs-ins") ? esl_opt_GetInteger(go, "--cs-ins")   : 0;
	  mindiff          = esl_opt_IsOn(go, "--cx-ins") ? 1. - esl_opt_GetReal(go, "--cx-ins") : 0;
	}
	/* print out the id matrix if nec */
	if( esl_opt_IsOn(go, "--c-mx")) { 
	  FILE *mxfp;
	  int i, j;
	  if ((mxfp = fopen(esl_opt_GetString(go, "--c-mx"), "w")) == NULL) esl_fatal("Failed to open --c-mx output file %s\n", esl_opt_GetString(go, "--c-mx"));
	  for(i = 0; i < msa->nseq; i++) { 
	    for(j = 0; j < msa->nseq; j++) { 
	      fprintf(mxfp, "%5d  %5d  %-30s  %-30s  %.5f\n", i, j, msa->sqname[i], msa->sqname[j], 1. - D->mx[i][j]);
	    }
	  }	  
	  fclose(mxfp);
	}
	  
	if((status = MSADivide(msa, D, do_cmindiff, do_ctarget_nc, do_ctarget_nsize, mindiff, nc, nsize, &nmsa, &cmsa, &xsize, errbuf)) != eslOK) goto ERROR;
	esl_msa_Destroy(msa); 
	msa = NULL;
	nmin = esl_opt_IsOn(go, "--c-nmin") ? esl_opt_GetInteger(go, "--c-nmin") : 1;
	for(m = 0; m < nmsa; m++) { 
	  if(cmsa[m]->nseq >= nmin) { 
	    status = esl_msa_Write(ofp, cmsa[m], (esl_opt_GetBoolean(go, "-1") ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM));
	    if      (status == eslEMEM) esl_fatal("Memory error when outputting alignment\n");
	    else if (status != eslOK)   esl_fatal("Writing alignment file failed with error %d\n", status);
	  }
	  esl_msa_Destroy(cmsa[m]);
	}
	write_ali = FALSE;
	free(cmsa);
      }
      else if ( esl_opt_IsOn(go, "--c-mx")) esl_fatal("--c-mx option requires at least one of: --cn-id, --cs-id, --cx-id, --cn-ins, --cs-ins, --cx-ins"); 

      /* remove GC annotation, if nec */
      if( esl_opt_IsOn(go, "--rm-gc")) {
	if((status = remove_gc_markup(msa, errbuf, esl_opt_GetString(go, "--rm-gc")) != eslOK)) goto ERROR;
	write_ali = TRUE;
      }

      /* write out list of sequences, if nec */
      if( esl_opt_IsOn(go, "--list")) {
	if ((listfp = fopen(esl_opt_GetString(go, "--list"), "w")) == NULL) 
	  esl_fatal("Failed to open --list output file %s\n", esl_opt_GetString(go, "--list"));
	int i;
	for(i = 0; i < msa->nseq; i++) fprintf(listfp, "%s\n", msa->sqname[i]);
	fclose(listfp);
     } 

      /* write out alignment, if nec */
      if((write_ali || esl_opt_GetBoolean(go, "-1")) && msa != NULL) {
	status = esl_msa_Write(ofp, msa, (esl_opt_GetBoolean(go, "-1") ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM));
	if      (status == eslEMEM) esl_fatal("Memory error when outputting alignment\n");
	else if (status != eslOK)   esl_fatal("Writing alignment file failed with error %d\n", status);
      }

      /* if nec, print #=GC RF annotation as a 1/0 mask (single line) to a file */
      if(esl_opt_GetString(go, "--omask") != NULL)
	{
	  if ((omaskfp = fopen(esl_opt_GetString(go, "--omask"), "w")) == NULL) 
	    esl_fatal("Failed to open --omask output file %s\n", esl_opt_GetString(go, "--omask"));
	  if((status = output_rf_as_mask(omaskfp, errbuf, msa)) != eslOK) goto ERROR;
	}
      if(msa != NULL)      esl_msa_Destroy(msa);
      if(othermsa != NULL) esl_msa_Destroy(othermsa);
    }

  /* If an msa read failed, we drop out to here with an informative status code. 
   */
  if      (status == eslEFORMAT) 
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)
    esl_fatal("No alignments found in file %s\n", alifile);

  if(esl_opt_GetString(go, "--omask") != NULL) fclose(omaskfp);
  
  /* Cleanup, normal return
   */
 END:
  if(otherafp != NULL) esl_msafile_Close(otherafp);
  if((esl_opt_GetString(go, "--morph") != NULL) && othermsa != NULL) esl_msa_Destroy(othermsa);

  if(esl_opt_GetString(go, "-o") != NULL) { fclose(ofp); }
  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  
  return 0;
  
 ERROR:
  if(afp != NULL) esl_msafile_Close(afp);
  if(go  != NULL) esl_getopts_Destroy(go);
  if(msa != NULL) esl_msa_Destroy(msa);

  esl_fatal(errbuf);
  return 1; /* never reached */
}


/* keep_or_remove_rf_gaps
 *                   
 * Given an MSA with #=GC RF markup, either remove or keep
 * all non-gap RF columns.
 */
static int
keep_or_remove_rf_gaps(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int keep_flag, int remove_flag)
{
  int     status;
  int    *useme;
  int64_t apos;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment.");
  if(keep_flag == TRUE && remove_flag == TRUE) ESL_XFAIL(eslEINVAL, errbuf, "in keep_or_remove_rf_gaps, keep_flag and remove_flag both TRUE.");
  if(keep_flag == FALSE && remove_flag == FALSE) ESL_XFAIL(eslEINVAL, errbuf, "in keep_or_remove_rf_gaps, keep_flag and remove_flag both FALSE.");

  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  if(keep_flag)
  {
    for(apos = 0; apos < msa->alen; apos++)    
      useme[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos]) ? FALSE : TRUE);
  }
  else if(remove_flag)
  {
    for(apos = 0; apos < msa->alen; apos++)    
      useme[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos]) ? TRUE : FALSE);
  }
  else ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "In keep_or_remove_rf_gaps, but neither -r nor -k enabled.");
  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  return eslEMEM;
}

/* keep_contiguous_column_block
 *                   
 * Keep only columns in range --start-all..--end-all, or --start-rf..--end-rf 
 */
static int
keep_contiguous_column_block(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int     status;
  int    *useme;
  int64_t apos;
  int     rf_mode;
  int     all_mode;
  int    *c2a_map = NULL;       /* msa map of consensus columns (non-gap RF residues) to alignment columns */
  int    clen;
  int    astart, aend;

  rf_mode  = ( esl_opt_IsOn(go, "--start-rf")  && esl_opt_IsOn(go, "--end-rf"))  ? TRUE : FALSE;
  all_mode = ( esl_opt_IsOn(go, "--start-all") && esl_opt_IsOn(go, "--end-all")) ? TRUE : FALSE;
  if((!rf_mode) && (!all_mode)) ESL_XFAIL(eslEINVAL, errbuf, "Entered keep_contiguous_column_block, but neither (--start-rf & --end-rf) nor (--start-all & --end-all) combination invoked.");
  
  /* contract check */
  if(rf_mode && msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "--start-rf and --end-rf required #=GC RF markup in alignment, but none exists.");

  if(rf_mode) { 
    if((status = map_cpos_to_apos(msa, &c2a_map, &clen))   != eslOK) goto ERROR;
    if(esl_opt_GetInteger(go, "--start-rf") < 1)    ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-rf must be > 1.");
    if(esl_opt_GetInteger(go, "--end-rf")   > clen) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --end-rf must be <= %d (which is the number of non-gap RF columns in the MSA).", clen);
    astart = c2a_map[esl_opt_GetInteger(go, "--start-rf")];
    aend   = c2a_map[esl_opt_GetInteger(go, "--end-rf")];
    if(astart > aend) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-rf <n> must be lower than <n> from --end-rf.");
  }
  else { 
    if(esl_opt_GetInteger(go, "--start-all") < 1)         ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-all must be > 1.");
    if(esl_opt_GetInteger(go, "--end-all")   > msa->alen) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --end-all must be <= %" PRId64 " (which is the number of columns in the MSA).", msa->alen);
    astart = esl_opt_GetInteger(go, "--start-all");
    aend   = esl_opt_GetInteger(go, "--end-all");
    if(astart > aend) ESL_XFAIL(eslEINVAL, errbuf, "<n> from --start-all <n> must be lower than <n> from --end-all.");
  }

  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  esl_vec_ISet(useme, msa->alen, FALSE);
  for(apos = astart-1; apos < aend; apos++) useme[apos] = TRUE;
  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  if(c2a_map != NULL) free(c2a_map);
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  return eslEMEM;
}

/* write_rf_gapthresh
 *                   
 * Given an MSA write/rewrite RF based on fraction
 * of gaps in each column. If fraction > gapthresh RF is an 'x',
 * otherwise it's a '.' (gap).
 */
static int
write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, float gapthresh)
{
  int      status;
  int64_t  apos;
  int64_t  gaps;
  int      i;

  if(msa->rf == NULL) ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (gaps = 0, i = 0; i < msa->nseq; i++)
	if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) gaps++;
      msa->rf[(apos-1)] = ((double) gaps / (double) msa->nseq > gapthresh) ? '.' : 'x';
    }
  msa->rf[msa->alen] = '\0';
  return eslOK;
 ERROR:
  return status;
}

/* write_rf_given_alen
 *                   
 * Given an MSA and a char string of 1s and 0s (a lanemask) of length
 * msa->alen, write/rewrite RF based as 'x' (non-gap) for 1, '.' (gap) for 0.
 */
static int
write_rf_given_alen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *amask)
{
  int      status;
  int64_t  apos;
  int64_t  mask_len;

  /* contract check, rfgiven_mask must be exact length of msa */
  if(amask == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask is NULL in write_rf_given, this shouldn't happen.\n");
  mask_len = strlen(amask);
  if(mask_len != msa->alen) 
    ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask length: %" PRId64 " is not equal to the MSA length (%" PRId64 ")\n", 
	     mask_len, msa->alen); 

  if(msa->rf == NULL) ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));

  for (apos = 1; apos <= msa->alen; apos++) {
      if     (amask[(apos-1)] == '0') msa->rf[(apos-1)] = '.';
      else if(amask[(apos-1)] == '1') msa->rf[(apos-1)] = 'x';
      else    ESL_FAIL(eslEINVAL, errbuf, "--mask-all mask char number %" PRId64 " is not a 1 nor a 0, but a %c\n", apos, amask[(apos-1)]);
  }
  msa->rf[msa->alen] = '\0';
  return eslOK;
 ERROR:
  return status;
}

/* write_rf_given_rflen
 *
 * Given an MSA and a char string of 1s and 0s (a lanemask) that is
 * the same length as the non-gap RF annotation in msa, rewrite msa
 * RF based as 'x' (non-gap) for 1, '.' (gap) for 0. 1s indicate which
 * non-gap RF columns to keep as 'x', and 0s indicate which non-gap
 * RF columns to make gaps '.'.
 */
static int
write_rf_given_rflen(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, char *rfmask)
{
  int64_t  apos, cpos;
  int64_t  mask_len;

  /* contract check, mask must be exact length of msa */
  if(rfmask == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-rf mask is NULL in write_rf_given, this shouldn't happen.\n");
  if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mask-rf mask requires RF annotation in MSA (try -g)\n");
  mask_len = strlen(rfmask);

  cpos = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) {
      cpos++;
      if     (rfmask[(cpos-1)] == '0') msa->rf[(apos-1)] = '.';
      else if(rfmask[(cpos-1)] == '1') msa->rf[(apos-1)] = 'x';
    }
    else msa->rf[(apos-1)] = '.'; 
  }
  msa->rf[msa->alen] = '\0';
  return eslOK;
}


/* write_rf_given_useme
 *
 * Given an MSA and a integer array <useme> of size msa->alen, set
 * msa->rf column [0..alen-1] i as 'x' if useme[i] == TRUE, and as
 * '.' if useme[i] == FALSE.
 */
static int
write_rf_given_useme(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int *useme)
{
  int     status;
  int64_t apos;

  if(msa->rf == NULL) ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  for (apos = 0; apos < msa->alen; apos++) msa->rf[apos] = useme[apos] ? 'x' : '.';
  msa->rf[msa->alen] = '\0';
  return eslOK;
 ERROR:
  return status;
}

/* individualize_consensus
 *                   
 * Given an MSA with a consensus structure impose it to create
 * individual secondary structures. Simple rule, for consensus
 * bp i,j if seq positions i and j are both non-gaps seq i,j are 
 * paired, if >= 1 is a gap, they're not paired.
 */
static int
individualize_consensus(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int     status;
  int64_t apos;
  int  i;
  int *cct;		/* 0..alen-1 base pair partners array for consensus        */
  int *ct;		/* 0..alen-1 base pair partners array for current sequence */
  char *ss;             /* individual secondary structure we've built              */
  char *ss_cons_nopseudo; /* no-pseudoknot version of consensus structure */

  if(msa->ss_cons == NULL)                                ESL_FAIL(eslEINVAL, errbuf, "--sindi requires MSA to have consensus structure annotation.\n");
  if(! (msa->flags & eslMSA_DIGITAL))                     ESL_FAIL(eslEINVAL, errbuf, "individualize_consensus() MSA is not digitized.\n");
    
  ESL_ALLOC(cct, sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ct,  sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ss,  sizeof(char) * (msa->alen+1));
  ESL_ALLOC(ss_cons_nopseudo, sizeof(char) * (msa->alen+1));

  esl_wuss_nopseudo(msa->ss_cons, ss_cons_nopseudo);
  if (esl_wuss2ct(ss_cons_nopseudo, msa->alen, cct) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");

  /* go through each position of each sequence, 
     if it's a gap and it is part of a base pair, remove that base pair */
  for (i = 0; i < msa->nseq; i++)
    {
      esl_vec_ICopy(cct, (msa->alen+1), ct);
      for (apos = 1; apos <= msa->alen; apos++)
	if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos]))
	  { 
	    if (ct[apos] != 0)  ct[ct[apos]] = 0;
	    ct[apos] = 0;
	  }
      /* convert to WUSS SS string and append to MSA */
      if (esl_ct2wuss(ct, msa->alen, ss) != eslOK) ESL_FAIL(status, errbuf, "Unexpected error converting de-knotted bp ct array to wuss notation.");
      esl_msa_AppendGR(msa, "SS", i, ss);
    }
  free(cct);
  free(ct);
  free(ss);
  free(ss_cons_nopseudo);
  return eslOK;
 ERROR:
  return status;
}

/* merge_msa
 *                   
 * Use the RF line as denoting consensus columns to merge
 * msa1 and msa2. msa1 is rewritten with merged msa. 
 * Important: msa1 will only contain sequence data from msa2.
 */
static int
merge_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_merged_msa)
{
  int status;
  int *agaps1 = NULL;
  int *agaps2 = NULL;
  int *c2a_map1 = NULL;       /* msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *c2a_map2 = NULL;       /* msa2 map of consensus columns (non-gap RF residues) to alignment columns */


  int *new_c2a_map1 = NULL;   /* merged msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *new_c2a_map2 = NULL;   /* merged msa2 map of consensus columns (non-gap RF residues) to alignment columns */
  int *aadd1 = NULL;          /* [1..apos..msa1->alen] number of columns to add after column apos to merge msa1 */
  int *aadd2 = NULL;          /* [1..apos..msa2->alen] number of columns to add after column apos to merge msa2 */
  int apos1, apos2;           /* counters over alignment positions of msa1, msa2 */
  int cpos = 0;               /* counter over consensus positions */
  int clen, clen2, new_clen1, new_clen2; /* consensus lengths */
  int tmp_ngaps;              /* temp var for number of gaps */
  int cur_apos1, nxt_apos1, cur_apos2, nxt_apos2; /* impt alignment positions */
  int astart2;                /* impt alignment positions */
  int ngaps1, ngaps2;         /* [1..apos..msa->alen] number of gaps in column apos of msa1, msa2 */
  int *msa2_cols_to_keep   = NULL; /* temp array for picking columns to keep in msa2 */
  int nadd1, nadd2;           /* temp vars, number of columns to keep, add */
  int radd = 0;               /* number of residues in msa2 columns corresponding to all 100% gap columns added to msa1 */
  int i, ip;                  /* sequence index counters */
  int x;                      /* general counter */
  int orig_msa1_nseq;         /* number of sequences in original msa1 */

  ESL_MSA *new_msa1;
  ESL_MSA *new_msa2;

  /* contract check */
  if(msa1->abc->type  != msa2->abc->type) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have same alphabet.");
  if(msa1->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have RF annotation.");
  if(msa2->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have RF annotation.");
  
  /* Determine number of gaps in each column of msa1 and msa2 */
  if((status = get_gaps_per_column(msa1, &agaps1)) != eslOK) goto ERROR;
  if((status = get_gaps_per_column(msa2, &agaps2)) != eslOK) goto ERROR;

  /* Map consensus columns to alignment positions */
  if((status = map_cpos_to_apos(msa1, &c2a_map1, &clen))   != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(msa2, &c2a_map2, &clen2))  != eslOK) goto ERROR;
  if(clen != clen2)
    ESL_XFAIL(eslEINVAL, errbuf, "With --merge both MSAs must have same consensus (non-gap RF) length.");
  
  /* Fill 'aadd1' and 'aadd2' arrays, these are the number columns of 100% gaps
   * we have to add after each 'msa1' and 'msa2' column respectively to make
   * the alignments the same size with identical non-gap RF lines.
   * Identical non-gap RF lines means for each aligned position i:
   *   isgap(msa1->rf[i]) == isgap(msa2->rf[i])
   */
  ESL_ALLOC(aadd1,  sizeof(int) * (msa1->alen+1));
  ESL_ALLOC(aadd2,  sizeof(int) * (msa2->alen+1));
  esl_vec_ISet(aadd1,  msa1->alen+1, 0);
  esl_vec_ISet(aadd2,  msa1->alen+1, 0);
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(cpos > 0) {
	cur_apos1 = c2a_map1[cpos];
	cur_apos2 = c2a_map2[cpos];
      }
      else cur_apos1 = cur_apos2 = 1;
      if(cpos < clen) {
	nxt_apos1 = c2a_map1[(cpos+1)];
	nxt_apos2 = c2a_map2[(cpos+1)];
      }
      else {
	nxt_apos1 = msa1->alen + 1;
	nxt_apos2 = msa2->alen + 1;
      }
      ngaps1 = nxt_apos1 - cur_apos1 - 1;
      ngaps2 = nxt_apos2 - cur_apos2 - 1;

      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* we don't have to add any columns to either msa (okay if 0) */
	{
	  /* do nothing */
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns to msa1 */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tmsa1 add     %4d all gap columns\n", (ngaps2-ngaps1)); 
	  nadd1 = ngaps2 - ngaps1;
	  /* determine where to put the gaps */
	  if(nxt_apos1 == (cur_apos1 + 1)) /* no choice, we have to put 100% gaps after cur_apos1 */
	    {
	      if(cpos == 0) aadd1[0] += nadd1;
	      else aadd1[c2a_map1[cpos]] += nadd1;
	    }
	  else
	    {
	      if(cpos == 0) 
		{ apos1 = astart2 = 0; }
	      else 
		{ 
		  apos1 = c2a_map1[cpos] + 1;
		  astart2 = cur_apos2+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps2, astart2, nxt_apos2-1, nadd1, &(msa2_cols_to_keep));
	      radd += (msa2->nseq * nadd1) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd1) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd1[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error!");
	      free(msa2_cols_to_keep);
	    }
	}
      else if(ngaps1 >  ngaps2) /* we need to add some new 100% gap columns to msa 2 */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tmsa2 add     %4d all gap columns\n", (ngaps1 - ngaps2));
	  nadd2 = ngaps1 - ngaps2;
	  /* determine where to put the gaps */
	  if(nxt_apos2 == (cur_apos2 + 1)) /* no choice, we have to put 100% gaps after cur_apos2 */
	    {
	      if(cpos == 0) aadd2[0] += nadd2;
	      else aadd2[c2a_map2[cpos]] += nadd2;
	    }
	  /*else
	    {
	      if(cpos == 0) 
		{ apos2 = astart1 = 0; }
	      else 
		{ 
		  apos2 = c2a_map2[cpos] + 1;
		  astart1 = cur_apos1+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps1, astart1, nxt_apos1-1, nadd, &(msa1_cols_to_keep));
	      radd += (msa2->nseq * nadd) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd1[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error!");
	      free(msa2_cols_to_keep);
	      }*/


	  /*hereherehere*/
	  /*if(cpos == 0) astart1 = 0;
	  else astart1 = cur_apos1+1;
	  if(ngaps2 == 0)
	    {
	      for(apos1 = astart1; apos1 < nxt_apos1; apos1++) akeep[apos1] = FALSE;
	      }*/
	  //else /* determine if it's likely flush left, or flush right first */
	  /*{
	      if(is_flush_left(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;         apos1 < (astart1 + nkeep); apos1++) 
		    akeep[apos1] = TRUE;
		  for(apos1 = (astart1 + nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = FALSE;
		}		  
	      else if(is_flush_right(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;             apos1 < (nxt_apos1 - nkeep); apos1++) 
		    akeep[apos1] = FALSE;
		  for(apos1 = (nxt_apos1 - nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = TRUE;
		    }*/
	  //else /* not flush left or flush right, pick least gappy columns to keep */
	  /*{
		  pick_gappiest_columns(agaps1, astart1, (nxt_apos1-1), (ngaps1 - nkeep), &(msa1_cols_to_remove));
		  for(apos1 = astart1; apos1 < nxt_apos1; apos1++) 
		    akeep[apos1] = (msa1_cols_to_remove[apos1 - astart1] == TRUE) ? FALSE : TRUE; 
		  free(msa1_cols_to_remove);
		}		
		}*/
	}
    }

  nadd1 = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      nadd1 += aadd1[apos1];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos1, aadd1[apos1]); }
    }
  nadd1 += aadd1[0];
  if(esl_opt_GetBoolean(go, "--verbose")) printf("Adding  %d columns to msa 1\n", nadd1);

  nadd2 = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa2 alignment column:\n"); }
  for(apos2 = 1; apos2 <= msa2->alen; apos2++)
    {
      nadd2 += aadd2[apos2];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos2, aadd2[apos2]); }
    }
  nadd2 += aadd2[0];
  if(esl_opt_GetBoolean(go, "--verbose")) printf("Adding  %d columns to msa 2\n", nadd2);

  /* add the 100% gap columns to msa1 and msa2 */
  status = add_gap_columns_to_msa(errbuf, msa1, aadd1, &new_msa1, TRUE);
  status = add_gap_columns_to_msa(errbuf, msa2, aadd2, &new_msa2, TRUE);

  /* Make new_c2a_map1 and new_c2a_map2, they should be identical */
  if((status = map_cpos_to_apos(new_msa1, &new_c2a_map1, &new_clen1))  != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(new_msa2, &new_c2a_map2, &new_clen2))  != eslOK) goto ERROR;
  if(new_clen1 != new_clen2) 
    ESL_XFAIL(eslEINVAL, errbuf, "Coding error, during alignment merge, after adding gaps, MSA lengths differ.");

  if(esl_opt_GetBoolean(go, "--verbose")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(new_c2a_map1[cpos] != new_c2a_map2[cpos]) 
	esl_fatal("Coding error. Alignments to merge do not have same consensus position map\n");
      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d %4d %4d\n", cpos, new_c2a_map1[cpos], new_c2a_map2[cpos]);
    }

  /* merge msa2 into msa1 */


  /* first make sure all the info that should be the same is the same */
  if(new_msa1->alen  != new_msa2->alen)  esl_fatal("Coding error. Alignments to merge do not have same lengths.\n");
  if(new_msa1->flags != new_msa2->flags) esl_fatal("Alignments to merge do not have flags (this *could* be worked around, implement it if you want).\n");
  if(new_msa1->abc->type != new_msa2->abc->type) esl_fatal("Alignments to merge do not have same alphabet.\n");
  for(x = 0; x < eslMSA_NCUTS; x++) 
    {
      if     ( new_msa1->cutset[x] && !new_msa2->cutset[x]) esl_fatal("Alignments to merge do not have same cutoff info.\n");
      else if(!new_msa1->cutset[x] &&  new_msa2->cutset[x]) esl_fatal("Alignments to merge do not have same cutoff info.\n");
      else if( new_msa1->cutset[x] &&  new_msa2->cutset[x])
	if(fabs(new_msa1->cutoff[x] - new_msa2->cutoff[x]) > 0.0001)
	  esl_fatal("Alignments to merge do not have same cutoff info.\n");
    }

  /* now merge new_msa1 and new_msa2, by expanding new_msa1, and swapping ptrs to data in new_msa2, 
   * new_msa1 becomes merged alignment
   * new_msa2 becomes pathetic shell of an alignment
   *
   * to expand a MSA, the alen must be 0 (flag for esl_msa_Expand()) I 
   * reset it to 0 here and then back again. This may be ill advised 
   */
  new_msa1->alen = 0;
  while(new_msa1->sqalloc < (new_msa1->nseq + new_msa2->nseq)) 
    esl_msa_Expand(new_msa1);
  new_msa1->alen = new_msa2->alen;
  orig_msa1_nseq = new_msa1->nseq;
  
  if((new_msa1->ss_cons == NULL && new_msa2->ss_cons != NULL) ||
     (new_msa1->ss_cons != NULL && new_msa2->ss_cons == NULL) ||
     ((new_msa1->ss_cons != NULL && new_msa2->ss_cons != NULL) && 
      (strcmp(new_msa1->ss_cons, new_msa2->ss_cons)   != 0))) esl_fatal("Alignments to merge do not have same consensus structure.\n");
  if((new_msa1->sa_cons == NULL && new_msa2->sa_cons != NULL) ||
     (new_msa1->sa_cons != NULL && new_msa2->sa_cons == NULL) ||
     ((new_msa1->sa_cons != NULL && new_msa2->sa_cons != NULL) && 
      (strcmp(new_msa1->sa_cons, new_msa2->sa_cons)   != 0))) esl_fatal("Alignments to merge do not have same consensus structure.\n");
  if((new_msa1->pp_cons == NULL && new_msa2->pp_cons != NULL) ||
     (new_msa1->pp_cons != NULL && new_msa2->pp_cons == NULL) ||
     ((new_msa1->pp_cons != NULL && new_msa2->pp_cons != NULL) && 
      (strcmp(new_msa1->pp_cons, new_msa2->pp_cons)   != 0))) esl_fatal("Alignments to merge do not have same consensus posteriors.\n");
  if((new_msa1->aseq == NULL && new_msa2->aseq != NULL) ||
     (new_msa1->aseq != NULL && new_msa2->aseq == NULL))  esl_fatal("Alignments to merge aseqs null/non-null mismatch.\n");
#ifdef eslAUGMENT_ALPHABET
  if((new_msa1->ax == NULL && new_msa2->ax != NULL) ||
     (new_msa1->ax != NULL && new_msa2->ax == NULL))  esl_fatal("Alignments to merge ax null/non-null mismatch.\n");
#endif /*eslAUGMENT_ALPHABET*/
  if((new_msa1->sqacc == NULL && new_msa2->sqacc != NULL) ||
     (new_msa1->sqacc != NULL && new_msa2->sqacc == NULL))  esl_fatal("Alignments to merge sqacc null/non-null mismatch.\n");
  if((new_msa1->sqdesc == NULL && new_msa2->sqdesc != NULL) ||
     (new_msa1->sqdesc != NULL && new_msa2->sqdesc == NULL))  esl_fatal("Alignments to merge sqdesc null/non-null mismatch.\n");
  if((new_msa1->ss == NULL && new_msa2->ss != NULL) ||
     (new_msa1->ss != NULL && new_msa2->ss == NULL))  esl_fatal("Alignments to merge ss null/non-null mismatch.\n");
  if((new_msa1->sa == NULL && new_msa2->sa != NULL) ||
     (new_msa1->sa != NULL && new_msa2->sa == NULL))  esl_fatal("Alignments to merge sa null/non-null mismatch.\n");
  if((new_msa1->pp == NULL && new_msa2->pp != NULL) ||
     (new_msa1->pp != NULL && new_msa2->pp == NULL))  esl_fatal("Alignments to merge pp null/non-null mismatch.\n");

     /* rf lines were already indirectly checked, they must be equal */

  for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
    {
      ip = i - orig_msa1_nseq;
      if(new_msa1->aseq != NULL) new_msa1->aseq[i]   = new_msa2->aseq[ip];
#ifdef eslAUGMENT_ALPHABET
      if(new_msa1->ax   != NULL) new_msa1->ax[i]     = new_msa2->ax[ip];
#endif /*eslAUGMENT_ALPHABET*/
      new_msa1->sqname[i] = new_msa2->sqname[ip];
      new_msa1->wgt[i]    = new_msa2->wgt[ip];
      new_msa1->nseq++;

      if(new_msa1->sqacc  != NULL) new_msa1->sqacc[i]  = new_msa2->sqacc[ip];
      if(new_msa1->sqdesc != NULL) new_msa1->sqdesc[i] = new_msa2->sqdesc[ip];
      if(new_msa1->ss     != NULL) new_msa1->ss[i]     = new_msa2->ss[ip];
      if(new_msa1->sa     != NULL) new_msa1->sa[i]     = new_msa2->sa[ip];
      if(new_msa1->pp     != NULL) new_msa1->pp[i]     = new_msa2->pp[ip];

      /* new_msa1->name,desc,acc,au untouched (is this unwise?) */

      if(new_msa1->sqlen  != NULL) new_msa1->sqlen[i]  = new_msa2->sqlen[ip];

      if(new_msa1->sslen != NULL) new_msa1->sslen[i]  = new_msa2->sslen[ip];
      if(new_msa1->salen != NULL) new_msa1->salen[i]  = new_msa2->salen[ip];
      if(new_msa1->pplen != NULL) new_msa1->pplen[i]  = new_msa2->pplen[ip];
      /* lastidx not touched, should be unimportant */
    }
  /* copy and free comments (no need to swap pointers thanks to convenient esl_msa_AddComment() function */
  for(x = 0; x < new_msa2->ncomment; x++) {
    esl_msa_AddComment(new_msa1, new_msa2->comment[x]);
    free(new_msa2->comment[x]);
  }
  /* copy and free GF markup */
  for(x = 0; x < new_msa2->ngf; x++) {
    esl_msa_AddGF(new_msa1, new_msa2->gf_tag[x], new_msa2->gf[x]);
    free(new_msa2->gf_tag[x]);
    free(new_msa2->gf[x]);
  }
  /* copy and free GS markup */
  for(x = 0; x < new_msa2->ngs; x++) {
    for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
      {
	ip = i - orig_msa1_nseq;
	esl_msa_AddGS(new_msa1, new_msa2->gs_tag[x], i, new_msa2->gs[x][ip]);
	free(new_msa2->gs[x][ip]);
      }
    free(new_msa2->gs_tag[x]);
  }
  /* don't touch GC (per column) annotation, is this unwise? */

  /* copy and free GR markup */
  for(x = 0; x < new_msa2->ngr; x++) {
    for(i = orig_msa1_nseq; i < (orig_msa1_nseq + new_msa2->nseq); i++)
      {
	ip = i - orig_msa1_nseq;
	esl_msa_AppendGR(new_msa1, new_msa2->gr_tag[x], i, new_msa2->gr[x][ip]);
	free(new_msa2->gr[x][ip]);
      }
    free(new_msa2->gr_tag[x]);
  }
  /* don't touch keyhashes or SSI offset, shouldn't be a problem since we're just
   * printing the alignment. 
   */

  *ret_merged_msa = new_msa1;
  free(new_msa2); /* the guts are still valid, being pointed to by new_msa1 */

  free(agaps1);
  free(agaps2);
  free(c2a_map1);
  free(c2a_map2);
  free(new_c2a_map1);
  free(aadd1);
  free(aadd2);
  return eslOK;

 ERROR:
  if(agaps1       != NULL) free(agaps1);
  if(agaps2       != NULL) free(agaps2);
  if(c2a_map1     != NULL) free(c2a_map1);
  if(c2a_map2     != NULL) free(c2a_map2);
  if(new_c2a_map1 != NULL) free(new_c2a_map1);
  if(aadd1        != NULL) free(aadd1);
  if(aadd2        != NULL) free(aadd2);
  return status;
}

/* morph_msa
 *                   
 * Use the RF line as denoting consensus columns to morph
 * msa1 into msa2's gap structure. This may require removing
 * some columns from msa1, and adding some 100% gap columns
 * to msa1.
 */
static int
morph_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **new_msa1)
{
  int status;
  int *agaps1 = NULL;
  int *agaps2 = NULL;
  int *c2a_map1 = NULL;       /* msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *c2a_map2 = NULL;       /* msa2 map of consensus columns (non-gap RF residues) to alignment columns */
  int *new_c2a_map1 = NULL;   /* morphed msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int *akeep = NULL;          /* [1..apos..msa1->alen] TRUE to keep column apos in morphed msa1, FALSE not to */
  int *aadd = NULL;           /* [1..apos..msa1->alen] number of columns to add after column apos to morphed msa1 */
  int apos1, apos2;           /* counters over alignment positions of msa1, msa2 */
  int cpos = 0;               /* counter over consensus positions */
  int clen, clen2, new_clen1; /* consensus lengths */
  int tmp_ngaps;              /* temp var for number of gaps */
  int cur_apos1, nxt_apos1, cur_apos2, nxt_apos2; /* impt alignment positions */
  int astart1, astart2;       /* impt alignment positions */
  int ngaps1, ngaps2;         /* [1..apos..msa->alen] number of gaps in column apos of msa1, msa2 */
  int *msa1_cols_to_remove = NULL; /* temp array for picking columns to remove from msa1 */
  int *msa2_cols_to_keep   = NULL; /* temp array for picking columns to keep in msa2 */
  int nkeep, nadd;            /* temp vars, number of columns to keep, add */
  int radd = 0;               /* number of residues in msa2 columns corresponding to all 100% gap columns added to msa1 */
  int delete = 0;             /* number of residues in msa1 we have to delete during morph */

  /* contract check */
  if(msa1->abc->type  != msa2->abc->type) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same alphabet.");
  if(msa1->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have RF annotation.");
  if(msa2->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have RF annotation.");
  
  /* Determine number of gaps in each column of msa1 and msa2 */
  if((status = get_gaps_per_column(msa1, &agaps1)) != eslOK) goto ERROR;
  if((status = get_gaps_per_column(msa2, &agaps2)) != eslOK) goto ERROR;

  /* Map consensus columns to alignment positions */
  if((status = map_cpos_to_apos(msa1, &c2a_map1, &clen))   != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(msa2, &c2a_map2, &clen2))  != eslOK) goto ERROR;
  if(clen != clen2)
    ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same consensus (non-gap RF) length.");
  
  /* Fill the 'akeep' array [1..msa1->alen] and 'aadd' array [0..msa1->alen] which 
   * tells us which columns in msa1 we'll keep, and how many columns of 100% gaps
   * to add after each msa1 column. 
   */
  ESL_ALLOC(akeep, sizeof(int) * (msa1->alen+1));
  ESL_ALLOC(aadd,  sizeof(int) * (msa1->alen+1));
  esl_vec_ISet(akeep, msa1->alen+1, FALSE);
  esl_vec_ISet(aadd,  msa1->alen+1, 0);
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(cpos > 0) {
	cur_apos1 = c2a_map1[cpos];
	cur_apos2 = c2a_map2[cpos];
      }
      else cur_apos1 = cur_apos2 = 1;
      if(cpos < clen) {
	nxt_apos1 = c2a_map1[(cpos+1)];
	nxt_apos2 = c2a_map2[(cpos+1)];
      }
      else {
	nxt_apos1 = msa1->alen + 1;
	nxt_apos2 = msa2->alen + 1;
      }
      akeep[cur_apos1] = TRUE; /* keep the consensus column */
      ngaps1 = nxt_apos1 - cur_apos1 - 1;
      ngaps2 = nxt_apos2 - cur_apos2 - 1;

      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* keep all columns in between (okay if 0) */
	{
	  for(apos1 = cur_apos1+1; apos1 < nxt_apos1; apos1++) akeep[apos1] = TRUE; 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tadd     %4d all gap columns\n", (ngaps2-ngaps1)); 
	  nadd = ngaps2 - ngaps1;
	  /* keep all the inserts we have in msa1 */
	  for(apos1 = cur_apos1+1; apos1 < nxt_apos1; apos1++) akeep[apos1] = TRUE;
	  if(nxt_apos1 == (cur_apos1 + 1)) /* no choice, we have to put 100% gaps after cur_apos1 */
	    {
	      if(cpos == 0) aadd[0] += nadd;
	      else aadd[c2a_map1[cpos]] += nadd;
	    }
	  else
	    {
	      if(cpos == 0) 
		{ apos1 = astart2 = 0; }
	      else 
		{ 
		  apos1 = c2a_map1[cpos] + 1;
		  astart2 = cur_apos2+1; 
		}
	      tmp_ngaps = pick_gappiest_columns(agaps2, astart2, nxt_apos2-1, nadd, &(msa2_cols_to_keep));
	      radd += (msa2->nseq * nadd) - tmp_ngaps;
	      if(esl_opt_GetBoolean(go, "--verbose")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
	      for(apos2 = astart2; apos2 < nxt_apos2; apos2++) 
		{
		  if(msa2_cols_to_keep[(apos2 - astart2)] == TRUE) aadd[apos1]++;
		  else apos1++;
		}
	      if(apos1 != nxt_apos1) 
		esl_fatal("Coding error 10.");
	      free(msa2_cols_to_keep);
	    }
	}
      else if(ngaps1 >  ngaps2) /* we need to delete some of our msa1 columns */
	{ 
	  if(esl_opt_GetBoolean(go, "--verbose")) printf("\tdelete  %4d/%4d    columns\n", (ngaps1 - ngaps2), (ngaps1));  
	  nkeep = ngaps2;
	  if(cpos == 0) astart1 = 0;
	  else astart1 = cur_apos1+1;
	  if(ngaps2 == 0)
	    {
	      for(apos1 = astart1; apos1 < nxt_apos1; apos1++) akeep[apos1] = FALSE;
	    }
	  else /* determine if it's likely flush left, or flush right first */
	    {
	      if(is_flush_left(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;         apos1 < (astart1 + nkeep); apos1++) 
		    akeep[apos1] = TRUE;
		  for(apos1 = (astart1 + nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = FALSE;
		}		  
	      else if(is_flush_right(agaps1, astart1, nxt_apos1-1))
		{
		  for(apos1 = astart1;             apos1 < (nxt_apos1 - nkeep); apos1++) 
		    akeep[apos1] = FALSE;
		  for(apos1 = (nxt_apos1 - nkeep); apos1 <  nxt_apos1;              apos1++)
		    akeep[apos1] = TRUE;
		}
	      else /* not flush left or flush right, pick least gappy columns to keep */
		{
		  pick_gappiest_columns(agaps1, astart1, (nxt_apos1-1), (ngaps1 - nkeep), &(msa1_cols_to_remove));
		  for(apos1 = astart1; apos1 < nxt_apos1; apos1++) 
		    akeep[apos1] = (msa1_cols_to_remove[apos1 - astart1] == TRUE) ? FALSE : TRUE; 
		  free(msa1_cols_to_remove);
		}		
	    }
	}
    }

  nadd = 0;
  nkeep = 0;
  if(esl_opt_GetBoolean(go, "--verbose")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      if(akeep[apos1]) nkeep++;
      else delete += (msa1->nseq - agaps1[apos1]);
      nadd += aadd[apos1];
      if(esl_opt_GetBoolean(go, "--verbose")) { printf("%5d %5d\n", apos1, aadd[apos1]); }
    }
  nadd += aadd[0];
  printf("\n\nKeeping %d columns, deleting %d residues.\n", nkeep, delete);
  printf("Adding  %d columns, which have %d total non-gaps in MSA2.\n", nadd, radd);

  /* Rewrite the msa->rf line so that we can call keep_or_remove_rf_gaps to remove
   * the columns we don't want. Then restore the rf line. Do this by making a new
   * #=GC ORIGRF markup line. This is a serious violation of easel conventions.
   */
  char *origrf;
  esl_strdup(msa1->rf, msa1->alen, &origrf);
  esl_msa_AppendGC(msa1, "ORIGRF", origrf);
  /* overwrite RF temporarily with 'x's for any column we're keeping, '.' for any we're losing */
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    msa1->rf[(apos1-1)] = akeep[apos1] == FALSE ? '.' : 'x';
  
  /* add the 100% gap columns */
  status = add_gap_columns_to_msa(errbuf, msa1, aadd, new_msa1, FALSE);

  /* remove unwanted columns */
  keep_or_remove_rf_gaps(go, errbuf, *new_msa1, TRUE, FALSE); 

  /* restore RF line */
  free((*new_msa1)->rf);
  esl_strdup((*new_msa1)->gc[(*new_msa1)->ngc-1], (*new_msa1)->alen, &((*new_msa1)->rf));
  free(origrf);
  free((*new_msa1)->gc_tag[((*new_msa1)->ngc-1)]);
  free((*new_msa1)->gc[((*new_msa1)->ngc-1)]);
  (*new_msa1)->ngc--;

  /* Make new new_c2a_map1, it should be identical to c2a_map2. */
  if((status = map_cpos_to_apos((*new_msa1), &new_c2a_map1, &new_clen1))  != eslOK) goto ERROR;
  if(new_clen1 != clen) 
    ESL_XFAIL(eslEINVAL, errbuf, "With --morph both MSAs must have same consensus (non-gap RF) length.");

  if(esl_opt_GetBoolean(go, "--verbose")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(c2a_map2[cpos] != new_c2a_map1[cpos]) 
	esl_fatal("Coding error. Morphed alignment does not have same consensus position map as %s\n", esl_opt_GetString(go, "--morph"));
      if(esl_opt_GetBoolean(go, "--verbose")) printf("%4d %4d %4d %4d\n", cpos, c2a_map2[cpos], new_c2a_map1[cpos], (c2a_map2[cpos] - new_c2a_map1[cpos]));
    }
  
  free(agaps1);
  free(agaps2);
  free(c2a_map1);
  free(c2a_map2);
  free(new_c2a_map1);
  free(akeep);
  free(aadd);
  return eslOK;

 ERROR:
  if(agaps1       != NULL) free(agaps1);
  if(agaps2       != NULL) free(agaps2);
  if(c2a_map1     != NULL) free(c2a_map1);
  if(c2a_map2     != NULL) free(c2a_map2);
  if(new_c2a_map1 != NULL) free(new_c2a_map1);
  if(akeep        != NULL) free(akeep);
  if(aadd         != NULL) free(aadd);
  return status;
}

/* add_gap_columns_to_msa
 *                   
 * Given an MSA and an array specifying a number
 * of all gap columns to add after each column,
 * add them. Reallocate all arrays as necessary.
 * if(do_treat_as_rf_gap) make new column a gap
 * in the RF line, else make it an 'x'.
 *
 * toadd is numbered 1..alen.
 */
static int
add_gap_columns_to_msa(char *errbuf, ESL_MSA *msa, int *toadd, ESL_MSA **ret_msa, int do_treat_as_rf_gap)
{
  int status;
  int i,j;
  int apos;
  int nnew = 0;
  ESL_ALPHABET *abc;
  char *newstr;
  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in add_gap_columns_to_msa(), msa must be digitized.");
  for(apos = 0; apos <= msa->alen; apos++)
    nnew += toadd[apos];

  /* Textize the alignment */
  abc = msa->abc;
  esl_msa_Textize(msa);

  ESL_MSA *newmsa;
  /*printf("msa->nseq: %d\n", msa->nseq);
    printf("msa->alen: %d\n", msa->alen);*/
  newmsa = esl_msa_Create(msa->nseq, (msa->alen+nnew));

  /* Copy and add gaps to all valid data that is [0..(alen-1)] or [1..alen] */ 
  if(msa->ss_cons != NULL) 
    {
      ESL_ALLOC(newmsa->ss_cons, sizeof(char) * (msa->alen+nnew+1));
      if((status = cp_and_add_gaps_to_aseq(newmsa->ss_cons, msa->ss_cons, msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
    }
  if(msa->sa_cons != NULL) 
    {
      ESL_ALLOC(newmsa->sa_cons, sizeof(char) * (msa->alen+nnew+1));
      if((status = cp_and_add_gaps_to_aseq(newmsa->sa_cons, msa->sa_cons, msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
    }
  if(msa->pp_cons != NULL) 
    {
      ESL_ALLOC(newmsa->pp_cons, sizeof(char) * (msa->alen+nnew+1));
      if((status = cp_and_add_gaps_to_aseq(newmsa->pp_cons, msa->pp_cons, msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
    }
  if(msa->rf != NULL)
    {
      ESL_ALLOC(newmsa->rf, sizeof(char) * (msa->alen+nnew+1));
      if(do_treat_as_rf_gap)
	{
	  if((status = cp_and_add_gaps_to_aseq(newmsa->rf,      msa->rf,      msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	}
      else if((status = cp_and_add_gaps_to_aseq(newmsa->rf,      msa->rf,      msa->alen, toadd, nnew, 'x') != eslOK)) goto ERROR;
    }

  if(msa->ss != NULL)
    {
      ESL_ALLOC(newmsa->ss, sizeof(char *) * msa->nseq);
      for(i = 0; i < msa->nseq; i++)
      {
	if(msa->ss[i] != NULL)
	  {
	    ESL_ALLOC(newmsa->ss[i], sizeof(char) * (msa->alen+nnew+1));
	    if((status = cp_and_add_gaps_to_aseq(newmsa->ss[i], msa->ss[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	  }
      }
    }

  if(msa->sa != NULL)
    {
      for(i = 0; i < msa->nseq; i++)
      {
	if(msa->sa[i] != NULL)
	  {
	    ESL_ALLOC(newmsa->sa[i], sizeof(char) * (msa->alen+nnew+1));
	    if((status = cp_and_add_gaps_to_aseq(newmsa->sa[i], msa->sa[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	  }
      }
    }  

  if(msa->pp != NULL)
    {
      for(i = 0; i < msa->nseq; i++)
      {
	if(msa->pp[i] != NULL)
	  {
	    ESL_ALLOC(newmsa->pp[i], sizeof(char) * (msa->alen+nnew+1));
	    if((status = cp_and_add_gaps_to_aseq(newmsa->pp[i], msa->pp[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	  }
      }
    }  

  if(msa->ncomment > 0)
    {
      for(j = 0; j < msa->ncomment; j++)
	{
	  if(msa->comment[j] != NULL) 
	    esl_msa_AddComment(newmsa, msa->comment[j]);
	}
    }

  if(msa->ngf > 0)
    {
      for(i = 0; i < msa->ngf; i++)
	if(msa->gf[i] != NULL) 
	    esl_msa_AddGF(newmsa, msa->gf_tag[i], msa->gf[i]);
    }

  if(msa->ngs > 0)
    {
      for(j = 0; j < msa->ngs; j++)
	{
	  for(i = 0; i < msa->nseq; i++)
	    if(msa->gs[j][i] != NULL) 
	      esl_msa_AddGS(newmsa, msa->gs_tag[j], i, msa->gs[j][i]);
	}
    }

  if(msa->ngc > 0)
    {
      for(i = 0; i < msa->ngc; i++)
	{
	  if(msa->gc[i] != NULL) 
	    {
	      ESL_ALLOC(newstr, sizeof(char) * (msa->alen+nnew+1));
	      if((status = cp_and_add_gaps_to_aseq(newstr, msa->gc[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
	      esl_msa_AppendGC(newmsa, msa->gc_tag[i], newstr);
	      free(newstr);
	    }
	}
    }

  if(msa->gr != NULL)
    {  
      for(j = 0; j < msa->ngr; j++)
	{
	  for(i = 0; i < msa->nseq; i++)
	    {
	      if(msa->gr[j][i] != NULL) 
		{
		  ESL_ALLOC(newstr, sizeof(char) * (msa->alen+nnew+1));
		  if((status = cp_and_add_gaps_to_aseq(newstr, msa->gr[j][i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
		  esl_msa_AppendGC(newmsa, msa->gc_tag[i], newstr);
		  free(newstr);
		}
	    }
	}
    }
    
  /* copy the aseqs, free as we go to save memory */
  for(i = 0; i < msa->nseq; i++)
    {
      esl_strdup(msa->sqname[i], -1, &(newmsa->sqname[i]));
      if((status = cp_and_add_gaps_to_aseq(newmsa->aseq[i], msa->aseq[i], msa->alen, toadd, nnew, '.') != eslOK)) goto ERROR;
      free(msa->aseq[i]);
      msa->aseq[i] = NULL;
    }    
  newmsa->abc = abc;
  esl_msa_Digitize(newmsa->abc, newmsa);
  esl_msa_Destroy(msa);
  *ret_msa = newmsa;
      
  return eslOK;
  
 ERROR:
  return status;
}

/*cp_and_add_gaps_to_aseq
 *                   
 * Given an aligned [0..alen-1] original text string,
 * add toadd[apos-1] gaps after each residue. 
 * new_aseq must be already allocated. 
 *
 * toadd is numbered 1..alen.
 */
static int cp_and_add_gaps_to_aseq(char *new_aseq, char *orig_aseq, int alen, int *toadd, int nnew, char gapchar)
{
  int orig_apos = 0;
  int new_apos  = 0;
  int i;

  for(i = 0; i < toadd[0]; i++)
    new_aseq[new_apos++] = gapchar;
  for(orig_apos = 0; orig_apos < alen; orig_apos++)
    {
      new_aseq[new_apos++] = orig_aseq[orig_apos];
      for(i = 0; i < toadd[(orig_apos+1)]; i++)
	new_aseq[new_apos++] = gapchar;
    }
  new_aseq[new_apos] = '\0';
  return eslOK;
}

/* is_flush_left
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend,
 * return TRUE if the residues in this interval appear to 
 * leftflushed inserts, FALSE otherwise
 */
static int is_flush_left(int *ngaps, int astart, int aend)
{
  if(astart == -1 || aend == -1) esl_fatal("is_flush_left invalid column positions.");
  
  int i;
  int gaps = ngaps[astart];
  for(i = astart+1; i <= aend; i++)
    {
      if(ngaps[i] < gaps) return FALSE;
      gaps = ngaps[i];
    }
  return TRUE;
}

/* is_flush_right
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend,
 * return TRUE if the residues in this interval appear to 
 * rightflushed inserts, FALSE otherwise
 */
static int is_flush_right(int *ngaps, int astart, int aend)
{
  if(astart == -1 || aend == -1) esl_fatal("is_flush_right invalid column positions.");
  
  int i;
  int gaps = ngaps[astart];
  for(i = astart+1; i <= aend; i++)
    {
      if(ngaps[i] > gaps) return FALSE;
      gaps = ngaps[i];
    }
  return TRUE;
}

/* pick_gappiest_columns
 *                   
 * Given an array with number of gaps in each column
 * of an alignment, and an interval of columns astart..aend.
 * Pick the npick gappiest columns and store that info
 * in ret_cols_to_pick.
 *
 * Returns total number of gaps the npick columns picked.
 */
static int pick_gappiest_columns(int *ngaps, int astart, int aend, int npick, int **ret_cols_to_pick)
{
  if(astart == -1 || aend == -1) esl_fatal("pick_gappiest_columns invalid column positions.");
  if((aend-astart+1) < npick)    esl_fatal("pick_gappiest_columns number to pick (%d) exceeds number of possibilities (%d).", npick, (aend-astart+1));

  int status;
  int i,c;
  int *tmp_ngaps;
  int *cols_to_pick;
  int topick;
  int total_gaps = 0;

  ESL_ALLOC(tmp_ngaps,    sizeof(int) * (aend-astart+1));
  ESL_ALLOC(cols_to_pick, sizeof(int) * (aend-astart+1));

  esl_vec_ISet(cols_to_pick, (aend-astart+1), FALSE);
  for(i = astart; i <= aend; i++)
    tmp_ngaps[(i-astart)] = ngaps[astart];
  for(c = 0; c < npick; c++)
    {
      topick               = esl_vec_IArgMax(tmp_ngaps, (aend-astart+1));
      cols_to_pick[topick] = TRUE;
      total_gaps          += tmp_ngaps[topick];
      tmp_ngaps[topick]    = -1;
    }
  free(tmp_ngaps);
  *ret_cols_to_pick = cols_to_pick;
  return total_gaps;

 ERROR:
  esl_fatal("Memory allocation error.");
  return -1;
}

/* get_gaps_per_column 
 *                   
 * Given an MSA, determine the number of gaps per
 * column, and return a newly allocated array with this
 * into in *ret_ngaps. 
 */
static int get_gaps_per_column(ESL_MSA *msa, int **ret_ngaps)
{
  int status;
  int i, apos;
  int *ngaps = NULL;
  /* contract check */
  if(! msa->flags & eslMSA_DIGITAL) { status = eslEINVAL; goto ERROR; }

  ESL_ALLOC(ngaps, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(ngaps, msa->alen+1, 0);
  for(i = 0; i < msa->nseq; i++) {
    for(apos = 1; apos <= msa->alen; apos++)
      ngaps[apos] += esl_abc_XIsGap(msa->abc, msa->ax[i][apos]);
  }
  *ret_ngaps = ngaps;
  return eslOK;

 ERROR:
  if(ngaps != NULL) free(ngaps);
  return status;
}

/* map_cpos_to_apos
 *                   
 * Given an MSA, determine the alignment position each
 * consensus position refers to. 
 */
static int map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen)
{
  int status;
  int clen = 0;
  int *c2a_map = NULL;
  int cpos = 0;
  int apos = 0;
  /* contract check */
  if(msa->rf == NULL) { status = eslEINVAL; goto ERROR; }

  /* count consensus columns */
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  /* build map */
  ESL_ALLOC(c2a_map, sizeof(int) * (clen+1));
  c2a_map[0] = -1;
  for(apos = 1; apos <= msa->alen; apos++) 
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) c2a_map[++cpos] = apos;

  *ret_c2a_map = c2a_map;
  *ret_clen    = clen;
  return eslOK;

 ERROR:
  if(c2a_map != NULL) free(c2a_map);
  return status;
}


/* read_sqfile
 *                   
 * Read all seqs in a sequence file and return them. Originally
 * written for --trim option.
 */
static int read_sqfile(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc, int nseq, ESL_SQ ***ret_sq)
{
  int status;
  ESL_SQ **sq; 
  int i;
  
  /* get seqs from sqfile */
  ESL_ALLOC(sq, sizeof(ESL_SQ *) * (nseq + 1)); /* +1 for the last guy we allocate but don't use */
  i = 0;
  sq[i] = esl_sq_CreateDigital(abc);
  while ((status = esl_sqio_Read(sqfp, sq[i])) == eslOK) { 
    i++;
    if(i > nseq) esl_fatal("With --trim, sequence file must have same number seqs as in <msafile>\n"); 
    sq[i] = esl_sq_CreateDigital(abc);
  }
  if (i != nseq) esl_fatal("With --trim, sequence file must have same number seqs as in <msafile>\n"); 
  /* status should be eslEOF on normal end; if it isn't, deal w/ error */
  esl_sq_Destroy(sq[i]); /* destroy final allocated but unused seq */

  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
					    sqfp->filename, sqfp->linenumber, sqfp->errbuf);     
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    status, sqfp->filename);
  esl_sqfile_Close(sqfp);
  *ret_sq = sq;

  return eslOK;

 ERROR:
  esl_fatal("Memory allocation error.");
  return status; /* NEVERREACHED */
}


/* trim_msa
 *                   
 * Given an MSA and unaligned 'trimmed' versions (subsequences) of all seqs in that MSA, 
 * replace all chars that have been trimmed away (not in subsequences) with gaps in the MSA.
 */
static int trim_msa(ESL_MSA *msa, ESL_SQ **sq, char *errbuf)
{
  int status;
  int i;
  int apos, uapos;
  int astart,  aend;
  int uastart, uaend;
  char *offset;
  char *aseq;
  char *uaseq;
  char *uasubseq;
  int *a2ua_map;
  int *ua2a_map;
  int ualen;

  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), msa must be digitized.");

  ESL_ALLOC(aseq,  sizeof(char) * (msa->alen+1));

  for(i = 0; i < msa->nseq; i++)
    {
      if (sq[i]->dsq == NULL) ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq's must be digitized.");
      if (sq[i]->n   == 0)    ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq[%d] is zero-length\n", i);

      ESL_ALLOC(a2ua_map, sizeof(int) * (msa->alen+1));
      esl_vec_ISet(a2ua_map, (msa->alen+1), -1);
      uapos = apos = 1;
      while(apos <= msa->alen)
	{
	  while(apos <= msa->alen && esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) apos++;
	  if(apos <= msa->alen) a2ua_map[apos] = uapos++;
	  apos++;
	}
      ualen = uapos;
      ESL_ALLOC(ua2a_map, sizeof(int) * (ualen+1));
      ua2a_map[0] = -1;
      for(apos = 1; apos <= msa->alen; apos++)
	if(a2ua_map[apos] != -1)
	  ua2a_map[a2ua_map[apos]] = apos;

      ESL_ALLOC(uasubseq, sizeof(char) * (sq[i]->n+1));
      esl_abc_Textize(msa->abc, sq[i]->dsq, sq[i]->n, uasubseq);
      esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, aseq);

      esl_strdup(aseq, -1, &(uaseq));
      esl_strdealign(uaseq, uaseq, "-_.", NULL);
      offset = strstr(uaseq, uasubseq);
      if(offset == NULL) ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq[%d] is not a subseq of msa seq %d\n", i, i);
      uastart = offset  - uaseq + 1;
      uaend   = uastart + strlen(uasubseq) - 1;
      astart  = ua2a_map[uastart];
      aend    = ua2a_map[uaend];
      free(ua2a_map);
      free(a2ua_map);

      for(apos = 1;        apos <  astart;    apos++) msa->ax[i][apos] = msa->abc->K; /* make it a gap */
      for(apos = aend + 1; apos <= msa->alen; apos++) msa->ax[i][apos] = msa->abc->K; /* make it a gap */
      free(uaseq);
      free(uasubseq);
    }

  for(i = 0; i < msa->nseq; i++)
    esl_sq_Destroy(sq[i]);
  free(sq);
      
  free(aseq);
  return eslOK;

 ERROR:
  return status;
}

/* dump_insert_info
 *                   
 * Given an MSA with RF annotation, print out information about how many 'insertions' come
 * after each non-gap RF column (consensus column). 
 */
static int dump_insert_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int **ict;
  int *total_ict, *med_ict;
  int i, l;
  int clen;
  int nseq;
  int *len;
  int *isize;
  float *ifract;
  int imax = 0;
  int ntotal = 0;
  int nins = 0;
  float cumulative = 0.;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_insert_info(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --iplot.");

  ESL_ALLOC(total_ict,  sizeof(int) * (msa->alen+2));
  ESL_ALLOC(med_ict,  sizeof(int) * (msa->alen+2));
  esl_vec_ISet(total_ict, (msa->alen+2), 0);
  esl_vec_ISet(med_ict, (msa->alen+2), 0);

  ESL_ALLOC(ict,  sizeof(int *) * (msa->alen+2));
  for(i = 0; i <= msa->alen; i++)
    {
      ESL_ALLOC(ict[i],  sizeof(int) * (msa->nseq));
      esl_vec_ISet(ict[i], (msa->nseq), 0);
    }

  fprintf(fp, "# %8s  %10s  %8s  %8s  %8s\n", "cons col", "nseq w/ins",  "freq ins", "avg len",  "med len");
  fprintf(fp, "# %8s  %10s  %8s  %8s  %8s\n", "--------", "----------", "--------", "--------", "--------");

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) 
	cpos++;
      else
	for(i = 0; i < msa->nseq; i++)
	  if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	    ict[cpos][i]++;
	    total_ict[cpos]++;
	  }	  
    }
  clen = cpos;

  /* determine avg median length for each insertion */
  for(cpos = 0; cpos <= clen; cpos++)
    {
      if(total_ict[cpos] > 0) { 
	nseq = 0;
	for(i = 0; i < msa->nseq; i++) { 
	  if(ict[cpos][i] >= 1) nseq++;
	}
	ESL_ALLOC(len, sizeof(int) * nseq);
	l = 0;
	for(i = 0; i < msa->nseq; i++) { 
	  if(ict[cpos][i] >= 1)
	    len[l++] = ict[cpos][i];
	}
	qsort(len, nseq, sizeof(int), compare_ints);
	med_ict[cpos] = len[nseq / 2];
	free(len);
      }      
    }
  for(cpos = 0; cpos <= clen; cpos++)
    {
      nseq = 0;
      for(i = 0; i < msa->nseq; i++) if(ict[cpos][i] >= 1) nseq++;
      if(nseq > 0) 
	fprintf(fp, "  %8d  %10d  %8.6f  %8.3f  %8d\n", cpos, nseq, (float) nseq / (float) msa->nseq, ((float) total_ict[cpos] / (float) nseq), med_ict[cpos]);
    }

  /* Possibly temporary: print plot of distribution of insert sizes, x value is insert size, 1..max ins length,
   * y is fraction of inserts accounted for by sizes <= x.
   */
  ESL_ALLOC(isize, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ifract, sizeof(float) * (msa->alen+1));
  esl_vec_ISet(isize, msa->alen+1, 0);
  for(cpos = 0; cpos <= clen; cpos++) {
    for(i = 0; i < msa->nseq; i++) { 
      if(ict[cpos][i] > 0) { 
	isize[ict[cpos][i]]++;
	imax = ESL_MAX(imax, ict[cpos][i]);
	nins++;
      }
    }
  }
  ntotal = esl_vec_ISum(total_ict, (msa->alen+2));

  for(i = 0; i <= msa->alen; i++) {
    ifract[i] = ((float) isize[i] * (float) i) / (float) ntotal;
  }

  printf("\n\n");
  printf("%d total inserted residues\n", ntotal);
  printf("%d total inserts (avg: %.3f)\n", (nins), (float) ntotal / (float) nins);
  for(i = 1; i <= imax; i++) {
    if(isize[i] != 0) { 
      cumulative += ifract[i];
      printf("%5d %5d %.5f %.8f %.8f\n", i, isize[i], (float) isize[i] / (float) nins, ifract[i], cumulative);
    }  
  }
  free(isize);
  free(ifract);

  for(i = 0; i <= msa->alen; i++)
    free(ict[i]);
  free(ict);
  free(total_ict);
  free(med_ict);

  return eslOK;

 ERROR:
  return status;
}


/* dump_residue_info
 *                   
 * Given an MSA, print out the number of sequences with
 * a non-gap residue in each column of the alignment.
 */
static int dump_residue_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int rct;
  int i;
  int has_rf;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_residue_info(), msa must be digitized.");
  has_rf = (msa->rf == NULL) ? FALSE : TRUE;

  if(has_rf) { 
    fprintf(fp, "# %8s  %7s  %8s  %8s\n", "cons col", "aln col", "num res",  "freq res");
    fprintf(fp, "# %8s  %7s  %8s  %8s\n", "--------", "-------", "--------", "--------");
  }  
  else { 
    fprintf(fp, "# %7s  %8s  %8s\n", "aln col", "num res",  "freq res");
    fprintf(fp, "# %7s  %8s  %8s\n", "-------", "--------", "--------");
  }
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    rct = 0;
    if(has_rf && (! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)]))) cpos++;
    for(i = 0; i < msa->nseq; i++) { 
      if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) rct++; 
    }

    if(has_rf) fprintf(fp, "  %8d  %7d  %8d  %8.6f\n", cpos, apos, rct, (float) rct / (float) msa->nseq);
    else       fprintf(fp, "  %7d  %8d  %8.6f\n", apos, rct, (float) rct / (float) msa->nseq);
  }

  return eslOK;

 ERROR:
  return status;
}


/* dump_cres_info
 *                   
 * Given an MSA, for each sequence <x>, print out the number of 
 * columns for which sequence <x> is the ONLY sequence with 
 * a residue in the column.
 * Originally written to pick out which sequences are solely responsible
 * for the most columns. 
 */
static int dump_cres_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos;
  int rct;
  int i;
  int lasti = -1;
  int *cres;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_cres_info(), msa must be digitized.");

  fprintf(fp, "# %-30s  %7s\n", "seq name", "ncols");
  fprintf(fp, "# %-30s  %7s\n", "------------------------------", "-------");

  ESL_ALLOC(cres, sizeof(int) * msa->nseq);
  esl_vec_ISet(cres, msa->nseq, 0);
  for(apos = 1; apos <= msa->alen; apos++) {
    rct = 0;
    for(i = 0; i < msa->nseq; i++) { 
      if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	rct++; 
	lasti = i;
      }
    }
    if(rct == 1) { /* special case, lasti will be the only i of the only seq with a residue in apos */
      cres[lasti]++;
    }
  }
  for(i = 0; i < msa->nseq; i++) { 
    if(cres[i] > 0) { 
      fprintf(fp, "  %-30s  %7d\n", msa->sqname[i], cres[i]);
    }
  }
  free(cres);

  return eslOK;

 ERROR:
  return status;
}


/* dump_delete_info
 *                   
 * Given an MSA, print out the number of sequences with
 * gaps residue in each consensus column of the alignment.
 */
static int dump_delete_info(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int dct;
  int i;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in dump_residue_info(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --dinfo.");

  fprintf(fp, "# Number of sequences in file: %d\n", msa->nseq);
  fprintf(fp, "# Only non-gap RF columns with > 0 deletes are listed.\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# %8s  %7s  %8s  %8s\n", "cons col", "aln col", "num del",  "freq del");
  fprintf(fp, "# %8s  %7s  %8s  %8s\n", "--------", "-------", "--------", "--------");
  
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    dct = 0;
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
      cpos++;
      for(i = 0; i < msa->nseq; i++) { 
	if(esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) dct++; 
      }
      if(dct > 0) fprintf(fp, "  %8d  %7d  %8d  %8.6f\n", cpos, apos, dct, (float) dct / (float) msa->nseq);
    }
  }
  return eslOK;

 ERROR:
  return status;
}

/* plot_inserts
 *                   
 * Given an MSA with RF annotation, print a postscript heatmap of how
 * many insertions are after each non-gap RF column (consensus column)
 * in each sequence.
 */
static int plot_inserts(FILE *fp, ESL_MSA *msa, int do_log, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  ESL_DMATRIX *I;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --iplot.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in plot_inserts(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  I = esl_dmatrix_Create(msa->nseq, (clen+1));
  esl_dmatrix_SetZero(I);

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) 
      cpos++;
    else
      for(i = 0; i < msa->nseq; i++)
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) I->mx[i][cpos] += 1.;
  }
  if(do_log) {
    for(i = 0; i < msa->nseq; i++)
      for(cpos = 0; cpos <= clen; cpos++)
	if(I->mx[i][cpos] > 0) 
	  I->mx[i][cpos] = log(I->mx[i][cpos]);
	else 
	  I->mx[i][cpos] = -1; /* don't want 0s to change to -inf, for coloring scheme */
  }
  else { 
    for(i = 0; i < msa->nseq; i++)
      for(cpos = 0; cpos <= clen; cpos++)
	if(I->mx[i][cpos] == 0) I->mx[i][cpos] = (-1 * esl_dmx_Max(I)) / 2; /* for better resolution on heatmap */
  }

  /* dmx_Visualize(fp, I, esl_dmx_Min(I), esl_dmx_Max(I)); */
  dmx_Visualize(fp, I, (-1 * esl_dmx_Max(I)), esl_dmx_Max(I));
  /* esl_dmatrix_Dump(stdout, I, NULL, NULL); */
  esl_dmatrix_Destroy(I);
  return eslOK;

 ERROR:
  return status;
}


/* plot_gaps
 *                   
 * Given an MSA with RF annotation, print a postscript checkboard grid 
 * showing which sequences have gaps in each non-gap RF column. 
 */
static int plot_gaps(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  ESL_DMATRIX *G;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --gplot.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in plot_gaps(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  G = esl_dmatrix_Create(msa->nseq, (clen+1));
  esl_dmatrix_SetZero(G);

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) {
      cpos++;
      for(i = 0; i < msa->nseq; i++)
	if(esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) G->mx[i][cpos] += 1.;
    }
  }
  /* dmx_Visualize(fp, G, esl_dmx_Min(G), esl_dmx_Max(G)); */
  dmx_Visualize(fp, G, -1, 1);
  esl_dmatrix_Destroy(G);
  return eslOK;

 ERROR:
  return status;
}

/* get_tree_order
 *                   
 * Given a tree, determine the branching order of the sequences
 * it represents by traversing it preorder.
 */
static int get_tree_order(ESL_TREE *T, char *errbuf, int **ret_order)
{
  int status;
  int opos = 0;
  int nd;
  int *order; 
  ESL_STACK *pda;
  ESL_ALLOC(order, sizeof(int) * T->N);

  opos = 0;
  pda  = esl_stack_ICreate();
  esl_stack_IPush(pda, T->right[0]);
  esl_stack_IPush(pda, T->left[0]);
  while (esl_stack_IPop(pda, &nd) != eslEOD)
    {
      if (nd > 0) { /* a node */
	esl_stack_IPush(pda, T->right[nd]); /* index for right child */
	esl_stack_IPush(pda, T->left[nd]);  /* index for left child */
      }
      else /* nd <= 0, a child */
	order[opos++] = nd * -1;
    }
  *ret_order = order;
  esl_stack_Destroy(pda);
  return eslOK;

 ERROR:
  return status;
}

/* reorder_msa
 *                   
 * Given an array specifying a new order for the sequences in
 * the MSA, reorder it by swapping pointers.
 */
static int
reorder_msa(ESL_MSA *msa, int *order, char *errbuf)
{
  int status;
  char **tmp; 
  ESL_ALLOC(tmp, sizeof(char *) * msa->nseq);
  int i, a;

  /* contract check */
  /* 'order' must be have nseq elements, elements must be in range [0..nseq-1], no duplicates  */
  int *covered;
  ESL_ALLOC(covered, sizeof(int) * msa->nseq);
  esl_vec_ISet(covered, msa->nseq, 0);
  for(i = 0; i < msa->nseq; i++) { 
    if(covered[order[i]]) ESL_FAIL(eslEINVAL, errbuf, "reorder_msa() order array has duplicate entries for i: %d\n", i);
    covered[order[i]] = 1;
  }
  free(covered);

  /* swap aseq or ax (one or the other must be non-NULL) */
  if(msa->flags & eslMSA_DIGITAL) { /* digital MSA */
    ESL_DSQ **tmp_dsq; 
    ESL_ALLOC(tmp_dsq, sizeof(ESL_DSQ *) * msa->nseq);
    for(i = 0; i < msa->nseq; i++) tmp_dsq[i] = msa->ax[i];
    for(i = 0; i < msa->nseq; i++) msa->ax[i] = tmp_dsq[order[i]];
    free(tmp_dsq);
  }
  else { /* text MSA */
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->aseq[i];
    for(i = 0; i < msa->nseq; i++) msa->aseq[i] = tmp[order[i]];
  }

  /* swap sqnames (mandatory) */
  for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqname[i];
  for(i = 0; i < msa->nseq; i++) msa->sqname[i] = tmp[order[i]];

  /* swap sqacc, if they exist */
  if(msa->sqacc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqacc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqacc[i] = tmp[order[i]];
  }

  /* swap sqdesc, if they exist */
  if(msa->sqdesc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqdesc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqdesc[i] = tmp[order[i]];
  }

  /* swap ss, if they exist */
  if(msa->ss != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->ss[i];
    for(i = 0; i < msa->nseq; i++) msa->ss[i] = tmp[order[i]];
  }

  /* swap sa, if they exist */
  if(msa->sa != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sa[i];
    for(i = 0; i < msa->nseq; i++) msa->sa[i] = tmp[order[i]];
  }

  /* swap pp, if they exist */
  if(msa->pp != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->pp[i];
    for(i = 0; i < msa->nseq; i++) msa->pp[i] = tmp[order[i]];
  }

  /* swap gs annotation, if it exists */
  for(a = 0; a < msa->ngs; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gs[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gs[a][i] = tmp[order[i]];
  }

  /* swap gr annotation, if it exists */
  for(a = 0; a < msa->ngr; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gr[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gr[a][i] = tmp[order[i]];
  }
  free(tmp);
  return eslOK;

 ERROR: 
  return status;
}

/****************************************************************
 * Stolen from hmmer/h3/heatmap.c SVN revision 2171
 * as dmx_Visualize. Then modified so that the full 
 * matrix is printed (not half split diagonally).
 */
/* dmx_Visualize()
 * Incept:    SRE, Wed Jan 24 11:58:21 2007 [Janelia]
 *
 * Purpose:   
 *            
 *            Color scheme roughly follows Tufte, Envisioning
 *            Information, p.91, where he shows a beautiful
 *            bathymetric chart. The CMYK values conjoin two
 *            recommendations from ColorBrewer (Cindy Brewer
 *            and Mark Harrower) 
 *            [http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer.html],
 *            specifically the 9-class sequential2 Blues and
 *            9-class sequential YlOrBr.
 * 
 *            Might eventually become part of Easel, once mature?
 *           
 * Note:      Binning rules basically follow same convention as
 *            esl_histogram. nb = xmax-xmin/w, so w = xmax-xmin/nb; 
 *            picking bin is (int) ceil((x - xmin)/w) - 1. (xref
 *            esl_histogram_Score2Bin()). This makes bin b contain
 *            values bw+min < x <= (b+1)w+min. (Which means that 
 *            min itself falls in bin -1, whoops - but we catch
 *            all bin<0 and bin>=nshades and put them in the extremes.
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max)
 {
   int    nshades   = 18;
   double cyan[]    = { 1.00, 1.00, 0.90, 0.75, 0.57, 0.38, 0.24, 0.13, 0.03,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.07, 0.20, 0.40, 0.60};
   double magenta[] = { 0.55, 0.45, 0.34, 0.22, 0.14, 0.08, 0.06, 0.03, 0.01,
			0.00, 0.03, 0.11, 0.23, 0.40, 0.55, 0.67, 0.75, 0.80};
   double yellow[]  = { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.10, 0.25, 0.40, 0.65, 0.80, 0.90, 1.00, 1.00, 1.00};
   double black[]   = { 0.30, 0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
   double w;			
   int    i,j;
   int    bin;
   int    boxsize;		/* box size in points */
   int    xcoord, ycoord;	/* postscript coords in points */
   int    leftmargin, rightmargin;
   int    bottommargin, topmargin;
   float  fboxsize;		/* box size in fractional points */

   /* Set some defaults that might become arguments later.
    */
   leftmargin   = rightmargin = 20;
   bottommargin = topmargin   = 20;

   /* Determine some working parameters 
    */
   w = (max-min) / (double) nshades; /* w = bin size for assigning values->colors*/
   boxsize = ESL_MAX(1, (ESL_MIN((792 - bottommargin) / D->n, 
				 (612 - leftmargin)   / D->m)));
   fboxsize= ESL_MIN( (792. - ((float) bottommargin + topmargin))   / (float) D->n, 
		      (612. - ((float) leftmargin   + rightmargin)) / (float) D->m);


   fprintf(fp, "%.4f %.4f scale\n", (fboxsize/(float) boxsize), (fboxsize/(float) boxsize));
   /* printf("n: %d\nm: %d\n", D->n, D->m); */
   for (i = 0; i < D->n; i++) {
     /* printf("\n"); */
     /* for (j = i; j < D->n; j++) */
     for (j = 0; j < D->m; j++)
       {
	 /* printf("i: %4d j: %4d %5.1f\n", i, j, D->mx[i][j]); */
	 xcoord = j * boxsize + leftmargin;
	 ycoord = (D->m-(i+1)) * boxsize + bottommargin; /* difference w/heatmap.c: (D->m-i+1) */
	 
	 if      (D->mx[i][j] == -eslINFINITY) bin = 0;
	 else if (D->mx[i][j] ==  eslINFINITY) bin = nshades-1;
	 else {
	   bin    = (int) ceil((D->mx[i][j] - min) / w) - 1;
	   if (bin < 0)        bin = 0;
	   if (bin >= nshades) bin = nshades-1;
	 }
	 
	fprintf(fp, "newpath\n");
	fprintf(fp, "  %d %d moveto\n", xcoord, ycoord);
	fprintf(fp, "  0  %d rlineto\n", boxsize);
	fprintf(fp, "  %d 0  rlineto\n", boxsize);
	fprintf(fp, "  0 -%d rlineto\n", boxsize);
	fprintf(fp, "  closepath\n");
  	fprintf(fp, " %.2f %.2f %.2f %.2f setcmykcolor\n",
		cyan[bin], magenta[bin], yellow[bin], black[bin]);
	fprintf(fp, "  fill\n");
      }
   }
  fprintf(fp, "showpage\n");
  return eslOK;
}


/* read_mask_file
 *
 * Given an open file pointer, read the first token of the
 * file and return it as *ret_mask.
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;

  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));
  for(n = 0; n < toklen; n++) mask[n] = tok[n];
  mask[n] = '\0';
  *ret_mask = mask;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}


/* map_msas
 *                   
 * For each non-gap RF column in msa1, determine the corresponding column
 * in msa2. This implementation requires:
 *  - msa1 and msa2 contain exactly the same sequences in the same order
 *  - msa non-gap RF len is < msa2->alen.
 * Note: the seqs in msa1 and msa2 do not have to have the same names.
 *
 * Uses a DP algorithm similar to Needleman-Wunsch, but simpler because
 * we require that all non-gap RF columns from msa1 must map to exactly 1
 * column in msa2, i.e. from a Needleman-Wunsch perspective, we can't have
 * gaps in one of our 'sequences'.
 */
static int
map_msas(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, char **ret_msa1_to_msa2_map)
{
  int status;
  int *c2a_map1 = NULL;       /* msa1 map of consensus columns (non-gap RF residues) to alignment columns */
  int clen1;                  /* consensus (non-gap RF) length of msa1 */
  int **one2two;              /* [0..c..clen1][0..a..msa2->alen] number of residues from non-gap RF column c of msa1
			       * aligned in column a of msa 2 */
  int apos1, apos2;           /* counters over alignment position in msa1, msa2 respectively */
  int cpos1;                  /* counter over non-gap RF (consensus) position in msa1 */
  int diagonal;               /* score for diagonal move in dp recursion */
  int vertical;               /* score for vertical move in dp recursion */
  int **mx;                   /* [0..c..clen1][0..a..msa2->alen] dp matrix, score of max scoring aln 
			       * from 1..c in msa1 and 1..a in msa 2 */
  int **tb;                   /* [0..c..clen1][0..a..msa2->alen] traceback ptrs, 0 for diagonal, 1 for vertical */
  char *seq1, *seq2;          /* temporary strings for ensuring dealigned sequences in msa1 and msa2 are identical */
  int64_t len1, len2;         /* length of seq1, seq2 */
  int isgap1, isgap2;         /* is this residue a gap in msa1, msa2? */
  int i;                      /* counter over sequences */
  int *res1_per_cpos;         /* [0..c..clen1] number of residues in non-gap RF column c of msa1 */
  int max_sc;                 /* max score of full path (alignment) through dp mx */
  int max_apos2;              /* temp val for finding endpoint of alignment */
  int tb_sc;                  /* score of traceback, should equal max_sc */
  int *one_rf2two_map;        /* [0..c..clen1] the alignment, msa2 column that non-gap RF column c in msa1 maps to */
  int total_cres1 = 0;        /* total number of residues in msa1 that are in non-gap RF columns */
  float coverage;             /* fraction of total_cres1 that are within mapped msa2 columns from one_rf2two_map, 
			       * this is tb_sc / total_cres1 */
  char *msa1_to_msa2_map;     /* map from msa1 to msa2, this is the optimal alignment found by dp */
  int total_msa1_res = 0;     /* total number of residues in MSA1, we use -1 * this value to initialize dp matrix */
  int be_verbose = esl_opt_GetBoolean(go, "--verbose");

  /* contract check */
  if(msa1->rf == NULL)                 ESL_FAIL(eslEINVAL, errbuf, "with --map %s must have RF annotation.", esl_opt_GetArg(go, 1));
  if(! (msa1->flags & eslMSA_DIGITAL)) ESL_FAIL(eslEINVAL, errbuf, "in map_msas() msa1 (%s) not digitized.\n", esl_opt_GetArg(go, 1));
  if(! (msa2->flags & eslMSA_DIGITAL)) ESL_FAIL(eslEINVAL, errbuf, "in map_msas() msa2 (%s) not digitized.\n", esl_opt_GetString(go, "--map"));
  
  /* Map msa1 non-gap RF (consensus) columns to alignment positions */
  if((status = map_cpos_to_apos(msa1, &c2a_map1, &clen1))   != eslOK) goto ERROR;
  if(clen1 > msa2->alen) ESL_XFAIL(eslEINVAL, errbuf, "non-gap RF length of msa in <msafile> %s (%d) is greater than --map alignment length of %s (%" PRId64 ").", esl_opt_GetArg(go, 1), clen1, esl_opt_GetString(go, "--map"), msa2->alen);
  if(be_verbose) { 
    printf("%25s non-gap RF (consensus) length: %d\n",          esl_opt_GetArg(go, 1), clen1);
    printf("%25s alignment length:              %" PRId64 "\n", esl_opt_GetString(go, "--map"), msa2->alen);
  }
  /* collect counts in one2two[i][j]: number of sequences for which residue aligned in msa1 non-gap column i
   * is aligned in msa2 alignment column j.
   */
  ESL_ALLOC(seq1, sizeof(char) * (msa1->alen+1));
  ESL_ALLOC(seq2, sizeof(char) * (msa2->alen+1));
  ESL_ALLOC(one2two, sizeof(int *) * (msa1->alen+1));
  for(apos1 = 0; apos1 <= msa1->alen; apos1++) { 
    ESL_ALLOC(one2two[apos1], sizeof(int) * (msa2->alen+1));
    esl_vec_ISet(one2two[apos1], (msa2->alen+1), 0);
  }

  total_msa1_res = 0;
  for(i = 0; i < msa1->nseq; i++) { 
    /* ensure raw (unaligned) seq i in the 2 msas is the same */
    esl_abc_Textize(msa1->abc, msa1->ax[i], msa1->alen, seq1); 
    esl_abc_Textize(msa1->abc, msa2->ax[i], msa2->alen, seq2); /* note: msa*1*->abc used on purpose, allows DNA/RNA to peacefully coexist in this func */
    esl_strdealign(seq1, seq1, "-_.", &len1);
    esl_strdealign(seq2, seq2, "-_.", &len2);

    if(len1 != len2) { 
      ESL_FAIL(eslEINVAL, errbuf, "--map error: unaligned seq number %d (msa1: %s, msa2: %s) differs in length %s (%" PRId64 ") and %s (%" PRId64 "), those files must contain identical raw seqs\n",
	       i, msa1->sqname[i], msa2->sqname[i], esl_opt_GetArg(go, 1), len1, esl_opt_GetString(go, "--map"), len2);
    }
    if(strncmp(seq1, seq2, len1) != 0)  ESL_FAIL(eslEINVAL, errbuf, "--map error: unaligned seq number %d differs between %s and %s, those files must contain identical raw seqs\n", i, esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--map"));
    total_msa1_res += len1;
    
    apos1 = apos2 = 1;
    while((apos1 <= msa1->alen) || (apos2 <= msa2->alen)) {
      isgap1 = esl_abc_XIsGap(msa1->abc, msa1->ax[i][apos1]);
      isgap2 = esl_abc_XIsGap(msa2->abc, msa2->ax[i][apos2]);
      if      ( isgap1 &&  isgap2) { apos1++; apos2++; }
      else if ( isgap1 && !isgap2) { apos1++;          }
      else if (!isgap1 &&  isgap2) {          apos2++; }
      else if ( msa1->ax[i][apos1] == msa2->ax[i][apos2]) { 
	one2two[apos1++][apos2++]++;
	/* two2one[apos2][apos1]++; */
      }
    }
  }

  /******************************************************************
   * DP alignment of msa1 to msa2
   * dp matrix: mx[cpos1][apos2] cpos1=1..clen1, apos2=1..msa2->alen (cpos1=0 || apos2=0 is invalid)
   * mx[cpos1][apos2] = score of maximal alignment for cpos1'=1..cpos1, apos2'=1..apos2 INCLUDING
   *                    cpos1 and apos2. Score is number of residues from msa1 consensus columns
   *                    1..cpos1 that exist in their respective aligned columns in msa2 (the growing
   *                    maximally scoring alignment).
   */

  /******************************************************************
   * initialization 
   */
  ESL_ALLOC(mx, sizeof(int *) * (clen1+1));
  ESL_ALLOC(tb, sizeof(int *) * (clen1+1));
  for(cpos1 = 0; cpos1 <= clen1; cpos1++) { 
    ESL_ALLOC(mx[cpos1], sizeof(int) * (msa2->alen+1));
    ESL_ALLOC(tb[cpos1], sizeof(int) * (msa2->alen+1));
    esl_vec_ISet(mx[cpos1], (msa2->alen+1), -1 * (total_msa1_res + 1)); /* initialize to worst possible score we could have, - 1,
									 * this was a bug before, if we init to 0, we can
									 * get alignments that don't go all the back to cpos1 = 1,
									 * but stop at say cpos1 = 2 because cpos1[2][1] was set as
									 * 0, even though it should be impossible. */
    esl_vec_ISet(tb[cpos1], (msa2->alen+1), -2); /* -2 is a bogus value, if we see it during traceback, there's a problem */
  }
  ESL_ALLOC(res1_per_cpos, sizeof(int) * (clen1+1));
  esl_vec_ISet(res1_per_cpos, (clen1+1), 0);

  mx[1][1] = one2two[c2a_map1[1]][1];
  tb[1][1] = -1; /* last cell, special value */

  /* initialize on cpos = 1, no choice, must come from vertical move */
  cpos1 = 1;
  apos1 = c2a_map1[cpos1];
  res1_per_cpos[cpos1] = one2two[apos1][1];
  for(apos2 = 2; apos2 <= msa2->alen; apos2++) {
    mx[cpos1][apos2] = mx[cpos1][(apos2-1)] - one2two[apos1][(apos2-1)] + one2two[apos1][apos2];
    tb[cpos1][apos2] = 1; /* vertical move */
    res1_per_cpos[cpos1] += one2two[apos1][apos2];
  }
  /*****************************************************************
   * recursion
   */
  for(cpos1 = 2; cpos1 <= clen1; cpos1++) {
    apos1 = c2a_map1[cpos1];
    res1_per_cpos[cpos1] = one2two[apos1][1];
    for(apos2 = 2; apos2 <= msa2->alen; apos2++) {
      /* only one msa2 column apos2 can align to each msa1 consensus column, 
       * if we take the vertical step, we're saying it wasn't apos2-1, it's apos2
       * that maps to cpos1, so we have to subtract out the score that apos2-1 
       * mapped to cpos that was added into
       * mx[cpos1][g-1] on previous recursion step. 
       */
      vertical = mx[ cpos1   ][(apos2-1)] - one2two[apos1][(apos2-1)]; 
      diagonal = mx[(cpos1-1)][(apos2-1)];
      if(diagonal >= vertical) {
	mx[cpos1][apos2] = diagonal;
	tb[cpos1][apos2] = 0; /* diagonal move */
      } 
      else {
	mx[cpos1][apos2] = vertical;
	tb[cpos1][apos2] = 1; /* vertical move */
      }
      mx[cpos1][apos2]     += one2two[apos1][apos2];
      res1_per_cpos[cpos1] += one2two[apos1][apos2];
    }
  }

  /*****************************************************************
   * traceback 
   */
  
  /* need to find end point, cpos1=clen1, apos2=argmax_apos2 dp[cpos1][apos2] + one2two[c2a_map1[cpos1]][apos2] */ 
  max_sc    = mx[clen1][1];
  max_apos2 = 1;
  apos1 = c2a_map1[clen1];
  for(apos2 = 2; apos2 <= msa2->alen; apos2++) {
    if((mx[clen1][apos2]) > max_sc) { 
      max_sc    = mx[clen1][apos2];
      max_apos2 = apos2;
    }
  }
  if(be_verbose) printf("max score %d\nmax apos2 %d\n", max_sc, max_apos2);
  ESL_ALLOC(one_rf2two_map, sizeof(int) * (clen1+1));

  /* traceback, and build one_rf2two_map[] */
  apos2 = max_apos2;
  cpos1 = clen1;
  one_rf2two_map[cpos1] = apos2;
  tb_sc = one2two[apos1][apos2];
  if     (be_verbose && res1_per_cpos[cpos1] == 0) printf("1 cc %4d --> 2 %4d %5d / %5d (%.4f)\n", cpos1, apos2, one2two[apos1][apos2], res1_per_cpos[cpos1], (float) 0.);
  else if(be_verbose)                              printf("1 cc %4d --> 2 %4d %5d / %5d (%.4f)\n", cpos1, apos2, one2two[apos1][apos2], res1_per_cpos[cpos1], ((float) one2two[apos1][apos2] / (float) res1_per_cpos[cpos1])); 

  total_cres1 = 0;
  apos1 = c2a_map1[cpos1];
  while(tb[cpos1][apos2] != -1) {
    if(tb[cpos1][apos2] == 0) { /* diagonal move */
      if(tb[cpos1][apos2] != -1) { 
	cpos1--; apos2--;
	apos1 = c2a_map1[cpos1];
	one_rf2two_map[cpos1] = apos2;
	if(be_verbose && res1_per_cpos[cpos1] == 0) printf ("1 cc %4d --> 2 %4d %5d / %5d (0.0000)\n", cpos1, apos2, one2two[apos1][apos2], res1_per_cpos[cpos1]);
	else {
	  if(be_verbose) printf("1 cc %4d --> 2 %4d %5d / %5d (%.4f)\n", cpos1, apos2, one2two[apos1][apos2], res1_per_cpos[cpos1], ((float) one2two[apos1][apos2] / (float) res1_per_cpos[cpos1])); 
	  total_cres1 += res1_per_cpos[cpos1];
	}
	tb_sc += one2two[apos1][apos2];
      }
    }
    else if(tb[cpos1][apos2] == 1)  apos2--; /* vertical move */
    else if(tb[cpos1][apos2] != -1) /* shouldn't happen */
      ESL_FAIL(eslEINVAL, errbuf, "--map error: in dp traceback, tb[cpos1: %d][apos2: %d] %d\n", cpos1, apos2, tb[cpos1][apos2]);
  }
  total_cres1 += res1_per_cpos[cpos1];
  /* done DP code 
   **********************************/

  if(be_verbose) printf("Total trace back sc: %d\n", tb_sc);
  if(tb_sc != max_sc) ESL_FAIL(eslEINVAL, errbuf, "--map error: in dp traceback, tb_sc (%d) != max_sc (%d)\n", tb_sc, max_sc);
  coverage = (float) tb_sc / (float) total_cres1;
  printf("Coverage: %6d / %6d (%.4f)\nCoverage is fraction of consensus residues from %s in optimally mapped columns in %s\n", tb_sc, total_cres1, coverage, esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--map"));

  /* create 1/0 mask */
  ESL_ALLOC(msa1_to_msa2_map, sizeof(char) * (msa2->alen+1));
  apos2 = 1;
  for(cpos1 = 1; cpos1 <= clen1; cpos1++) {
    while(apos2 < one_rf2two_map[cpos1]) { msa1_to_msa2_map[(apos2-1)] = '0'; apos2++; }
    msa1_to_msa2_map[(apos2-1)] = '1'; 
    apos2++; 
  }
  while(apos2 <= msa2->alen) { msa1_to_msa2_map[(apos2-1)] = '0'; apos2++; }

  msa1_to_msa2_map[msa2->alen] = '\0';
  *ret_msa1_to_msa2_map = msa1_to_msa2_map;

  /* clean up and return */
  for(cpos1 = 0; cpos1 <= clen1; cpos1++) { 
    free(mx[cpos1]);
    free(tb[cpos1]);
  }
  free(mx);
  free(tb);

  for(apos1 = 0; apos1 <= msa1->alen; apos1++) free(one2two[apos1]);
  free(one2two);
  free(one_rf2two_map);
  free(res1_per_cpos);
  free(c2a_map1);

  free(seq1);
  free(seq2);
  return eslOK;
  
 ERROR: 
  return status;
}


/* map_sub_msas
 *                   
 * msa1 and msa2 contain the same named sequences, msa1 contains a superset 
 * of the columns in msa2. Determine which of the msa1 columns the msa2
 * columns correspond to.
 */
static int
map_sub_msas(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa1, ESL_MSA *msa2, char **ret_msa1_to_msa2_mask)
{
  int status;
  int  apos1, apos2;          /* counters over alignment position in msa1, msa2 respectively */
  int i;
  int *msa1_to_msa2_map;    /* [0..apos1..msa1->alen] msa2 alignment position that apos1 corresponds to */
  char *mask;

  /* contract check */
  if(! (msa1->flags & eslMSA_DIGITAL)) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas() msa1 (%s) not digitized.\n", esl_opt_GetArg(go, 1));
  if(! (msa2->flags & eslMSA_DIGITAL)) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas() msa2 (%s) not digitized.\n", esl_opt_GetString(go, "--submap"));
  if(msa1->alen <= msa2->alen) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas() alignment length for msa1 (%" PRId64 "d) <= length for msa2 (%" PRId64 ")\n", msa1->alen, msa2->alen);
  
  ESL_ALLOC(mask, sizeof(char) * (msa1->alen+1));
  for(apos1 = 0; apos1 < msa1->alen; apos1++) mask[apos1] = '0';
  mask[msa1->alen] = '\0';

  ESL_ALLOC(msa1_to_msa2_map, sizeof(int) * msa1->alen+1);
  esl_vec_ISet(msa1_to_msa2_map, (msa1->alen+1), -1);

  /* both alignments must have same 'named' sequences in same order */
  if(msa1->nseq != msa2->nseq) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas() msa1 has %d sequences, msa2 has %d sequences\n", msa1->nseq, msa2->nseq);
  for(i = 0; i < msa1->nseq; i++) { 
    if(strcmp(msa1->sqname[i], msa2->sqname[i]) != 0) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas() msa1 seq %d is named %s, msa2 seq %d is named %s\n", i, msa1->sqname[i], i, msa2->sqname[i]);
  }

  apos1 = 1;
  apos2 = 1;
  while((apos2 <= msa2->alen) || (apos1 <= msa1->alen)) { /* determine which apos1 (alignment column in msa1), apos2 (alignment column in msa2) corresponds to */
    for(i = 0; i < msa1->nseq; i++) { 
      if(msa1->ax[i][apos1] != msa2->ax[i][apos2]) { 
	apos1++; 
	break; /* try next apos1 */ 
      }
    }	
    if(i == msa1->nseq) { /* found a match */
      msa1_to_msa2_map[apos1] = apos2;
      mask[(apos1-1)] = '1';
      apos1++;
      apos2++;
    }
  }
  if((apos1 != (msa1->alen+1)) || (apos2 != (msa2->alen+1))) ESL_FAIL(eslEINVAL, errbuf, "in map_sub_msas(), failure mapping alignments, end of loop apos1-1 = %d (msa1->alen: %" PRId64 ") and apos2-1 = %d (msa2->alen: %" PRId64 ")\n", apos1-1, msa1->alen, apos2-1, msa2->alen);

  free(msa1_to_msa2_map);
  *ret_msa1_to_msa2_mask = mask;
  return eslOK;
  
 ERROR: 
  return status;
}

/* handle_post_opts
 *                   
 * Read "#=GR POST" annotation into a 2D matrix, each sequence
 * is a row, each residue is a column. Handle any command line
 * options that use the posterior info.
 *
 */      
static int handle_post_opts(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int    status;
  int    s,c;                 /* counters over sequences, columns of MSA */
  int   *nongap_c, *nongap_s; /* number of non-gap posterior values for each column/sequence respectively */
  float *sum_c, *sum_s;       /* sum of non-gap posterior values for each column/sequence respectively */
  float *min_c, *min_s;       /* min of non-gap posterior values for each column/sequence respectively */
  float *avg_c, *avg_s;       /* average non-gap posterior values for each column/sequence respectively */
  int   *athresh_c;           /* [0..c..msa->alen-1] number of sequences with residue in this column with post value >= pthresh */
  float *athresh_fract_c;     /* [0..c..msa->alen-1] fraction of non-gap sequences with residue in this column with post value >= pthresh */
  int    ridx1, ridx2;
  int    r;
  float  p;
  int    do_pfract = esl_opt_IsOn(go, "--pfract");
  int    do_prf    = esl_opt_IsOn(go, "--p-rf");
  int    do_pinfo  = esl_opt_IsOn(go, "--pinfo");
  float  pfract;
  float  pthresh   = esl_opt_GetReal(go, "--pthresh"); /* default is 0.95 */
  int  *useme; 
  int *c2a_map;
  int clen;
  int cpos;
  int nkept;
  int ir1, ir2;
  int ndigits;
  int nongap_total = 0;
  int nongap_total_rf = 0;
  int sum_total = 0;
  int sum_total_rf = 0;
  FILE *pinfofp = NULL;  /* output file for --pinfo */

  if((!do_pfract) && (!do_pinfo)) ESL_FAIL(eslEINVAL, errbuf, "handle_post_opts(): --pinfo nor --pfract options selected, shouldn't be in this function.");

  /* Find out which #=GR line is the POST, Post, or post line (if more than one exist, last one is chosen) */
  ridx1 = -1;
  ridx2 = -1;
  ndigits = 0;
  for (r = 0; r < msa->ngr; r++) { 
    if (strcmp(msa->gr_tag[r], "POST")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "Post")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "post")   == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "POSTX.") == 0) { ridx1 = r; ndigits = 1; }
    if (strcmp(msa->gr_tag[r], "POST.X") == 0) { ridx2 = r; ndigits = 2; }
  }
  if(ndigits == 1 && ridx1 == -1) { 
    if(do_pfract) ESL_FAIL(eslEINVAL, errbuf, "--pfract requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
    if(do_pinfo)  ESL_FAIL(eslEINVAL, errbuf, "--pinfo  requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  if(ndigits == 2 && (ridx1 == -1 || ridx2 == -1)) { 
    if(do_pfract) ESL_FAIL(eslEINVAL, errbuf, "--pfract requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
    if(do_pinfo)  ESL_FAIL(eslEINVAL, errbuf, "--pinfo  requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  if(msa->rf == NULL) { 
    if(do_prf)    ESL_FAIL(eslEINVAL, errbuf, "--p-rf requires \"#=GC RF\" annotation in %s.\n", esl_opt_GetArg(go,1));
  }
  /*ESL_ALLOC(post_sc, sizeof(int *) * msa->nseq);
    for(i = 0; i < msa->nseq; i++) 
    ESL_ALLOC(post_sc, sizeof(int) * msa->alen);
    ESL_ALLOC(post_cs, sizeof(int *) * msa->alen);
    for(i = 0; i < msa->alen; i++) 
    ESL_ALLOC(post_cs, sizeof(int) * msa->nseq);
  */
     
  ESL_ALLOC(nongap_c, sizeof(int) * msa->alen);
  ESL_ALLOC(sum_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(min_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(avg_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(athresh_c,sizeof(int)   * msa->alen);
  ESL_ALLOC(athresh_fract_c,sizeof(float) * msa->alen);
  esl_vec_ISet(nongap_c, msa->alen, 0);
  esl_vec_FSet(sum_c,    msa->alen, 0.);
  esl_vec_FSet(min_c,    msa->alen, 10.);
  esl_vec_ISet(athresh_c,msa->alen, 0);

  ESL_ALLOC(nongap_s, sizeof(int) * msa->nseq);
  ESL_ALLOC(sum_s, sizeof(float) * msa->nseq);
  ESL_ALLOC(min_s, sizeof(float) * msa->nseq);
  ESL_ALLOC(avg_s, sizeof(float) * msa->nseq);
  esl_vec_ISet(nongap_s, msa->nseq, 0);
  esl_vec_FSet(sum_s,    msa->nseq, 0.);
  esl_vec_FSet(min_s,    msa->nseq, 10.);

  if(ndigits == 1) {    
    for(s = 0; s < msa->nseq; s++) { 
      for(c = 0; c < msa->alen; c++) { 
	if(! esl_abc_CIsGap(msa->abc, msa->gr[ridx1][s][c])) {
	  switch(msa->gr[ridx1][s][c]) { 
	  case '*': p = 1.0; break;
	  case '9': p = 0.9; break;
	  case '8': p = 0.8; break;
	  case '7': p = 0.7; break;
	  case '6': p = 0.6; break;
	  case '5': p = 0.5; break;
	  case '4': p = 0.4; break;
	  case '3': p = 0.3; break;
	  case '2': p = 0.2; break;
	  case '1': p = 0.1; break;
	  case '0': p = 0.0; break;
	  default: 
	    ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, unrecognized residue: %c\n", s, c, msa->gr[ridx1][s][c]);
	  }
	  sum_c[c] += p;
	  sum_s[s] += p;
	  nongap_c[c]++;
	  nongap_s[s]++;
	  min_c[c] = ESL_MIN(min_c[c], p);
	  min_s[s] = ESL_MIN(min_s[s], p);
	  if(p >= pthresh) athresh_c[c]++;
	}
	else p = -1; /* gap */
	
	/*post_sc[s][c] = p;
	  post_cs[c][s] = p;*/
      }
    }
  }
  if(ndigits == 2) { 
    for(s = 0; s < msa->nseq; s++) { 
      for(c = 0; c < msa->alen; c++) { 
	if(! esl_abc_CIsGap(msa->abc, msa->gr[ridx1][s][c])) {
	  if(esl_abc_CIsGap(msa->abc, msa->gr[ridx2][s][c])) ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value non-gap but post 'ones' value is gap.\n", s, c);
	  if(msa->gr[ridx1][s][c] == '*') {
	    if(msa->gr[ridx2][s][c] != '*') ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value '*' but post 'ones' value != '*'.\n", s, c);
	    p = 1.0;
	  }
	  else {
	    ir1 = (int) (msa->gr[ridx1][s][c] - '0');
	    ir2 = (int) (msa->gr[ridx2][s][c] - '0');
	    p = ((float) ir1 * 10. + ir2) * .01;
	    /* printf("r1: %c %d r2: %c %d p: %.2f\n", msa->gr[ridx1][s][c], ir1, msa->gr[ridx2][s][c], ir2, p);*/
	  }
	  sum_c[c] += p;
	  sum_s[s] += p;
	  nongap_c[c]++;
	  nongap_s[s]++;
	  min_c[c] = ESL_MIN(min_c[c], p);
	  min_s[s] = ESL_MIN(min_s[s], p);
	  if(p >= pthresh) athresh_c[c]++;
	}
	else p = -1; /* gap */
	
	/*post_sc[s][c] = p;
	  post_cs[c][s] = p;*/
      }
    }
  }

  if(msa->rf != NULL) map_cpos_to_apos(msa, &c2a_map, &clen);
  else c2a_map = NULL;

  /* get averages */
  for(s = 0; s < msa->nseq; s++) { 
    avg_s[s]  =  (float) sum_s[s] / (float) nongap_s[s];
  }
  cpos = 1;
  for(c = 0; c < msa->alen; c++) { 
    avg_c[c]  = (float) sum_c[c] / (float) nongap_c[c];
    sum_total += sum_c[c];
    nongap_total += nongap_c[c];
    if(c2a_map != NULL) {
      if(c2a_map[cpos] == (c+1)) {  /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	cpos++; 
	sum_total_rf += sum_c[c];
	nongap_total_rf += nongap_c[c];
      }
    }
  }
  /* determine the fraction of sequences in each column with POSTs that exceed pthresh */
  for(c = 0; c < msa->alen; c++) 
    athresh_fract_c[c] = (nongap_c[c] > 0) ? ((float) athresh_c[c] / (float) nongap_c[c]) : 0;


  printf("\nAverage posterior value:                            %.5f (%d non-gap residues)\n", (float) sum_total / (float) nongap_total, nongap_total);
  if(c2a_map != NULL) 
    printf("Average posterior value in non-gap #=GC RF columns: %.5f (%d non-gap RF residues)\n", (float) sum_total_rf / (float) nongap_total_rf, nongap_total_rf);
  printf("\n");


  /* if nec, print posterior info */
  cpos = 1;
  if(do_pinfo) { 
    if ((pinfofp = fopen(esl_opt_GetString(go, "--pinfo"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --pinfo output file %s\n", esl_opt_GetString(go, "--pinfo"));
    fprintf(pinfofp, "# Posterior stats per column:\n");
    if(msa->rf != NULL) { 
      fprintf(pinfofp, "# %5s %5s %6s %6s %6s > %5.3f\n", "rfcol", "col", "nongap", "avg", "min", pthresh);
      fprintf(pinfofp, "# %5s %5s %6s %6s %6s %7s\n", "-----", "-----", "------", "------", "------", "-------");
      for(c = 0; c < msa->alen; c++) { 
	if(c2a_map[cpos] == (c+1)) { 
	  fprintf(pinfofp, "  %5d ", cpos);
	  cpos++; /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	}
	else fprintf(pinfofp, "  %5s ", "");
	if(nongap_c[c] == 0) fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)),      0.0,      0.0, athresh_fract_c[c]);
	else                 fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)), avg_c[c], min_c[c], athresh_fract_c[c]);
      }
    }
    else { /* msa->rf is NULL, we can't indicate the non-gap RF columns */
      fprintf(pinfofp, "%5s %6s %6s %6s > %5.3f\n", "col", "nongap", "avg", "min", pthresh);
      fprintf(pinfofp, "%5s %6s %6s %6s %7s\n", "-----", "------", "------", "------", "-------");
      for(c = 0; c < msa->alen; c++) 
	fprintf(pinfofp, "%5d %6.3f %6.3f %6.1f %7.3f\n", c+1, ((float) (nongap_c[c]) / ((float) msa->nseq)), avg_c[c], min_c[c], athresh_fract_c[c]);
    }
    fprintf(pinfofp, "\n\n");

    fprintf(pinfofp, "# Posterior stats per sequence:\n");
    fprintf(pinfofp, "# %5s %-60s %6s %6s %6s\n", "idx",   "seq name", "nongap", "avg", "min");
    fprintf(pinfofp, "# %5s %-60s %6s %6s %6s\n", "-----", "------------------------------------------------------------", "------", "------", "------");
    for(s = 0; s < msa->nseq; s++) { 
      fprintf(pinfofp, "  %5d %-60s %6.3f %6.3f %6.2f\n", s+1, msa->sqname[s], ((float) (nongap_s[s]) / ((float) msa->alen)), avg_s[s], min_s[s]); 
    }
    fclose(pinfofp);
  }

  /* optionally, add/rewrite msa->rf if --pfract enabled */
  if(do_pfract) { 
    pfract = esl_opt_GetReal(go, "--pfract");  
    ESL_ALLOC(useme, sizeof(int) * (msa->alen+1));
    if(do_prf) { /* only look at consensus columns */
      cpos = 1;
      for(c = 0; c < msa->alen; c++) { 
	if(c2a_map[cpos] == (c+1)) { /* off-by-one, c2a_map is 1..clen, c is 0..alen */
	  cpos++;
	  if(athresh_fract_c[c] >= pfract) useme[c] = 1; 
	  else                             useme[c] = 0; 
	}
	else useme[c] = 0; 
      }    
    }
    else { /* look at all columns */
      for(c = 0; c < msa->alen; c++) {  
	if(athresh_fract_c[c] >= pfract) useme[c] = 1; 
	else                             useme[c] = 0; 
      }
    }
    useme[msa->alen] = '\0';

    nkept = 0;
    write_rf_given_useme(go, errbuf, msa, useme);
    for(c = 0; c < msa->alen; c++) nkept += useme[c];
    free(useme);
    if(do_prf)  printf("\n%d of %d RF columns (%.3f) pass threshold\n\n", nkept, clen, (float) nkept / (float) clen);
    else        printf("\n%d of %" PRId64 " columns (%.3f) pass threshold\n\n", nkept, msa->alen, (float) nkept / (float) msa->alen);
  }

  /*for(s = 0; s < msa->nseq; s++) free(post_sc[s]);
    free(post_sc);
    for(c = 0; c < msa->alen; c++) free(post_cs[c]);
    free(post_cs);*/
  
  free(athresh_fract_c);
  free(athresh_c);
  free(nongap_s);
  free(nongap_c);
  free(min_c);
  free(min_s);
  free(avg_c);
  free(avg_s);
  free(sum_s);
  free(sum_c);
  if(c2a_map != NULL) free(c2a_map);

  return eslOK;

 ERROR:
  return status;
}


/* output_rf_as_mask
 *
 * Given an MSA with rf annotation, convert it to a lanemask of 1s and 0s.
 * 1s for non-gap RF columns, 0s for gap RF columns.
 */
static int
output_rf_as_mask(FILE *fp, char *errbuf, ESL_MSA *msa)
{
  int status;
  int  apos;
  char *mask;

  if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "msa->rf is NULL, and we're trying to convert it to a 1/0 mask.");
  ESL_ALLOC(mask, sizeof(char) * (msa->alen+1));

  for (apos = 0; apos < msa->alen; apos++) 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos]))  mask[apos] = '0';
    else                                         mask[apos] = '1';
  mask[msa->alen] = '\0';

  fprintf(fp, "%s\n", mask);
  free(mask);

  return eslOK;
 ERROR:
  return status;
}


/* expand_msa2mask
 *
 * Given an MSA <msa> and a lanemask <xmask> with exactly msa->alen 1s in it.
 * Add 100% gap columns in between each column as dictated by <xmask>.
 *
 * For example if lanemask is 100101, msa->alen is 3, we add 2 100% gap
 * columns after column 1, and 1 100% gap column after column 2, to make
 * the msa length = length(xmask) = 6.
 */
static int
expand_msa2mask(char *errbuf, ESL_MSA *msa, char *xmask, ESL_MSA **newmsa)
{
  int status;
  int  mpos;
  int  masklen;
  int *nzeroesA;
  int  nones = 0;

  if(xmask == NULL) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), xmask is NULL.");

  masklen = strlen(xmask);
  /* count 1s in xmask */
  for (mpos = 0; mpos < masklen; mpos++) { 
    if     (xmask[mpos] == '1') nones++;
    else if(xmask[mpos] == '0') ; /* do nothing */
    else    ESL_FAIL(eslEINVAL, errbuf, "--xmask mask char number %d is not a 1 nor a 0, but a %c\n", mpos+1, xmask[mpos]);
  }
  if(nones != msa->alen) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), number of 1s in --xmask file: %d != msa->alen: %" PRId64 ", they must be equal.", nones, msa->alen);

  /* determine number of 0s after each consensus column */
  nones = 0;
  ESL_ALLOC(nzeroesA, sizeof(int) * masklen+1);
  esl_vec_ISet(nzeroesA, (masklen+1), 0);
  for (mpos = 0; mpos < masklen; mpos++) { 
    if     (xmask[mpos] == '1') nones++;
    else if(xmask[mpos] == '0') nzeroesA[nones]++;
    else    ESL_FAIL(eslEINVAL, errbuf, "--xmask mask char number %d is not a 1 nor a 0, but a %c\n", mpos+1, xmask[mpos]);
  }
  
  /*int i;
  for (i = 0; i <= nones; i++) { 
    printf("nzeroes[%3d]: %3d\n", i, nzeroesA[i]);
    }*/

  /* add the 100% gap columns */
  if((status = add_gap_columns_to_msa(errbuf, msa, nzeroesA, newmsa, TRUE)) != eslOK) return status ;
  /* new alen should equal masklen */
  if((*newmsa)->alen != masklen) ESL_FAIL(eslEINVAL, errbuf, "expand_msa2mask(), new msa->alen: (%" PRId64 ") != length of mask (%d), this shouldn't happen.", (*newmsa)->alen, masklen);
  free(nzeroesA);

  return eslOK;
 ERROR:
  return status;
}

/* Function: compare_ints()
 * 
 * Purpose:  Comparison function for qsort(). Used 
 *           by msa_median_length().
 */ 
static int 
compare_ints(const void *el1, const void *el2)
{
  if      ((* ((int *) el1)) > (* ((int *) el2)))  return 1;
  else if ((* ((int *) el1)) < (* ((int *) el2)))  return 1;
  return 0;
}

/* Function: msa_median_length()
 * 
 * Purpose:  Returns the median (unaligned) length of 
 *           the sequences in an alignment.
 */
static int
msa_median_length(ESL_MSA *msa)
{
  int  status;
  int *len;
  int  i;
  int  median;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(len, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) {
    esl_sq_GetFromMSA(msa, i, sq);
    len[i] = sq->n;
    esl_sq_Reuse(sq);
    /*printf("i: %d len: %d\n", i, len[i]);*/
  }

  qsort(len, msa->nseq, sizeof(int), compare_ints);

  median = len[msa->nseq / 2];
  free(len);

  esl_sq_Destroy(sq);
  return median;

 ERROR:
  esl_fatal("msa_median_length() memory allocation error.");
  return 0.; /* NEVERREACHED */
}

/* Function: msa_remove_seqs_below_minlen()
 * 
 * Purpose:  Remove sequences in MSA whose dealigned length is less than a minimum length.
 */
static int
msa_remove_seqs_below_minlen(ESL_MSA *msa, float minlen, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i;

  ESL_MSA *new_msa;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) {
    esl_sq_GetFromMSA(msa, i, sq);
    useme[i] = ((float) sq->n >= minlen) ? TRUE : FALSE;
    /*printf("useme[i:%d]: %d\n", i, useme[i]);*/
    esl_sq_Reuse(sq);
  }

  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  free(useme);
  esl_sq_Destroy(sq);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  esl_fatal("msa_remove_seqs_below_minlen() memory allocation error.");
  return eslOK; /* NEVERREACHED */
}

/* Function: msa_remove_truncated_seqs()
 * 
 * Purpose:  Remove sequences in MSA that have all gaps in the first <ntrunc> 5' leading 
 *           non-gap RF columns OR the last <ntrunc> 3' leading non-gap RF columns
 */
static int
msa_remove_truncated_seqs(ESL_MSA *msa, char *errbuf, int ntrunc, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i;
  int  leading_okay, trailing_okay;
  int  apos, cpos_ct;
  int  nused = 0;
  ESL_MSA *new_msa;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in msa_remove_truncated_seqs(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --detrunc.");

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);

  for(i = 0; i < msa->nseq; i++) { 
    /* if ALL of the first 5' <ntrunc> non-gap RF columns are gaps in this seq, we'll remove it */
    leading_okay  = FALSE;
    cpos_ct = 0; 
    apos = 1;
    while(!leading_okay && (cpos_ct < ntrunc) && (apos <= msa->alen)) { 
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos_ct++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) leading_okay = TRUE;
      }
      apos++;
    }

    trailing_okay = FALSE;
    cpos_ct = 0;
    apos = msa->alen;
    while(!trailing_okay && (cpos_ct < ntrunc) && (apos >= 1)) { 
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos_ct++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) trailing_okay = TRUE;
      }
      apos--;
    }
    useme[i] = (leading_okay && trailing_okay) ? TRUE : FALSE;
    if(useme[i]) nused++;
  }
  if(nused == 0) ESL_FAIL(eslEINVAL, errbuf, "--detrunc removed ALL sequences!");
  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  free(useme);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "msa_remove_truncated_seqs(): memory allocation error.");
  return eslOK; /* NEVERREACHED */
}

/* dump_infocontent
 *                   
 * Given an MSA with RF annotation, print a postscript heatmap of the 
 * information content of each non-gap RF column (consensus column).
 */
static int dump_infocontent(FILE *fp, ESL_MSA *msa, char *errbuf)
{
  int status;
  int apos, cpos;
  int i;
  int clen;
  double *obs, *ent, *bg;

  /* contract check */
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --icinfo.");
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in dump_infocontent(), msa must be digitized.");

  clen = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  ESL_ALLOC(ent, sizeof(double) * clen);
  ESL_ALLOC(obs, sizeof(double) * msa->abc->K);
  ESL_ALLOC(bg, sizeof(double) * msa->abc->K);
  esl_vec_DSet(bg, msa->abc->K, 1./(msa->abc->K));

  cpos = 0;
  fprintf(fp, "# %4s  %5s\n", "cpos", "info");
  fprintf(fp, "# %4s  %5s\n", "----", "-----");
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { /* a consensus position */
      esl_vec_DSet(obs, msa->abc->K, 0.);
      for(i = 0; i < msa->nseq; i++)
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	  esl_abc_DCount(msa->abc, obs, msa->ax[i][apos], 1.);
	}
      esl_vec_DNorm(obs, msa->abc->K);
      ent[cpos] = esl_vec_DEntropy(bg, msa->abc->K) - esl_vec_DEntropy(obs, msa->abc->K);
      fprintf(fp, " %4d  %5.3f\n", cpos, ent[cpos]);
      cpos++;
    }
  }

  free(ent);
  free(obs);
  free(bg);
  return eslOK;

 ERROR:
  return status;
}

/* number_columns
 *                   
 * Add annotation to an MSA numbering the columns, either all
 * the columns (if <do_all>) or just non-gap #=GC RF columns.
 */
static int
number_columns(ESL_MSA *msa, int do_all, char *errbuf)
{
  int  status;
  int i;
  char *numstring;
  char *tag;
  int alen_ndigits;
  int tagwidth;
  int a,b,apos;
  int bmin;
  int pos2print;
  /* contract check */
  if(!do_all && msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment.");

  alen_ndigits = int_ndigits(msa->alen);
  tagwidth = do_all ? (3+alen_ndigits) : (5+alen_ndigits); /* "COL.X" or RFCOL.X" */

  ESL_ALLOC(tag, sizeof(char) * (tagwidth+1));
  ESL_ALLOC(numstring, sizeof(char) * (msa->alen+1));
  numstring[msa->alen] = '\0';
  tag[tagwidth] = '\0';
  if(do_all) { 
    bmin = 3;
    tag[0] = 'C';
    tag[1] = 'O';
    tag[2] = 'L';
  }
  else { 
    bmin = 5;
    tag[0] = 'R';
    tag[1] = 'F';
    tag[2] = 'C';
    tag[3] = 'O';
    tag[4] = 'L';
  }

  for(a = 0; a < alen_ndigits; a++) { 
    for(b = 0; b < alen_ndigits; b++) tag[b+bmin] = (a == b) ? 'X' : '.';
    pos2print = 1;
    for(apos = 1; apos <= msa->alen; apos++) { 
      if(!do_all && (esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)]))) numstring[(apos-1)] = '.';
      else numstring[(apos-1)] = get_char_digit_x_from_int(pos2print++, (alen_ndigits-a));
	/*printf("called get_char_digit_x_from_int(%d, %d)\n",apos, (alen_ndigits-a));*/
    }
    esl_msa_AppendGC(msa, tag, numstring);
  }

  ESL_ALLOC(numstring, sizeof(char) * (msa->alen + 1));
  for(i = 0; i < msa->alen; i++) { 
    numstring[i] = digit_to_char(i);
  }
  numstring[msa->alen] = '\0';
  free(numstring);
  return eslOK;

 ERROR:
  return eslEMEM;
}


/* digit_to_char
 *                   
 * Given a digit (0-9) return the character reprentation of it.
 * There must be a better way to do this; oh well.
 */
static char
digit_to_char(int digit) 
{
  if(digit == 0) return '0';
  if(digit == 1) return '1';
  if(digit == 2) return '2';
  if(digit == 3) return '3';
  if(digit == 4) return '4';
  if(digit == 5) return '5';
  if(digit == 6) return '6';
  if(digit == 7) return '7';
  if(digit == 8) return '8';
  if(digit == 9) return '9';
  else return '?';
}

/* Function: int_ndigits
 * Returns: The number of digits in <i>.
 */
static int
int_ndigits(int i)
{
  int n   = 0;
  while(i > 0) { i/=10; n++; }
  return n;
}

/* get_char_digit_x_from_int
 *                   
 * Given two integers <i> and <place> return the 
 * character version of the <place>'th digit in <i>.
 * Example <i> = 14378 <place> = 4 would return 7.
 */
static char
get_char_digit_x_from_int(int i, int place)
{
  int n,a,divisor;
  n = int_ndigits(i);

  if(n < place) return digit_to_char(0);

  divisor = 1;
  for(a = 0; a < (place-1); a++) divisor *= 10;
  /* subtract leading digits before the one we care about */
  i %= (divisor*10);
  return digit_to_char (i / divisor);
}

/* Function: read_seq_name_file
 * Date:     EPN, Thu Jun  5 13:21:36 2008
 * 
 * Read a file listing sequence names to remove or keep.
 * Store sequences in *ret_seqlist and return it.
 * Each white-space delimited token is considered a 
 * different sequence name. No checking is done in this 
 * function, but rather in subsequent functions. Each sequence name is 
 * 
 * Returns eslOK on success.
 */
int
read_seq_name_file(char *filename, char *errbuf, char ***ret_seqlist, int *ret_seqlist_n)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  int nalloc     = 10;
  int chunksize  = 10;
  char **seqlist = NULL;
  int n = 0;
  int i;
  void *tmp;

  ESL_ALLOC(seqlist, sizeof(char *) * nalloc);
  if (esl_fileparser_Open(filename, NULL,  &efp) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "failed to open %s in read_seq_name_file\n", filename);
  
  while((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslEOF) {
    if(n == nalloc) { nalloc += chunksize; ESL_RALLOC(seqlist, tmp, sizeof(char *) * nalloc); }
    if((status = esl_strdup(tok, -1, &(seqlist[n++]))) != eslOK) ESL_FAIL(status, errbuf, "error in esl_strdup.");
  }
  esl_fileparser_Close(efp);
  *ret_seqlist = seqlist;
  *ret_seqlist_n = n;
  return eslOK;

 ERROR:
  if(seqlist != NULL) {
    for(i = 0; i < n; i++) free(seqlist[i]); 
    free(seqlist);
  }
  return status;
}


/* Function: msa_keep_or_remove_seqs()
 * 
 * Purpose:  Given a list of <seqlist_n> sequences in <seqlist>, either remove those
 *           sequences from msa, or remove all other sequences besides those from msa.
 *           Create and return the new msa with only the specified seqs in <ret_new_msa>.
 * 
 * Note:     Terribly inefficient, does a linear search for each seq, no sorting or anything.
 *
 * Returns: eslOK on success, eslEINVAL if a sequence name in seqlist does not exist in the msa.
 * 
 */
static int
msa_keep_or_remove_seqs(ESL_MSA *msa, char *errbuf, char **seqlist, int seqlist_n, int do_keep, ESL_MSA **ret_new_msa)
{
  int  status;
  int *useme;
  int  i, ip, n;
  int *order_all, *order_new;

  ESL_MSA *new_msa;
  ESL_SQ *sq;
  sq = esl_sq_CreateDigital(msa->abc);

  ESL_ALLOC(useme,     sizeof(int) * msa->nseq);
  ESL_ALLOC(order_all, sizeof(int) * msa->nseq);
  ESL_ALLOC(order_new, sizeof(int) * seqlist_n);
  esl_vec_ISet(order_all, msa->nseq, -1);
  if(do_keep) esl_vec_ISet(useme, msa->nseq, FALSE);
  else        esl_vec_ISet(useme, msa->nseq, TRUE); 
  for(n = 0; n < seqlist_n; n++) { 
    for (i = 0; i < msa->nseq; i++) {
      if(strcmp(seqlist[n], msa->sqname[i]) == 0) { 
	useme[i] = do_keep ? TRUE : FALSE;
	order_all[i] = n;
	break;
      }
      if(i == (msa->nseq-1)) ESL_FAIL(eslEINVAL, errbuf, "ERROR sequence %s does not exist in the MSA!", seqlist[n]);
    }
  }      

  if((status = esl_msa_SequenceSubset(msa, useme, &new_msa)) != eslOK) esl_fatal("esl_msa_SequenceSubset() had a problem.");
  /* if do_keep, reorder to order of names in the list file */
  if(do_keep) { 
    ip = 0;
    for(i = 0; i < msa->nseq; i++) { 
      if(order_all[i] != -1) order_new[order_all[i]] = ip++;
    }
    if((status = reorder_msa(new_msa, order_new, errbuf)) != eslOK) return status;
  }

  free(useme);
  free(order_all);
  free(order_new);
  *ret_new_msa = new_msa;
  return eslOK;

 ERROR:
  esl_fatal("msa_keep_or_remove_seqs() memory allocation error.");
  return eslOK; /* NEVERREACHED */
}


/* Function:  insert_x_pair_shared()
 * Synopsis:  Calculate the fraction of inserts shared between of two aligned digital seqs.
 * Incept:    EPN, Wed Jun 25 10:33:23 2008
 *
 * Purpose:   Returns the fraction of the total
 *            number of inserts in both sequences that are shared.
 *            An 'insert' is present in sequence s for consensus column 
 *            (non-gap RF column) c if at least 1 residue exists between 
 *            consensus column c and c+1. If sequence t also has an insert
 *            between c and c+1, they share that insert. If that were the
 *            only insert in either of the two sequences, then they would
 *            share 1.0 fraction of inserts.
 *            
 * Args:      msa          - MSA, digitized, with RF annotation
 *            i            - index of seq 1
 *            j            - indes of seq 2
 *            cfirst       - first consensus position to consider
 *            clast        - last consensus position to consider
 *            opt_pshared  - optRETURN: pairwise insert identity, 0<=x<=1
 *            opt_nshared  - optRETURN: # of inserts shared
 *            opt_nins     - optRETURN: nins
 *
 * Returns:   <eslOK> on success. <opt_distance>, <opt_nid>, <opt_n>
 *            contain the answers, for any of these that were passed
 *            non-<NULL> pointers.
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 */
int
insert_x_pair_shared(ESL_MSA *msa, int i, int j, int cfirst, int clast, double *opt_pshared, int *opt_nshared, int *opt_nins)
{
  int     shared;               /* total shared inserts */
  int     nins;                 /* number of inserts in either sequence */
  int     apos;                 /* position in aligned seqs   */
  int     cpos;
  int     insi, insj;
  int     seen_insert = FALSE;
  shared = nins = 0;
  
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos++;
	seen_insert = FALSE;
      }
      else { /* not a consensus column, an insert column */
	insi = (! esl_abc_XIsGap(msa->abc, msa->ax[i][apos]));
	insj = (! esl_abc_XIsGap(msa->abc, msa->ax[j][apos]));
	if(cpos >= cfirst && cpos <= clast) { 
	  if(insi && insj   && !seen_insert)   shared++;
	  if((insi || insj) && !seen_insert) { nins++; seen_insert = TRUE; }
	}
      }
    }
  /*if (opt_pshared  != NULL)  *opt_pshared  = ( nins==0 ? 0. : (double) shared / (double) nins );*/
  if (opt_pshared  != NULL)  *opt_pshared  = ( nins==0 ? 1. : (double) shared / (double) nins );
  if (opt_nshared  != NULL)  *opt_nshared  = shared;
  if (opt_nins     != NULL)  *opt_nins     = nins;
  return eslOK;
}


/* Function:  insert_x_pair_shared_length()
 * Synopsis:  Calculate the fraction of inserts shared between of two aligned digital seqs,
 *            weighted by the length of the inserts.
 * Incept:    EPN, Wed Jun 25 10:33:23 2008
 *
 * Purpose:   Returns the weighted fraction of the total
 *            number of inserts in both sequences that are shared.
 *            
 * Args:      msa          - MSA, digitized, with RF annotation
 *            i            - index of seq 1
 *            j            - indes of seq 2
 *            cfirst       - first consensus position to consider
 *            clast        - last consensus position to consider
 *            opt_pshared  - optRETURN: pairwise insert identity, 0<=x<=1
 *            opt_nshared  - optRETURN: weighted # inserts shared
 *            opt_nins     - optRETURN: nins, number of columns with an insert
 *
 * Returns:   <eslOK> on success. <opt_distance>, <opt_nid>, <opt_n>
 *            contain the answers, for any of these that were passed
 *            non-<NULL> pointers.
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 */
int
insert_x_pair_shared_length(ESL_MSA *msa, int i, int j, int cfirst, int clast, double *opt_pshared, double *opt_nshared, int *opt_nins)
{
  double  shared;               /* weighted shared inserts */
  int     nins;                 /* number of inserts in either sequence */
  int     apos;                 /* position in aligned seqs   */
  int     cpos;
  int     leni, lenj;
  shared = nins = leni = lenj = 0;
  
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { 
	cpos++;
	if(cpos >= cfirst && cpos <= clast) { 
	  if((leni + lenj) > 0) { 
	    nins++;
	    if(leni >= lenj) shared += (double) lenj / (double) leni;
	    else             shared += (double) leni / (double) lenj;
	  }
	  leni = lenj = 0;
	}
      }
      else { /* not a consensus column, an insert column */
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) leni++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[j][apos])) lenj++;
      }
    }
  /*if (opt_pshared  != NULL)  *opt_pshared  = ( nins==0 ? 0. : (double) shared / (double) nins );*/
  if (opt_pshared  != NULL)  *opt_pshared  = ( nins==0 ? 1. : (double) shared / (double) nins );
  if (opt_nshared  != NULL)  *opt_nshared  = shared;
  if (opt_nins     != NULL)  *opt_nins     = nins;
  return eslOK;
}

/* Function:  insert_x_diffmx()
 * Synopsis:  NxN insert difference matrix for N aligned digital seqs.         
 * Incept:    EPN, Wed Jun 25 10:25:01 2008
 *
 * Purpose:   For each pair of sequences calculates the fraction of number
 *            of inserts that are different between the two sequences.
 *            An 'insert' is present in sequence s for consensus column 
 *            (non-gap RF column) c if at least 1 residue exists between 
 *            consensus column c and c+1. If sequence t also has an insert
 *            between c and c+1, they share that insert. If that were the
 *            only insert in either of the two sequences, then they would
 *            share 1.0 fractional insert id, and 1.0 - 1.0 = 0.0 fractional 
 *            insert difference, thus the insert diff mx entry between 
 *            seq s and t would be 0.0.
 *
 * Args:      go - command-line options
 *            errbuf - for printing error messages
 *            msa   - aligned dsq's, [0..N-1][1..alen]                  
 *            do_length_weight - weight insert similarity by length of inserts
 *            do_only_internal_inserts - TRUE to only count inserts that are at positions
 *                                       internal to both seq i, j (don't count those truncated off in either seq)
 *            ret_D - RETURN: NxN matrix of fractional differences
 *            
 * Returns:   <eslOK> on success, and <ret_D> contains the difference
 *            matrix; caller is obligated to free <D> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
insert_x_diffmx(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int do_length_weight, int do_only_internal_inserts, ESL_DMATRIX **ret_D)
{
  int status;
  ESL_DMATRIX *D = NULL;
  int i,j;
  int N = msa->nseq;
  int nshared;
  double nshared_len;
  int nins;
  int *firstA, *lastA; /* [0..i..nseq-1] first and last consensus column occupied by seq i, only used if do_only_internal_inserts == TRUE */
  int  clen;
  int ifirst, ilast, jfirst, jlast;

  if(msa->rf == NULL)                  ESL_FAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment.");
  if(! (msa->flags & eslMSA_DIGITAL))  ESL_FAIL(eslEINVAL, errbuf, "insert_x_diffmx() MSA is not digitized.\n");

  if (( D = esl_dmatrix_Create(N,N) ) == NULL) goto ERROR;
  if ((status = determine_first_last_consensus_columns(msa, errbuf, &firstA, &lastA, &clen)) != eslOK) return status;

  /* TEMP  for (i = 0; i < N; i++) printf("i: %4d %4d %4d\n", i, firstA[i], lastA[i]); */

  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      ifirst = do_only_internal_inserts ? firstA[i] : 0;
      ilast  = do_only_internal_inserts ? lastA[i]  : clen;
      for (j = i+1; j < N; j++)
	{
	  jfirst = do_only_internal_inserts ? firstA[j] : 0;
	  jlast  = do_only_internal_inserts ? lastA[j]  : clen;
	  if(do_length_weight) { 
	    status = insert_x_pair_shared_length(msa, i, j, ESL_MAX(ifirst, jfirst), ESL_MIN(ilast, jlast), &(D->mx[i][j]), &nshared_len, &nins);
	    if(esl_opt_GetBoolean(go, "--verbose")) printf("D %4d %4d %.3f %8.3f of %4d\n", i, j, 1. - D->mx[i][j], nshared_len, nins);
	  }
	  else { 
	    status = insert_x_pair_shared(msa, i, j, ESL_MAX(ifirst, jfirst), ESL_MIN(ilast, jlast), &(D->mx[i][j]), &nshared, &nins);
	    if(esl_opt_GetBoolean(go, "--verbose")) printf("D %4d %4d %.3f %4d of %4d\n", i, j, 1. - D->mx[i][j], nshared, nins);
	  }
	  D->mx[i][j] = 1. - D->mx[i][j]; /* convert from id to distance */
	  if (status != eslOK) ESL_XEXCEPTION(status, "Pairwise insert identity calculation failed at seqs %d,%d\n", i,j);
	  D->mx[j][i] =  D->mx[i][j];
	}
      if(esl_opt_GetBoolean(go, "--verbose")) printf("\n");
    }
  if (ret_D != NULL) *ret_D = D; else esl_dmatrix_Destroy(D);
  return eslOK;

 ERROR:
  if (D     != NULL)  esl_dmatrix_Destroy(D);
  if (ret_D != NULL) *ret_D = NULL;
  return status;
}

/* Function: MSADivide()
 * From Infernal's cmbuild.c
 * EPN, Wed Mar 21 17:26:39 2007
 * 
 * Purpose:  Given an MSA and a distance matrix, divide the MSA it into 
 *           multiple MSAs, each with a different cluster of the original 
 *           sequences. Where clusters are defined by single linkage
 *           clustering based on the distance matrix.
 *
 *           Different modes:
 *           
 *        1. if(do_mindiff): define clusters
 *           such that we maximize the number of clusters while
 *           satisfying: minimum fractional difference b/t any 
 *           2 seqs in different clusters >= 'mindiff'. 
 *           The contract states that mindiff > 0. in this case.
 *           
 *        2. if(do_nc): define clusters 
 *           such that we have exactly 'target_nc' clusters by
 *           searching for the 'mindiff' that gives exactly
 *           'target_nc' clusters. (We guarantee we can do this
 *           by rounding 'diff' fractional difference values b/t
 *           seqs to nearest 0.001). 
 *
 *        3. if(do_nsize): define clusters 
 *           such that we have 1 cluster that has at least nsize
 *           sequences in it by searching for the 'mindiff' that 
 *           achieves that.
 *
 * Args:    
 * ESL_MSA *mmsa        - the master MSA, we cluster the seqs in this guy
 *                        and build a new MSA from each cluster
 * ESL_DMATRIX *D;      - the distance matrix
 * int     do_mindiff   - TRUE (mode 1): satisfy clusters are at least mindiff different
 * int     do_nc        - TRUE (mode 2): set mindiff such that we get exactly target_nc clusters
 * int     do_nsize     - TRUE (mode 3): set mindiff such that we get 1 cluster with nsize seqs
 * float   mindiff      - the minimum fractional difference allowed between 2 seqs of different clusters
 *                        (0. indicates mode 2 or 3) 
 * int     target_nc    - if(do_nc) number of clusters to define, else irrelevant
 * int     target_nsize - if(do_nsize) min size of largets cluster, else irrelevant
 * int     *ret_num_msa - number of MSAs in ret_MSA
 * ESL_MSA  ***ret_cmsa - new MSAs, one for each cluster
 * ESL_MSA  *ret_xsize  - max size of a cluster
 * char     *errbuf     - buffer for error messages
 *           
 * Return: ret_cmsa (alloc'ed here) and ret_num_msa
 */
int 
MSADivide(ESL_MSA *mmsa, ESL_DMATRIX *D, int do_mindiff, int do_nc, int do_nsize, float mindiff, int target_nc,
	  int target_nsize, int *ret_num_msa, ESL_MSA ***ret_cmsa, int *ret_xsize, char *errbuf)
{
  int   status;        /* Easel status code */
  ESL_MSA **cmsa = NULL;/* the new MSAs we're creating from clusters of seqs in mmsa */
  int   i;             /* counter over sequences */
  int   m;             /* counter over new MSAs */
  int   n;             /* counter over tree nodes */
  ESL_TREE    *T = NULL;/* the tree, created by Single-Linkage Clustering */
  double *diff = NULL; /* [0..T->N-2], diff[n]= min distance between any leaf in right and
		        * left subtree of node n of tree T */
  double *minld = NULL;/* [0..T->N-2], min dist from node to any taxa in left  subtree */
  double *minrd = NULL;/* [0..T->N-2], min dist from node to any taxa in right subtree */
  int     nc;          /* number of clusters/MSAs  */
  int    *clust = NULL;/* [0..T->N-1], cluster number (0..nc-1) this seq is in */
  int    *csize = NULL;/* [0..nc-1], size of each cluster */
  int   **useme = NULL;/* [0.m.nc-1][0.i.N] TRUE to use seq i in new MSA m, FALSE not to */
  int     best;        /* 'best' node, returned by select_node() */
  int     xsize;       /* size of cluster under 'best' node (largest cluster) */

  /* Contract check */
  if((do_nc + do_mindiff + do_nsize) != 1) ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() exactly 1 of do_nc, do_mindiff, do_nsize must be TRUE.");
  if( do_nc && target_nc == 0)             ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() target_nc is 0 but do_nc is TRUE!");
  if( do_nsize && target_nsize == 0)       ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() target_nsize is 0 but do_nsize is TRUE!");
  if( do_mindiff && mindiff <= 0.)         ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() mindiff is <= 0. but do_mindiff is TRUE!");
  if( do_mindiff && target_nc != 0)        ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() do_mindiff is TRUE, but target_nc != 0");
  if( do_mindiff && target_nsize != 0)     ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() do_mindiff is TRUE, but target_nsize != 0");
  /* mmsa must be digital */
  if(!(mmsa->flags & eslMSA_DIGITAL))                 ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() MSA is not digital.");

  if(do_nc || do_nsize) mindiff = 0.;

  /* Infer tree by single linkage clustering given the distance matrix */
  if((status = esl_tree_SingleLinkage(D, &T)) != eslOK)                        ESL_FAIL(status, errbuf, "esl_tree_SingleLinkage() error, status: %d", status);
  if((status = esl_tree_SetTaxaParents(T)) != eslOK)                           ESL_FAIL(status, errbuf, "esl_tree_SetTaxaParentse() error, status: %d", status);
  /*esl_tree_WriteNewick(stdout, T);*/
  if((status = esl_tree_Validate(T, errbuf) != eslOK)) return status;
  
  /* determine the diff values: 
   * (use: n_child > n, unless n's children are taxa)
   * diff[n] is minimum distance between any taxa (leaf) in left subtree of 
   * n to any taxa in right subtree of n. 
   */
  ESL_ALLOC(diff,  (sizeof(double) * (T->N - 1)));  /* one for each node */
  ESL_ALLOC(minld, (sizeof(double) * (T->N - 1))); 
  ESL_ALLOC(minrd, (sizeof(double) * (T->N - 1))); 
  for (n = (T->N-2); n >= 0; n--) {
    minld[n] = T->ld[n] + ((T->left[n]  > 0) ? (minld[T->left[n]])  : 0);
    minrd[n] = T->rd[n] + ((T->right[n] > 0) ? (minrd[T->right[n]]) : 0);
    diff[n]  = minld[n] + minrd[n];
    diff[n] *= 1000.; 
    diff[n]  = (float) ((int) diff[n]);
    diff[n] /= 1000.; 
    /*printf("diff[n:%d]: %f\n", n, diff[n]);*/
  }
  free(minld); minld = NULL;
  free(minrd); minrd = NULL;
  /*for (n = 0; n < (T->N-1); n++)
    printf("diff[n:%d]: %f\n", n, diff[n]);
    for (n = 0; n < (T->N-1); n++)
    printf("left[n:%d]: %d right[n:%d]: %d\n", n, T->left[n], n, T->right[n]);*/
  
  if(do_mindiff) { /* Mode 1 */
    /* Define clusters that are at least mindiff different
     * from each other. */
    if((status = select_node(T, diff, mindiff, &clust, &nc, &xsize, &best, errbuf)) != eslOK) return status;
    printf("# Alignment split into %d clusters\n", nc);
    printf("# Maximum identity b/t any 2 seqs in different clusters: %.2f\n", (1.-mindiff));
    printf("# Largest cluster contains %d sequences.\n", xsize);
    printf("#\n");
  }
  else if (do_nc) { /* Mode 2, do_nc == TRUE, mindiff was set to 0.0 above */
    /* Find the minimum fractional difference (mindiff) that 
     * gives exactly target_nc clusters, also define clusters
     * based on that mindiff, this is all done with find_mindiff(),
     * which does a binary search for mindiff, we're guaranteed to 
     * find exactly target_nc clusters b/c diff values are rounded
     * to nearest 0.001. */
    if(target_nc > (T->N)) target_nc = T->N; /* max num clusters is num seqs */
    if((status = find_mindiff(T, diff, FALSE, target_nc, &clust, &nc, &xsize, &best, &mindiff, errbuf)) != eslOK) return status;
    printf("# Alignment split into %d clusters.\n", nc);
    printf("# Maximum identity b/t any 2 seqs in different clusters: %.2f\n", (1.-mindiff));
    printf("# Largest cluster contains %d sequences.\n", xsize);
    printf("#\n");
  }
  else { /* Mode 3, do_nsize == TRUE, mindiff was set to 0.0 above */
    /* Find the minimum fractional difference (mindiff) that 
     * gives 1 cluster with size of at least <target_nsize> sequences
     * based on that mindiff, this is all done with find_mindiff(),
     * which does a binary search for mindiff.
     */
    if(target_nsize > (T->N)) target_nsize = T->N; /* max num clusters is num seqs */
    if((status = find_mindiff(T, diff, TRUE, target_nsize, &clust, &nc, &xsize, &best, &mindiff, errbuf)) != eslOK) return status;
    printf("# Alignment split into %d clusters.\n", nc);
    printf("# Largets cluster contains %d sequences.\n", xsize);
    printf("# Maximum identity b/t any 2 seqs in different clusters: %.2f\n", (1.-mindiff));
    printf("#\n");
  }
  /* Determine the size of each cluster */
  ESL_ALLOC(csize, (sizeof(int) * (nc)));
  esl_vec_ISet(csize, nc, 0);
  for(i = 0; i < mmsa->nseq; i++)
    csize[clust[i]]++;
  
  /* Create one new MSA for each cluster,
   * if(do_orig): keep the original MSA as cmsa[nc] */
  ESL_ALLOC(cmsa, (sizeof(ESL_MSA *) * (nc)));
  for(m = 0; m < nc; m++) cmsa[m] = NULL;

  ESL_ALLOC(useme, (sizeof(int *) * (nc+1)));
  for(m = 0; m <= nc; m++) {
    ESL_ALLOC(useme[m], (sizeof(int)) * mmsa->nseq);
    if(m < nc) esl_vec_ISet(useme[m], mmsa->nseq, FALSE);
    else       esl_vec_ISet(useme[m], mmsa->nseq, TRUE); /* keep all seqs in cmsa[nc]*/
  }
  
  for(i = 0; i < mmsa->nseq; i++)
    if(clust[i] != -1) 
      useme[clust[i]][i] = TRUE;
  printf("#   idx    nseq\n");
  printf("#  ----  ------\n");
  for(m = 0; m < nc; m++) {
    if((status = esl_msa_SequenceSubset(mmsa, useme[m], &(cmsa[m]))) != eslOK) ESL_FAIL(status, errbuf, "MSADivide(), esl_msa_SequenceSubset error, status: %d.", status);
    printf("   %4d  %6d\n", m+1, cmsa[m]->nseq);
    free(useme[m]);
  }
  printf("\n");

  free(useme[nc]);
  free(useme);
  
  *ret_num_msa = nc;
  *ret_cmsa = cmsa;
  *ret_xsize = xsize;
  
  esl_tree_Destroy(T);
  free(diff);
  diff = NULL;

  if(clust != NULL)  free(clust);
  if(csize != NULL)  free(csize);
  return eslOK;
  
 ERROR: 
  if(diff  != NULL) free(diff);
  if(minld != NULL) free(minld);
  if(minrd != NULL) free(minrd);
  if(clust != NULL) free(clust);
  if(csize != NULL) free(csize);
  if(cmsa  != NULL) {
    for(m = 0; m < nc; m++)
      if(cmsa[m] != NULL) esl_msa_Destroy(cmsa[m]);
    free(cmsa);
  }
  return status;
}

/* Function: select_node()
 * EPN, Fri Mar 23 08:48:37 2007 
 * Adapted from SRE's select_node() in maketestset.c originally written
 * for the PROFMARK HMMER benchmark.
 * 
 * 
 * Purpose:  Define clusters of the taxa (seqs) in the tree such
 *           that minimum disparity b/t any 2 seqs in different 
 *           clusters is greater than <mindiff> and the number of
 *           clusters is maximized. <ret_best> is the index of the node
 *           of the tree under which the largest cluster belongs.
 *           <ret_nc> is the number of clusters after clustering, 
 *           <ret_clust> is an array [0..T->N-1] specifying which
 *           cluster each taxa belongs to.
 *           
 *           For high disparities, this cluster may contain all
 *           the sequences, and we'll return the root node (0).
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * double    mindiff  - (see description above)
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 * int      *ret_xsize- size of largest cluster
 * int      *ret_best - RETURN: index of node of tree under which largest cluster belongs (see Purpose).
 * char     *errbuf   - buffer for error messages
 *
 * Returns: node index (as explained in Purpose)
 */
static int
select_node(ESL_TREE *T, double *diff, double mindiff, int **ret_clust, int *ret_nc, int *ret_xsize, int *ret_best, char *errbuf)
{
  int status;     /* Easel status code */
  ESL_STACK *ns1; /* stack for traversing tree */
  ESL_STACK *ns2; /* another stack for traversing tree */
  int c;	  /* counter for clusters */
  int best;       /* index of current best node */
  int maxsize;    /* size of cluster for best node */
  int n, np;      /* counters over tree nodes */
  int *clust;     /* [1..T->N-1] cluster number this seq is in */

  /*printf("in selec_node mindiff: %f T->N: %d\n", mindiff, T->N);*/
  /* set tree cladesizes if not already set */
  if(T->cladesize == NULL) 
    if((status = esl_tree_SetCladesizes(T)) != eslOK) ESL_FAIL(status, errbuf, "select_node(), esl_tree_SetCladeSizes error, status: %d.", status);

  ESL_ALLOC(clust, (sizeof(int) * T->N));
  esl_vec_ISet(clust, T->N, 0);

  if((ns1 = esl_stack_ICreate()) == NULL) ESL_FAIL(status, errbuf, "select_node(), failed to create a stack, probably out of memory, status: %d.", status);
  if((ns2 = esl_stack_ICreate()) == NULL) ESL_FAIL(status, errbuf, "select_node(), failed to create a stack, probably out of memory, status: %d.", status);

  /* push root on stack to start */
  if((status = esl_stack_IPush(ns1, 0)) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status);	
  maxsize  = 0;
  best     = 0;
  c        = 0;
  while (esl_stack_IPop(ns1, &n) != eslEOD) {
    if ((n == 0 || diff[T->parent[n]] > mindiff) &&
	diff[n] <= mindiff) { /* we're at a cluster */
      if (T->cladesize[n] > maxsize) {
	maxsize = T->cladesize[n];
	best = n;
      }
      /* determine all taxa in the clade rooted at n*/
      esl_stack_IPush(ns2, n);	
      while (esl_stack_IPop(ns2, &np) != eslEOD) {
	/*printf("np: %d T->left[np]: %d\n", np, T->left[np]);*/
	if(T->left[np]  <= 0) clust[(-1*T->left[np])]  = c;
	else { if((status = esl_stack_IPush(ns2, T->left[np])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
	if(T->right[np] <= 0) clust[(-1*T->right[np])]  = c;
	else { if((status = esl_stack_IPush(ns2, T->right[np])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
      }
      c++;
    }
    else {		/* we're not a cluster, keep traversing */
      /*printf("n: %d T->left[n]: %d\n", n, T->left[n]);*/
      if(T->left[n]  <= 0) clust[(-1*T->left[n])]  = c++; /* single seq with its own cluster */
      else { if((status = esl_stack_IPush(ns1, T->left[n])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
      if(T->right[n] <= 0) clust[(-1*T->right[n])] = c++; /* single seq with its own cluster */
      else { if((status = esl_stack_IPush(ns1, T->right[n])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
    }
  }
  esl_stack_Destroy(ns1);
  esl_stack_Destroy(ns2);
  *ret_nc = c;
  *ret_clust = clust;
  *ret_xsize = maxsize;
  *ret_best  = best;
  /*printf("nc: %d(%d) best: %d maxsize: %d nc: %d mindiff: %.3f\n\n", *ret_nc, c, best, maxsize, c, mindiff);
    for(n = 0; n < T->N; n++) printf("clust[%d]: %d\n", n, clust[n]);*/
  return eslOK;

 ERROR: 
  if(clust != NULL) free(clust);
  ESL_FAIL(status, errbuf, "select_node(), memory allocation error, status: %d.", status); 
}


/* Function: find_mindiff()
 * EPN, Fri Mar 23 18:59:42 2007
 * 
 * Purpose:  Given a tree resulting from single linkage clustering,
 *           find the min fractional difference (mindiff) that when used to
 *           define clusters (such that no seq in cluster A is less
 *           than mindiff different than any seq in cluster B), 
 *           gives either (A) (if(do_nc)) number of clusters >= target or
 *           (B) (if(do_nsize)) >= 1 cluster with >= target sequences
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * int      do_nsize  - TRUE  to find mindiff that gives >= 1 cluster with <target> seqs
 *                      FALSE to find mindiff that gives <target> clusters
 * int      target    - number of clusters (! if(do_nsize)) we want, or min size of
 *                      biggest cluster (if(do_nsize))
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 * int      *ret_xsize- size of largest cluster
 * int      *ret_best - cluster idx of largest cluster
 * int      *ret_mindiff - mindiff that achieves target
 * char     *errbuf   - buffer for error messages
 *
 * Returns: fractional difference (as explained in Purpose)
 */
static float
find_mindiff(ESL_TREE *T, double *diff, int do_nsize, int target, int **ret_clust, int *ret_nc, int *ret_xsize, int *ret_best, float *ret_mindiff, char *errbuf)
{
  int   status;
  float high_diff  = 1.0;
  float low_diff   = 0.0;
  int   high       = 0;
  int   low        = 0;
  float mindiff    = 0.5;
  int   curr_nc    = -1;
  int   curr_xsize = -1;
  int   curr       = -1;
  int   curr_best  = -1;
  int   keep_going = TRUE;
  float thresh     = 0.001;
  int  *clust      = NULL;

  /* Contract check */
  if(target > T->N) ESL_FAIL(eslEINCOMPAT, errbuf, "find_mindiff(), desired target is greater than number of seqs in the tree");

  while(keep_going) {
    if(clust != NULL) free(clust);
    if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_xsize, &curr_best, errbuf)) != eslOK) return status;
    curr = do_nsize ? curr_xsize : curr_nc;
    if(((!do_nsize) && (curr < target)) || ((do_nsize) && (curr >= target))) {
      high_diff  = mindiff;
      high       = curr;
      /*printf("LOWER   curr: %d mindiff: %f low: %f (%d) high: %f (%d)\n", curr, mindiff, low_diff, low, high_diff, high);*/
      mindiff   -= (mindiff - low_diff) / 2.;
      if((fabs(high_diff-0.) < thresh) && (fabs(low-0.) < thresh))  keep_going = FALSE; 
      /* stop, high and low have converged at 0. */
    }
    else {/* if(do_nsize) curr_nc > target_nc, else if(!do_nsize) curr_nc >= target_nc */
      low_diff   = mindiff;
      low        = curr;
      /*printf("GREATER curr: %d mindiff: %f low: %f (%d) high: %f (%d)\n", curr, mindiff, low_diff, low, high_diff, high);*/
      mindiff   += (high_diff - mindiff) / 2.;
      if(fabs(high_diff-low_diff) < thresh)  keep_going = FALSE; /* stop, high and low have converged */
    }
  }
  if(do_nsize) { 
    if(curr < target) { /* we couldn't reach our target in search due to precision */
      if(high >= target) { /* we could reach it at high */
	mindiff = high_diff;
	if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_xsize, &curr_best, errbuf)) != eslOK) return status;
      }
      else { /* we couldn't reach our threshold, this shouldn't happen */
	ESL_FAIL(eslEINVAL, errbuf,"Error in find_mindiff(), even with mindiff of %.5f can't produce cluster of size: %d\n", mindiff, target);
      }
    }
  }
  else { /* ! do_nsize, trying to achieve <target> clusters */
    /* it's possible we can't reach our target, if so, set mindiff as minimum value that gives 
     * less than target clusters. */
    if(curr != target) {
      /*printf("targ: %d curr: %d low: %d (%f) high: %d (%f)\n", target, curr, low, low_diff, high, high_diff);*/
      if(high < target) {
	mindiff = high;
	if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_xsize, &curr_best, errbuf)) != eslOK) return status;
      }
      else
	while(high > target) {
	  high += thresh;
	  if(high > 1.0)  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "find_mindiff(), mindiff has risen above 1.0");
	  mindiff = high;
	  if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_xsize, &curr_best, errbuf)) != eslOK) return status;
	  high = curr_nc;
	}
    }
  }
  /*printf("FINAL mindiff: %f\n", mindiff);  */
  *ret_nc      = curr_nc;
  *ret_clust   = clust;
  *ret_xsize   = curr_xsize;
  *ret_best    = curr_best;
  *ret_mindiff = mindiff;

  return eslOK;
}


/* determine_first_last_consensus_columns
 *                   
 * Given an MSA, determine the first and last consensus columns
 * occupied by each sequence
 */
static int determine_first_last_consensus_columns(ESL_MSA *msa, char *errbuf, int **ret_fA, int **ret_lA, int *ret_clen)
{
  int status;
  int clen = 0;
  int *fA = NULL;
  int *lA = NULL;
  int cpos = 0;
  int apos = 0;
  int i;

  /* contract check */
  if(msa->rf == NULL) { status = eslEINVAL; goto ERROR; }

  /* determine the first and last occupied consensus position in each sequence */
  ESL_ALLOC(fA, sizeof(int) * msa->nseq);
  ESL_ALLOC(lA, sizeof(int) * msa->nseq);

  /* count consensus columns */
  for(apos = 1; apos <= msa->alen; apos++) if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  esl_vec_ISet(lA, msa->nseq, 0);
  esl_vec_ISet(fA, msa->nseq, clen);
  /* this could be more efficient */
  for(i = 0; i < msa->nseq; i++) { 
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* apos is a consensus position */
	cpos++;
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][(apos+1)])) { /* cpos for seq i is not a gap */
	  fA[i] = ESL_MIN(fA[i], cpos);
	  lA[i] = ESL_MAX(lA[i], cpos);
	}
      }
    }
  }
  *ret_fA   = fA;
  *ret_lA   = lA;
  *ret_clen = clen;
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "determine_first_last_consensus_columns(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function:  dst_nongap_XPairId()
 * Synopsis:  Pairwise identity of two aligned digital seqs.
 *            Differs from esl_dst_XPairId() in that denominator is 
 *            seq identity is number of columns that are non-gap in both 
 *            sequences (instead of length of the shorter of the two seqs).
 *            
 * Incept:    EPN, Fri Jun 27 15:07:44 2008
 *
 * Args:      abc          - digital alphabet in use
 *            ax1          - aligned digital seq 1
 *            ax2          - aligned digital seq 2
 *            opt_pid      - optRETURN: pairwise identity, 0<=x<=1
 *            opt_nid      - optRETURN: # of identities
 *            opt_n        - optRETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <opt_distance>, <opt_nid>, <opt_n>
 *            contain the answers, for any of these that were passed
 *            non-<NULL> pointers.
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 */
int
dst_nongap_XPairId(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2, 
		   double *opt_distance, int *opt_nid, int *opt_n)
{
  int     status;
  int     idents;               /* total identical positions  */
  int     len;                  /* number of nongap colns in both seqs */
  int     i;                    /* position in aligned seqs   */

  idents = len = 0;
  for (i = 1; ax1[i] != eslDSQ_SENTINEL && ax2[i] != eslDSQ_SENTINEL; i++) 
    {
      if (esl_abc_XIsCanonical(abc, ax1[i]) && esl_abc_XIsCanonical(abc, ax2[i])) { 
	len++;
	if(ax1[i] == ax2[i]) idents++;
      }
    }

  if (ax1[i] != eslDSQ_SENTINEL || ax2[i] != eslDSQ_SENTINEL) 
    ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");

  if (opt_distance != NULL)  *opt_distance = ( len==0 ? 0. : (double) idents / (double) len );
  if (opt_nid      != NULL)  *opt_nid      = idents;
  if (opt_n        != NULL)  *opt_n        = len;
  return eslOK;

 ERROR:
  if (opt_distance != NULL)  *opt_distance = 0.;
  if (opt_nid      != NULL)  *opt_nid      = 0;
  if (opt_n        != NULL)  *opt_n        = 0;
  return status;
}


/* Function:  dst_nongap_XDiffMx()
 * Synopsis:  NxN difference matrix for N aligned digital seqs.         
 *            Differs from esl_dst_XDiffMx() in that denominator for
 *            seq identity is number of columns that are non-gap in both 
 *            sequences (instead of length of the shorter of the two seqs).
 *            
 * Incept:    EPN, Fri Jun 27 15:10:53 2008
 *
 * Args:      abc   - digital alphabet in use
 *            ax    - aligned dsq's, [0..N-1][1..alen]                  
 *            N     - number of aligned sequences
 *            ret_D - RETURN: NxN matrix of fractional differences
 *            
 * Returns:   <eslOK> on success, and <ret_D> contains the difference
 *            matrix; caller is obligated to free <D> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
dst_nongap_XDiffMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D)
{
  int status;
  ESL_DMATRIX *D = NULL;
  int i,j;

  if (( D = esl_dmatrix_Create(N,N) ) == NULL) goto ERROR;
  
  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      for (j = i+1; j < N; j++)
	{
	  status = dst_nongap_XPairId(abc, ax[i], ax[j], &(D->mx[i][j]), NULL, NULL);
	  if (status != eslOK) ESL_XEXCEPTION(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
	  D->mx[i][j] =  1.0 - D->mx[i][j];
	  D->mx[j][i] =  D->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D; else esl_dmatrix_Destroy(D);
  return eslOK;

 ERROR:
  if (D     != NULL)  esl_dmatrix_Destroy(D);
  if (ret_D != NULL) *ret_D = NULL;
  return status;
}



/* find_seqs_with_given_insert
 *                   
 * Given an MSA with RF annotation, determine which sequences have inserts after column <target>
 * of at least size <min> and at most size <max>. Fill an array <useme> of size msa->nseq with 
 * TRUE seq i has the insert, FALSE if it doesn't.
 */
static int find_seqs_with_given_insert(ESL_MSA *msa, char *errbuf, int target, int min, int max, int **ret_useme)
{
  int status;
  int apos, cpos, clen;
  int **ict;
  int *useme;
  int i;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "in find_seqs_with_given_insert(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --seq-ins.");

  ESL_ALLOC(useme,sizeof(int) * (msa->nseq));
  ESL_ALLOC(ict,  sizeof(int *) * (msa->alen+2));
  for(i = 0; i <= msa->alen; i++)
    {
      ESL_ALLOC(ict[i],  sizeof(int) * (msa->nseq));
      esl_vec_ISet(ict[i], (msa->nseq), 0);
    }

  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++)
    {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) 
	cpos++;
      else 
	for(i = 0; i < msa->nseq; i++)
	  if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) { 
	    ict[cpos][i]++;
	  }	  
    }
  clen = cpos;
  if(target > clen) ESL_XFAIL(eslEINVAL, errbuf, "--seq-ins <n> enabled with <n> = %d, but non-gap RF length of alignment is only %d columns.", target, clen);

  for(i = 0; i < msa->nseq; i++) { 
    useme[i] = ((ict[target][i] >= min) && (ict[target][i] <= max)) ?  TRUE : FALSE;
  }

  *ret_useme = useme;
  for(i = 0; i <= msa->alen; i++)
    free(ict[i]);
  free(ict);
  return eslOK;

 ERROR:
  return status;
}

/* minorize_msa()
 *                   
 * Given an MSA with #=GS <seq name> <<tag>> <minor set name>, make a new msa that
 * includes all seqs that have the same <minor set name> #=GS annotation for the
 * same <<tag>> (string value of <tag> is set at command line and passed into 
 * this function)
 * Also set the #=GC RF markup for each minor subset as either: 
 * (A) the #=GF <x> <y> with <x> equal to <minor set name> and <y>
 *     as the RF line.
 *
 * or, if no such #=GF markup exists define the consensus with the 
 * gap fraction rule, any column with <= <x> (from --gapthresh <x>)
 * gaps becomes an 'x' in the RF seq, all other columns are '.'.
 */
static int
minorize_msa(const ESL_GETOPTS *go, ESL_MSA *msa, char *errbuf, FILE *fp, char *tag)
{
  int    status;
  int   *useme = NULL;
  int    nalloc = 1;
  int    nmin = 0;
  int   *which_minor = NULL; /* [0..i..msa->nseq-1] which minor subset each sequence belongs to, -1 if none */
  char **minorA = NULL;      /* [0..m..nmin-1]      ptr to minor subset name in #=GS markup, only a ptr, don't free the char strings */
  int    f, g, i, m, mt;
  int    gt = -1;
  void  *tmp;
  char  *rf;
  int    apos;
  int    ip;
  int   *order;

  /* contract check */
  if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "-M requires #=GC RF markup in alignment.");
  if(msa->gs == NULL) ESL_FAIL(eslEINVAL, errbuf, "-M requires #=GS markup in alignment denoting minor subsets.");

  /* determine which tag matches <tag> */
  for(g = 0; g < msa->ngs; g++) { 
    if(strcmp(msa->gs_tag[g], tag) == 0) { 
      gt = g;
      break;
    }
  }
  if(gt == -1) ESL_FAIL(eslEINVAL, errbuf, "No #=GS markup has tag: %s\n", tag);

  /* determine which minor set each seq belongs to, reallocate minorA as we go and see new minor subsets */
  ESL_ALLOC(which_minor, sizeof(int) * msa->nseq);
  ESL_ALLOC(minorA, sizeof(char *) * nalloc);
  esl_vec_ISet(which_minor, msa->nseq, -1);

  for(i = 0; i < msa->nseq; i++) { 
    if(msa->gs[gt][i] != NULL) { 
      mt = -1;
      for(m = 0; m < nmin; m++) { 
	if(strcmp(minorA[m], msa->gs[gt][i]) == 0) { 
	  mt = m;
	  break;
	}
      }
      if(mt == -1) { 
	if((nmin+1) == nalloc) { 
	  nalloc++;
	  ESL_RALLOC(minorA, tmp, sizeof(char *) * nalloc);
	}
	minorA[nmin] = msa->gs[gt][i]; 
	mt = nmin++;
      }
      which_minor[i] = mt;
    }
  }
  for(i = 0; i < msa->nseq; i++) if(which_minor[i] == -1) ESL_FAIL(eslEINVAL, errbuf, "-M requires ALL sequences have #=GS markup with user supplied tag %s. Seq %d (%s) has none.\n", esl_opt_GetString(go, "-M"), i, msa->sqname[i]);

  /* Now make the minor alignments by keeping only the relevant subset of seqs.
   * We do not call esl_msa_MinimGaps() b/c we want the alignment length to be 
   * identical b/t all minor msas, and the RF also so cmalign knows how to 
   * map the minor alignments to the major alignment */
  ESL_MSA **minor_msaA;
  ESL_ALLOC(minor_msaA, sizeof(ESL_MSA *) * nmin);
  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for(m = 0; m < nmin; m++) { 
    for(i = 0; i < msa->nseq; i++) useme[i] = (which_minor[i] == m)  ? TRUE : FALSE;
    if((status = esl_msa_SequenceSubset(msa, useme, &(minor_msaA[m]))) != eslOK) ESL_FAIL(status, errbuf, "Error taking subset for minor subset %d with name: %s\n", m, minorA[m]);

    /* set name */
    esl_msa_SetName(minor_msaA[m], minorA[m]);

    /* unless --M-rf free RF annotation and set new annotation (--M-rf tells us to keep the initial RF annotation for all minor alignments */
    if(! esl_opt_GetBoolean(go, "--M-rf")) { 
      if(minor_msaA[m]->rf == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "Error creating minor alignment %d, RF is NULL.", m);
      free(minor_msaA[m]->rf);
      minor_msaA[m]->rf = NULL;
      /* Check if we have #=GF markup denoting the RF line for each minor subset */
      rf = NULL;
      for(f = 0; f < msa->ngf; f++)  
	if(strcmp(minorA[m], msa->gf_tag[f]) == 0) { rf = msa->gf[f]; break; }
      if(rf != NULL) { /* ensure the RF annotation is the length of the alignment */
	if(strlen(rf) != msa->alen) ESL_FAIL(eslEINCOMPAT, errbuf, "'#=GF %s <RF sequence>' markup is of length %d but it must be equal to aln length (%" PRId64 ").", msa->gf_tag[f], (int) strlen(rf), msa->alen);
	if((status = esl_strdup(rf, msa->alen, &(minor_msaA[m]->rf))) != eslOK) ESL_FAIL(status, errbuf, "Error duplicating RF for minor alignment %d\n", m);
	/* make sure minor_msaA[m]->rf does not have any non-gap columns where msa->rf has a gap (cmalign -M demands all minor consensus columns are also major consensus columns) */
	for (apos = 0; apos < minor_msaA[m]->alen; apos++) { 
	  if (! (esl_abc_CIsGap(minor_msaA[m]->abc, minor_msaA[m]->rf[apos]))) { 
	    if (esl_abc_CIsGap(msa->abc, msa->rf[apos])) ESL_FAIL(eslEINCOMPAT, errbuf, "'#=GF %s <RF sequence>' markup has a non-gap (%c char) at aln position %d, but the major alignment has a gap there! cmalign will choke on this.\n", msa->gf_tag[f], minor_msaA[m]->rf[(apos-1)], apos);
	  }
	}
      }
      else { /* no #=GF markup denoting RF line for alignment m existed in the input file, define it based on gaps */
	if((status = write_rf_gapthresh(go, errbuf, minor_msaA[m], esl_opt_GetReal(go, "--gapthresh"))) != eslOK) return status;
	/* careful, remember with cmalign -M, the minor alignments can only have training alignment column c defined as a consensus column 
	 * if it is also consensus in the major alignment, so we have to remove any minor_msaA[m]->rf 'x' columns that are gaps in <msa>
	 */
	for (apos = 0; apos < minor_msaA[m]->alen; apos++)
	  if (esl_abc_CIsGap(msa->abc, msa->rf[apos])) minor_msaA[m]->rf[apos] = '.';
      }
    }
  }

  /* Print out the alignments, first major, then minors */
  /* first, reorder major alignment so that it contains the minor alignment seqs in order */
  ip = 0;
  ESL_ALLOC(order, sizeof(int) * msa->nseq);
  for(m = 0; m < nmin; m++) { 
    for(i = 0; i < msa->nseq; i++) {
      if(which_minor[i] == m) order[i] = ip++;
    }
  }
  if((status = reorder_msa(msa, order, errbuf)) != eslOK) return status;

  esl_msa_Write(fp, msa, (esl_opt_GetBoolean(go, "-1")) ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM);
  for(m = 0; m < nmin; m++) { 
    esl_msa_Write(fp, minor_msaA[m], (esl_opt_GetBoolean(go, "-1")) ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM);
    esl_msa_Destroy(minor_msaA[m]);
  }
  free(minor_msaA);
  
  free(order);
  free(useme);
  free(which_minor);
  free(minorA);
  return eslOK;

 ERROR:
  if(which_minor != NULL) free(which_minor);
  if(minorA != NULL)      free(minorA);
  if(useme != NULL)       free(useme);
  return eslEMEM;
}


/* remove_gc_markup()
 *                   
 * Given a GC tag <tag>, remove that markup from an MSA.
 * Return eslEINVAL if <tag> does not exist.
 */
static int
remove_gc_markup(ESL_MSA *msa, char *errbuf, char *tag)
{
  int    does_not_exist = FALSE;

  /* Currently, we can only handle removal of parsed GC markup, RF, SS_cons, SA_cons, PP_cons 
   * (the main reason is b/c I didn't know how to deal with possibility of the ESL_KEYHASH in msa->gc_idx).
   */
  if (strcmp(tag, "RF") == 0) { 
    if   (msa->rf == NULL) does_not_exist = TRUE;
    else { free(msa->rf); msa->rf = NULL; }
  }
  else if(strcmp(tag, "SS_cons") == 0) { 
    if   (msa->ss_cons == NULL) does_not_exist = TRUE; 
    else { free(msa->ss_cons); msa->ss_cons = NULL; }
  }
  else if (strcmp(tag, "SA_cons") == 0) { 
    if   (msa->sa_cons == NULL) does_not_exist = TRUE;
    else { free(msa->sa_cons); msa->sa_cons = NULL; }
  }
  else if (strcmp(tag, "PP_cons") == 0) { 
    if   (msa->pp_cons == NULL) does_not_exist = TRUE;
    else { free(msa->pp_cons); msa->pp_cons = NULL; }
  }
  else { 
    ESL_FAIL(eslEINVAL, errbuf, "--rm-gc <s> only works if <s> is \'RF\', \'SS_cons\', \'SA_cons\', or \'PP_cons\'");
  }
  if(does_not_exist) { 
    ESL_FAIL(eslEINVAL, errbuf, "--rm-gc %s enabled but %s GC markup exists in the MSA.", tag, tag);
  }
  return eslOK;
}
