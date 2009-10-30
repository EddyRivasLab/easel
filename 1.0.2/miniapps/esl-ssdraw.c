/* Draw secondary structure diagrams given a postscript SS template. 
 * Initial development of this program was for SSU rRNA structures
 * with templates derived from Gutell's CRW. 
 *
 * EPN, Mon Jun 23 14:46:05 2008
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

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

#define ERRBUFSIZE 1024

#define ALIMODE 0
#define INDIMODE 1
#define SIMPLEMASKMODE 2
#define DRAWFILEMODE 3

#define RAINBOWRHSCHEME 0
#define RAINBOWRLSCHEME 1
#define NRAINBOWRHSCHEME 11
#define NRAINBOWRLSCHEME 11

#define RBSIXRHSCHEME 2
#define RBSIXRLSCHEME 3
#define NRBSIXRHSCHEME 6
#define NRBSIXRLSCHEME 6

#define NOC 9
#define CYANOC 0
#define MAGENTAOC 1
#define YELLOWOC 2
#define BLACKOC 3
#define LIGHTGREYOC 4
#define DARKGREYOC 5
#define REDOC 6
#define PURPLEOC 7
#define ORANGEOC 8

#define LEGTEXTNCHARS 60
#define NCMYK 4
#define ICYAN 0
#define IMAGENTA 1
#define IYELLOW 2
#define IBLACK 3

/* define the color for blank cells where no value is appropriate */
#define BLANKCYAN    0.0
#define BLANKMAGENTA 0.0
#define BLANKYELLOW  0.0
#define BLANKBLACK   0.5
/*if I extend to allow RED, BLUE, GREEN with combos of CMYK, but it's a pain right now due to implementation 
  #define IRED 4
  #define IBLUE 5
  #define IGREEN 6
  #define NCOLORS 7
*/
#define LEG_NBOXES  11
#define LEG_BOXSIZE 24.
#define LEG_MINFONTSIZE 10
#define LEGX_OFFSET 24.
#define LEGY_OFFSET -24.
#define SPECIAL_FONT "Courier-BoldOblique"
#define LEG_FONT "Courier-Bold"
#define LEG_EXTRA_COLUMNS 12 /* how many extra columns we need for printing stats in the legend */

#define DEFAULT_FONT "Courier-Bold"
#define RESIDUES_FONT "Helvetica-Bold"
#define HUNDREDS_FONT "Helvetica"

#define SS_BOXSIZE 8.

#define RESIDUES_FONTSIZE 8.
#define HUNDREDS_FONTSIZE 8.
#define LEG_FONTSIZE_UNSCALED 9.6
/*#define HEADER_FONTSIZE_UNSCALED 14.4*/
#define HEADER_FONTSIZE_UNSCALED 12
#define HEADER_MODELNAME_MAXCHARS 20
#define TICKS_LINEWIDTH 2.
#define BP_LINEWIDTH 1.

#define POSTSCRIPT_PAGEWIDTH 612.
#define POSTSCRIPT_PAGEHEIGHT 792.
#define PAGE_TOPBUF 12.
#define PAGE_SIDEBUF 12.
#define PAGE_BOTBUF 12.
#define COURIER_HEIGHT_WIDTH_RATIO 1.65

/* Structure: scheme_color_legend
 * Incept:    EPN, Thu Jun 25 20:20:38 2009
 *
 * Parameters describing a one-dimensional legend of colors
 * from a preset scheme for use in a SSPostscript_t data structure.
 */
typedef struct scheme_color_legend_s {
  int    scheme;            /* preset color scheme index */
  int    nbins;             /* number of colors (bins) in this scheme */
  char   *text1;            /* first line of text for legend, a single string */
  char   *text2;            /* second line of text for legend, a single string */
  float *limits;            /* [nbins+1] limits for each bin, limits[0] is min value we would expect to see, limits[nbins] is max */
  int *counts;              /* [nbins] number of cells we've painted each color */
  int *counts_masked;       /* [nbins] number of cells within mask ('1's) that we've painted each color */
} SchemeColorLegend_t;

/* Structure: onecell_color_legend
 * Incept:    EPN, Tue Sep 30 13:06:15 2008
 *
 * Parameters describing a single colored cell legend for a
 * SSPostscript_t data structure.
 */
typedef struct onecell_color_legend_s {
  float  col[NCMYK];        /* [CMYK] color value for the cell */
  char   *text;             /* text for legend */
  int    nres;              /* number of residues colored by the color in col[NCMYK] */
  int    nres_masked;       /* number of residues within a mask colored by the color in col[NCMYK] */
} OneCellColorLegend_t;

/* Structure: ss_postscript
 * Incept:    EPN, Mon Jun 23 15:50:30 2008
 *
 * A clumsy data structure for storing the information that will
 * become a postscript secondary structure diagram based on a 
 * template created by Robin Gutell and colleagues.
 *
 */
typedef struct ss_postscript_s {
  int     npage;        /* number of pages in eventual postscript */
  char   *modelname;    /* name of model, read from template file */
  int    *modeA;        /* [0..npage-1] page mode, ALIMODE, INDIMODE, or SIMPLEMASKMODE */
  char  **descA;        /* [0..npage-1] description for each page */
  int     desc_max_chars; /* max num characters for a page description */
  float   headerx;      /* x coordinate (bottom left corner) of header area */
  float   headery;      /* y coordinate (bottom left corner) of header area */
  float   headerx_charsize;/* size of a character in x-dimension in the header */
  float   headery_charsize;/* size of a character in y-dimension in the header */
  float   headerx_desc; /* x coordinate (bottom left corner) of header area */
  float   legx;         /* x coordinate (bottom left corner) of legend area */
  float   legy;         /* y coordinate (bottom left corner) of legend area */
  float   cur_legy;     /* y coordinate of current line in legend */
  float   legx_charsize;/* size of a character in x-dimension in the legend */
  float   legy_charsize;/* size of a character in y-dimension in the legend */
  int     legx_max_chars; /* max num residues in x direction we can print in legend before running off page */
  int     legy_max_chars; /* max num residues in y direction we can print in legend before running off page */
  int     legx_stats;   /* x position for printing stats in the legend */
  float   pagex_max;    /* max x position on page */
  float   pagey_max;    /* max y position on page */
  float   scale;        /* scale parameter, read from template file */
  char  **regurgA;      /* [0..nregurg-1][] lines from the template file to regurgitate, these are unchanged. */
  int     nregurg;      /* number of lines (char *'s) in the regurg_textAA 2D array */
  float  *hundredsxA;   /* [0..nhundreds-1] x value for hundreds (el 0 is for '100', 1 is for '200', etc.) */
  float  *hundredsyA;   /* [0..nhundreds-1] y value for hundreds (el 0 is for '100', 1 is for '200', etc.) */
  int     nhundreds;    /* number of elements in hundredsx and hundredsy */
  float  *ticksx1A;     /* [0..nticks-1] x begin value for ticks */
  float  *ticksx2A;     /* [0..nticks-1] x end   value for ticks */
  float  *ticksy1A;     /* [0..nticks-1] y begin value for ticks */
  float  *ticksy2A;     /* [0..nticks-1] x end   value for ticks */
  int     nticks;       /* number of ticks */
  float  *bpx1A;        /* [0..nbp-1] x begin value for bp connect line */
  float  *bpx2A;        /* [0..nbp-1] x end   value for bp connect line */
  float  *bpy1A;        /* [0..nbp-1] y begin value for bp connect line */
  float  *bpy2A;        /* [0..nbp-1] x end   value for bp connect line */
  int     nbp;          /* number of bp */
  float  *rxA;          /* [0..clen-1] x coordinate for each residue in the eventual postscript */
  float  *ryA;          /* [0..clen-1] y coordinate for each residue in the eventual postscript */
  int     clen;         /* the number of residues in the template file */
  char  **rrAA;         /* [0..npage-1][0..clen-1] residue character in the eventual postscript */
  float ***rcolAAA;     /* [0..npage-1][0..clen-1][0..3] color for block on page p, position c, CMYK in the eventual postscript */
  OneCellColorLegend_t ***occlAAA;/* [0..npage-1][0..l..nocclA[p]  ptr to one cell color legend l for page p */
  int     *nocclA;      /* [0..npage-1] number of one cell color legends for each page */
  SchemeColorLegend_t  **sclAA;/* [0..npage-1]  ptr to scheme color legend l for page p, NULL if none */
  char    *mask;        /* mask for this postscript, columns which are '0' get drawn differently */
  int      nalloc;      /* number of elements to add to arrays when reallocating */
  int     *msa_ct;      /* [1..ps->clen] CT array for msa this postscript corresponds to, 
			 * msa_ct[i] is the position that consensus residue i base pairs to, or 0 if i is unpaired. */
  int      msa_nbp;     /* number of bps read from current MSA (in msa_ct), should equal nbp, but only if bps read from template file */
  float    msa_avgid;   /* average id b/t all pairs of seqs in the MSA */
  float    msa_avglen;  /* average length of dealigned seqs in the MSA */
  int     *uaseqlenA;   /* [0..ps->msa->nseq-1] unaligned sequence length for all sequences in the MSA, only computed if --indi */
  int     *seqidxA;     /* [0..ps->npage-1] the sequence index in the MSA each page corresponds to, only valid if --indi */
  ESL_MSA *msa;         /* pointer to MSA this object corresponds to */
} SSPostscript_t;

static SSPostscript_t *create_sspostscript();
static int  setup_sspostscript(SSPostscript_t *ps, char *errbuf);
static OneCellColorLegend_t *create_onecell_colorlegend(float *cmykA, int nres, int nres_masked);
static SchemeColorLegend_t *create_scheme_colorlegend(int scheme, int ncols, float *limits);
static int  add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text, int legx_max_chars, char *errbuf);
static int  add_text_to_onecell_colorlegend(SSPostscript_t *ps, OneCellColorLegend_t *occl, char *text, int legx_max_chars, char *errbuf);
static int  add_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *text, char *errbuf);
static int  add_diffmask_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *mask1, char *mask2, char *errbuf);
static int  draw_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps);
static int  draw_legend_column_headers(FILE *fp, SSPostscript_t *ps, char *errbuf);
static int  draw_onecell_colorlegend(FILE *fp, OneCellColorLegend_t *occl, SSPostscript_t *ps, int occl_idx);
static int  draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page);
static void free_sspostscript(SSPostscript_t *ps);
static int  add_pages_sspostscript(SSPostscript_t *ps, int ntoadd, int page_mode);
static int  map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int **ret_a2c_map, int *ret_clen);
static int  parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, int msa_clen, SSPostscript_t **ret_ps);
static int  parse_template_page(ESL_FILEPARSER *efp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps);
static int  parse_modelname_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_scale_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf, int *ret_read_showpage);
static int  parse_regurgitate_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_text_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_lines_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  validate_justread_sspostscript(SSPostscript_t *ps, char *errbuf);
static int  validate_and_update_sspostscript_given_msa(const ESL_GETOPTS *go, SSPostscript_t *ps, ESL_MSA *msa, char *errbuf);
static int  individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  rf_seq_sspostscript (const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int ss_idx, int zerores_idx);
static int  delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float **hc_onecell, int incmask_idx, int excmask_idx);
static int  diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2, float **hc_onecell, int incboth_idx, int inc1_idx, int inc2_idx, int excboth_idx);
static int  compare_ints(const void *el1, const void *el2);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen, int *ret_mask_has_internal_zeroes);
static int  drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps);
static void PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, float wt);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  get_date(char *errbuf, char **ret_date);
static int  set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl, int within_mask);
static int  set_onecell_values(char *errbuf, float *vec, int ncolvals, float *onecolor);
static int  add_mask_to_ss_postscript(SSPostscript_t *ps, char *mask);
static int  draw_masked_block(FILE *fp, float x, float y, float *colvec, int do_circle_mask, int do_square_mask, int do_x_mask, int do_border, float boxsize);
static int  draw_header_and_footer(FILE *fp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int page, int pageidx2print);

static char banner[] = "draw postscript secondary structure diagrams.";
static char usage[]  = "[options] <msafile> <SS postscript template> <output postscript file name>\n\
The <msafile> must be in Stockholm format.";

#define MASKTYPEOPTS "-d,-c,-x" /* exclusive choice for mask types */
#define INCOMPATWITHSINGLEOPTS "--prob,--ins,--dall,--dint,--struct,--indi,--all,--dfile" /* exclusive choice for mask types */
#define INCOMPATWITHDFILEOPTS "-q,--prob,--ins,--dall,--dint,--struct,--indi,--all,--mask-col,--mask-diff" /* exclusive choice for mask types */

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },
  { "-q",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "do not draw info content diagram (or RF sequence if --indi)", 1 },
  { "--mask",   eslARG_INFILE, NULL, NULL, NULL, NULL,NULL,NULL,             "for all diagrams, mark masked ('0') columns from mask in <f>", 1 },
  { "--prob",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "draw posterior probability diagram(s)", 1 },

  { "--ins",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw insert diagram", 2 },
  { "--dall",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw delete diagram w/all deletions (incl. terminal deletes)", 2 },
  { "--dint",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw delete diagram w/only internal (non-terminal) deletions", 2 },
  { "--struct", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw structural information content diagram", 2 },

  { "--indi",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "draw diagrams for individual sequences instead of the aln", 3 },
  { "--all",    eslARG_NONE,  FALSE, NULL, NULL, NULL,"--indi", NULL,        "with --indi, draw individual diagrams of all sequences", 3 },

  { "--mask-u", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as squares", 4 },
  { "--mask-x", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as x's", 4 },
  { "--mask-a", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask-u or --mask-x, draw alternative mask style", 4 },

  { "--mask-col",eslARG_NONE, NULL, NULL, NULL, NULL,"--mask",  INCOMPATWITHSINGLEOPTS, "w/--mask draw black/cyan diagram denoting masked columns", 5 },
  { "--mask-diff",eslARG_INFILE,NULL, NULL, NULL, NULL,"--mask",INCOMPATWITHSINGLEOPTS, "with --mask-col <f1>, compare mask in <f1> to mask in <f>", 5 },

  { "--dfile",  eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, INCOMPATWITHDFILEOPTS, "read 'draw' file specifying >=1 diagrams", 6 },

  { "--no-leg", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw legend", 7 },
  { "--no-head",eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw header", 7 },
  { "--no-foot",eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw footer", 7 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  ESL_ALPHABET *abc     = NULL;	/* biological alphabet             */
  char         *alifile = NULL;	/* alignment file name             */
  char         *outfile = NULL;	/* output ps file name             */
  char         *templatefile = NULL; /* template file, specifying >= 1 SS diagrams 
				      * (each must have a unique consensus length) */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *msa     = NULL;	/* one multiple sequence alignment */
  int           status;		/* easel return code               */
  int           clen;           /* non-gap RF (consensus) length of each alignment */
  int           read_msa = FALSE; /* set to true when we read an alignment */
  int           apos;
  char          errbuf[ERRBUFSIZE];
  FILE           *ofp       = NULL;    /* output file for postscript */
  SSPostscript_t *ps        = NULL;    /* the postscript data structure we create */
  int            *hc_nbins  = NULL;
  float        ***hc_scheme = NULL;
  float         **hc_onecell = NULL;
  int           z;
  char          *command = NULL;
  char          *date = NULL;
  int master_mode;
  char *mask = NULL;
  int masklen;
  char *mask2 = NULL;
  int mask2len;
  int mask_has_internal_zeroes = FALSE;
  int mask2_has_internal_zeroes = FALSE;
  
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

  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for alignment summary diagrams (incompatible with --indi):");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for individual mode (require --indi):");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions controlling style of masked positions:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for drawing simple block-only diagrams of masks:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noption for reading in a file dictating colors:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      puts("\noptions for omitting parts of the diagram:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 3) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  alifile      = esl_opt_GetArg(go, 1);
  templatefile = esl_opt_GetArg(go, 2);
  outfile      = esl_opt_GetArg(go, 3);

  if((status = get_command(go, errbuf, &command)) != eslOK) esl_fatal(errbuf);
  if((status = get_date(errbuf, &date))           != eslOK) esl_fatal(errbuf);

  /****************/
  /* Premlinaries */
  /****************/
  /* allocate and fill predefined one-cell colors, these are hardcoded */
  ESL_ALLOC(hc_onecell, sizeof(float *) * NOC);
  for(z = 0; z < NOC; z++) hc_onecell[z] = NULL;
  for(z = 0; z < NOC; z++) { ESL_ALLOC(hc_onecell[z], sizeof(float) * NCMYK); }
  
  hc_onecell[CYANOC][0] = 1.0;
  hc_onecell[CYANOC][1] = 0.0;
  hc_onecell[CYANOC][2] = 0.0;
  hc_onecell[CYANOC][3] = 0.0;

  hc_onecell[MAGENTAOC][0] = 0.0;
  hc_onecell[MAGENTAOC][1] = 1.0;
  hc_onecell[MAGENTAOC][2] = 0.0;
  hc_onecell[MAGENTAOC][3] = 0.0;

  hc_onecell[YELLOWOC][0] = 0.0;
  hc_onecell[YELLOWOC][1] = 0.0;
  hc_onecell[YELLOWOC][2] = 1.0;
  hc_onecell[YELLOWOC][3] = 0.0;

  hc_onecell[BLACKOC][0] = 0.0;
  hc_onecell[BLACKOC][1] = 0.0;
  hc_onecell[BLACKOC][2] = 0.0;
  hc_onecell[BLACKOC][3] = 1.0;

  hc_onecell[LIGHTGREYOC][0] = 0.0;
  hc_onecell[LIGHTGREYOC][1] = 0.0;
  hc_onecell[LIGHTGREYOC][2] = 0.0;
  hc_onecell[LIGHTGREYOC][3] = 0.2;

  hc_onecell[DARKGREYOC][0] = 0.0;
  hc_onecell[DARKGREYOC][1] = 0.0;
  hc_onecell[DARKGREYOC][2] = 0.0;
  hc_onecell[DARKGREYOC][3] = 0.5;

  hc_onecell[REDOC][0] = 0.0;
  hc_onecell[REDOC][1] = 1.0;
  hc_onecell[REDOC][2] = 1.0;
  hc_onecell[REDOC][3] = 0.0;

  hc_onecell[PURPLEOC][0] = 1.0;
  hc_onecell[PURPLEOC][1] = 1.0;
  hc_onecell[PURPLEOC][2] = 0.0;
  hc_onecell[PURPLEOC][3] = 0.0;

  hc_onecell[ORANGEOC][0] = 0.0;
  hc_onecell[ORANGEOC][1] = 0.5;
  hc_onecell[ORANGEOC][2] = 1.0;
  hc_onecell[ORANGEOC][3] = 0.0;

  /***********************************/
  /* allocate and fill predefined color schemes, these are hardcoded */
  ESL_ALLOC(hc_scheme, sizeof(float **) * 4);
  for (z = 0; z < 4; z++) hc_scheme[z] = NULL;
  ESL_ALLOC(hc_scheme[0], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) hc_scheme[0][z] = NULL;
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[0][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[1], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) hc_scheme[1][z] = NULL;
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[1][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[2], sizeof(float *) * 6); 
  for(z = 0; z < 6; z++) hc_scheme[2][z] = NULL;
  for(z = 0; z < 6; z++) { ESL_ALLOC(hc_scheme[2][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[3], sizeof(float *) * 6); 
  for(z = 0; z < 6; z++) hc_scheme[3][z] = NULL;
  for(z = 0; z < 6; z++) { ESL_ALLOC(hc_scheme[3][z], sizeof(float) * NCMYK); }

  ESL_ALLOC(hc_nbins, sizeof(int) * 4);
  hc_nbins[0] = NRAINBOWRHSCHEME;
  hc_nbins[1] = NRAINBOWRLSCHEME;
  hc_nbins[2] = NRBSIXRHSCHEME;
  hc_nbins[3] = NRBSIXRLSCHEME;

  /***********************************/
  /*Scheme 0 and 1: Rainbow(red high) 11 is 0, Rainbow (red low) 11 is 1 */
  hc_scheme[0][ 0][0] = 0.92; hc_scheme[0][ 0][1] = 0.84; hc_scheme[0][ 0][2] = 0.00; hc_scheme[0][ 0][3] = 0.08; /*blue*/
  hc_scheme[1][10][0] = 0.92; hc_scheme[1][10][1] = 0.84; hc_scheme[1][10][2] = 0.00; hc_scheme[1][10][3] = 0.08; /*blue*/

  hc_scheme[0][ 1][0] = 0.78; hc_scheme[0][ 1][1] = 0.56; hc_scheme[0][ 1][2] = 0.00; hc_scheme[0][ 1][3] = 0.22;
  hc_scheme[1][ 9][0] = 0.78; hc_scheme[1][ 9][1] = 0.56; hc_scheme[1][ 9][2] = 0.00; hc_scheme[1][ 9][3] = 0.22;

  hc_scheme[0][ 2][0] = 0.50; hc_scheme[0][ 2][1] = 0.00; hc_scheme[0][ 2][2] = 0.00; hc_scheme[0][ 2][3] = 0.50;
  hc_scheme[1][ 8][0] = 0.50; hc_scheme[1][ 8][1] = 0.00; hc_scheme[1][ 8][2] = 0.00; hc_scheme[1][ 8][3] = 0.50;

  hc_scheme[0][ 3][0] = 0.61; hc_scheme[0][ 3][1] = 0.00; hc_scheme[0][ 3][2] = 0.56; hc_scheme[0][ 3][3] = 0.22;
  hc_scheme[1][ 7][0] = 0.61; hc_scheme[1][ 7][1] = 0.00; hc_scheme[1][ 7][2] = 0.56; hc_scheme[1][ 7][3] = 0.22;

  hc_scheme[0][ 4][0] = 0.42; hc_scheme[0][ 4][1] = 0.00; hc_scheme[0][ 4][2] = 1.00; hc_scheme[0][ 4][3] = 0.00;
  hc_scheme[1][ 6][0] = 0.42; hc_scheme[1][ 6][1] = 0.00; hc_scheme[1][ 6][2] = 1.00; hc_scheme[1][ 6][3] = 0.00;

  hc_scheme[0][ 5][0] = 0.00; hc_scheme[0][ 5][1] = 0.00; hc_scheme[0][ 5][2] = 1.00; hc_scheme[0][ 5][3] = 0.00;
  hc_scheme[1][ 5][0] = 0.00; hc_scheme[1][ 5][1] = 0.00; hc_scheme[1][ 5][2] = 1.00; hc_scheme[1][ 5][3] = 0.00;

  hc_scheme[0][ 6][0] = 0.00; hc_scheme[0][ 6][1] = 0.21; hc_scheme[0][ 6][2] = 1.00; hc_scheme[0][ 6][3] = 0.00;
  hc_scheme[1][ 4][0] = 0.00; hc_scheme[1][ 4][1] = 0.21; hc_scheme[1][ 4][2] = 1.00; hc_scheme[1][ 4][3] = 0.00;

  hc_scheme[0][ 7][0] = 0.00; hc_scheme[0][ 7][1] = 0.42; hc_scheme[0][ 7][2] = 1.00; hc_scheme[0][ 7][3] = 0.00;
  hc_scheme[1][ 3][0] = 0.00; hc_scheme[1][ 3][1] = 0.42; hc_scheme[1][ 3][2] = 1.00; hc_scheme[1][ 3][3] = 0.00;

  hc_scheme[0][ 8][0] = 0.00; hc_scheme[0][ 8][1] = 0.63; hc_scheme[0][ 8][2] = 1.00; hc_scheme[0][ 8][3] = 0.00;
  hc_scheme[1][ 2][0] = 0.00; hc_scheme[1][ 2][1] = 0.63; hc_scheme[1][ 2][2] = 1.00; hc_scheme[1][ 2][3] = 0.00;

  hc_scheme[0][ 9][0] = 0.00; hc_scheme[0][ 9][1] = 0.84; hc_scheme[0][ 9][2] = 1.00; hc_scheme[0][ 9][3] = 0.00;
  hc_scheme[1][ 1][0] = 0.00; hc_scheme[1][ 1][1] = 0.84; hc_scheme[1][ 1][2] = 1.00; hc_scheme[1][ 1][3] = 0.00;

  hc_scheme[0][10][0] = 0.00; hc_scheme[0][10][1] = 0.94; hc_scheme[0][10][2] = 1.00; hc_scheme[0][10][3] = 0.00; /*red*/
  hc_scheme[1][ 0][0] = 0.00; hc_scheme[1][ 0][1] = 0.94; hc_scheme[1][ 0][2] = 1.00; hc_scheme[1][ 0][3] = 0.00; /*red*/
  /***********************************/
  /*Scheme 0 and 1: Rainbow(red high) 11 is 0, Rainbow (red low) 11 is 1 */
  hc_scheme[2][0][0] = 0.92; hc_scheme[2][0][1] = 0.84; hc_scheme[2][0][2] = 0.00; hc_scheme[2][0][3] = 0.08; /*blue*/
  hc_scheme[3][5][0] = 0.92; hc_scheme[3][5][1] = 0.84; hc_scheme[3][5][2] = 0.00; hc_scheme[3][5][3] = 0.08; /*blue*/

  hc_scheme[2][1][0] = 0.50; hc_scheme[2][1][1] = 0.00; hc_scheme[2][1][2] = 0.00; hc_scheme[2][1][3] = 0.50;
  hc_scheme[3][4][0] = 0.50; hc_scheme[3][4][1] = 0.00; hc_scheme[3][4][2] = 0.00; hc_scheme[3][4][3] = 0.50;

  hc_scheme[2][2][0] = 0.42; hc_scheme[2][2][1] = 0.00; hc_scheme[2][2][2] = 1.00; hc_scheme[2][2][3] = 0.00;
  hc_scheme[3][3][0] = 0.42; hc_scheme[3][3][1] = 0.00; hc_scheme[3][3][2] = 1.00; hc_scheme[3][3][3] = 0.00;

  hc_scheme[2][3][0] = 0.00; hc_scheme[2][3][1] = 0.21; hc_scheme[2][3][2] = 1.00; hc_scheme[2][3][3] = 0.00;
  hc_scheme[3][2][0] = 0.00; hc_scheme[3][2][1] = 0.21; hc_scheme[3][2][2] = 1.00; hc_scheme[3][2][3] = 0.00;

  hc_scheme[2][4][0] = 0.00; hc_scheme[2][4][1] = 0.63; hc_scheme[2][4][2] = 1.00; hc_scheme[2][4][3] = 0.00;
  hc_scheme[3][1][0] = 0.00; hc_scheme[3][1][1] = 0.63; hc_scheme[3][1][2] = 1.00; hc_scheme[3][1][3] = 0.00;

  hc_scheme[2][5][0] = 0.00; hc_scheme[2][5][1] = 0.94; hc_scheme[2][5][2] = 1.00; hc_scheme[2][5][3] = 0.00; /*red*/
  hc_scheme[3][0][0] = 0.00; hc_scheme[3][0][1] = 0.94; hc_scheme[3][0][2] = 1.00; hc_scheme[3][0][3] = 0.00; /*red*/
  /***************************************************************/
  master_mode = ALIMODE;
  if(esl_opt_GetBoolean(go, "--indi"))       master_mode = INDIMODE;
  if(esl_opt_GetBoolean(go, "--mask-col"))   master_mode = SIMPLEMASKMODE;
  if(! esl_opt_IsDefault(go, "--mask-diff")) master_mode = SIMPLEMASKMODE;
  if(! esl_opt_IsDefault(go, "--dfile"))     master_mode = DRAWFILEMODE;

  /***********************************************
   * Open the MSA file; determine alphabet; set for digital input
   ***********************************************/
  
  fmt = eslMSAFILE_STOCKHOLM;
  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed with error %d\n", status);

  /* open postscript output file for writing */
  if ((ofp = fopen(outfile, "w")) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open output postscript file %s\n", esl_opt_GetArg(go, 2));

  /* Assert RNA, it's the ribosome */
  abc = esl_alphabet_Create(eslRNA);
  afp->abc = abc;

  /* Read the mask files, if nec */
  if(! esl_opt_IsDefault(go, "--mask")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask"), errbuf, &mask, &masklen, &mask_has_internal_zeroes)) != eslOK) esl_fatal(errbuf);
  }
  if(! esl_opt_IsDefault(go, "--mask-diff")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-diff"), errbuf, &mask2, &mask2len, &mask2_has_internal_zeroes)) != eslOK) esl_fatal(errbuf);
    if(masklen != mask2len) esl_fatal("Mask in %f length (%d) differs from mask in %f (%d)!", esl_opt_GetString(go, "--mask"), masklen, esl_opt_GetString(go, "--mask-diff"), mask2len);
  }

  /* Read the alignment */
  if((status = esl_msa_Read(afp, &msa)) == eslOK) { 
    read_msa = TRUE;
    msa->abc = abc;
    if(msa->rf == NULL) esl_fatal("First MSA in %s does not have RF annotation.", alifile);
    clen = 0;
    for(apos = 0; apos < msa->alen; apos++) if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) clen++;
    
    /* We've read the alignment, now read the template postscript file (we do this second b/c the RF len of the alignment tells us which postscript template to use) */
    if((status = parse_template_file(templatefile, go, errbuf, clen, &ps) != eslOK)) esl_fatal(errbuf);
    /* determine position for header and legend */
    if((status = setup_sspostscript(ps, errbuf) != eslOK)) esl_fatal(errbuf);

    if(ps->clen == 0)    esl_fatal("MSA has consensus (non-gap RF) length of %d which != template file consensus length of %d.", clen, ps->clen);
    if(clen != ps->clen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != template file consensus length of %d.", clen, ps->clen);

    /* add the mask if there is one */
    if(mask != NULL && (master_mode != SIMPLEMASKMODE)) add_mask_to_ss_postscript(ps, mask);
    if(mask != NULL && ps->clen != masklen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != lane mask length of %d from mask file %s.", clen, masklen, esl_opt_GetString(go, "--mask"));

    if((status = validate_and_update_sspostscript_given_msa(go, ps, msa, errbuf)) != eslOK) esl_fatal(errbuf);
    
    if(master_mode == ALIMODE) { 
      if(! esl_opt_GetBoolean(go, "-q")) { 
	if((status = infocontent_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--struct")) { 
	if((status = structural_infocontent_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, DARKGREYOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--ins")) { /* make a new postscript page marking insertions */
	if((status = insert_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dall")) { /* make a new postscript page marking all deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, TRUE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dint")) { /* make a new postscript page marking internal deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, FALSE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--prob")) { 
	if((status = posteriors_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
    }
    else if(master_mode == INDIMODE) { 
      if(! esl_opt_GetBoolean(go, "-q")) { 
	if((status = rf_seq_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
      }
      /* determine which sequences we are going to draw diagrams for */
      if(esl_opt_GetBoolean(go, "--all")) { 
	if((status = individual_seqs_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
	if(esl_opt_GetBoolean(go, "--prob")) { 
	  if((status = posteriors_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
	}
      }
    }
    else if(master_mode == SIMPLEMASKMODE) { 
      if(esl_opt_GetBoolean(go, "--mask-col")) { 
	if(ps->clen != masklen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != lane mask length of %d.", clen, masklen);
	/* Paint positions excluded by the mask magenta, unless the mask has zero internal exclusions.
	 * Such a mask is a 'truncating' mask, that excludes only a 5' contiguous set of columns, and a 3' contiguous set of columns, in this case paint excluded positions light grey. */
	if((status = colormask_sspostscript(go, errbuf, ps, msa, mask, hc_onecell, BLACKOC, (mask_has_internal_zeroes ? MAGENTAOC : LIGHTGREYOC))) != eslOK) esl_fatal(errbuf);
      }
      if(! esl_opt_IsDefault(go, "--mask-diff")) { 
	if((status = diffmask_sspostscript(go, errbuf, ps, msa, mask, mask2, hc_onecell, BLACKOC, CYANOC, MAGENTAOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
    }
    else if(master_mode == DRAWFILEMODE) { 
      if((status = drawfile2sspostscript(go, errbuf, ps)) != eslOK) esl_fatal(errbuf);
    }

    if((status = draw_sspostscript(ofp, go, errbuf, command, date, hc_scheme, ps)) != eslOK) esl_fatal(errbuf);
    free(command);
    fclose(ofp);
    esl_msa_Destroy(msa);
  }
  else { 
    /* If the msa read failed, we drop out to here with an informative status code. 
     */
    if (status == eslEFORMAT) 
      esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
		afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
    else if (status != eslEOF)
      esl_fatal("Alignment file read failed with error code %d\n", status);
  }
  if (!read_msa) { 
    esl_fatal("No alignments found in file %s\n", alifile);
  }
  /* Cleanup, normal return
   */
  if(mask != NULL) free(mask);
  if(date != NULL) free(date);
  free_sspostscript(ps);
  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  free(hc_nbins);
  for(z = 0; z < NOC; z++) free(hc_onecell[z]);
  free(hc_onecell);
  for(z = 0; z < 11; z++) free(hc_scheme[0][z]);
  for(z = 0; z < 11; z++) free(hc_scheme[1][z]);
  for(z = 0; z < 6; z++) free(hc_scheme[2][z]);
  for(z = 0; z < 6; z++) free(hc_scheme[3][z]);
  free(hc_scheme[0]);
  free(hc_scheme[1]);
  free(hc_scheme[2]);
  free(hc_scheme[3]);
  free(hc_scheme);
  return 0;

 ERROR: 
  esl_fatal("Memory allocation error in main().");
}

/* Function: create_sspostscript()
 * 
 * Purpose:  Create and initialize a SS postscript data structure.
 * Return:   ps
 */
SSPostscript_t *
create_sspostscript()
{
  int status;
  SSPostscript_t *ps;

  ESL_ALLOC(ps, sizeof(SSPostscript_t));

  ps->npage    = 0;
  ps->modelname = NULL;
  ps->modeA = NULL;
  ps->descA = NULL;
  ps->headerx = 0.;
  ps->headery = 0.;
  ps->headerx_desc = 0.;
  ps->headerx_charsize = 0.;
  ps->headery_charsize = 0.;
  ps->desc_max_chars = 0;
  ps->legx = 0.;
  ps->legy = 0.;
  ps->cur_legy = 0.;
  ps->legx_charsize = 0.;
  ps->legy_charsize = 0.;
  ps->legx_max_chars = 0;
  ps->legx_stats = 0.;
  ps->pagex_max = 0.;
  ps->pagey_max = 0.;
  ps->scale = 0.;
  ps->regurgA  = NULL;
  ps->nregurg  = 0;
  ps->hundredsxA = ps->hundredsyA = NULL;
  ps->nhundreds = 0;
  ps->ticksx1A = ps->ticksx2A = ps->ticksy1A = ps->ticksy2A = NULL;
  ps->nticks = 0;
  ps->bpx1A = ps->bpx2A = ps->bpy1A = ps->bpy2A = NULL;
  ps->nbp = 0;
  ps->rxA = ps->ryA = NULL;
  ps->clen = 0;
  ps->rrAA        = NULL;
  ps->rcolAAA     = NULL;
  ps->occlAAA     = NULL;
  ps->nocclA      = NULL;
  ps->sclAA       = NULL;
  ps->mask        = NULL;
  ps->nalloc      = 50;
  ps->msa_ct      = NULL;
  ps->msa_nbp     = 0;
  ps->msa_avglen  = 0.;
  ps->msa_avgid   = 0.;
  ps->uaseqlenA   = NULL;
  ps->seqidxA     = NULL;
  ps->msa         = NULL;
  return ps;

 ERROR: esl_fatal("create_sspostscript(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: setup_sspostscript()
 * 
 * Purpose:  Determine positions for header and legend in a SSPostscript_t()
 * Return:   eslOK
 */
int
setup_sspostscript(SSPostscript_t *ps, char *errbuf)
{
  if(ps->clen == 0) ESL_FAIL(eslEINVAL, errbuf, "Failed to ready any residues in template file.");

  /* set up legx, legy, this is a hack (takes advantage of position of 3' residue in all SSU models) */
  ps->legx = ps->rxA[ps->clen-1] + LEGX_OFFSET;
  ps->legy = ps->ryA[ps->clen-1] + LEGY_OFFSET;
  ps->cur_legy = ps->legy;

  ps->pagex_max = POSTSCRIPT_PAGEWIDTH / ps->scale;
  ps->pagey_max = POSTSCRIPT_PAGEHEIGHT / ps->scale;

  ps->headerx = 0. + PAGE_SIDEBUF;
  ps->headery = ps->pagey_max - PAGE_TOPBUF - ((HEADER_FONTSIZE_UNSCALED) / ps->scale);

  /* determine max number of residues we can print before we run off the page in the legend section */
  float xroom, yroom;
  xroom  = ps->pagex_max - ps->legx;
  yroom  = (ps->legy - ps->pagey_max) * -1;
  ps->legx_charsize = (LEG_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  ps->legy_charsize = (LEG_FONTSIZE_UNSCALED) / ps->scale; 
  ps->legx_max_chars = (int) (xroom / ps->legx_charsize);
  ps->legy_max_chars = (int) (yroom / ps->legy_charsize);
  ps->legx_stats     = ps->pagex_max - PAGE_SIDEBUF - (LEG_EXTRA_COLUMNS * ps->legx_charsize);

  /* determine max size of description that will fit in header */
  float header_fontwidth, header_max_chars;
  header_fontwidth      = (HEADER_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  ps->headerx_charsize  = (HEADER_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  header_max_chars      = (int) (ps->pagex_max / ps->headerx_charsize) - 2;
  ps->headery_charsize  = (HEADER_FONTSIZE_UNSCALED) / ps->scale; 
  ps->desc_max_chars    = header_max_chars - (HEADER_MODELNAME_MAXCHARS + 6 + 6 + 8 +2); /*6,6,8 for #res,#bps,#seq plus 2 spaces each, plus 2 for after name */
  ps->headerx_desc      = ps->pagex_max - PAGE_SIDEBUF - (ps->desc_max_chars * ps->headerx_charsize);

  return eslOK;
}

/* Function: free_sspostscript()
 * 
 * Purpose:  Free a SS postscript data structure.
 * Return:   (void)
 */
void
free_sspostscript(SSPostscript_t *ps)
{
  int i, p, c, l;

  if(ps->modelname != NULL) free(ps->modelname);

  if(ps->modeA != NULL)  free(ps->modeA);

  if(ps->descA != NULL) {
    for(i = 0; i < ps->npage; i++) { 
      if(ps->descA[i] != NULL) {
	free(ps->descA[i]);
      }
    }
    free(ps->descA);
  }

  if(ps->regurgA != NULL) {
    for(i = 0; i < ps->nregurg; i++) { 
      if(ps->regurgA[i] != NULL) {
	free(ps->regurgA[i]);
      }
    }
    free(ps->regurgA);
  }

  if(ps->hundredsxA != NULL) free(ps->hundredsxA);
  if(ps->hundredsyA != NULL) free(ps->hundredsyA);
  if(ps->ticksx1A != NULL) free(ps->ticksx1A);
  if(ps->ticksy1A != NULL) free(ps->ticksy1A);
  if(ps->ticksx2A != NULL) free(ps->ticksx2A);
  if(ps->ticksy2A != NULL) free(ps->ticksy2A);
  if(ps->bpx1A != NULL) free(ps->bpx1A);
  if(ps->bpy1A != NULL) free(ps->bpy1A);
  if(ps->bpx2A != NULL) free(ps->bpx2A);
  if(ps->bpy2A != NULL) free(ps->bpy2A);
  if(ps->rxA != NULL) free(ps->rxA);
  if(ps->ryA != NULL) free(ps->ryA);

  if(ps->rrAA != NULL) { 
    for(p = 0; p < ps->npage; p++) 
      if(ps->rrAA[p] != NULL) free(ps->rrAA[p]);
    free(ps->rrAA);
  }

  if(ps->rcolAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->rcolAAA[p] != NULL) { 
	for(c = 0; c < ps->clen; c++) free(ps->rcolAAA[p][c]);
	free(ps->rcolAAA[p]); 
      }
    }
    free(ps->rcolAAA); 
  }

  if(ps->occlAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->occlAAA[p] != NULL) { 
	for(l = 0; l < ps->nocclA[p]; l++) { 
	  if(ps->occlAAA[p][l]->text != NULL) free(ps->occlAAA[p][l]->text);
	  free(ps->occlAAA[p][l]); /* rest is statically allocated memory */
	}
      free(ps->occlAAA[p]);
      }
    }
    free(ps->occlAAA);
  }

  if(ps->sclAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->sclAA[p] != NULL) { 
	if(ps->sclAA[p]->limits != NULL) free(ps->sclAA[p]->limits);
	if(ps->sclAA[p]->counts != NULL) free(ps->sclAA[p]->counts);
	if(ps->sclAA[p]->counts_masked != NULL) free(ps->sclAA[p]->counts_masked);
	if(ps->sclAA[p]->text1 != NULL) free(ps->sclAA[p]->text1);
	if(ps->sclAA[p]->text2 != NULL) free(ps->sclAA[p]->text2);
	free(ps->sclAA[p]); /* statically allocated memory */
      }
    }
    free(ps->sclAA);
  }

  if(ps->nocclA != NULL) free(ps->nocclA);
  if(ps->msa_ct != NULL) free(ps->msa_ct);
  if(ps->mask != NULL) free(ps->mask);
  if(ps->seqidxA != NULL) free(ps->seqidxA);
  if(ps->uaseqlenA != NULL) free(ps->uaseqlenA);

  free(ps);
  return;
}

/* Function: create_onecell_colorlegend()
 * 
 * Purpose:  Create and initialize a one cell color legend data structure.
 * Return:   occl
 */
OneCellColorLegend_t *
create_onecell_colorlegend(float *col, int nres, int nres_masked)
{
  int status;
  OneCellColorLegend_t *occl;

  ESL_ALLOC(occl, sizeof(OneCellColorLegend_t));

  /* initialize */
  esl_vec_FSet(occl->col, NCMYK, 0.);
  occl->text = NULL;

  /* set caller specified values */
  esl_vec_FCopy(col, NCMYK, occl->col);

  occl->nres = nres;
  occl->nres_masked = nres_masked;
  return occl;

 ERROR: esl_fatal("create_onecell_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: create_scheme_colorlegend()
 * 
 * Purpose:  Create and initialize a scheme color legend data structure.
 * Return:   scl
 */
SchemeColorLegend_t *
create_scheme_colorlegend(int scheme, int nbins, float *limits)
{
  int status;
  SchemeColorLegend_t *scl;
  int i;
  int *counts;
  int *counts_masked;
  ESL_ALLOC(scl, sizeof(SchemeColorLegend_t));

  /* initialize */
  scl->text1 = NULL;
  scl->text2 = NULL;
  
  /* set caller specified values */
  scl->scheme = scheme;
  scl->nbins = nbins;
  ESL_ALLOC(scl->limits, sizeof(float) * (nbins+1));
  for(i = 0; i <= nbins; i++) { scl->limits[i] = limits[i]; }

  ESL_ALLOC(counts, sizeof(int) * nbins);
  ESL_ALLOC(counts_masked, sizeof(int) * nbins);
  esl_vec_ISet(counts, nbins, 0);
  esl_vec_ISet(counts_masked, nbins, 0);
  scl->counts = counts;
  scl->counts_masked = counts_masked;

  return scl;

 ERROR: esl_fatal("create_scheme_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: add_text_to_scheme_colorlegend()
 * 
 * Purpose:  Add text to an existing scheme color legend data structure.
 * Throws:   Exception if the text is too long.
 */
int
add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text, int legx_max_chars, char *errbuf)
{
  int status;
  int i, idx;
  int max_chars_per_line;

  if(scl->text1 != NULL) esl_fatal("add_text_to_scheme_colorlegend(), text already exists!\n"); 
  if(scl->text2 != NULL) esl_fatal("add_text_to_scheme_colorlegend(), text already exists!\n"); 
  if(text == NULL) esl_fatal("add_text_to_scheme_colorlegend(), passed in text is NULL!\n"); 

  max_chars_per_line = legx_max_chars - LEG_EXTRA_COLUMNS -2;
  if(((int) strlen(text)) <= max_chars_per_line) {
    /* case 1, entire text can fit in one line */
    if((status = esl_strdup(text, -1, &(scl->text1))) != eslOK) esl_fatal("add_text_to_scheme_colorlegend(), error copying text");
    return eslOK;
  }
  else if(((int) strlen(text)) > ((2 * max_chars_per_line) - 6)) { 
    /* case 2, entire text can't even fit in two lines, 
     * (this is inexact doesn't account for size of words which is an issue b/c we break at newline,
     * this is why I have the extra '- 6');
     */
    ESL_FAIL(eslEINVAL, errbuf, "add_text_to_scheme_colorlegend(), text is %d chars, max allowed is %d (%s)\n", (int) strlen(text), ((2 * max_chars_per_line) - 6), text);
  }
  else { /* split it up into two lines */
    idx = max_chars_per_line - 1;
    while(text[idx] != ' ') { 
      idx--;
      if(idx < 0) ESL_FAIL(eslEINVAL, errbuf, "add_text_to_scheme_colorlegend(), couldn't find a breakpoint for splitting the string (%s)\n", text);
    }
    /* copy first string into text1 */
    ESL_ALLOC(scl->text1, sizeof(char) * (idx+1));
    for(i = 0; i < idx; i++) scl->text1[i] = text[i];
    scl->text1[idx] = '\0';
    /* copy remainder into text2 */
    int len = (int) strlen(text);
    idx++;
    ESL_ALLOC(scl->text2, sizeof(char) * (len - idx+1));
    for(i = idx; i < len; i++) scl->text2[i-idx] = text[i];
    scl->text2[len-idx] = '\0';
  }
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Error adding text to scheme_colorlegend probably out of memory");
}


/* Function: add_text_to_onecell_colorlegend()
 * 
 * Purpose:  Add text to an existing one cell color legend data structure.
 * Throws:   Exception if the text is too long.
 */
int
add_text_to_onecell_colorlegend(SSPostscript_t *ps, OneCellColorLegend_t *occl, char *text, int legx_max_chars, char *errbuf)
{
  int status;
  int max_chars_per_line;
  if(occl->text != NULL) esl_fatal("add_text_to_onecell_colorlegend(), text already exists!\n"); 
  if(text == NULL) esl_fatal("add_text_to_onecell_colorlegend(), passed in text is NULL!\n"); 

  max_chars_per_line = legx_max_chars - LEG_EXTRA_COLUMNS - 2 - ((int) ((LEG_BOXSIZE * 1.5) / ps->legx_charsize));
  /*printf("max: %d cur: %d text: %s\n", max_chars_per_line, (int) strlen(text), text);*/
  if(((int) strlen(text)) > (max_chars_per_line)) { 
    ESL_FAIL(eslEINVAL, errbuf, "add_text_to_onecell_colorlegend(), text is %d chars, max allowed is %d (%s)\n", (int) strlen(text), max_chars_per_line, text);
  }
  if((status = esl_strdup(text, -1, &(occl->text))) != eslOK) esl_fatal("add_text_to_onecell_colorlegend(), error copying text");
  return eslOK;
}

/* Function: add_page_desc_to_sspostscript()
 * 
 * Purpose:  Add text describing a particular page of a postscript object.
 * Throws:   Exception if the text is too long.
 */
int
add_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *text, char *errbuf)
{
  int status;
  int i, j;
  int max_both_lines;

  if(ps->descA[page] != NULL) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), description for page %d already exists!\n", page); 
  if(text == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), passed in text is NULL!\n"); 

  max_both_lines =(2. * ps->desc_max_chars);
  if(ps->modeA[page] == INDIMODE || ps->modeA[page] == SIMPLEMASKMODE) 
    max_both_lines--; /* b/c we have to add a '-' to split up the larger strings onto 2 lines */

  /* check to see if we can fit the text onto two lines of max width ps->desc_max_chars */
  int textlen = (int) strlen(text);
  if(textlen <= ps->desc_max_chars) { /* fine, this will fit on one line */
    if((status = esl_strdup(text, -1, &(ps->descA[page]))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), error copying text");
  }
  else if (textlen <= max_both_lines) { 
    if(ps->modeA[page] == ALIMODE) { 
      /* maybe fine, make sure there's a break point (space ' ') that will break this string into two strings of length <= ps->desc_max_chars */
      i = ps->desc_max_chars;
      while(text[i] != ' ') { 
	i--; 
	if(i < 0) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), first word of text (%s) is more than max allowed of %d chars", text, ps->desc_max_chars);
      }
      /* found last space before max width, make sure it breaks the line up into two chunks of valid size */
      if((textlen - (i+1)) <= ps->desc_max_chars) { /* we're good */
	if((status = esl_strdup(text, -1, &(ps->descA[page]))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), error copying text");
	ps->descA[page][i] = '\n'; /* so we can remember where the break is */
      }
      else ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), couldn't find break point (' ') for partitioning text into two valid size chunks (%s)", text);
    }
    else { /* INDIMODE or SIMPLEMASKMODE, sequence/mask name bigger than 1 line, but not 2, we put a '-' in it at the end of line 1 and add a '\n' so we remember where it was */
      ESL_ALLOC(ps->descA[page], sizeof(char) * (textlen + 3)); /* +3 so we have space for the extra '-' and '\n' */
      for(i = 0; i < ps->desc_max_chars; i++) ps->descA[page][i] = text[i];
      i = ps->desc_max_chars;
      ps->descA[page][i] = '-';
      i++; 
      ps->descA[page][i] = '\n';
      for(i = ps->desc_max_chars; i < textlen; i++) ps->descA[page][i+2] = text[i];
      ps->descA[page][textlen+2] = '\0';
    }
  }
  else /* the text won't fit on two lines */
    if(ps->modeA[page] != INDIMODE) { /* not fine, this won't fit on two lines */
      ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), text is %d chars, max allowed is %d (%s)\n", textlen, 2 * ps->desc_max_chars, text);
    }
    else { /* INDIMODE or SIMPLEMASKMODE, sequence/mask name exceeds max, we put a '-' in it at the end of line 1 and truncate it */
      ESL_ALLOC(ps->descA[page], sizeof(char) * (max_both_lines + 2)); /* +2 so we have space for the extra '\n' (which we won't print), and the '\0' */
      /* first look to see if there's a space before desc_max_chars, this will never happen for a seq name, but may for a mask description */
      j = ps->desc_max_chars;
      while(text[j] != ' ' && j > 0) { 
	j--; 
      }
      if(j == 0) { j = ps->desc_max_chars; } /* no space */
      for(i = 0; i < j; i++) ps->descA[page][i] = text[i];
      i = j;
      if(j == ps->desc_max_chars && text[j] != ' ') { 
	ps->descA[page][i] = '-';
      }
      i++; 
      ps->descA[page][i] = '\n';
      for(i = j; i < j + ps->desc_max_chars; i++) ps->descA[page][i+2] = text[i];
      ps->descA[page][i+2] = '\0';
  }
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "add_page_desc_to_sspostscript() error, probably out of memory.");
}


/* Function: add_diffmask_page_desc_to_sspostscript()
 * 
 * Purpose:  Add text describing a diff mask page of a postscript object.
 * Throws:   Exception if the text is too long.
 */
int
add_diffmask_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *mask1, char *mask2, char *errbuf)
{
  int status;
  int i;
  char *mask1desc = NULL;
  char *mask2desc = NULL;
  int len2copy;
  int mask1len;
  int mask2len;
  void *tmp;

  if(ps->descA[page] != NULL) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), description for page %d already exists!\n", page); 
  if(mask1 == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), passed in mask1 is NULL!\n"); 
  if(mask2 == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), passed in mask2 is NULL!\n"); 

  /* check to see if we can fit the text onto two lines of max width ps->desc_max_chars */
  mask1len = (int) strlen(mask1);
  mask2len = (int) strlen(mask2);
  if((status = esl_strcat(&(mask1desc), -1, "mask 1: ", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((mask1len + 8) <= ps->desc_max_chars) { /* this will fit on one line */
    if((status = esl_strcat(&(mask1desc), -1, mask1, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  }
  else { /* won't fit on one line, include as much as we can */
    len2copy = ps->desc_max_chars - 8 - 3; /* 8 is for "mask 1: " at beginning, 3 is for "..." at end */
    ESL_RALLOC(mask1desc, tmp, sizeof(char) * ps->desc_max_chars+1);
    for(i = 0; i < len2copy; i++) { 
      mask1desc[8+i] = mask1[i];
    }
    mask1desc[8+len2copy]   = '.';
    mask1desc[8+len2copy+1] = '.';
    mask1desc[8+len2copy+2] = '.';
    mask1desc[8+len2copy+3] = '\0'; 
  }

  /* repeat for mask 2 */
  if((status = esl_strcat(&(mask2desc), -1, "mask 2: ", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((mask2len + 8) <= ps->desc_max_chars) { /* this will fit on one line */
    if((status = esl_strcat(&(mask2desc), -1, mask2, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  }
  else { /* won't fit on one line, include as much as we can */
    len2copy = ps->desc_max_chars - 8 - 3; /* 8 is for "mask 1: " at beginning, 3 is for "..." at end */
    ESL_RALLOC(mask2desc, tmp, sizeof(char) * ps->desc_max_chars+1);
    for(i = 0; i < len2copy; i++) { 
      mask2desc[8+i] = mask2[i];
    }
    mask2desc[8+len2copy]   = '.';
    mask2desc[8+len2copy+1] = '.';
    mask2desc[8+len2copy+2] = '.';
    mask2desc[8+len2copy+3] = '\0'; 
  }

  /* concatenate them in ps->descA[page] */
  if((status = esl_strcat(&(ps->descA[page]), -1, mask1desc, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((status = esl_strcat(&(ps->descA[page]), -1, "\n", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((status = esl_strcat(&(ps->descA[page]), -1, mask2desc, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");

  free(mask1desc);
  free(mask2desc);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "add_page_desc_to_sspostscript() error, probably out of memory.");
}

/* Function: add_mask_to_sspostscript
 * 
 * Purpose:  Add a mask to a sspostscript object.
 */
int
add_mask_to_ss_postscript(SSPostscript_t *ps, char *mask)
{
  int status;
  if(ps->mask != NULL) { esl_fatal("add_mask_to_ss_postscript(), mask is non-null!\n"); }
  if(mask == NULL) esl_fatal("add_mask_to_ss_postscript(), passed in mask is NULL!\n"); 
  if((status = esl_strdup(mask, -1, &(ps->mask))) != eslOK) esl_fatal("add_mask_to_ss_postscript(), error copying mask");
  return eslOK;
}


/* Function: draw_legend_column_headers()
 * 
 * Purpose:  Draw the legend column headers.
 * Return:   eslOK
 */
int 
draw_legend_column_headers(FILE *fp, SSPostscript_t *ps, char *errbuf)
{
  int status;
  int i;
  float x, y, legend_fontsize;
  char *cur_string;
  int cur_width = 0;

  legend_fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;

  x = ps->legx;
  y = ps->cur_legy;
  if(ps->mask != NULL) { 
    y -= 0.625 * LEG_BOXSIZE;
  }
  /*fprintf(fp, "/%s findfont %f scalefont setfont\n", SPECIAL_FONT, legend_fontsize);*/
  fprintf(fp, "(%s) %.4f %.4f moveto show\n", "LEGEND", x, (y + (LEG_BOXSIZE * .25)));
  /*fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, legend_fontsize);*/

  x = ps->legx_stats;
  y = ps->cur_legy;
  cur_width = ps->legx_max_chars - LEG_EXTRA_COLUMNS -2;

  ESL_ALLOC(cur_string, sizeof(char) * (cur_width+1));
  for(i = 0; i < cur_width; i++) cur_string[i] = '-'; 
  cur_string[cur_width] = '\0';

  if(ps->mask != NULL) { 
    fprintf(fp, "(%4s  %4s) %.4f %.4f moveto show\n", "", " in ", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%4s  %4s) %.4f %.4f moveto show\n", "all", "mask", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", cur_string, ps->legx, (y + (LEG_BOXSIZE * .25)));
    fprintf(fp, "(----  ----) %.4f %.4f moveto show\n", x, (y + (LEG_BOXSIZE * .25)));
  }
  else { 
    fprintf(fp, "(%5s) %.4f %.4f moveto show\n", "count", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", cur_string, ps->legx, (y + (LEG_BOXSIZE * .25)));
    fprintf(fp, "(-----) %.4f %.4f moveto show\n", x, (y + (LEG_BOXSIZE * .25)));
  }
  ps->cur_legy = y - (1.0 * LEG_BOXSIZE);
  
  free(cur_string);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "ERROR drawing legend column headers, probably out of memory.");
}

/* Function: draw_onecell_colorlegend()
 * 
 * Purpose:  Print a one cell color legend to an open file.
 * Return:   eslOK
 */
int 
draw_onecell_colorlegend(FILE *fp, OneCellColorLegend_t *occl, SSPostscript_t *ps, int occl_idx)
{
  float x, y;
  int cp;
  float fontsize;

  /* object is valid, print it */
  x = ps->legx;
  y = ps->cur_legy;

  fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;

  /* print cell */
  fprintf(fp, "newpath\n");
  fprintf(fp, "  %.2f %.2f moveto", x, y);
  fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
  fprintf(fp, "  ");
  for(cp = 0; cp < NCMYK; cp++) fprintf(fp, "%.4f ", occl->col[cp]);
  fprintf(fp, "setcmykcolor\n");
  fprintf(fp, "  fill\n");
  
  x += LEG_BOXSIZE * 1.5;

  /* print text for this legend */
  if(occl->text != NULL) { 
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, fontsize);
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", occl->text, x, y + (LEG_BOXSIZE * 0.25));

    /* print stats */
    x = ps->legx_stats;
    if(ps->mask != NULL) { 
      fprintf(fp, "(%4d  %4d) %.4f %.4f moveto show\n", occl->nres, occl->nres_masked, x, y + (LEG_BOXSIZE * 0.25));
    }
    else { 
      fprintf(fp, "(%5d) %.4f %.4f moveto show\n", occl->nres, x, y + (LEG_BOXSIZE * 0.25));
    }
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  y -= LEG_BOXSIZE * 1.5;
  ps->cur_legy = y;
  
  return eslOK;
}


/* Function: draw_scheme_colorlegend()
 * 
 * Purpose:  Print a scheme color legend to an open file.
 * Return:   eslOK
 */
int 
draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page)
{
  float x, y;
  int cp;
  int c,i,n1s;
  float fontsize;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  int do_mask;
  float old_x;

  do_mask = (ps->mask == NULL) ? FALSE : TRUE;
  do_border = (!esl_opt_GetBoolean(go, "--mask-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "--mask-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "--mask-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  x = ps->legx;
  y = ps->cur_legy;
  //y = ps->legy - (ps->nocclA[page] * (LEG_BOXSIZE * 1.5));
  fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;
  fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, fontsize);
  fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");


  float colvec[NCMYK];
  colvec[0] = colvec[1] = colvec[2] = 0.;
  colvec[3] = 1.0;
  if(do_mask) { /* print cells showing difference between masked and unmasked */
    /*x -= LEG_BOXSIZE;*/
    /*y -= LEG_BOXSIZE;*/
    fprintf(fp, "%.1f setlinewidth\n", LEG_BOXSIZE/4.);
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", colvec[cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.625;
    fprintf(fp, "(included by mask) %.4f %.4f moveto show\n", x, y);
    y -= LEG_BOXSIZE * 0.625;
    fprintf(fp, "((all colors)) %.4f %.4f moveto show\n", x, y);
    x -= LEG_BOXSIZE * 1.5;

    /* print stats for included by mask */
    old_x = x;
    n1s = 0;
    for(i = 0; i < ps->clen; i++) if(ps->mask[i] == '1') n1s++; 
    x = ps->legx_stats;
    y += LEG_BOXSIZE * 0.3125;
    fprintf(fp, "(%4s  %4d) %.4f %.4f moveto show\n", "-", n1s, x, y);
    y -= LEG_BOXSIZE * 0.3125;

    x = old_x;
    y -= LEG_BOXSIZE * 1.5;
    draw_masked_block(fp, x, y, colvec, do_circle_mask, do_square_mask, do_x_mask, do_border, LEG_BOXSIZE);

    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.625;
    fprintf(fp, "(excluded by mask) %.4f %.4f moveto show\n", x, y);
    y -= LEG_BOXSIZE * 0.625;
    fprintf(fp, "((all colors)) %.4f %.4f moveto show\n", x, y);

    /* print stats for excluded by mask */
    old_x = x;
    x = ps->legx_stats;
    y += LEG_BOXSIZE * 0.3125;
    fprintf(fp, "(%4s  %4d) %.4f %.4f moveto show\n", "-", ps->clen-n1s, x, y);

    y -= LEG_BOXSIZE * 1.8125;
    x = ps->legx;
  }

  /* print text for this legend */
  if(scl->text1 != NULL) { 
    if(scl->text2 == NULL) { 
      fprintf(fp, "(%s:) %.4f %.4f moveto show\n", scl->text1, x, (y + (LEG_BOXSIZE * .25)));
    }
    else { 
      fprintf(fp, "(%s) %.4f %.4f moveto show\n", scl->text1, x, (y + (LEG_BOXSIZE * .25)));
      y -= LEG_BOXSIZE * 0.625;
      fprintf(fp, "(%s:) %.4f %.4f moveto show\n", scl->text2, x, (y + (LEG_BOXSIZE * .25)));
    }
  }
  y -= LEG_BOXSIZE;
  
  /* print masked scheme color cells */
  /*if(do_mask) { 
    fprintf(fp, "%.1f setlinewidth\n", LEG_BOXSIZE/4.);
    for(c = 0; c < scl->nbins; c++) { 
    draw_masked_block(fp, x, y, hc_scheme[c], do_circle_mask, do_square_mask, do_x_mask, do_border, LEG_BOXSIZE);
    y -= LEG_BOXSIZE;
    }
    y += (LEG_BOXSIZE * scl->nbins);
    x += 1.5 * LEG_BOXSIZE;
    fprintf(fp, "1.0 setlinewidth\n");
  }
  */

  /* print scheme color cells and labels next to them */
  for(c = 0; c < scl->nbins; c++) { 
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", hc_scheme[c][cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.25;
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    if(c == scl->nbins-1) fprintf(fp, "(\\[%.3f-%.3f\\]) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);
    else                  fprintf(fp, "(\\[%.3f-%.3f\\)) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);

    /* print stats */
    old_x = x;
    x = ps->legx_stats;
    if(ps->mask != NULL) { 
      fprintf(fp, "(%4d  %4d) %.4f %.4f moveto show\n", scl->counts[c], scl->counts_masked[c], x, y);
    }
    else { 
      fprintf(fp, "(%5d) %.4f %.4f moveto show\n", scl->counts[c], x, y);
    }

    x = old_x - LEG_BOXSIZE * 1.5;
    y -= LEG_BOXSIZE * 0.25;
    y -= LEG_BOXSIZE;
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  
  ps->cur_legy = y;
  return eslOK;
}

/* Function: draw_sspostscript()
 * 
 * Purpose:  Print a SS postscript data structure.
 * Return:   eslOK on success;
 *           eslEINCOMPAT if ps->npage == 0
 */
int
draw_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps)
{
  int status;
  int p, pi, i, c, l, si;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  int *page_orderA;
  int rfoffset;

  do_border = (!esl_opt_GetBoolean(go, "--mask-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "--mask-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "--mask-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  if(ps->npage == 0) ESL_FAIL(eslEINCOMPAT, errbuf, "draw_sspostscript, ps->npage == 0\n");

  /* determine print order or pages, this is just 0..npage-1 UNLESS:
   *  --indi and --all and --prob all enabled, in which case, we put the 
   *    posteriors immediately after the sequences
   */
  ESL_ALLOC(page_orderA, sizeof(int) * ps->npage);
  if(esl_opt_GetBoolean(go, "--indi") && 
     esl_opt_GetBoolean(go, "--all") && 
     esl_opt_GetBoolean(go, "--prob")) { 
    pi = 0;
    rfoffset = 0;
    if(! esl_opt_GetBoolean(go, "-q")) { 
      page_orderA[0] = 0; /* print consensus sequence first */ 
      pi++;
      rfoffset = 1;
    } /* otherwise, we didn't print the consensus sequence */
    for(si = 0; si < ps->msa->nseq; si++) { 
      page_orderA[pi] = si + rfoffset; /* the indi sequence page */
      pi++;
      page_orderA[pi] = (si+ps->msa->nseq) + rfoffset; /* the posterior page */
      pi++;
    }
  }
  else { 
    for(pi = 0; pi < ps->npage; pi++) page_orderA[pi] = pi;
  }

  /* draw the pages */
  for(pi = 0; pi < ps->npage; pi++) { 
    p = page_orderA[pi];
    ps->cur_legy = ps->legy;

    /* scale section */
    fprintf(fp, "%% begin scale\n");
    fprintf(fp, "%.2f %.2f scale\n", ps->scale, ps->scale);
    fprintf(fp, "%% end scale\n\n");
      
    /* header section */
    if((status = draw_header_and_footer(fp, go, errbuf, ps, p, pi+1)) != eslOK) return status;

    /* regurgitated section */
    if(ps->regurgA != NULL) {
      fprintf(fp, "%% begin regurgitate\n");
      for(i = 0; i < ps->nregurg; i++)
	fprintf(fp, "%s", ps->regurgA[i]);
      fprintf(fp, "%% end regurgitate\n\n");
    }

    /* 'text hundreds' section */
    for(i = 0; i < ps->nhundreds; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin text hundreds\n");
	fprintf(fp, "/%s findfont %.2f scalefont setfont\n", HUNDREDS_FONT, HUNDREDS_FONTSIZE);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "(%d) %.2f %.2f moveto show\n", (i+1) * 100, ps->hundredsxA[i], ps->hundredsyA[i]); 
      if(i == (ps->nhundreds-1)) { 
	fprintf(fp, "%% end text hundreds\n\n");
      }
    }
    
    /* 'lines ticks' section */
    for(i = 0; i < ps->nticks; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin lines ticks\n");
	fprintf(fp, "%.2f setlinewidth\n", TICKS_LINEWIDTH);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", ps->ticksx1A[i], ps->ticksy1A[i], ps->ticksx2A[i], ps->ticksy2A[i]);
      if(i == (ps->nticks-1)) { 
	fprintf(fp, "%% end lines ticks\n\n");
      }
    }

    /* 'lines bpconnects' section */
    for(i = 0; i < ps->nbp; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin lines bpconnects\n");
	fprintf(fp, "%.2f setlinewidth\n", BP_LINEWIDTH);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", ps->bpx1A[i], ps->bpy1A[i], ps->bpx2A[i], ps->bpy2A[i]);
      if(i == (ps->nbp-1)) { 
	fprintf(fp, "%% end lines bpconnects\n\n");
      }
    }

    /* 'text residues' section */
    /* NOTE: I only print this out so that this file could possibly be used as a template */
    fprintf(fp, "%% begin text residues\n");
    fprintf(fp, "/%s findfont %.2f scalefont setfont\n", RESIDUES_FONT, RESIDUES_FONTSIZE);
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
    for(i = 0; i < ps->clen; i++) { 
      fprintf(fp, "() %.2f %.2f moveto show\n", ps->rxA[i], ps->ryA[i]);
    }
    fprintf(fp, "%% end text residues\n");

    /* the rest of the text will be ignored by ssu-draw if the output
     * file we're creating is read in as a template file later on
     */
    fprintf(fp, "%% begin ignore\n");
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* set to black */
    fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, LEG_FONTSIZE_UNSCALED / ps->scale);

    /* draw legend headers, if we have a legend */
    if((ps->nocclA[p] > 0) || (ps->sclAA != NULL && ps->sclAA[p] != NULL)) { 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	if((status = draw_legend_column_headers(fp, ps, errbuf)) != eslOK) return status;
      }
    }

    /* print one cell color legends, if any */
    if(ps->occlAAA != NULL && ps->occlAAA[p] != NULL) { 
      for(l = 0; l < ps->nocclA[p]; l++) 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	draw_onecell_colorlegend(fp, ps->occlAAA[p][l], ps, l);
      }
    }
    /* print scheme color legends, if any */
    if(ps->sclAA != NULL && ps->sclAA[p] != NULL) { 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	draw_scheme_colorlegend(go, fp, ps->sclAA[p], hc_scheme[ps->sclAA[p]->scheme], ps, p);
      }
    }

    if(ps->rcolAAA != NULL && ps->rcolAAA[p] != NULL) { 
      if(ps->mask != NULL) { 
	fprintf(fp, "2.0 setlinewidth\n");
	if(do_border && do_x_mask)      { fprintf(fp, "1.0 setlinewidth\n"); }
	if(do_border && do_square_mask) { fprintf(fp, "2.0 setlinewidth\n"); }
	if(do_border && do_circle_mask) { fprintf(fp, "2.5 setlinewidth\n"); }
	for(c = 0; c < ps->clen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  if(ps->mask[c] == '0') { 
	    draw_masked_block(fp, ps->rxA[c]-1., ps->ryA[c]-1., ps->rcolAAA[p][c], do_circle_mask, do_square_mask, do_x_mask, do_border, SS_BOXSIZE);
	  }
	  else { /* cell is within mask, ps->mask[c] == '1' */
	    fprintf(fp, "newpath\n");
	    fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1., ps->ryA[c] -1.);
	    fprintf(fp, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto closepath\n");
	    fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);
	    fprintf(fp, "  fill\n"); 
	  }
	}
	fprintf(fp, "1.00 setlinewidth\n");
      }
      else { /* no mask, all cells are printed the same */
	for(c = 0; c < ps->clen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  fprintf(fp, "newpath\n");
	  fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1., ps->ryA[c] -1.);
	  fprintf(fp, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto closepath\n");
	  /*fprintf(fp, "  %.2f %.2f %.2f %.2f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);*/
	  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);
	  fprintf(fp, "  fill\n");
	}
      }
      /* back to black */
      fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    }

    if(ps->rrAA[p] != NULL) { 
      fprintf(fp, "/%s findfont %f scalefont setfont\n", RESIDUES_FONT, RESIDUES_FONTSIZE);
      for(c = 0; c < ps->clen; c++) { 
	fprintf(fp, "(%c) %.2f %.2f moveto show\n", ps->rrAA[p][c], ps->rxA[c], ps->ryA[c]);
      }
    }
    fprintf(fp, "grestore\nshowpage\n");
    fprintf(fp, "%% end ignore\n\n");
  }
  free(page_orderA);
  return eslOK;

 ERROR: ESL_FAIL(eslEINVAL, errbuf, "draw_sspostscript() error, probably out of memory.");
}


/* parse_template_file
 *
 * Read secondary structure templates from a postscript 
 * template file derived from the Gutell CRW website until
 * the one that corresponds to our alignment is found, 
 * (defined as the first template structure that has the
 * same consensus length as the passed in <msa_clen>)
 * or we run out structure templates. This function
 * repeatedly calls parse_template() which actually parses
 * the templates.
 * 
 * If we run out of structure templates of anything is 
 * invalid in any of the templates in the file, we 
 * return a non-eslOK status code to caller and fill
 * errbuf with the error message.
 * 
 * Returns:  eslOK on success.
 */
int
parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, int msa_clen, SSPostscript_t **ret_ps)
{
  int             status;
  ESL_FILEPARSER *efp;
  SSPostscript_t *ps;
  int             found_match = FALSE;

  if (esl_fileparser_Open(filename, &efp) != eslOK) esl_fatal("ERROR, failed to open template file %s in parse_template_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');

  status = eslOK;
  while((found_match == FALSE) && (status == eslOK)) {
    status = parse_template_page(efp, go, errbuf, &ps);
    if((status != eslOK) && (status != eslEOF)) { 
      esl_fileparser_Close(efp);
      return status;
    }
    if(ps->clen == msa_clen) { found_match = TRUE; }
    else                     { free_sspostscript(ps); }
  }
  if(found_match == FALSE) { 
    esl_fileparser_Close(efp);
    esl_fatal("ERROR, did not find template structure to match alignment consensus length of %d in:\n%s\n", msa_clen, filename);
  }

  /* if we get here, we've found a match */

  /* validate the template we just read */
  if((status = validate_justread_sspostscript(ps, errbuf)) != eslOK) return status;

  esl_fileparser_Close(efp);
  *ret_ps = ps;

  return eslOK;
}
  

/* parse_template_page
 *
 * Read a single secondary structure template page for a
 * single CM from a postscript template file derived from the 
 * Gutell CRW website. This function is called repeatedly
 * on the same file to read multiple structure templates.
 * The logic is to keep reading until we find one with the
 * same consensus length as our input alignment. 
 * 
 * The structure template is read in sections.
 * Each section begins with a line like this: 
 * % begin <type1> <type2> 
 * 
 * list of valid tokens for <type1>:
 * modelname
 * scale
 * regurgitate
 * ignore 
 * lines
 * text
 * 
 * if <type1> is lines or text, then <type2> is read, 
 * valid tokens for <type2> if <type1> is 'text'
 * hundreds
 * residues
 * 
 * valid tokens for <type2> if <type1> is 'lines'
 * ticks
 * bpconnects
 * 
 * The 'regurgitate' lines are stored, but never changed.
 * The 'ignore' lines are not even stored.
 * All other <type1> lines are stored in data structures that
 * can be manipulated (though not all are at this point).
 * 
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_template_page(ESL_FILEPARSER *efp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps)
{
  int             status;
  char           *tok;
  int             toklen;
  SSPostscript_t *ps = NULL;
  int            read_showpage = FALSE;
  int            reached_eof = FALSE;

  /* Create the postscript object */
  ps = create_sspostscript();

  while ((read_showpage == FALSE) && ((status = esl_fileparser_GetToken(efp, &tok, &toklen))  == eslOK)) {
    if(strcmp(tok, "%") == 0) { 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
	if(strcmp(tok, "begin") == 0) { 
	  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
	    if(strcmp(tok, "modelname") == 0) { 
	      if((status = parse_modelname_section(efp, errbuf, ps)) != eslOK) return status;
	    }
	    else if(strcmp(tok, "scale") == 0) { 
	      /*printf("parsing scale\n");*/
	      if((status = parse_scale_section(efp, errbuf, ps)) != eslOK) return status;
	    }
	    else if(strcmp(tok, "ignore") == 0) { 
	      /*printf("parsing ignore\n");*/
	      if((status = parse_ignore_section(efp, errbuf, &read_showpage)) != eslOK) return status;
	    }	
	    else if(strcmp(tok, "regurgitate") == 0) { 
	      /*printf("parsing regurgitate\n");*/
	      if((status = parse_regurgitate_section(efp, errbuf, ps)) != eslOK) return status;
	    }	
	    else if(strcmp(tok, "text") == 0) { 
	      /*printf("parsing text\n");*/
	      if((status = parse_text_section(efp, errbuf, ps)) != eslOK) return status;		   
	    }
	    else if(strcmp(tok, "lines") == 0) { 
	      /*printf("parsing lines\n");*/
	      if((status = parse_lines_section(efp, errbuf, ps)) != eslOK) return status;		   
	    }	
	    else { 
	      ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), error, unknown section type %s.", tok);
	    }
	  }
	  else { 
	    ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), error last read line number %d.", efp->linenumber);
	  }
	}
	else { 
	  ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), expected line beginning with %% begin, but read tok: %s instead of begin, last read line number %d.", tok, efp->linenumber);
	}
      }
      else { 
	ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), ran out of tokens early, error last read line number %d.", efp->linenumber);
      }
    }
    else { 
      ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), expected line beginning with %%, read tok: %s, last read line number %d.", tok, efp->linenumber);
    }
  }
  if(read_showpage == FALSE && status != eslEOF) { 
    ESL_FAIL(status, errbuf, "parse_template_page(), error, ran out of tokens, but not at end of file?, last read line number %d.", efp->linenumber);
  }
  if(status == eslEOF) reached_eof = TRUE;

  *ret_ps = ps;

  if(reached_eof) return eslEOF;
  else            return eslOK;
}

/* parse_modelname_section
 *
 * Parse the modelname section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_modelname_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  char *curstr = NULL;
  int   curlen = 0;
  int  ntok = 0;
  /* this section should be exactly 3 lines, one of which we've already read,
   * here's an example, next token should be the first 0.65
   * % begin modelname
   * % archaea
   * % end scale
   */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, middle line token 1 should be a percent sign but it's %s", tok); 
  /* read remainder of line, this is the model name */
  while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen))  == eslOK) { 
    if(ntok > 0) {
      if((status = esl_strcat(&curstr, curlen, " ",  1))      != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
      curlen += 1;
    }
    if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
    curlen += toklen;
    ntok++;
  }
  ESL_ALLOC(ps->modelname, sizeof(char) * (curlen+1));
  strcpy(ps->modelname, curstr);

  /* next line should be '% end modelname' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 1 of 3 should be a percent sign but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 2 of 3"); 
  if (strcmp(tok, "end") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 2 of 3 should be 'end' but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 3 of 3"); 
  if (strcmp(tok, "modelname") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 3 of 3 should be 'modelname' but it's %s", tok); 

  if(curstr != NULL) free(curstr);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Error, parsing modelname section, memory error?");
}


/* parse_scale_section
 *
 * Parse the scale section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_scale_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;

  /* this section should be exactly 3 lines, one of which we've already read,
   * here's an example, next token should be the first 0.65
   * % begin scale
   * 0.65 0.65 scale
   * % end scale
   */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 1 of 3"); 
  ps->scale = atof(tok);
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 2 of 3"); 
  if(esl_FCompare(ps->scale, atof(tok), eslSMALLX1) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, x and y scales are not equal %.2f != %.2f", ps->scale, atof(tok)); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 3 of 3"); 
  if (strcmp(tok, "scale") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, token 3 of 3 should be 'scale' but it's %s", tok); 

  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 1 of 3 should be a percent sign but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 2 of 3"); 
  if (strcmp(tok, "end") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 2 of 3 should be 'end' but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 3 of 3"); 
  if (strcmp(tok, "scale") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 3 of 3 should be 'scale' but it's %s", tok); 

  return eslOK;
}


/* parse_ignore_section
 *
 * Parse an ignore section of a template postscript file.
 * We ignore this data. This function's purpose is to read
 * tokens until we see the "% end ignore" line signalling
 * the end of the ignore section.
 * 
 * Returns:  eslOK on success.
 */
int
parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf, int *ret_read_showpage)
{
  int status;
  char *tok;
  int   toklen;
  int   keep_reading = TRUE;
  int   read_showpage = FALSE;
  while((keep_reading) && (status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
    if (strcmp(tok, "%") == 0) { 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "ignore") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      keep_reading = FALSE;
      status = eslOK;
    }
    else if(strcmp(tok, "showpage") == 0) { 
      read_showpage = TRUE;
    }
  }
  if(status == eslEOF) ESL_FAIL(status, errbuf, "Error, parsing ignore section, finished file looking for '%% end ignore' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing ignore section, last line number read %d", efp->linenumber);

  *ret_read_showpage = read_showpage;
  return eslOK;
}


/* parse_regurgitate_section
 *
 * Parse a regurgitate section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_regurgitate_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  char *curstr = NULL;
  char *newstr;
  int   nalloc = ps->nregurg;
  int   curlen;
  void *tmp;
  int   ntok = 0;

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    curlen = ntok = 0;
    while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen))  == eslOK) { 
      if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
	if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if (strcmp(tok, "regurgitate") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	seen_end = TRUE;
	break;
      }
      else { 
	if(ntok > 0) { 
	  if((status = esl_strcat(&curstr, curlen, " ",  1))  != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (2).");
	  curlen += 1;
	}
	if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (1).");
	curlen += toklen;
	ntok++;
      }
    }
    if(seen_end) break;
    if((status = esl_strcat(&curstr, curlen, "\n", 1)) != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (3).");
    curlen += 1;
    ESL_ALLOC(newstr, sizeof(char) * (curlen+1));
    strcpy(newstr, curstr);
    if(ps->nregurg == nalloc) { 
      nalloc += ps->nalloc; ESL_RALLOC(ps->regurgA, tmp, sizeof(char *) * nalloc); 
    }
    ps->regurgA[ps->nregurg++] = newstr;
    free(curstr);
    curstr = NULL;
  }
  if(status == eslEOF) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, finished file looking for '%% end regurgitate' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, last line number read %d", efp->linenumber);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing regurgitate section");
}


/* parse_text_section
 *
 * Parse a text section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_text_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  int   nalloc;
  void *tmp;
  int do_hundreds = FALSE;
  int do_residues = FALSE;

  /* find out which section we're in, 'hundreds' or 'residues' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if      (strcmp(tok, "hundreds") == 0) { do_hundreds = TRUE; nalloc = ps->nhundreds; }
  else if (strcmp(tok, "residues") == 0) { do_residues = TRUE; nalloc = ps->clen; }

  /* read the first two special lines, should be a 5-token line ending with setfont and a 5-token line ending with setcmykcolor,
   * we don't store these, but we require that they're there. */
  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if(strcmp(tok, "setfont") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");

  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if(strcmp(tok, "setcmykcolor") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if (strcmp(tok, "text") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if(do_hundreds) { 
	if (strcmp(tok, "hundreds") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text hundreds' after it"); 
      }
      if(do_residues) {
	if (strcmp(tok, "residues") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text residues' after it"); 
      }
      seen_end = TRUE;
      break;
    }
    /* if we get here, we haven't seen the end, we're reading a normal line, tok is the string, we discard this */
    if(do_hundreds && ps->nhundreds == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->hundredsxA, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->hundredsyA, tmp, sizeof(float) * nalloc); 
    }
    if(do_residues && ps->clen == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->rxA, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ryA, tmp, sizeof(float) * nalloc); 
    }
    /* get x */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if(do_hundreds) ps->hundredsxA[ps->nhundreds] = atof(tok);
    if(do_residues) ps->rxA[ps->clen] = atof(tok);
    /* get y */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if(do_hundreds) ps->hundredsyA[ps->nhundreds] = atof(tok);
    if(do_residues) ps->ryA[ps->clen] = atof(tok);
    
    /* verify moveto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "moveto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text main section, fourth token should be 'moveto', line %d", efp->linenumber);
    /* verify show */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "show") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text main section, fifth token should be 'show', line %d", efp->linenumber);
    
    if(do_hundreds) ps->nhundreds++;
    if(do_residues) ps->clen++;
  }
  if(!seen_end) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, didn't see end! line: %d\n", efp->linenumber);
  if(status == eslEOF && do_hundreds) 
    ESL_FAIL(status, errbuf, "Error, parsing text section, finished file looking for '%% end text hundreds' line");
  if(status == eslEOF && do_residues) 
    ESL_FAIL(status, errbuf, "Error, parsing text section, finished file looking for '%% end text residues' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing text section, last line number read %d", efp->linenumber);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing text section");
}

/* parse_lines_section
 *
 * Parse a lines section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_lines_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  int   nalloc;
  void *tmp;
  int do_ticks = FALSE;
  int do_bpconnects = FALSE;

  /* find out which section we're in, 'ticks' or 'bpconnects' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if      (strcmp(tok, "ticks") == 0)      { do_ticks = TRUE; nalloc = ps->nticks; }
  else if (strcmp(tok, "bpconnects") == 0) { do_bpconnects = TRUE; nalloc = ps->nbp; }

  /* read the first two special lines, should be a 2-token line ending with setlinewidth and a 5-token line ending with setcmykcolor,
   * we don't store these, but we require that they're there. */
  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");
  if(strcmp(tok, "setlinewidth") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");

  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if(strcmp(tok, "setcmykcolor") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if (strcmp(tok, "lines") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if(do_ticks) { 
	if (strcmp(tok, "ticks") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines ticks' after it"); 
      }
      if(do_bpconnects) {
	if (strcmp(tok, "bpconnects") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines bpconnects' after it"); 
      }
      seen_end = TRUE;
      break;
    }
    /* if we get here, we haven't seen the end, we're reading a normal line, tok is the first x coord, we record this */
    /* first we expand our arrays if nec */
    if(do_ticks && ps->nticks == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->ticksx1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksy1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksx2A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksy2A, tmp, sizeof(float) * nalloc); 
    }
    if(do_bpconnects && ps->nbp == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->bpx1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpy1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpx2A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpy2A, tmp, sizeof(float) * nalloc); 
    }
    /* store x1 */
    if(do_ticks) ps->ticksx1A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpx1A[ps->nbp] = atof(tok);
    /* get y1 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksy1A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpy1A[ps->nbp] = atof(tok);
    /* get x2 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksx2A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpx2A[ps->nbp] = atof(tok);
    /* get y2 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksy2A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpy2A[ps->nbp] = atof(tok);
    
    /* verify newpath */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "newpath") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, fifth token should be 'newpath', line %d", efp->linenumber);
    /* verify moveto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "moveto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, sixth token should be 'moveto', line %d", efp->linenumber);
    /* verify lineto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "lineto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, seventh token should be 'lineto', line %d", efp->linenumber);
    /* verify stroke */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "stroke") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, eigth token should be 'stroke', line %d", efp->linenumber);
    
    if(do_ticks) ps->nticks++;
    if(do_bpconnects) ps->nbp++;
  }
  if(!seen_end) ESL_FAIL(status, errbuf, "Error, parsing lines section, didn't see end! line: %d\n", efp->linenumber);
  if(status == eslEOF && do_ticks) 
    ESL_FAIL(status, errbuf, "Error, parsing lines section, finished file looking for '%% end lines ticks' line");
  if(status == eslEOF && do_bpconnects) 
    ESL_FAIL(status, errbuf, "Error, parsing lines section, finished file looking for '%% end lines bpconnects' line");

  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing lines section, last line number read %d", efp->linenumber);

  return eslOK;
 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing lines section");
}

/* Function: individual_seqs_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with info for each seq in the MSA 
 * Return:   eslOK on success.
 */
int
individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa)
{
  int status;
  int p, i, pp;
  int cpos, apos;
  int orig_npage = ps->npage;

  if((status = add_pages_sspostscript(ps, msa->nseq, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
  }

  /* fill ps->rrAA with residues and gaps */
  for(i = 0; i < msa->nseq; i++) {
    pp = orig_npage + i;
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
	ps->rrAA[pp][cpos] = (char) msa->aseq[i][apos];
	/* printf("ps->rrAA[%3d][%4d]: %c\n", pp, cpos, ps->rrAA[pp][cpos]); */
	cpos++;
      }
    }
    ps->rrAA[pp][cpos] = '\0';
    ps->seqidxA[pp] = i;
    /* add description to ps */
    if((status = add_page_desc_to_sspostscript(ps, pp, msa->sqname[i], errbuf)) != eslOK) return status;
  }
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "individual_seqs_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: rf_seq_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, the RF sequence.
 * Return:   eslOK on success.
 */
int
rf_seq_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa)
{
  int status;
  int p, pp;
  int cpos, apos;
  int orig_npage = ps->npage;

  if((status = add_pages_sspostscript(ps, 1, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
  }

  /* fill ps->rrAA with residues and gaps for RF sequence */
  pp = orig_npage;
  cpos = 0;
  for(apos = 1; apos <= msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) { /* a consensus position */
      ps->rrAA[pp][cpos] = msa->rf[(apos-1)];
      /* printf("ps->rrAA[%3d][%4d]: %c\n", pp, cpos, ps->rrAA[pp][cpos]); */
      cpos++;
    }
  }

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "*CONSENSUS*", errbuf)) != eslOK) return status;

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "rf_seq_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: infocontent_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, colored squares indicating
 *           the information content of each consensus column.
 *           
 * Return:   eslOK on success.
 */
int
infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  double **obs   = NULL;
  double  *ent   = NULL;
  double  *bg    = NULL;
  float   *limits = NULL;
  int zero_obs;
  int nonecell = 0;
  int nonecell_masked = 0;
  int within_mask;

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->sclAA[p]   = NULL;
    ps->occlAAA[p] = NULL;
  }
  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->clen; c++) 
      ps->rcolAAA[p][c] = NULL;
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(ent, sizeof(double) * ps->clen);
  ESL_ALLOC(obs, sizeof(double *) * ps->clen);
  ESL_ALLOC(bg, sizeof(double) * msa->abc->K);
  esl_vec_DSet(bg, msa->abc->K, 1./(msa->abc->K));

  for(cpos = 0; cpos < ps->clen; cpos++) 
    obs[cpos] = NULL;
  
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    ESL_ALLOC(obs[cpos], sizeof(double) * msa->abc->K);
    esl_vec_DSet(obs[cpos], msa->abc->K, 0.);
  }

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.4;
  limits[2] = 0.8;
  limits[3] = 1.2;
  limits[4] = 1.6;
  limits[5] = 1.99;
  limits[6] = 2.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);

  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */
  for(i = 0; i < msa->nseq; i++) { 
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { /* seq i is not a gap at cpos */
	  esl_abc_DCount(msa->abc, obs[cpos], esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][apos]), 1.);
	}
	cpos++;
      }
    }
  }
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    zero_obs = (esl_DCompare(esl_vec_DSum(obs[cpos], msa->abc->K), 0., eslSMALLX1) == eslOK) ? TRUE : FALSE;
    esl_vec_DNorm(obs[cpos], msa->abc->K);
    ent[cpos] = esl_vec_DEntropy(bg, msa->abc->K) - esl_vec_DEntropy(obs[cpos], msa->abc->K);

    if(zero_obs) { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell++;
      if(ps->mask != NULL && ps->mask[cpos] == '1') nonecell_masked++; 
    }
    else { 
      within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent[cpos], ps->sclAA[pp], within_mask)) != eslOK) return status;
    }

    ps->rrAA[pp][cpos] = ' ';
    /*ps->rrAA[pp][cpos] = (esl_FCompare(ent[cpos], 0., eslSMALLX1) == eslOK) ? '-' : ' ';*/
  }

  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "100% gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 1;

  /* add text to legend */
  /*sprintf(text, "information content (bits) (total: %.2f bits)", esl_vec_DSum(ent, ps->clen));*/
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "information content (bits)", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "information content per position", errbuf)) != eslOK) return status;

  free(ent);
  for(cpos = 0; cpos < ps->clen; cpos++) free(obs[cpos]);
  free(obs);
  free(bg);
  free(limits);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "infocontent_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: delete_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with a new page w/colored squares indicating
 *           the number of sequences with gaps (deletions) at each consensus column. 
 *           If do_all is TRUE the page shows all deletions. If false only 'internal' deletions,
 *           those that come after the first occupied consensus column of each sequence and
 *           before the final occupied consensus column for each sequence, are shown.
 *           
 * Return:   eslOK on success.
 */
int
delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  int *dct;
  int *dct_internal;
  int *fA, *lA;
  int nonecell = 0;
  int nonecell_masked = 0;
  int within_mask;
  float *limits = NULL;

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(int *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(int) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(dct,          sizeof(int) * ps->clen);
  ESL_ALLOC(dct_internal, sizeof(int) * ps->clen);
  esl_vec_ISet(dct,          ps->clen, 0);
  esl_vec_ISet(dct_internal, ps->clen, 0);

  /* determine the first and last occupied consensus position in each sequence */
  ESL_ALLOC(fA, sizeof(int) * msa->nseq);
  ESL_ALLOC(lA, sizeof(int) * msa->nseq);
  esl_vec_ISet(lA, msa->nseq, 0);
  esl_vec_ISet(fA, msa->nseq, ps->clen-1);
  /* this could be more efficient */
  for(i = 0; i < msa->nseq; i++) { 
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* apos is a consensus position */
	cpos++;
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { /* cpos for seq i is not a gap */
	  fA[i] = ESL_MIN(fA[i], cpos);
	  lA[i] = ESL_MAX(lA[i], cpos);
	}
      }
    }
  }

  for(i = 0; i < msa->nseq; i++) { 
    /*printf("fA[%4d] %4d lA[%4d] %4d\n", i, fA[i], i, lA[i]);*/
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
	cpos++;
	if(esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { 
	  dct[(cpos-1)]++; 
	  if(cpos >= fA[i] && cpos <= lA[i]) dct_internal[(cpos-1)]++;
	}
      }
    }
  }

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);

  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */
  if(do_all) { 
    /* draw delete page with all deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if(dct[cpos] == 0) { 
	if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status; 
	nonecell++;
	if(ps->mask != NULL && ps->mask[cpos] == '1') nonecell_masked++; 
      }
      else {
	within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
	if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct[cpos]) / (float) (msa->nseq), ps->sclAA[pp], within_mask)) != eslOK) return status;
      }
    }
  }
  else { /* do_all is FALSE, draw delete page with only internal deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if(dct_internal[cpos] == 0) { 
	if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status; 
	nonecell++;
	if(ps->mask != NULL && ps->mask[cpos] == '1') nonecell_masked++; 
      }
      else {
	within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
	if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct_internal[cpos]) / (float) (msa->nseq), ps->sclAA[pp], within_mask)) != eslOK) return status;
      }
    }
  }

  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  ps->nocclA[pp] = 1;

  if(do_all) { 
    if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero deletions", ps->legx_max_chars, errbuf)) != eslOK) return status;
  }
  else {
    if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero internal deletions", ps->legx_max_chars, errbuf)) != eslOK) return status; 
  }

  /* add color legend and description */
  if(do_all) { 
    /*sprintf(text, "fraction seqs w/deletes ('-'=0 deletes; avg/seq: %.2f)", (float) esl_vec_ISum(dct, ps->clen) / (float) msa->nseq);*/
    if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs with deletes", ps->legx_max_chars, errbuf)) != eslOK) return status;
    if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of deletions at each position", errbuf)) != eslOK) return status;
  }
  else { /* !do_all, only internal deletes counted */
    /*sprintf(text, "fraction seqs w/internal deletes ('-'=0; avg/seq: %.2f)", (float) esl_vec_ISum(dct_internal, ps->clen) / (float) msa->nseq);*/
    if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs w/internal deletions", ps->legx_max_chars, errbuf)) != eslOK) return status;
    if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of internal (non-terminal) deletions in each position", errbuf)) != eslOK) return status;
  }

  free(dct);
  free(dct_internal);
  free(fA);
  free(lA);
  free(limits);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "delete_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: insert_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, with colors in log 
 *           scale indicating the fraction of seqs with inserts after each
 *           position, and numbers indicating the median length of inserts in those
 *           sequences that have inserts at each position. Positions with 0 inserts
 *           in all sequences are marked '-' with no color. Positions with median
 *           length 10 or greater are marked with '*'. 
 *           
 * Return:   eslOK on success.
 */
int
insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  int **ict;
  int *total_ict, *med_ict, *nseq_ict;
  int imed;
  float col;
  char res;
  int nonecell = 0;
  int nonecell_masked = 0;
  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(ict,  sizeof(int *) * (ps->clen+1));
  for(cpos = 0; cpos <= ps->clen; cpos++) { 
    ESL_ALLOC(ict[cpos],  sizeof(int) * (msa->nseq));
    esl_vec_ISet(ict[cpos], (msa->nseq), 0);
  }
  
  ESL_ALLOC(total_ict,  sizeof(int) * (ps->clen+1));
  ESL_ALLOC(nseq_ict,  sizeof(int) * (ps->clen+1));
  ESL_ALLOC(med_ict,  sizeof(int) * (ps->clen+1));
  esl_vec_ISet(total_ict, (ps->clen+1), 0);
  esl_vec_ISet(nseq_ict, (ps->clen+1), 0);
  esl_vec_ISet(med_ict, (ps->clen+1), 0);

  cpos = 0;
  for(apos = 0; apos < msa->alen; apos++) { 
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) cpos++;
    else { 
      for(i = 0; i < msa->nseq; i++)
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { 
	  total_ict[cpos]++;
	  ict[cpos][i]++;
	}	  
    }
  }

  int nseq;
  int *len;
  int l;
  /* determine avg median length for each insertion */
  for(cpos = 0; cpos <= ps->clen; cpos++) { 
    if(total_ict[cpos] > 0) { 
      nseq = 0;
      for(i = 0; i < msa->nseq; i++) { 
	if(ict[cpos][i] >= 1) nseq_ict[cpos]++;
      }
      ESL_ALLOC(len, sizeof(int) * nseq_ict[cpos]);
      l = 0;
      for(i = 0; i < msa->nseq; i++) { 
	if(ict[cpos][i] >= 1) { 
	  len[l++] = ict[cpos][i];
	}
      }
      qsort(len, nseq, sizeof(int), compare_ints);
      med_ict[cpos] = len[nseq / 2];
      free(len);
    }
  }

  pp = orig_npage;

  /* add color legend */
  float *limits;
  int within_mask;
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);

  for(cpos = 1; cpos <= ps->clen; cpos++) { 
    if(nseq_ict[cpos] == 0) { 
      res = '-';
      col = 0.0;
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][(cpos-1)], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell++;
      if(ps->mask != NULL && ps->mask[cpos-1] == '1') nonecell_masked++; 
    }
    else {
      imed = (int) med_ict[cpos];
      switch (imed) { 
      case 0: res = '0';
	break;
      case 1: res = '1';
	break;
      case 2: res = '2';
	break;
      case 3: res = '3';
	break;
      case 4: res = '4';
	break;
      case 5: res = '5';
	break;
      case 6: res = '6';
	break;
      case 7: res = '7';
	break;
      case 8: res = '8';
	break;
      case 9: res = '9';
	break;
      default: res = '*';
	break;
      }
      /*col = 1. / (1. - log((float) nseq_ict[cpos] / (float) msa->nseq)); */
      col = (float) nseq_ict[cpos] / (float) msa->nseq;
      within_mask = (ps->mask != NULL && ps->mask[(cpos-1)] == '1') ? TRUE : FALSE;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][(cpos-1)], NCMYK, hc_scheme[hc_scheme_idx], col, ps->sclAA[pp], within_mask)) != eslOK) return status;
    }
    /*ps->rrAA[pp][(cpos-1)] = res;*/
    ps->rrAA[pp][(cpos-1)] = ' ';
  }

  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  ps->nocclA[pp] = 1;
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero insertions", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add color legend */
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs w/insertions", ps->legx_max_chars, errbuf)) != eslOK) return status;
  if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of insertions after each position", errbuf)) != eslOK) return status;

  for(i = 0; i < ps->clen; i++) free(ict[i]);
  free(ict);
  free(total_ict);
  free(nseq_ict);
  free(med_ict);
  free(limits);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "insert_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: posteriors_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with info on posterior probabilities in the MSA.
 * Return:   eslOK on success.
 */
static int
posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int    status;
  int    s,c;           /* counters over sequences, columns of MSA */
  int   *nongap_c;      /* number of non-gap posterior values for each column */
  int   *nongaprf_c;    /* number of non-gap RF non-gap posterior values for each column */
  float *sum_c;         /* sum of non-gap posterior values for each column */
  float *sumrf_c;       /* sum of non-gap RF non-gap posterior values for each column */
  int    nongap_s;      /* number of non-gap posterior values for current sequence */
  int    nongaprf_s;    /* number of non-gap RF non-gap posterior values for current sequence */
  float  sum_s;         /* sum of non-gap posterior values for current sequence */
  float  sumrf_s;       /* sum of non-gap RF non-gap posterior values for current sequence */
  float  avgrf_c;       /* avg non-gap RF non-gap posterior value for current column */
  float  avg_s;         /* avg non-gap posterior values for current sequence */
  float  avgrf_s;       /* avg non-gap RF non-gap posterior values for current sequence */
  int    p;
  float  prob;
  int *c2a_map;
  int *a2c_map;
  int clen;
  int cpos;
  int orig_npage = ps->npage;
  int new_npage = 0;
  int ir1, ir2;
  int pp;
  int do_avg = FALSE;
  int do_indi = FALSE;
  int nfirst_indi_page = -1;
  int navg_page = -1;
  int ridx1, ridx2, r;
  float *limits; /* bin limits for the color scheme */
  int nonecell_avg = 0;
  int nonecell_avg_masked = 0;
  int nonecell_seq = 0;
  int nonecell_seq_masked = 0;
  int within_mask;

  if(ps->mask == NULL) { nonecell_avg_masked = nonecell_seq_masked = -1; } /* special flag */

  if(msa->rf == NULL) esl_fatal("No RF annotation in alignment");

  if(esl_opt_GetBoolean(go, "--indi")) { do_indi = TRUE; nfirst_indi_page = orig_npage + new_npage; new_npage += msa->nseq; }
  else                                 { do_avg = TRUE;  new_npage += 1; navg_page = orig_npage; }

  if(do_indi) {
    if((status = add_pages_sspostscript(ps, new_npage, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");
  }
  else {
    if((status = add_pages_sspostscript(ps, new_npage, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");
  }

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  /* Find out which #=GR line is the POST, Post, or post line (if more than one exist, last one is chosen) */
  ridx1 = ridx2 = -1;
  for (r = 0; r < msa->ngr; r++) { 
    if (strcmp(msa->gr_tag[r], "POSTX.") == 0) { ridx1 = r; }
    if (strcmp(msa->gr_tag[r], "POST.X") == 0) { ridx2 = r; }
  }
  if((ridx1 == -1) || (ridx2 == -1)) { 
    ESL_FAIL(eslEINVAL, errbuf, "--prob requires \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s (from infernal v1.x\'s cmalign).\n", esl_opt_GetArg(go,1));
  }
  
  /* map consensus columns to alignment positions */
  map_cpos_to_apos(msa, &c2a_map, &a2c_map, &clen);

  /* per column stats */
  ESL_ALLOC(nongap_c, sizeof(int) * msa->alen);
  ESL_ALLOC(sum_c,    sizeof(float) * msa->alen);
  ESL_ALLOC(nongaprf_c, sizeof(int) * msa->alen);
  ESL_ALLOC(sumrf_c,    sizeof(float) * msa->alen);
  esl_vec_ISet(nongap_c, msa->alen, 0);
  esl_vec_FSet(sum_c,    msa->alen, 0.);
  esl_vec_ISet(nongaprf_c, msa->alen, 0);
  esl_vec_FSet(sumrf_c,    msa->alen, 0.);

  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.8;
  limits[2] = 0.9;
  limits[3] = 0.925;
  limits[4] = 0.95;
  limits[5] = 0.975;
  limits[6] = 1.00;

  /* step through each sequence and each column, collecting stats */
  pp = nfirst_indi_page;
  for(s = 0; s < msa->nseq; s++) { 
    nonecell_seq = 0;
    nonecell_seq_masked = (ps->mask == NULL) ? -1 : 0;
    nongap_s = nongaprf_s = 0;
    sum_s    = sumrf_s = 0.;
    if(do_indi) { /* add color legend for this sequence */
      ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);
    }
    for(c = 0; c < msa->alen; c++) { 
      if(! esl_abc_CIsGap(msa->abc, msa->gr[ridx1][s][c])) {
	if(esl_abc_CIsGap(msa->abc, msa->gr[ridx2][s][c])) ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value non-gap but post 'ones' value is gap.\n", s, c);
	if(msa->gr[ridx1][s][c] == '*') {
	  if(msa->gr[ridx2][s][c] != '*') ESL_FAIL(eslEINVAL, errbuf, "reading post annotation for seq: %d aln column: %d, post 'tens' value '*' but post 'ones' value != '*'.\n", s, c);
	  prob = 1.0;
	}
	else {
	  ir1 = (int) (msa->gr[ridx1][s][c] - '0');
	  ir2 = (int) (msa->gr[ridx2][s][c] - '0');
	  prob = ((float) ir1 * 10. + ir2) * .01;
	  /*printf("c: %d r1: %c %d r2: %c %d p: %.2f\n", c, msa->gr[ridx1][s][c], ir1, msa->gr[ridx2][s][c], ir2, prob);*/
	}
	sum_c[c] += prob;
	nongap_c[c]++;
	sum_s += prob;
	nongap_s++;
	if(a2c_map[c] != -1) { /* consensus position */
	  cpos = a2c_map[c];
	  sumrf_c[c] += prob;
	  nongaprf_c[c]++;
	  sumrf_s += prob;
	  nongaprf_s++;
	  if(do_indi) { /* compute color for this column */
	    within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
	    if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], prob, ps->sclAA[pp], within_mask)) != eslOK) return status;
	    ps->rrAA[pp][cpos] = ' ';
	  }
	}
	//if(prob >= pthresh) athresh_c[c]++;
      }
      else if (do_indi) { /* gap, if it's a consensus column, draw blank square */
	if(a2c_map[c] != -1) { /* consensus position */
	  cpos = a2c_map[c];
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
	  nonecell_seq++;
	  if(ps->mask != NULL && ps->mask[cpos] == '1') nonecell_seq_masked++; 
	  ps->rrAA[pp][cpos] = ' ';
	}
      }
    } /* done with this sequence */
    if(do_indi) { 
      avg_s   =  (float) sum_s / (float) nongap_s;
      avgrf_s =  (float) sumrf_s / (float) nongaprf_s;

      /* add one-cell color legend */
      ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell_seq, nonecell_seq_masked);
      if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "gap", ps->legx_max_chars, errbuf)) != eslOK) return status;
      ps->nocclA[pp] = 1;

      if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "posterior probability (alnment confidence)", ps->legx_max_chars, errbuf)) != eslOK) return status;
      ps->seqidxA[pp] = s;
      if((status = add_page_desc_to_sspostscript(ps, pp, msa->sqname[s], errbuf)) != eslOK) return status;
      pp++;
    }
  } /* done with all sequences */

  if(do_avg) { /* add average colors */
    pp = navg_page;
    ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);
    for(c = 0; c < msa->alen; c++) { 
      if(a2c_map[c] != -1) { /* cons position */
	cpos = a2c_map[c]; 
	if(nongap_c[c] > 0) { 
	  avgrf_c = sum_c[c] /= (float) nongap_c[c];
	  within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], avgrf_c, ps->sclAA[pp], within_mask)) != eslOK) return status;
	}
	else { 
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
	  nonecell_avg++;
	  if(ps->mask != NULL && ps->mask[cpos] == '1') nonecell_avg_masked++; 
	}
	ps->rrAA[pp][cpos] = ' ';
      }
    }

    /* add one-cell color legend */
    ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell_avg, nonecell_avg_masked);
    if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "100% gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;
    ps->nocclA[pp] = 1;
    
    /* add color legend */
    if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "average posterior probability (alnment confidence)", ps->legx_max_chars,errbuf)) != eslOK) return status;

    /* add description to ps */
    if((status = add_page_desc_to_sspostscript(ps, pp, "average posterior probability (confidence) per position", errbuf)) != eslOK) return status;
  }


  free(nongap_c);
  free(nongaprf_c);
  free(sum_c);
  free(sumrf_c);
  free(c2a_map);
  free(a2c_map);
  free(limits);

  return eslOK;

 ERROR:
  return status;
}

/* Function: colormask_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page based on a lanemask, each column
 *           is either black (if included, a '1' in the mask) or pink (not included, a '0' in the
 *           mask.
 *
 * Return:   eslOK on success.
 */
int
colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float **hc_onecell, int incmask_idx, int excmask_idx)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_inside_mask = 0;
  int ncols_outside_mask = 0;
  char *mask_desc = NULL;

  if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  for(cpos = 0; cpos < ps->clen; cpos++) { 
    if(mask[cpos] == '1') { /* included */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[incmask_idx])) != eslOK) return status; 
      ncols_inside_mask++;
    }
    else if(mask[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[excmask_idx])) != eslOK) return status; 
      ncols_outside_mask++;
    }
    else ESL_FAIL(eslEINVAL, errbuf, "--mask mask char number %d is not a 1 nor a 0, but a %c\n", cpos, mask[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incmask_idx], ncols_inside_mask, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "columns included by mask", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[excmask_idx], ncols_outside_mask, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "columns excluded by mask", ps->legx_max_chars, errbuf)) != eslOK) return status;

  if((status = esl_strcat(&(mask_desc), -1, "mask file: ", -1)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");;
  if((status = esl_strcat(&(mask_desc), -1, esl_opt_GetString(go, "--mask"), -1)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");;
  if((status = add_page_desc_to_sspostscript(ps, pp, mask_desc, errbuf)) != eslOK) return status;

  ps->nocclA[pp] = 2;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "colormask_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}



/* Function: diffmask_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page based on a comparison between
 *           two masks, each column is either black (if included (a '1') in both masks), red
 *           (if a '1' in mask 1 and a '0' in mask 2), cyan (a '0' in mask 1 and a '1' in mask 2)
 *           or grey (a '0' in both masks).
 *
 * Return:   eslOK on success.
 */
int
diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2, float **hc_onecell, int incboth_idx, int inc1_idx, int inc2_idx, int excboth_idx)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_in_both = 0;
  int ncols_out_both = 0;
  int ncols_in_1_out_2 = 0;
  int ncols_out_1_in_2 = 0;

  if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->occlAAA[p],  sizeof(OneCellColorLegend_t **) * 4);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  for(cpos = 0; cpos < ps->clen; cpos++) { 
    if(mask1[cpos] == '1' && mask2[cpos] == '1') { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[incboth_idx])) != eslOK) return status; 
      ncols_in_both++;
    }
    else if(mask1[cpos] == '1' && mask2[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[inc1_idx])) != eslOK) return status; 
      ncols_in_1_out_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '1') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[inc2_idx])) != eslOK) return status; 
      ncols_out_1_in_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[excboth_idx])) != eslOK) return status; 
      ncols_out_both++;
    }
    else if(mask1[cpos] != '0' && mask1[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-col char number %d is not a 1 nor a 0, but a %c\n", cpos, mask1[cpos]);
    else if(mask2[cpos] != '0' && mask2[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-diff char number %d is not a 1 nor a 0, but a %c\n", cpos, mask2[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incboth_idx], ncols_in_both, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "included by both masks", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[inc1_idx], ncols_in_1_out_2, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "incl. mask 1, excl. mask 2", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][2] = create_onecell_colorlegend(hc_onecell[inc2_idx], ncols_out_1_in_2, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][2], "excl. mask 1, incl. mask 1", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][3] = create_onecell_colorlegend(hc_onecell[excboth_idx], ncols_out_both, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][3], "excluded by both masks", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 4;

  if((status = add_diffmask_page_desc_to_sspostscript(ps, pp, esl_opt_GetString(go, "--mask"), esl_opt_GetString(go, "--mask-diff"), errbuf)) != eslOK) return status;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "diffmask_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}
/* Function: add_pages_sspostscript()
 * 
 * Purpose:  Add and initialize blank pages to a postscript object.
 */ 
static int 
add_pages_sspostscript(SSPostscript_t *ps, int ntoadd, int page_mode)
{
  int status;
  void *tmp;
  int p;

  if(ps->npage == 0) { 
    assert(ps->rrAA    == NULL);
    assert(ps->rcolAAA == NULL);
    ESL_ALLOC(ps->rrAA,    sizeof(char *)   * (ps->npage + ntoadd));
    ESL_ALLOC(ps->rcolAAA, sizeof(float **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->occlAAA, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->nocclA,  sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->sclAA,   sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->descA,   sizeof(char *) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->modeA,   sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->seqidxA, sizeof(int) * (ps->npage + ntoadd));
  }
  else { 
    assert(ps->rrAA    != NULL);
    assert(ps->rcolAAA != NULL);
    ESL_RALLOC(ps->rrAA,   tmp, sizeof(char *)   * (ps->npage + ntoadd));
    ESL_RALLOC(ps->rcolAAA,tmp, sizeof(float **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->occlAAA,tmp, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->nocclA, tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->sclAA,  tmp, sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->descA,  tmp, sizeof(char *) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->modeA,  tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->seqidxA,tmp, sizeof(int) * (ps->npage + ntoadd));
  }
  for(p = ps->npage; p < (ps->npage + ntoadd); p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->occlAAA[p] = NULL;
    ps->nocclA[p]  = 0;
    ps->sclAA[p]   = NULL;
    ps->descA[p]   = NULL;
    ps->modeA[p]   = page_mode;
    ps->seqidxA[p] = -1;
  }
  ps->npage += ntoadd;

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

/* map_cpos_to_apos
 *                   
 * Given an MSA, determine the alignment position each
 * consensus (#=GC RF) position refers to. 
 * Both maps that are returned are indexed starting from 0.
 * c2a_map[0..clen-1]
 * a2c_map[0..alen-1]
 */
static int map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int **ret_a2c_map, int *ret_clen)
{
  int status;
  int clen = 0;
  int *c2a_map = NULL;
  int *a2c_map = NULL;
  int cpos = 0;
  int apos = 0;
  /* contract check */
  if(msa->rf == NULL) { status = eslEINVAL; goto ERROR; }

  /* count consensus columns */
  for(apos = 0; apos < msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) clen++;

  /* build map */
  ESL_ALLOC(c2a_map, sizeof(int) * clen);
  ESL_ALLOC(a2c_map, sizeof(int) * msa->alen);
  esl_vec_ISet(a2c_map, msa->alen, -1);

  for(apos = 0; apos < msa->alen; apos++) { 
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      a2c_map[apos] = cpos;
      c2a_map[cpos] = apos;
      cpos++;
    }
  }

  *ret_c2a_map = c2a_map;
  *ret_a2c_map = a2c_map;
  *ret_clen    = clen;
  return eslOK;

 ERROR:
  if(c2a_map != NULL) free(c2a_map);
  return status;
}

/* read_mask_file
 *
 * Given an open file pointer, read the first token of the
 * file and return it as *ret_mask. Also return length of
 * mask in *ret_masklen, and a flag indicating whether or
 * not the mask has any internal zeroes in *ret_mask_has_internal_zeroes.
 * An internal '0' is one that occurs between at least one
 * 5' and one 3' '1'
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen, int *ret_mask_has_internal_zeroes)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;
  /* for determining if we have internal zeroes */
  int             seen_1 = FALSE;                /* becomes TRUE when we see first (5'-most) '1' */
  int             seen_1_then_0 = FALSE;         /* becomes TRUE when we see first (5'-most) '0' that is 3' of first '1' */
  int             seen_1_then_0_then_1 = FALSE;  /* becomes TRUE when we see first (5'-most) '1' that is 3' of first '0' that is 3' of first '1' */

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));

  for(n = 0; n < toklen; n++) { 
    mask[n] = tok[n];
    if(mask[n] == '0') { 
      if((seen_1) && (!seen_1_then_0)) { seen_1_then_0 = TRUE; }
    }
    else if (mask[n] == '1') { 
      if(!seen_1) { seen_1 = TRUE; }
      if((seen_1) && (seen_1_then_0) && (!seen_1_then_0_then_1)) { seen_1_then_0_then_1 = TRUE; }
    }
    else { ESL_FAIL(eslEINVAL, errbuf, "character %d of mask file is invalid: %c (must be a '1' or a '0')\n", n, mask[n]); }

    mask[n] = tok[n];
  }
  mask[n] = '\0';

  *ret_mask = mask;
  *ret_masklen= toklen;
  *ret_mask_has_internal_zeroes = seen_1_then_0_then_1;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}

/* Function: drawfile2sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with >= 1 new page(s), with colors described
 *           in an input 'draw' file, with >= 1 sets of <x> lines of data, each set 
 *           is separated by a line with only "//". <x> must be equal to the consensus
 *           ps->clen. Each line has at least 4 floats explaining 
 *           the CMYK values for the color to use at each position of the SS diagram,
 *           and optionally contains an extra single character which is the residue
 *           to put at that position.
 *           
 * Return:   eslOK on success.
 */
int
drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps)
{
  int status;
  int p, pp;
  int cpos, c;
  int orig_npage = ps->npage;
  ESL_FILEPARSER *efp;
  char           *s;
  char *dfile = esl_opt_GetString(go, "--dfile");
  if (esl_fileparser_Open(dfile, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in draw_file2sspostscript\n", dfile);
  esl_fileparser_SetCommentChar(efp, '#');

  pp = orig_npage - 1;
  cpos = 0;

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      /* example line without residue markup:
       * 0.000 0.000 0.000 0.500
       *
       * example line with residue markup:
       * 0.000 0.000 0.000 0.500 A
       */

      cpos++;
      if(cpos == 1) { /* add a new page */
	if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");
	
	for(p = (ps->npage-1); p < ps->npage; p++) { 
	  ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
	  ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
	  for(c = 0; c < ps->clen; c++) { 
	    ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
	  }
	}
	pp++; /* if first page, pp == orig_npage now */
      }
      if(cpos == (ps->clen+1)) { /* should be a single token, a "\\" on this line */ 
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read a final token at the end of description of draw page %d on line %d of drawfile %s\n", (pp - orig_npage + 1), efp->linenumber, dfile);
	if (strcmp(s, "//") != 0) 
	  esl_fatal("Failed to read a final \"//\" token (read %s) at the end of description of draw page %d on line %d of drawfile %s\n", s, (pp - orig_npage + 1), efp->linenumber, dfile);
	cpos = 0;
      }
      else { 
	/* get C value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read C of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][0] = atof(s);

	/* get M value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read M of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][1] = atof(s);

	/* get Y value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	esl_fatal("Failed to read Y of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][2] = atof(s);

	/* get K value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read K of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][3] = atof(s);

	/* optionally read a residue value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) == eslOK) {
	  if(((int) strlen(s)) != 1) esl_fatal("Read multi-character string (%s) for consensus residue %d on line %d of drawfile %s\n", s, cpos, efp->linenumber, dfile);
	  ps->rrAA[pp][(cpos-1)] = s[0];
	}
	else ps->rrAA[pp][(cpos-1)] = ' ';
      }
    }
  if(pp == (orig_npage - 1)) { /* no new pages were read, this is an error */
    esl_fatal("Failed to read a single page from drawfile %s\n", dfile);
  }

  esl_fileparser_Close(efp);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "drawfile2sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: structural_infocontent_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, colored squares indicating
 *           the structural information content of each base paired consensus column.
 *           Structural information content is the extra information gained from modelling
 *           the pair together (info of vector of bps, size 16) versus separately (sum
 *           of info of the two independent vector of singlets, size 4).
 *           
 * Return:   eslOK on success.
 */
int
structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int ss_idx, int zerores_idx)
{
  int status;
  int p, pp, c, i;
  int cpos, apos, rcpos, lapos, rapos;
  int orig_npage = ps->npage;
  double **obs,   *ent,   *bg;
  double **obs_p, *ent_p, *bg_p;
  double tmp_bg, tmp_bg_p;
  int *ct;
  ESL_DSQ ldsq;
  ESL_DSQ rdsq;
  int *nres;
  int nss = 0;
  int nzerores = 0;
  int nss_masked = 0;
  int nzerores_masked = 0;
  int within_mask;
  if(ps->mask == NULL) nss_masked = nzerores_masked = -1; /* special flag */

  if(msa->ss_cons == NULL) ESL_FAIL(status, errbuf, "--struct requires #=GC SS_cons annotation in the alignment.");
  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  int *c2a_map, *a2c_map, clen;
  map_cpos_to_apos(msa, &c2a_map, &a2c_map, &clen);

  ESL_ALLOC(ent, sizeof(double) * ps->clen);
  ESL_ALLOC(obs, sizeof(double *) * ps->clen);
  ESL_ALLOC(bg,  sizeof(double) * msa->abc->K);
  esl_vec_DSet(bg, msa->abc->K, 1./(msa->abc->K));

  ESL_ALLOC(ent_p, sizeof(double) * ps->clen);
  ESL_ALLOC(obs_p, sizeof(double *) * ps->clen);
  ESL_ALLOC(bg_p,  sizeof(double) * msa->abc->K * msa->abc->K);
  esl_vec_DSet(bg_p, (msa->abc->K * msa->abc->K), 1./(msa->abc->K * msa->abc->K));
  
  ESL_ALLOC(nres, sizeof(int) * ps->clen);
  esl_vec_ISet(nres, ps->clen, 0);

  ESL_ALLOC(ct, sizeof(int) * (msa->alen+1));
  if (esl_wuss2ct(msa->ss_cons, msa->alen, ct) != eslOK) ESL_FAIL(status, errbuf, "structural_infocontent_sspostscript problem getting ct from SS_cons.");

  for(cpos = 0; cpos < ps->clen; cpos++) { 
    ESL_ALLOC(obs[cpos], sizeof(double) * msa->abc->K);
    esl_vec_DSet(obs[cpos], msa->abc->K, 0.);
    ESL_ALLOC(obs_p[cpos], sizeof(double) * (msa->abc->K * msa->abc->K));
    esl_vec_DSet(obs_p[cpos], (msa->abc->K * msa->abc->K), 0.);
  }
  pp = orig_npage;

  /* add color legend */
  float *limits;
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.16;
  limits[2] = 0.33;
  limits[3] = 0.50;
  limits[4] = 0.66;
  limits[5] = 0.83;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits);

  /* get observed residues at each cpos */
  for(i = 0; i < msa->nseq; i++) { 
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { /* seq i is not a gap at cpos */
	  nres[cpos]++;
	  /* only count base paired positions for which both left and right half are not gaps, 
	   * check if we're base paired */
	  if(ct[apos+1] != 0) { 
	    if(ct[apos+1] > (apos+1)) { /* cpos is left half of base pair */
	      /* check if right half is a gap */
	      rapos = ct[apos+1]-1; 
	      if((! esl_abc_CIsGap(msa->abc, msa->aseq[i][rapos])) && 
		  (! esl_abc_CIsGap(msa->abc, msa->rf[rapos]))) { /* seq i is not a gap at right half, and is a cons posn at right half */
		esl_abc_DCount(msa->abc, obs[cpos], esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][apos]), 1.);
		rcpos = a2c_map[rapos];
		assert(rcpos != -1);
		ldsq = esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][apos]);
		rdsq = esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][rapos]);
		PairCount(msa->abc, obs_p[cpos],  ldsq, rdsq, 1.);
		PairCount(msa->abc, obs_p[rcpos], ldsq, rdsq, 1.);
	      }
	    }
	    else { /* cpos is right half of base pair */
	      /* check if left half is a gap */
	      lapos = ct[apos+1]-1; 
	      if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][lapos])) { /* seq i is not a gap at left half */
		esl_abc_DCount(msa->abc, obs[cpos], esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][apos]), 1.);
	      }
	    }
	  }
	}
	cpos++;
      }
    }
  }

  /* determine entropy of each singlet */
  tmp_bg = esl_vec_DEntropy(bg, msa->abc->K);
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    esl_vec_DNorm(obs[cpos], msa->abc->K);
    ent[cpos] = tmp_bg - esl_vec_DEntropy(obs[cpos], msa->abc->K);
  }

  /* determine entropy of each pair */
  tmp_bg_p = esl_vec_DEntropy(bg_p, msa->abc->K * msa->abc->K);
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    apos  = c2a_map[cpos];
    if(ct[apos+1] != 0) { 
      esl_vec_DNorm(obs_p[cpos], msa->abc->K * msa->abc->K);

      rapos = ct[apos+1]-1; 
      rcpos = a2c_map[rapos];

      ent_p[cpos] = tmp_bg_p - esl_vec_DEntropy(obs_p[cpos], msa->abc->K*msa->abc->K);

      /*printf("lpos: %5d  rpos: %5d  entP: %8.3f  entL: %8.3f  entR: %8.3f  ", 
	cpos, rcpos, ent_p[cpos], ent[cpos], ent[rcpos]);*/
      ent_p[cpos] -= (ent[cpos] + ent[rcpos]);
      ent_p[cpos] /= 2.;
      /*printf("final: %8.3f\n", ent_p[cpos]);*/
      if(ent_p[cpos] < (-1. * eslSMALLX1)) { 
	ESL_FAIL(eslEINCONCEIVABLE, errbuf, "pair information < 0.: %f (lpos: %d rpos: %d)\n", ent_p[cpos], cpos, rcpos);
      }
    }
    else ent_p[cpos] = -1.0;
  }

  for(cpos = 0; cpos < ps->clen; cpos++) { 
    if(ent_p[cpos] < (-1. * eslSMALLX1)) { 
      nss++;
      if(ps->mask != NULL && ps->mask[cpos] == '1') nss_masked++; 
    }
    if(nres[cpos] == 0) { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[zerores_idx])) != eslOK) return status;
      ent_p[cpos] = 0.; /* impt to set to 0., so esl_vec_DSum(ent_p... call below to calc total struct info is accurate */
      nzerores++;
      if(ps->mask != NULL && ps->mask[cpos] == '1') nzerores_masked++; 
    }
    else if(ent_p[cpos] < (-1. * eslSMALLX1)) { /* single stranded base, paint grey */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[ss_idx])) != eslOK) return status;
      ent_p[cpos] = 0.; /* impt to set to 0., so esl_vec_DSum(ent_p... call below to calc total struct info is accurate */
    }
    else { 
      within_mask = (ps->mask != NULL && ps->mask[cpos] == '1') ? TRUE : FALSE;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent_p[cpos], ps->sclAA[pp], within_mask)) != eslOK) return status;
    }
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add text to the one cell legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[ss_idx], nss, nss_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "single-stranded", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add text to the second one cell legend */
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[zerores_idx], nzerores, nzerores_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "100% gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 2;

  /* add text to the scheme legend */
  /*sprintf(text,  "structural info content per basepaired posn", esl_vec_DSum(ent_p, ps->clen) * 2.);*/
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "extra information from structure (bits)", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "extra information from structure per basepaired position", errbuf)) != eslOK) return status;

  free(nres);
  free(ent);
  free(ent_p);
  for(cpos = 0; cpos < ps->clen; cpos++) { free(obs[cpos]); free(obs_p[cpos]); }
  free(obs);
  free(obs_p);
  free(bg);
  free(bg_p);
  free(c2a_map);
  free(a2c_map);
  free(ct);
  free(limits);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "infocontent_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: PairCount()
 * Date:     SRE, Tue Aug  1 10:34:20 2000 [St. Louis]
 *
 * Purpose:  Given a possibly degenerate symbol code for left
 *           and right symbols in a pair, increment a symbol
 *           counter array appropriately.
 *           
 * Args:     abc      - pointer to the internal alphabet
 *           counters - vector to count into [0..abc->K^2-1]
 *           syml     - index of left symbol  [0..abc->sym_iupac-1]
 *           symr     - index of right symbol [0..abc->sym_iupac-1]
 *           wt       - weight to use for the count (often 1.0).          
 *
 * Returns:  void
 */
void
PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, float wt)
{
  int status;
  if (syml < abc->K && symr < abc->K) {
    counters[(int) (syml * abc->K + symr)] += wt;
    return;
  }
  else {
    float *left = NULL;
    float *right = NULL;
    ESL_ALLOC(left,  sizeof(float) * abc->K);
    ESL_ALLOC(right, sizeof(float) * abc->K);

    int   l,r;
    
    esl_vec_FSet(left,  abc->K, 0.);
    esl_vec_FSet(right, abc->K, 0.);
    esl_abc_FCount(abc, left,  syml, wt);
    esl_abc_FCount(abc, right, symr, wt);

    for (l = 0; l < abc->K; l++)
      for (r = 0; r < abc->K; r++)
	counters[l*abc->K +r] += left[l] * right[r];
    free(left);
    free(right);
  }
  return;

 ERROR:
  esl_fatal("Memory error");
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call ssu-draw
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

/* Function: get_date
 * Date:     EPN, Fri Jan 25 13:59:22 2008
 *
 * Purpose:  Return a string that gives the current date.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_date(char *errbuf, char **ret_date)
{
  int    status;
  time_t date = time(NULL);
  char  *sdate = NULL;

  if((status = esl_strdup(ctime(&date), -1, &sdate)) != eslOK) goto ERROR;
  esl_strchop(sdate, -1); /* doesn't return anything but eslOK */

  *ret_date = sdate;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_date() error status: %d, probably out of memory.", status);
  return status; 
}

/* Function: set_scheme_values()
 * 
 * Purpose:  Set color values from a predefined scheme given min, max, 
 *           value and number of colors.
 */ 
static int 
set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl, int within_mask)
{
  float min, max;
  int ci, bi;
  min = scl->limits[0];
  max = scl->limits[scl->nbins];
  if((min-val) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f < min: %.4f\n", val, min); }
  if((val-max) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f > max: %.4f\n", val, max); }

  bi = 0;
  while((val > scl->limits[bi+1]) && (bi <= (scl->nbins-1))) { bi++; }    
  /*printf("%.3f %d (%.3f)\n", val, bi, scl->limits[bi+1]);*/
  scl->counts[bi]++;
  if(within_mask) scl->counts_masked[bi]++;
  for(ci = 0; ci < ncolvals; ci++) { 
    vec[ci] = scheme[bi][ci]; 
  }
  return eslOK;
}

/* Function: set_onecell_values()
 * 
 * Purpose:  Set color values as a  predefined single
 *           color.
 */ 
static int 
set_onecell_values(char *errbuf, float *vec, int ncolvals, float *onecolor)
{
  int ci;
  for(ci = 0; ci < ncolvals; ci++) { vec[ci] = onecolor[ci]; }
  return eslOK;
}
  


/* Function: draw_masked_block()
 * 
 * Purpose:  Given coords, color, and mask style options draw a masked block.
 */ 
static int 
draw_masked_block(FILE *fp, float x, float y, float *colvec, int do_circle_mask, int do_square_mask, int do_x_mask, int do_border, float boxsize)
{
  if (do_circle_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, " %.2f %.2f %.1f 0 360 arc closepath\n", x + (boxsize/2.), y + (boxsize/2.), boxsize * (3./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  stroke\n"); 
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, " %.2f %.2f %.1f 0 360 arc closepath\n", x + (boxsize/2.), y + (boxsize/2.), boxsize * (3./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
    }
  }
  else if(do_square_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x +1., y +1.);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize*0.75, boxsize*0.75, boxsize*0.75);
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  stroke\n");
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x + 1.5, y + 1.5);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize*(5./8.), boxsize*(5./8.), boxsize*(5./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
    }
  }
  else if (do_x_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize, boxsize, boxsize);
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
      
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 0.);
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  %.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
      fprintf(fp, "  %.2f %.2f moveto", x + boxsize, y);
      fprintf(fp, "  -%.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  %.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x + boxsize, y);
      fprintf(fp, "  -%.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
    }
  }
  return eslOK;
}

/* Function: validate_justread_sspostscript()
 * 
 * Purpose:  Validate a sspostscript just created by parsing
 *           a template file. Nothing fancy here, just make
 *           sure all we've written everything we expect.
 */ 
static int 
validate_justread_sspostscript(SSPostscript_t *ps, char *errbuf)
{
  if(ps->modelname == NULL) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read modelname from template file.");
  if(ps->nbp == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines bpconnects' section from template file.");
  if(ps->clen == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'text residues' section from template file.");

  /* Stuff we don't currently require, but we may want to eventually */
  /*if(ps->nhundreds == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'text hundreds' section from template file.");*/
  /*if(ps->nticks == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines ticks' section from template file.");*/
  /*if(ps->nbp == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines bpconnects' section from template file.");*/

  return eslOK;
}


/* Function: validate_and_update_sspostscript_given_msa()
 * 
 * Purpose:  Validate that a sspostscript works with a MSA.
 */ 
static int 
validate_and_update_sspostscript_given_msa(const ESL_GETOPTS *go, SSPostscript_t *ps, ESL_MSA *msa, char *errbuf)
{
  int status;
  int *msa_ct;
  int msa_nbp = 0;
  int *tmp_ct;
  int msa_clen;
  int apos, cpos;
  int i;

  ps->msa = msa;

  /* get the CT array for this msa */
  ESL_ALLOC(tmp_ct, sizeof(int) * (msa->alen+1));
  if (esl_wuss2ct(msa->ss_cons, msa->alen, tmp_ct) != eslOK) ESL_FAIL(status, errbuf, "Problem getting ct from SS_cons, does first alignment of MSA file have SS_cons annotation?");
  
  msa_clen = 0;
  for(apos = 0; apos < msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) msa_clen++;
  }
  /* convert tmp_ct which is in alignment coords [1..alen] to consensus coords [0..clen-1]*/
  ESL_ALLOC(msa_ct, sizeof(int) * (msa_clen));
  cpos = 0;
  for(apos = 0; apos < msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
      if((tmp_ct[(apos+1)] > (apos+1)) && (! esl_abc_CIsGap(msa->abc, msa->rf[tmp_ct[apos+1]-1]))) { /* a consensus position paired to another consensus posn */
	msa_nbp++;
      }
      msa_ct[cpos++] = tmp_ct[(apos+1)];
    }
  }
  free(tmp_ct);

  if(ps->msa_ct != NULL) { free(ps->msa_ct); ps->msa_ct = NULL;}
  ps->msa_ct = msa_ct;
  ps->msa_nbp = msa_nbp;

  if(ps->clen != msa_clen) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected consensus length of %d in MSA, but read %d\n", ps->clen, msa_clen);
  if(ps->nbp != 0 && ps->nbp != msa_nbp) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected %d basepairs in MSA's SS_cons, but read %d\n", ps->nbp, msa_nbp);

  /* allocate the uaseqlenA data structure to store unaligned seq lengths, only compute them if --indi though */
  if(ps->uaseqlenA != NULL) { free(ps->uaseqlenA); ps->uaseqlenA = NULL; }
  ESL_ALLOC(ps->uaseqlenA, sizeof(int) * msa->nseq);
  esl_vec_ISet(ps->uaseqlenA, msa->nseq, 0);
  if(esl_opt_GetBoolean(go, "--indi")) { /* individual mode, get each sequence's length */
    for(i = 0; i < msa->nseq; i++) { 
      for(apos = 0; apos < msa->alen; apos++) {
	if(! esl_abc_CIsGap(ps->msa->abc, ps->msa->aseq[i][apos])) 
	  ps->uaseqlenA[i]++;
      }
    }
  }
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "validate_and_update_sspostscript_given_msa(), error status %d, probably out of memory.\n", status);
  return status; 
  }


/* Function: draw_header_and_footer()
 * 
 * Purpose:  Draw header for a page
 */ 
static int 
draw_header_and_footer(FILE *fp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int page, int pageidx2print)
{
  int status;
  int i, split_idx;
  float x, y;
  float header_fontsize;
  int model_width, desc_width, desc_column_width;
  char *model_dashes, *desc_dashes, *desc2print;
  char *desc_string = NULL;
  float xmodel;
  char *model2print = NULL;

  header_fontsize = HEADER_FONTSIZE_UNSCALED / ps->scale; 

  fprintf(fp, "%% begin ignore\n");
  fprintf(fp, "/%s findfont %.2f scalefont setfont\n", DEFAULT_FONT, header_fontsize);
  fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */

  if(! (esl_opt_GetBoolean(go, "--no-head"))) { 
    model_width = ESL_MAX(strlen("model"), (int) strlen(ps->modelname));
    if(model_width > HEADER_MODELNAME_MAXCHARS) { 
      ESL_ALLOC(model2print, sizeof(char) * (HEADER_MODELNAME_MAXCHARS+1));
      for(i = 0; i < (HEADER_MODELNAME_MAXCHARS-3); i++) model2print[i] = ps->modelname[i];
      model2print[HEADER_MODELNAME_MAXCHARS-3] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS-2] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS-1] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS] = '\0';
    }
    else { 
      if((status = esl_strdup(ps->modelname, -1, &(model2print))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying modelname");
    }
    
    model_width = ESL_MIN(model_width, HEADER_MODELNAME_MAXCHARS);
    ESL_ALLOC(model_dashes, sizeof(char) * (model_width+1));
    for(i = 0; i < model_width; i++) model_dashes[i] = '-'; 
    model_dashes[model_width] = '\0';
    
    if(ps->modeA[page] == ALIMODE || ps->modeA[page] == SIMPLEMASKMODE) { 
      if((status = esl_strdup("description", -1, &(desc_string))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying description");
    }
    else if(ps->modeA[page] == INDIMODE) { 
      if((status = esl_strdup("sequence name", -1, &(desc_string))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying description");
    }
    
    xmodel = ps->headerx_desc - (ps->headerx_charsize * (model_width  + 6 + 6 + 8 + 2)); /*6,6,8 for #res,#bps,#seq|seqlen) plus 2 spaces each, first 2 for 2 spaces before desc*/
    x = xmodel;
    y = ps->headery;
    
    fprintf(fp, "(%-*s  %4s  %4s) %.2f %.2f moveto show\n", model_width, "model", "#res", "#bps", x, y);
    y -= header_fontsize * 0.75;
    fprintf(fp, "(%-*s  %4s  %4s) %.2f %.2f moveto show\n", model_width, model_dashes, "----", "----", x, y);
    y -= header_fontsize * 0.75;
    fprintf(fp, "(%-*s  %4d  %4d) %.2f %.2f moveto show", model_width, model2print, ps->clen, ps->msa_nbp, x, y);
    free(model_dashes);
    x += (model_width + 6 + 6 +2) * ps->headerx_charsize;
    
    if(ps->modeA[page] == ALIMODE) { 
      y+= header_fontsize * 1.5;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "#seqs", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "------", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6d) %.2f %.2f moveto show", ps->msa->nseq, x, y);
    }
    else if(ps->modeA[page] == INDIMODE && (ps->seqidxA[page] != -1)) { /* ps->seqidxA[page] == -1 if we're printing the consensus sequence */
      y+= header_fontsize * 1.5;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "seqlen", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "------", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6d) %.2f %.2f moveto show", ps->uaseqlenA[ps->seqidxA[page]], x, y);
    }
    
    if(ps->descA[page] != NULL) { 
      x =  ps->headerx_desc;
      y += 2. * header_fontsize * 0.75;
      desc_width = ESL_MAX((int) strlen(desc_string), (int) strlen(ps->descA[page]));
      if(desc_width > ps->desc_max_chars) { 
	/* split into two lines, the add_page_desc_to_sspostscript() function added a '\n' where the split should be */
	i = 0; 
	while(ps->descA[page][i] != '\n') { 
	  i++; 
	  if(i >= desc_width) ESL_FAIL(eslEINVAL, errbuf, "drawing header, failed to find split point from add_page_desc_to_() in two-line description (%s)", ps->descA[page]);
	}
	split_idx = i;
	desc_column_width = split_idx;
      }
      else desc_column_width = desc_width;
      
      ESL_ALLOC(desc_dashes, sizeof(char) * (desc_column_width+1));
      for(i = 0; i < desc_column_width; i++) desc_dashes[i] = '-'; 
      desc_dashes[desc_column_width] = '\0';
      
      fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc_string, x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc_dashes, x, y);
      y -= header_fontsize * 0.75;
      free(desc_dashes);
      
      if(desc_width > ps->desc_max_chars) {
	ESL_ALLOC(desc2print, sizeof(char) * (split_idx+1));
	for(i = 0; i < split_idx; i++) desc2print[i] = ps->descA[page][i];
	desc2print[split_idx] = '\0';
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc2print, x, y);
	free(desc2print);
	
	x = ps->headerx_desc;
	y-= ps->headery_charsize * 1;
	ESL_ALLOC(desc2print, sizeof(char) * ((desc_width - split_idx) -1+1));
	for(i = split_idx+1; i < desc_width; i++) desc2print[(i-(split_idx+1))] = ps->descA[page][i];
	desc2print[(desc_width-(split_idx+1))] = '\0';
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc2print, x, y);
	free(desc2print);
      }
      else { 
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_width, ps->descA[page], x, y);
      }
    }
    /* masked row of header goes here if desired */
  }

  /* draw footer */
  float footer_fontsize, footerx_charsize;
  footer_fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;
  footerx_charsize = ps->legx_charsize;
  
  fprintf(fp, "/%s findfont %.2f scalefont setfont\n", DEFAULT_FONT, footer_fontsize);
  if(! (esl_opt_GetBoolean(go, "--no-foot"))) { 
    /* draw alignment file name in lower left hand corner */
    if(ps->mask != NULL) { 
      if(esl_opt_GetString(go, "--mask-diff") != NULL) { 
	fprintf(fp, "(alifile: %s; mask 1 file: %s; mask 2 file: %s;) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--mask"), esl_opt_GetString(go, "--mask-diff"), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
      }
      else { 
	fprintf(fp, "(alifile: %s; mask file: %s;) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--mask"), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
      }
    }
    else { 
      fprintf(fp, "(alifile: %s) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
    }

    /* put page number */
    /* determine ndigits */
    int tmp, ndigits;
    tmp = pageidx2print;
    ndigits = 1;
    while(tmp > 10) { tmp /= 10; ndigits++; }
    x = ps->pagex_max - (PAGE_SIDEBUF) - (footerx_charsize * (5 + ndigits)); 
    fprintf(fp, "(page %d) %.2f %.2f moveto show\n", pageidx2print, x, PAGE_BOTBUF);
    
  }    
  fprintf(fp, "(structure diagram derived from CRW database: http://www.rna.ccbb.utexas.edu/) %.2f %.2f moveto show\n", PAGE_SIDEBUF , PAGE_BOTBUF);
  fprintf(fp, "%% end ignore\n");

  if(model2print != NULL) free(model2print);
  return eslOK;

 ERROR: ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), memory error.");
}


