/* Draw SSU secondary structure diagrams given a Gutell SS template 
 * and an SSU alignment. 
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

#define RAINBOWRHSCHEME 0
#define RAINBOWRLSCHEME 1
#define NRAINBOWRHSCHEME 11
#define NRAINBOWRLSCHEME 11

#define RBSIXRHSCHEME 2
#define RBSIXRLSCHEME 3
#define NRBSIXRHSCHEME 6
#define NRBSIXRLSCHEME 6

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
#define LEG_ONED_NBOXES  11
#define LEG_ONED_BOXSIZE 24.
#define LEG_MINTEXTSIZE 10
/*#define TITLE_FONTSIZE 18*/
#define TITLE_FONTSIZE 12

/* Structure: color_legend
 * Incept:    EPN, Tue Aug 26 06:50:01 2008
 *
 * Parameters describing a legend for colors for an
 * SSPostscript_t data structure.
 */
typedef struct color_legend_s {
  int    which_color[NCMYK];/* [CMYK] TRUE or FALSE, if that color is to be used on the color legend */  
  float  min[NCMYK];        /* [CMYK] min value to start at for each color (irrelevant if which_color[] is FALSE */  
  float  max[NCMYK];        /* [CMYK] max value to end at for each color (irrelevant if which_color[] is FALSE */  
  float  x;                 /* x coordinate to start legend at */
  float  y;                 /* y coordinate to start legend at */
  char   *text[NCMYK];      /* text for legend, one string for each color */
  float  boxsize;           /* size of box for each residue */
  int    nboxes;            /* number of boxes */ 
} ColorLegend_t;

/* Structure: scheme_color_legend
 * Incept:    EPN, Thu Jun 25 20:20:38 2009
 *
 * Parameters describing a one-dimensional legend of colors
 * from a preset scheme for use in a SSPostscript_t data structure.
 */
typedef struct scheme_color_legend_s {
  int    scheme;            /* preset color scheme index */
  int    nbins;             /* number of colors (bins) in this scheme */
  //  int    min;               /* min value to start at for coloring */
  //  int    max;               /* max value to end at for coloring */
  //float  scale;             /* scale for values, if 100., a value of 31 really equals 31/100.=0.31 */
  float  x;                 /* x coordinate to start legend at */
  float  y;                 /* y coordinate to start legend at */
  char   *text;             /* text for legend, a single string */
  float  boxsize;           /* size of box for each residue */
  float *limits;            /* [nbins+1] limits for each bin, limits[0] is min value we would expect to see, limits[nbins] is max */
} SchemeColorLegend_t;

/* Structure: one_cell_color_legend
 * Incept:    EPN, Tue Sep 30 13:06:15 2008
 *
 * Parameters describing a single colored cell legend for a
 * SSPostscript_t data structure.
 */
typedef struct one_cell_color_legend_s {
  float  col[NCMYK];        /* [CMYK] color value for the cell */
  float  x;                 /* x coordinate to start legend at */
  float  y;                 /* y coordinate to start legend at */
  char   *text;             /* text for legend */
  float  boxsize;           /* size of box for each residue */
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
  char  **regurgAA;     /* [0..nregurg-1][] lines from the original Gutell file to regurgitate, these are unchanged. */
  int     nregurg;      /* number of lines (char *'s) in the regurgAA 2D array */
  int     npage;        /* number of pages in eventual postscript */
  int     clen;         /* the number of residues in the Gutell template file */
  int     title_begin;  /* line on which the title information that we'll rewrite begins */
  int     title_nlines; /* number of lines of title information */
  float   titlex;       /* x coordinate (upper left corner) of title area */
  float   titley;       /* y coordinate (upper left corner) of title area */
  float   legx;         /* x coordinate (upper left corner) of legend area */
  float   legy;         /* y coordinate (upper left corner) of legend area */
  float   cur_legx;     /* current x coordinate (upper left corner) of legend area */
  float   cur_legy;     /* current y coordinate (upper left corner) of legend area */
  float  *rxA;          /* [0..clen-1] x coordinate for each residue in the eventual postscript */
  float  *ryA;          /* [0..clen-1] y coordinate for each residue in the eventual postscript */
  char  **rrAA;         /* [0..npage-1][0..clen-1] residue character in the eventual postscript */
  float ***rcolAAA;     /* [0..npage-1][0..clen-1][0..3] color for block on page p, position c, CMYK in the eventual postscript */
  ColorLegend_t ***clAAA;/* [0..npage-1][0..l..nclA[p]  ptr to color legend l for page p */
  int     *nclA;        /* [0..npage-1] number of color legends for each page */
  OneCellColorLegend_t ***occlAAA;/* [0..npage-1][0..l..nocclA[p]  ptr to one cell color legend l for page p */
  int     *nocclA;      /* [0..npage-1] number of one cell color legends for each page */
  SchemeColorLegend_t  **sclAA;/* [0..npage-1]  ptr to scheme color legend l for page p, NULL if none */
  char   **maskAA;      /* [0..npage-1][0..clen-1] mask, columns which are '0' get drawn differently */
} SSPostscript_t;


static SSPostscript_t *create_sspostscript();
static ColorLegend_t *create_colorlegend(SSPostscript_t *ps, int *for_which_colors, float *min, float *max, float boxsize, int nboxes, char **text);
static ColorLegend_t *create_one_dim_colorlegend(SSPostscript_t *ps, int color_idx, float min, float max, float boxsize, int nboxes, char *text);
static OneCellColorLegend_t *create_one_cell_colorlegend(SSPostscript_t *ps, float *cmykA, float boxsize, char *text);
static SchemeColorLegend_t *create_scheme_colorlegend(SSPostscript_t *ps, int scheme, int ncols, float boxsize, char *text, float *limits);
static int  add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text);
static int  print_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps);
static int  print_colorlegend(FILE *fp, ColorLegend_t *cl);
static int  print_onecellcolorlegend(FILE *fp, OneCellColorLegend_t *occl);
static int  print_scheme_colorlegend(FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme);
static void free_sspostscript(SSPostscript_t *ps);
static int  addpages_sspostscript(SSPostscript_t *ps, int ntoadd);
static int  map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int **ret_a2c_map, int *ret_clen);
static int  read_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps);
static int  individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  rf_seq_sspostscript (const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins);
static int  delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins);
static int  insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins);
static int  compare_ints(const void *el1, const void *el2);
static int  posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen);
static int  drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps);
static void PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, float wt);
static int  structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_bins);
static int  phylosignal_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask);
static int  colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask);
static int  diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  get_date(char *errbuf, char **ret_date);
static int  set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl);
static int  add_mask_to_ss_postscript(SSPostscript_t *ps, int page, char *mask);

static char banner[] = "draw Gutell based postscript SSU secondary structure diagrams.";
static char usage[]  = "[options] <msafile> <Gutell SS postscript template> <output postscript file name>\n\
The <msafile> must be in Stockholm format.";

#define MASKOPTS "--mask,--mask-col" /* exclusive choice for masking */
#define MASKTYPEOPTS "-d,-c,-x" /* exclusive choice for mask types */

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-q",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "DO NOT create SS info content diagram (on by default)", 0 },
  { "-s",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create SS diagram for each sequence in the alignment",    0 },
  { "-i",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create insert SS diagram",    1 },
  { "-d",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create insert SS diagram",    1 },
  { "-c",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create insert SS diagram",    1 },
  { "-x",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create insert SS diagram",    1 },
  { "--rf",     eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create SS diagram for RF sequence",    1 },
  { "--struct", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create structural info content SS diagram",    1 },
  { "--p-avg",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create average posterior probability SS diagram",1 },
  { "--p-indi", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create posterior probability diagram for each sequence",1 },
  { "--p-min",  eslARG_REAL,  "0.90",NULL, "0.09999<x<=1.",NULL,NULL,NULL,       "set minimum posterior probability to color to <x>",1 },
  { "--phy",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create SS diagram displaying phylogenetic signal per position", 1},
  { "--dall",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create delete diagram w/all deletions (incl. terminal deletes)",    1 },
  { "--dint",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create delete diagram w/only internal (non-terminal) deletions", 1 },
  { "--mask",    eslARG_INFILE, NULL, NULL, NULL, NULL,NULL,MASKOPTS,        "for all diagrams, mark masked columns from mask in <f>", 1 },
  { "--mask-col",eslARG_INFILE,NULL, NULL, NULL, NULL,NULL, MASKOPTS,        "create black/pink colored SS diagram denoting masked columns", 1 },
  { "--mask-diff",eslARG_INFILE,NULL, NULL, NULL, NULL,"--mask-col","--mask","with --mask-col <f1>, compare mask in <f1> to mask in <f>", 1},
  { "--dfile",  eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, NULL,            "read 'draw' file specifying >=1 SS diagram drawings", 1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  ESL_ALPHABET *abc     = NULL;	/* biological alphabet             */
  char         *alifile = NULL;	/* alignment file name             */
  char         *outfile = NULL;	/* output ps file name             */
  char         *templatefile = NULL; /* Gutell template file       */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *msa     = NULL;	/* one multiple sequence alignment */
  int           status;		/* easel return code               */
  int           clen;           /* non-gap RF (consensus) length of each alignment */
  int           nali = 0;	/* number of alignments read       */
  int           apos;
  char          errbuf[eslERRBUFSIZE];
  FILE         *ofp;		/* output file for postscript */
  
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
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 3) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  alifile = esl_opt_GetArg(go, 1);
  templatefile = esl_opt_GetArg(go, 2);
  outfile = esl_opt_GetArg(go, 3);

  char *command;
  char *date;
  if((status = get_command(go, errbuf, &command)) != eslOK) esl_fatal(errbuf);
  if((status = get_date(errbuf, &date))           != eslOK) esl_fatal(errbuf);

  SSPostscript_t *ps;
  if((status = read_template_file(templatefile, go, errbuf, &ps) != eslOK)) esl_fatal(errbuf);

  /***********************************************
   * Open the MSA file; determine alphabet; set for digital input
   ***********************************************/

  fmt = eslMSAFILE_STOCKHOLM;
  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed with error %d\n", status);

  /* open PS output file for writing */
  if ((ofp = fopen(outfile, "w")) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open output postscript file %s\n", esl_opt_GetArg(go, 3));

  /* Assert RNA, it's the ribosome */
  abc = esl_alphabet_Create(eslRNA);
  afp->abc = abc;

  char *mask = NULL;
  int masklen;
  char *mask2 = NULL;
  int masklen2;
  if(! esl_opt_IsDefault(go, "--mask")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask"), errbuf, &mask, &masklen)) != eslOK) esl_fatal(errbuf);
  }
  if(! esl_opt_IsDefault(go, "--mask-col")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-col"), errbuf, &mask, &masklen)) != eslOK) esl_fatal(errbuf);
    if(! esl_opt_IsDefault(go, "--mask-diff")) { 
      if((status = read_mask_file(esl_opt_GetString(go, "--mask-diff"), errbuf, &mask2, &masklen2)) != eslOK) esl_fatal(errbuf);
      if(masklen != masklen2) esl_fatal("Mask in %f length (%d) differs from mask in %f (%d)!", esl_opt_GetString(go, "--mask-col"), masklen, esl_opt_GetString(go, "--mask-diff"), masklen2);
    }
  }

  /***********************************/
  /* allocate and fill predefined color schemes, these are hardcoded */
  int     *hc_nbins;
  float ***hc_scheme;
  int z;
  ESL_ALLOC(hc_scheme, sizeof(float **) * 4);
  ESL_ALLOC(hc_scheme[0], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[0][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[1], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[1][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[2], sizeof(float *) * 6); 
  for(z = 0; z < 6; z++) { ESL_ALLOC(hc_scheme[2][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[3], sizeof(float *) * 6); 
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
  hc_scheme[3][0][0] = 0.00; hc_scheme[3][0][1] = 0.94; hc_scheme[3][0][2] = 1.00; hc_scheme[3][1][3] = 0.00; /*red*/
  /***************************************************************/
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      msa->abc = abc;
      if(msa->rf == NULL) esl_fatal("MSA number: %d in %s does not have RF annotation.", nali, alifile);
      clen = 0;
      for(apos = 0; apos < msa->alen; apos++) if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) clen++;
      if(ps->clen == 0)    esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != template file consensus length of %d. Did you add the 'residue_start' line?", nali, clen, ps->clen);
      if(clen != ps->clen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != template file consensus length of %d.", nali, clen, ps->clen);
      if(mask != NULL && ps->clen != masklen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != lane mask length of %d from mask file %s.", nali, clen, masklen, esl_opt_GetString(go, "--mask"));

      if(! esl_opt_GetBoolean(go, "-q")) { 
	if((status = infocontent_sspostscript(go, errbuf, ps, msa, mask, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--struct")) { 
	if((status = structural_infocontent_sspostscript(go, errbuf, ps, msa, NULL, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "-i")) { /* make a new postscript page marking insertions */
	if((status = insert_sspostscript(go, errbuf, ps, msa, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dall")) { /* make a new postscript page marking all deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, TRUE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dint")) { /* make a new postscript page marking internal deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, FALSE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--rf")) { /* make a new postscript page for the RF sequence in the alignment */
	if((status = rf_seq_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
      }
      int do_post = (esl_opt_GetBoolean(go, "--p-avg")) ? TRUE : FALSE;
      if(do_post) { 
	if((status = posteriors_sspostscript(go, errbuf, ps, msa, NULL, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME])) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "-s")) { /* make a new postscript page for each sequence in the alignment */
	if((status = individual_seqs_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
      }
      if(! esl_opt_IsDefault(go, "--dfile")) { 
	if((status = drawfile2sspostscript(go, errbuf, ps)) != eslOK) esl_fatal(errbuf);
      }
      if(! esl_opt_IsDefault(go, "--mask-col")) { 
	if(! esl_opt_IsDefault(go, "--mask-diff")) { 
	  if((status = diffmask_sspostscript(go, errbuf, ps, msa, mask, mask2)) != eslOK) esl_fatal(errbuf);
	}
	else {
	  if(ps->clen != masklen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != lane mask length of %d.", nali, clen, masklen);
	  if((status = colormask_sspostscript(go, errbuf, ps, msa, mask)) != eslOK) esl_fatal(errbuf);
	}
      }

      if((status = print_sspostscript(ofp, go, errbuf, command, date, hc_scheme, ps)) != eslOK) esl_fatal(errbuf);
      fclose(ofp);
      esl_msa_Destroy(msa);
    }

  /* If an msa read failed, we drop out to here with an informative status code. 
   */
  if (status == eslEFORMAT) 
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)
    esl_fatal("No alignments found in file %s\n", alifile);

  /* Cleanup, normal return
   */
  if(mask != NULL) free(mask);
  free_sspostscript(ps);
  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
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

  ps->regurgAA = NULL;
  ps->nregurg  = 0;
  ps->npage    = 0;
  ps->clen     = 100;
  ps->title_nlines = 0;
  ESL_ALLOC(ps->rxA,  sizeof(float) * ps->clen);
  ESL_ALLOC(ps->ryA,  sizeof(float) * ps->clen);
  ps->rrAA        = NULL;
  ps->rcolAAA     = NULL;
  ps->clAAA       = NULL;
  ps->occlAAA     = NULL;
  ps->sclAA       = NULL;
  ps->maskAA      = NULL;
  ps->legx = ps->legy = ps->cur_legx = ps->cur_legy = 0.;

  return ps;

 ERROR: esl_fatal("create_sspostscript(): memory allocation error.");
  return NULL; /* NEVERREACHED */
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

  if(ps->regurgAA != NULL) {
    for(i = 0; i < ps->nregurg; i++) { 
      if(ps->regurgAA[i] != NULL) {
	free(ps->regurgAA[i]);
      }
    }
    free(ps->regurgAA);
  }

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

  if(ps->clAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->clAAA[p] != NULL) { 
	for(l = 0; l < ps->nclA[p]; l++) { 
	  free(ps->clAAA[p][l]); /* all statically allocated memory */
	}
      }
    }
  }

  if(ps->occlAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->occlAAA[p] != NULL) { 
	for(l = 0; l < ps->nocclA[p]; l++) { 
	  free(ps->occlAAA[p][l]); /* all statically allocated memory */
	}
      }
    }
  }

  if(ps->sclAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->sclAA[p] != NULL) { 
	if(ps->sclAA[p]->limits != NULL) free(ps->sclAA[p]->limits);
	free(ps->sclAA[p]); /* statically allocated memory */
      }
    }
  }

  if(ps->sclAA != NULL) free(ps->sclAA);
  if(ps->nclA != NULL) free(ps->nclA);
  if(ps->nocclA != NULL) free(ps->nocclA);

  free(ps);
  return;
}

/* Function: create_colorlegend()
 * 
 * Purpose:  Create and initialize a color legend data structure.
 * Return:   cl
 */
ColorLegend_t *
create_colorlegend(SSPostscript_t *ps, int *which_color, float *min, float *max, float boxsize, int nboxes, char **text)
{
  int status;
  ColorLegend_t *cl;
  int c;

  ESL_ALLOC(cl, sizeof(ColorLegend_t));

  /* initialize */
  esl_vec_ISet(cl->which_color, NCMYK, FALSE);
  esl_vec_FSet(cl->min, NCMYK, 0.);
  esl_vec_FSet(cl->max, NCMYK, 0.);
  cl->x = ps->cur_legx;
  cl->y = ps->cur_legy;
  for(c = 0; c < NCMYK; c++) cl->text[c] = NULL;

  /* set caller specified values */
  esl_vec_ICopy(which_color, NCMYK, cl->which_color);
  esl_vec_FCopy(min, NCMYK, cl->min);
  esl_vec_FCopy(max, NCMYK, cl->max);
  cl->boxsize = boxsize;
  cl->nboxes  = nboxes;
  for(c = 0; c < NCMYK; c++) cl->text[c] = text[c];

  /* update postscripts legx, legy */
  ps->cur_legx += 0.;
  ps->cur_legy -= 3 * boxsize;

  return cl;

 ERROR: esl_fatal("create_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: create_one_dim_colorlegend()
 * 
 * Purpose:  Create and initialize a color legend data structure.
 * Return:   cl
 */
ColorLegend_t *
create_one_dim_colorlegend(SSPostscript_t *ps, int color_idx, float min, float max, float boxsize, int nboxes, char *text)
{
  int status;
  int which_color[NCMYK];
  float minA[NCMYK];
  float maxA[NCMYK];
  char **textAA;
  int c;

  if(color_idx >= NCMYK) esl_fatal("create_one_dim_colorlegend(), color_idx %d invalid (must be < %d)", color_idx, NCMYK);
  if(min > max)          esl_fatal("create_one_dim_colorlegend(), min (%f) > max (%f)", min, max);

  esl_vec_ISet(which_color, NCMYK, FALSE);
  esl_vec_FSet(minA, NCMYK, 0.);
  esl_vec_FSet(maxA, NCMYK, 0.);
  ESL_ALLOC(textAA, sizeof(char *) * NCMYK);
  for(c = 0; c < NCMYK; c++) textAA[c] = NULL;

  which_color[color_idx] = TRUE;
  minA[color_idx] = min;
  maxA[color_idx] = max;
  if((status = esl_strdup(text, -1, &(textAA[color_idx]))) != eslOK) esl_fatal("create_one_dim_colorlegend(), error copying text");

  return create_colorlegend(ps, which_color, minA, maxA, boxsize, nboxes, textAA);
									      
 ERROR:
  esl_fatal("create_one_dim_colorlegend(), memory error.");
  return NULL; /* NEVERREACHED */
}  

/* Function: create_one_cell_colorlegend()
 * 
 * Purpose:  Create and initialize a one cell color legend data structure.
 * Return:   occl
 */
OneCellColorLegend_t *
create_one_cell_colorlegend(SSPostscript_t *ps, float *col, float boxsize, char *text)
{
  int status;
  OneCellColorLegend_t *occl;

  ESL_ALLOC(occl, sizeof(OneCellColorLegend_t));

  /* initialize */
  esl_vec_FSet(occl->col, NCMYK, 0.);
  occl->x = ps->cur_legx;
  occl->y = ps->cur_legy;
  occl->text = NULL;
  
  /* set caller specified values */
  esl_vec_FCopy(col, NCMYK, occl->col);
  occl->boxsize = boxsize;
  if((status = esl_strdup(text, -1, &(occl->text))) != eslOK) esl_fatal("create_one_cell_colorlegend(), error copying text");
  /* update postscripts legx, legy */
  ps->cur_legx += 0.;
  ps->cur_legy -= 2 * boxsize;

  return occl;

 ERROR: esl_fatal("create_one_cell_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: create_scheme_colorlegend()
 * 
 * Purpose:  Create and initialize a scheme color legend data structure.
 * Return:   scl
 */
SchemeColorLegend_t *
create_scheme_colorlegend(SSPostscript_t *ps, int scheme, int nbins, float boxsize, char *text, float *limits)
{
  int status;
  SchemeColorLegend_t *scl;
  int i;
  ESL_ALLOC(scl, sizeof(SchemeColorLegend_t));

  /* initialize */
  scl->x = ps->cur_legx;
  scl->y = ps->cur_legy;
  scl->text = NULL;
  
  /* set caller specified values */
  scl->scheme = scheme;
  scl->nbins = nbins;
  ESL_ALLOC(scl->limits, sizeof(float) * (nbins+1));
  for(i = 0; i <= nbins; i++) { scl->limits[i] = limits[i]; printf("LIMIT[%d]: %.3f\n", i, scl->limits[i]); }
  scl->boxsize = boxsize;
  if((status = esl_strdup(text, -1, &(scl->text))) != eslOK) esl_fatal("create_scheme_colorlegend(), error copying text");
  /* update postscripts legx, legy */
  ps->cur_legx += 0.;
  ps->cur_legy -= 2 * boxsize;

  return scl;

 ERROR: esl_fatal("create_scheme_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: add_text_to_scheme_colorlegend()
 * 
 * Purpose:  Add text to an existing scheme color legend data structure.
 * Return:   scl
 */
int
add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text)
{
  int status;
  if(scl->text != NULL) esl_fatal("add_text_to_scheme_colorlegend(), text already exists!\n"); 
  if(text == NULL) esl_fatal("add_text_to_scheme_colorlegend(), passed in text is NULL!\n"); 
  if((status = esl_strdup(text, -1, &(scl->text))) != eslOK) esl_fatal("add_text_to_scheme_colorlegend(), error copying text");
  return eslOK;
}

/* Function: add_mask_to_sspostscript
 * 
 * Purpose:  Add a mask to a sspostscript object.
 */
int
add_mask_to_ss_postscript(SSPostscript_t *ps, int page, char *mask)
{
  int status;
  if(ps->maskAA[page] != NULL) { esl_fatal("add_mask_to_ss_postscript(), mask for page: %d is non-null!\n", page); }
  if(mask == NULL) esl_fatal("add_mask_to_ss_postscript(), passed in mask is NULL!\n"); 
  if((status = esl_strdup(mask, -1, &(ps->maskAA[page]))) != eslOK) esl_fatal("add_mask_to_ss_postscript(), error copying mask");
  return eslOK;
}

/* Function: print_colorlegend()
 * 
 * Purpose:  Print a color legend to an open file.
 * Return:   eslOK
 */
int 
print_colorlegend(FILE *fp, ColorLegend_t *cl)
{
  int ndims = 0;
  float cur[NCMYK];     /* CMYK value for current color we're printing */
  float colstep[NCMYK]; /* step for each color */
  float x, y;
  int c, b, cp;
  int textlen;
  float textsize;

  /* verify we have a valid color legend object */
  for(c = 0; c < NCMYK; c++) { 
    if(cl->which_color[c] == TRUE) { 
      if(cl->max[c] < cl->min[c]) esl_fatal("print_colorlegend(): colorlegend object is invalid, max[%d]: %f > min[%d]: %f\n", c, cl->max[c], cl->min[c]);
      ndims++;
      colstep[c] = (cl->max[c] - cl->min[c]) / (cl->nboxes - 1);
    }
    else { 
      if(esl_FCompare(cl->min[c], 0., eslSMALLX1) != eslOK) esl_fatal("print_colorlegend(): colorlegend object is invalid, which_color[%d] is FALSE, but min[%d] is non-zero (%f).", c, c, cl->min[c]);
      if(esl_FCompare(cl->max[c], 0., eslSMALLX1) != eslOK) esl_fatal("print_colorlegend(): colorlegend object is invalid, which_color[%d] is FALSE, but max[%d] is non-zero (%f).", c, c, cl->max[c]);
    }
  }
  if(ndims == 0) esl_fatal("print_colorlegend(): colorlegend object is invalid, which_color[] is FALSE for all colors.");
  if(ndims >  2) esl_fatal("print_colorlegend(): colorlegend object is invalid, want to print %d dimensions, but max allowed is 2.", ndims);

  /* object is valid, print it */
  fprintf(fp, "%slegstart\n", "%");
  x = cl->x;
  y = cl->y;
  if(ndims == 1) { 
    /* print text for this legend */
    for(c = 0; c < NCMYK; c++) { 
      if(cl->which_color[c] == FALSE) continue;
      if(cl->text[c] != NULL) { 
	textlen = strlen(cl->text[c]);
	textsize = ESL_MIN(2. * (cl->boxsize * cl->nboxes) / (float) textlen, cl->boxsize);
	textsize = ESL_MAX(textsize, LEG_MINTEXTSIZE);
	/* back to black */
	fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
	fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", textsize);
	fprintf(fp, "(%s) %.4f %.4f lwstring\n", cl->text[c], x, y);
	/* reset font size to 8 */
	fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
      }
    }

    /* print heatmap */
    esl_vec_FSet(cur, NCMYK, 0.);
    for(c = 0; c < NCMYK; c++) { 
      if(cl->which_color[c] == FALSE) continue;
      cur[c] = cl->min[c];
      x = cl->x;
      y -= cl->boxsize * 1.5;
      for(b = 0; b < cl->nboxes; b++) { 
	fprintf(fp, "newpath\n");
	fprintf(fp, "  %.2f %.2f moveto", x, y);
	fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", cl->boxsize, cl->boxsize, (-1 * cl->boxsize));
	fprintf(fp, "  ");
	for(cp = 0;   cp < c;     cp++) fprintf(fp, "%.4f ", 0.);
	fprintf(fp, "%.4f ", (cur[c] - cl->min[c]) / (cl->max[c] - cl->min[c])); 
	for(cp = c+1; cp < NCMYK; cp++) fprintf(fp, "%.4f ", 0.);
	fprintf(fp, "setcmykcolor\n");
	fprintf(fp, "  fill\n");
	cur[c] += colstep[c];
	x += cl->boxsize;
      }
    }
    /* reset color to black */ 
    fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);

    /* print labels underneath heatmap */
    for(c = 0; c < NCMYK; c++) { 
      if(cl->which_color[c] == FALSE) continue;
      cur[c] = cl->min[c];
      x = cl->x;
      y -= cl->boxsize * 0.5;
      fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", ((float) cl->boxsize)/ 2.5);
      for(b = 0; b < cl->nboxes; b++) { 
	fprintf(fp, "(%3.2f) %.4f %.4f lwstring\n", cur[c], x, y);
	cur[c] += colstep[c];
	x += cl->boxsize;
      }
      /* reset font size to 8 */
      fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
    }
  }
  if(ndims == 2) { esl_fatal("print_colorlegend with 2 dimensions is not yet implemented."); }
  return eslOK;
}


/* Function: print_onecellcolorlegend()
 * 
 * Purpose:  Print a one cell color legend to an open file.
 * Return:   eslOK
 */
int 
print_onecellcolorlegend(FILE *fp, OneCellColorLegend_t *occl)
{
  float x, y;
  int cp;
  float textsize;

  /* object is valid, print it */
  fprintf(fp, "%sone cell legstart\n", "%");
  x = occl->x;
  y = occl->y;

  /* print cell */
  fprintf(fp, "newpath\n");
  fprintf(fp, "  %.2f %.2f moveto", x, y);
  fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", occl->boxsize, occl->boxsize, (-1 * occl->boxsize));
  fprintf(fp, "  ");
  for(cp = 0; cp < NCMYK; cp++) fprintf(fp, "%.4f ", occl->col[cp]);
  fprintf(fp, "setcmykcolor\n");
  fprintf(fp, "  fill\n");
  
  x += occl->boxsize * 1.5;

  /* print text for this legend */
  if(occl->text != NULL) { 
    textsize = 12;
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", textsize);
    fprintf(fp, "(%s) %.4f %.4f lwstring\n", occl->text, x, (y + (occl->boxsize * .25)));
    /* reset font size to 8 */
    fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  
  /* reset font size to 8 */
  fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
  return eslOK;
}


/* Function: print_scheme_colorlegend()
 * 
 * Purpose:  Print a scheme color legend to an open file.
 * Return:   eslOK
 */
int 
print_scheme_colorlegend(FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme)
{
  float x, y;
  int cp;
  int c;
  float textsize;

  /* object is valid, print it */
  fprintf(fp, "%sone cell legstart\n", "%");
  x = scl->x;
  y = scl->y;
  textsize = 16;

  /* print text for this legend */
  if(scl->text != NULL) { 
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", textsize);
    fprintf(fp, "(%s) %.4f %.4f lwstring\n", scl->text, x, (y + (scl->boxsize * .25)));
  }
  y -= scl->boxsize;

  fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", textsize);
  /* print cells and labels next to them */
  for(c = 0; c < scl->nbins; c++) { 
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", scl->boxsize, scl->boxsize, (-1 * scl->boxsize));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) fprintf(fp, "%.4f ", hc_scheme[c][cp]);
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += scl->boxsize * 1.5;
    y += scl->boxsize * 0.25;
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    if(c == scl->nbins-1) fprintf(fp, "(\\[%.3f-%.3f\\]) %.4f %.4f lwstring\n", scl->limits[c], scl->limits[c+1], x, y);
    else                  fprintf(fp, "(\\[%.3f-%.3f\\)) %.4f %.4f lwstring\n", scl->limits[c], scl->limits[c+1], x, y);
    x -= scl->boxsize * 1.5;
    y -= scl->boxsize * 0.25;
    y -= scl->boxsize;
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  
  /* reset font size to 8 */
  fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
  return eslOK;
}

/* Function: print_sspostscript()
 * 
 * Purpose:  Print a SS postscript data structure.
 * Return:   eslOK on success;
 *           eslEINCOMPAT if ps->npage == 0
 */
int
print_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps)
{
  int p, i, c, l;
  float title_fontsize;
  int do_square_mask, do_diamond_mask, do_x_mask, do_circle_mask;
  do_square_mask = do_diamond_mask = do_x_mask = do_circle_mask = FALSE;
  if(esl_opt_GetBoolean(go, "-d")) { do_diamond_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "-x")) { do_x_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "-c")) { do_circle_mask = TRUE; }
  else do_square_mask = TRUE;

  title_fontsize = TITLE_FONTSIZE;

  if(ps->npage == 0) ESL_FAIL(eslEINCOMPAT, errbuf, "print_sspostscript, ps->npage == 0\n");

  for(p = 0; p < ps->npage; p++) { 
    if(ps->regurgAA != NULL) {
      /* print from beginning up to title section */
      for(i = 0; i < ps->title_begin; i++)
	fprintf(fp, "%s", ps->regurgAA[i]);

      /* print title section */
      i = ps->title_begin;
      /* back to black */
      fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
      fprintf(fp, "/Helvetica findfont %.2f scalefont setfont\n", title_fontsize);
      /* to print with gutell's coords */
      fprintf(fp, "(\"%s\" page %d/%d) %s\n", command, p+1, ps->npage, ps->regurgAA[i]);
      /* to print with preset coords */
      /* WHAT I SHOULD DO: use Gutells x coord -200.00 units (x coord is first coord), but this
       * will require further parsing of the gutell line. */
      /*fprintf(fp, "(\"%s\" page %d/%d) -360.00 -200.00 lwstring\n", command, p+1, ps->npage);*/
      /* reset font size */
      fprintf(fp, "/Helvetica findfont 12.00 scalefont setfont\n");
      fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
      i++;
      /*fprintf(fp, "/Helvetica findfont 24.00 scalefont setfont\n");
	fprintf(fp, "(%s) %s\n", date, ps->regurgAA[i]);
	i++;
      */
      /*for(i = ps->title_begin; i < ps->title_begin + ps->title_nlines; i++) { 
	fprintf(fp, "TEST ");
	fprintf(fp, "%s", ps->regurgAA[i]);
	}*/

      /* print from title section to residue section */
      for(i = ps->title_begin + ps->title_nlines; i <ps->nregurg; i++)
	fprintf(fp, "%s", ps->regurgAA[i]);
    }
    
    /* print color legend, if any */
    if(ps->clAAA != NULL && ps->clAAA[p] != NULL) { 
      for(l = 0; l < ps->nclA[p]; l++) print_colorlegend(fp, ps->clAAA[p][l]);
    }

    /* print one cell color legends, if any */
    if(ps->occlAAA != NULL && ps->occlAAA[p] != NULL) { 
      for(l = 0; l < ps->nocclA[p]; l++) print_onecellcolorlegend(fp, ps->occlAAA[p][l]);
    }

    /* print scheme color legends, if any */
    if(ps->sclAA != NULL && ps->sclAA[p] != NULL) { 
      print_scheme_colorlegend(fp, ps->sclAA[p], hc_scheme[ps->sclAA[p]->scheme]);
    }

    if(ps->rcolAAA != NULL && ps->rcolAAA[p] != NULL) { 
      if(ps->maskAA[p] != NULL) { 
	for(c = 0; c < ps->clen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  fprintf(fp, "newpath\n");
	  if(ps->maskAA[p][c] == '0') { 
	    if(do_square_mask) { 
	      fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1. + 2., ps->ryA[c] -1. + 2.);
	      fprintf(fp, "  0 4 rlineto 4 0 rlineto 0 -4 rlineto closepath\n");
	    }
	    else if (do_diamond_mask) { 
	      fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1. + 4., ps->ryA[c] -1.);
	      fprintf(fp, "  -4 4 rlineto 4 4 rlineto 4 -4 rlineto closepath\n");
	    }
	    else if (do_circle_mask) { 
	      fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1. + 4., ps->ryA[c] -1.);
	      fprintf(fp, "  0 4 rlineto 4 0 rlineto 0 -4 rlineto closepath\n");
	    }
	    else if (do_x_mask) { 
	      fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1. + 4., ps->ryA[c] -1.);
	      fprintf(fp, "  0 4 rlineto 4 0 rlineto 0 -4 rlineto closepath\n");
	    }
	  }
	  else { 
	    fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1., ps->ryA[c] -1.);
	    fprintf(fp, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto closepath\n");
	  }
	  /*fprintf(fp, "  %.2f %.2f %.2f %.2f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);*/
	  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);
	  fprintf(fp, "  fill\n");
	}
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
      for(c = 0; c < ps->clen; c++) { 
	fprintf(fp, "(%c) %.2f %.2f lwstring\n", ps->rrAA[p][c], ps->rxA[c], ps->ryA[c]);
      }
    }
    fprintf(fp, "stroke\ngrestore\nm4showpage\n");
  }
  return eslOK;
}

/* read_template_file
 *
 * Read a Gutell postscript template file until we see a line
 * reading "%start". We'll regurgitate the pre-start info
 * from this file in all of the SS figs we draw.
 *
 * Returns:  eslOK on success.
 */
int
read_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  char          **regurgAA;
  int             nlines = 0;
  int             nalloc = 50;
  char           *newstr = NULL;
  char           *curstr = NULL;
  int             curlen = 0;
  void           *tmp;
  SSPostscript_t *ps = NULL;
  int             seen_residue_start = FALSE;
  int             seen_residue_end = FALSE;
  int             c = 0;
  float           x, y;
  int             ignore_flag = FALSE;
  int             in_title = FALSE;
  int             title_begin = 0; 
  int             title_end   = 0; 
  int             title_ntok = 0;
  /* Create the postscript object */
  ps = create_sspostscript();

  ESL_ALLOC(regurgAA, sizeof(char *) * nalloc);

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_template_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK && !seen_residue_start)
  {
    curlen = 0;
    if(in_title) { ignore_flag = TRUE; title_ntok = 0; } 
    while (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)  == eslOK) { 
      if(in_title && tok[(strlen(tok)-1)] == ')') ignore_flag = FALSE;
      if(strcmp(tok, "%title_start") == 0) { in_title = TRUE;  ignore_flag = TRUE;  title_begin = nlines; }
      if(strncmp(tok, "(5')", 4) == 0)     { in_title = FALSE; ignore_flag = FALSE; title_end   = nlines; }
      if(strncmp(tok, "%residue_start", 14) == 0) { seen_residue_start = TRUE; break; }
      if(!(in_title && ignore_flag)) { /* we're going to regurgitate this */
	if(in_title && title_ntok == 0) { title_ntok++; } /* skip final token of title lines */
	else { 
	  if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "read_template_file(), error (1) reading header template file.");
	  curlen += toklen;
	  if((status = esl_strcat(&curstr, curlen, " ",  1))      != eslOK) ESL_FAIL(status, errbuf, "read_template_file(), error (2) reading header template file.");
	  curlen += 1;
	  if(in_title) title_ntok++;
	}
      }
      /*if(in_title) printf("tok: %d line: %d tok: %s\n", title_ntok, (nlines - title_begin), tok);*/
      if(in_title && title_ntok == 2 && ((nlines - title_begin) == 0)) ps->titlex = atof(tok); /* titlex coord */
      if(in_title && title_ntok == 3 && ((nlines - title_begin) == 0)) ps->titley = atof(tok); /* titlex coord */
      if(in_title && title_ntok == 2 && ((nlines - title_begin) == 1)) ps->legx = atof(tok); /* legx coord */
      if(in_title && title_ntok == 3 && ((nlines - title_begin) == 1)) ps->legy = atof(tok); /* legx coord */
    }
    if(seen_residue_start) break; 
    if(!(in_title && ignore_flag)) { /* we're going to regurgitate this */
      if((status = esl_strcat(&curstr, curlen, "\n", 1)) != eslOK) ESL_FAIL(status, errbuf, "read_template_file(), error (3) reading header template file.");
      curlen += 1;

      ESL_ALLOC(newstr, sizeof(char) * (curlen+1));
      strcpy(newstr, curstr);
      regurgAA[nlines++] = newstr;
      if(nlines == nalloc) { nalloc += 50; ESL_RALLOC(regurgAA, tmp, sizeof(char *) * nalloc); }
      free(curstr);
      curstr = NULL;
    }
  }
  ps->regurgAA = regurgAA;
  ps->nregurg = nlines;

  /* We're done reading the header information that we'll regurgitate, now we need
   * the sequence residue coordinates */
  while (esl_fileparser_NextLine(efp) == eslOK && !seen_residue_end)
  {
    if(c == ps->clen) { 
      ps->clen += 100;
      ESL_RALLOC(ps->rxA, tmp, sizeof(float) * ps->clen);
      ESL_RALLOC(ps->ryA, tmp, sizeof(float) * ps->clen);
      assert(ps->rrAA == NULL);
    }
    /* example line:
     *(A) 61.30 -831.00 lwstring
     */
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) esl_fatal("Failed to read residue on line %d of postscript template file %s\n", efp->linenumber, filename);
    if(strncmp(tok, "stroke", 5) == 0) { seen_residue_end = TRUE; break; }
    /* residue is not used b/c we overwrite it (ex. (A)) */
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) esl_fatal("Failed to read x coord on line %d of postscript template file %s\n", efp->linenumber, filename);
    x = atof(tok);
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) esl_fatal("Failed to read y coord on line %d of postscript template file %s\n", efp->linenumber, filename);
    y = atof(tok);
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) esl_fatal("Failed to read 'lwstring' on line %d of postscript template file %s\n", efp->linenumber, filename);
    
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslEOL) esl_fatal("Failed to read EOL on line %d of postscript template file %s\n", efp->linenumber, filename);
    
    ps->rxA[c] = x;
    ps->ryA[c] = y;
    c++;
  }      
  ps->clen = c;
  if(title_begin == 0) esl_fatal("Failed to read title section in postscript template file %s. Add \"%stitle_start\" line before \"/Helvetica findfont 24.00 scalefont setfont\" line.", filename, "%");
  ps->title_begin  = title_begin;
  ps->title_nlines = title_end - title_begin;
  /*printf("begin: %d nlines: %d\n", ps->title_begin, ps->title_nlines);
    printf("titlex: %f titley: %f\n", ps->titlex, ps->titley);
    printf("legx: %f legy: %f\n", ps->legx, ps->legy);*/
  esl_fileparser_Close(efp);

  *ret_ps = ps;
  return eslOK;
  
 ERROR:
  return eslEMEM;
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

  if((status = addpages_sspostscript(ps, msa->nseq)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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
infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  double **obs, *ent, *bg;

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    /*ESL_ALLOC(ps->clAAA[p],    sizeof(ColorLegend_t **) * 1);*/
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(ent, sizeof(double) * ps->clen);
  ESL_ALLOC(obs, sizeof(double *) * ps->clen);
  ESL_ALLOC(bg, sizeof(double) * msa->abc->K);
  esl_vec_DSet(bg, msa->abc->K, 1./(msa->abc->K));
  
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    ESL_ALLOC(obs[cpos], sizeof(double) * msa->abc->K);
    esl_vec_DSet(obs[cpos], msa->abc->K, 0.);
  }

  pp = orig_npage;

  /* add color legend */
  float *limits;
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.4;
  limits[2] = 0.8;
  limits[3] = 1.2;
  limits[4] = 1.6;
  limits[5] = 1.99;
  limits[6] = 2.00;
  ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);

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
    esl_vec_DNorm(obs[cpos], msa->abc->K);
    ent[cpos] = esl_vec_DEntropy(bg, msa->abc->K) - esl_vec_DEntropy(obs[cpos], msa->abc->K);

    if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent[cpos], ps->sclAA[pp])) != eslOK) return status;

    ps->rrAA[pp][cpos] = ' ';
    //ps->rrAA[pp][cpos] = (esl_FCompare(ent[cpos], 0., eslSMALLX1) == eslOK) ? '-' : ' ';
  }

  /* add text to legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * strlen("information content (bits) (total: 1000000.00 bits)"));
  sprintf(text, "information content (bits) (total: %.2f bits)", esl_vec_DSum(ent, ps->clen));
  add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  free(text);

  if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);

  free(ent);
  for(cpos = 0; cpos < ps->clen; cpos++) free(obs[cpos]);
  free(obs);
  free(bg);
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
delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  int *dct;
  int *dct_internal;
  int *fA, *lA;

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(int *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
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
  float *limits;
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.01;
  limits[2] = 0.2;
  limits[3] = 0.4;
  limits[4] = 0.6;
  limits[5] = 0.8;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);

  if(do_all) { 
    /* draw delete page with all deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct[cpos]) / (float) (msa->nseq), ps->sclAA[pp])) != eslOK) return status;
    }
  }
  else { /* do_all is FALSE, draw delete page with only internal deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct_internal[cpos]) / (float) (msa->nseq), ps->sclAA[pp])) != eslOK) return status;
    }
  }

  /* add color legend */
  char *text;
  if(do_all) { 
    ESL_ALLOC(text, sizeof(char) * strlen("fraction seqs w/deletes ('-'=0 deletes; avg/seq: 1000.00)"));
    sprintf(text, "fraction seqs w/deletes ('-'=0 deletes; avg/seq: %.2f)", (float) esl_vec_ISum(dct, ps->clen) / (float) msa->nseq);
    add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  }
  else { /* !do_all, only internal deletes counted */
    ESL_ALLOC(text, sizeof(char) * strlen("fraction seqs w/internal deletes ('-'=0; avg/seq: 1000.00)"));
    sprintf(text, "fraction seqs w/internal deletes ('-'=0; avg/seq: %.2f)", (float) esl_vec_ISum(dct_internal, ps->clen) / (float) msa->nseq);
    add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  }
  free(text);
  free(dct);
  free(dct_internal);
  free(fA);
  free(lA);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "delete_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: insert_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, with colors in log 
 *           scale indicating the fraction of sequences with inserts after each
 *           position, and numbers indicating the median length of inserts in those
 *           sequences that have inserts at each position. Positions with 0 inserts
 *           in all sequences are marked '-' with no color. Positions with median
 *           length 10 or greater are marked with '*'. 
 *           
 * Return:   eslOK on success.
 */
int
insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, float ***hc_scheme, int hc_scheme_idx, int hc_nbins)
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

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
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
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.01;
  limits[2] = 0.2;
  limits[3] = 0.4;
  limits[4] = 0.6;
  limits[5] = 0.8;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);

  for(cpos = 1; cpos <= ps->clen; cpos++) { 
    if(nseq_ict[cpos] == 0) { 
      res = '-';
      col = 0.0;
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
    }
    ps->rrAA[pp][(cpos-1)] = res;
    if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][(cpos-1)], NCMYK, hc_scheme[hc_scheme_idx], col, ps->sclAA[pp])) != eslOK) return status;
  }

  /* add color legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * strlen("fraction seqs w/inserts; 'N' = median size, if N=*, N > 10; avg/seq: 1000.00"));
  sprintf(text, "fraction seqs w/inserts; 'N' = median size, if N=*, N > 10; avg/seq: %.2f", (float) esl_vec_ISum(total_ict, ps->clen+1) / (float) msa->nseq);
  add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  free(text);
  
  for(i = 0; i < ps->clen; i++) free(ict[i]);
  free(ict);
  free(total_ict);
  free(nseq_ict);
  free(med_ict);

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "insert_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: addpages_sspostscript()
 * 
 * Purpose:  Add and initialize blank pages to a postscript object.
 */ 
static int 
addpages_sspostscript(SSPostscript_t *ps, int ntoadd)
{
  int status;
  void *tmp;
  int p;

  if(ps->npage == 0) { 
    assert(ps->rrAA    == NULL);
    assert(ps->rcolAAA == NULL);
    ESL_ALLOC(ps->rrAA,    sizeof(char *)   * (ps->npage + ntoadd));
    ESL_ALLOC(ps->rcolAAA, sizeof(float **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->clAAA,   sizeof(ColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->nclA,    sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->occlAAA, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->nocclA,  sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->sclAA,   sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->maskAA,  sizeof(char **) * (ps->npage + ntoadd));
  }
  else { 
    assert(ps->rrAA    != NULL);
    assert(ps->rcolAAA != NULL);
    ESL_RALLOC(ps->rrAA,    tmp, sizeof(char *)   * (ps->npage + ntoadd));
    ESL_RALLOC(ps->rcolAAA, tmp, sizeof(float **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->clAAA,   tmp, sizeof(ColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->nclA,    tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->occlAAA, tmp, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->nocclA,  tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->sclAA,   tmp, sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->maskAA,  tmp, sizeof(char **) * (ps->npage + ntoadd));
  }
  for(p = ps->npage; p < (ps->npage + ntoadd); p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->clAAA[p]   = NULL;
    ps->nclA[p]    = 0;
    ps->occlAAA[p] = NULL;
    ps->nocclA[p]  = 0;
    ps->sclAA[p]   = NULL;
    ps->maskAA[p]  = NULL;
  }
  ps->npage += ntoadd;
  ps->cur_legx   = ps->legx;
  ps->cur_legy   = ps->legy;
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


/* Function: posteriors_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with info on posterior probabilities in the MSA.
 * Return:   eslOK on success.
 */
static int
posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins)
{
  if(mask != NULL) esl_fatal("posteriors with mask not yet upgraded.");

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
  char *text;    /* text for color legends */

  if(msa->rf == NULL) esl_fatal("No RF annotation in alignment");

  if(esl_opt_GetBoolean(go, "--p-avg"))  { do_avg = TRUE;  new_npage += 1; navg_page = orig_npage; }
  if(esl_opt_GetBoolean(go, "--p-indi")) { do_indi = TRUE; nfirst_indi_page = orig_npage + new_npage; new_npage += msa->nseq; }

  if((status = addpages_sspostscript(ps, new_npage)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    /*if(mask == NULL) ESL_ALLOC(ps->clAAA[p],    sizeof(ColorLegend_t **) * 1);
      else             ESL_ALLOC(ps->clAAA[p],    sizeof(ColorLegend_t **) * 2);*/
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
    ESL_FAIL(eslEINVAL, errbuf, "--p-avg and --p-indi require \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s (from infernal v1.x\'s cmalign).\n", esl_opt_GetArg(go,1));
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
    nongap_s = nongaprf_s = 0;
    sum_s    = sumrf_s = 0.;
    if(do_indi) { /* add color legend for this sequence */
      ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);
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
	    if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], prob, ps->sclAA[pp])) != eslOK) return status;
	    ps->rrAA[pp][cpos] = ' ';
	  }
	}
	//if(prob >= pthresh) athresh_c[c]++;
      }
      else if (do_indi) { /* gap, if it's a consensus column, draw blank square */
	if(a2c_map[c] != -1) { /* consensus position */
	  cpos = a2c_map[c];
	  ps->rcolAAA[pp][cpos][0] = BLANKCYAN; 
	  ps->rcolAAA[pp][cpos][1] = BLANKMAGENTA; 
	  ps->rcolAAA[pp][cpos][2] = BLANKYELLOW; 
	  ps->rcolAAA[pp][cpos][3] = BLANKBLACK; 
	  ps->rrAA[pp][cpos] = ' ';
	}
      }
    } /* done with this sequence */
    if(do_indi) { 
      avg_s   =  (float) sum_s / (float) nongap_s;
      avgrf_s =  (float) sumrf_s / (float) nongaprf_s;
      ESL_ALLOC(text, sizeof(char) * strlen("avg posterior probability; 1.000 (RF) 1.000 (all)"));
      sprintf(text, "avg posterior probability; %.3f (RF) %.3f (all)", avgrf_s, avg_s);
      add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
      free(text);
      pp++;
    }
  } /* done with all sequences */

  if(do_avg) { /* add average colors */
    pp = navg_page;
    ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);
    for(c = 0; c < msa->alen; c++) { 
      if(a2c_map[c] != -1) { /* cons position */
	cpos = a2c_map[c]; 
	if(nongap_c[c] > 0) { 
	  avgrf_c = sum_c[c] /= (float) nongap_c[c];
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], avgrf_c, ps->sclAA[pp])) != eslOK) return status;
	}
	else { 
	  ps->rcolAAA[pp][cpos][0] = BLANKCYAN; 
	  ps->rcolAAA[pp][cpos][1] = BLANKMAGENTA; 
	  ps->rcolAAA[pp][cpos][2] = BLANKYELLOW; 
	  ps->rcolAAA[pp][cpos][3] = BLANKBLACK; 
	}
	ps->rrAA[pp][cpos] = ' ';
      }
    }
    /* add color legend */
    ESL_ALLOC(text, sizeof(char) * strlen("avg posterior probability; 1.000 (RF) 1.000 (all)"));
    sprintf(text, "avg posterior probability; %.3f (RF) %.3f (all)", 
	    (esl_vec_FSum(sumrf_c, msa->alen) / (float) (esl_vec_ISum(nongaprf_c, msa->alen))),
	    (esl_vec_FSum(sum_c, msa->alen)   / (float) (esl_vec_ISum(nongap_c, msa->alen))));
    add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
    free(text);
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
 * file and return it as *ret_mask.
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));
  for(n = 0; n < toklen; n++) mask[n] = tok[n];
  mask[n] = '\0';

  *ret_mask = mask;
  *ret_masklen= toklen;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}

/* Function: draw_file2sspostscript()
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
	if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");
	
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
	  if(strlen(s) != 1) esl_fatal("Read multi-character string (%s) for consensus residue %d on line %d of drawfile %s\n", s, cpos, efp->linenumber, dfile);
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
structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins)
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
  

  if(msa->ss_cons == NULL) ESL_FAIL(status, errbuf, "--struct requires #=GC SS_cons annotation in the alignment.");
  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
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
  ps->sclAA[pp] = create_scheme_colorlegend(ps, hc_scheme_idx, hc_nbins, LEG_ONED_BOXSIZE, NULL, limits);

  /* get observed residues at each cpos */
  for(i = 0; i < msa->nseq; i++) { 
    cpos = 0;
    for(apos = 0; apos < msa->alen; apos++) {
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { /* seq i is not a gap at cpos */
	  /* only count base paired positions for which both left and right half are not gaps, 
	   * check if we're base paired */
	  if(ct[apos+1] != 0) { 
	    if(ct[apos+1] > (apos+1)) { /* cpos is left half of base pair */
	      /* check if right half is a gap */
	      rapos = ct[apos+1]-1; 
	      if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][rapos])) { /* seq i is not a gap at right half */
		esl_abc_DCount(msa->abc, obs[cpos], esl_abc_DigitizeSymbol(msa->abc, msa->aseq[i][apos]), 1.);
		rcpos = a2c_map[rapos+1]-1;
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
    apos  = c2a_map[cpos+1]-1;
    if(ct[apos+1] != 0) { 
      esl_vec_DNorm(obs_p[cpos], msa->abc->K * msa->abc->K);

      rapos = ct[apos+1]-1; 
      rcpos = a2c_map[rapos+1]-1;

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
    if(ent_p[cpos] < (-1. * eslSMALLX1)) { /* single stranded base, paint grey */
      ps->rcolAAA[pp][cpos][0] = BLANKCYAN; 
      ps->rcolAAA[pp][cpos][1] = BLANKMAGENTA; 
      ps->rcolAAA[pp][cpos][2] = BLANKYELLOW; 
      ps->rcolAAA[pp][cpos][3] = BLANKBLACK; 
      ent_p[cpos] = 0.; /* impt to set to 0., so esl_vec_DSum(ent_p... call below to calc total struct info is accurate */
    }
    else if(mask == NULL) { 
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent_p[cpos], ps->sclAA[pp])) != eslOK) return status;
    }
    else { 
      if(mask[cpos] == '0') { 
	ps->rcolAAA[pp][cpos][0] = 0.0;
	ps->rcolAAA[pp][cpos][1] = ent_p[cpos];
	ps->rcolAAA[pp][cpos][2] = ent_p[cpos];
	ps->rcolAAA[pp][cpos][3] = 0.0;
      }
      else if(mask[cpos] == '1') { 
	ps->rcolAAA[pp][cpos][0] = ent_p[cpos];
	ps->rcolAAA[pp][cpos][1] = ent_p[cpos];
	ps->rcolAAA[pp][cpos][2] = 0.0; 
	ps->rcolAAA[pp][cpos][3] = 0.0; 
      }
      else ESL_FAIL(eslEINVAL, errbuf, "--mask mask char number %d is not a 1 nor a 0, but a %c\n", cpos, mask[cpos]);
    }
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add text to the legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * strlen("structural info content per basepaired posn (total: 1000.00 bits)"));
  sprintf(text,  "structural info content per basepaired posn (total: %.2f bits)", esl_vec_DSum(ent_p, ps->clen) * 2.);
  add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  free(text);

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


/* Function: phylosignal_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, greyscale (CMYK K value) squares indicating
 *           the bits of phylogenetic signal/information per column. 
 *         
 *           If mask is non-NULL, columns that are *outside* the mask are colored red instead of 
 *           black.
 *
 * Return:   eslOK on success.
 */
int
phylosignal_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  double **obs;
  double ent;
  double summed_ent;
  double inmask_summed_ent;

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    if(mask == NULL) ESL_ALLOC(ps->clAAA[p],    sizeof(ColorLegend_t **) * 1);
    else             ESL_ALLOC(ps->clAAA[p],    sizeof(ColorLegend_t **) * 2);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(obs, sizeof(double *) * ps->clen);
  
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    ESL_ALLOC(obs[cpos], sizeof(double) * msa->abc->K);
    esl_vec_DSet(obs[cpos], msa->abc->K, 0.);
  }

  pp = orig_npage;
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
  summed_ent = inmask_summed_ent;
  for(cpos = 0; cpos < ps->clen; cpos++) { 
    esl_vec_DNorm(obs[cpos], msa->abc->K);
    ent = esl_vec_DEntropy(obs[cpos], msa->abc->K);

    if(mask == NULL || mask[cpos] == '1') { 
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 0.0;
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = ent/2.;
      inmask_summed_ent += ent;
    }
    else if (mask[cpos] == '0') { /* mask != NULL && mask[cpos] != '0' */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = ent/2.; 
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = 0.0; 
    }
    else ESL_FAIL(eslEINVAL, errbuf, "--mask mask char number %d is not a 1 nor a 0, but a %c\n", cpos, mask[cpos]);
    summed_ent += ent;
    ps->rrAA[pp][cpos] = (esl_FCompare(ent, 0., eslSMALLX1) == eslOK) ? '-' : ' '; /* mark blank squares with a '-' */
  }
  printf("       Consensus columns: %d\n", ps->clen);
  printf("phylogenetic information: %.2f bits\n", summed_ent);

  if(mask != NULL) { 
    printf("             within mask: %.2f bits (%.4f)\n", inmask_summed_ent, inmask_summed_ent / summed_ent);
  }

  /* add color legend */
  char *text;
  if(mask == NULL) { 
    ESL_ALLOC(text, sizeof(char) * strlen("entropy (phylogenetic signal) in bits (total: 10000.00)"));
    sprintf(text, "entropy (phylogenetic signal) in bits (total: %.2f)", summed_ent);
    ps->clAAA[pp][0] = create_one_dim_colorlegend(ps, IBLACK, 0., 2., LEG_ONED_BOXSIZE, LEG_ONED_NBOXES, text);
    ps->nclA[pp] = 1;
    free(text);
  }
  else { 
    ESL_ALLOC(text, sizeof(char) * strlen("within mask  entropy (phylogenetic signal) in bits (total: 10000.00 (1000))"));
    sprintf(text, "within mask  entropy (phylogenetic signal) in bits (total: %.2f (%.2f%s))", inmask_summed_ent, 100 * inmask_summed_ent / summed_ent, "%");
    ps->clAAA[pp][0] = create_one_dim_colorlegend(ps, IBLACK,   0., 2., LEG_ONED_BOXSIZE, LEG_ONED_NBOXES, text);
    free(text);

    ESL_ALLOC(text, sizeof(char) * strlen("outside mask entropy (phylogenetic signal) in bits (total: 10000.00 (1000))"));
    sprintf(text, "outside mask entropy (phylogenetic signal) in bits (total: %.2f (%.2f%s))", summed_ent - inmask_summed_ent, 100. - (100 * inmask_summed_ent / summed_ent), "%");
    ps->clAAA[pp][1] = create_one_dim_colorlegend(ps, IMAGENTA,   0., 2., LEG_ONED_BOXSIZE, LEG_ONED_NBOXES, text);
    free(text);

    ps->nclA[pp] = 2;
  }

  for(cpos = 0; cpos < ps->clen; cpos++) free(obs[cpos]);
  free(obs);
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "phylosignal_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
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
colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_inside_mask = 0;
  int ncols_outside_mask = 0;
  
  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->occlAAA[p],    sizeof(ColorLegend_t **) * 2);
    for(c = 0; c < ps->clen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  for(cpos = 0; cpos < ps->clen; cpos++) { 
    if(mask[cpos] == '1') { 
      /* green, 100% cyan, 100% yellow */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 0.0;
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = 1.0;
      ncols_inside_mask++;
    }
    else if(mask[cpos] == '0') {
      /* red, 100% magenta, 100% yellow */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 1.0;
      ps->rcolAAA[pp][cpos][2] = 1.0;
      ps->rcolAAA[pp][cpos][3] = 0.0; 
      ncols_outside_mask++;
    }
    else ESL_FAIL(eslEINVAL, errbuf, "--mask mask char number %d is not a 1 nor a 0, but a %c\n", cpos, mask[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  char *text;
  float col[NCMYK];
  ESL_ALLOC(text, sizeof(char) * strlen("columns included within mask (1000 of 1000 (1.000))"));
  sprintf(text, "columns included within mask (%4d of %4d (%.3f))", ncols_inside_mask, ps->clen, (float) ncols_inside_mask / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[IBLACK] = 1.0;
  ps->occlAAA[pp][0] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * strlen("columns excluded from  mask (1000 of 1000 (1.000))"));
  sprintf(text, "columns excluded from  mask (%4d of %4d (%.3f))", ncols_outside_mask, ps->clen, (float) ncols_outside_mask / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[IMAGENTA] = 1.0;
  col[IYELLOW] = 1.0;
  ps->occlAAA[pp][1] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

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
diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_in_both = 0;
  int ncols_out_both = 0;
  int ncols_in_1_out_2 = 0;
  int ncols_out_1_in_2 = 0;
  int textsize = 0;

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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
      /* black */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 0.0;
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = 1.0;
      ncols_in_both++;
    }
    else if(mask1[cpos] == '1' && mask2[cpos] == '0') {
      /* red */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 1.0;
      ps->rcolAAA[pp][cpos][2] = 1.0;
      ps->rcolAAA[pp][cpos][3] = 0.0; 
      ncols_in_1_out_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '1') {
      /* cyan */
      ps->rcolAAA[pp][cpos][0] = 1.0;
      ps->rcolAAA[pp][cpos][1] = 0.0;
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = 0.0; 
      ncols_out_1_in_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '0') {
      /* cyan */
      ps->rcolAAA[pp][cpos][0] = 0.0;
      ps->rcolAAA[pp][cpos][1] = 0.0;
      ps->rcolAAA[pp][cpos][2] = 0.0;
      ps->rcolAAA[pp][cpos][3] = 0.2; 
      ncols_out_both++;
    }
    else if(mask1[cpos] != '0' && mask1[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-col char number %d is not a 1 nor a 0, but a %c\n", cpos, mask1[cpos]);
    else if(mask2[cpos] != '0' && mask2[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-diff char number %d is not a 1 nor a 0, but a %c\n", cpos, mask2[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  char *text;
  float col[NCMYK];
  ESL_ALLOC(text, sizeof(char) * (strlen("columns included within both masks (1000 of 1000 (1.000))")+1));
  sprintf(text, "columns included within both masks (%4d of %4d (%.3f))", ncols_in_both, ps->clen, (float) ncols_in_both / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[IBLACK] = 1.0;
  ps->occlAAA[pp][0] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

  textsize  = strlen("columns incl. in --mask-col mask but not --mask-diff mask (1000 of 1000 (1.000))");
  /*textsize += strlen(esl_opt_GetString(go, "--mask-col"));
    textsize += strlen(esl_opt_GetString(go, "--mask-diff"));*/
  ESL_ALLOC(text, sizeof(char) * (textsize+1));
  /*sprintf(text, "columns included in mask %s but not mask %s (%4d of %4d (%.3f))", esl_opt_GetString(go, "--mask-col"), esl_opt_GetString(go, "--mask-diff"), ncols_in_1_out_2, ps->clen, (float) ncols_in_1_out_2 / (float) ps->clen);*/
  sprintf(text, "columns incl. in --mask-col mask but not mask --mask-diff mask (%4d of %4d (%.3f))", ncols_in_1_out_2, ps->clen, (float) ncols_in_1_out_2 / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[IMAGENTA] = 1.0;
  col[IYELLOW] = 1.0;
  ps->occlAAA[pp][1] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

  textsize  = strlen("columns incl. in --mask-diff mask but not --mask-col mask (1000 of 1000 (1.000))");
  textsize  += 22; /* len of --mask-col plus --mask-diff plus 1*/
  /*textsize += strlen(esl_opt_GetString(go, "--mask-col"));
    textsize += strlen(esl_opt_GetString(go, "--mask-diff"));*/
  ESL_ALLOC(text, sizeof(char) * (textsize+1));
  /*sprintf(text, "columns included in mask %s but not mask %s (%4d of %4d (%.3f))", esl_opt_GetString(go, "--mask-diff"), esl_opt_GetString(go, "--mask-col"), ncols_out_1_in_2, ps->clen, (float) ncols_out_1_in_2 / (float) ps->clen);*/
  sprintf(text, "columns included in --mask-diff mask but not --mask-col mask (%4d of %4d (%.3f))", ncols_out_1_in_2, ps->clen, (float) ncols_out_1_in_2 / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[ICYAN] = 1.0;
  ps->occlAAA[pp][2] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * (strlen("columns excluded from both masks (1000 of 1000 (1.000))")+1));
  sprintf(text, "columns excluded from both masks (%4d of %4d (%.3f))", ncols_out_both, ps->clen, (float) ncols_out_both / (float) ps->clen);
  esl_vec_FSet(col, NCMYK, 0.);
  col[IBLACK] = 0.2;
  ps->occlAAA[pp][3] = create_one_cell_colorlegend(ps, col, LEG_ONED_BOXSIZE, text);
  free(text);

  ps->nocclA[pp] = 4;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "diffmask_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call esl-ssudraw
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
set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl)
{
  float min, max;
  int ci, bi;
  min = scl->limits[0];
  max = scl->limits[scl->nbins];
  if((min-val) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f < min: %.4f\n", val, min); }
  if((val-max) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f > max: %.4f\n", val, max); }

  bi = 0;
  while((val > scl->limits[bi+1]) && (bi <= (scl->nbins-1))) { bi++; }    
  printf("%.3f %d (%.3f)\n", val, bi, scl->limits[bi+1]);
  for(ci = 0; ci < ncolvals; ci++) { vec[ci] = scheme[bi][ci]; }
  return eslOK;
}
  
