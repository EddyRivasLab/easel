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
#define LEG_MINTEXTSIZE 10
#define LEGX_OFFSET 24.
#define LEGY_OFFSET -24.
#define LEG_FONT "Courier-Bold"

#define DEFAULT_FONT "Courier-Bold"
#define RESIDUE_FONT "Helvetica-Bold"

#define SS_BOXSIZE 8.

#define RESIDUES_DEFAULT_FONTSIZE 8.
#define HUNDREDS_DEFAULT_FONTSIZE 8.
#define TITLE_DEFAULT_FONTSIZE 24.
#define TICKS_DEFAULT_LINEWIDTH 2.
#define BP_DEFAULT_LINEWIDTH 1.

#define POSTSCRIPT_PAGEWIDTH 612.
#define POSTSCRIPT_PAGEHEIGHT 792.
#define PAGE_TOPBUF 18.
#define PAGE_SIDEBUF 18.
#define PAGE_BOTBUF 18.

/* Structure: scheme_color_legend
 * Incept:    EPN, Thu Jun 25 20:20:38 2009
 *
 * Parameters describing a one-dimensional legend of colors
 * from a preset scheme for use in a SSPostscript_t data structure.
 */
typedef struct scheme_color_legend_s {
  int    scheme;            /* preset color scheme index */
  int    nbins;             /* number of colors (bins) in this scheme */
  char   *text;             /* text for legend, a single string */
  float  boxsize;           /* size of box for each residue */
  float *limits;            /* [nbins+1] limits for each bin, limits[0] is min value we would expect to see, limits[nbins] is max */
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
  int     npage;        /* number of pages in eventual postscript */
  char   *modelname;    /* name of model, read from template file */
  char  **titleA;       /* text for the generic title that will appear */
  int     ntitle;       /* number of lines of title information */
  float   titlex;       /* x coordinate (bottom left corner) of title area */
  float   titley;       /* y coordinate (bottom left corner) of title area */
  float   legx;         /* x coordinate (bottom left corner) of legend area */
  float   legy;         /* y coordinate (bottom left corner) of legend area */
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
  char   **maskAA;      /* [0..npage-1][0..clen-1] mask, columns which are '0' get drawn differently */
  int      nalloc;      /* number of elements to add to arrays when reallocating */
  int     *msa_ct;      /* [1..ps->clen] CT array for msa this postscript corresponds to, 
			 * msa_ct[i] is the position that consensus residue i base pairs to, or 0 if i is unpaired. */
  int      msa_nbp;     /* number of bps read from current MSA (in msa_ct), should equal nbp, but only if bps read from template file */
  int      msa_idx;     /* msa index we're currently on in MSA file */
} SSPostscript_t;

static SSPostscript_t *create_sspostscript();
static int  setup_sspostscript(SSPostscript_t *ps, char *errbuf);
static OneCellColorLegend_t *create_onecell_colorlegend(float *cmykA, float boxsize, char *text);
static SchemeColorLegend_t *create_scheme_colorlegend(int scheme, int ncols, float boxsize, char *text, float *limits);
static int  add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text);
static int  draw_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps);
static int  draw_onecell_colorlegend(FILE *fp, OneCellColorLegend_t *occl, SSPostscript_t *ps, int occl_idx);
static int  draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page, int do_mask);
static void free_sspostscript(SSPostscript_t *ps);
static int  addpages_sspostscript(SSPostscript_t *ps, int ntoadd);
static int  map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int **ret_a2c_map, int *ret_clen);
static int  parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps);
static int  parse_modelname_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_scale_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf);
static int  parse_regurgitate_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_text_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_lines_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  validate_justread_sspostscript(SSPostscript_t *ps, char *errbuf);
static int  validate_and_update_sspostscript_given_msa(SSPostscript_t *ps, ESL_MSA *msa, char *errbuf, int msa_idx);
static int  individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  rf_seq_sspostscript (const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_bins, float **hc_onecell, int ss_idx, int zerores_idx);
static int  delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float **hc_onecell, int incmask_idx, int excmask_idx);
static int  diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2, float **hc_onecell, int incboth_idx, int inc1_idx, int inc2_idx, int excboth_idx);
static int  compare_ints(const void *el1, const void *el2);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen);
static int  drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps);
static void PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, float wt);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  get_date(char *errbuf, char **ret_date);
static int  set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl);
static int  set_onecell_values(char *errbuf, float *vec, int ncolvals, float *onecolor);
static int  add_mask_to_ss_postscript(SSPostscript_t *ps, int page, char *mask);
static int  draw_masked_block(FILE *fp, float x, float y, float *colvec, int do_circle_mask, int do_square_mask, int do_x_mask, int do_border, float boxsize);

static char banner[] = "draw Gutell based postscript SSU secondary structure diagrams.";
static char usage[]  = "[options] <msafile> <Gutell SS postscript template> <output postscript file name>\n\
The <msafile> must be in Stockholm format.";

#define MASKTYPEOPTS "-d,-c,-x" /* exclusive choice for mask types */

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-q",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "DO NOT create SS info content diagram (on by default)", 0 },
  { "-s",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create SS diagram for each sequence in the alignment",    0 },
  { "-u",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as squares", 1 },
  { "-x",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as x's", 1 },
  { "-a",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask and -u or -x, draw alternative mask style", 1 },
  { "--rf",     eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create SS diagram for RF sequence",    1 },
  { "--struct", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create structural info content SS diagram",    1 },
  { "--p-avg",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create average posterior probability SS diagram",1 },
  { "--p-indi", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create posterior probability diagram for each sequence",1 },
  { "--ins",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,         "create insert SS diagram",    1 },
  { "--dall",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create delete diagram w/all deletions (incl. terminal deletes)",    1 },
  { "--dint",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "create delete diagram w/only internal (non-terminal) deletions", 1 },
  { "--mask",    eslARG_INFILE, NULL, NULL, NULL, NULL,NULL,NULL,            "for all diagrams, mark masked columns from mask in <f>", 1 },
  { "--mask-col",eslARG_NONE, NULL, NULL, NULL, NULL,"--mask", NULL,        "w/--mask create black/orange diagram denoting masked columns", 1 },
  { "--mask-diff",eslARG_INFILE,NULL, NULL, NULL, NULL,"--mask",NULL,        "with --mask-col <f1>, compare mask in <f1> to mask in <f>", 1},
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
  if((status = parse_template_file(templatefile, go, errbuf, &ps) != eslOK)) esl_fatal(errbuf);
  /* determine position for title and legend */
  if((status = setup_sspostscript(ps, errbuf) != eslOK)) esl_fatal(errbuf);
  
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
  if(! esl_opt_IsDefault(go, "--mask-diff")) { 
      if((status = read_mask_file(esl_opt_GetString(go, "--mask-diff"), errbuf, &mask2, &masklen2)) != eslOK) esl_fatal(errbuf);
      if(masklen != masklen2) esl_fatal("Mask in %f length (%d) differs from mask in %f (%d)!", esl_opt_GetString(go, "--mask"), masklen, esl_opt_GetString(go, "--mask-diff"), masklen2);
  }

  /***********************************/
  /* allocate and fill predefined one-cell colors, these are hardcoded */
  float **hc_onecell;
  int z;
  ESL_ALLOC(hc_onecell, sizeof(float *) * NOC);
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
  int     *hc_nbins;
  float ***hc_scheme;
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
  hc_scheme[3][0][0] = 0.00; hc_scheme[3][0][1] = 0.94; hc_scheme[3][0][2] = 1.00; hc_scheme[3][0][3] = 0.00; /*red*/
  /***************************************************************/
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      msa->abc = abc;
      if(msa->rf == NULL) esl_fatal("MSA number: %d in %s does not have RF annotation.", nali, alifile);
      clen = 0;
      for(apos = 0; apos < msa->alen; apos++) if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) clen++;
      if(ps->clen == 0)    esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != template file consensus length of %d.", nali, clen, ps->clen);
      if(clen != ps->clen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != template file consensus length of %d.", nali, clen, ps->clen);
      if(mask != NULL && ps->clen != masklen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != lane mask length of %d from mask file %s.", nali, clen, masklen, esl_opt_GetString(go, "--mask"));

      if((status = validate_and_update_sspostscript_given_msa(ps, msa, errbuf, nali)) != eslOK) esl_fatal(errbuf);

      if(! esl_opt_GetBoolean(go, "-q")) { 
	if((status = infocontent_sspostscript(go, errbuf, ps, msa, mask, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--struct")) { 
	if((status = structural_infocontent_sspostscript(go, errbuf, ps, msa, mask, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, DARKGREYOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--ins")) { /* make a new postscript page marking insertions */
	if((status = insert_sspostscript(go, errbuf, ps, msa, mask, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dall")) { /* make a new postscript page marking all deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, mask, TRUE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--dint")) { /* make a new postscript page marking internal deletes */
	if((status = delete_sspostscript(go, errbuf, ps, msa, mask, FALSE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--rf")) { /* make a new postscript page for the RF sequence in the alignment */
	if((status = rf_seq_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
      }
      int do_post = (esl_opt_GetBoolean(go, "--p-avg")) ? TRUE : FALSE;
      if(do_post) { 
	if((status = posteriors_sspostscript(go, errbuf, ps, msa, mask, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "-s")) { /* make a new postscript page for each sequence in the alignment */
	if((status = individual_seqs_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
      }
      if(! esl_opt_IsDefault(go, "--dfile")) { 
	if((status = drawfile2sspostscript(go, errbuf, ps)) != eslOK) esl_fatal(errbuf);
      }
      if(esl_opt_GetBoolean(go, "--mask-col")) { 
	if(ps->clen != masklen) esl_fatal("MSA number: %d has consensus (non-gap RF) length of %d which != lane mask length of %d.", nali, clen, masklen);
	if((status = colormask_sspostscript(go, errbuf, ps, msa, mask, hc_onecell, BLACKOC, CYANOC)) != eslOK) esl_fatal(errbuf);
      }
      if(! esl_opt_IsDefault(go, "--mask-diff")) { 
	if((status = diffmask_sspostscript(go, errbuf, ps, msa, mask, mask2, hc_onecell, BLACKOC, CYANOC, MAGENTAOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }

      if((status = draw_sspostscript(ofp, go, errbuf, command, date, hc_scheme, ps)) != eslOK) esl_fatal(errbuf);
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

  ps->npage    = 0;
  ps->modelname = NULL;
  ps->titleA = NULL;
  ps->ntitle = 0;
  ps->titlex = 0.;
  ps->titley = 0.;
  ps->legx = 0.;
  ps->legy = 0.;
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
  ps->sclAA       = NULL;
  ps->maskAA      = NULL;
  ps->nalloc      = 50;
  ps->msa_ct      = NULL;
  ps->msa_nbp     = 0;
  return ps;

 ERROR: esl_fatal("create_sspostscript(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: setup_sspostscript()
 * 
 * Purpose:  Determine positions for title and legend in a SSPostscript_t()
 * Return:   eslOK
 */
int
setup_sspostscript(SSPostscript_t *ps, char *errbuf)
{
  int status;
  float pagex;
  float pagey;

  if(ps->clen == 0) ESL_FAIL(eslEINVAL, errbuf, "Failed to ready any residues in template file.");

  /* set up legx, legy, this is a hack (takes advantage of position of 3' residue in all SSU models) */
  ps->legx = ps->rxA[ps->clen-1] + LEGX_OFFSET;
  ps->legy = ps->ryA[ps->clen-1] + LEGY_OFFSET;

  pagex = POSTSCRIPT_PAGEWIDTH / ps->scale;
  pagey = POSTSCRIPT_PAGEHEIGHT / ps->scale;

  ps->titlex = pagex/2.;
  ps->titley = pagey - PAGE_TOPBUF - TITLE_DEFAULT_FONTSIZE;

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

  if(ps->titleA != NULL) {
    for(i = 0; i < ps->ntitle; i++) { 
      if(ps->titleA[i] != NULL) {
	free(ps->titleA[i]);
      }
    }
    free(ps->titleA);
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
  if(ps->nocclA != NULL) free(ps->nocclA);
  if(ps->msa_ct != NULL) free(ps->msa_ct);

  free(ps);
  return;
}

/* Function: create_onecell_colorlegend()
 * 
 * Purpose:  Create and initialize a one cell color legend data structure.
 * Return:   occl
 */
OneCellColorLegend_t *
create_onecell_colorlegend(float *col, float boxsize, char *text)
{
  int status;
  OneCellColorLegend_t *occl;

  ESL_ALLOC(occl, sizeof(OneCellColorLegend_t));

  /* initialize */
  esl_vec_FSet(occl->col, NCMYK, 0.);
  occl->text = NULL;
  
  /* set caller specified values */
  esl_vec_FCopy(col, NCMYK, occl->col);
  occl->boxsize = boxsize;
  if((status = esl_strdup(text, -1, &(occl->text))) != eslOK) esl_fatal("create_onecell_colorlegend(), error copying text");

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
create_scheme_colorlegend(int scheme, int nbins, float boxsize, char *text, float *limits)
{
  int status;
  SchemeColorLegend_t *scl;
  int i;
  ESL_ALLOC(scl, sizeof(SchemeColorLegend_t));

  /* initialize */
  scl->text = NULL;
  
  /* set caller specified values */
  scl->scheme = scheme;
  scl->nbins = nbins;
  ESL_ALLOC(scl->limits, sizeof(float) * (nbins+1));
  for(i = 0; i <= nbins; i++) { scl->limits[i] = limits[i]; }
  scl->boxsize = boxsize;
  if(text != NULL) if((status = esl_strdup(text, -1, &(scl->text))) != eslOK) esl_fatal("create_scheme_colorlegend(), error copying text");

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
  float textsize;

  /* object is valid, print it */
  fprintf(fp, "%sone cell legstart\n", "%");
  x = ps->legx;
  y = ps->legy - (LEG_BOXSIZE * 1.5 * (occl_idx)) ;
  textsize = 16;

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
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", textsize);
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", occl->text, x, (y + (occl->boxsize * .25)));
    /* reset font size to 8 */
    fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  
  /* reset font size to 8 */
  fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
  return eslOK;
}


/* Function: draw_scheme_colorlegend()
 * 
 * Purpose:  Print a scheme color legend to an open file.
 * Return:   eslOK
 */
int 
draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page, int do_mask)
{
  float x, y;
  int cp;
  int c;
  float textsize;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  do_border = (!esl_opt_GetBoolean(go, "-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  /* object is valid, print it */
  fprintf(fp, "%sone cell legstart\n", "%");
  x = ps->legx;
  y = ps->legy - (ps->nocclA[page] * (LEG_BOXSIZE * 1.5));
  textsize = 16;

  /* print text for this legend */
  if(scl->text != NULL) { 
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, textsize);
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", scl->text, x, (y + (scl->boxsize * .25)));
  }
  y -= scl->boxsize;

  fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, textsize);
  /* print masked scheme color cells */
  if(do_mask) { 
    fprintf(fp, "%.1f setlinewidth\n", scl->boxsize/4.);
    for(c = 0; c < scl->nbins; c++) { 
      draw_masked_block(fp, x, y, hc_scheme[c], do_circle_mask, do_square_mask, do_x_mask, do_border, scl->boxsize);
      y -= scl->boxsize;
    }
    y += (scl->boxsize * scl->nbins);
    x += 1.5 * scl->boxsize;
    fprintf(fp, "1.0 setlinewidth\n");
  }
	
  /* print scheme color cells and labels next to them */
  for(c = 0; c < scl->nbins; c++) { 
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", scl->boxsize, scl->boxsize, (-1 * scl->boxsize));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", hc_scheme[c][cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += scl->boxsize * 1.5;
    y += scl->boxsize * 0.25;
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    if(c == scl->nbins-1) fprintf(fp, "(\\[%.3f-%.3f\\]) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);
    else                  fprintf(fp, "(\\[%.3f-%.3f\\)) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);
    x -= scl->boxsize * 1.5;
    y -= scl->boxsize * 0.25;
    y -= scl->boxsize;
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);

  float colvec[NCMYK];
  colvec[0] = colvec[1] = colvec[2] = 0.;
  colvec[3] = 1.0;
  if(do_mask) { /* print cells showing difference between masked and unmasked */
    x -= scl->boxsize * 1.5;
    y -= scl->boxsize;
    fprintf(fp, "%.1f setlinewidth\n", scl->boxsize/4.);
    draw_masked_block(fp, x, y, colvec, do_circle_mask, do_square_mask, do_x_mask, do_border, scl->boxsize);

    /* print label */
    x += scl->boxsize * 1.5;
    y += scl->boxsize * 0.25;
    fprintf(fp, "(positions excluded by mask (all colors)) %.4f %.4f moveto show\n", x, y);
    x -= scl->boxsize * 1.5;
    y -= scl->boxsize * 0.25;
   
    y -= scl->boxsize * 1.5;
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", scl->boxsize, scl->boxsize, (-1 * scl->boxsize));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", colvec[cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    x += scl->boxsize * 1.5;
    y += scl->boxsize * 0.25;
    fprintf(fp, "(positions included by mask (all colors)) %.4f %.4f moveto show\n", x, y);
    x -= scl->boxsize * 1.5;
    y -= scl->boxsize * 0.25;
  }
  
  /* reset font size to 8 */
  fprintf(fp, "/Helvetica findfont 8.00 scalefont setfont\n");
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
  int p, i, c, l;
  float title_fontsize;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  do_border = (!esl_opt_GetBoolean(go, "-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  if(ps->npage == 0) ESL_FAIL(eslEINCOMPAT, errbuf, "draw_sspostscript, ps->npage == 0\n");

  for(p = 0; p < ps->npage; p++) { 
    /* scale section */
    fprintf(fp, "%% begin scale\n");
    fprintf(fp, "%.2f %.2f scale\n", ps->scale, ps->scale);
    fprintf(fp, "%% end scale\n\n");
      
    /* title section */
    fprintf(fp, "%% begin ignore\n");
    fprintf(fp, "/%s findfont %.2f scalefont setfont\n", DEFAULT_FONT, TITLE_DEFAULT_FONTSIZE);
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
    fprintf(fp, "(%s: %d residues; %d basepairs) %.2f %.2f moveto show\n", ps->modelname, ps->clen, ps->msa_nbp, ps->titlex, ps->titley);
    fprintf(fp, "%% end ignore\n");

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
	fprintf(fp, "/Helvetica findfont %.2f scalefont setfont\n", HUNDREDS_DEFAULT_FONTSIZE);
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
	fprintf(fp, "%.2f setlinewidth\n", TICKS_DEFAULT_LINEWIDTH);
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
	fprintf(fp, "%.2f setlinewidth\n", BP_DEFAULT_LINEWIDTH);
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
    fprintf(fp, "/Helvetica findfont %.2f scalefont setfont\n", RESIDUES_DEFAULT_FONTSIZE);
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
    for(i = 0; i < ps->clen; i++) { 
      fprintf(fp, "() %.2f %.2f moveto show\n", ps->rxA[i], ps->ryA[i]);
    }
    fprintf(fp, "%% end text residues\n");

    /* the rest of the text will be ignored by esl-ssudraw if the output
     * file we're creating is read in as a template file later on
     */
    fprintf(fp, "%% begin ignore\n");
    /* print one cell color legends, if any */
    if(ps->occlAAA != NULL && ps->occlAAA[p] != NULL) { 
      for(l = 0; l < ps->nocclA[p]; l++) 
	draw_onecell_colorlegend(fp, ps->occlAAA[p][l], ps, l);
    }
    /* print scheme color legends, if any */
    if(ps->sclAA != NULL && ps->sclAA[p] != NULL) { 
      draw_scheme_colorlegend(go, fp, ps->sclAA[p], hc_scheme[ps->sclAA[p]->scheme], ps, p, ((ps->maskAA[p] == NULL) ? FALSE : TRUE));
    }

    if(ps->rcolAAA != NULL && ps->rcolAAA[p] != NULL) { 
      if(ps->maskAA[p] != NULL) { 
	fprintf(fp, "2.0 setlinewidth\n");
	if(do_border && do_x_mask)      { fprintf(fp, "1.0 setlinewidth\n"); }
	if(do_border && do_square_mask) { fprintf(fp, "2.0 setlinewidth\n"); }
	if(do_border && do_circle_mask) { fprintf(fp, "2.5 setlinewidth\n"); }
	for(c = 0; c < ps->clen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  if(ps->maskAA[p][c] == '0') { 
	    draw_masked_block(fp, ps->rxA[c]-1., ps->ryA[c]-1., ps->rcolAAA[p][c], do_circle_mask, do_square_mask, do_x_mask, do_border, SS_BOXSIZE);
	  }
	  else { /* cell is within mask, ps->maskAA[p][c] == '1' */
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
      for(c = 0; c < ps->clen; c++) { 
	fprintf(fp, "(%c) %.2f %.2f moveto show\n", ps->rrAA[p][c], ps->rxA[c], ps->ryA[c]);
      }
    }
    fprintf(fp, "grestore\nshowpage\n");
    fprintf(fp, "%% end ignore\n\n");
  }
  return eslOK;
}

/* parse_template_file
 *
 * Read a postscript template file derived from the 
 * Gutell CRW website. The file is read in sections.
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
parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  SSPostscript_t *ps = NULL;
  /* Create the postscript object */
  ps = create_sspostscript();

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in parse_template_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  while ((status = esl_fileparser_GetToken(efp, &tok, &toklen))  == eslOK) { 
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
	      if((status = parse_ignore_section(efp, errbuf)) != eslOK) return status;
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
	      ESL_FAIL(eslEINVAL, errbuf, "parse_template_file(), error, unknown section type %s.", tok);
	    }
	  }
	  else { 
	    ESL_FAIL(eslEINVAL, errbuf, "parse_template_file(), error last read line number %d.", efp->linenumber);
	  }
	}
	else { 
	  ESL_FAIL(eslEINVAL, errbuf, "parse_template_file(), expected line beginning with %% begin, but read tok: %s instead of begin, last read line number %d.", tok, efp->linenumber);
	}
      }
      else { 
	ESL_FAIL(eslEINVAL, errbuf, "parse_template_file(), ran out of tokens early, error last read line number %d.", efp->linenumber);
      }
    }
    else { 
      ESL_FAIL(eslEINVAL, errbuf, "parse_template_file(), expected line beginning with %%, read tok: %s, last read line number %d.", tok, efp->linenumber);
    }
  }
  if(status != eslEOF) { 
    ESL_FAIL(status, errbuf, "parse_template_file(), error, ran out of tokens, but not at end of file?, last read line number %d.", efp->linenumber);
  }
  esl_fileparser_Close(efp);

  /* validate the file we just read */
  if((status = validate_justread_sspostscript(ps, errbuf)) != eslOK) return status;
  
  *ret_ps = ps;
  return eslOK;
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
  char *newstr;
  int   curlen = 0;

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
    if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
    curlen += toklen;
    if((status = esl_strcat(&curstr, curlen, " ",  1))      != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
    curlen += 1;
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
parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf)
{
  int status;
  char *tok;
  int   toklen;
  int   keep_reading = TRUE;
  while((keep_reading) && (status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
    if (strcmp(tok, "%") == 0) { 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "ignore") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      keep_reading = FALSE;
      status = eslOK;
    }
  }
  if(status == eslEOF) ESL_FAIL(status, errbuf, "Error, parsing ignore section, finished file looking for '%% end ignore' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing ignore section, last line number read %d", efp->linenumber);

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
  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    curlen = 0;
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
	if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (1).");
	curlen += toklen;
	if((status = esl_strcat(&curstr, curlen, " ",  1))      != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (2).");
	curlen += 1;
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
infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int status;
  int p, pp, c, i;
  int cpos, apos;
  int orig_npage = ps->npage;
  double **obs, *ent, *bg;
  int zero_obs;

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->clen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
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
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);

  int nonecell = 0;
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
    }
    else { 
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent[cpos], ps->sclAA[pp])) != eslOK) return status;
    }

    ps->rrAA[pp][cpos] = ' ';
    /*ps->rrAA[pp][cpos] = (esl_FCompare(ent[cpos], 0., eslSMALLX1) == eslOK) ? '-' : ' ';*/
  }

  /* add one-cell color legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero residues (all gaps) (1000/1000)")+1));
  sprintf(text, "positions with zero residues (all gaps) (%4d/%4d)", nonecell, ps->clen);
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], LEG_BOXSIZE, text);
  ps->nocclA[pp] = 1;
  free(text);

  /* add text to legend */
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
delete_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
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
  float *limits;
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);

  int nonecell = 0;
  if(do_all) { 
    /* draw delete page with all deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if(dct[cpos] == 0) { 
	if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status; 
	nonecell++;
      }
      else {
	if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct[cpos]) / (float) (msa->nseq), ps->sclAA[pp])) != eslOK) return status;
      }
    }
  }
  else { /* do_all is FALSE, draw delete page with only internal deletes */
    for(cpos = 0; cpos < ps->clen; cpos++) { 
      ps->rrAA[pp][cpos] = ' ';
      if(dct_internal[cpos] == 0) { 
	if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status; 
	nonecell++;
      }
      else {
	if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], (float) (dct_internal[cpos]) / (float) (msa->nseq), ps->sclAA[pp])) != eslOK) return status;
      }
    }
  }


  /* add one-cell color legend */
  char *text;
  if(do_all) { 
    ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero deletions (1000/1000)")+1));
    sprintf(text, "positions with zero deletions (%4d/%4d)", nonecell, ps->clen);
  }
  else { 
    ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero internal deletions (1000/1000)")+1));
    sprintf(text, "positions with zero internal deletions (%4d/%4d)", nonecell, ps->clen);
  }
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], LEG_BOXSIZE, text);
  ps->nocclA[pp] = 1;
  free(text);

  /* add color legend */
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

  if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);

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
insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
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

  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->clen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->clen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
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
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);

  for(cpos = 1; cpos <= ps->clen; cpos++) { 
    if(nseq_ict[cpos] == 0) { 
      res = '-';
      col = 0.0;
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][(cpos-1)], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell++;
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
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][(cpos-1)], NCMYK, hc_scheme[hc_scheme_idx], col, ps->sclAA[pp])) != eslOK) return status;
    }
    /*ps->rrAA[pp][(cpos-1)] = res;*/
    ps->rrAA[pp][(cpos-1)] = ' ';
  }

  /* add one-cell color legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero inserts (1000/1000)")+1));
  sprintf(text, "positions with zero inserts (%4d/%4d)", nonecell, ps->clen);
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], LEG_BOXSIZE, text);
  ps->nocclA[pp] = 1;
  free(text);

  /* add color legend */
  ESL_ALLOC(text, sizeof(char) * (strlen("fraction of sequences with inserts:")));
  sprintf(text, "fraction of sequences with inserts:");
  add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  free(text);

  if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);
  
  for(i = 0; i < ps->clen; i++) free(ict[i]);
  free(ict);
  free(total_ict);
  free(nseq_ict);
  free(med_ict);

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
posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
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
  char *text;    /* text for color legends */
  int nonecell_avg = 0;
  int nonecell_seq = 0;

  if(msa->rf == NULL) esl_fatal("No RF annotation in alignment");

  if(esl_opt_GetBoolean(go, "--p-avg"))  { do_avg = TRUE;  new_npage += 1; navg_page = orig_npage; }
  if(esl_opt_GetBoolean(go, "--p-indi")) { do_indi = TRUE; nfirst_indi_page = orig_npage + new_npage; new_npage += msa->nseq; }

  if((status = addpages_sspostscript(ps, new_npage)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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
    nonecell_seq = 0;
    nongap_s = nongaprf_s = 0;
    sum_s    = sumrf_s = 0.;
    if(do_indi) { /* add color legend for this sequence */
      ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);
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
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
	  nonecell_seq++;
	  ps->rrAA[pp][cpos] = ' ';
	}
      }
    } /* done with this sequence */
    if(do_indi) { 
      avg_s   =  (float) sum_s / (float) nongap_s;
      avgrf_s =  (float) sumrf_s / (float) nongaprf_s;

      /* add one-cell color legend */
      ESL_ALLOC(text, sizeof(char) * (strlen("gap positions (1000/1000)")+1));
      sprintf(text, "gap positions (%4d/%4d)", nonecell_seq, ps->clen);
      ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], LEG_BOXSIZE, text);
      ps->nocclA[pp] = 1;
      free(text);

      ESL_ALLOC(text, sizeof(char) * (strlen("posterior probability; 1.000 (RF) 1.000 (all)")+1));
      sprintf(text, "posterior probability; %.3f (RF) %.3f (all)", avgrf_s, avg_s);
      add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
      free(text);
      if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);
      pp++;
    }
  } /* done with all sequences */

  if(do_avg) { /* add average colors */
    pp = navg_page;
    ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);
    for(c = 0; c < msa->alen; c++) { 
      if(a2c_map[c] != -1) { /* cons position */
	cpos = a2c_map[c]; 
	if(nongap_c[c] > 0) { 
	  avgrf_c = sum_c[c] /= (float) nongap_c[c];
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], avgrf_c, ps->sclAA[pp])) != eslOK) return status;
	}
	else { 
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
	  nonecell_avg++;
	}
	ps->rrAA[pp][cpos] = ' ';
      }
    }

    /* add one-cell color legend */
    ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero residues (all gaps) (1000/1000)")+1));
    sprintf(text, "positions with zero residues (all gaps) (%4d/%4d)", nonecell_avg, ps->clen);
    ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], LEG_BOXSIZE, text);
    ps->nocclA[pp] = 1;
    free(text);
    
    /* add color legend */
    ESL_ALLOC(text, sizeof(char) * (strlen("avg posterior probability; 1.000 (RF) 1.000 (all)")+1));
    sprintf(text, "avg posterior probability; %.3f (RF) %.3f (all)", 
	    (esl_vec_FSum(sumrf_c, msa->alen) / (float) (esl_vec_ISum(nongaprf_c, msa->alen))),
	    (esl_vec_FSum(sum_c, msa->alen)   / (float) (esl_vec_ISum(nongap_c, msa->alen))));
    add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
    free(text);
    if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);
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
  
  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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
  char *text;
  ESL_ALLOC(text, sizeof(char) * (strlen("columns included within mask (1000 of 1000 (1.000))")+1));
  sprintf(text, "columns included within mask (%4d of %4d (%.3f))", ncols_inside_mask, ps->clen, (float) ncols_inside_mask / (float) ps->clen);
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incmask_idx], LEG_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * (strlen("columns excluded from  mask (1000 of 1000 (1.000))")+1));
  sprintf(text, "columns excluded from  mask (%4d of %4d (%.3f))", ncols_outside_mask, ps->clen, (float) ncols_outside_mask / (float) ps->clen);
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[excmask_idx], LEG_BOXSIZE, text);
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
  char *text;
  ESL_ALLOC(text, sizeof(char) * (strlen("included by both masks (1000 of 1000 (1.000))")+1));
  sprintf(text, "included by both masks (%4d of %4d (%.3f))", ncols_in_both, ps->clen, (float) ncols_in_both / (float) ps->clen);
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incboth_idx], LEG_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * (strlen("included by mask 1 but not mask 2 (1000 of 1000 (1.000))")+1));
  sprintf(text, "included by mask 1 but not mask 2 (%4d of %4d (%.3f))", ncols_in_1_out_2, ps->clen, (float) ncols_in_1_out_2 / (float) ps->clen);
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[inc1_idx], LEG_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * (strlen("included by mask 2 but not mask 1 (1000 of 1000 (1.000))")+1));
  sprintf(text, "included by mask 2 but not mask 1 (%4d of %4d (%.3f))", ncols_out_1_in_2, ps->clen, (float) ncols_out_1_in_2 / (float) ps->clen);
  ps->occlAAA[pp][2] = create_onecell_colorlegend(hc_onecell[inc2_idx], LEG_BOXSIZE, text);
  free(text);

  ESL_ALLOC(text, sizeof(char) * (strlen("excluded by both masks (1000 of 1000 (1.000))")+1));
  sprintf(text, "excluded by both masks (%4d of %4d (%.3f))", ncols_out_both, ps->clen, (float) ncols_out_both / (float) ps->clen);
  ps->occlAAA[pp][3] = create_onecell_colorlegend(hc_onecell[excboth_idx], LEG_BOXSIZE, text);
  free(text);

  ps->nocclA[pp] = 4;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "diffmask_sspostscript(): memory allocation error.");
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
    ESL_RALLOC(ps->occlAAA, tmp, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->nocclA,  tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->sclAA,   tmp, sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->maskAA,  tmp, sizeof(char **) * (ps->npage + ntoadd));
  }
  for(p = ps->npage; p < (ps->npage + ntoadd); p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->occlAAA[p] = NULL;
    ps->nocclA[p]  = 0;
    ps->sclAA[p]   = NULL;
    ps->maskAA[p]  = NULL;
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
structural_infocontent_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int ss_idx, int zerores_idx)
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

  if(msa->ss_cons == NULL) ESL_FAIL(status, errbuf, "--struct requires #=GC SS_cons annotation in the alignment.");
  if((status = addpages_sspostscript(ps, 1)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

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
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, LEG_BOXSIZE, NULL, limits);

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
	      if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][rapos])) { /* seq i is not a gap at right half */
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
    if(ent_p[cpos] < (-1. * eslSMALLX1)) nss++;
    if(nres[cpos] == 0) { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[zerores_idx])) != eslOK) return status;
      ent_p[cpos] = 0.; /* impt to set to 0., so esl_vec_DSum(ent_p... call below to calc total struct info is accurate */
      nzerores++;
    }
    else if(ent_p[cpos] < (-1. * eslSMALLX1)) { /* single stranded base, paint grey */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[ss_idx])) != eslOK) return status;
      ent_p[cpos] = 0.; /* impt to set to 0., so esl_vec_DSum(ent_p... call below to calc total struct info is accurate */
    }
    else { 
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_scheme[hc_scheme_idx], ent_p[cpos], ps->sclAA[pp])) != eslOK) return status;
    }
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add text to the one cell legend */
  char *text;
  ESL_ALLOC(text, sizeof(char) * (strlen("single-stranded positions (1000/1000)")+1));
  sprintf(text, "single-stranded positions (%4d/%4d)", nss, ps->clen);
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[ss_idx], LEG_BOXSIZE, text);
  ps->nocclA[pp] = 1;
  free(text);

  /* add text to the second one cell legend */
  ESL_ALLOC(text, sizeof(char) * (strlen("positions with zero residues (all gaps) (1000/1000)")+1));
  sprintf(text, "positions with zero residues (all gaps) (%4d/%4d)", nzerores, ps->clen);
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[zerores_idx], LEG_BOXSIZE, text);
  ps->nocclA[pp] = 2;
  free(text);

  /* add text to the scheme legend */
  ESL_ALLOC(text, (sizeof(char) * strlen("structural info content per basepaired posn (total: 1000.00 bits)")+1));
  sprintf(text,  "structural info content per basepaired posn (total: %.2f bits)", esl_vec_DSum(ent_p, ps->clen) * 2.);
  add_text_to_scheme_colorlegend(ps->sclAA[pp], text);
  free(text);

  if(mask != NULL) add_mask_to_ss_postscript(ps, pp, mask);

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
  /*printf("%.3f %d (%.3f)\n", val, bi, scl->limits[bi+1]);*/
  for(ci = 0; ci < ncolvals; ci++) { vec[ci] = scheme[bi][ci]; }
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
  int status;
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
validate_and_update_sspostscript_given_msa(SSPostscript_t *ps, ESL_MSA *msa, char *errbuf, int msa_idx)
{
  int status;
  int *msa_ct;
  int msa_nbp = 0;
  int *tmp_ct;
  int *c2a_map, *a2c_map, msa_clen;
  int apos, cpos;

  ps->msa_idx = msa_idx;

  /* get the CT array for this msa */
  ESL_ALLOC(tmp_ct, sizeof(int) * (msa->alen+1));
  if (esl_wuss2ct(msa->ss_cons, msa->alen, tmp_ct) != eslOK) ESL_FAIL(status, errbuf, "Problem getting ct from SS_cons, does alignment %d of MSA file have SS_cons annotation?", msa_idx);
  /* map cpos to apos */
  map_cpos_to_apos(msa, &c2a_map, &a2c_map, &msa_clen);
  /* convert tmp_ct which is in alignment coords [1..alen] to consensus coords [0..clen-1]*/
  cpos = 0;
  for(apos = 0; apos < msa->alen; apos++) {
    if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) { /* a consensus position */
      if(tmp_ct[(apos+1)] != 0) msa_nbp++;
      msa_ct[cpos++] = tmp_ct[(apos+1)];
    }
  }
  free(tmp_ct);

  if(ps->msa_ct != NULL) { free(ps->msa_ct); ps->msa_ct = NULL;}
  ps->msa_ct = msa_ct;
  ps->msa_nbp = msa_nbp;

  if(ps->clen != msa_clen) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected consensus length of %d in MSA, but read %d\n", ps->clen, msa_clen);
  if(ps->nbp != 0 && ps->nbp != msa_nbp) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected %d basepairs in MSA's SS_cons, but read %d\n", ps->nbp, msa_nbp);

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "validate_and_update_sspostscript_given_msa(), error status %d, probably out of memory.\n", status);
  return status; 
}
