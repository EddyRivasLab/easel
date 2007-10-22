/* Manipulate a multiple sequence alignment in some useful ways.
 *
 * From easel's esl-alistat which was from squid's alistat (1995)
 * EPN, Fri Aug 10 08:52:30 2007
 * SVN $Id$
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <easel.h>
#include <esl_distance.h>
#include <esl_getopts.h>
#include <esl_sqio.h>
#include <esl_msa.h>
#include <esl_distance.h>
#include <esl_dmatrix.h>
#include <esl_vectorops.h>
#include <esl_stack.h>
#include <esl_tree.h>
#include <esl_wuss.h>

static char banner[] = "manipulate a multiple sequence alignment file";
static char usage[]  = "[options] <msafile>\n\
The <msafile> must be in Stockholm format.";

static int  keep_or_remove_rf_gaps(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int keep_flag, int remove_flag);
static int  write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);
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
static int  plot_inserts(FILE *fp, ESL_MSA *msa, char *errbuf);
static int  get_tree_order(ESL_TREE *T, char *errbuf, int **ret_order);
static int  reorder_msa(ESL_MSA *msa, int *order, char *errbuf);
static int  dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);

static ESL_OPTIONS options[] = {
  /* name          type        default  env   range      togs reqs  incomp                      help                                                       docgroup */
  { "-h",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "help; show brief info on version and usage",                     1 },
  { "-o",          eslARG_OUTFILE,NULL, NULL, NULL,      NULL,NULL, NULL,                       "output the alignment to file <f>, not stdout",                   1 },
  { "-s",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                        "print statistics (esl-alistat behavior)",                 0 },
  { "-i",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, "-g,-k,-r,--morph,--amino", "annotate individual secondary structures by imposing consensus", 1 },
  { "-g",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "add/rewrite #=GC RF markup marking consensus columns",           1 },
  { "--gapthresh", eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"-g", NULL,                       "with -g, fraction of gaps to allow in a consensus column",       1 },
  { "-k",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "keep  only columns w/(possibly post -g) non-gap #=GC RF markup", 1 },
  { "-r",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "remove all columns w/(possibly post -g) non-gap #=GC RF markup", 1 },
  { "-v",          eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL, NULL,                       "be verbose (usually with --morph or --merge)",                   1 },
  { "--merge",     eslARG_STRING,FALSE, NULL, NULL,      NULL,NULL,"--morph,-g,-k,-r",          "merge msa in <msafile> with msa in <f>",                         2 },
  { "--morph",     eslARG_STRING,FALSE, NULL, NULL,      NULL,NULL,"--merge",                   "morph msa in <msafile> to msa in <f>'s gap structure",           2 },
  { "--trim",      eslARG_STRING,FALSE, NULL, NULL,      NULL,NULL,"--merge,--morph",           "trim aligned seqs in <msafile> to subseqs in <f>",               2 },
  { "--iinfo",     eslARG_STRING, NULL, NULL, NULL,      NULL,NULL,"--merge,--morph",           "print info on # of insertions b/t all non-gap RF cols to <f>",   2 },
  { "--iplot",     eslARG_STRING, NULL, NULL, NULL,      NULL,NULL,"--merge,--morph",           "plot heatmap of # of insertions b/t all non-gap RF cols to <f>", 2 },
  { "--tree",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--merge,--morph",           "reorder MSA to tree ordering following single linkage clustering", 2 },
  { "--amino",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--dna,--rna",               "<msafile> contains protein alignments",                          3 },
  { "--dna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--rna",             "<msafile> contains DNA alignments",                              3 },
  { "--rna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--dna",             "<msafile> contains RNA alignments",                              3 },
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
  int           i;		/* counter over seqs               */
  int           rlen;		/* a raw (unaligned) seq length    */
  int           small, large;	/* smallest, largest sequence      */
  uint64_t      nres;		/* total # of residues in msa      */
  double        avgid;		/* average fractional pair id      */
  int    max_comparisons;       /* maximum # comparisons for avg id */
  FILE         *ofp;		/* output file (default is stdout) */
  char          errbuf[eslERRBUFSIZE];
  int           write_ali = FALSE; /* set to TRUE if we should print a new MSA */

  /* --merge and --morph related vars */
  ESL_MSAFILE  *otherafp = NULL;	/* other input alignment file (with --morph) */
  ESL_MSA      *othermsa = NULL;	/* other input alignment      (with --morph) */

  /* --trim related vars */
  ESL_SQFILE   *trimfp = NULL;  /* sequence file with subsequences for --trim */

  /* --iinfo, --iplot related vars */
  FILE *iinfofp = NULL;  /* output file for --iinfo */
  FILE *iplotfp = NULL;  /* output file for --iplot */

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
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\nexpert miscellaneous options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
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
  max_comparisons = 1000;

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
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
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
  if((esl_opt_GetBoolean(go, "-i")) && (abc->type != eslRNA && abc->type != eslDNA))
    esl_fatal("-i option pertains to base pairs and only makes sense with DNA or RNA alphabets.");

  /* optionally, open --morph or --merge msa file for reading, --merge and --morph are incompatible,
   * so we'll never try to do both. */
  if(esl_opt_GetString(go, "--morph") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--morph"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--morph alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --morph alignment %s\n", 
					      esl_opt_GetString(go, "--morph"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
    }
  if(esl_opt_GetString(go, "--merge") != NULL)
    {
      status = esl_msafile_OpenDigital(abc, esl_opt_GetString(go, "--merge"), eslMSAFILE_STOCKHOLM, NULL, &otherafp);
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--merge alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --merge alignment %s\n", 
					      esl_opt_GetString(go, "--merge"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
    }

  /***********************************************
   * Read MSAs one at a time.
   ***********************************************/

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      if(esl_opt_GetBoolean(go, "-s")) /* print stats */
	{
	  nres = 0;
	  small = large = -1;
	  for (i = 0; i < msa->nseq; i++)
	    {
	      rlen  = esl_abc_dsqrlen(msa->abc, msa->ax[i]);
	      nres += rlen;
	      if (small == -1 || rlen < small) small = rlen;
	      if (large == -1 || rlen > large) large = rlen;
	    }
	  
	  esl_dst_XAverageId(abc, msa->ax, msa->nseq, max_comparisons, &avgid);
	  printf("Alignment number:    %d\n",     nali);
	  if (msa->name != NULL)
	    printf("Alignment name:      %s\n",   msa->name); 
	  printf("Format:              %s\n",     esl_msa_DescribeFormat(afp->format));
	  printf("Number of sequences: %d\n",     msa->nseq);
	  printf("Alignment length:    %d\n",     msa->alen);
	  printf("Total # residues:    %llu\n",   (unsigned long long) nres);
	  printf("Smallest:            %d\n",     small);
	  printf("Largest:             %d\n",     large);
	  printf("Average length:      %.1f\n",   (double) nres / (double) msa->nseq);
	  /*printf("Average identity:    %.0f\n",   100.*avgid); */
	  printf("//\n");
	}


      /* read other msa if --morph or --merge enabled */
      if((esl_opt_GetString(go, "--morph") != NULL) || (esl_opt_GetString(go, "--merge") != NULL))
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
	if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--trim"));
	else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
	else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
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
      if(esl_opt_GetBoolean(go, "-g"))
	{
	  if((status = write_rf_gapthresh(go, errbuf, msa)) != eslOK) goto ERROR;
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

      /* impose consensus structure to get individual secondary structures, if nec */
      if(esl_opt_GetBoolean(go, "-i"))
	{
	  if((status = individualize_consensus(go, errbuf, msa) != eslOK)) goto ERROR;
	  write_ali = TRUE;
	}

      /* handle the --tree option, if enabled */
      if(! esl_opt_IsDefault(go, "--tree"))
	{
	  ESL_TREE    *T = NULL;/* the tree, created by Single-Linkage Clustering */
	  ESL_DMATRIX *D = NULL;/* the distance matrix */
	  
	  /* Create distance matrix and infer tree by single linkage clustering */
	  esl_dst_XDiffMx(msa->abc, msa->ax, msa->nseq, &D);
	  esl_tree_SingleLinkage(D, &T);
	  esl_tree_SetTaxaParents(T);
	  /* esl_tree_WriteNewick(stdout, T); */
	  esl_tree_Validate(T, NULL);

	  /* Get new order for seqs in the MSA based on the tree */
	  int *order;
	  if((status = get_tree_order(T, errbuf, &order)) != eslOK) goto ERROR;
	  /*for(i = 0; i < msa->nseq; i++) {
	    printf("order[%3d]: %3d\n", i, order[i]);
	    }*/
	  esl_tree_Destroy(T);
	  esl_dmatrix_Destroy(D);
	  if((status = reorder_msa(msa, order, errbuf)) != eslOK) goto ERROR;
	  write_ali = TRUE;
	  free(order);
	}	  

      /* handle the --iinfo option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--iinfo")) {
	if ((iinfofp = fopen(esl_opt_GetString(go, "--iinfo"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --iinfo output file %s\n", esl_opt_GetString(go, "--iinfo"));
	if((status = dump_insert_info(iinfofp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(iinfofp);
      }

      /* handle the --iplot option, if enabled, do this after all MSA has been manipulated due to other options */
      if(! esl_opt_IsDefault(go, "--iplot")) {
	if ((iplotfp = fopen(esl_opt_GetString(go, "--iplot"), "w")) == NULL) 
	  ESL_FAIL(eslFAIL, errbuf, "Failed to open --iplot output file %s\n", esl_opt_GetString(go, "--iplot"));
	if((status = plot_inserts(iplotfp, msa, errbuf) != eslOK)) goto ERROR;
	fclose(iplotfp);
      }

      /* write out alignment, if nec */
      if(write_ali) {
	status = esl_msa_Write(ofp, msa, eslMSAFILE_STOCKHOLM);
	if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
	else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);
      }
      esl_msa_Destroy(msa);
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

  /* Cleanup, normal return
   */
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
  int status;
  int          *useme;
  int           apos;

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
  esl_msa_ColumnSubset(msa, useme);
  free(useme);
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
write_rf_gapthresh(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa)
{
  int  status;
  int  apos;
  int  gaps;
  int  i;
  double gapthresh;

  if(msa->rf == NULL) ESL_ALLOC(msa->rf, sizeof(char) * msa->alen);
  gapthresh = esl_opt_GetReal(go, "--gapthresh");
  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (gaps = 0, i = 0; i < msa->nseq; i++)
	if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) gaps++;
      msa->rf[(apos-1)] = ((double) gaps / (double) msa->nseq > gapthresh) ? '.' : 'x';
    }
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
  int  status;
  int  apos;
  int  i;
  int *cct;		/* 0..alen-1 base pair partners array for consensus        */
  int *ct;		/* 0..alen-1 base pair partners array for current sequence */
  char *ss;             /* individual secondary structure we've built              */

  if(msa->ss_cons == NULL)                                ESL_FAIL(status, errbuf, "-i requires MSA to have consensus structure annotation.\n");
  if(! (msa->flags & eslMSA_DIGITAL))                     ESL_FAIL(status, errbuf, "individualize_consensus() MSA is not digitized.\n");
    
  ESL_ALLOC(cct, sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ct,  sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ss,  sizeof(char) * (msa->alen+1));

  if (esl_wuss2ct(msa->ss_cons, msa->alen, cct) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");
  
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
      if (esl_ct2wuss(ct, msa->alen, ss) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string had pseudoknots, we can't handle this yet.");
      esl_msa_AppendGR(msa, "SS", i, ss);
    }
  free(cct);
  free(ct);
  free(ss);
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

      if(esl_opt_GetBoolean(go, "-v")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* we don't have to add any columns to either msa (okay if 0) */
	{
	  /* do nothing */
	  if(esl_opt_GetBoolean(go, "-v")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns to msa1 */
	{ 
	  if(esl_opt_GetBoolean(go, "-v")) printf("\tmsa1 add     %4d all gap columns\n", (ngaps2-ngaps1)); 
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
	      if(esl_opt_GetBoolean(go, "-v")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd1) - tmp_ngaps), radd);
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
	  if(esl_opt_GetBoolean(go, "-v")) printf("\tmsa2 add     %4d all gap columns\n", (ngaps1 - ngaps2));
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
	      if(esl_opt_GetBoolean(go, "-v")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
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
  if(esl_opt_GetBoolean(go, "-v")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      nadd1 += aadd1[apos1];
      if(esl_opt_GetBoolean(go, "-v")) { printf("%5d %5d\n", apos1, aadd1[apos1]); }
    }
  nadd1 += aadd1[0];
  if(esl_opt_GetBoolean(go, "-v")) printf("Adding  %d columns to msa 1\n", nadd1);

  nadd2 = 0;
  if(esl_opt_GetBoolean(go, "-v")) { printf("Printing number of all gap columns to add after each msa2 alignment column:\n"); }
  for(apos2 = 1; apos2 <= msa2->alen; apos2++)
    {
      nadd2 += aadd2[apos2];
      if(esl_opt_GetBoolean(go, "-v")) { printf("%5d %5d\n", apos2, aadd2[apos2]); }
    }
  nadd2 += aadd2[0];
  if(esl_opt_GetBoolean(go, "-v")) printf("Adding  %d columns to msa 2\n", nadd2);

  /* add the 100% gap columns to msa1 and msa2 */
  status = add_gap_columns_to_msa(errbuf, msa1, aadd1, &new_msa1, TRUE);
  status = add_gap_columns_to_msa(errbuf, msa2, aadd2, &new_msa2, TRUE);

  /* Make new_c2a_map1 and new_c2a_map2, they should be identical */
  if((status = map_cpos_to_apos(new_msa1, &new_c2a_map1, &new_clen1))  != eslOK) goto ERROR;
  if((status = map_cpos_to_apos(new_msa2, &new_c2a_map2, &new_clen2))  != eslOK) goto ERROR;
  if(new_clen1 != new_clen2) 
    ESL_XFAIL(eslEINVAL, errbuf, "Coding error, during alignment merge, after adding gaps, MSA lengths differ.");

  if(esl_opt_GetBoolean(go, "-v")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(new_c2a_map1[cpos] != new_c2a_map2[cpos]) 
	esl_fatal("Coding error. Alignments to merge do not have same consensus position map\n");
      if(esl_opt_GetBoolean(go, "-v")) printf("%4d %4d %4d\n", cpos, new_c2a_map1[cpos], new_c2a_map2[cpos]);
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

      /* new_msa1->name,desc,acc,au untouched (is this unwise?) */

      if(new_msa1->sqlen  != NULL) new_msa1->sqlen[i]  = new_msa2->sqlen[ip];

      if(new_msa1->sslen != NULL) new_msa1->sslen[i]  = new_msa2->sslen[ip];
      if(new_msa1->salen != NULL) new_msa1->salen[i]  = new_msa2->salen[ip];
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

      if(esl_opt_GetBoolean(go, "-v")) printf("%4d: ", cpos); 
      if(ngaps1 == ngaps2) /* keep all columns in between (okay if 0) */
	{
	  for(apos1 = cur_apos1+1; apos1 < nxt_apos1; apos1++) akeep[apos1] = TRUE; 
	  if(esl_opt_GetBoolean(go, "-v")) printf("\n");      
	}
      else if(ngaps1 <  ngaps2) /* we need to add some new 100% gap columns */
	{ 
	  if(esl_opt_GetBoolean(go, "-v")) printf("\tadd     %4d all gap columns\n", (ngaps2-ngaps1)); 
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
	      if(esl_opt_GetBoolean(go, "-v")) printf("\t\tresidues added: %d (%d)\n", ((msa2->nseq * nadd) - tmp_ngaps), radd);
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
	  if(esl_opt_GetBoolean(go, "-v")) printf("\tdelete  %4d/%4d    columns\n", (ngaps1 - ngaps2), (ngaps1));  
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
  if(esl_opt_GetBoolean(go, "-v")) { printf("Printing number of all gap columns to add after each msa1 alignment column:\n"); }
  for(apos1 = 1; apos1 <= msa1->alen; apos1++)
    {
      if(akeep[apos1]) nkeep++;
      else delete += (msa1->nseq - agaps1[apos1]);
      nadd += aadd[apos1];
      if(esl_opt_GetBoolean(go, "-v")) { printf("%5d %5d\n", apos1, aadd[apos1]); }
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

  if(esl_opt_GetBoolean(go, "-v")) printf("printing final test\n\n");
  for(cpos = 1; cpos <= clen; cpos++) 
    {
      if(c2a_map2[cpos] != new_c2a_map1[cpos]) 
	esl_fatal("Coding error. Morphed alignment does not have same consensus position map as %s\n", esl_opt_GetString(go, "--morph"));
      if(esl_opt_GetBoolean(go, "-v")) printf("%4d %4d %4d %4d\n", cpos, c2a_map2[cpos], new_c2a_map1[cpos], (c2a_map2[cpos] - new_c2a_map1[cpos]));
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
    if (status == eslEFORMAT)
      esl_fatal("\
Sequence file parse error, line %d of file %s:\n\
%s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);
    else if (status != eslEOF)
      esl_fatal("Sequence file %s read failed with error code %d\n",
		sqfp->filename, status);
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
      if(! (sq[i]->flags & eslSQ_DIGITAL))  ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq's must be digitized.");
      if(sq[i]->n == 0) ESL_XFAIL(eslEINVAL, errbuf, "in trim_msa(), sq[%d] is zero-length\n", i);

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
      esl_sq_Dealign(uaseq, uaseq, "-_.", msa->alen);
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
  int i;
  int clen;
  int nseq;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in dump_insert_info(), msa must be digitized.");
  if(msa->rf == NULL) ESL_XFAIL(eslEINVAL, errbuf, "No #=GC RF markup in alignment, it is needed for --iplot.");

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
	  if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) ict[cpos][i]++;
    }
  clen = cpos;
  for(cpos = 0; cpos <= clen; cpos++)
    {
      nseq = 0;
      for(i = 0; i < msa->nseq; i++)
	if(ict[cpos][i] >= 1) nseq++;
      if(nseq > 0) 
	printf("%5d %5d\n", cpos, nseq);
    }

  for(i = 0; i <= msa->alen; i++)
    free(ict[i]);
  free(ict);
      
  return eslOK;

 ERROR:
  return status;
}

/* plot_inserts
 *                   
 * Given an MSA with RF annotation, print a postscript heatmap of how
 * many insertions are afer each non-gap RF column (consensus column)
 * in each sequence.
 */
static int plot_inserts(FILE *fp, ESL_MSA *msa, char *errbuf)
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

  if(! (msa->flags & eslMSA_DIGITAL))
    ESL_XFAIL(eslEINVAL, errbuf, "in dump_insert_info(), msa must be digitized.");

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
  dmx_Visualize(fp, I, esl_dmx_Min(I), esl_dmx_Max(I));
  esl_dmatrix_Destroy(I);
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
    for(i = 0; i < msa->nseq; i++) msa->sa[i] = tmp[order[i]];
  }

  /* swap sa, if they exist */
  if(msa->sa != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sa[i];
    for(i = 0; i < msa->nseq; i++) msa->sa[i] = tmp[order[i]];
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
   boxsize = ESL_MIN( (792 - bottommargin) / D->n, 
		      (612 - leftmargin)   / D->m);
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

