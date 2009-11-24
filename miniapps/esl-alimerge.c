/* Merge alignments into a single alignment based on their reference (RF) annotation.
 * 
 * EPN, Fri Nov 20 16:28:59 2009
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_distance.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_vectorops.h"

static char banner[] = "merge alignments based on their reference (RF) annotation";
static char usage1[]  = "[options] <alignment file 1> <alignment file 2>";
static char usage2[]  = "[options] --list <file listing n > 1 ali files to merge>\n\
\n\
  Input alignments must be in Stockholm or Pfam format.\n\
  Ouput format choices\n\
  --------------------\n\
  stockholm [default]\n\
  pfam\n\
  a2m\n\
  psiblast\n\
  afa";

static void read_list_file(char *listfile, char ***ret_alifile_list, int *ret_nalifile);
static void update_maxinsert(ESL_MSA *msa, int clen, int *maxinsert);
static int  validate_and_copy_msa_annotation(const ESL_GETOPTS *go, int outfmt, ESL_MSA *mmsa, ESL_MSA **msaA, int clen, int nmsa, int alen_merged, int *maxinsert, char *errbuf);
static int  add_msa(ESL_MSA *mmsa, ESL_MSA *msa_to_add, int *maxinsert, int clen, int alen_merged, char *errbuf);
static int  gapize_string(char *src_str, int64_t src_len, int64_t dst_len, int *ngapA, char gapchar, char **ret_dst_str);
static int  validate_no_nongaps_in_rf_gaps(const ESL_ALPHABET *abc, char *rf_str, char *other_str, int64_t len);
static int  determine_gap_columns_to_add(ESL_MSA *msa, int *maxinsert, int clen, int **ret_ngapA, char *errbuf);
 
static ESL_OPTIONS options[] = {
  /* name         type          default  env   range togs reqs  incomp           help                                                             docgroup */
  { "--list",     eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,            "command-line argument is a file that lists ali files to merge",  99 },
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",                     1 },
  { "-o",         eslARG_OUTFILE,  NULL, NULL, NULL, NULL,NULL, NULL,            "output the final alignment to file <f>, not stdout",             1 },
  { "-v",         eslARG_NONE,    FALSE, NULL, NULL, NULL,"-o", NULL,            "print info on merge to stdout; requires -o",                     1 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL,NULL, NULL,            "NOT YET DISPLAYED",                                              99 },
  { "--outformat",eslARG_STRING,  FALSE, NULL, NULL, NULL,NULL, NULL,            "specify that output aln be format <s> (see choices above)",      1 },
  { "--rna",      eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "alignments to merge are RNA alignments",                         1 },
  { "--dna",      eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "alignments to merge are DNA alignments",                         1 },
  { "--amino",    eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "alignments to merge are protein alignments",                     1 },
  { "--stall",    eslARG_NONE,    FALSE, NULL, NULL, NULL,NULL, NULL,            "arrest after start: for debugging under gdb",                    99 },  
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  int           status;		               /* easel return code               */
  ESL_GETOPTS  *go      = NULL;	               /* application configuration       */
  ESL_ALPHABET *abc     = NULL;      	       /* biological alphabet             */
  char         *msafile1 = NULL;	       /* msa file 1 (stays NULL if --list) */
  char         *msafile2 = NULL;	       /* msa file 2 (stays NULL if --list) */
  char         *listfile = NULL;	       /* list file name (stays NULL unless --list) */
  int           infmt   = eslMSAFILE_UNKNOWN;  /* format code for input alifiles  */
  int           outfmt  = eslMSAFILE_UNKNOWN;  /* format code for output ali      */
  ESL_MSAFILE  *afp     = NULL;	               /* open alignment file             */
  FILE         *ofp;		               /* output file (default is stdout) */
  char        **alifile_list = NULL;           /* list of alignment files to merge */
  int           nalifile;                      /* size of alifile_list             */
  int           do_stall;                      /* used to stall when debugging     */
  int           abctype;
  int           fi;                            /* counter over alignment files */
  int           ai;                            /* counter over alignments */
  int           nali_cur;                      /* number of sequences in this alignment */
  int           nali_tot;                      /* number of sequences in this alignment */
  int           nseq_tot;                      /* number of sequences in all alignments */
  ESL_MSA     **msaA = NULL;                   /* [0..nali_tot-1] all msas read from all files */
  int          *maxinsert = NULL;              /* [0..cpos..rflen+1] max number of inserts 
						* before each consensus position in all alignments */
  int           nalloc = 0;                    /* current size of msaA */
  int           chunksize = 10;                /* size to increase nalloc by when realloc'ing */
  void         *tmp;                           /* for ESL_RALLOC() */
  int           clen;                          /* consensus length (non-gap #=GC RF length) of all alignments */
  int           cur_clen;                      /* consensus length (non-gap #=GC RF length) of current alignments */
  int           apos;                          /* alignment position */
  ESL_MSA      *mmsa = NULL;                   /* the merged alignment created by merging all alignments in msaA */
  int           alen_mmsa;                     /* number of columns in merged MSA */
  char          errbuf[eslERRBUFSIZE];         /* buffer for error messages */
  char         *tmpstr;                        /* used if -v, for printing file names */

  /* output formatting, only relevant if -v */
  char      *namedashes = NULL;                /* string of dashes, an underline */
  int        ni;                               /* counter                        */
  int        namewidth;                        /* max width of file name         */
  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage1);
      esl_usage (stdout, argv[0], usage2);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      exit(0);
    }

  if(((! esl_opt_GetBoolean(go, "--list")) && (esl_opt_ArgNumber(go) != 2)) ||
     ((  esl_opt_GetBoolean(go, "--list")) && (esl_opt_ArgNumber(go) != 1))) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage1);
      esl_usage(stdout, argv[0], usage2);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if(esl_opt_GetBoolean(go, "--list")) { 
    listfile = esl_opt_GetArg(go, 1);
  }
  else { 
    msafile1 = esl_opt_GetArg(go, 1);
    msafile2 = esl_opt_GetArg(go, 2);
  }

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      esl_fatal("Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else ofp = stdout;

  if (esl_opt_IsOn(go, "--informat")) {
    infmt = esl_msa_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (infmt == eslMSAFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat", esl_opt_GetString(go, "--informat")); 
  }

  if (esl_opt_IsOn(go, "--outformat")) {
    outfmt = esl_msa_EncodeFormat(esl_opt_GetString(go, "--outformat"));
    if (outfmt == eslMSAFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --outformat", esl_opt_GetString(go, "--outformat")); 
  }
  else outfmt = eslMSAFILE_STOCKHOLM;

  do_stall = esl_opt_GetBoolean(go, "--stall"); /* a stall point for attaching gdb */
  while (do_stall); 

  /* determine file names to merge */
  if(listfile != NULL) { /* read list file */
    read_list_file(listfile, &alifile_list, &nalifile);
    if(nalifile == 0) esl_fatal("Failed to read a single alignment file name from %s\n", listfile);
  }
  else { /* we're merging two alignment files from command-line */
    nalifile = 2;
    ESL_ALLOC(alifile_list, sizeof(char *) * nalifile);
    if((status = esl_strdup(msafile1, -1, &(alifile_list[0]))) != eslOK) esl_fatal("Error storing alignment file name %s, error status: %d\n", msafile1, status);
    if((status = esl_strdup(msafile2, -1, &(alifile_list[1]))) != eslOK) esl_fatal("Error storing alignment file name %s, error status: %d\n", msafile2, status);
  }

  if(esl_opt_GetBoolean(go, "-v")) { 
      /* determine the longest file name in alifile_list */
      namewidth = 9; /* length of 'file name' */
      for(fi = 0; fi < nalifile; fi++) { 
	if((status = esl_FileTail(alifile_list[fi], FALSE, &tmpstr)) != eslOK) esl_fatal("Memory allocation error.");
	namewidth = ESL_MAX(namewidth, strlen(tmpstr));
	free(tmpstr);
      }

      ESL_ALLOC(namedashes, sizeof(char) * (namewidth+1));
      namedashes[namewidth] = '\0';
      for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';
      fprintf(stdout, "# Reading %d alignment files...\n", nalifile);
      fprintf(stdout, "#\n");
      fprintf(stdout, "# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "",        namewidth,"",          "",        "",          "",            "",               "ncols");
      fprintf(stdout, "# %7s  %-*s  %7s  %9s  %9s  %13s  %8s\n", "file #",  namewidth,"file name", "ali #",   "#seq/ali",  "ncols/ali",   "# seq total",    "required");
      fprintf(stdout, "# %7s  %*s  %7s  %9s  %9s  %13s  %8s\n", "-------", namewidth, namedashes,  "-------", "---------", "---------",   "-------------", "--------");
    }

  /* Allocate and initialize */
  nalloc = chunksize;
  ESL_ALLOC(msaA, sizeof(ESL_MSA *) * nalloc);

  /****************************************************************
   *  Read alignments one at a time, storing them all, separately *
   ****************************************************************/

  nali_tot = 0;
  nseq_tot = 0;
  for(fi = 0; fi < nalifile; fi++) { 
    status = esl_msafile_Open(alifile_list[fi], eslMSAFILE_UNKNOWN, NULL, &afp);
    if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile_list[fi]);
    else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile_list[fi]);
    else if (status != eslOK)        esl_fatal("Alignment file %s open failed with error %d\n", alifile_list[fi], status);

    if(abc == NULL) { /* this will only be true of first alignment of first file */
      if      (esl_opt_GetBoolean(go, "--amino"))   abc = esl_alphabet_Create(eslAMINO);
      else if (esl_opt_GetBoolean(go, "--dna"))     abc = esl_alphabet_Create(eslDNA);
      else if (esl_opt_GetBoolean(go, "--rna"))     abc = esl_alphabet_Create(eslRNA);
      else {
	status = esl_msafile_GuessAlphabet(afp, &abctype);
	if (status == eslEAMBIGUOUS)    esl_fatal("Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", alifile_list[fi]);
	else if (status == eslEFORMAT)  esl_fatal("Alignment file parse failed: %s\n", afp->errbuf);
	else if (status == eslENODATA)  esl_fatal("Alignment file %s is empty\n", alifile_list[fi]);
	else if (status != eslOK)       esl_fatal("Failed to read alignment file %s\n", alifile_list[fi]);
	abc = esl_alphabet_Create(abctype);
      }
    }
    nali_cur = 0;

    while((status = esl_msa_Read(afp, &(msaA[nali_tot]))) == eslOK) { 
      nali_cur++;
      nali_tot++;
      nseq_tot  += msaA[(nali_tot-1)]->nseq;

      if(nali_tot == nalloc) { /* we need more msa ptrs, reallocate */
	nalloc += chunksize; 
	ESL_RALLOC(msaA, tmp, sizeof(ESL_MSA *) * nalloc+chunksize); 
	for(ai = nali_tot; ai < nalloc; ai++) msaA[ai] = NULL;
      }

      if(msaA[(nali_tot-1)]->rf == NULL) esl_fatal("Error, all alignments must have #=GC RF annotation; alignment %d of file %d does not (%s)\n", nali_cur, (fi+1), alifile_list[fi]); 
      msaA[(nali_tot-1)]->abc = abc;

      /* get current consensus (non-gap RF) length) */
      cur_clen = 0;
      for(apos = 0; apos < (int) msaA[(nali_tot-1)]->alen; apos++) { 
	if(! esl_abc_CIsGap(msaA[(nali_tot-1)]->abc, msaA[(nali_tot-1)]->rf[apos])) cur_clen++;
      }
      if(nali_tot == 1) { /* first alignment, store clen, allocate maxinsert */
	clen = cur_clen;
	ESL_ALLOC(maxinsert, sizeof(int) * (clen+1)); 
	esl_vec_ISet(maxinsert, (clen+1), 0);
      }
      else if(cur_clen != clen) { 
	esl_fatal("Error, all alignments must have identical non-gap #=GC RF lengths; expected (RF length of first ali read): %d,\nalignment %d of file %d length is %d (%s))\n", clen, nali_cur, (fi+1), cur_clen, alifile_list[fi]); 
      }
      update_maxinsert(msaA[(nali_tot-1)], clen, maxinsert);

      if(esl_opt_GetBoolean(go, "-v")) { 
	if((status = esl_FileTail(alifile_list[fi], FALSE, &tmpstr)) != eslOK) esl_fatal("Memory allocation error.");
	fprintf(stdout, "  %7d  %-*s  %7d  %9d  %9" PRId64 "  %13d  %8d\n", (fi+1),  namewidth, tmpstr, nali_tot, msaA[(nali_tot-1)]->nseq, msaA[(nali_tot-1)]->alen, nseq_tot, (clen+esl_vec_ISum(maxinsert, (clen+1))));
	free(tmpstr);
      }
    } /* end of while esl_msa_Read() loop */

    if      (status == eslEFORMAT) esl_fatal("Alignment file %s, parse error:\n%s\n", alifile_list[fi], afp->errbuf);
    else if (status == eslEINVAL)  esl_fatal("Alignment file %s, parse error:\n%s\n", alifile_list[fi], afp->errbuf);
    else if (status != eslEOF)     esl_fatal("Alignment file %s, read failed with error code %d\n", alifile_list[fi], status);
    if(nali_cur == 0)              esl_fatal("Failed to read any alignments from file %s\n", alifile_list[fi]);
    esl_msafile_Close(afp);
  } /* end of for (fi=0; fi < nalifile; fi++) */

  /*********************************************
   *  Merge all alignments into the merged MSA *
   *********************************************/

  /* We allocate space for all sequences, but left sequences as NULL (-1), this saves us from 
   * needing the full amount of memory now; we'll copy aligned sequences from each alignment,
   * freeing them in the orignal msaA alignments as we go to reduce amount of memory we require
   * to do the merge.
   */     
  mmsa = esl_msa_Create(nseq_tot, -1); 
  alen_mmsa = clen + esl_vec_ISum(maxinsert, (clen+1)); 

  /* Determine what annotation from the input alignments 
   * we will include in the merged MSA.
   * See comments in header of validate_and_copy_msa_annotation()
   * for rules on what we include.
   */
  if((status = validate_and_copy_msa_annotation(go, outfmt, mmsa, msaA, nali_tot, clen, alen_mmsa, maxinsert, errbuf)) != eslOK)
    esl_fatal("Error while checking and copying individual MSA annotation to merged MSA:%s\n", errbuf);
  
  if(esl_opt_GetBoolean(go, "-v")) { 
    fprintf(stdout, "#\n");
    fprintf(stdout, "# Merging alignments ... ");
    fflush(stdout);
  }

  for(ai = 0; ai < nali_tot; ai++) { 
    if((status = add_msa(mmsa, msaA[ai], maxinsert, clen, alen_mmsa, errbuf)) != eslOK) esl_fatal("Error, merging alignment %d of %d:\n%s.", (ai+1), nali_tot, errbuf); 
    esl_msa_Destroy(msaA[ai]); /* note: the aligned sequences will have already been freed in add_msa() */
    msaA[ai] = NULL;
  }
  mmsa->alen = alen_mmsa; /* it was -1, b/c we filled in each seq as we marched through each msaA[] alignment */
  if(esl_opt_GetBoolean(go, "-v")) { 
    fprintf(stdout, "done.\n");
    fprintf(stdout, "#\n");
  }

  if(ofp != stdout) {
    fprintf(stdout, "# Saving alignment to file %s ... ", esl_opt_GetString(go, "-o"));
    fflush(stdout);
  }

  status = esl_msa_Write(ofp, mmsa, outfmt);
  if(status != eslOK) esl_fatal("Error, during alignment output; status code: %d\n", status);

  if(ofp != stdout) { 
    fprintf(stdout, "done.\n");
    fclose(ofp);
  }

  /* clean up and exit */
  if(alifile_list != NULL) { 
    for(fi = 0; fi < nalifile; fi++) { 
      if(alifile_list[fi] != NULL) free(alifile_list[fi]); 
    }
    free(alifile_list);
  }

  if(namedashes != NULL) free(namedashes);
  if(msaA != NULL)       free(msaA);
  if(maxinsert != NULL)  free(maxinsert);
  if(mmsa != NULL)       esl_msa_Destroy(mmsa);
  if(abc  != NULL)       esl_alphabet_Destroy(abc);
    
  esl_getopts_Destroy(go);
  return 0;

 ERROR: 
  esl_fatal("Out of memory.");
  return eslEMEM; /*NEVERREACHED*/
}

/* Function: read_list_file
 * Date:     EPN, Fri Nov 20 16:41:32 2009
 * 
 * Read a file listing alignment files to merge.
 * Store file names in *ret_alifile_list and return it,
 * return number of files in ret_nalifile and return it.
 * Each white-space delimited token is considered a 
 * different alignment name. 
 * 
 * Dies if we encounter an error.
 * 
 * Returns: void.
 */
void
read_list_file(char *listfile, char ***ret_alifile_list, int *ret_nalifile)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int nalloc     = 10;
  int chunksize  = 10;
  char **alifile_list = NULL;
  int n = 0;
  void *tmp;

  ESL_ALLOC(alifile_list, sizeof(char *) * nalloc);
  status = esl_fileparser_Open(listfile, NULL,  &efp);
  if     (status == eslENOTFOUND) esl_fatal("List file %s does not exist or is not readable\n", listfile);
  else if(status == eslEMEM)      esl_fatal("Ran out of memory when opening list file %s\n", listfile);
  else if(status != eslOK)        esl_fatal("Error opening list file %s\n", listfile);
  
  while((status = esl_fileparser_GetToken(efp, &tok, NULL)) != eslEOF) {
    if(n == nalloc) { nalloc += chunksize; ESL_RALLOC(alifile_list, tmp, sizeof(char *) * nalloc); }
    if((status = esl_strdup(tok, -1, &(alifile_list[n++]))) != eslOK) {
      esl_fatal("Error storing alignment file name while reading list file %s, error status: %d\n", listfile, status);
    }
  }
  esl_fileparser_Close(efp);
  *ret_alifile_list = alifile_list;
  *ret_nalifile = n;

  return;

 ERROR:
  esl_fatal("Out of memory.");
  return; /*NOTREACHED*/
}

/* Function: update_maxinsert
 * Date:     EPN, Sun Nov 22 09:40:48 2009
 * 
 * Update maxinsert[] an array that keeps track of 
 * the max number of inserted (gap #=GC RF) columns
 * before each cpos (consensus (non-gap #=GC RF) column).
 *
 * Consensus columns are index [0..cpos..(clen-1)].
 * 
 * maxinsert[0]      is number of inserts before 1st cpos.
 * maxinsert[clen-1] is number of inserts before final cpos.
 * maxinsert[clen]   is number of inserts after  final cpos.
 * 
 * Caller has already checked that msa->rf != NULL
 * and its non-gap length is clen. If we find either
 * of these is not true, we die (but this shouldn't happen).
 * 
 * Returns: void.
 */
void
update_maxinsert(ESL_MSA *msa, int clen, int *maxinsert) 
{
  int apos;
  int cpos = 0;
  int nins = 0;

  for(apos = 0; apos < msa->alen; apos++) { 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
    }
    else { 
      maxinsert[cpos] = ESL_MAX(maxinsert[cpos], nins);
      cpos++;
      nins = 0;
    }
  }
      
  /* update final value, maxinsert[clen+1], the number of inserts
   * after the final consensus position */
  maxinsert[cpos] = ESL_MAX(maxinsert[cpos], nins);
  if(cpos != clen) esl_fatal("Unexpected error in update_maxinsert(), expected clen (%d) not equal to actual clen (%d).\n", clen, cpos);

  return;
}

/* Function: validate_and_copy_msa_annotation
 * Date:     EPN, Tue Nov 24 05:35:50 2009
 * 
 * Decide what individual MSA annotation from
 * the input alignments in msaA[], if any, will be 
 * included in the merged alignment (mmsa) and 
 * add that info to it.
 * 
 * Rules for what to include in mmsa:
 *
 * msaA[]->name, msaA[]->desc, msaA[]->acc annotation is not 
 * included in the merged alignment since, presumably, they 
 * should be specific to the individual alignments.
 *
 * We include author annotation in merged alignment if it is 
 * identical in all msaA[] input alignments.
 *
 * We include comments and per-file (GF) annotation if they 
 * are present and identical in all input msaA[] alignments.
 *
 * We include per-column (GC) annotation if it is present and
 * identical *with-respect-to* #=GC RF annotation AND all 
 * the annotation in gap  #=GC RF columns in all msaA[] are gaps 
 * ('.'). This also pertains to the following parsed per-column 
 * annotation: ss_cons, sa_cons, pp_cons, and rf. With rf,
 * de-gapped rf annotation must be identical in all input 
 * alignments period, if it is not we'll die with an error message. 
 *
 * Per-sequence information and per-residue information is always
 * included in merged alignment. This is done by add_msa() function.
 * 
 * Returns: eslOK on success.
 *          eslEINCONCEIVABLE if input/msa is corrupt in some way (example: ngf>0 but gf_tag[0] == NULL)
 *          eslEMEM on memory error
 *          if !eslOK, errbuf is filled before return
 */
int
validate_and_copy_msa_annotation(const ESL_GETOPTS *go, int outfmt, ESL_MSA *mmsa, ESL_MSA **msaA, int nmsa, int clen, int alen_merged, int *maxinsert, char *errbuf)
{
  int status;
  int *ngapA = NULL;
  int j;                   /* counter over alignment annotations */
  int j2;                  /* counter over alignment annotations */
  int ai;                  /* counter over alignments */
  char *dealigned  = NULL; /* a temporary, dealigned string */
  char *dealigned2 = NULL; /* another temporary, dealigned string */
  char *gapped_out = NULL; /* a temporary string with gaps added to fit into merged aln */
  int do_add;
  int found_tag;
  int be_verbose = FALSE;

  /* we only print info about annotation if -v AND we'll actually
   * output it (as stockholm or pfam) 
   * (we actually don't even need to be in this function if we're
   * not output in stockholm or pfam...)
   */
  if((esl_opt_GetBoolean(go, "-v")) && 
     (outfmt == eslMSAFILE_STOCKHOLM || outfmt == eslMSAFILE_PFAM)) 
    { be_verbose = TRUE; }

  if(be_verbose) fprintf(stdout, "#\n");

  if(nmsa == 0) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "in validate_and_copy_msa_annotation(): zero child alignments.");

  /* First, determine how many all gap columns to insert after each alignment position
   * of the first child msa, so we can (possibly) gap out GC,SS_cons,SA_cons,PP_cons annotation 
   * to appropriate length when adding it to the merged MSA. */
  if((status = determine_gap_columns_to_add(msaA[0], maxinsert, clen, &(ngapA), errbuf)) != eslOK) 
    return status;

  /* Note: esl_strcmp() can handle NULL strings (they are not identical to non-NULL strings) */

  /*********************************************************************/
  /* Check if author annotation is identical in all alignments */
  do_add = TRUE; /* until proven otherwise */
  if(msaA[0]->au != NULL) { 
    for(ai = 1; ai < nmsa; ai++) { 
      if(esl_strcmp(msaA[0]->au, msaA[ai]->au) != 0) { do_add = FALSE; break; }
    }
    if(do_add) { 
      if(be_verbose) fprintf(stdout, "# Identical author annotation from all alignments transferred to merged alignment.\n"); 
      if((status = esl_strdup(msaA[0]->au, -1, &(mmsa->au))) != eslOK) goto ERROR;
    }
    else if(be_verbose) fprintf(stdout, "# Author annotation is not identical in all alignments; not included in merged alignment.\n"); 
  }
  else if(be_verbose) fprintf(stdout, "# Author annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check per-file (GF) annotation, must be present and identical in all msaA[] alignments to be included */
  if(msaA[0]->ngf > 0) { 
    for(j = 0; j < msaA[0]->ngf; j++) { 
      do_add = TRUE; /* until proven otherwise */
      /* verify that what we think is true is true */
      if(msaA[0]->gf_tag[j] == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpectedly, GF tag %d of msaA[0] is NULL, but msaA[0]->ngf is %d.\n", j, msaA[0]->ngf);
      if(msaA[0]->gf[j]    == NULL)  ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpectedly, GF annotation %d of msaA[0] is NULL, but msaA[0]->ngf is %d.\n", j, msaA[0]->ngf);
      for(ai = 1; ai < nmsa; ai++) { 
	found_tag = FALSE;
	for(j2 = 0; j2 < msaA[ai]->ngf; j2++) { 
	  if(esl_strcmp(msaA[0]->gf_tag[j], msaA[ai]->gf_tag[j2]) == 0) {
	    found_tag = TRUE;
	    if(esl_strcmp(msaA[0]->gf[j], msaA[ai]->gf_tag[j2]) != 0) { 
	      do_add = FALSE; 
	    }
	    break; /* if we found a match, do_add remains TRUE */
	  }
	}
	if(found_tag && do_add) { 
	  if(be_verbose) fprintf(stdout, "# Identical GF tag %s annotation from all alignments transferred to merged alignment.\n", msaA[0]->gf_tag[j]);
	  if((status = esl_msa_AddGF(mmsa, msaA[0]->gf_tag[j], msaA[0]->gf[j])) != eslOK) goto ERROR;
	}
	else { 
	  if(be_verbose) fprintf(stdout, "# GF tag %s annotation from first alignment absent from >= 1 other alignments; not included in merged alignment.\n", msaA[0]->gf_tag[j]);
	}
      }
    }
  }
  else if(be_verbose) fprintf(stdout, "# Unparsed GF annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check comments, all must be identically ordered and identical in all msaA[] aligments to include them */
  if(msaA[0]->ncomment > 0) { 
    do_add = TRUE; /* until proven otherwise */
    /* make sure all alignments have same number of comments */
    for(ai = 1; ai < nmsa; ai++) { 
      if(msaA[ai]->ncomment != msaA[0]->ncomment) { 
	do_add = FALSE;
	break;
      }
    }
    if(do_add) { 
      /* make sure all alignments have identical comments */
      for(j = 0; j < msaA[0]->ncomment; j++) { 
	/* verify that what we think is true is true */
	if(msaA[0]->comment[j] == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpectedly, comment %d of msaA[0] is NULL, but msaA[0]->ncomment is %d.\n", j, msaA[0]->ncomment);
	for(ai = 1; ai < nmsa; ai++) { 
	  if(esl_strcmp(msaA[0]->comment[j], msaA[ai]->comment[j]) != 0) { /* comment doesn't match */
	    do_add = FALSE;
	    break;
	  }
	}
      }
    }
    if(do_add) { 
      for(j = 0; j < msaA[0]->ncomment; j++) { 
	if((status = esl_msa_AddComment(mmsa, msaA[0]->comment[j]))!= eslOK) goto ERROR;
      }
      if(be_verbose) fprintf(stdout, "# All alignments have identical comments in the same order. These were transferred to merged alignment.\n"); 
    }
    else { 
      if(be_verbose) fprintf(stdout, "# Comments are not identical in all alignments; not included in merged alignment.\n"); 
    }
  }
  else if(be_verbose) fprintf(stdout, "# No comments in (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check unparsed per-column (GC) annotation, it must include all gaps in gap RF columns and 
   * be identical once gap RF columns are removed in all msaA[] alignments to be included. */

  if(msaA[0]->ngc > 0) { 
    for(j = 0; j < msaA[0]->ngc; j++) { 
      do_add = TRUE; /* until proven otherwise */
      /* verify that what we think is true is true */
      if(msaA[0]->gc_tag[j] == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpectedly, GC tag %d of msaA[0] is NULL, but msaA[0]->ngf is %d.\n", j, msaA[0]->ngc);
      if(msaA[0]->gc[j]     == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpectedly, GC annotation %d of msaA[0] is NULL, but msaA[0]->ngf is %d.\n", j, msaA[0]->ngc);

      /* ensure it does not have non-gaps in gap RF columns */
      if(validate_no_nongaps_in_rf_gaps(msaA[0]->abc, msaA[0]->rf, msaA[0]->gc[j], msaA[0]->alen)) { /* returns TRUE if gc[j] has 0 non-gap characters in gap columns of RF annotation */
	/* dealign gc line */
	if((status = esl_strdup(msaA[0]->gc[j], msaA[0]->alen, &(dealigned))) != eslOK) goto ERROR;
	if((status = esl_strdealign(dealigned, msaA[0]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning GC tag %s of msaA[0]", msaA[0]->gc_tag[j]);

	for(ai = 1; ai < nmsa; ai++) { 
	  found_tag = FALSE;
	  for(j2 = 0; j2 < msaA[ai]->ngc; j2++) { 
	    if(esl_strcmp(msaA[0]->gc_tag[j], msaA[ai]->gc_tag[j2]) == 0) {
	      found_tag = TRUE;
	      /* ensure it does not have non-gaps in gap RF columns */
	      if(validate_no_nongaps_in_rf_gaps(msaA[ai]->abc, msaA[ai]->rf, msaA[ai]->gc[j2], msaA[ai]->alen)) { /* returns TRUE if gc[j2] has 0 non-gap characters in gap columns of RF annotation */
		/* dealign */
		if((status = esl_strdup(msaA[ai]->gc[j2], msaA[ai]->alen, &(dealigned2))) != eslOK) goto ERROR;
		if((status = esl_strdealign(dealigned2, msaA[ai]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning GC tag %s of msaA[%d]", msaA[ai]->gc_tag[j2], ai);
		/* check identity */
		if(esl_strcmp(dealigned, dealigned2) != 0) { do_add = FALSE; }
		free(dealigned2); 
		dealigned2 = NULL; 
		break; /* if we matched, do_add remains TRUE */
	      }
	    }
	  }
	} /* end of (for(ai = 1...)) */
	if(dealigned != NULL) { free(dealigned); dealigned = NULL; }
	if(found_tag && do_add) { 
	  /* gap out the the GC annotation to fit in merged alignment */
	  if((status = gapize_string(msaA[0]->gc[j], msaA[0]->alen, alen_merged, ngapA, '.', &(gapped_out))) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding gaps to create GC tag %s annotation for merged alignment.", msaA[0]->gc_tag[j]);
	  if((status = esl_msa_AppendGC(mmsa, msaA[0]->gc_tag[j], gapped_out)) != eslOK) goto ERROR;
	  free(gapped_out);
	  gapped_out = NULL;
	  if(be_verbose) fprintf(stdout, "# Identical GC tag %s annotation from all alignments transferred to merged alignment.\n", msaA[0]->gc_tag[j]); 
	}
	else { 
	  if(be_verbose) fprintf(stdout, "# GC tag %s annotation from first alignment absent from or different in >= 1 other alignments; not included in merged alignment.\n", msaA[0]->gc_tag[j]);
	}
      } 
    } /* end of for(j = 0 j < msaA[0]->ngc... */
  } /* end of if(msaA[0]->ngc > 0) */
  else if(be_verbose) fprintf(stdout, "# Unparsed GC annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check ss_cons: it must include all gaps in gap RF columns and be identical once gap RF columns are removed in all 
   * msaA[] alignments to be included. (Same requirements as unparsed GC annotation, so code block below is analogous to one above). */
  if(msaA[0]->ss_cons != NULL) { 
    do_add = TRUE; /* until proven otherwise */
    /* ensure it does not have non-gaps in gap RF columns */
    if(validate_no_nongaps_in_rf_gaps(msaA[0]->abc, msaA[0]->rf, msaA[0]->ss_cons, msaA[0]->alen)) { /* returns TRUE if ss_cons has 0 non-gap characters in gap columns of RF annotation */
      /* dealign ss_cons */
      if((status = esl_strdup(msaA[0]->ss_cons, msaA[0]->alen, &(dealigned))) != eslOK) goto ERROR;
      if((status = esl_strdealign(dealigned, msaA[0]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning ss_cons of msaA[0]");
      for(ai = 1; ai < nmsa; ai++) { 
	if(msaA[ai]->ss_cons == NULL) { 
	  do_add = FALSE; 
	  break;
	}
	/* ss_cons != NULL, ensure it does not have non-gaps in gap RF columns */
	if(validate_no_nongaps_in_rf_gaps(msaA[ai]->abc, msaA[ai]->rf, msaA[ai]->ss_cons, msaA[ai]->alen)) { /* returns TRUE if ss_cons has 0 non-gap characters in gap columns of RF annotation */
	  /* dealign */
	  if((status = esl_strdup(msaA[ai]->ss_cons, msaA[ai]->alen, &(dealigned2))) != eslOK) goto ERROR;
	  if((status = esl_strdealign(dealigned2, msaA[ai]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning ss_cons of msaA[%d]", ai);
	  /* check identity */
	  if(esl_strcmp(dealigned, dealigned2) != 0) { do_add = FALSE; }
	  free(dealigned2); 
	  dealigned2 = NULL; 
	  break; /* if we matched, do_add remains TRUE */
	}
      } /* end of (for(ai = 1...)) */
      if(dealigned != NULL) { free(dealigned); dealigned = NULL; }
      if(do_add) { 
	/* gap out the the ss_cons to fit in merged alignment */
	if((status = gapize_string(msaA[0]->ss_cons, msaA[0]->alen, alen_merged, ngapA, '.', &(gapped_out))) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding gaps to create SS_cons annotation for merged alignment.");
	if(mmsa->ss_cons != NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding SS_cons to merged alignment, it is already non-NULL.");
	if((status = esl_strdup(gapped_out, alen_merged, &(mmsa->ss_cons))) != eslOK) goto ERROR;
	free(gapped_out);
	gapped_out = NULL;
	if(be_verbose) fprintf(stdout, "# Identical SS_cons annotation from all alignments transferred to merged alignment.\n");
      }
      else { 
	if(be_verbose) fprintf(stdout, "# SS_cons annotation from first alignment absent from or different in >= 1 other alignments; not included in merged alignment.\n");
      }
    }
  } /* end of if(msaA[0]->ss_cons != NULL) */
  else if(be_verbose) fprintf(stdout, "# SS_cons annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check sa_cons: it must include all gaps in gap RF columns and be identical once gap RF columns are removed in all 
   * msaA[] alignments to be included. (Same requirements as unparsed GC annotation, so code block below is analogous to one above). */
  if(msaA[0]->sa_cons != NULL) { 
    do_add = TRUE; /* until proven otherwise */
    /* ensure it does not have non-gaps in gap RF columns */
    if(validate_no_nongaps_in_rf_gaps(msaA[0]->abc, msaA[0]->rf, msaA[0]->sa_cons, msaA[0]->alen)) { /* returns TRUE if sa_cons has 0 non-gap characters in gap columns of RF annotation */
      /* dealign sa_cons */
      if((status = esl_strdup(msaA[0]->sa_cons, msaA[0]->alen, &(dealigned))) != eslOK) goto ERROR;
      if((status = esl_strdealign(dealigned, msaA[0]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning sa_cons of msaA[0]");
      for(ai = 1; ai < nmsa; ai++) { 
	if(msaA[ai]->sa_cons == NULL) { 
	  do_add = FALSE; 
	  break;
	}
	/* sa_cons != NULL, ensure it does not have non-gaps in gap RF columns */
	if(validate_no_nongaps_in_rf_gaps(msaA[ai]->abc, msaA[ai]->rf, msaA[ai]->sa_cons, msaA[ai]->alen)) { /* returns TRUE if sa_cons has 0 non-gap characters in gap columns of RF annotation */
	  /* dealign */
	  if((status = esl_strdup(msaA[ai]->sa_cons, msaA[ai]->alen, &(dealigned2))) != eslOK) goto ERROR;
	  if((status = esl_strdealign(dealigned2, msaA[ai]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning sa_cons of msaA[%d]", ai);
	  /* check identity */
	  if(esl_strcmp(dealigned, dealigned2) != 0) { do_add = FALSE; }
	  free(dealigned2); 
	  dealigned2 = NULL; 
	  break; /* if we matched, do_add remains TRUE */
	}
      } /* end of (for(ai = 1...)) */
      if(dealigned != NULL) { free(dealigned); dealigned = NULL; }
      if(do_add) { 
	/* gap out the the sa_cons to fit in merged alignment */
	if((status = gapize_string(msaA[0]->sa_cons, msaA[0]->alen, alen_merged, ngapA, '.', &(gapped_out))) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding gaps to create SA_cons annotation for merged alignment.");
	if(mmsa->sa_cons != NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding SA_cons to merged alignment, it is already non-NULL.");
	if((status = esl_strdup(gapped_out, alen_merged, &(mmsa->sa_cons))) != eslOK) goto ERROR;
	free(gapped_out);
	gapped_out = NULL;
	if(be_verbose) fprintf(stdout, "# Identical SA_cons annotation from all alignments transferred to merged alignment.\n");
      }
      else { 
	if(be_verbose) fprintf(stdout, "# SA_cons annotation from first alignment absent from or different in >= 1 other alignments; not included in merged alignment.\n");
      }
    }
  } /* end of if(msaA[0]->sa_cons != NULL) */
  else if(be_verbose) fprintf(stdout, "# SA_cons annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Check pp_cons: it must include all gaps in gap RF columns and be identical once gap RF columns are removed in all 
   * msaA[] alignments to be included. (Same requirements as unparsed GC annotation, so code block below is analogous to one above). */
  if(msaA[0]->pp_cons != NULL) { 
    do_add = TRUE; /* until proven otherwise */
    /* ensure it does not have non-gaps in gap RF columns */
    if(validate_no_nongaps_in_rf_gaps(msaA[0]->abc, msaA[0]->rf, msaA[0]->pp_cons, msaA[0]->alen)) { /* returns TRUE if pp_cons has 0 non-gap characters in gap columns of RF annotation */
      /* dealign pp_cons */
      if((status = esl_strdup(msaA[0]->pp_cons, msaA[0]->alen, &(dealigned))) != eslOK) goto ERROR;
      if((status = esl_strdealign(dealigned, msaA[0]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning pp_cons of msaA[0]");
      for(ai = 1; ai < nmsa; ai++) { 
	if(msaA[ai]->pp_cons == NULL) { 
	  do_add = FALSE; 
	  break;
	}
	/* pp_cons != NULL, ensure it does not have non-gaps in gap RF columns */
	if(validate_no_nongaps_in_rf_gaps(msaA[ai]->abc, msaA[ai]->rf, msaA[ai]->pp_cons, msaA[ai]->alen)) { /* returns TRUE if pp_cons has 0 non-gap characters in gap columns of RF annotation */
	  /* dealign */
	  if((status = esl_strdup(msaA[ai]->pp_cons, msaA[ai]->alen, &(dealigned2))) != eslOK) goto ERROR;
	  if((status = esl_strdealign(dealigned2, msaA[ai]->rf, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning pp_cons of msaA[%d]", ai);
	  /* check identity */
	  if(esl_strcmp(dealigned, dealigned2) != 0) { do_add = FALSE; }
	  free(dealigned2); 
	  dealigned2 = NULL; 
	  break; /* if we matched, do_add remains TRUE */
	}
      } /* end of (for(ai = 1...)) */
      if(dealigned != NULL) { free(dealigned); dealigned = NULL; }
      if(do_add) { 
	/* gap out the the pp_cons to fit in merged alignment */
	if((status = gapize_string(msaA[0]->pp_cons, msaA[0]->alen, alen_merged, ngapA, '.', &(gapped_out))) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding gaps to create PP_cons annotation for merged alignment.");
	if(mmsa->pp_cons != NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding PP_cons to merged alignment, it is already non-NULL.");
	if((status = esl_strdup(gapped_out, alen_merged, &(mmsa->pp_cons))) != eslOK) goto ERROR;
	free(gapped_out);
	gapped_out = NULL;
	if(be_verbose) fprintf(stdout, "# Identical PP_cons annotation from all alignments transferred to merged alignment.\n");
      }
      else { 
	if(be_verbose) fprintf(stdout, "# PP_cons annotation from first alignment absent from or different in >= 1 other alignments; not included in merged alignment.\n");
      }
    }
  } /* end of if(msaA[0]->pp_cons != NULL) */
  else if(be_verbose) fprintf(stdout, "# PP_cons annotation absent from (at least) first alignment; not included in merged alignment.\n"); 

  /*********************************************************************/
  /* Finally, validate that RF annotation is identical in all alignments after removing gaps. */

  if(msaA[0]->rf == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "All alignments must have #= GC RF annotation."); 
  /* dealign rf */
  if((status = esl_strdup(msaA[0]->rf, msaA[0]->alen, &(dealigned))) != eslOK) goto ERROR;
  if((status = esl_strdealign(dealigned, dealigned, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning RF of msaA[0]");
  for(ai = 1; ai < nmsa; ai++) { 
    if(msaA[ai]->rf == NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "All alignments must have #= GC RF annotation."); 
    /* dealign */
    if((status = esl_strdup(msaA[ai]->rf, msaA[ai]->alen, &(dealigned2))) != eslOK) goto ERROR;
    if((status = esl_strdealign(dealigned2, dealigned2, "-_.", NULL)) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "unexpected error dealigning RF of msaA[%d]", ai);
    /* check identity */
    if(esl_strcmp(dealigned, dealigned2) != 0) { 
      ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "All alignments must have identical #=GC RF annotation, once gaps (\".\",\"-\",\"_\") are removed.\nAlignment %d de-gapped RF annotation differs from that of alignment 1.\n%s\n%s", ai+1, dealigned, dealigned2);
    }
    if(dealigned2 != NULL) { free(dealigned2); dealigned2 = NULL; }
  }
  if(dealigned  != NULL) { free(dealigned);  dealigned = NULL; }
  /* gap out the the RF to fit in merged alignment */
  if((status = gapize_string(msaA[0]->rf, msaA[0]->alen, alen_merged, ngapA, '.', &(gapped_out))) != eslOK) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding gaps to create RF annotation for merged alignment.");
  if(mmsa->rf != NULL) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "Error adding RF to merged alignment, it is already non-NULL.");
  if((status = esl_strdup(gapped_out, alen_merged, &(mmsa->rf))) != eslOK) goto ERROR;
  free(gapped_out);
  gapped_out = NULL;
  if(be_verbose) fprintf(stdout, "# Identical RF annotation from all alignments transferred to merged alignment.\n");

  if(dealigned  != NULL) free(dealigned);
  if(dealigned2 != NULL) free(dealigned2);
  if(gapped_out != NULL) free(gapped_out);
  if(ngapA != NULL) free(ngapA);
  return eslOK;
  
 ERROR:
  if(dealigned  != NULL) free(dealigned);
  if(dealigned2 != NULL) free(dealigned2);
  if(gapped_out != NULL) free(gapped_out);
  if(ngapA      != NULL) free(ngapA);
  return status;
}

/* Function: add_msa
 * Date:     EPN, Mon Nov 23 05:54:37 2009
 * 
 * Add a "child" MSA we read from a file to the merged
 * MSA - the merged alignment that we'll eventually
 * output. We free each string in the child as soon
 * as we've added it to the merged, to save memory.
 * 
 * We add all sequence data (aseq), and per sequence
 * annotation, including sqname, sqdesc, sqacc, pp, ss,
 * sa, as well as non-parsed GS and GR annotation.
 * 
 * <maxinsert>[0..clen] is an array specifying the 
 * number of inserted columns necessary between
 * each consensus position.
 * 
 * maxinsert[0]      is number of inserts before 1st cpos.
 * maxinsert[clen-1] is number of inserts before final cpos.
 * maxinsert[clen]   is number of inserts after  final cpos.
 * 
 * <alen_merged> is the number of columns in the merged
 * alignment. This is the non-gap RF length plus the
 * sum of the maxinsert vector.
 * 
 * Returns: eslOK on success.
 *          eslEMEM on memory allocation failure.
 */
int
add_msa(ESL_MSA *mmsa, ESL_MSA *msa_to_add, int *maxinsert, int clen, int alen_merged, char *errbuf)
{
  int status;
  int i;              /* counter over sequences in msa_to_add */
  int j;              /* counter over alignment annotations */
  int mi;             /* counter over sequences in mmsa */
  void *tmp;          /* for reallocations */
  char *tmpstr;       /* used for copying GR annotation */
  int nseq_existing;  /* number of sequences already added to mmsa, by previous calls of this function */
  int *ngapA = NULL;  /* number of all gap columns to add after each alignment position to fill
		       * it out to width of final, merged alignment */

  nseq_existing = mmsa->nseq;

  /* determine how many all gap columns to insert after each alignment position
   * of the child msa when copying it to the merged msa */
  if((status = determine_gap_columns_to_add(msa_to_add, maxinsert, clen, &(ngapA), errbuf)) != eslOK) 
    return status;

  nseq_existing = mmsa->nseq; 

  /* Append msa_to_add's sequence data and per-sequence annotation to mmsa after adding necessary gaps */
  /* sequence names and aligned sequence data (digitized) */
  for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      esl_strdup(msa_to_add->sqname[i], -1, &(mmsa->sqname[mi]));

      if((status = gapize_string(msa_to_add->aseq[i], msa_to_add->alen, alen_merged, ngapA, '.', &(mmsa->aseq[mi]))) != eslOK) 
	ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d.\n", i+1);
      free(msa_to_add->aseq[i]); /* free immediately */
      msa_to_add->aseq[i] = NULL;
  }

  /* parsed annotation that is optional */
  /* sqacc */
  if(msa_to_add->sqacc != NULL) { 
    if(mmsa->sqacc == NULL) { /* allocate for all sequences, even ones added in previous calls to add_msa() */
      ESL_ALLOC(mmsa->sqacc, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
      for(mi = 0; mi < nseq_existing; mi++) { mmsa->sqacc[mi] = NULL; }
    }
    else { /* reallocate; to add space for new seqs */
      ESL_RALLOC(mmsa->sqacc, tmp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
    }
    for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      if(msa_to_add->sqacc[i] != NULL) { 
	if((status = esl_strdup(msa_to_add->sqacc[i], -1, &(mmsa->sqacc[mi]))) != eslOK)  
	  ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence accession number %d.\n", i+1);
	free(msa_to_add->sqacc[i]); /* free immediately */
	msa_to_add->sqacc[i] = NULL;
      }
      else { 
	mmsa->sqacc[mi] = NULL; 
      }
    }
  }
  /* sqdesc */
  if(msa_to_add->sqdesc != NULL) { 
    if(mmsa->sqdesc == NULL) { /* allocate for all sequences, even ones added in previous calls to add_msa() */
      ESL_ALLOC(mmsa->sqdesc, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
      for(mi = 0; mi < nseq_existing; mi++) { mmsa->sqdesc[mi] = NULL; }
    }
    else { /* reallocate; to add space for new seqs */
      ESL_RALLOC(mmsa->sqdesc, tmp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
    }
    for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      if(msa_to_add->sqdesc[i] != NULL) { 
	if((status = esl_strdup(msa_to_add->sqdesc[i], -1, &(mmsa->sqdesc[mi]))) != eslOK)
	  ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence description number %d.\n", i+1);
	free(msa_to_add->sqdesc[i]); /* free immediately */
	msa_to_add->sqdesc[i] = NULL;
      }
      else { 
	mmsa->sqdesc[mi] = NULL; 
      }
    }
  }
  /* per-seq posterior probabilities */
  if(msa_to_add->pp != NULL) { 
    if(mmsa->pp == NULL) { /* allocate for all sequences, even ones added in previous calls to add_msa() */
      ESL_ALLOC(mmsa->pp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
      for(mi = 0; mi < nseq_existing; mi++) { mmsa->pp[mi] = NULL; }
    }
    else { /* reallocate; to add space for new seqs */
      ESL_RALLOC(mmsa->pp, tmp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
    }
    for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      if(msa_to_add->pp[i] != NULL) { 
	if((status = gapize_string(msa_to_add->pp[i], msa_to_add->alen, alen_merged, ngapA, '.', &(mmsa->pp[mi]))) != eslOK) 
	  ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d posterior probabilities.\n", i+1);
	free(msa_to_add->pp[i]); /* free immediately */
	msa_to_add->pp[i] = NULL;
      }
      else { 
	mmsa->pp[mi] = NULL; 
      }
    }
  }
  /* per-seq secondary structure */
  if(msa_to_add->ss != NULL) { 
    if(mmsa->ss == NULL) { /* allocate for all sequences, even ones added in previous calls to add_msa() */
      ESL_ALLOC(mmsa->ss, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
      for(mi = 0; mi < nseq_existing; mi++) { mmsa->ss[mi] = NULL; }
    }
    else { /* reallocate; to add space for new seqs */
      ESL_RALLOC(mmsa->ss, tmp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
    }
    for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      if(msa_to_add->ss[i] != NULL) { 
	if((status = gapize_string(msa_to_add->ss[i], msa_to_add->alen, alen_merged, ngapA, '.', &(mmsa->ss[mi]))) != eslOK) 
	  ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d secondary structure.\n", i+1);
	free(msa_to_add->ss[i]); /* free immediately */
	msa_to_add->ss[i] = NULL;
      }
      else { 
	mmsa->ss[mi] = NULL; 
      }
    }
  }
  /* per-seq surface accessibility */
  if(msa_to_add->sa != NULL) { 
    if(mmsa->sa == NULL) { /* allocate for all sequences, even ones added in previous calls to add_msa() */
      ESL_ALLOC(mmsa->sa, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
      for(mi = 0; mi < nseq_existing; mi++) { mmsa->sa[mi] = NULL; }
    }
    else { /* reallocate; to add space for new seqs */
      ESL_RALLOC(mmsa->sa, tmp, sizeof(char *) * (nseq_existing + msa_to_add->nseq));
    }
    for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
      if(msa_to_add->sa[i] != NULL) { 
	if((status = gapize_string(msa_to_add->sa[i], msa_to_add->alen, alen_merged, ngapA, '.', &(mmsa->sa[mi]))) != eslOK) 
	  ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d surface accessibility.\n", i+1);
	free(msa_to_add->sa[i]); /* free immediately */
	msa_to_add->sa[i] = NULL;
      }
      else { 
	mmsa->sa[mi] = NULL; 
      }
    }
  }
  /* Unparsed per-sequence (GS) annotation */
  if(msa_to_add->ngs > 0) { 
    for(j = 0; j < msa_to_add->ngs; j++) { 
      for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) {
	if(msa_to_add->gs[j][i] != NULL) 
	  if((status =esl_msa_AddGS(mmsa, msa_to_add->gs_tag[j], mi, msa_to_add->gs[j][i])) != eslOK) 
	    ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d GS annotation.\n", i+1);
      }
      free(msa_to_add->gs[j][i]); /* free immediately */
      msa_to_add->gs[j][i] = NULL;
    }
  }
  /* caller will free the rest of gs via esl_msa_Destroy() */
  
  /* Unparsed per-residue (GR) annotation */
  if(msa_to_add->gr != NULL) { 
    for(j = 0; j < msa_to_add->ngr; j++) { 
      for(i = 0, mi = nseq_existing; i < msa_to_add->nseq; i++, mi++) { 
	if(msa_to_add->gr[j][i] != NULL) {
	  if((status = gapize_string(msa_to_add->gr[j][i], msa_to_add->alen, alen_merged, ngapA, '.', &(tmpstr))) != eslOK) 
	    ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d GR annotation.\n", i+1);
	  if((status = esl_msa_AppendGR(mmsa, msa_to_add->gr_tag[j], mi, tmpstr)) != eslOK)
	    ESL_XFAIL(status, errbuf, "Memory allocation error when copying sequence number %d GR annotation.\n", i+1);
	  free(tmpstr);
	  free(msa_to_add->gr[j][i]); /* free immediately */
	  msa_to_add->gr[j][i] = NULL;
	}
      }
    }
  }
  /* caller will free the rest of gr in esl_msa_Destroy() */

  /* msa_to_add should be destroyed by caller */
    
  /* update nseq in mmsa */
  mmsa->nseq += msa_to_add->nseq;

  if(ngapA != NULL) free(ngapA);
  return eslOK;
  
 ERROR:
  if(ngapA != NULL) free(ngapA);
  return status;
}

/* gapize_string
 *                   
 * Given a string, create a new one that is a copy of it, 
 * but with gaps added before each position (apos) as specified 
 * by ngapA[0..apos..len]. <gapchar> specifies the gap character
 * to add.
 * 
 * ngapA[0]       - number of gaps to add before first posn
 * ngapA[apos]    - number of gaps to add before posn apos
 * ngapA[src_len] - number of gaps to add after  final posn
 * 
 * ret_str is allocated here.
 *
 * Returns eslOK on success.
 *         eslEMEM on memory error.
 */
int 
gapize_string(char *src_str, int64_t src_len, int64_t dst_len, int *ngapA, char gapchar, char **ret_dst_str)
{
  int status;
  int src_apos = 0;
  int dst_apos  = 0;
  int i;
  char *dst_str;

  ESL_ALLOC(dst_str, sizeof(char) * (dst_len+1));
  dst_str[dst_len] = '\0';

  /* add gaps before first position */
  for(i = 0; i < ngapA[0]; i++) dst_str[dst_apos++] = gapchar;

  /* add gaps after every position */
  for(src_apos = 0; src_apos < src_len; src_apos++) { 
    dst_str[dst_apos++] = src_str[src_apos];
    for(i = 0; i < ngapA[(src_apos+1)]; i++) dst_str[dst_apos++] = gapchar;
  }

  *ret_dst_str = dst_str;
  return eslOK;

 ERROR: 
  return eslEMEM;
}

/* validate_no_nongaps_in_rf_gaps
 *                   
 * Given an RF string with gaps defined as by alphabet <abc>
 * and another string of same length. Make sure none of the 
 * positions that are gaps in the RF string are non-gaps in the
 * other string. Return TRUE if none are. Return FALSE if at
 * least one is.
 * 
 * Returns TRUE if 0 characters in <other_str> in same position
 *         as a gap in <rf_str> are non-gaps. FALSE otherwise.
 */
int 
validate_no_nongaps_in_rf_gaps(const ESL_ALPHABET *abc, char *rf_str, char *other_str, int64_t len) 
{
  int64_t i;
  for(i = 0; i < len; i++) { 
    if((esl_abc_CIsGap(abc, rf_str[i])) && (! esl_abc_CIsGap(abc, other_str[i]))) return FALSE;
  }
  return TRUE;
}


/* determine_gap_columns_to_add
 *                   
 * Given <maxinsert>, an array of the number of gap RF (inserts) 
 * positions between each non-gap RF (consensus) position in the 
 * eventual merged alignment, calculate how many inserts we need to
 * add at each position of <msa> to expand it out to the
 * appropriate size of the eventual merged alignment.
 * 
 * maxinsert[0]      is number of inserts before 1st cpos in merged aln
 * maxinsert[cpos]   is number of inserts after  final cpos in merged aln
 *                   for cpos = 1..clen 
 * clen is the number of non-gap RF positions in msa (and in eventual merged msa).             
 * 
 * We allocate fill and return ret_ngapA[0..msa->alen] here.
 *
 * ret_ngapA[0]      is number of inserts to add before 1st position of msa 
 * ret_ngapA[apos]   is number of inserts to add after alignment position apos
 *                   for apos = 1..msa->alen
 * 
 * Returns eslOK on success.
 *         eslEMEM on memory alloaction error 
 *         eslERANGE if a value exceeds what we expected (based on earlier
 *                   checks before this function was entered).
 *         if !eslOK, errbuf if filled.
 */
int
determine_gap_columns_to_add(ESL_MSA *msa, int *maxinsert, int clen, int **ret_ngapA, char *errbuf)
{
  int status;
  int apos;
  int apos_for_inserts;
  int prv_apos = 0;  /* alignment position corresponding to consensus position cpos-1 */
  int cpos = 0;
  int nins = 0;
  int *ngapA = NULL;

  ESL_ALLOC(ngapA, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(ngapA, (msa->alen+1), 0);
  
  for(apos = 0; apos < msa->alen; apos++) { 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
    }
    else { 
      if(maxinsert[cpos] < nins) { 
	if(ngapA != NULL) free(ngapA);
	ESL_XFAIL(status, errbuf, "%d inserts before cpos %d greater than max expected (%d).\n", nins, cpos, maxinsert[cpos]); 
      }

      if (cpos == 0)  
	apos_for_inserts = prv_apos; /* inserts before first position: flush right (so add all-gap columns after leftmost column) */
      else 
	apos_for_inserts = apos - ((apos - prv_apos) / 2); /* all other positions: split inserts */

      ngapA[apos_for_inserts] = maxinsert[cpos] - nins;
      cpos++;
      prv_apos = apos;
      nins = 0;
    }
  }
  /* deal with final consensus position */
  apos_for_inserts = apos; /* inserts after final position: flush left (so add all-gap columns after rightmost column) */
  ngapA[apos_for_inserts] = maxinsert[cpos] - nins;

  /* validate that clen is what it should be */
  if(cpos != clen) { 
    if(ngapA != NULL) free(ngapA);
    ESL_FAIL(status, errbuf, "consensus length (%d) is not the expected length (%d).", cpos, clen);
  }

  *ret_ngapA = ngapA;

  return eslOK;

 ERROR: 
  if(ngapA != NULL) free(ngapA);
  ESL_FAIL(status, errbuf, "Memory allocation error.");
  return status; /*NEVERREACHED*/
}
