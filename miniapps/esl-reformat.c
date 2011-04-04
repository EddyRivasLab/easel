/* Convert sequence file formats
 *             
 * SRE, Sun Feb 27 08:24:33 2005
 * from squid's sreformat (1993).
 * SVN $Id$            
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_sq.h"
#include "esl_msa.h"
#include "esl_wuss.h"

static char banner[] = "convert sequence file formats";

static char usage[] = "[-options] <format> <seqfile>\n\
  Output format choices: Unaligned      Aligned\n\
                         -----------    -------\n\
                         fasta          stockholm\n\
                                        pfam\n\
                                        a2m\n\
                                        psiblast\n\
                                        afa\n\
\n";

#define INCOMPATWITHSMALLOPT "--mingap,--nogap,--ignore,--acceptx"

static ESL_OPTIONS options[] = {
   /* name          type        default env   range togs  reqs        incompat                     help                                      docgroup */
  { "-d",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-r",                  "convert to DNA alphabet (U->T)",                     0 },
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       NULL,                  "help; print brief info on version and usage",        0 },
  { "-l",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-u",                  "convert to lower case",                              0 },
  { "-n",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-x",                  "remove DNA IUPAC codes; convert ambig chars to N",   0 },
  { "-o",         eslARG_STRING,  NULL, NULL, NULL, NULL, NULL,       NULL,                  "send output to file <f>, not stdout",                0 },
  { "-r",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-d",                  "convert to RNA alphabet (T->U)",                     0 }, 
  { "-u",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-l",                  "convert to upper case",                              0 },
  { "-x",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "-n",                  "convert non-IUPAC chars (e.g. X) in DNA to N",       0 },
  { "--gapsym",   eslARG_STRING,  NULL, NULL, NULL, NULL, NULL,       "--mingap,--nogap",    "convert all gaps to character <c>",                  0 },
  { "--informat", eslARG_STRING,  NULL, NULL, NULL, NULL, NULL,       NULL,                  "input sequence file is in format <s>",               0 },
  { "--mingap",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "--nogap",             "remove columns containing all gaps (seqfile=MSA)",   0 },
  { "--keeprf",   eslARG_NONE,   FALSE, NULL, NULL, NULL, "--mingap", NULL,                  "with --mingap, keep all nongap #=GC RF columns",     0 },
  { "--nogap",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "--mingap,--gapsym",   "remove columns containing any gaps (seqfile=MSA)",   0 },
  { "--wussify",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "--dewuss,--fullwuss", "convert old RNA structure markup lines to WUSS",     0 },
  { "--dewuss",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "--wussify,--fullwuss","convert WUSS RNA structure markup to old format",    0 },
  { "--fullwuss", eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       "--wussify,--dewuss",  "convert simple WUSS notation to full (output) WUSS", 0 },
  { "--ignore",   eslARG_STRING, FALSE, NULL, NULL, NULL, NULL,       NULL,                  "ignore input seq characters listed in string <s>",   0 },
  { "--acceptx",  eslARG_STRING, FALSE, NULL, NULL, NULL, NULL,       NULL,                  "accept input seq chars in string <s> as X",          0 },
  { "--rename",   eslARG_STRING, FALSE, NULL, NULL, NULL, NULL,       NULL,                  "rename and number each sequence <s>.<n>",            0 },
  { "--small",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,       INCOMPATWITHSMALLOPT,  "use minimal RAM, input must be pfam, output must be afa or pfam",0 },
  { 0,0,0,0,0,0,0,0 },
};

static void symconvert(char *s, char *oldsyms, char *newsyms);
static void regurgitate_pfam_as_afa(ESL_MSAFILE *afp, FILE *ofp, char *alifile, char *gapsym, int force_lower, int force_upper, int force_rna, int force_dna, int iupac_to_n, int x_is_bad, char *rename, int *ret_reached_eof);
static int  regurgitate_pfam_as_pfam(ESL_MSAFILE *afp, FILE *ofp, char *gapsym, int force_lower, int force_upper, int force_rna, int force_dna, int iupac_to_n, int x_is_bad, int wussify, int dewuss, int fullwuss);
  
int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;	        /* application configuration               */
  char        *outformat;	/* output format as a string               */
  char        *infile;          /* name of input sequence file             */
  int          infmt;		/* input format as a code; eslSQFILE_FASTA */
  int          outfmt;		/* output format as a code                 */
  int          status;		/* return code from an Easel call          */
  FILE        *ofp;		/* output stream                           */

  char  *informat;		/* input format as string; "fasta"           */
  char  *outfile;		/* output file, or NULL                      */
  int    force_rna;		/* TRUE to force RNA alphabet                */
  int    force_dna;		/* TRUE to force DNA alphabet                */
  int    force_lower;		/* TRUE to force lower case                  */
  int    force_upper;		/* TRUE to force upper case                  */
  int    iupac_to_n;            /* TRUE to convert ambiguities all to N's    */
  int    x_is_bad;		/* TRUE to convert X to N                    */
  int    do_mingap;		/* TRUE to remove cols containing all gaps   */
  int    do_nogap;		/* TRUE to remove cols containing any gaps   */
  char  *gapsym;		/* NULL if unset; else, char for gaps        */
  int    wussify;		/* TRUE to convert old KH SS markup to WUSS  */
  int    dewuss;		/* TRUE to convert WUSS back to old KH       */
  int    fullwuss;		/* TRUE to convert simple WUSS to full WUSS  */
  char  *rename; 		/* if non-NULL rename seqs to <s>.<n>        */
  int    do_small;		/* TRUE to operate in small memory mode      */
  int    reached_eof;           /* reached EOF? used only in small mem mode  */
  int    idx;                   /* counter over sequences                    */
  int    nali;                  /* number of alignments read                 */
  char   errbuf[eslERRBUFSIZE]; /* for error messages                        */
  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "-h"))
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("  where options are:\n");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0= group; 2 = indentation; 80=textwidth*/
      exit(EXIT_SUCCESS);
    }

  force_dna   = esl_opt_GetBoolean(go, "-d");
  force_lower = esl_opt_GetBoolean(go, "-l");
  iupac_to_n  = esl_opt_GetBoolean(go, "-n");
  outfile     = esl_opt_GetString (go, "-o");
  force_rna   = esl_opt_GetBoolean(go, "-r");
  force_upper = esl_opt_GetBoolean(go, "-u");
  x_is_bad    = esl_opt_GetBoolean(go, "-x");
  gapsym      = esl_opt_GetString( go, "--gapsym");
  informat    = esl_opt_GetString( go, "--informat");
  do_mingap   = esl_opt_GetBoolean(go, "--mingap");
  do_nogap    = esl_opt_GetBoolean(go, "--nogap");
  wussify     = esl_opt_GetBoolean(go, "--wussify");
  dewuss      = esl_opt_GetBoolean(go, "--dewuss");
  fullwuss    = esl_opt_GetBoolean(go, "--fullwuss");
  rename      = esl_opt_GetString (go, "--rename");
  do_small    = esl_opt_GetBoolean(go, "--small");

  if (esl_opt_ArgNumber(go) != 2) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  outformat = esl_opt_GetArg(go, 1);
  infile    = esl_opt_GetArg(go, 2);

  infmt = eslSQFILE_UNKNOWN;
  if (informat != NULL)
    {
      infmt = esl_sqio_EncodeFormat(informat);
      if (infmt == eslSQFILE_UNKNOWN)
	esl_fatal("%s is not a recognized input seqfile format\n", informat);
    }

  outfmt = esl_sqio_EncodeFormat(outformat);
  if (outfmt == eslSQFILE_UNKNOWN)
    esl_fatal("%s is not a recognized output seqfile format\n", outformat);

  /* if --small, make sure infmt == pfam and (outfmt == afa || outfmt == pfam) */
  if(do_small && (infmt != eslMSAFILE_PFAM || (outfmt != eslMSAFILE_AFA && outfmt != eslMSAFILE_PFAM)))  
    esl_fatal("--small requires '--informat pfam' and output format of either 'afa' or 'pfam'");

  if (gapsym != NULL && strlen(gapsym) != 1)
    esl_fatal("Argument to --gapsym must be a single character.");
  
  if (outfile == NULL) ofp = stdout;
  else if ((ofp = fopen(outfile, "w")) == NULL)
    esl_fatal("Failed to open output file %s\n", outfile);

  /***********************************************
   * Reformat the file, printing to stdout.
   ***********************************************/

  /* If the output format is an alignment, then the input format
   * has to be an alignment.
   */
  if (esl_sqio_IsAlignment(outfmt))
    {
      ESL_MSAFILE *afp;
      ESL_MSA     *msa;

      status = esl_msafile_Open(infile, infmt, NULL, &afp);
      if (status == eslENOTFOUND)
	esl_fatal("Alignment file %s not readable\n", infile);
      else if (status == eslEFORMAT) 
	esl_fatal("Couldn't determine format of alignment %s\n", infile);
      else if (status == eslEINVAL)
	esl_fatal("Can't autodetect format of stdin or .gz; use --informat\n");
      else if (status != eslOK) 
	esl_fatal("Alignment file open failed with error %d\n", status);

      if ( esl_opt_IsOn(go, "--ignore"))  esl_fatal("The --ignore option is unimplemented for alignment reformatting.");
      if ( esl_opt_IsOn(go, "--acceptx")) esl_fatal("The --acceptx option is unimplemented for alignment reformatting.");

      nali = 0;

      if (do_small) { 
	if(infmt == eslMSAFILE_PFAM && outfmt == eslMSAFILE_AFA) {
	  if(afp->do_stdin) esl_fatal("--small with afa out format and stdin input is unimplemented.");
	  regurgitate_pfam_as_afa(afp, ofp, infile, gapsym, force_lower, force_upper, force_rna, force_dna, iupac_to_n, x_is_bad, rename, &reached_eof);
	  if(! reached_eof) esl_fatal("Input file contains >1 alignments, but afa formatted output file can only contain 1");
	}
	else if (infmt == eslMSAFILE_PFAM && outfmt == eslMSAFILE_PFAM) {
	  if(rename != NULL) esl_fatal("--rename is unimplemented for combination of --small and output format pfam"); 
	  while((status = regurgitate_pfam_as_pfam(afp, ofp, gapsym, force_lower, force_upper, force_rna, force_dna, iupac_to_n, x_is_bad, wussify, dewuss, fullwuss)) != eslEOF) { 
	    if      (status == eslEFORMAT) esl_fatal("--small alignment file parse error:\n%s\n", afp->errbuf);
	    else if (status == eslEINVAL)  esl_fatal("--small alignment file parse error:\n%s\n", afp->errbuf);
	    else if (status != eslOK)      esl_fatal("--small alignment file read failed with error code %d\n", status);
	  }
	  esl_msafile_Close(afp);
	}
	else { /* do_small enabled, but neither (infmt==pfam && outfmt=afa) nor (infmt==pfam && outfmt==pfam) */
	  esl_fatal("--small requires '--informat pfam' and output format of either 'afa' or 'pfam'");
	}
      }
      else { /* normal mode, --small not enabled */
	while ((status = esl_msa_Read(afp, &msa)) != eslEOF)
	  {
	    nali++;
	    if      (status == eslEFORMAT) esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
	    else if (status == eslEINVAL)  esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
	    else if (status != eslOK)      esl_fatal("Alignment file read failed with error code %d\n", status);

	    if (nali > 1 && outfmt == eslMSAFILE_AFA)      esl_fatal("Input file contains >1 alignments, but afa formatted output file can only contain 1");
	    if (nali > 1 && outfmt == eslMSAFILE_A2M)      esl_fatal("Input file contains >1 alignments, but a2m formatted output file can only contain 1");
	    if (nali > 1 && outfmt == eslMSAFILE_PSIBLAST) esl_fatal("Inputt file contains >1 alignments, but psiblast formatted output file can only contain 1");

	    if (do_mingap)    if((status = esl_msa_MinimGaps(msa, errbuf, "-_.~", esl_opt_GetBoolean(go, "--keeprf"))) != eslOK) esl_fatal(errbuf);
	    if (do_nogap)     if((status = esl_msa_NoGaps   (msa, errbuf, "-_.~"))                                     != eslOK) esl_fatal(errbuf);
	    if (gapsym!=NULL) esl_msa_SymConvert(msa, "-_.", gapsym);
	    if (force_lower)  esl_msa_SymConvert(msa,
						 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
						 "abcdefghijklmnopqrstuvwxyz");
	    if (force_upper)  esl_msa_SymConvert(msa,
						 "abcdefghijklmnopqrstuvwxyz",
						 "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	    if (force_rna)    esl_msa_SymConvert(msa, "Tt", "Uu");
	    if (force_dna)    esl_msa_SymConvert(msa, "Uu", "Tt");
	    if (iupac_to_n)   esl_msa_SymConvert(msa, 
						 "RYMKSWHBVDrymkswhbvd",
						 "NNNNNNNNNNnnnnnnnnnn");
	    if (x_is_bad)     esl_msa_SymConvert(msa, "Xx", "Nn");
	  
	    if (rename)
	      {
		for (idx = 0; idx < msa->nseq; idx++)
		  esl_msa_FormatSeqName(msa, idx, "%s.%d", rename, idx+1);
	      }

	    if (wussify)
	      {
		if (msa->ss_cons != NULL) 
		  esl_kh2wuss(msa->ss_cons, msa->ss_cons);
		if (msa->ss != NULL)
		  for (idx = 0; idx < msa->nseq; idx++)
		    if (msa->ss[idx] != NULL)
		      esl_kh2wuss(msa->ss[idx], msa->ss[idx]);
	      }

	    if (dewuss)
	      {
		if (msa->ss_cons != NULL)
		  esl_wuss2kh(msa->ss_cons, msa->ss_cons);
		if (msa->ss != NULL)
		  for (idx = 0; idx < msa->nseq; idx++)
		    if (msa->ss[idx] != NULL)
		      esl_wuss2kh(msa->ss[idx], msa->ss[idx]);
	      }

	    if (fullwuss)
	      {
		if (msa->ss_cons != NULL)
		  {
		    status = esl_wuss_full(msa->ss_cons, msa->ss_cons);
		    if (status == eslESYNTAX) 
		      esl_fatal("Bad consensus SS: not in WUSS format\n");
		    else if (status != eslOK)
		      esl_fatal("Conversion of SS_cons failed, code %d\n", status);
		  }
		if (msa->ss != NULL)
		  for (idx = 0; idx < msa->nseq; idx++)
		    if (msa->ss[idx] != NULL)
		      {
			status = esl_wuss_full(msa->ss[idx], msa->ss[idx]);
			if (status == eslESYNTAX) 
			  esl_fatal("Bad SS for %s: not in WUSS format\n",
				    msa->sqname[idx]);
			else if (status != eslOK)
			  esl_fatal("Conversion of SS for %s failed, code %d\n", 
				    msa->sqname[idx], status);
		      }
	      }

	    esl_msa_Write(ofp, msa, outfmt);
	    esl_msa_Destroy(msa);
	  }
	esl_msafile_Close(afp);
      } 
    } /* end of alignment->alignment conversion */
  else
    { /* else: conversion to unaligned file formats */
      ESL_SQFILE  *sqfp;	/* open input sequence file                */
      ESL_SQ      *sq;		/* an input sequence                       */

      status = esl_sqfile_Open(infile, infmt, NULL, &sqfp);
      if (status == eslENOTFOUND)
	esl_fatal("Couldn't open seqfile %s\n", infile);
      else if (status == eslEFORMAT)
	esl_fatal("Couldn't determine format of seqfile %s\n", infile);      
      else if (status == eslEINVAL)
	esl_fatal("Can't autodetect format of stdin or .gz; use --informat\n");
      else if (status != eslOK)
	esl_fatal("Open of seqfile %s failed, code %d\n", infile, status);
      
      if ( esl_opt_IsOn(go, "--ignore"))  esl_sqio_Ignore  (sqfp, esl_opt_GetString(go, "--ignore"));
      if ( esl_opt_IsOn(go, "--acceptx")) esl_sqio_AcceptAs(sqfp, esl_opt_GetString(go, "--acceptx"), 'X');

      sq  = esl_sq_Create();
      idx = 0;
      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  if (force_lower) symconvert(sq->seq, 
				      "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
				      "abcdefghijklmnopqrstuvwxyz");
	  if (force_upper) symconvert(sq->seq, 
				      "abcdefghijklmnopqrstuvwxyz",
				      "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	  if (force_rna)   symconvert(sq->seq, "Tt", "Uu");
	  if (force_dna)   symconvert(sq->seq, "Uu", "Tt");
	  if (iupac_to_n)  symconvert(sq->seq, 
				      "RYMKSWHBVDrymkswhbvd",
				      "NNNNNNNNNNnnnnnnnnnn");
	  if (x_is_bad)    symconvert(sq->seq, "Xx", "Nn");
	  
	  if (wussify && sq->ss != NULL) esl_kh2wuss(sq->ss, sq->ss);	    
	  if (dewuss  && sq->ss != NULL) esl_wuss2kh(sq->ss, sq->ss);	    

	  if (fullwuss && sq->ss != NULL)
	    {
	      status = esl_wuss_full(sq->ss, sq->ss);
	      if (status == eslESYNTAX) 
		esl_fatal("Bad SS for %s: not in WUSS format\n", sq->name);
	      else if (status != eslOK)
		esl_fatal("Conversion of SS for %s failed, code %d\n", 
			  sq->name, status);
	    }

	  if (rename) esl_sq_FormatName(sq, "%s.%d", rename, idx+1);

	  esl_sqio_Write(ofp, sq, outfmt, FALSE);
	  esl_sq_Reuse(sq);
	  idx++;
	}
      /* status should be eslEOF on normal end; if it isn't, deal w/ error */
      if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					       sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					       status, sqfp->filename);
      
      esl_sq_Destroy(sq);
      esl_sqfile_Close(sqfp);
    } /* end of unaligned seq conversion */

  if (ofp != stdout) fclose(ofp);
  esl_getopts_Destroy(go);
  exit(0);
}

/* symconvert()
 * 
 * single seq version of esl_msa_SymConvert(); see
 * documentation there.
 * 
 * no reason yet to include in sqio API, but that may change.
 * 
 * inefficient to use this for upper/lower case conversion,
 * prob by an order of magnitude (because of the strchr() call,
 * which could be replaced by a range test), but I bet it's
 * unnoticeable.
 */
static void
symconvert(char *s, char *oldsyms, char *newsyms)
{
  int   pos;
  char *sptr;
  int   special;

  special = (strlen(newsyms) == 1 ? TRUE : FALSE);

  for (pos = 0; s[pos] != '\0'; pos++)
    if ((sptr = strchr(oldsyms, s[pos])) != NULL)
      s[pos] = (special ? *newsyms : newsyms[sptr-oldsyms]);
}


/* msafile_getline():
 * load the next line of <afp> into <afp->buf>. 
 * Returns eslOK on success, eslEOF on normal eof.
 * Throws eslEMEM on alloc failure.
 */
static int
msafile_getline(ESL_MSAFILE *afp)
{
  int status;
  status = esl_fgets(&(afp->buf), &(afp->buflen), afp->f);
  afp->linenumber++;
  return status;
}


/* regurgitate_pfam_as_afa()
 * 
 * Given an open Pfam formatted msafile, read the next alignment and 
 * regurgitate it in aligned FASTA (AFA) format without storing
 * it in a esl_msa data structure.
 * 
 * We need to do two passes through the file because in Pfam
 * sequence accessions (#=GS <seqname> AC) and sequence descriptions
 * (#=GS <seqname> DE) appear altogether before any aligned sequence
 * data, while in AFA they appear on the same line as the sequence
 * name (accession, then description).
 *
 * Example: 
 * # STOCKHOLM 1.0
 * #=GS tRNA1 AC RF00005-1
 * #=GS tRNA2 AC RF00005-2
 * #=GS tRNA1 DE first tRNA
 * #=GS tRNA2 DE second tRNA
 * 
 * tRNA1 GCGGAUUUAGCUCAGUUGGG.AGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
 * tRNA2 UCCGAUAUAGUGUAAC.GGCUAUCACAUCACGCUUUCACCGUGGAGA.CCGGGGUUCGACUCCCCGUAUCGGAG
 * 
 * converts to AFA:
 * >tRNA1 RF00005-1 first tRNA
 * GCGGAUUUAGCUCAGUUGGG.AGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAU
 * CCACAGAAUUCGCA
 * >tRNA2 RF00005-2 second tRNA
 * UCCGAUAUAGUGUAAC.GGCUAUCACAUCACGCUUUCACCGUGGAGA.CCGGGGUUCGAC
 * UCCCCGUAUCGGAG
 * 
 * In the first pass, output the sequence names and accessions we find
 * as '#=GS <seqname> AC' lines in the Pfam alignment to an accession
 * tmpfile, and output sequence names and descriptions we find as 
 * as '#=GS <seqname> DE' lines in the Pfam alignment to a description
 * tmpfile.
 *
 * In the second pass, rewind all (up to 3) files: <ac_tmpfile>,
 * <de_tmpfile> and the Pfam alignment file and start reading them
 * again.  As we're reading them, output the accessions, descriptions
 * and aligned sequence data in the proper order to an aligned FASTA
 * file.
 * 
 * Set <ret_reached_eof> as TRUE if the alignment read and reformatted
 * appears to be the only one remaining in afp.  Set <ret_reached_eof>
 * as FALSE if afp appears to include at least one more alignment.
 * 
 * Returns void. Dies upon any input error.
 */
static void
regurgitate_pfam_as_afa(ESL_MSAFILE *afp, FILE *ofp, char *alifile, char *gapsym, int force_lower, int force_upper, int force_rna, int force_dna, int iupac_to_n, int x_is_bad, char *rename, int *ret_reached_eof)
{
  char      *s = NULL;
  int        status;
  int        status2;
  char      *seqname = NULL;
  char      *first_seqname = NULL;
  char      *text = NULL;
  char      *tag = NULL;
  char      *gs = NULL;
  int        nseq_read = 0;
  int        reached_eof;
  /* variables related to reading accessions */
  char       ac_tmpfile[16] = "esltmpXXXXXX";
  FILE      *ac_fp = NULL; /* file ptr for accession tmpfile */
  char      *ac_buf = NULL;	/* buffer for line input w/ sre_fgets()      */
  int        ac_buflen = 0;	/* current allocated length for buf          */
  char      *ac_s = NULL;	        
  char      *ac_seqname = NULL;
  char      *ac = NULL;
  int        have_ac = FALSE;
  /* variables related to reading descriptions */
  char       de_tmpfile[16] = "esltmpXXXXXX";
  FILE      *de_fp = NULL; /* file ptr for description tmpfile */
  char      *de_buf = NULL;	/* buffer for line input w/ sre_fgets()      */
  int        de_buflen = 0;	/* current allocated length for buf          */
  char      *de_s = NULL;	        
  char      *de_seqname = NULL;
  char      *de = NULL;
  int        have_de = FALSE;
  /* variables related to printing out sequences */
  char      *aseq = NULL;
  int        aseqlen = 0;
  int64_t    pos;
  char       aseqbuf[61];
  int        acpl;       /* actual number of character per line */
  
  if (feof(afp->f))  { esl_fatal("--small parse error, no alignments read"); }
  afp->errbuf[0] = '\0';
   
  /**************************************************************************************************
   * First pass, go through each line of the Pfam file and output all GS DE annotation to a tmpfile *
   **************************************************************************************************/

  /* Check the magic Stockholm header line, allowing blank lines */
  do {
    status = msafile_getline(afp);
    if     (status == eslEOF) return;
    else if(status != eslOK)  esl_fatal("--small parse error. problem reading line %d of msafile", afp->linenumber);
  } while (esl_str_IsBlank(afp->buf));
    
  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0) esl_fatal("--small parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);
  while ((status2 = msafile_getline(afp)) == eslOK) {
    s = afp->buf;
    while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */
    
    if (*s == '#') { /* only lines we'll regurgitate are '#=GS DE' lines, we don't even check other lines for validity */
      if (strncmp(s, "#=GS", 4) == 0) {
	/* we don't validate the sequence exists, this would require storing all seqnames */
	s = afp->buf;
	if (esl_strtok(&s, " \t\n\r", &gs)      != eslOK) esl_fatal("small mem parse failed (line %d): bad #=GS line", afp->linenumber);
	if (esl_strtok(&s, " \t\n\r", &seqname) != eslOK) esl_fatal("small mem parse failed (line %d): bad #=GS line", afp->linenumber);
	if (esl_strtok(&s, " \t\n\r", &tag)     != eslOK) esl_fatal("small mem parse failed (line %d): bad #=GS line", afp->linenumber);
	if (esl_strtok(&s, "\n\r",    &text)    != eslOK) esl_fatal("small mem parse failed (line %d): bad #=GS line", afp->linenumber);
	if(strcmp(tag , "AC") == 0) { 
	  if(ac_fp == NULL) { /* first AC line, open tmpfile */
	    if (esl_tmpfile(ac_tmpfile, &ac_fp) != eslOK) esl_fatal("small mem parse failed, unable to open accession tmpfile");
	  }
	  fprintf(ac_fp, "%s %s\n", seqname, text);
	}
	if(strcmp(tag , "DE") == 0) { 
	  if(de_fp == NULL) { /* first DE line, open tmpfile */
	    if (esl_tmpfile(de_tmpfile, &de_fp) != eslOK) esl_fatal("small mem parse failed, unable to open description tmpfile");
	  }
	  fprintf(de_fp, "%s %s\n", seqname, text);
	}
      }
    } /* end of 'if (*s == '#')' */ 
    else if (strncmp(s, "//",   2) == 0) {
      reached_eof = FALSE;
      /* End of alignment. Normal way out, but check to see if there are any more alignments first */
      do {
	status = msafile_getline(afp);
	if (status == eslEOF) { 
	  reached_eof = TRUE;
	  break;
	}
	else if(status != eslOK)  esl_fatal("--small parse error. problem reading line %d of msafile", afp->linenumber);
      } while (esl_str_IsBlank(afp->buf));
      if(! reached_eof && strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0) esl_fatal("--small parse failed (line %d) unexpected lines after the end of first alignment", afp->linenumber);
      /* else reached_eof == FALSE: more alignments exist, go ahead and regurgitate the first alignment, inform caller we didn't read final alignment by setting <ret_reached_eof> as FALSE */
      break; /* normal way out */
    }
  }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) esl_fatal("--small parse failed (line %d): didn't find // at end of alignment", afp->linenumber);
  /* done with pass 1, close alifile and reopen it */
  esl_msafile_Close(afp);
  if((status = esl_msafile_Open(alifile, eslMSAFILE_PFAM, NULL, &afp)) != eslOK) esl_fatal("--small, second pass, unable to open file %s for reading", alifile);
  
  /******************************************************************************************
   * Pass 2, rewind the tmpfiles and the alifile, step through each, outputting appropriately
   ******************************************************************************************/
  if(ac_fp != NULL) { /* open the tmpfile with the seq accessions */
    rewind(ac_fp);
    if((status = esl_fgets(&(ac_buf), &(ac_buflen), ac_fp)) != eslOK) esl_fatal("--small accession tmpfile parse failed");
    ac_s = ac_buf;
    if (esl_strtok_adv(&ac_s, " \t\n\r", &ac_seqname, NULL, NULL) != eslOK) esl_fatal("--small accession tmpfile parse failed");
    if (esl_strtok_adv(&ac_s, "\n\r",    &ac,         NULL, NULL) != eslOK) esl_fatal("--small accession tmpfile parse failed");
  }
  if(de_fp != NULL) { /* open the tmpfile with the seq descriptions */
    rewind(de_fp);
    if((status = esl_fgets(&(de_buf), &(de_buflen), de_fp)) != eslOK) esl_fatal("--small description tmpfile parse failed");
    de_s = de_buf;
    if (esl_strtok_adv(&de_s, " \t\n\r", &de_seqname, NULL, NULL) != eslOK) esl_fatal("--small description tmpfile parse failed");
    if (esl_strtok_adv(&de_s, "\n\r",    &de,         NULL, NULL) != eslOK) esl_fatal("--small description tmpfile parse failed");
  }
  /* Check the magic Stockholm header line, allowing blank lines */
  do {
    status = msafile_getline(afp);
    if     (status == eslEOF) return;
    else if(status != eslOK)  esl_fatal("--small parse error pass 2. problem reading line %d of msafile", afp->linenumber);
  } while (esl_str_IsBlank(afp->buf));
  
  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0) esl_fatal("--small parse pass 2 failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);
  while ((status2 = msafile_getline(afp)) == eslOK) {
    s = afp->buf;
    while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */
    if (*s == '#' || *s == '\n' || *s == '\r') { /* comment or blank line, do nothing */ ; } 
    else if (strncmp(s, "//",   2) == 0)    { break; /* normal way out */ }
    else { /* sequence line */
      /* parse line into temporary strings */
      s = afp->buf;
      if (esl_strtok_adv(&s, " \t\n\r", &seqname, NULL,  NULL) != eslOK) esl_fatal("--small parse pass 2 failed (line %d): bad sequence line", afp->linenumber);
      if (esl_strtok_adv(&s, " \t\n\r", &aseq, &aseqlen, NULL) != eslOK) esl_fatal("--small parse pass 2 failed (line %d): bad sequence line", afp->linenumber);
      /* make sure we haven't just read a second line of the first sequence in file (we must be in Pfam 1 line/seq file) */
      if(nseq_read == 0) { if ((status = esl_strdup(seqname, -1, &(first_seqname))) != eslOK) esl_fatal("small mem parse failed unable to copy seqname"); }
      else if(strcmp(first_seqname, seqname) == 0) { esl_fatal("--small parse pass 2 failed (line %d): two seqs named %s. Alignment appears to be in interleaved Stockholm (not Pfam) format.", afp->linenumber, seqname); }
      nseq_read++;
      /* determine if we have an accesion and/or a description for this sequence */
      have_de = have_ac = FALSE;
      if(ac_seqname != NULL && (strcmp(ac_seqname, seqname) == 0)) have_ac = TRUE;
      if(de_seqname != NULL && (strcmp(de_seqname, seqname) == 0)) have_de = TRUE;
      if(rename != NULL) { 
	fprintf(ofp, ">%s.%d%s%s%s%s\n", rename, nseq_read, (have_ac ? " " : "") , (have_ac ? ac : ""), (have_de ? " " : "") , (have_de ? de : "")); 
      }
      else {
	fprintf(ofp, ">%s%s%s%s%s\n", seqname, (have_ac ? " " : "") , (have_ac ? ac : ""), (have_de ? " " : "") , (have_de ? de : "")); 
      }
      if(have_ac) {
	status = esl_fgets(&(ac_buf), &(ac_buflen), ac_fp);
	if(status == eslEOF) { 
	  ac_seqname = NULL;
	}
	else if (status == eslOK) { 
	  ac_s = ac_buf;
	  if (esl_strtok_adv(&ac_s, " \t\n\r", &ac_seqname, NULL, NULL) != eslOK) esl_fatal("--small accession tmpfile parse failed");
	  if (esl_strtok_adv(&ac_s, "\n\r",    &ac,         NULL, NULL) != eslOK) esl_fatal("--small accession tmpfile parse failed");
	}
      }
      if(have_de) {
	status = esl_fgets(&(de_buf), &(de_buflen), de_fp);
	if(status == eslEOF) { 
	  de_seqname = NULL;
	}
	else if (status == eslOK) { 
	  de_s = de_buf;
	  if (esl_strtok_adv(&de_s, " \t\n\r", &de_seqname, NULL, NULL) != eslOK) esl_fatal("--small description tmpfile parse failed");
	  if (esl_strtok_adv(&de_s, "\n\r",    &de,         NULL, NULL) != eslOK) esl_fatal("--small description tmpfile parse failed");
	}
      }

      /* now print sequence, after converting symbols as nec */
      if (gapsym != NULL) symconvert(aseq, "-_.", gapsym);
      if (force_lower)    symconvert(aseq,
				           "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
					   "abcdefghijklmnopqrstuvwxyz");
      if (force_upper)    symconvert(aseq,
				  	   "abcdefghijklmnopqrstuvwxyz",
					   "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      if (force_rna)      symconvert(aseq, "Tt", "Uu");
      if (force_dna)      symconvert(aseq, "Uu", "Tt");
      if (iupac_to_n)     symconvert(aseq, 
					   "RYMKSWHBVDrymkswhbvd",
					   "NNNNNNNNNNnnnnnnnnnn");
      if (x_is_bad)     symconvert(aseq,   "Xx", "Nn");

      pos = 0;
      while (pos < aseqlen) { 
	acpl = (aseqlen - pos > 60)? 60 : aseqlen - pos;
	strncpy(aseqbuf, aseq + pos, acpl);
	aseqbuf[acpl] = '\0';
	fprintf(ofp, "%s\n", aseqbuf);	      
	pos += 60;
      }
    }
  }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK)   esl_fatal("--small parse pass 2 failed (line %d): didn't find // at end of alignment", afp->linenumber);
  if (de_seqname != NULL) esl_fatal("--small parse pass 2 failed, sequence %s with #=GS DE line does not exist in alignment or is in different order.", de_seqname);

  if(ac_fp != NULL) fclose(ac_fp);
  if(de_fp != NULL) fclose(de_fp);
  esl_msafile_Close(afp);

  if(first_seqname != NULL) free(first_seqname);
  if(ac_buf != NULL)        free(ac_buf);
  if(de_buf != NULL)        free(de_buf);

  *ret_reached_eof = reached_eof;
  return;
}

/* determine_spacelen
 *                   
 * Determine number of consecutive ' ' characters 
 * in the string pointed to by s.
 */
static int
determine_spacelen(char *s)
{
  int spacelen = 0;
  while (*s == ' ') { spacelen++; s++; } 
  return spacelen;
}


/* regurgitate_pfam_as_pfam()
 * 
 * Given an open Pfam formatted msafile, read the next alignment and
 * regurgitate it, after modifying it as necessary (change dna to rna,
 * wussify SS, etc) in Pfam format.
 * 
 * Returns <eslOK> on success. 
 * Returns <eslEOF> if there are no more alignments in <afp>.
 * Returns <eslEFORMAT> if parse fails because of a file format
 * problem, in which case afp->errbuf is set to contain a formatted
 * message that indicates the cause of the problem.
 */
static int
regurgitate_pfam_as_pfam(ESL_MSAFILE *afp, FILE *ofp, char *gapsym, int force_lower, int force_upper, int force_rna, int force_dna, int iupac_to_n, int x_is_bad, int wussify, int dewuss, int fullwuss)
{
  char      *s = NULL;
  int        status;
  int        status2;
  char      *gc = NULL;
  char      *gr = NULL;
  char      *seqname = NULL;
  char      *first_seqname = NULL;
  int        namelen;
  char      *tag = NULL;
  int        taglen;
  char      *text = NULL;
  int        textlen = 0;
  int        nseq_read = 0;
  int        parse_gc_and_gr;
  int        spacelen, spacelen2;
  int        flushpoint = 10000;
  int        exp_alen = -1;

  parse_gc_and_gr = (wussify || dewuss || fullwuss) ? TRUE : FALSE; /* should we parse out GR/GC lines and check if they're SS lines? */

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';
   
  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR; /* includes EOF  */
  } while (esl_str_IsBlank(afp->buf));
  
  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);
  fprintf(ofp, "%s", afp->buf);

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) { 
    if(afp->linenumber % flushpoint == 0) fflush(ofp);
    s = afp->buf;
    while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */
    
    if (*s == '#') {
      if (parse_gc_and_gr && strncmp(s, "#=GC", 4) == 0) { 
	/* parse line into temporary strings */
	s = afp->buf;
	if (esl_strtok    (&s, " \t\n\r", &gc)                  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	if (esl_strtok_adv(&s, " \t\n\r", &tag,  &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	spacelen = determine_spacelen(s);
	if (esl_strtok_adv(&s, " \t\n\r", &text, &textlen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	
	/* verify alignment length */
	if(exp_alen == -1) exp_alen = textlen;
	else if(exp_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line, len %d, expected %d", afp->linenumber, textlen, exp_alen);
	
	if(strcmp(tag, "SS_cons") == 0) {
	  if     (wussify)  esl_kh2wuss(text, text);
	  else if(dewuss)   esl_wuss2kh(text, text);
	  else if(fullwuss) { 
	    status = esl_wuss_full(text, text);
	    if      (status == eslESYNTAX) esl_fatal("Bad SS_cons line: not in WUSS format, alifile line: %d", afp->linenumber);
	    else if (status != eslOK)      esl_fatal("Conversion of SS_cons line failed, code %d, alifile line: %d", status, afp->linenumber);
	  }
	}		  
	fprintf(ofp, "#=GC %-*s %s\n", taglen+spacelen, tag, text);
      }
      else if (parse_gc_and_gr && strncmp(s, "#=GR", 4) == 0) { 
	/* parse line into temporary strings */
	s = afp->buf;
	if (esl_strtok    (&s, " \t\n\r", &gr)                      != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad #=GR line", afp->linenumber);
	if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad #=GR line", afp->linenumber);
	spacelen = determine_spacelen(s);
	if (esl_strtok_adv(&s, " \t\n\r", &tag,    &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad #=GR line", afp->linenumber);
	spacelen2 = determine_spacelen(s);
	if (esl_strtok_adv(&s, " \t\n\r", &text,  &textlen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad #=GR line", afp->linenumber);
	/* verify alignment length */
	if      (exp_alen == -1)      exp_alen = textlen;
	else if (exp_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad seq line, len %d, expected %d", afp->linenumber, textlen, exp_alen);
	
	if(strcmp(tag, "SS") == 0) {
	  if     (wussify)  esl_kh2wuss(text, text);
	  else if(dewuss)   esl_wuss2kh(text, text);
	  else if(fullwuss) { 
	    status = esl_wuss_full(text, text);
	    if      (status == eslESYNTAX) esl_fatal("Bad SS line: not in WUSS format, alifile line: %d", afp->linenumber);
	    else if (status != eslOK)      esl_fatal("Conversion of SS line failed, code %d, alifile line: %d", status, afp->linenumber);
	  }
	}		  
	fprintf(ofp, "#=GR %-*s %-*s %s\n", namelen+spacelen, seqname, taglen+spacelen2, tag, text);
      }
      else { /* '#' prefixed line that is not #=GR (or it is #=GR and wussify,dewuss,fullwuss are all FALSE) */
	fprintf(ofp, "%s", afp->buf); /* print the line */
      }
    } /* end of 'if (*s == '#')' */ 
    else if (strncmp(s, "//",   2) == 0)   { break; /* normal way out */ }
    else if (*s == '\n' || *s == '\r')     { fprintf(ofp, "%s", afp->buf); continue; }
    else { /* sequence line */
      /* parse line into temporary strings */
      s = afp->buf;
      if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad sequence line", afp->linenumber);
      spacelen = determine_spacelen(s);
      if (esl_strtok_adv(&s, " \t\n\r", &text,    &textlen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): bad sequence line", afp->linenumber);
      /* verify alignment length */
      if     (exp_alen == -1)      exp_alen = textlen;
      else if(exp_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad seq line, len %d, expected %d", afp->linenumber, textlen, exp_alen);
      
      /* make sure we haven't just read a second line of the first sequence in file (we must be in Pfam 1 line/seq file) */
      if(nseq_read == 0) { if ((status = esl_strdup(seqname, -1, &(first_seqname))) != eslOK) goto ERROR; }
      else if(strcmp(first_seqname, seqname) == 0) { ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): two seqs named %s. Alignment appears to be in Stockholm format. Reformat to Pfam with esl-reformat.", afp->linenumber, seqname); }
      nseq_read++;
      
      /* make adjustments as necessary */
      if (gapsym != NULL) symconvert(text, "-_.", gapsym);
      else                symconvert(text, "-_.", "-"); /* AFA default is to have gap characters as '-' */
      if (force_lower)    symconvert(text,
				     "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
				     "abcdefghijklmnopqrstuvwxyz");
      if (force_upper)    symconvert(text,
				     "abcdefghijklmnopqrstuvwxyz",
				     "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      if (force_rna)      symconvert(text, "Tt", "Uu");
      if (force_dna)      symconvert(text, "Uu", "Tt");
      if (iupac_to_n)     symconvert(text, 
				     "RYMKSWHBVDrymkswhbvd",
				     "NNNNNNNNNNnnnnnnnnnn");
      if (x_is_bad)     symconvert(text,   "Xx", "Nn");
      /* print it out */
      fprintf(ofp, "%-*s %s\n", namelen+spacelen, seqname, text);
    }
  }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "--small parse failed (line %d): didn't find // at end of alignment", afp->linenumber);
  if (first_seqname     != NULL) free(first_seqname);
  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

