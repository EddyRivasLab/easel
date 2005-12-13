/* sreformat - convert sequence file formats
 *             
 * SRE, Sun Feb 27 08:24:33 2005
 * from squid's sreformat (1993).
 * SVN $Id$            
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <easel.h>
#include <esl_getopts.h>
#include <esl_sqio.h>
#include <esl_msa.h>
#include <esl_wuss.h>



static char banner[] = "\
sreformat :: convert sequence file formats";

static char usage[] = "\
Usage: sreformat [-options] <format> <seqfile>\n\
  Output format choices: Unaligned      Aligned\n\
                         -----------    -------\n\
                         fasta          stockholm\n\
                                        pfam\n\
\n\
  Available options are:\n\
  -h     : help; print brief info on version and usage\n\
  -o <f> : send output to file <f>, not stdout\n\
  -d     : force DNA alphabet for nucleic acid sequence\n\
  -r     : force RNA alphabet for nucleic acid sequence\n\
  -l     : force lower case\n\
  -u     : force upper case\n\
  -x     : convert non-IUPAC chars (i.e. X's) in DNA to N's\n\
  -n     : remove IUPAC codes; convert all ambig chars in DNA to N's\n\
";

static char experts[] = "\
  Expert options:\n\
    --informat <s>: input sequence file is in format <s>\n\
    --mingap      : remove columns containing all gaps (seqfile=alignment)\n\
    --nogap       : remove columns containing any gaps (seqfile=alignment)\n\
    --gapsym <c>  : convert all gaps to character '<c>'\n\
    --wussify     : convert old format RNA structure markup lines to WUSS\n\
    --dewuss      : convert WUSS notation RNA structure markup to old format\n\
    --fullwuss    : convert simple WUSS notation to full (output) WUSS\n\
";

static ESL_OPTIONS options[] = {
   /* name          type        default env   range togs  reqs  incompat */
  { "-d",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-r" },
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL },
  { "-l",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-u" },
  { "-n",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-x" },
  { "-o",         eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, NULL },
  { "-r",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-d" },
  { "-u",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-l" },
  { "-x",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-n" },
  { "--gapsym",   eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, "--mingap,--nogap" },
  { "--informat", eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, NULL },
  { "--mingap",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--nogap", },
  { "--nogap",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--mingap,--gapsym" },
  { "--wussify",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--dewuss,--fullwuss"  },
  { "--dewuss",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--wussify,--fullwuss" },
  { "--fullwuss", eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--wussify,--dewuss" },
  { 0,0,0,0,0,0,0,0 },
};

static void symconvert(char *s, char *oldsyms, char *newsyms);

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;		/* application configuration               */
  char        *outformat;	/* output format as a string               */
  char        *infile;          /* name of input sequence file             */
  int          infmt;		/* input format as a code; eslSQFILE_FASTA */
  int          outfmt;		/* output format as a code                 */
  int          status;		/* return code from an Easel call          */
  FILE        *ofp;		/* output stream                           */

  char  *informat;		/* input format as string; "fasta"           */
  char  *outfile;		/* output file, or NULL                      */
  int    show_help;		/* TRUE to show help/usage                   */
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

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  esl_opt_GetBooleanOption(go, "-d",         &force_dna);
  esl_opt_GetBooleanOption(go, "-h",         &show_help);
  esl_opt_GetBooleanOption(go, "-l",         &force_lower);
  esl_opt_GetBooleanOption(go, "-n",         &iupac_to_n);
  esl_opt_GetStringOption( go, "-o",         &outfile);
  esl_opt_GetBooleanOption(go, "-r",         &force_rna);
  esl_opt_GetBooleanOption(go, "-u",         &force_upper);
  esl_opt_GetBooleanOption(go, "-x",         &x_is_bad);
  esl_opt_GetStringOption( go, "--gapsym",   &gapsym);
  esl_opt_GetStringOption( go, "--informat", &informat);
  esl_opt_GetBooleanOption(go, "--mingap",   &do_mingap);
  esl_opt_GetBooleanOption(go, "--nogap",    &do_nogap);
  esl_opt_GetBooleanOption(go, "--wussify",  &wussify);
  esl_opt_GetBooleanOption(go, "--dewuss",   &dewuss);
  esl_opt_GetBooleanOption(go, "--fullwuss", &fullwuss);
  
  if (show_help) 
    {
      esl_banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }

  if (esl_opt_ArgNumber(go) != 2) 
    esl_fatal("Incorrect number of command line arguments.\n%s\n", usage); 

  outformat = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);
  infile    = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);

  infmt = eslSQFILE_UNKNOWN;
  if (informat != NULL)
    {
      infmt = esl_sqfile_FormatCode(informat);
      if (infmt == eslSQFILE_UNKNOWN)
	esl_fatal("%s is not a recognized input seqfile format\n");
    }
    
  outfmt = esl_sqfile_FormatCode(outformat);
  if (outfmt == eslSQFILE_UNKNOWN)
    esl_fatal("%s is not a recognized output seqfile format\n");

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
  if (esl_sqfile_IsAlignment(outfmt))
    {
      ESL_MSAFILE *afp;
      ESL_MSA     *msa;
      int          idx;

      status = esl_msafile_Open(infile, infmt, NULL, &afp);
      if (status == eslENOTFOUND)
	esl_fatal("Alignment file %s not readable\n", infile);
      else if (status == eslEFORMAT) 
	esl_fatal("Couldn't determine format of alignment %s\n", infile);
      else if (status == eslEINVAL)
	esl_fatal("Can't autodetect format of stdin or .gz; use --informat\n");
      else if (status != eslOK) 
	esl_fatal("Alignment file open failed with error %d\n", status);

      while ((status = esl_msa_Read(afp, &msa)) == eslOK)
	{
	  if (do_mingap)    esl_msa_MinimGaps(msa, "-_.");
	  if (do_nogap)     esl_msa_NoGaps(msa, "-_.");
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
		    esl_fatal("Bad consensus SS: not in WUSS format\n",);
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
      
      /* Alignment input should end with normal status eslEOF;
       * if it didn't, deal with the problem.
       */
      if (status == eslEFORMAT)
	esl_fatal("\
Alignment file parse error, line %d of file %s:\n\
%s\n\
Offending line is:\n\
%s\n", afp->linenumber, afp->fname, afp->errbuf, afp->buf);
      else if (status != eslEOF)
	esl_fatal("Alignment file read failed with error code %d\n", status);

      esl_msafile_Close(afp);
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
      
      sq = esl_sq_Create();
      while ((status = esl_sq_Read(sqfp, sq)) == eslOK)
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

	  esl_sq_Write(ofp, sq, outfmt);
	  esl_sq_Reuse(sq);
	}

      /* status should be eslEOF on normal end; if it isn't, deal w/ error */
      if (status == eslEFORMAT)
	esl_fatal("\
Sequence file parse error, line %d of file %s:\n\
%s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);
      else if (status != eslEOF)
	esl_fatal("Sequence file %s read failed with error code %d\n",
		  sqfp->filename, status);
      
      esl_sq_Destroy(sq);
      esl_sqfile_Close(sqfp);
    } /* end of unaligned seq conversion */

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




/*****************************************************************
 * @LICENSE@
 *****************************************************************/
