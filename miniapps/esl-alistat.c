/* Show statistics about a multiple sequence alignment file or MSA database.
 * 
 * From squid's alistat (1995)
 * SRE, Wed May 16 08:23:23 2007 [Janelia] [Philip Glass, The Fog of War]
 * SVN $Id$
 */
#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_msa.h"
#include "esl_distance.h"

static char banner[] = "show summary statistics for a multiple sequence alignment file";
static char usage[]  = "[options] <msafile>\n\
The <msafile> must be in Stockholm format.";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-1",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "use tabular output, one line per alignment",              0 },
  { "--amino",  eslARG_NONE,"default",NULL,NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   0 },
  { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       0 },
  { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       0 },
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
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed with error %d\n", status);

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

  /***********************************************
   * Read MSAs one at a time.
   ***********************************************/

  if (esl_opt_GetBoolean(go, "-1")) {
    puts("#");
    printf("# %-4s %-20s %10s %7s %7s %12s %6s %6s %10s %3s\n", "idx", "name", "format", "nseq", "alen", "nres", "small", "large", "avlen", "%id");
    printf("# %-4s %-20s %10s %7s %7s %12s %6s %6s %10s %3s\n", "----", "--------------------", "----------", "-------", "-------", "------------", "------", "------", "----------", "---");
  }

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

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
      
      if (esl_opt_GetBoolean(go, "-1")) 
	{
	  printf("%-6d %-20s %10s %7d %7d %12llu %6d %6d %10.1f %3.0f\n",
		 nali, 
		 msa->name,
		 esl_msa_DescribeFormat(afp->format),
		 msa->nseq,
		 msa->alen,
		 (unsigned long long) nres,
		 small,
		 large,
		 (double) nres / (double) msa->nseq,
		 100.*avgid);
	}
      else
	{
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
	  printf("Average identity:    %.0f%%\n", 100.*avgid); 
	  printf("//\n");
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
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  return 0;
}


