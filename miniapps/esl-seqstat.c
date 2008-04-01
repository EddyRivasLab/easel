/* Simple statistics on a sequence file
 * 
 * SRE, Sun Feb 24 15:33:53 2008 [UA5315 to St. Louis]
 * SVN $Id$  
 * from squid's seqstat (1994)
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

static char banner[] = "show simple statistics on a sequence file";
static char usage1[] = "   [options] <seqfile>";

#define ALPH_OPTS "--rna,--dna,--amino" /* toggle group, alphabet type options          */

static ESL_OPTIONS options[] = {
  /* name         type           default   env range togs  reqs  incomp      help                                      docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL,      NULL, "help; show brief info on version and usage",          1 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL, NULL,      NULL, "specify that input file is in format <s>",            1 },
  { "--rna",      eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains RNA sequence",        1 },
  { "--dna",      eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains DNA sequence",        1 },
  { "--amino",    eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains protein sequence",    1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage1);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage1);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = NULL;
  char           *seqfile   = NULL;
  ESL_SQFILE     *sqfp      = NULL;
  int             infmt     = eslSQFILE_UNKNOWN;
  int             alphatype = eslUNKNOWN;
  ESL_ALPHABET   *abc       = NULL;
  ESL_SQ         *sq        = NULL;
  long long       nseq      = 0;   
  long long       nres      = 0;
  long long       small     = 0;
  long long       large     = 0;
  double         *monoc     = NULL; /* monoresidue composition per sequence  */
  double         *monoc_all = NULL; /* monoresidue composition over all seqs */
  int             status    = eslOK;
  int             i;
  int             x;

  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) 
    cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK)
    cmdline_failure(argv[0], "Error in app configuration: %s\n",   go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )
    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go) != 1) 
    cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");

  seqfile = esl_opt_GetArg(go, 1);

  if (esl_opt_GetString(go, "--informat") != NULL) {
    infmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }

  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file %s", seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", seqfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
    if      (status == eslEAMBIGUOUS) esl_fatal("Couldn't guess alphabet from first sequence in %s", seqfile);
    else if (status == eslEFORMAT)    esl_fatal("Sequence file parse error, line %d of file %s:\n%s\n",
					       sqfp->linenumber, seqfile, sqfp->errbuf);
    else if (status == eslENODATA)    esl_fatal("Sequence file %s contains no data?", seqfile);
    else if (status != eslOK)         esl_fatal("Failed to guess alphabet (error code %d)\n", status);
  }
  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);

  ESL_ALLOC(monoc,     (abc->Kp) * sizeof(double));  
  ESL_ALLOC(monoc_all, (abc->Kp) * sizeof(double));  
  esl_vec_DSet(monoc_all, abc->Kp, 0.0);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      esl_vec_DSet(monoc, abc->Kp, 0.0);
      for (i = 1; i <= sq->n; i++) 
	monoc[sq->dsq[i]]++;
      esl_vec_DAdd(monoc_all, monoc, abc->Kp);

      if (nseq == 0) { small = large = sq->n; }
      else {
	small = ESL_MIN(small, sq->n);
	large = ESL_MAX(large, sq->n);
      }

      nres += sq->n;
      nseq++;
      esl_sq_Reuse(sq);
    }

  printf("Format:              %s\n",   esl_sqio_DescribeFormat(sqfp->format));
  printf("Alphabet type:       %s\n",   esl_abc_DescribeType(abc->type));
  printf("Number of sequences: %lld\n", nseq);
  printf("Total # residues:    %lld\n", nres);
  printf("Smallest:            %lld\n", small);
  printf("Largest:             %lld\n", large);
  printf("Average length:      %.1f\n", (float) nres / (float) nseq);


  printf("\nResidue composition:\n");
  for (x = 0; x < abc->Kp; x++)
    if (x < abc->K || monoc_all[x] > 0)
      printf("residue: %c   %10.0f  %.4f\n", abc->sym[x], monoc_all[x], monoc_all[x] / (double) nres);

  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  return status;
}
