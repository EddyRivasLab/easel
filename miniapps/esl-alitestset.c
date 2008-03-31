/* Construct a training alignment/test sequences set from an MSA.
 * 
 * This procedure is used in constructing our internal RMARK and PMARK
 * benchmarks.
 * 
 * SRE, Thu Mar 27 11:05:38 2008 [Janelia]
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_vectorops.h"


static char banner[] = "show summary statistics for a multiple sequence alignment file";
static char usage[]  = "[options] <msafile>\n\
The <msafile> must be in Stockholm format; it can be a multi-MSA file.";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-F",       eslARG_REAL, "0.70", NULL,"0<x<=1.0",NULL,NULL,NULL,         "filter out seqs <x*average length",                       0 },
  { "--amino",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   0 },
  { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       0 },
  { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}


int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r     = NULL;	/* random number generator         */
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  char         *alifile = NULL;	/* alignment file name             */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_ALPHABET *abc     = NULL;	/* biological alphabet             */
  ESL_MSA      *origmsa = NULL;	/* one multiple sequence alignment */
  ESL_MSA      *msa     = NULL;	/* MSA after frags are removed     */
  int           status;		/* easel return code               */
  int           nali;		/* number of alignments read       */
  int           i;		/* counter over seqs               */
  int           rlen;		/* a raw (unaligned) seq length    */
  int           small, large;	/* smallest, largest sequence      */
  uint64_t      nres;		/* total # of residues in msa      */
  double        avgid;		/* average fractional pair id      */
  int    max_comparisons;       /* maximum # comparisons for avg id */
  
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go)                  != 1)     cmdline_failure(argv[0], "Incorrect number of command line arguments\n");

  alifile = esl_opt_GetArg(go, 1);
  fmt     = eslMSAFILE_STOCKHOLM;

  r = esl_randomness_Create(42);

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

  while ((status = esl_msa_Read(afp, &origmsa)) == eslOK)
    {
      remove_fragments(go, origmsa, &msa);
      separate_sets(r, msa, idthresh1, idthresh2, &trainmsa, &teststack);


      esl_msa_Destroy(trainmsa);
      esl_msa_Destroy(origmsa);
      esl_msa_Destroy(msa);
    }
  if      (status == eslEFORMAT)  esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
					    afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)      esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)           esl_fatal("No alignments found in file %s\n", alifile);

  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  return 0;
}
      
  
/* Step 1. Label all sequence fragments < fragfrac of average raw length */
static int
remove_fragments(ESL_GETOPTS *go, ESL_MSA *msa, ESL_MSA **ret_filteredmsa)
{
  int     *useme    = NULL;
  double   len      = 0.0;
  double   fragfrac = esl_opt_GetReal(go, "-F");
  int      status;

  for (i = 0; i < msa->nseq; i++) 
    len += esl_abc_dsqrlen(msa->abc, msa->ax[i]);
  len *= fragfrac / (double) msa->nseq;

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) 
    useme[i] = (esl_abc_dsqrlen(msa->abc, msa->ax[i]) < len) ? 0 : 1;

  return esl_msa_SequenceSubset(msa, useme, ret_filteredmsa);

 ERROR:
  if (useme != NULL) free(useme); 
  *ret_filteredmsa = NULL;
  return status;
}

/* Step 2. Extract the training set and test set.
 */
static int
separate_sets(ESL_RANDOMNESS *r, ESL_MSA *msa, double idthresh1, double idthresh2,
	      ESL_MSA *ret_trainmsa, ESL_STACK *ret_teststack)
{      
  ESL_MSA   *trainmsa  = NULL;
  ESL_MSA   *test_msa  = NULL;
  ESL_STACK *teststack = NULL;
  int *assignment = NULL;
  int *nin        = NULL;
  int *useme      = NULL;
  int  nc         = 0;
  int  c;
  int  ctrain;			/* index of the cluster that becomes the training alignment */
  int  ntrain;			/* number of seqs in the training alignment */
  int  nskip;
  int  i;
  int  status;

  if ((teststack = esl_stack_PCreate()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(useme, sizeof(int) * msa->nseq);

  if ((status = esl_msacluster_SingleLinkage(msa, idthresh1, &assignment, &nin, &nc)) != eslOK) goto ERROR;
  ctrain = esl_vec_IArgMax(nin, nc);
  ntrain = esl_vec_IMax(nin, nc);

  for (i = 0; i < msa->nseq; i++) useme[i] = (!is_frag[i] && assignment[i] == ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &trainmsa)) != eslOK) goto ERROR;

  free(nin);         nin        = NULL;
  free(assignment);  assignment = NULL;

  /* Put all the other sequences into an MSA of their own; from these, we'll
   * choose test sequences.
   */
  for (i = 0; i < msa->nseq; i++) useme[i] = (!is_frag[i] && assignment[i] != ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &test_msa))                        != eslOK) goto ERROR;
  if ((status = esl_msacluster_SingleLinkage(msa, idthresh2, &assignment, &nin, &nc)) != eslOK) goto ERROR;

  for (c = 0; c < nc; c++)
    {
      nskip = esl_rnd_Choose(r, nin[c]); /* pick a random seq in this cluster to be the test. */
      for (i=0; nskip >= 0 && i < msa->nseq; i++)
	if (assignment[i] == c) nskip--;
      
      esl_sq_FetchFromMSA(msa, i, &sq);
      esl_stack_PPush(teststack, (void *) sq);
    }

  esl_msa_Destroy(test_msa);
  free(useme);
  free(nin);
  free(assignment);

  *ret_trainmsa  = msa;
  *ret_teststack = teststack;
  return eslOK;

 ERROR:
  if (useme      != NULL) free(useme);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  esl_msa_Destroy(trainmsa); 
  esl_msa_Destroy(test_msa); 
  while (esl_stack_PPop(teststack, &sq) == eslOK) esl_sq_Destroy(sq);
  esl_stack_Destroy(teststack);
  *ret_trainmsa  = NULL;
  *ret_teststack = NULL;
  return status;
}

