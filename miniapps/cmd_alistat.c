#include "esl_config.h"

#include <stdio.h>
#include <sys/stat.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_distance.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_subcmd.h"

static ESL_OPTIONS cmd_options[] = {
  /* name             type        default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL,  "show brief help on version and usage",                 0 },
  { "-1",          eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL,  "use tabular output, one line per alignment",           0 },
  { "--dna",       eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL,  "use DNA alphabet (don't autodetect)",                  0 },
  { "--rna",       eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL,  "use RNA alphabet (don't autodetect)",                  0 },
  { "--amino",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL,  "use protein alphabet (don't autodetect)",              0 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static void alistat_default(const char *msafile, ESL_MSAFILE *afp);
static void alistat_oneline(const char *msafile, ESL_MSAFILE *afp);

int
esl_cmd_alistat(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv);
  ESL_ALPHABET   *abc     = NULL;
  char           *msafile = esl_opt_GetArg(go, 1);
  ESL_MSAFILE    *afp     = NULL;
  int             fmt     = eslMSAFILE_UNKNOWN;
  int             status;
  
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  if (( status = esl_msafile_Open(&abc, msafile, /*env=*/NULL, fmt, /*fmtd=*/NULL, &afp)) != eslOK)
    esl_msafile_OpenFailure(afp, status);

  if (esl_opt_GetBoolean(go, "-1")) alistat_oneline(msafile, afp);
  else                              alistat_default(msafile, afp);
  
  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;
}


static void
alistat_oneline(const char *msafile, ESL_MSAFILE *afp)
{
  ESL_MSA    *msa         = NULL;
  FILE       *fp          = NULL;
  int         nali        = 0;
  esl_pos_t   last_offset = -1;
  esl_pos_t   totsize;
  int64_t     rlen, smallest, largest, nres;
  double      avgid;
  int         max_comparisons = 1000;
  struct stat fileinfo;
  int64_t     recsize;
  float       ratio;
  int         i;
  int         status;

  /* Get the total file size, in bytes */
  if (( fp = fopen(msafile, "r")) == NULL)  esl_fatal("Failed to open %s as a file\n", msafile);
  fstat(fileno(fp), &fileinfo);
  totsize = fileinfo.st_size;
  fclose(fp);

  esl_dataheader(stdout,
		 -6,  "idx",
		 -20, "name",
		 -10, "format",
		 10,  "nseq",
		 10,  "alen",
		 12,  "nres",
		 6,   "small",
		 6,   "large",
		 8,   "avglen",
		 3,   "%id",
		 12,  "recsize",
		 10,  "size/nres",
		 0);  // 0 is needed to signal arglist termination

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      nali++;

      /* disk record size stats -- for *previous* msa, delayed off-by-one output */
      if (last_offset != -1)   
	{
	  recsize = msa->offset - last_offset;
	  ratio   = (float) recsize / (float) nres;  // <nres> is from the _previous_ MSA, previous loop iteration
	  printf("%12" PRId64 " %10.2f\n", recsize, ratio);
	}

      /* raw sequence length stats */
      nres = 0;
      smallest = largest = -1;
      for (i = 0; i < msa->nseq; i++)
	{
	  rlen  = esl_abc_dsqrlen(msa->abc, msa->ax[i]); 
	  nres += rlen;  // <nres> output is deferred to next time around the loop
	  if (smallest == -1 || rlen < smallest) smallest = rlen;
	  if (largest  == -1 || rlen > largest)  largest  = rlen;
	}

      /* percent identity stats */
      esl_dst_XAverageId(msa->abc, msa->ax, msa->nseq, max_comparisons, &avgid);

      printf("%-6d %-20s %10s %10d %10" PRId64 " %12" PRId64 " %6" PRId64 " %6" PRId64 " %8.1f %3.0f ",
	     nali,
	     msa->name,
	     esl_msafile_DecodeFormat(afp->format),
	     msa->nseq,
	     msa->alen,
	     nres,
	     smallest,
	     largest,
	     (double) nres / (double) msa->nseq,
	     100. * avgid);

      last_offset = msa->offset;
      esl_msa_Destroy(msa);
    }
  if (nali == 0 || status != eslEOF) esl_msafile_ReadFailure(afp, status); 

  // and for the very last msa in the file... 
  if (last_offset != -1)   
    {
      recsize = totsize - last_offset;
      ratio   = (float) recsize / (float) nres;
      printf("%12" PRId64 " %10.2f\n", recsize, ratio);
    }
}


static void
alistat_default(const char *msafile, ESL_MSAFILE *afp)
{
  ESL_MSA    *msa             = NULL;
  int         nali            = 0;
  int         max_comparisons = 1000;
  int64_t     rlen, smallest, largest, nres;
  double      avgid;
  int         i;
  int         status;

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      /* raw seq length stats */
      nres = 0;
      smallest = largest = -1;
      for (i = 0; i < msa->nseq; i++)
	{
	  rlen  = esl_abc_dsqrlen(msa->abc, msa->ax[i]); 
	  nres += rlen; 
	  if (smallest == -1 || rlen < smallest) smallest = rlen;
	  if (largest  == -1 || rlen > largest)  largest  = rlen;
	}

      /* percent identity stats */
      esl_dst_XAverageId(msa->abc, msa->ax, msa->nseq, max_comparisons, &avgid);

      printf("Alignment name:      %s\n",          msa->name);
      printf("Format:              %s\n",          esl_msafile_DecodeFormat(afp->format));
      printf("Alphabet:            %s\n",          esl_abc_DecodeType(msa->abc->type));
      printf("Number of sequences: %d\n",          msa->nseq);
      printf("Alignment length:    %" PRId64 "\n", msa->alen);
      printf("Total # residues:    %" PRId64 "\n", nres);
      printf("Smallest:            %" PRId64 "\n", smallest);
      printf("Largest:             %" PRId64 "\n", largest);
      printf("Average length:      %.1f\n",        (double) nres / (double) msa->nseq);
      printf("Average identity:    %.0f%%\n",      100.*avgid);
      printf("//\n");

      esl_msa_Destroy(msa);
      nali++;
    }
  if (nali == 0 || status != eslEOF) esl_msafile_ReadFailure(afp, status); 

}
