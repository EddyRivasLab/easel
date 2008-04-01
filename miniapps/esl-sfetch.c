/* Fetch a sequence (or part of one) from a sequence flatfile.
 * 
 * From squid's sfetch and ffetch
 * SRE, Mon Mar 31 16:12:50 2008 [Janelia] 
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_fileparser.h"
#include "esl_keyhash.h"
#include "esl_ssi.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static char banner[] = "retrieve sequence(s) from a file";
static char usage1[] = "[options] <sqfile> <name>          (retrieves one sequence named <name>)";
static char usage2[] = "[options] -f <sqfile> <namefile>   (retrieves all sequences named in <namefile>)";
static char usage3[] = "[options] --index <sqfile>         (index <sqfile>)";

static void
cmdline_failure(char *argv0, char *format, ...) 
{
  va_list argp;
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage1);
  esl_usage(stdout, argv0, usage2);
  esl_usage(stdout, argv0, usage3);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage1);
  esl_usage (stdout, argv0, usage2);
  esl_usage (stdout, argv0, usage3);
  puts("\n where options are:");
  esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
  exit(0);
}

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL,          "help; show brief info on version and usage",        0 },
  { "-f",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL,"--index",      "second cmdline arg is a file of names to retrieve", 0 },
  { "-o",       eslARG_OUTFILE,FALSE,NULL, NULL, NULL, NULL,"-O,--index",   "output sequences to file <f> instead of stdout",    0 },
  { "-O",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL,"-o,-f,--index","output sequence to file named <key>",               0 },
  { "--index",  eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL,          "index <sqfile>, creating <sqfile>.ssi",             0 },
  { "--informat",eslARG_STRING,FALSE,NULL, NULL, NULL, NULL,  NULL,         "specify that input file is in format <s>",          0 },
 { 0,0,0,0,0,0,0,0,0,0 },
};

static void create_ssi_index(ESL_GETOPTS *go, ESL_SQFILE *sqfp);
static void multifetch(ESL_GETOPTS *go, FILE *ofp, char *keyfile, ESL_SQFILE *sqfp);
static void onefetch(ESL_GETOPTS *go, FILE *ofp, char *key, ESL_SQFILE *sqfp);
static void regurgitate_one_uniprot_entry(FILE *ofp, ESL_SQFILE *sqfp);


int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	                        /* application configuration       */
  char         *seqfile = NULL;	                        /* sequence file name              */
  int           format  = eslSQFILE_UNKNOWN;		/* format code for seqfile         */
  ESL_SQFILE   *sqfp    = NULL;                         /* open sequence file              */
  FILE         *ofp     = NULL;	                        /* output stream for sequences     */
  int           status;		                        /* easel return code               */

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv)             != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)                           != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )                               cmdline_help   (argv[0], go);
  if (esl_opt_ArgNumber(go) < 1)                                   cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        


  /* Open the sequence file */
  seqfile = esl_opt_GetArg(go, 1);
  if (esl_opt_GetString(go, "--informat") != NULL) {
    infmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  /* Open the output file, if any */
  if (esl_opt_GetBoolean(go, "-O")) 
    {
      if ((ofp = fopen(esl_opt_GetArg(go, 2), "w")) == NULL)
	esl_fatal("Failed to open output file %s\n", esl_opt_GetArg(go, 2));
    }
  else if (esl_opt_GetString(go, "-o") != NULL)
    {
      if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
	esl_fatal("Failed to open output file %s\n", esl_opt_GetString(go, "-o"));
    }
  else ofp = stdout;

  if (esl_opt_GetBoolean(go, "--index")) 
    {
      if (esl_opt_ArgNumber(go) != 1) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      create_ssi_index(go, sqfp);
    }
  else if (esl_opt_GetBoolean(go, "-f"))
    {
      if (esl_opt_ArgNumber(go) != 2) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      multifetch(go, ofp, esl_opt_GetArg(go, 2), sqfp);
    }
  else 
    {
      if (esl_opt_ArgNumber(go) != 2) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      onefetch(go, ofp, esl_opt_GetArg(go, 2), sqfp);
      if (ofp != stdout) printf("\n\nRetrieved sequence %s.\n",  esl_opt_GetArg(go, 2));
    }

  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);
  return 0;
}


/* Create an SSI index file for open sequence file <sqfp>.
 * Both name and accession of sequences are stored as keys.
 */
static void
create_ssi_index(ESL_GETOPTS *go, ESL_SQFILE *sqfp)
{
  int         status;
  ESL_NEWSSI *ns      = esl_newssi_Create();
  ESL_SQ     *sq      = NULL;
  int         nali    = 0;
  char       *ssifile = NULL;
  FILE       *sfp     = NULL;
  uint16_t    fh;

  if (afp->ssi != NULL) 
    esl_fatal("Alignment file %s already has an SSI index. Delete or move it first.\n", afp->fname);
  if (esl_FileNewSuffix(afp->fname, "ssi", &ssifile) != eslOK)
    esl_fatal("Failed to name SSI file for %s\n", afp->fname);
  if ((sfp = fopen(ssifile, "wb")) == NULL)
    esl_fatal("Failed to open SSI file %s\n", ssifile);
  if (esl_newssi_AddFile(ns, afp->fname, afp->format, &fh) != eslOK)
    esl_fatal("Failed to add MSA file %s to new SSI index\n", afp->fname);

  printf("Working...    "); 
  fflush(stdout);
  
  while ((status = esl_msa_Read(afp, &msa)) != eslEOF)
    {
      if (status == eslEFORMAT)
	esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", afp->linenumber, afp->fname, afp->errbuf, afp->buf);
      else if (status != eslOK)
	esl_fatal("Alignment file read failed with error code %d\n", status);
      nali++;

      if (msa->name == NULL) 
	esl_fatal("Every alignment in file must have a name to be indexed. Failed to find name of alignment #%d\n", nali);

      if (esl_newssi_AddKey(ns, msa->name, fh, msa->offset, 0, 0) != eslOK)
	esl_fatal("Failed to add key %s to SSI index", msa->name);

      if (msa->acc != NULL) {
	if (esl_newssi_AddAlias(ns, msa->acc, msa->name) != eslOK)
	  esl_fatal("Failed to add secondary key %s to SSI index", msa->acc);
      }
      esl_msa_Destroy(msa);
    }
  
  if (esl_newssi_Write(sfp, ns) != eslOK) 
    esl_fatal("Failed to write keys to ssi file %s\n", ssifile);

  printf("done.\n");
  if (ns->nsecondary > 0) 
    printf("Indexed %d alignments (%ld names and %ld accessions).\n", nali, (long) ns->nprimary, (long) ns->nsecondary);
  else 
    printf("Indexed %d alignments (%ld names).\n", nali, (long) ns->nprimary);
  printf("SSI index written to file %s\n", ssifile);

  fclose(sfp);
  free(ssifile);
  esl_newssi_Destroy(ns);
  return;
}  
