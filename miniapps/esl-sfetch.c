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

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	                        /* application configuration       */
  char         *seqfile = NULL;	                        /* sequence file name              */
  int           infmt   = eslSQFILE_UNKNOWN;		/* format code for seqfile         */
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
  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) cmdline_failure(argv[0], "Sequence file %s not found.\n",     seqfile);
  else if (status == eslEFORMAT)   cmdline_failure(argv[0], "Format of file %s unrecognized.\n", seqfile);
  else if (status == eslEINVAL)    cmdline_failure(argv[0], "Can't autodetect stdin or .gz.\n");
  else if (status != eslOK)        cmdline_failure(argv[0], "Open failed, code %d.\n", status);

  /* Open the output file, if any */
  if (esl_opt_GetBoolean(go, "-O")) 
    {
      if ((ofp = fopen(esl_opt_GetArg(go, 2), "w")) == NULL)
	cmdline_failure(argv[0], "Failed to open output file %s\n", esl_opt_GetArg(go, 2));
    }
  else if (esl_opt_GetString(go, "-o") != NULL)
    {
      if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
	cmdline_failure(argv[0], "Failed to open output file %s\n", esl_opt_GetString(go, "-o"));
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
      if (! esl_sqio_IsAlignment(sqfp->format)) esl_sqfile_OpenSSI(sqfp, NULL);
      multifetch(go, ofp, esl_opt_GetArg(go, 2), sqfp);
    }
  else 
    {
      if (esl_opt_ArgNumber(go) != 2) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      if (! esl_sqio_IsAlignment(sqfp->format)) esl_sqfile_OpenSSI(sqfp, NULL);
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
  ESL_SQ     *sq      = esl_sq_Create();
  int         nseq    = 0;
  char       *ssifile = NULL;
  FILE       *sfp     = NULL;
  uint16_t    fh;

  esl_strdup(sqfp->filename, -1, &ssifile);
  esl_strcat(&ssifile, -1, ".ssi", 4);
  if (esl_FileExists(ssifile)                            == TRUE)   esl_fatal("SSI file %s already exists; delete or rename", ssifile);
  if ((sfp = fopen(ssifile, "wb"))                       == NULL)   esl_fatal("Failed to open SSI file %s for writing\n",     ssifile);

  if (esl_newssi_AddFile(ns, sqfp->filename, sqfp->format, &fh) != eslOK)
    esl_fatal("Failed to add sequence file %s to new SSI index\n", sqfp->filename);

  printf("Working...    "); 
  fflush(stdout);
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      nseq++;
      if (sq->name == NULL) esl_fatal("Every sequence must have a name to be indexed. Failed to find name of seq #%d\n", nseq);

      if (esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->n) != eslOK)
	esl_fatal("Failed to add key %s to SSI index", sq->name);

      if (sq->acc[0] != '\0') {
	if (esl_newssi_AddAlias(ns, sq->acc, sq->name) != eslOK)
	  esl_fatal("Failed to add secondary key %s to SSI index", sq->acc);
      }
      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal("Parse failed, line %d, file %s: %s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  
  /* Determine if the file was suitable for fast subseq lookup. */
  if (sqfp->bpl > 0 && sqfp->rpl > 0) {
    if ((status = esl_newssi_SetSubseq(ns, fh, sqfp->bpl, sqfp->rpl)) != eslOK) 
      esl_fatal("Failed to set %s for fast subseq lookup.");
  }

  /* Save the SSI file to disk */
  if (esl_newssi_Write(sfp, ns) != eslOK)  esl_fatal("Failed to write keys to ssi file %s\n", ssifile);

  /* Done - output and exit. */
  printf("done.\n");
  if (ns->nsecondary > 0) 
    printf("Indexed %d sequences (%ld names and %ld accessions).\n", nseq, (long) ns->nprimary, (long) ns->nsecondary);
  else 
    printf("Indexed %d sequences (%ld names).\n", nseq, (long) ns->nprimary);
  printf("SSI index written to file %s\n", ssifile);

  fclose(sfp);
  free(ssifile);
  esl_sq_Destroy(sq);
  esl_newssi_Destroy(ns);
  return;
}

/* multifetch:
 * given a file containing lines with one name or key per line;
 * parse the file line-by-line;
 * if we have an SSI index available, retrieve the seqs by key
 * as we see each line;
 * else, without an SSI index, store the keys in a hash, then
 * read the entire seq file in a single pass, outputting seqs
 * that are in our keylist. 
 * 
 * Note that with an SSI index, you get the seqs in the order they
 * appear in the <keyfile>, but without an SSI index, you get seqs in
 * the order they occur in the seq file.
 */
static void
multifetch(ESL_GETOPTS *go, FILE *ofp, char *keyfile, ESL_SQFILE *sqfp)
{
  ESL_KEYHASH    *keys   = esl_keyhash_Create();
  ESL_FILEPARSER *efp    = NULL;
  int             nseq   = 0;
  int             nkeys  = 0;
  char           *key;
  int             keylen;
  int             keyidx;
  int             status;

  
  if (esl_fileparser_Open(keyfile, &efp) != eslOK)  esl_fatal("Failed to open key file %s\n", keyfile);
  esl_fileparser_SetCommentChar(efp, '#');

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &key, &keylen) != eslOK)
	esl_fatal("Failed to read seq name on line %d of file %s\n", efp->linenumber, keyfile);
      
      status = esl_key_Store(keys, key, &keyidx);
      if (status == eslEDUP) esl_fatal("seq key %s occurs more than once in file %s\n", key, keyfile);
	
      /* if we have an SSI index, just fetch them as we go. */
      if (sqfp->ssi != NULL) { onefetch(go, ofp, key, sqfp);  nseq++; }
      nkeys++;
    }

  /* If we don't have an SSI index, we haven't fetched anything yet; do it now. */
  if (sqfp->ssi == NULL) 
    {
      ESL_SQ *sq     = esl_sq_Create();

      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  if ( (sq->name[0] != '\0' && esl_key_Lookup(keys, sq->name, NULL) == eslOK) ||
	       (sq->acc[0]  != '\0' && esl_key_Lookup(keys, sq->acc,  NULL) == eslOK))
	    {
	      esl_sqio_Write(ofp, sq, eslSQFILE_FASTA);
	      nseq++;
	    }
	  esl_sq_Reuse(sq);
	}
      if (status != eslEOF) esl_fatal("Parse failed, line %d, file %s: %s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);
      esl_sq_Destroy(sq);
    }
  
  if (nkeys != nseq) esl_fatal("Tried to retrieve %d keys, but only retrieved %d sequences\n", nkeys, nseq);

  if (ofp != stdout) printf("\nRetrieved %d sequences.\n", nseq);

  esl_keyhash_Destroy(keys);
  esl_fileparser_Close(efp);
  return;
}
  


/* onefetch():
 * Given one <key> (a seq name or accession), retrieve the corresponding sequence.
 * In SSI mode, we can do this quickly by positioning the file, then regurgitating
 * every line until the end-of-record marker; we don't even have to parse.
 * Without an SSI index, we have to parse the file sequentially 'til we find
 * the one we're after.
 */
static void
onefetch(ESL_GETOPTS *go, FILE *ofp, char *key, ESL_SQFILE *sqfp)
{
  int status;

  if (sqfp->ssi != NULL)	/* the fast way */
    {
      status = esl_sqfile_PositionByKey(sqfp, key);
      if      (status == eslENOTFOUND) esl_fatal("seq %s not found in SSI index for file %s\n", key, sqfp->filename);
      else if (status == eslEFORMAT)   esl_fatal("Failed to parse SSI index for %s\n", sqfp->filename);
      else if (status != eslOK)        esl_fatal("Failed to look up location of seq %s in SSI index of file %s\n", key, sqfp->filename);
      
      if ((status = esl_sqio_Echo(ofp, sqfp)) != eslOK) esl_fatal("Echo failed: %s\n", sqfp->errbuf);
    }
  else				/* the slow way */
    {
      ESL_SQ  *sq   = esl_sq_Create();

      while ((status = esl_sqio_Read(sqfp, sq)) != eslEOF)
	{
	  if (status != eslOK) esl_fatal("Parse failed, line %d, file %s: %s\n", sqfp->linenumber, sqfp->filename, sqfp->errbuf);

	  if (strcmp(key, sq->name) == 0 || strcmp(key, sq->acc) == 0) 
	    {
	      esl_sqio_Write(ofp, sq, eslSQFILE_FASTA);
	      break;
	    }
	  esl_sq_Reuse(sq);
	}
      esl_sq_Destroy(sq);
    }
}

