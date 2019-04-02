#include "esl_config.h"

#include <string.h>

#include "easel.h"
#include "esl_subcmd.h"

extern int esl_cmd_alistat   (const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);  // cmd_alistat.c
extern int esl_cmd_downsample(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);  // cmd_downsample.c
extern int esl_cmd_filter    (const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);  // cmd_filter.c

ESL_SUBCMD subcommands[] = {
  /* function            subcmd_name  nargs        arg_description               help_line */
  { esl_cmd_alistat,    "alistat",       1, "[-options] <msafile>",         "summary statistics for a multiple seq alignment file"     },
  { esl_cmd_downsample, "downsample",    2, "[-options] <m> <infile>",      "downsample <m> things from larger <infile> of n things"   },
  { esl_cmd_filter,     "filter",        2, "[-options] <maxid> <msafile>", "remove seqs >= <maxid> fractional identity from MSA"      },
};

static ESL_OPTIONS top_options[] = {
   /* name         type          default  env  range tog's   reqs incomp  help                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show overall brief help summary", 1  },
  { "--version",  eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show version number",             1  },
  { "--help",     eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show overall brief help summary", 99 },  // accept --help as an undocumented special case
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static int
top_usage(const char *topcmd)
{
  char *lastslash = strrchr(topcmd, '/');
  if (lastslash) topcmd = lastslash+1;

  if (printf("Usage:\n")                                                                < 0) esl_fatal("printf failed");
  if (printf("  %s -h                : show overall brief help summary\n",      topcmd) < 0) esl_fatal("printf failed");
  if (printf("  %s --version         : show version number\n",                  topcmd) < 0) esl_fatal("printf failed");
  if (printf("  %s <cmd> -h          : show brief help for an Easel command\n", topcmd) < 0) esl_fatal("printf failed");
  if (printf("  %s <cmd> [<args>...] : run an Easel command\n",                 topcmd) < 0) esl_fatal("printf failed");
  return eslOK;
}

static int
top_help(const char *topcmd)
{
  int   ncmds     =  sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int   i;
  int   status;

  if ( printf("easel: little utilities for biological sequence analysis\n")        < 0) esl_fatal("printf failed");
  if ( printf("version %s (%s): %s\n\n", EASEL_VERSION, EASEL_DATE, EASEL_URL)     < 0) esl_fatal("printf failed");
  if (( status = top_usage(topcmd)) != eslOK) return status;
  if ( printf("\navailable commands:\n")                                           < 0) esl_fatal("printf failed");
  for (i = 0; i < ncmds; i++)
    if ( printf("  %-12s %s\n", subcommands[i].subcmd, subcommands[i].description) < 0) esl_fatal("printf failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(top_options);
  int ncmds = sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int idx;
  int status;
 
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  
  if (esl_opt_GetBoolean(go, "--version") == TRUE) { printf("%s\n", EASEL_VERSION); status = eslOK; goto DONE; }
  if (esl_opt_GetBoolean(go, "--help")    == TRUE) { status = top_help(argv[0]);    goto DONE; }
  if (esl_opt_GetBoolean(go, "-h")        == TRUE) { status = top_help(argv[0]);    goto DONE; }
  if (argc - go->optind == 0)                      { status = top_help(argv[0]);    goto DONE; }

  for (idx = 0; idx < ncmds; idx++)
    if (strcmp(go->argv[go->optind], subcommands[idx].subcmd) == 0) break;
  if (idx == ncmds) { status = top_usage(argv[0]); goto DONE; }

  status = subcommands[idx].func(argv[0], &subcommands[idx], argc-go->optind, argv+go->optind);
  
 DONE:
  esl_getopts_Destroy(go);
  return status;
}

