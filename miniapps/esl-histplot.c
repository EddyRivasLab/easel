/* esl-histplot  - collate data into a histogram and output xmgrace datafile.
 *                  
 * SRE, Tue Feb 21 14:18:05 2006                  
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>

#include <easel.h>
#include <esl_getopts.h>
#include <esl_histogram.h>

static char banner[] = "\
esl-histplot :: collate data histogram, output xmgrace datafile";

static char usage[] = "\
Usage: esl-histplot [-options] <datafile>\n\
  Available options are:\n\
  -h     : help; print brief info on version and usage\n\
  -f <n> : use field <n> as data, 1..N (default=1, first field)\n\
  -o <f> : output xmgrace xy datafile to file <f>\n\
";

static char experts[] = "\
  Expert options:\n\
  [none]\n\
";

static ESL_OPTIONS options[] = {
   /* name          type        default env   range  togs  reqs  incompat */
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL, NULL },
  { "-f",         eslARG_INT,      "1", NULL, "n>0", NULL, NULL, NULL },
  { "-o",         eslARG_STRING,  NULL, NULL, NULL,  NULL, NULL, NULL },
  { 0,0,0,0,0,0,0,0 },
};


int
main(int argc, char **argv)
{
  ESL_HISTOGRAM  *h;		/* full histogram w/ collated values          */
  ESL_GETOPTS    *go;		/* application configuration                  */
  char           *datafile;	/* input data; "-" means stdin                */
  FILE           *ifp;		/* input stream                               */
  char           *buf;		/* ptr to line buffer, for esl_fgets()        */
  int             nbuf;		/* allocated line lengths, for esl_fgets()    */
  char           *s;		/* ptr to line, for esl_strtok()              */
  char           *tok;		/* ptr to a data field on the line            */
  int             i;		/* counter over fields                        */
  int             toklen;	/* string length of a field                   */
  double          x;		/* value of field, after conversion to double */

  int             show_help;	/* TRUE to show usage and exit                  */
  int             which_field;	/* which field to use as data, 1..nf (default 1)*/
  char           *outfile;	/* output xmgrace xy data file                  */
  FILE           *ofp;		/* output data stream                           */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  esl_opt_GetBooleanOption(go, "-h",         &show_help);
  esl_opt_GetIntegerOption(go, "-f",         &which_field);
  esl_opt_GetStringOption (go, "-o",         &outfile);

  if (show_help) 
    {
      esl_banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    esl_fatal("Incorrect number of command line arguments.\n%s\n", usage); 


  datafile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);
  esl_getopts_Destroy(go);

  /*****************************************************************
   * Open the input and output datafiles, and init the histogram.
   *****************************************************************/

  if (strcmp(datafile, "-") == 0)
    {
      ifp = stdin;
    }
  else
    {
      ifp = fopen(datafile, "r");
      if (ifp == NULL)
	esl_fatal("Failed to open input data file %s\n", datafile);
    }

  if (outfile == NULL)
    {
      ofp = stdout;
    }
  else
    {
      ofp = fopen(outfile, "w");
      if (ofp == NULL)
	esl_fatal("Failed to open output xmgrace data file %s\n", outfile);
    }

  h = esl_histogram_CreateFull(-100,100,1.0);
  if (h == NULL) esl_fatal("Failed to create histogram");


  /*****************************************************************
   * Collect the data
   *****************************************************************/

  buf  = NULL;
  nbuf = 0;
  while (esl_fgets(&buf, &n, fp) == eslOK)
    {
      s = buf;
      for (i = 0; i < which_field; i++)
	{
	  esl_strtok(&s, " \t\n", &tok, &toklen);
	  if (tok == NULL) continue;
	}

      x = atof(tok);
      esl_histogram_Add(h, x);
    }
  
  /*****************************************************************
   * Output
   *****************************************************************/

  esl_histogram_PlotSurvival(ofp, h);

  if (outfile != NULL) fclose(ofp);
  if (strcmp(datafile, "-") != 0) fclose(ifp);
  esl_histogram_Destroy(h);
  if (buf != NULL) free(buf);
  return 0;
}
