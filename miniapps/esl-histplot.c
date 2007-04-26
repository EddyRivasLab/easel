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
#include <esl_exponential.h>
#include <esl_gumbel.h>


static char banner[] = "\
esl-histplot :: collate data histogram, output xmgrace datafile";

static char usage[] = "\
Usage: esl-histplot [-options] <datafile>";

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                     1 },
  {"-o",  eslARG_OUTFILE, NULL,  NULL, NULL, NULL, NULL, NULL, "output file for histogram",               1},

  {"-b",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "input file is binary, array of doubles",  2 },
  {"-f",  eslARG_INT,       "1", NULL, "n>0",NULL, NULL, "-b", "which field to read on text line (1..n)", 2 },

  {"-w",   eslARG_REAL,   "1.0", NULL, NULL, NULL, NULL, NULL, "bin size for histogram",                  3 },
  {"--min", eslARG_REAL, "-100.",NULL, NULL, NULL, NULL, NULL, "initial lower bound of histogram",        3 },
  {"--max", eslARG_REAL,  "100.",NULL, NULL, NULL, NULL, NULL, "initial upper bound of histogram",        3 },
  {"--surv", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "output survival plot, not histogram",     3 },
   
  {"--gumbel",        eslARG_NONE,FALSE,NULL, NULL, NULL, NULL, NULL,   "fit data to a Gumbel distribution",                   4 }, 
  {"--gumbel-lambda", eslARG_REAL,"0.693", NULL, NULL, NULL, NULL, NULL, "fit data to a Gumbel distribution w/ known lambda",   4 }, 
  {"--exp-tail",      eslARG_REAL, "0.01", NULL, NULL, NULL, NULL, NULL,  "fit tail to an exponential distribution",             4 },

  { 0,0,0,0,0,0,0,0,0,0},
};


int
main(int argc, char **argv)
{
  ESL_HISTOGRAM  *h;		/* full histogram w/ collated values          */
  ESL_GETOPTS    *go;		/* application configuration                  */
  char           *datafile;	/* input data; "-" means stdin                */
  FILE           *ifp;		/* input stream                               */
  double          x;		/* value of field, after conversion to double */

  int             which_field;	/* which field to use as data, 1..nf (default 1)*/
  char           *outfile;	/* output xmgrace xy data file                  */
  FILE           *ofp;		/* output data stream                           */
  double          hmin;		/* initial histogram lower bound                */
  double          hmax;		/* initial histogram upper bound                */
  double          hbinsize;	/* histogram's bin size                         */
  int             is_binary;	/* datafile is a binary file of doubles         */

  double *xv;
  int     n;
  double  params[2];

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\nGeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 0= group; 2 = indentation; 80=textwidth*/
    puts("\nOptions that control how to read the input file:");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
    puts("\nOptions that control how to display the output XY file:");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 2= group; 2 = indentation; 80=textwidth*/
    puts("\nOptional ML fitting of the data to distributions:");
    esl_opt_DisplayHelp(stdout, go, 4, 2, 80); /* 3= group; 2 = indentation; 80=textwidth*/
    return eslOK;
  }
  esl_opt_GetIntegerOption(go, "-f",         &which_field);
  esl_opt_GetStringOption (go, "-o",         &outfile);
  esl_opt_GetDoubleOption (go, "-w",         &hbinsize);
  esl_opt_GetDoubleOption (go, "--min",      &hmin);
  esl_opt_GetDoubleOption (go, "--max",      &hmax);

  if (esl_opt_ArgNumber(go) != 1) 
    esl_fatal("Incorrect number of command line arguments.\n%s\n", usage); 

  datafile = esl_opt_GetCmdlineArg(go, eslARG_INFILE, NULL);


  /*****************************************************************
   * Open the input and output datafiles, and init the histogram.
   *****************************************************************/

  if (strcmp(datafile, "-") == 0) ifp = stdin;
  else {
    ifp = fopen(datafile, "rb");
    if (ifp == NULL) esl_fatal("Failed to open input data file %s\n", datafile);
  }

  if (outfile == NULL) ofp = stdout;
  else {
    ofp = fopen(outfile, "w");
    if (ofp == NULL) esl_fatal("Failed to open output xmgrace data file %s\n", outfile);
  }

  h = esl_histogram_CreateFull(hmin,hmax,hbinsize);
  if (h == NULL) esl_fatal("Failed to create histogram");


  /*****************************************************************
   * Collect the data
   *****************************************************************/

  if (esl_opt_IsSet(go, "-b"))
    {
      while (fread(&x, sizeof(double), 1, ifp) == 1)
	esl_histogram_Add(h, x);
    } 
  else 
    {
      char           *buf;		/* ptr to line buffer, for esl_fgets()        */
      int             nbuf;		/* allocated line lengths, for esl_fgets()    */
      char           *s;		/* ptr to line, for esl_strtok()              */
      char           *tok;		/* ptr to a data field on the line            */
      int             i;		/* counter over fields                        */
      int             toklen;	/* string length of a field                   */

      buf  = NULL;
      nbuf = 0;
      while (esl_fgets(&buf, &nbuf, ifp) == eslOK)
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
      free(buf);
    }
  
  /*****************************************************************
   * Optionally, fit the data
   *****************************************************************/

  if (esl_opt_IsSet(go, "--gumbel"))
    {
      esl_histogram_GetData(h, &xv, &n);
      esl_gumbel_FitComplete(xv, n, &(params[0]), &(params[1]));
      esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &params);
    }
  else if (esl_opt_IsSet(go, "--gumbel-lambda"))
    {
      esl_opt_GetDoubleOption(go, "--gumbel-lambda", &(params[1]));
      esl_histogram_GetData(h, &xv, &n);
      esl_gumbel_FitCompleteLoc(xv, n, params[1], &(params[0]));
      esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &params);
    }
  else if (esl_opt_IsSet(go, "--exp-tail"))
    {
      double tailp;
      esl_opt_GetDoubleOption(go, "--exp-tail", &tailp);
      esl_histogram_GetTailByMass(h, tailp, &xv, &n, NULL);
      esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
      esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);
    }
      
  /*****************************************************************
   * Output
   *****************************************************************/

  if   (esl_opt_IsSet(go, "--surv"))  esl_histogram_PlotSurvival(ofp, h);
  else                                esl_histogram_Plot(ofp, h);


  /*****************************************************************
   * Cleanup
   *****************************************************************/

  if (outfile != NULL)            fclose(ofp);
  if (strcmp(datafile, "-") != 0) fclose(ifp);
  esl_histogram_Destroy(h);
  esl_getopts_Destroy(go);
  return 0;
}
