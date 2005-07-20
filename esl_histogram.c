/* esl_histogram.c
 * Collecting and displaying histograms.
 * 
 * SRE, Fri Jul  1 13:21:45 2005 [St. Louis]
 * SVN $Id$
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <easel.h>
#include <esl_histogram.h>
#ifdef eslAUGMENT_GUMBEL
#include <esl_gumbel.h>
#endif
#ifdef eslAUGMENT_GEV
#include <esl_gev.h>
#endif
#ifdef eslAUGMENT_STATS
#include <esl_stats.h>
#endif

/* Function:  esl_histogram_Create()
 * Incept:    SRE, Fri Jul  1 13:40:26 2005 [St. Louis]
 *
 * Purpose:   Creates and returns a new histogram object, initially
 *            allocated to count scores $>=$ <xmin> and $<$ <xmax> into
 *            bins of width <w>. Thus, a total of <xmax>-<xmin>/<w> bins
 *            are initially created. 
 *            
 *            The bounds <xmin> and <xmax> only need to be initial
 *            guesses.  The histogram object will reallocate itself
 *            dynamically as needed to accomodate scores that exceed
 *            current bounds.
 *
 *            For example, <esl_histogram_Create(-100, 100, 0.5)> 
 *            would init the object to collect scores into 400 bins 
 *            $[-100<=x<-99.5],[-99.5<=x<-99.0]...[99.5<=x<100.0]$.
 *  
 * Args:      xmin - caller guesses that minimum score will be >= xmin
 *            xmax - caller guesses that max score will be < xmax
 *            w    - size of bins (1.0, for example)
 *            
 * Returns:   ptr to new <ESL_HISTOGRAM> object, which caller is responsible
 *            for free'ing with <esl_histogram_Destroy()>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_HISTOGRAM *
esl_histogram_Create(double xmin, double xmax, double w)
{
  ESL_HISTOGRAM *h = NULL;
  int i;

  if ((h = malloc(sizeof(ESL_HISTOGRAM))) == NULL)    
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  h->bin   = NULL;
  h->nbins = 1 + (int)((xmax-xmin)/w);
  h->w     = w;
  h->xmin  = xmin;
  h->xmax  = xmax;
  h->imin  = h->nbins;
  h->imax  = -1;
  h->total = 0;
  h->expect= NULL;

  if ((h->bin = malloc(sizeof(int) * h->nbins)) == NULL)
    { esl_histogram_Destroy(h); ESL_ERROR_NULL(eslEMEM, "malloc failed");}
  for (i = 0;  i < h->nbins; i++)
    h->bin[i] = 0;

  return h;
}

/* Function:  esl_histogram_Destroy()
 * Incept:    SRE, Sat Jul  2 19:41:17 2005 [St. Louis]
 *
 * Purpose:   Frees an <ESL_HISTOGRAM> object <h>.
 */
void
esl_histogram_Destroy(ESL_HISTOGRAM *h)
{
  if (h ==  NULL) return;
  if (h->bin    != NULL) free(h->bin); 
  if (h->expect != NULL) free(h->expect);
  free(h);
  return;
}


/* Function:  esl_histogram_Add()
 * Incept:    SRE, Sat Jul  2 19:41:45 2005 [St. Louis]
 *
 * Purpose:   Adds score <x> to a histogram <h>.
 *           
 *            The histogram will be automatically reallocated as
 *            needed if the score is smaller or larger than the
 *            current allocated bounds.  
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
esl_histogram_Add(ESL_HISTOGRAM *h, double x)
{
  int i,j;			/* what bin we're in                       */
  int nnew;			/* # of new bins created by a reallocation */

  i = (int) floor((x - h->xmin)/h->w);	
  
  /* Reallocate below?
   */
  if (i < 0) 
    {				
      nnew = -i*2;	/* overallocate by 2x */
      h->bin = realloc(h->bin, sizeof(int) * (nnew+ h->nbins));
      if (h->bin == NULL) { ESL_ERROR(eslEMEM, "reallocation failed"); }
      
      memmove(h->bin+nnew, h->bin, sizeof(int) * h->nbins);
      h->nbins += nnew;
      i        += nnew;
      h->xmin  -= nnew*h->w;
      h->imin  += nnew;
      if (h->imax > -1) h->imax += nnew;
      for (j = 0; j < nnew; j++) h->bin[j] = 0;
    }
  /* Reallocate above?
   */
  else if (i >= h->nbins)
    {
      nnew = (i-h->nbins+1) * 2; /* 2x overalloc */
      
      h->bin = realloc(h->bin, sizeof(int) * (nnew+ h->nbins));
      if (h->bin == NULL) { ESL_ERROR(eslEMEM, "reallocation failed"); }

      for (j = h->nbins; j < h->nbins+nnew; j++)
	h->bin[j] = 0;      
      if (h->imin == h->nbins) h->imin+=nnew;
      h->xmax  += nnew*h->w;
      h->nbins += nnew;
    }

  /* Bump the bin counter, and whatever other data we're keeping.
   */
  h->bin[i]++;
  h->total++;
  if (i > h->imax) h->imax = i;
  if (i < h->imin) h->imin = i;
  return eslOK;
}
  

/* Function:  esl_histogram_Print() 
 * Incept:    SRE, Sat Jul  2 16:03:37 2005 [St. Louis]
 *
 * Purpose:   Print a "prettified" histogram <h> to a file pointer <fp>.
 *            Deliberately a look-and-feel clone of Bill Pearson's 
 *            excellent FASTA output.
 * 
 * Args:      fp     - open file to print to (stdout works)
 *            h      - histogram to print
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_histogram_Print(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;
  int maxbar;
  int imode;
  int units;
  int num;
  char buffer[81];		/* output line buffer */
  int  pos;			/* position in output line buffer */
  int  ilowbound, lowcount;	/* cutoffs on the low side  */
  int  ihighbound, highcount;	/* cutoffs on the high side */
  int  emptybins = 3;

  /* Find out how we'll scale the histogram.  We have 59 characters to
   * play with on a standard 80-column terminal display: leading "%5d
   * %6d %6d|" occupies 20 chars.  Save the peak position, we'll use
   * it later.
   */
  maxbar = 0;
  for (i = 0; i < h->nbins; i++)
    if (h->bin[i] > maxbar) 
      {
	maxbar  = h->bin[i];     /* max height    */
	imode   = i;
      }

  /* Truncate histogram display on both sides, ad hoc fashion.
   * Start from the peak; then move out until we see <emptybins> empty bins,
   * and stop.
   */
  for (num = 0, ihighbound = imode; ihighbound < h->imax; ihighbound++)
    {
      if (h->bin[ihighbound] > 0) { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }
  for (num = 0, ilowbound = imode; ilowbound > h->imin; ilowbound--)
    {
      if (h->bin[ilowbound] > 0)  { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }

				/* collect counts outside of bounds */
  for (lowcount = 0, i = h->imin; i < ilowbound; i++)
    lowcount += h->bin[i];
  for (highcount = 0, i = h->imax; i > ihighbound; i--)
    highcount += h->bin[i];

				/* maxbar might need to be raised now; then set our units  */
  if (lowcount  > maxbar) maxbar = lowcount;
  if (highcount > maxbar) maxbar = highcount;
  units = ((maxbar-1)/ 59) + 1;

  /* Print the histogram
   */
  fprintf(fp, "%5s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  fprintf(fp, "%5s %6s %6s\n", "-----", "---", "---");
  buffer[80] = '\0';
  buffer[79] = '\n';
  for (i = h->imin; i <= h->imax; i++)
    {
      memset(buffer, ' ', 79 * sizeof(char));
      x = i*h->w + h->xmin;

      /* Deal with special cases at edges
       */
      if      (i < ilowbound)  continue;
      else if (i > ihighbound) continue;
      else if (i == ilowbound && i != h->imin) 
	{
	  sprintf(buffer, "<%4.0f %6d %6s|", x+h->w, lowcount, "-");
	  if (lowcount > 0) {
	    num = 1+(lowcount-1) / units;
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}
      else if (i == ihighbound && i != h->imax)
	{
	  sprintf(buffer, ">%4.0f %6d %6s|", x, highcount, "-");
	  if (highcount > 0) {
	    num = 1+(highcount-1) / units;
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}

      /* Deal with most cases
       */
      if (h->expect != NULL) 
	sprintf(buffer, "%5.0f %6d %6d|", 
		x, h->bin[i], (int) h->expect[i]);
      else
	sprintf(buffer, "%5.0f %6d %6s|", x, h->bin[i], "-");
      buffer[20] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->bin[i] > 0) {
	num = 1 + (h->bin[i]-1) / units;
	for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       */
      if (h->expect != NULL && (int) h->expect[i] > 0)
	{
	  pos = 20 + (int)(h->expect[i]-1) / units;
	  if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	  buffer[pos] = '*';
	}

      /* Print the line
       */
      fputs(buffer, fp);
    }

  return eslOK;
}
  

/* Function:  esl_histogram_Plot()
 * Incept:    SRE, Sat Jul  2 19:43:20 2005 [St. Louis]
 *
 * Purpose:   Print a histogram <h> to open file ptr <fp>, 
 *            as an XY file suitable for input to the xmgrace
 *            graphing program.
 *
 */
void
esl_histogram_Plot(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;

  /* First data set is the observed histogram
   */
  for (i = h->imin; i <= h->imax; i++)
    if (h->bin[i] > 0)
      {
	x = i*h->w + h->xmin;
	fprintf(fp, "%f %d\n", x, h->bin[i]);
      }
  fprintf(fp, "&\n");

  /* Second data set is the theoretical histogram
   */
  if (h->expect != NULL)
    {
      for (i = 0; i < h->nbins; i++)
	if (h->expect[i] > 0.)
	  {
	    x = i*h->w + h->xmin;
	    fprintf(fp, "%.2f %.2f\n", x, h->expect[i]);
	  }
      fprintf(fp, "&\n");
    }
}

/*****************************************************************
 * Functions for setting expected histogram frequencies
 *****************************************************************/
#ifdef eslAUGMENT_GUMBEL
int
esl_histogram_SetGumbel(ESL_HISTOGRAM *h, double mu, double lambda)
{
  int    i;
  double x1,x2;

  if (h->expect == NULL) 
    {
      h->expect = malloc(sizeof(double) * h->nbins);
      if (h->expect == NULL) ESL_ERROR(eslEMEM, "malloc failed");
    }

  for (i = 0; i < h->nbins; i++)
    {
      x1 = i * h->w + h->xmin;	   /* scores in bin i are > x1 */
      x2 = (i+1) * h->w + h->xmin; /*     ...and <= x2.        */

      h->expect[i] = h->total *
	(esl_gumbel_cdf(x2, mu, lambda) - esl_gumbel_cdf(x1, mu, lambda));
    }

  return eslOK;
}
#endif /*eslAUGMENT_GUMBEL*/

#ifdef eslAUGMENT_GEV
int
esl_histogram_SetGEV(ESL_HISTOGRAM *h, double mu, double lambda, double alpha)
{
  int    i;
  double x1,x2;

  if (h->expect == NULL) 
    {
      h->expect = malloc(sizeof(double) * h->nbins);
      if (h->expect == NULL) ESL_ERROR(eslEMEM, "malloc failed");
    }

  for (i = 0; i < h->nbins; i++)
    {
      x1 = i * h->w + h->xmin;	   /* scores in bin i are > x1 */
      x2 = (i+1) * h->w + h->xmin; /*     ...and <= x2.        */

      h->expect[i] = h->total *
	(esl_gev_cdf(x2, mu, lambda, alpha) - esl_gev_cdf(x1, mu, lambda, alpha));
    }
  return eslOK;
}
#endif /*eslAUGMENT_GEV*/


#ifdef eslAUGMENT_STATS
int
esl_histogram_GTestGoodness(ESL_HISTOGRAM *h, int ndeg, 
			    int *ret_nbins, double *ret_G, double *ret_p)
{
  int    i;
  int    nbins;			/* number of bins counted toward G */
  double G = 0.;		/* the G-statistic */
  double Gp;			/* P(test > G) by chi-square distribution */
  int    status;
  
  /* Calculate the G statistic = 2 * log likelihood ratio.
   */
  nbins = 0;
  for (i = 0; i < h->nbins; i++)
    if (h->bin[i] > 0)
      {
	G += (double) h->bin[i] * log((double) h->bin[i]/ h->expect[i]);
	nbins++;
      }
  G *= 2.;
  
  /* G is distributed approximately as \chi^2
   */
  if (nbins-1-ndeg >= 0) 
    {
      status = esl_stats_ChiSquaredTest(nbins-1-ndeg, G, &Gp);
      if (status != eslOK) return status;
    }
  else Gp = 0.;

  if (ret_nbins != NULL) *ret_nbins = nbins;
  if (ret_G     != NULL) *ret_G     = G;
  if (ret_p     != NULL) *ret_p     = Gp;
  return eslOK;
}

int
esl_histogram_ChiSquaredGoodness(ESL_HISTOGRAM *h, int ndeg, 
				 int *ret_nbins, double *ret_X2, double *ret_p)
{
  int    i;
  int    nbins;			/* number of bins counted toward X^2 */
  double chisq = 0.;		/* the X^2 statistic */
  double delta;			/* obs - exp in a bin */
  double chip;		        /* P(test > X^2) by chi-square distribution */
  int    status;
  
  /* Calculate the X^2 statistic = \sum_i (obs_i-exp_i)^2/exp, over
   * all bins containing some minimum size (arbitrarily 5)
   */
  nbins = 0;
  for (i = 0; i < h->nbins; i++)
    if (h->bin[i] >= 5 && h->expect[i] >= 5) 
      {
	delta = h->bin[i] - h->expect[i];
	chisq += delta * delta / h->expect[i];
	nbins++;
      }
  
  /* X^2 is distributed approximately as \chi^2
   */
  if (nbins-1-ndeg >= 0)
    {
      status = esl_stats_ChiSquaredTest(nbins-1-ndeg, chisq, &chip);
      if (status != eslOK) return status;
    }
  else chip = 0.;

  if (ret_nbins != NULL) *ret_nbins = nbins;
  if (ret_X2    != NULL) *ret_X2    = chisq;
  if (ret_p     != NULL) *ret_p     = chip;
  return eslOK;
}


#endif /*eslAUGMENT_STATS*/


/*****************************************************************
 * Example main()
 *****************************************************************/
#ifdef eslHISTOGRAM_EXAMPLE
/*::cexcerpt::histogram_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslHISTOGRAM_EXAMPLE esl_histogram.c esl_random.c easel.c
 * run:     ./example 
 */
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_Create(-100, 100, 1);
  int    nsamples    = 1000;
  double mean        = 20.0;
  double stddev      = 10.0;
  int    i;
  double x;

  for (i = 0; i < nsamples; i++) {
    x = esl_rnd_Gaussian(r, mean, stddev);
    esl_histogram_Add(h, x);
  }

  esl_histogram_Print(stdout, h);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example::end::*/
#endif /*eslHISTOGRAM_EXAMPLE*/
