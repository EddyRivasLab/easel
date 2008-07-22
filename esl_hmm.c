/* General hidden Markov models (discrete, of alphabetic strings
 * 
 * SRE, Fri Jul 18 09:00:14 2008 [Janelia]
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_hmm.h"
#include "esl_random.h"
#include "esl_vectorops.h"

/* Function:  esl_hmm_Create()
 * Synopsis:  Creates a new HMM.
 * Incept:    SRE, Fri Jul 18 09:01:54 2008 [Janelia]
 *
 * Purpose:   Creates a new HMM of <M> states for
 *            generating or modeling strings in the
 *            alphabet <abc>.
 *
 * Returns:   a pointer to the new HMM.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_HMM *
esl_hmm_Create(ESL_ALPHABET *abc, int M)
{
  ESL_HMM *hmm = NULL;
  int      k;
  int      status;

  ESL_ALLOC(hmm, sizeof(ESL_HMM));
  hmm->t = NULL;
  hmm->e = NULL;

  ESL_ALLOC(hmm->t, sizeof(float *) * M);
  ESL_ALLOC(hmm->e, sizeof(float *) * M);
  hmm->t[0] = NULL;
  hmm->e[0] = NULL;

  ESL_ALLOC(hmm->t[0], sizeof(float) * M * (M+1));  /* state M is the implicit end state */
  ESL_ALLOC(hmm->e[0], sizeof(float) * M * abc->K);
  ESL_ALLOC(hmm->pi,   sizeof(float) * (M+1));      /* initial transition to state M means a L=0 sequence */
  
  for (k = 1; k < M; k++)
    {
      hmm->t[k] = hmm->t[0] + k*(M+1);
      hmm->e[k] = hmm->e[0] + k*(abc->K);
    }
  
  hmm->M   = M;
  hmm->K   = abc->K;
  hmm->abc = abc;
  return hmm;

 ERROR:
  esl_hmm_Destroy(hmm);
  return NULL;
}

int
esl_hmm_SetDegeneracies(ESL_HMM *hmm)
{
  int Kp = hmm->abc->Kp;
  int K  = hmm->abc->K;
  int k, x,y;

  for (k = 0; k < hmm->M; k++)
    {
      hmm->e[k][K]    = 0.0;	/* gap char                      */
      hmm->e[k][Kp-1] = 1.0;	/* missing data (treated as N/X) */
      for (x = K+1; x <= Kp-2; x++) /* other degeneracies are summed over their residues */
	{
	  hmm->e[k][x] = 0.0;
	  for (y = 0; y < K; y++)
	    if (hmm->abc->degen[x][y]) hmm->e[k][x] += hmm->e[k][y];
	}
    }
  return eslOK;
}


/* Function:  esl_hmm_Destroy()
 * Synopsis:  Destroys an HMM.
 * Incept:    SRE, Fri Jul 18 09:06:22 2008 [Janelia]
 *
 * Purpose:   Frees an HMM.
 */
void
esl_hmm_Destroy(ESL_HMM *hmm)
{
  if (hmm == NULL) return;

  if (hmm->t != NULL) {
    if (hmm->t[0] != NULL) free(hmm->t[0]);
    free(hmm->t);
  }
  if (hmm->e != NULL) {
    if (hmm->e[0] != NULL) free(hmm->e[0]);
    free(hmm->e);
  }
  if (hmm->pi != NULL) free(hmm->pi);
  free(hmm);
  return;
}


ESL_HMX *
esl_hmx_Create(int allocL, int allocM)
{
  ESL_HMX *mx = NULL;
  int      i;
  int      status;
  
  ESL_ALLOC(mx, sizeof(ESL_HMX));
  mx->dp_mem = NULL;
  mx->dp     = NULL;
  mx->sc     = NULL;

  ESL_ALLOC(mx->dp_mem, sizeof(float) * (allocL+1) * allocM);
  mx->ncells = (allocL+1) * allocM;
  
  ESL_ALLOC(mx->dp, sizeof (float *) * (allocL+1));
  ESL_ALLOC(mx->sc, sizeof (float)   * (allocL+2));
  mx->allocL = allocL;

  for (i = 0; i <= allocL; i++)
    mx->dp[i] = mx->dp_mem + i*allocM;
  mx->validL = allocL;
  mx->allocM = allocM;

  mx->L = 0;
  mx->M = 0;
  return mx;

 ERROR:
  esl_hmx_Destroy(mx);
  return NULL;
}

void
esl_hmx_Destroy(ESL_HMX *mx)
{
  if (mx == NULL) return;

  if (mx->dp_mem != NULL) free(mx->dp_mem);
  if (mx->dp     != NULL) free(mx->dp);
  if (mx->sc     != NULL) free(mx->sc);
  free(mx);
  return;
}


/* Function:  esl_hmm_Emit()
 * Synopsis:  Emit a sequence from an HMM.
 * Incept:    SRE, Fri Jul 18 13:16:20 2008 [Janelia]
 *
 * Purpose:   Sample one sequence from an <hmm>, using random
 *            number generator <r>. Optionally return the sequence,
 *            the state path, and/or the length via <opt_dsq>, 
 *            <opt_path>, and <opt_L>.
 *            
 *            If <opt_dsq> or <opt_path> are requested, caller
 *            becomes responsible for free'ing their memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_hmm_Emit(ESL_RANDOMNESS *r, const ESL_HMM *hmm, ESL_DSQ **opt_dsq, int **opt_path, int *opt_L)
{
  int      k, L, allocL;
  int     *path = NULL;
  ESL_DSQ *dsq  = NULL;
  void    *tmp  = NULL;
  int      status;
  
  ESL_ALLOC(dsq,  sizeof(ESL_DSQ) * 256);
  ESL_ALLOC(path, sizeof(int)     * 256);
  allocL = 256;

  dsq[0]  = eslDSQ_SENTINEL;
  path[0] = -1;
  
  k = esl_rnd_FChoose(r, hmm->pi, hmm->M+1);
  L = 0;
  while (k != hmm->M)		/* M is the implicit end state */
    {
      L++;
      if (L >= allocL-1) {	/* Reallocate path and seq if needed */
	ESL_RALLOC(dsq,  tmp, sizeof(ESL_DSQ) * (allocL*2));
	ESL_RALLOC(path, tmp, sizeof(int)     * (allocL*2));
	allocL *= 2;
      }
	
      path[L] = k;
      dsq[L]  = esl_rnd_FChoose(r, hmm->e[k], hmm->abc->K);
      k       = esl_rnd_FChoose(r, hmm->t[k], hmm->M+1);
    }

  path[L+1] = hmm->M;		/* sentinel for "end state" */
  dsq[L+1]  = eslDSQ_SENTINEL;
  
  if (opt_dsq  != NULL) *opt_dsq  = dsq;   else free(dsq);
  if (opt_path != NULL) *opt_path = path;  else free(dsq);
  if (opt_L    != NULL) *opt_L    = L;     
  return eslOK;

 ERROR:
  if (path != NULL) free(path);
  if (dsq  != NULL) free(dsq);
  return status;
}


int
esl_hmm_Forward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *fwd, float *opt_sc)
{
  int   i, k, m;
  int   M     = hmm->M;
  float logsc = 0;
  float max;

  fwd->sc[0] = 0.0;

  if (L == 0) {
    fwd->sc[L+1] = logsc = log(hmm->pi[M]);
    if (opt_sc != NULL) *opt_sc = logsc;
    return eslOK;
  }

  max = 0.0;
  for (k = 0; k < M; k++) {
    fwd->dp[1][k] = hmm->e[k][dsq[1]] * hmm->pi[k];
    max = ESL_MAX(fwd->dp[1][k], max);
  }
  for (k = 0; k < M; k++) {
    fwd->dp[1][k] /= max;
  }
  fwd->sc[1] = log(max);

  for (i = 2; i <= L; i++)
    {
      max = 0.0;
      for (k = 0; k < M; k++)
	{
	  fwd->dp[i][k] = 0.0;
	  for (m = 0; m < M; m++)
	    fwd->dp[i][k] += fwd->dp[i-1][m] * hmm->t[m][k];

	  fwd->dp[i][k] *= hmm->e[k][dsq[i]];
	  
	  max = ESL_MAX(fwd->dp[i][k], max);
	}
      
      for (k = 0; k < M; k++)
	fwd->dp[i][k] /= max;
      fwd->sc[i] = log(max);
    }
	  
  
  fwd->sc[L+1] = 0.0;
  for (m = 0; m < M; m++) 
    fwd->sc[L+1] += fwd->dp[L][m] * hmm->t[m][M];
  fwd->sc[L+1] = log(fwd->sc[L+1]);

  logsc = 0.0;
  for (i = 1; i <= L+1; i++)
    logsc += fwd->sc[i];

  if (opt_sc != NULL) *opt_sc = logsc;
  return eslOK;
}


int
esl_hmm_Backward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *bck, float *opt_sc)
{
  int   i,k,m;
  int   M     = hmm->M;
  float logsc = 0.0;
  float max;
  
  bck->sc[L+1] = 0.0;

  if (L == 0) {
    bck->sc[0] = logsc = log(hmm->pi[M]);
    if (opt_sc != NULL) *opt_sc = logsc;
    return eslOK;
  }
  
  max = 0.0;
  for (k = 0; k < M; k++)
    {
      bck->dp[L][k] = hmm->t[k][M];
      max = ESL_MAX(bck->dp[L][k], max);
    }
  for (k = 0; k < M; k++)
    bck->dp[L][k] /= max;
  bck->sc[L] = log(max);

  for (i = L-1; i >= 1; i--)
    {
      max = 0.0;
      for (k = 0; k < M; k++)
	{
	  bck->dp[i][k] = 0.0;
	  for (m = 0; m < M; m++)
	    bck->dp[i][k] += bck->dp[i+1][m] * hmm->e[m][dsq[i+1]] * hmm->t[k][m];
	  
	  max = ESL_MAX(bck->dp[i][k], max);
	}

      for (k = 0; k < M; k++)
	bck->dp[i][k] /= max;
      bck->sc[i] = log(max);
    }

  bck->sc[0] = 0.0;
  for (m = 0; m < M; m++)
    bck->sc[0] += bck->dp[1][m] * hmm->e[m][dsq[1]] * hmm->pi[m];
  bck->sc[0] = log(bck->sc[0]);

  logsc = 0.0;
  for (i = 0; i <= L; i++) 
    logsc += bck->sc[i];

  if (opt_sc != NULL) *opt_sc = logsc;
  return eslOK;
}  
		   

int
esl_hmm_PosteriorDecoding(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *fwd, ESL_HMX *bck, ESL_HMX *pp)
{
  int i,k;

  for (i = 1; i <= L; i++)
    {
      for (k = 0; k < hmm->M; k++)
	pp->dp[i][k] = fwd->dp[i][k] * bck->dp[i][k];
      esl_vec_FNorm(pp->dp[i], hmm->M);
    }
  return eslOK;
}


int
esl_hmm_Expectation(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *fwd, ESL_HMX *bck, ESL_HMM *counts)
{
  int i, k;

  /* Expected emission counts */
  for (i = 1; i <= L; i++)
    {
      for (k = 0; k < hmm->M; k++)
	pp->dp[i][k] = fwd->dp[i][k] * bck->dp[i][k];
      esl_vec_FNorm(pp->dp[i], hmm->M);
      
      for (k = 0; k < hmm->M; k++)
      

    }



}


/*****************************************************************
 * x. Functions used in unit testing
 *****************************************************************/

static int
make_occasionally_dishonest_casino(ESL_HMM **ret_hmm, ESL_ALPHABET **ret_abc)
{
  ESL_ALPHABET *abc = esl_alphabet_Create(eslDICE);
  ESL_HMM      *hmm = esl_hmm_Create(abc, 2);
  int           x;

  /* State 0 = fair die */
  hmm->pi[0] = 1.0;
  hmm->pi[1] = 0.0;
  hmm->pi[2] = 0.0;		/* no L=0 seqs */

  hmm->t[0][0] = 0.96;
  hmm->t[0][1] = 0.03;
  hmm->t[0][2] = 0.01;		/* end from state 0; mean length 100 */

  for (x = 0; x < abc->K; x++)
    hmm->e[0][x] = 1.0 / (float) abc->K;

  /* State 1 = loaded die */
  hmm->t[1][0] = 0.05;
  hmm->t[1][1] = 0.95;
  hmm->t[1][2] = 0.0;		/* no end from state 1 */

  for (x = 0; x < abc->K-1; x++) hmm->e[1][x] = 0.5 / ((float) abc->K-1);
  hmm->e[1][abc->K-1] = 0.5;

  *ret_hmm = hmm;
  *ret_abc = abc;
  return eslOK;
}


  
/*****************************************************************
 * x. Test driver.
 *****************************************************************/
#ifdef eslHMM_TESTDRIVE
/* gcc -g -Wall -o hmm_utest -L. -I. -DeslHMM_TESTDRIVE esl_hmm.c -leasel -lm
 */
#include "esl_config.h"

#include <stdio.h>

#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_hmm.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  {"-r",  eslARG_NONE,     NULL, NULL, NULL, NULL, NULL, NULL, "use arbitrary random number seed",  0},
  {"-s",  eslARG_INT,      "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",     0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose output",               0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for hmm module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int             seed       = esl_opt_GetInteger(go, "-s");
  ESL_RANDOMNESS *r          = NULL;
  ESL_ALPHABET   *abc        = NULL;
  ESL_HMM        *hmm        = NULL;
  ESL_DSQ        *dsq        = NULL;
  int            *path       = NULL;
  ESL_HMX        *fwd        = NULL;
  ESL_HMX        *bck        = NULL;		
  ESL_HMX        *pp         = NULL;		
  float           fsc, bsc;
  int             L;
  int             i;
  float           fsum, bsum;

  if (esl_opt_GetBoolean(go, "-r")) r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(seed);

  make_occasionally_dishonest_casino(&hmm, &abc);

  esl_hmm_Emit(r, hmm, &dsq, &path, &L);

  fwd = esl_hmx_Create(L, hmm->M);
  bck = esl_hmx_Create(L, hmm->M);
  pp  = esl_hmx_Create(L, hmm->M);

  esl_hmm_Forward (dsq, L, hmm, fwd, &fsc);
  esl_hmm_Backward(dsq, L, hmm, bck, &bsc);
  esl_hmm_PosteriorDecoding(dsq, L, hmm, fwd, bck, pp);

  fsum = 0.0;
  bsum = bsc;

  fsum += fwd->sc[0];
  printf("%4d %c %s %8.3f %8.3f\n", 0, '-', "--", fwd->sc[0], bck->sc[0]);
  bsum -= bck->sc[0];

  for (i = 1; i <= L; i++)
    {
      fsum += fwd->sc[i];
      printf("%4d %c %s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
	     i, abc->sym[dsq[i]], path[i] == 0 ? "F " : " L", 
	     fwd->sc[i], bck->sc[i],
	     fsum, bsum, fsum+bsum,
	     pp->dp[i][0], pp->dp[i][1]);
      bsum -= fwd->sc[i];
    }

  printf("%4d %c %s %8.3f %8.3f\n", 0, '-', "--", fwd->sc[L+1], bck->sc[L+1]);
  printf("Forward score  = %f\n", fsc);
  printf("Backward score = %f\n", bsc);

  free(path);
  free(dsq);
  esl_hmx_Destroy(pp);
  esl_hmx_Destroy(bck);
  esl_hmx_Destroy(fwd);
  esl_alphabet_Destroy(abc);
  esl_hmm_Destroy(hmm);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslHMM_TESTDRIVE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
