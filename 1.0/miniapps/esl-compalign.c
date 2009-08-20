/* esl-compalign-rf - compare two sequence alignments
 *
 * EPN, Sun Aug  3 14:57:35 2008
 * From squid's compalign: Sean Eddy, Tue Nov  3 07:46:59 1992
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_msa.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

static char banner[] = "compare two multiple alignments";

static char usage[]  = "\
[-options] <trusted file> <test file>\n\
  Both files must be in Stockholm format with #=GC RF markup.\n\
  Sequences must occur in the same order in the two files.\n\
  Number of non-gap characters in #=GC RF markup must be identical.\n\
  Note: the scoring metric used is different from Squid\'s compalign.\n\
";

static int integerize_posterior_char(char c);
static int read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen);

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",      eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help; show brief info on version and usage",                     0 },
  { "-c",      eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "print per column statistics",                 0 },
  { "-p",      eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "print histogram of accuracy versus posterior probability",       0 },
  { "--c2dfile", eslARG_OUTFILE,NULL,NULL, NULL, NULL,"-c", NULL, "print per column stats to esl-ssudraw --dfile file <f>",         0 },
  { "--p2xm",    eslARG_OUTFILE,NULL,NULL, NULL, NULL,"-p", NULL, "print posterior stats to xmgrace file",         0 },
  { "--mask-p2xm", eslARG_OUTFILE,NULL,NULL, NULL, NULL,"--p2xm", NULL, "with --p2xm, only look at columns within mask in <f>",         0 },
  { "--amino",     eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--dna,--rna",               "<msafile> contains protein alignments",                         10 },
  { "--dna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--rna",             "<msafile> contains DNA alignments",                             10 },
  { "--rna",       eslARG_NONE,  FALSE, NULL, NULL,      NULL,NULL,"--amino,--dna",             "<msafile> contains RNA alignments",                             10 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;		/* application configuration       */
  int          kstatus, tstatus;/* return code from Easel routine  */
  int          fmt;		/* expected format of kfile, tfile */
  char        *kfile, *tfile;   /* known, test structure file      */
  ESL_MSAFILE *kfp, *tfp;       /* open kfile, tfile               */
  ESL_MSA     *ka,  *ta; 	/* known, trusted alignment        */
  int64_t      klen, tlen;	/* lengths of dealigned seqs       */
  int          i;		/* counter over sequences          */
  int          apos;		/* counter over alignment columns  */
  int          cpos;		/* counter over consensus (non-gap RF) columns  */
  int       is_cpos;            /* TRUE if current apos is a consensus pos, FALSE if not */
  int          uapos;		/* counter over unaligned residue positions */

  int        **kp;              /* [0..i..nseq-1][1..r..sq->n] = x known non-gap RF position of residue r in sequence i */
  int        **tp;              /* [0..i..nseq-1][1..r..sq->n] = x predicted non-gap RF position of residue r in sequence i */
  /* for both kp and pp, if x <= 0, residue r for seq i is not aligned to a non-gap RF position, but rather as an 'insert'
   * after non-gap RF position (x * -1) 
   */
  int        *km_pos;          /* [0..rflen] = x, in known aln,     number of residues aligned to non-gap RF column x; special case: mct[0] = 0 */
  int        *ki_pos;          /* [0..rflen] = x, in known aln,     number of residues inserted after non-gap RF column x */
  int        *tm_pos;          /* [0..rflen] = x, in predicted aln, number of residues aligned to non-gap RF column x; special case: mct[0] = 0 */
  int        *ti_pos;          /* [0..rflen] = x, in predicted aln, number of residues inserted after non-gap RF column x */
  int    *cor_tm_pos;          /* [0..rflen] = x, in predicted aln, number of correctly predicted residues aligned to non-gap RF column x; special case: mct[0] = 0 */
  int    *cor_ti_pos;          /* [0..rflen] = x, in predicted aln, number of correctly predicted residues inserted after non-gap RF column x */

  int        *km_seq;          /* [0..i..nseq-1] = x, in known aln,     number of residues aligned to non-gap RF columns in seq i; */
  int        *ki_seq;          /* [0..i..nseq-1] = x, in known aln,     number of residues inserted in seq i */
  int        *tm_seq;          /* [0..i..nseq-1] = x, in predicted aln, number of residues aligned to non-gap RF columns in seq i; */
  int        *ti_seq;          /* [0..i..nseq-1] = x, in predicted aln, number of residues inserted in seq i */
  int    *cor_tm_seq;          /* [0..i..nseq-1] = x, in predicted aln, number of correctly predicted residues aligned to non-gap RF columns in seq i */
  int    *cor_ti_seq;          /* [0..i..nseq-1] = x, in predicted aln, number of correctly predicted residues inserted in seq i */

  int     *seqlen;             /* [0..i..nseq-1] = x, unaligned seq i has length x */
  ESL_ALPHABET *abc;           /* alphabet for all alignments */
  int      rflen, t_rflen;     /* non-gap RF length (consensus lengths) */
  int   status;

  /* variables needed for -p and related options */
  int do_post = FALSE; /* TRUE if -p enabled */
  int do_post_for_this_cpos = FALSE; /* set for each consensus position, always TRUE unless --mask-p2xm */
  int p;               /* counter over integerized posteriors */
  int ridx1 = -1;      /* #=GR index for posterior value digit 2 */
  int ridx2 = -1;      /* #=GR index for posterior value digit 2 */
  int ndigits = 0;     /* number of posterior digits in alignment, 1 or 2 */
  int *ptm = NULL;     /* [0..p..100] number of total   matches with int posteriors of p */
  int *pti = NULL;     /* [0..p..100] number of total   inserts with int posteriors of p */
  int *cor_ptm = NULL; /* [0..p..100] number of correct matches with int posteriors of p */
  int *cor_pti = NULL; /* [0..p..100] number of correct inserts with int posteriors of p */
  int npostvals = 101; /* number of posterior int values 0..101 */
  int r;               /* counter over #=GR annotation */
  int pint;            /* integerized posterior */
  int cm_cor_ptm, cm_cor_pti, cm_ptm, cm_pti, cm_incor_ptm, cm_incor_pti; /* cumulative counts of posteriors */
  int tot_cor_ptm, tot_cor_pti, tot_ptm, tot_pti, tot_incor_ptm, tot_incor_pti; /* total counts of posteriors */
  char          errbuf[eslERRBUFSIZE];

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
      exit(EXIT_SUCCESS);
    }

  if (esl_opt_ArgNumber(go) != 2) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  kfile = esl_opt_GetArg(go, 1);
  tfile = esl_opt_GetArg(go, 2);
  
  fmt = eslMSAFILE_STOCKHOLM;

  /***********************************************
   * Open the two Stockholm files.
   ***********************************************/

  if (esl_msafile_Open(kfile, fmt, NULL, &kfp) != eslOK)
    esl_fatal("Failed to open trusted structure file %s for reading", kfile);
  if (esl_msafile_Open(tfile, fmt, NULL, &tfp) != eslOK)
    esl_fatal("Failed to open test structure file %s for reading", tfile);
  
  if      (esl_opt_GetBoolean(go, "--amino"))   abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     abc = esl_alphabet_Create(eslRNA);
  else {
    int type;
    status = esl_msafile_GuessAlphabet(kfp, &type);
    if (status == eslEAMBIGUOUS)    esl_fatal("Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", kfile);
    else if (status == eslEFORMAT)  esl_fatal("Alignment file parse failed: %s\n", kfp->errbuf);
    else if (status == eslENODATA)  esl_fatal("Alignment file %s is empty\n", kfile);
    else if (status != eslOK)       esl_fatal("Failed to read alignment file %s\n", kfile);
    abc = esl_alphabet_Create(type);
  }
  /* set both as same alphabet */
  esl_msafile_SetDigital(kfp, abc);
  esl_msafile_SetDigital(tfp, abc);

  do_post = esl_opt_GetBoolean(go, "-p");

  /* read the mask file if --mask-p2xm is enabled */
  char *mask = NULL;
  int masklen;
  if(! esl_opt_IsDefault(go, "--mask-p2xm")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-p2xm"), errbuf, &mask, &masklen)) != eslOK) esl_fatal(errbuf);
  }

  /***********************************************
   * Do alignment comparisons, one seq at a time;
   * this means looping over all seqs in all alignments.
   ***********************************************/

  while (1)
    {
      kstatus = esl_msa_Read(kfp, &ka);
      tstatus = esl_msa_Read(tfp, &ta);
      if (kstatus != eslOK || tstatus != eslOK) break; /* normal or errors. */

      /* Sanity check on alignment
       */
      if (ka->nseq != ta->nseq)
	esl_fatal("trusted, test alignments don't have same seq #\n");
      if (ka->rf == NULL)
	esl_fatal("trusted alignment has no reference annotation\n");
      if (ta->rf == NULL)
	esl_fatal("test alignment has no reference annotation\n");

      /* make sure the sequences are all identical */
      ESL_DSQ *ks;
      ESL_DSQ *ts;
      ESL_ALLOC(seqlen, sizeof(int) * ka->nseq);
      for(i = 0; i < ka->nseq; i++) { 
	if(strcmp(ka->sqname[i], ta->sqname[i]) != 0) esl_fatal("sequence i of trusted alignment %s has different name than seq i of predicted alignment %s\n", ka->sqname[i], ta->sqname[i]); 
	ESL_ALLOC(ks, sizeof(ESL_DSQ) * (ka->alen+2));
	memcpy(ks, ka->ax[i], (ka->alen+2) * sizeof(ESL_DSQ));
	esl_abc_XDealign(ka->abc, ks, ka->ax[i], &klen);

	ESL_ALLOC(ts, sizeof(ESL_DSQ) * (ta->alen+2));
	memcpy(ts, ta->ax[i], (ta->alen+2) * sizeof(ESL_DSQ));
	esl_abc_XDealign(ta->abc, ts, ta->ax[i], &tlen);

	if (tlen != klen)
	  esl_fatal("dealigned sequence mismatch, seq %d, when dealigned, is %d residues in the known alignment, but %d residues in the trusted alignment.", i, klen, tlen);

	if (memcmp(ks, ts, sizeof(ESL_DSQ) * klen) != 0) 
	  esl_fatal("dealigned sequence mismatch, seq %d %s, when dealigned, are not identical.", i, ka->sqname[i]);

	seqlen[i] = tlen;
	free(ks);
	free(ts);
      }

      /* determine non-gap RF length */
      rflen = 0;
      for(apos = 1; apos <= ka->alen; apos++) { 
	if(! (esl_abc_CIsGap(ka->abc, ka->rf[(apos-1)]))) rflen++;
      }
      t_rflen = 0;
      for(apos = 1; apos <= ta->alen; apos++) { 
	if(! (esl_abc_CIsGap(ta->abc, ta->rf[(apos-1)]))) t_rflen++;
      }
      if(t_rflen != rflen) esl_fatal("Trusted alignment non-gap RF length (%d) != predicted alignment non-gap RF length (%d).\n", rflen, t_rflen);

      /* if -p, make sure the test alignment has posterior probabilities, and allocate our counters for correct/incorrect per post value */
      if(do_post) { 
	if(! esl_opt_IsDefault(go, "--mask-p2xm")) {
	  if(masklen != rflen) { 
	    esl_fatal("Length of mask in %s (%d) not equal to non-gap RF len of alignments (%d)\n", esl_opt_GetString(go, "--mask-p2xm"), masklen, rflen);
	  }
	}
	for (r = 0; r < ta->ngr; r++) { 
	  if (strcmp(ta->gr_tag[r], "POST")   == 0) { ridx1 = r; ndigits = 1; }
	  if (strcmp(ta->gr_tag[r], "Post")   == 0) { ridx1 = r; ndigits = 1; }
	  if (strcmp(ta->gr_tag[r], "post")   == 0) { ridx1 = r; ndigits = 1; }
	  if (strcmp(ta->gr_tag[r], "POSTX.") == 0) { ridx1 = r; ndigits = 1; }
	  if (strcmp(ta->gr_tag[r], "POST.X") == 0) { ridx2 = r; ndigits = 2; }
	}
	if(ndigits == 0)                                 esl_fatal("-p requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", tfile);
	if(ndigits == 1 && ridx1 == -1)                  esl_fatal("-p requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", tfile);
	if(ndigits == 2 && (ridx1 == -1 || ridx2 == -1)) esl_fatal("-p requires \"#=GR POST\", \"#=GR Post\", \"#=GR post\", \"#=GR POSTX.\", or \"#=GR POSTX.\" and \"#=GR POST.X\" annotation in %s.\n", tfile);
	/* always allocate 0..100, if we only have 1 post value, only 0,10,20,30,40,50,60,70,80,90,100 will be filled with counts */
	ESL_ALLOC(ptm,     sizeof(int) * npostvals);
	ESL_ALLOC(pti,     sizeof(int) * npostvals);
	ESL_ALLOC(cor_ptm, sizeof(int) * npostvals);
	ESL_ALLOC(cor_pti, sizeof(int) * npostvals);
	esl_vec_ISet(ptm, npostvals, 0);
	esl_vec_ISet(pti, npostvals, 0);
	esl_vec_ISet(cor_ptm, npostvals, 0);
	esl_vec_ISet(cor_pti, npostvals, 0);
      }

      /* allocate and initialize our counters */
      ESL_ALLOC(kp, sizeof(int *) * ka->nseq);
      ESL_ALLOC(tp, sizeof(int *) * ta->nseq);
      for(i = 0; i < ka->nseq; i++) { 
	ESL_ALLOC(kp[i], sizeof(int) * (seqlen[i]+1));
	ESL_ALLOC(tp[i], sizeof(int) * (seqlen[i]+1));
	esl_vec_ISet(kp[i], seqlen[i]+1, -987654321);
	esl_vec_ISet(tp[i], seqlen[i]+1, -987654321);
      }

      ESL_ALLOC(km_pos, sizeof(int) * (rflen+1));
      ESL_ALLOC(ki_pos, sizeof(int) * (rflen+1));
      ESL_ALLOC(tm_pos, sizeof(int) * (rflen+1));
      ESL_ALLOC(ti_pos, sizeof(int) * (rflen+1));
      ESL_ALLOC(cor_tm_pos, sizeof(int) * (rflen+1));
      ESL_ALLOC(cor_ti_pos, sizeof(int) * (rflen+1));
      esl_vec_ISet(km_pos, rflen+1, 0);
      esl_vec_ISet(ki_pos, rflen+1, 0);
      esl_vec_ISet(tm_pos, rflen+1, 0);
      esl_vec_ISet(ti_pos, rflen+1, 0);
      esl_vec_ISet(cor_tm_pos, rflen+1, 0);
      esl_vec_ISet(cor_ti_pos, rflen+1, 0);

      ESL_ALLOC(km_seq, sizeof(int) * ka->nseq);
      ESL_ALLOC(ki_seq, sizeof(int) * ka->nseq);
      ESL_ALLOC(tm_seq, sizeof(int) * ka->nseq);
      ESL_ALLOC(ti_seq, sizeof(int) * ka->nseq);
      ESL_ALLOC(cor_tm_seq, sizeof(int) * ka->nseq);
      ESL_ALLOC(cor_ti_seq, sizeof(int) * ka->nseq);
      esl_vec_ISet(km_seq, ka->nseq, 0);
      esl_vec_ISet(ki_seq, ka->nseq, 0);
      esl_vec_ISet(tm_seq, ka->nseq, 0);
      esl_vec_ISet(ti_seq, ka->nseq, 0);
      esl_vec_ISet(cor_tm_seq, ka->nseq, 0);
      esl_vec_ISet(cor_ti_seq, ka->nseq, 0);

      /* determine non-gap RF location of each residue in known alignment */
      for(i = 0; i < ka->nseq; i++) { 
	uapos = cpos = 0;
	for(apos = 1; apos <= ka->alen; apos++) { 
	  is_cpos = FALSE;
	  if(! (esl_abc_CIsGap(ka->abc, ka->rf[(apos-1)]))) { 
	    cpos++; is_cpos = TRUE;
	  }
	  if(! esl_abc_XIsGap(ka->abc, ka->ax[i][apos])) { 
	    uapos++;
	    kp[i][uapos] = (is_cpos) ? cpos : (-1 * cpos);
	    if(is_cpos) { km_pos[cpos]++; km_seq[i]++; }
	    else        { ki_pos[cpos]++; ki_seq[i]++; }
	  }
	}
      }

      /* determine non-gap RF location of each residue in predicted alignment */
      for(i = 0; i < ta->nseq; i++) { 
	uapos = cpos = 0;
	for(apos = 1; apos <= ta->alen; apos++) { 
	  is_cpos = FALSE;
	  if(! (esl_abc_CIsGap(ta->abc, ta->rf[(apos-1)]))) { 
	    cpos++; is_cpos = TRUE;
	    if(do_post) { 
	      do_post_for_this_cpos = (mask != NULL && mask[cpos-1] == '0') ? FALSE : TRUE;
	    }
	  }
	  if(! esl_abc_XIsGap(ta->abc, ta->ax[i][apos])) { 
	    uapos++;
	    tp[i][uapos] = (is_cpos) ? cpos : (-1 * cpos);
	    if(do_post) { 
	      pint = 10 * integerize_posterior_char(ta->gr[ridx1][i][(apos-1)]);
	      if(ndigits == 2 && pint != 100) pint += integerize_posterior_char(ta->gr[ridx2][i][(apos-1)]);
	    }
	    if(is_cpos) { 
	      tm_pos[cpos]++; tm_seq[i]++; 
	      if(do_post_for_this_cpos) ptm[pint]++;
	    }
	    else { 
	      ti_pos[cpos]++; ti_seq[i]++; 
	      if(do_post) pti[pint]++;
	    }
	    if(kp[i][uapos] == tp[i][uapos]) { /* correctly predicted this residue */
	      if(is_cpos) { 
		cor_tm_seq[i]++; cor_tm_pos[cpos]++; 
		if(do_post_for_this_cpos) cor_ptm[pint]++;
	      } 
	      else {
		cor_ti_seq[i]++; cor_ti_pos[cpos]++; 
		if(do_post) cor_pti[pint]++;
	      } 
	    }
	  }
	}
      }
      if((! (esl_opt_GetBoolean(go, "-c"))) && (! esl_opt_GetBoolean(go, "-p"))) { 
	/* print per sequence statistics */
	char *namedashes;
	int ni;
	int namewidth = 8; /* length of 'seq name' */
	/* determine the longest name in msa */
	for(ni = 0; ni < ka->nseq; ni++) namewidth = ESL_MAX(namewidth, strlen(ka->sqname[ni]));
	ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
	namedashes[namewidth] = '\0';
	for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';
	
	printf("# %-*s  %5s  %20s  %20s  %20s\n", namewidth, "seq name", "len", "match columns", "insert columns", "all columns");
	printf("# %-*s  %5s  %20s  %20s  %20s\n", namewidth, namedashes, "-----", "--------------------", "--------------------", "--------------------");
	for(i = 0; i < ta->nseq; i++) { 
	  printf("  %-*s  %5d  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)\n", namewidth, ka->sqname[i], seqlen[i],
		 cor_tm_seq[i], km_seq[i], (km_seq[i] == 0) ? 0. : ((float) cor_tm_seq[i] / (float) km_seq[i]), 
		 cor_ti_seq[i], ki_seq[i], (ki_seq[i] == 0) ? 0. : ((float) cor_ti_seq[i] / (float) ki_seq[i]), 
		 (cor_tm_seq[i] + cor_ti_seq[i]), (km_seq[i] + ki_seq[i]), ((float) (cor_tm_seq[i] + cor_ti_seq[i]) / ((float) km_seq[i] + ki_seq[i]))); 
	}
	int cor_tm, cor_ti, km, ki;
	cor_tm = esl_vec_ISum(cor_tm_seq, ka->nseq);
	cor_ti = esl_vec_ISum(cor_ti_seq, ka->nseq);
	km = esl_vec_ISum(km_seq, ka->nseq);
	ki = esl_vec_ISum(ki_seq, ka->nseq);
	
	printf("# %-*s  %5s  %20s  %20s  %20s\n", namewidth, namedashes, "-----", "--------------------", "--------------------", "--------------------");
	printf("# %-*s  %5s  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)\n",
	       namewidth, "*all*", "-", 
	       cor_tm, km, ((float) cor_tm / (float) km), 
	       cor_ti, ki, ((float) cor_ti / (float) ki), 
	       (cor_tm+cor_ti), (km+ki), (((float) (cor_tm + cor_ti))/ ((float) (km + ki)))); 
	free(namedashes);
	for(i = 0; i < ka->nseq; i++) { 
	  free(kp[i]); 
	  free(tp[i]); 
	}
      }
      else if(esl_opt_GetBoolean(go, "-c")) { /* print per column statistics */
	printf("# %5s  %20s  %20s  %20s\n", "rfpos", "match", "insert", "both");
	printf("# %5s  %20s  %20s  %20s\n", "-----", "--------------------", "--------------------", "--------------------");
	for(cpos = 0; cpos <= rflen; cpos++) { 
	  printf("  %5d  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)  %4d / %4d  (%.3f)\n", cpos, 
		 
		 cor_tm_pos[cpos], km_pos[cpos], (km_pos[cpos] == 0) ? 0. : ((float) cor_tm_pos[cpos] / (float) km_pos[cpos]), 
		 cor_ti_pos[cpos], ki_pos[cpos], (ki_pos[cpos] == 0) ? 0. : ((float) cor_ti_pos[cpos] / (float) ki_pos[cpos]), 
		 (cor_tm_pos[cpos] + cor_ti_pos[cpos]), (km_pos[cpos] + ki_pos[cpos]), ((float) (cor_tm_pos[cpos] + cor_ti_pos[cpos]) / ((float) km_pos[cpos] + ki_pos[cpos]))); 
	}
      }
      else if(do_post) { /* do posterior output */
	FILE *pfp = NULL;
	  /* pfp will be an xmgrace file
	   * first series:
	   * x axis = posterior probability 
	   * y axis = fraction of match residues at >= pp x that are correct 
	   */
	if (esl_opt_GetString(go, "--p2xm") != NULL) {
	  if ((pfp = fopen(esl_opt_GetString(go, "--p2xm"), "w")) == NULL) 
	    esl_fatal("Failed to open --p2xm output file %s\n", esl_opt_GetString(go, "--p2xm"));
	}

	if(mask == NULL) printf("# %4s  %44s  %44s\n", "prob", "match columns             ", "insert columns             ");
	else             printf("# %4s  %44s  %44s\n", "prob", "match columns within mask ", "insert columns             ");
	printf("# %4s  %44s  %44s\n", "----", "--------------------------------------------", "---------------------------------------------");
	cm_ptm = cm_pti = cm_cor_ptm = cm_cor_pti = cm_incor_ptm = cm_incor_pti = 0;
	tot_ptm = esl_vec_ISum(ptm, npostvals);
	tot_pti = esl_vec_ISum(pti, npostvals);
	tot_cor_ptm = esl_vec_ISum(cor_ptm, npostvals);
	tot_cor_pti = esl_vec_ISum(cor_pti, npostvals);
	tot_incor_ptm = tot_ptm - tot_cor_ptm;
	tot_incor_pti = tot_pti - tot_cor_pti;
	for(p = (npostvals-1); p >= 0; p--) { 
	  cm_cor_ptm += cor_ptm[p];
	  cm_cor_pti += cor_pti[p];
	  cm_ptm     += ptm[p];
	  cm_pti     += pti[p];
	  cm_incor_ptm += ptm[p] - cor_ptm[p];
	  cm_incor_pti += pti[p] - cor_pti[p];
	  printf("  %4d %8d / %8d (%.5f) (%.5f) (%.5f)  %8d / %8d (%.5f) (%.5f) (%.5f)\n", 
		 p, cor_ptm[p], ptm[p], 
		 (ptm[p] == 0) ? 0. : (float) cor_ptm[p] / (float) ptm[p], 
		 (cm_ptm == 0) ? 0. : (float) cm_cor_ptm / (float) cm_ptm, 
		 (tot_incor_ptm == 0) ? 0. : (float) cm_incor_ptm / (float) tot_incor_ptm, 
		 cor_pti[p], pti[p], 
		 (pti[p] == 0) ? 0. : (float) cor_pti[p] / (float) pti[p],
		 (cm_pti == 0) ? 0. : (float) cm_cor_pti / (float) cm_pti,
		 (tot_incor_pti == 0) ? 0. : (float) cm_incor_pti / (float) tot_incor_pti);
	  if(pfp != NULL) fprintf(pfp, "%f %f\n", (float) p / 100., (cm_ptm == 0) ? 0. : (float) cm_cor_ptm / (float) cm_ptm);
	}
	if(pfp != NULL) { 
#if 0
	  /* x axis = posterior probability 
	   * y axis = fraction of match residues at (exactly) pp x that are correct 
	   */
	  fprintf(pfp, "&\n");
	  cm_ptm = cm_pti = cm_cor_ptm = cm_cor_pti = cm_incor_ptm = cm_incor_pti = 0;
	  for(p = (npostvals-1); p >= 0; p--) {
	    cm_cor_ptm += cor_ptm[p];
	    cm_cor_pti += cor_pti[p];
	    cm_ptm     += ptm[p];
	    cm_pti     += pti[p];
	    cm_incor_ptm += ptm[p] - cor_ptm[p];
	    cm_incor_pti += pti[p] - cor_pti[p];
	    fprintf(pfp, "%f %f\n", (float) p / 100., (ptm[p] == 0) ? 0. : (float) cor_ptm[p] / (float) ptm[p]);
	  }
	  fprintf(pfp, "&\n");
#endif
#if 0
	  /* x axis = posterior probability 
	   * y axis = fraction of all the incorrect residues with pp < x
	   */
	  cm_ptm = cm_pti = cm_cor_ptm = cm_cor_pti = cm_incor_ptm = cm_incor_pti = 0;
	  for(p = (npostvals-1); p >= 0; p--) {
	    cm_cor_ptm += cor_ptm[p];
	    cm_cor_pti += cor_pti[p];
	    cm_ptm     += ptm[p];
	    cm_pti     += pti[p];
	    cm_incor_ptm += ptm[p] - cor_ptm[p];
	    cm_incor_pti += pti[p] - cor_pti[p];
	    fprintf(pfp, "%f %f\n", (float) p / 100., (tot_incor_ptm == 0) ? 0. : (float) cm_incor_ptm / (float) tot_incor_ptm);
	  }
	  fprintf(pfp, "&\n");
#endif
#if 0
	  /* print the identity line for comparison */
	  for(p = (npostvals-1); p >= 0; p--)
 	    fprintf(pfp, "%f %f\n", (float) p / 100., (float) p / 100.);
	  fprintf(pfp, "&\n");
#endif
	  fclose(pfp);
	}
      }
      /* handle --c2dfile */
      FILE *dfp;
      if (esl_opt_GetString(go, "--c2dfile") != NULL) {
	if ((dfp = fopen(esl_opt_GetString(go, "--c2dfile"), "w")) == NULL) 
	  esl_fatal("Failed to open --c2dfile output file %s\n", esl_opt_GetString(go, "--c2dfile"));
	/* match stats, 4 fields, CMYK color values */
	for(cpos = 1; cpos <= rflen; cpos++) { 
	  if(km_pos[cpos] == 0) { /* special case, no known alignment residues, a blank position */
	    fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 0., 0., 0., 0.);
	  }
	  else { 
	    fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 
		    0., /* cyan */
		    1. - ((float) cor_tm_pos[cpos] / (float) km_pos[cpos]), /* magenta, fraction correct */
		    1. - ((float) km_pos[cpos] / ta->nseq), /* yellow, fraction of seqs with residue in column */
		    0.);
	  }		 
	}	
	fprintf(dfp, "//\n");
	/* insert stats, 4 fields, CMYK color values */
	cpos = 0; /* special case, combine insert posn 0 and 1 together */
	if(ki_pos[cpos] == 0) { /* special case, no known alignment residues, a blank position */
	  fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 0., 0., 0., 0.);
	}
	else { 
	  fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 
		  0., /* cyan */
		  1. - ((float) (cor_ti_pos[0] + cor_ti_pos[1]) / ((float) (ki_pos[0] + ki_pos[1]))), /* magenta, fraction correct */
		  0.,
		  0.);
	}
	/* insert stats posn 2..rflen */
	for(cpos = 2; cpos <= rflen; cpos++) { 
	  if(ki_pos[cpos] == 0) { /* special case, no known alignment residues, a blank position */
	    fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 0., 0., 0., 0.);
	  }
	  else { 
	    fprintf(dfp, "%.3f %.3f %.3f %.3f\n", 
		    0., /* cyan */
		    1. - ((float) cor_ti_pos[cpos] / (float) ki_pos[cpos]), /* magenta, fraction correct */
		    0.,
		    0.);
	  }
	} 
	fprintf(dfp, "//\n");
	fclose(dfp);
      }
      
      if(ptm != NULL) free(ptm);
      if(pti != NULL) free(pti);
      if(cor_ptm != NULL) free(cor_ptm);
      if(cor_ptm != NULL) free(cor_pti);
      free(kp);
      free(tp);
      free(km_seq);
      free(ki_seq);
      free(tm_seq);
      free(ti_seq);
      free(cor_tm_seq);
      free(cor_ti_seq);
      free(km_pos);
      free(ki_pos);
      free(tm_pos);
      free(ti_pos);
      free(cor_tm_pos);
      free(cor_ti_pos);
      free(seqlen);
      esl_msa_Destroy(ka);
      esl_msa_Destroy(ta);
    }
  
  /* At this point, we should have EOF status on both
   * alignment files; if we don't, there's an error we have to handle.
   */
  if (kstatus != eslEOF || tstatus != eslEOF)
    {
      if (kstatus == eslEFORMAT)
	esl_fatal("Parse error, line %d of trusted file %s:\n%s\n",
		  kfp->linenumber, kfp->fname, kfp->errbuf);
      if (tstatus == eslEFORMAT)
	esl_fatal("Parse error, line %d of test file %s:\n%s\n",
		  tfp->linenumber, tfp->fname, tfp->errbuf);
      if (kstatus == eslOK) 
	esl_fatal("Trusted file has more data than test file\n");
      if (tstatus == eslOK)
	esl_fatal("Test file has more data than trusted file\n");
      if (kstatus != eslEOF)
	esl_fatal("read error %d for trusted file\n", kstatus);
      if (tstatus != eslEOF)
	esl_fatal("read error %d for test file\n", tstatus);
    }

  if(mask != NULL) free(mask);
  esl_getopts_Destroy(go);
  esl_msafile_Close(tfp);
  esl_msafile_Close(kfp);
  return 0;

 ERROR:
  return status;
}

/* integerize_posterior_char
 *                   
 * Return a integer 0..10 that is the discretized integer form of
 * a posterior probability character 'c'
 * If the posterior annotation is a gap or otherwise bogus we die.
 */      
int
integerize_posterior_char(char c)
{
  if(c == '*') return 10;
  int i = (int) c;
  if(i >= 48 && i <= 57) return i-48; /* '0' is 48, '9' is 57 */
  else esl_fatal("Don't know what to do with posterior value: %c\n", c);
  return 0; /* NEVERREACHED */
}



/* read_mask_file
 *
 * Given an open file pointer, read the first token of the
 * file and return it as *ret_mask.
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;

  if (esl_fileparser_Open(filename, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));
  for(n = 0; n < toklen; n++) mask[n] = tok[n];
  mask[n] = '\0';

  *ret_mask = mask;
  *ret_masklen= toklen;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}
