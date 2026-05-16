/* cmd_printseq : format and print a DNA/RNA sequence
 */
#include <esl_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"


#define ALPH_OPTS "--dna,--rna"   /* toggle group: alphabet type options */

static ESL_OPTIONS cmd_options[] = {
  /* name           type         default   env  range toggles    reqs  incomp   help                                          docgroup */
  { "-h",         eslARG_NONE,    FALSE,   NULL, NULL,  NULL,      NULL,  NULL,  "show brief help on version and usage",              0 },
  { "-3",         eslARG_NONE,    FALSE,   NULL, NULL,  NULL,      NULL,  "-a",  "groups of 3, frame 0 translation in 3-letter code", 0 },
  { "-a",         eslARG_NONE,    FALSE,   NULL, NULL,  NULL,      NULL,  "-3",  "show translation in all three reading frames",      0 },
  { "-f",         eslARG_INT,       "1",   NULL,"n>=1", NULL,      NULL,  NULL,  "number first residue as <n>",                       0 },
  { "-g",         eslARG_INT,      "10",   NULL,"n>=0", NULL,      NULL,  NULL,  "print in groups of <n> (0 = no grouping)",          0 },
  { "-n",         eslARG_INT,      "60",   NULL,"n>=10",NULL,      NULL,  NULL,  "number of symbols per line",                        0 },
  { "-s",         eslARG_NONE,    FALSE,   NULL, NULL,  NULL,      NULL,  NULL,  "print single-stranded (no complement)",             0 },
  { "-t",         eslARG_NONE,    FALSE,   NULL, NULL,  NULL,      NULL,  NULL,  "print coordinates at top instead of left",          0 },
  { "--informat", eslARG_STRING,   NULL,   NULL, NULL,  NULL,      NULL,  NULL,  "specify that input file is in format <s>",          0 },
  { "--dna",      eslARG_NONE,    FALSE,   NULL, NULL,  ALPH_OPTS, NULL,  NULL,  "specify that <seqfile> contains DNA sequence",      0 },
  { "--rna",      eslARG_NONE,    FALSE,   NULL, NULL,  ALPH_OPTS, NULL,  NULL,  "specify that <seqfile> contains RNA sequence",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static const char *aa3(char c);
static int         coord_width(int64_t val);
static int64_t     frame_start(int64_t pos, int frame);
static void        print_top_coords(int64_t pos, int resoffset, int nt_perline, int group_by, int64_t L);
static int         print_seq_line (const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int64_t pos, int64_t L, int nt_perline, int group_by);
static int         print_comp_line(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int64_t pos, int64_t L, int nt_perline, int group_by);

/* esl_cmd_printseq():  implements `easel printseq`
 *
 * <topcmd> is original argv[0]: `easel` or `/path/to/easel`.
 *
 * <sub> is the <ESL_SUBCMD> data corresponding to this subcommand, passed
 * from the `easel` program:
 *    sub->func        = esl_cmd_printseq
 *    sub->subcmd      = "printseq"
 *    sub->nargs       = 1
 *    sub->usage       = usage string defined in miniapps/easel.c 
 *    sub->description = help string defined in miniapps/easel.c
 *
 * <argc> is the number of subcommand arguments, including "printseq" but
 * not including "easel" or any top command options.
 *
 * <argv> is the list of subcommand arguments, starting with argv[0] =
 * "printseq".
 */
int
esl_cmd_printseq(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS   *go           = esl_subcmd_CreateDefaultApp(topcmd, sub, cmd_options, argc, argv, /*custom opthelp?:*/NULL);
  char          *seqfile      = esl_opt_GetArg(go, 1);
  ESL_ALPHABET  *nt_abc       = NULL;
  ESL_ALPHABET  *aa_abc       = NULL;  // only allocated if -a or -3
  ESL_GENCODE   *gcode        = NULL;  // only allocated if -a or -3
  ESL_SQFILE    *sqfp         = NULL;
  ESL_SQ        *sq           = NULL;
  ESL_SQ        *sq2          = NULL;
  int            do_translate = esl_opt_GetBoolean(go, "-a");
  int            do_triplet   = esl_opt_GetBoolean(go, "-3");
  int            resoffset    = esl_opt_GetInteger(go, "-f") - 1;  // offset to add to all coords
  int            group_by     = esl_opt_GetInteger(go, "-g");      // 0: no groups (-a); default display: groups of 10; 3-letter code: 3
  int            nt_perline   = esl_opt_GetInteger(go, "-n");
  int            do_single    = esl_opt_GetBoolean(go, "-s");
  int            do_topcoord  = esl_opt_GetBoolean(go, "-t");
  int            infmt        = eslSQFILE_UNKNOWN;
  int            alphatype    = eslUNKNOWN;
  int64_t        L;
  int            cw;                // coordinate field width
  int64_t        pos;
  int            nprinted;
  int64_t        fpos;
  int            frame, i, aa;
  int            status;

  /* Enforce option constraints */
  if (do_translate)
    {
      if ( ! esl_opt_IsDefault(go, "-g") && group_by != 0) esl_fatal("-a requires -g 0 (no grouping)");
      group_by = 0;
    }
  if (do_triplet)
    {
      if ( ! esl_opt_IsDefault(go, "-g") && group_by != 3) esl_fatal("-3 requires -g 3");
      group_by = 3;
    }
  if (do_topcoord && group_by == 0)
    esl_fatal("-t requires a nonzero group size (-g)");

  if (group_by > 0 && nt_perline % group_by != 0)
    esl_fatal("number of bases per line (-n) must be an exact multiple of group size (-g)");


  /* --informat: assert the input file format */
  if (esl_opt_IsOn(go, "--informat") &&
      (infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslSQFILE_UNKNOWN)
    esl_fatal("%s is not a valid input sequence file format for --informat", esl_opt_GetString(go, "--informat"));

  /* Open in text mode first, so we can guess the alphabet if needed */
  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Couldn't open %s for reading",    seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of %s", seqfile);
  else if (status != eslOK)        esl_fatal("Open of %s failed, code %d",      seqfile, status);

  /* Alphabet: --dna|--rna assert; otherwise autodetect from the first sequence */
  if      (esl_opt_GetBoolean(go, "--dna")) alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--rna")) alphatype = eslRNA;
  else
    {
      status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
      if      (status == eslENOALPHABET) esl_fatal("Couldn't guess alphabet from first sequence in %s", seqfile);
      else if (status == eslEFORMAT)     esl_fatal("Parse failed (sequence file %s):\n  %s", seqfile, esl_sqfile_GetErrorBuf(sqfp));
      else if (status == eslENODATA)     esl_fatal("Sequence file %s contains no data?", seqfile);
      else if (status != eslOK)          esl_fatal("Failed to guess alphabet (error code %d)", status);
      if (alphatype == eslAMINO)
        esl_fatal("%s appears to contain protein sequence; easel printseq requires DNA or RNA", seqfile);
    }
  nt_abc = esl_alphabet_Create(alphatype);
  sq     = esl_sq_CreateDigital(nt_abc);
  esl_sqfile_SetDigital(sqfp, nt_abc);

  if (do_translate || do_triplet)
    {
      aa_abc = esl_alphabet_Create(eslAMINO);
      gcode  = esl_gencode_Create(nt_abc, aa_abc);
    }

  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Parse failed:\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status == eslEOF)     esl_fatal("No sequence found in %s", seqfile);
  else if (status != eslOK)      esl_fatal("Unexpected read error %d", status);

  L  = sq->n;
  cw = coord_width(L + resoffset) + 1;
  printf("# %s (%" PRId64 " nt)\n\n", sq->name, L);

  /* Main output loop */
  pos = 0;
  while (pos < L)
    {
      /* Top strand */
      if (do_topcoord)   print_top_coords(pos, resoffset, nt_perline, group_by, L);
      else               printf("%*" PRId64 ": ", cw, pos + 1 + resoffset);
      nprinted = print_seq_line(nt_abc, sq->dsq, pos, L, nt_perline, group_by);
      if (! do_topcoord) printf(" %" PRId64, pos + nprinted + resoffset);
      putchar('\n');

      /* Complement strand */
      if (! do_single)
        {
          if (! do_topcoord) printf("%*s: ", cw, "");
          print_comp_line(nt_abc, sq->dsq, pos, L, nt_perline, group_by);
          putchar('\n');
        }

      /* All-frames translation (-a) */
      if (do_translate)
        {
          for (frame = 0; frame < 3; frame++)
            {
              fpos = frame_start(pos, frame);

              if (! do_topcoord) printf("%*d : ", cw - 1, frame);

              for (i = 0; i < (fpos - pos); i++) putchar(' ');

              while (fpos < pos + nt_perline && fpos + 2 < L)
                {
                  aa = esl_gencode_GetTranslation(gcode, sq->dsq + fpos + 1);  // dsq is 1-based
                  printf("%c  ", aa_abc->sym[aa]);
                  fpos += 3;
                }
              putchar('\n');
            }
        }

      /* Frame-0 triplet translation (-3) */
      if (do_triplet)
        {
          fpos = frame_start(pos, 0);

          if (! do_topcoord) printf("%*s: ", cw, "");

          while (fpos < pos + nt_perline && fpos + 2 < L)
            {
              aa = esl_gencode_GetTranslation(gcode, sq->dsq + fpos + 1);
              printf("%s ", aa3(aa_abc->sym[aa]));
              fpos += 3;
            }
          putchar('\n');
        }

      putchar('\n');
      pos += nprinted;
    }

  /* Verify the file contains exactly one sequence; warn if it didn't.
   */
  {
    sq2    = esl_sq_CreateDigital(nt_abc);
    status = esl_sqio_ReadInfo(sqfp, sq2);
    if      (status == eslOK)      esl_fprintf(stderr, "# warning: %s contains more than one sequence; only the first is shown\n", seqfile);
    else if (status == eslEFORMAT) esl_fatal("Parse failed:\n  %s", esl_sqfile_GetErrorBuf(sqfp));
    else if (status != eslEOF)     esl_fatal("Unexpected read error %d", status);
    esl_sq_Destroy(sq2);
  }



  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_getopts_Destroy(go);
  return eslOK;
}




/* aa3()
 *
 * Returns the three-letter abbreviation for a one-letter amino acid
 * code <c>; for example, "Ala" for 'A'. Returns " * " for the stop
 * codon '*' so that all returned strings are three characters wide,
 * and "???" for any unrecognized symbol.
 *
 * <c> comes from esl_gencode_Translate(), which can only generate
 * the canonical 20 aas or *, but this function handles ambiguous
 * one-letter codes anyway, just to guard against other future uses.
 */
static const char *
aa3(char c)
{
  switch (c) {
  case 'A': return "Ala";
  case 'B': return "N|D";
  case 'C': return "Cys";
  case 'D': return "Asp";
  case 'E': return "Glu";
  case 'F': return "Phe";
  case 'G': return "Gly";
  case 'H': return "His";
  case 'I': return "Ile";
  case 'J': return "I|L";
  case 'K': return "Lys";
  case 'L': return "Leu";
  case 'M': return "Met";
  case 'N': return "Asn";
  case 'O': return "Pyl";
  case 'P': return "Pro";
  case 'Q': return "Gln";
  case 'R': return "Arg";
  case 'S': return "Ser";
  case 'T': return "Thr";
  case 'U': return "Sec";
  case 'V': return "Val";
  case 'W': return "Trp";
  case 'Y': return "Tyr";
  case 'Z': return "Q|E";
  case '*': return " * ";
  default:  return "???";
  }
}


/* coord_width()
 *
 * Returns the number of characters needed to print integer <val> in
 * decimal. Used to size the coordinate columns at the left margin so
 * they align across all lines of output.
 */
static int
coord_width(int64_t val)
{
  int w = 1;
  while (val >= 10) { val /= 10; w++; }
  return w;
}

/* frame_start()
 *
 * Returns the 0-based sequence position of the first codon in reading
 * frame <frame> (0, 1, or 2) that begins at or after position <pos>.
 * Used to align translation lines with the codons on the sequence
 * line above them when a line of output starts in the middle of a
 * codon.
 */
static int64_t
frame_start(int64_t pos, int frame)
{
  if (pos <= frame) return frame;
  return frame + 3 * ((pos - frame + 2) / 3);
}

/* print_top_coords()
 *
 * Prints the header line of coordinates that appears above each
 * sequence line with the `-t` option. Each coordinate is right-
 * justified to the end of its group, so the printed number marks the
 * last residue in that group (counting from <resoffset>+1). At most
 * <nt_perline> residues are accounted for, but the line stops early
 * if the end of the sequence (<L> residues) is reached.
 *
 * <pos> is the 0-based offset of the first residue on this line.
 * <group_by> is the grouping size; this routine assumes <group_by>
 * is nonzero (the `-t` option requires it).
 */
static void
print_top_coords(int64_t pos, int resoffset, int nt_perline, int group_by, int64_t L)
{
  int     nprinted = 0;
  int     gcount   = 0;
  int64_t coord;

  while (nprinted < nt_perline && (pos + nprinted) < L)
    {
      gcount++;
      nprinted++;
      if (gcount == group_by || nprinted == nt_perline || (pos + nprinted) == L)
        {
          coord = pos + nprinted + resoffset;
          printf("%*" PRId64 " ", gcount, coord);
          gcount = 0;
        }
    }
  putchar('\n');
}

/* print_seq_line()
 *
 * Prints up to <nt_perline> residues of the digital sequence <dsq>
 * (1-based, indexed by <pos>+1) starting at 0-based offset <pos>, in
 * uppercase. If <group_by> is nonzero, a single space is inserted
 * after every <group_by> residues, but no trailing space is emitted
 * at end of line or end of sequence. Stops at <L> residues.
 *
 * Returns the number of residues printed (which may be fewer than
 * <nt_perline> at the end of the sequence).
 */
static int
print_seq_line(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int64_t pos, int64_t L, int nt_perline, int group_by)
{
  int nprinted = 0;
  int gcount   = 0;

  while (nprinted < nt_perline && pos + nprinted < L)
    {
      putchar(toupper(abc->sym[dsq[pos + nprinted + 1]]));   // dsq is 1-based
      nprinted++;
      gcount++;
      if (group_by > 0 && gcount == group_by && nprinted < nt_perline && pos + nprinted < L)
        { putchar(' '); gcount = 0; }
    }
  return nprinted;
}

/* print_comp_line()
 *
 * Like print_seq_line(), but prints the complement of each residue,
 * producing the bottom strand of double-stranded output. Returns the
 * number of residues printed.
 */
static int
print_comp_line(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int64_t pos, int64_t L, int nt_perline, int group_by)
{
  int     nprinted = 0;
  int     gcount   = 0;
  ESL_DSQ x;

  while (nprinted < nt_perline && pos + nprinted < L)
    {
      x = abc->complement[dsq[pos + nprinted + 1]];
      putchar(toupper(abc->sym[x]));
      nprinted++;
      gcount++;
      if (group_by > 0 && gcount == group_by && nprinted < nt_perline && pos + nprinted < L)
        { putchar(' '); gcount = 0; }
    }
  return nprinted;
}


