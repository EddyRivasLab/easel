/* msa.h
 * Multiple sequence alignment file i/o.
 * 
 * SVN $Id$
 * SRE, Wed Jan 19 19:16:28 2005
 */
#ifndef ESL_MSA_INCLUDED
#define ESL_MSA_INCLUDED

#include <stdio.h>

/* The following constants define the Pfam/Rfam cutoff set we propagate
 * from Stockholm format msa's into HMMER and Infernal models.
 */
#define eslMSA_CUTOFF_TC1 0
#define eslMSA_CUTOFF_TC2 1
#define eslMSA_CUTOFF_GA1 2
#define eslMSA_CUTOFF_GA2 3
#define eslMSA_CUTOFF_NC1 4
#define eslMSA_CUTOFF_NC2 5
#define eslMSA_MAXCUTS    6

/* Object: ESL_MSA
 * 
 * A multiple sequence alignment.
 */
typedef struct {
  /* Mandatory information associated with the alignment.
   * (The important stuff.)
   */
  char **aseq;                  /* alignment itself, [0..nseq-1][0..alen-1] */
  char **sqname;                /* sequence names, [0..nseq-1][]            */
  float *wgt;	                /* sequence weights [0..nseq-1]             */
  int    alen;			/* length of alignment (columns)            */
  int    nseq;			/* number of seqs in alignment              */

  /* Optional information that we understand, and might have.
   * (The occasionally useful stuff.)
   */
  int    flags;			/* flags for what optional info is valid    */
  int    type;			/* eslRNA, eslDNA, eslAMINO, eslNONSTANDARD */
  char  *name;             	/* name of alignment, or NULL               */
  char  *desc;	                /* description of alignment, or NULL        */
  char  *acc;	                /* accession of alignment, or NULL          */
  char  *au;		        /* "author" information, or NULL            */
  char  *ss_cons;		/* consensus secondary structure, or NULL   */
  char  *sa_cons;               /* consensus surface accessibility, or NULL */
  char  *rf;                    /* reference coordinate system, or NULL     */
  char **sqacc;			/* accession numbers for sequences i        */
  char **sqdesc;		/* description lines for sequences i        */
  char **ss;                    /* per-seq secondary structures, or NULL    */
  char **sa;                    /* per-seq surface accessibilities, or NULL */
  float  cutoff[eslMSA_MAXCUTS];/* NC/TC/GA cutoffs propagated to Pfam/Rfam */
  int    cutset[eslMSA_MAXCUTS];/* TRUE if a cutoff is set; else FALSE      */

  /* Optional information that we don't understand.
   * (The stuff we're morally obligated to parse because it's
   *  in a file's markup tags, but all we know how to do
   *  is regurgitate it.)
   *
   * That is, we know what type of information it is, but it's
   * either (interpreted as) free-text comment, or it's Stockholm 
   * markup with unfamiliar tags.
   */
  char  **comment;              /* free text comments, or NULL      */
  int     ncomment;		/* number of comment lines          */
  int     alloc_ncomment;	/* number of comment lines alloc'ed */

  char  **gf_tag;               /* markup tags for unparsed #=GF lines  */
  char  **gf;                   /* annotations for unparsed #=GF lines  */
  int     ngf;			/* number of unparsed #=GF lines        */
  int     alloc_ngf;		/* number of gf lines alloc'ed          */

  char  **gs_tag;               /* markup tags for unparsed #=GS lines     */
  char ***gs;                   /* [0..ngs-1][0..nseq-1][free text] markup */
  GKI    *gs_idx;               /* hash of #=GS tag types                  */
  int     ngs;                  /* number of #=GS tag types                */
  
  char  **gc_tag;               /* markup tags for unparsed #=GC lines  */
  char  **gc;                   /* [0..ngc-1][0..alen-1] markup         */
  GKI    *gc_idx;               /* hash of #=GC tag types               */
  int     ngc;                  /* number of #=GC tag types             */

  char  **gr_tag;               /* markup tags for unparsed #=GR lines   */
  char ***gr;                   /* [0..ngr][0..nseq-1][0..alen-1] markup */
  GKI    *gr_idx;               /* hash of #=GR tag types                */
  int     ngr;			/* number of #=GR tag types              */

  /* Stuff we need for our own maintenance of the data structure
   */
  GKI   *index;		        /* name ->seqidx hash table */
  int    nseqalloc;		/* number of seqs currently allocated for   */
  int    nseqlump;		/* lump size for dynamic expansions of nseq */
  int   *sqlen;                 /* individual seq lengths during parsing    */
  int   *sslen;                 /* individual ss lengths during parsing     */
  int   *salen;                 /* individual sa lengths during parsing     */
  int    lastidx;		/* last index we saw; use for guessing next */
} MSA;

/* Flags for msa->flags
 */
#define eslMSA_SET_WGT  (1 << 0)  /* 1 if wgts were set, 0 if default 1.0's */


                                     
/* Object: MSAFILE
 * 
 * Defines an alignment file that's open for reading.
 */
typedef struct {
  FILE *f;                      /* open file pointer                         */
  char *fname;			/* name of file. used for diagnostic output  */
  int   linenumber;		/* what line are we on in the file           */

  char *buf;			/* buffer for line input w/ sre_fgets()      */
  int   buflen;			/* current allocated length for buf          */

  int   do_gzip;		/* TRUE if f is "gzip -dc |" (will pclose(f))*/
  int   do_stdin;		/* TRUE if f is stdin (won't close f)        */
  int   format;			/* format of alignment file we're reading    */

#ifdef ESL_SSI_INCLUDED		/* AUGMENTATION: SSI indexing of an MSA db   */
  SSIFILE *ssi;		        /* open SSI index file; or NULL, if none.    */
#endif
} MSAFILE;


/* Alignment file format codes.
 * Must coexist with sqio.c/squid.h unaligned file format codes.
 * Rules:
 *     - 0 is an unknown/unassigned format 
 *     - <100 reserved for unaligned formats
 *     - >100 reserved for aligned formats
 */
#define eslMSAFILE_UNKNOWN   0	  /* unknown format                          */
#define eslMSAFILE_STOCKHOLM 101  /* Pfam/Rfam Stockholm format              */
#define eslMSAFILE_SELEX     102  /* Obsolete(!): old SELEX format           */
#define eslMSAFILE_MSF	     103  /* GCG MSF format                          */
#define eslMSAFILE_CLUSTAL   104  /* Clustal V/W format                      */
#define eslMSAFILE_A2M	     105  /* aligned FASTA (A2M is UCSC terminology) */
#define eslMSAFILE_PHYLIP    106  /* Felsenstein's PHYLIP format             */
#define eslMSAFILE_EPS       107  /* Encapsulated PostScript (output only)   */

#define IsAlignmentFormat(fmt)  ((fmt) > 100)


/* from msa.c
 */
extern MSAFILE *MSAFileOpen(char *filename, int format, char *env);
extern MSA     *MSAFileRead(MSAFILE *afp);
extern void     MSAFileClose(MSAFILE *afp);
extern void     MSAFree(MSA *msa);
extern void     MSAFileWrite(FILE *fp, MSA *msa, int outfmt, int do_oneline);

extern int MSAFileRewind(MSAFILE *afp);
extern int MSAFilePositionByKey(MSAFILE *afp, char *key);
extern int MSAFilePositionByIndex(MSAFILE *afp, int idx);

extern int   MSAFileFormat(MSAFILE *afp);
extern MSA  *MSAAlloc(int nseq, int alen);
extern void  MSAExpand(MSA *msa);
extern char *MSAFileGetLine(MSAFILE *afp);
extern void  MSASetSeqAccession(MSA *msa, int seqidx, char *acc);
extern void  MSASetSeqDescription(MSA *msa, int seqidx, char *desc);
extern void  MSAAddComment(MSA *msa, char *s);
extern void  MSAAddGF(MSA *msa, char *tag, char *value);
extern void  MSAAddGS(MSA *msa, char *tag, int seqidx, char *value);
extern void  MSAAppendGC(MSA *msa, char *tag, char *value);
extern char *MSAGetGC(MSA *msa, char *tag);
extern void  MSAAppendGR(MSA *msa, char *tag, int seqidx, char *value);
extern void  MSAVerifyParse(MSA *msa);
extern int   MSAGetSeqidx(MSA *msa, char *name, int guess);

extern MSA  *MSAFromAINFO(char **aseq, AINFO *ainfo);   

extern void  MSAMingap(MSA *msa);
extern void  MSANogap(MSA *msa);
extern void  MSAShorterAlignment(MSA *msa, int *useme);
extern void  MSASmallerAlignment(MSA *msa, int *useme, MSA **ret_new);

extern char *MSAGetSeqAccession(MSA *msa, int idx);
extern char *MSAGetSeqDescription(MSA *msa, int idx);
extern char *MSAGetSeqSS(MSA *msa, int idx);
extern char *MSAGetSeqSA(MSA *msa, int idx);

extern float MSAAverageSequenceLength(MSA *msa);

/* from a2m.c
 */
extern MSA  *ReadA2M(MSAFILE *afp);
extern void  WriteA2M(FILE *fp, MSA *msa);

/* from clustal.c
 */
extern MSA  *ReadClustal(MSAFILE *afp);
extern void  WriteClustal(FILE *fp, MSA *msa);

/* from eps.c
 */
extern void EPSWriteSmallMSA(FILE *fp, MSA *msa);

/* from msf.c
 */
extern MSA  *ReadMSF(MSAFILE *afp);
extern void  WriteMSF(FILE *fp, MSA *msa);

/* from phylip.c
 */
extern MSA  *ReadPhylip(MSAFILE *afp);
extern void  WritePhylip(FILE *fp, MSA *msa);

/* from selex.c
 */
extern MSA  *ReadSELEX(MSAFILE *afp);
extern void  WriteSELEX(FILE *fp, MSA *msa);
extern void  WriteSELEXOneBlock(FILE *fp, MSA *msa);

/* from stockholm.c
 */
extern MSA  *ReadStockholm(MSAFILE *afp);
extern void  WriteStockholm(FILE *fp, MSA *msa);
extern void  WriteStockholmOneBlock(FILE *fp, MSA *msa);

#endif /*ESL_MSA_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
