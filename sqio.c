/* sqio.c
 * Sequence file i/o.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/parse.h>
#include <easel/sqio.h>

ESL_SQ *
esl_sq_Create(void)
{
  ESL_SQ *sq;
  int     status;

  sq = malloc(sizeof(ESL_SQ));
  if (sq == NULL) ESL_ERROR_VAL(NULL, ESL_EMEM, "malloc failed");
  status = esl_sq_Init(sq);
  if (status != ESL_OK) { free(sq); return NULL; }
  return sq;
}




int
esl_sq_Init(ESL_SQ *sq)
{
  sq->seq = malloc(sizeof(char) * ESL_SQ_CHUNKSIZE);
  if (sq->seq == NULL) ESL_ERROR(ESL_EMEM, "malloc failed");
  sq->ss  = NULL;
  sqfp->alloclen = ESL_SQ_CHUNKSIZE;
  esl_sq_Reuse(sq);
  return ESL_OK;
}

int
esl_sq_Reuse(ESL_SQ *sq)
{
  sqfp->name[0] = '\0';
  sqfp->acc[0]  = '\0';
  sqfp->desc[0] = '\0';
  if (sqfp->seq != NULL) sqfp->seq[0]  = '\0';
  if (sqfp->ss  != NULL) sqfp->ss[0] = '\0';
  sqfp->len = 0;
  return ESL_OK;
}

int 
esl_sq_Release(ESL_SQ *sq)
{
  if (sq->seq != NULL) free(sq->seq);
  if (sq->ss  != NULL) free(sq->ss);
  sq->alloclen = 0;
  return ESL_OK;
}

int
esl_sq_Destroy(ESL_SQ *sq)
{
  esl_sq_Fini(sq);
  free(sq);
  return ESL_OK;
}



int
esl_sqio_OpenFASTA(char *seqfile, ESL_SQFILE *sqfp)
{
  int n;
  int status;

  sqfp->fp = fopen(seqfile,"r");
  if (sqfp->fp == NULL) ESL_ERROR(ESL_ENOFILE, "file not found");

  n = strlen(seqfile);
  sqfp->filename = malloc(sizeof(char) * (n+1));
  if (sqfp->filename == NULL) ESL_ERROR(ESL_EMEM, "malloc failed");
  strcpy(sqfp->filename, seqfile);

  sqfp->format     = ESL_SQFORMAT_FASTA;
  sqfp->do_gzip    = FALSE;
  sqfp->do_stdin   = FALSE;

  sqfp->linenumber = 0;
  sqfp->buf        = NULL;
  sqfp->buflen     = 0;

  if ((status = esl_parse_fgets(&buf, &n, sqfp->fp)) != ESL_OK)
    {
      esl_sqio_CloseFASTA(sqfp);
      return status;
    }
}

int
esl_sqio_ReadFASTA(ESL_SQFILE *sqfp, ESL_SQ *s)
{
  char *tok, *sptr;
  int   len;
  int   status;
  
  /* The lookahead buffer should contain a FASTA descline.
   */
  if (*sqfp->buf != '>') ESL_ERROR(ESL_EFORMAT, "file not in FASTA format");
  
  /* Parse the descline.
   */
  /* First token after the > is the name. */
  sptr = sqfp->buf + 1;
  status = esl_parse_strtok(&sptr, " \t\n", &tok, &len);
  if (status != ESL_OK) return status;
  if (tok == NULL)       ESL_ERROR(ESL_EFORMAT, "no sequence name found");
  if (len >= SQ_MAXNAME) ESL_ERROR(ESL_EFORMAT, "sequence name too long");
  strcpy(s->name, tok);

  /* Remainder of the line is the description. */
  while (isspace(*sptr)) sptr++;
  if (*sptr != '\0') {
    strncpy(s->desc, sptr, SQ_MAXDESC-1);
    s->desc[SQ_MAXDESC-1] = '\0';
  }
  
  /* Parse the rest of the sequence, one line at a time.
   */
  s->n = 0;
  while ((status = sqfile_getnextline(sqfp)) == ESL_OK) 
    {
      if (*sqfp->buf == '>') break; /* we've reached the next sequence: so stop */
      addseq(s, sqfp->buf);
    }
  if (status == ESL_EMEM) return status; /* we're out of memory, couldn't get seq. */
  s->seq[s->len] = '\0';
  return status;
}


int
esl_sqio_CloseFASTA()
{}




static int
getnextline(ESL_SQFILE *sqfp)
{
  int status;

  sqfp->linenumber++;
  status = esl_parse_fgets(&sqfp->buf, &sqfp->buflen, sqfp->fp);
  return status; /* ESL_OK, ESL_EOF, ESL_EMEM */ 
}

static int
addseq(SQ *sq, char *line)
{
  char *s;
  int   n;

  n = sq->len;
  for (s = line; *s != '\0'; s++)
    {
      if (! isalpha(*s)) continue; /* accept any alphabetic character */
      
      sq->seq[n] = *s;		   /* store the character */
      n++;
      
      if (n == sq->alloclen)	   /* out of room? then realloc. */
	{
	  sq->alloclen += sq->allocchunk;
	  sq->seq = realloc(sq->seq, sizeof(char) * sq->alloclen);
	  if (seq->seq == NULL) ESL_ERROR(ESL_EMEM, "realloc failed");
	}
    }
  sq->len = n;
  /* note that caller is responsible for null-termination of the growing seq. */
  return ESL_OK;
}
