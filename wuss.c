/* wuss.c
 * RNA secondary structure markup in WUSS notation.
 * 
 * SRE, Tue Feb 15 08:43:23 2005
 * SVN $Id$
 * xref squid wuss.c.
 */

#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/stack.h>
#include <easel/wuss.h>

/* Function:  esl_wuss2ct()
 * Incept:    SRE, Tue Feb 15 08:44:54 2005 [St. Louis]
 *
 * Purpose:   Given a secondary structure string <ss>, <0..len-1>,
 *            in WUSS notation, convert it to a CT array, <1..len>,
 *            in <ct>. Caller provides a <ct> allocated for at least 
 *            <len+1> ints. <ct[i]> is the position that residue i
 *            base pairs to, or 0 if i is unpaired. <ct[0]> is undefined
 *            (but if you care: it is set to 0).
 *            
 *            WUSS notation is interpreted loosely here, as input
 *            WUSS.  Any matching bracket pair or upper/lower case
 *            alphabetic pair is interpreted as a base pair; any other
 *            WUSS annotation is interpreted as unpaired.
 *            
 * Returns:   <eslOK> on success. Returns <eslESYNTAX> if the WUSS
 *            string isn't valid.
 */
int 
esl_wuss2ct(char *ss, int len, int *ct)
{
  ESL_STACK *pda[27];     /* 1 secondary structure + up to 26 levels of pk's */
  int        i;
  int        pos, pair;
  int        status;      /* success or failure return status */

 /* Initialization: always initialize the main pda (0);
  * we'll init the pk pda's on demand.
  */
  pda[0] = esl_stack_ICreate();
  for (i = 1; i <= 26; i++) pda[i] = NULL;

  for (pos = 0; pos <= len; pos++) ct[pos] = 0;

  status = eslOK;
  for (pos = 1; pos <= len; pos++)
    {
      if (!isprint((int) ss[pos-1]))  /* armor against garbage */
	status = eslESYNTAX;

      /* left side of a pair: push position onto stack 0 (pos = 1..L) */
      else if (ss[pos-1] == '<' ||
	       ss[pos-1] == '(' ||
	       ss[pos-1] == '[' ||
	       ss[pos-1] == '{')
        esl_stack_IPush(pda[0], pos);
      
      /* right side of a pair; resolve pair; check for agreement */
      else if (ss[pos-1] == '>' || 
	       ss[pos-1] == ')' ||
	       ss[pos-1] == ']' ||
	       ss[pos-1] == '}')
        {
          if (esl_stack_IPop(pda[0], &pair) == eslEOD)
            { status = eslESYNTAX; }	/* no closing bracket */
          else if ((ss[pair-1] == '<' && ss[pos-1] != '>') ||
		   (ss[pair-1] == '(' && ss[pos-1] != ')') ||
		   (ss[pair-1] == '[' && ss[pos-1] != ']') ||
		   (ss[pair-1] == '{' && ss[pos-1] != '}'))
	    { status = eslESYNTAX; }	/* brackets don't match */
	  else
	    {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
                                /* same stuff for pseudoknots */
      else if (isupper((int) ss[pos-1])) 
	{
	  /* Create the PK stacks on demand.
	   */
	  i = ss[pos-1] - 'A' + 1;
	  if (pda[i] == NULL) pda[i] = esl_stack_ICreate();
	  esl_stack_IPush(pda[i], pos);
	}
      else if (islower((int) ss[pos-1])) 
	{
	  i = ss[pos-1] - 'a' + 1;
	  if (pda[i] == NULL || 
	      esl_stack_IPop(pda[i], &pair) == eslEOD)
            { status = eslESYNTAX; }
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
	}
      else if (strchr(":,_-.~", ss[pos-1]) == NULL)
	status = eslESYNTAX; /* bogus character */
    }
                                
  for (i = 0; i <= 26; i++)
    if (pda[i] != NULL) 
      { /* nothing should be left on stacks */
	if (esl_stack_ObjectCount(pda[i]) != 0)
	  status = eslESYNTAX;
	esl_stack_Destroy(pda[i]);
      }

  return status;
}


/* Function:  esl_wuss2kh()
 * Incept:    SRE, Tue Feb 15 10:05:35 2005 [St. Louis]
 *
 * Purpose:   Converts a secondary structure string <ss> in 
 *            WUSS notation back to old KHS format in <kh>.
 *            <kh> must be allocated for at least as much
 *            space as <ss>. <kh> may be the same as <ss>,
 *            in which case the conversion is done in-place.
 *
 * Note:      Left bp chars  are converted to >   (left base of base pairs)
 *            Right bp chars are converted to <   (right base of base pairs)
 *            Characters _-,:~ are converted to . (unpaired bases)
 *            Character  .     is untouched       (unpaired)
 *            Everything else is untouched, including any pseudoknot notation.
 * 
 * Returns:   <eslOK> on success.
 */
int
esl_wuss2kh(char *ss, char *kh)
{
  while (*ss != '\0')
    {
      if       (*ss == '<') *kh = '>';
      else if  (*ss == '(') *kh = '>';
      else if  (*ss == '[') *kh = '>';
      else if  (*ss == '{') *kh = '>';
      else if  (*ss == '>') *kh = '<';
      else if  (*ss == ')') *kh = '<';
      else if  (*ss == ']') *kh = '<';
      else if  (*ss == '}') *kh = '<';
      else if  (*ss == '_') *kh = '.';
      else if  (*ss == '-') *kh = '.';
      else if  (*ss == ',') *kh = '.';
      else if  (*ss == ':') *kh = '.';
      else if  (*ss == '~') *kh = '.';
      else *kh = *ss;
      ss++;
      kh++;
    }
  *kh = '\0';
  return eslOK;
}


/* Function:  esl_kh2wuss()
 * Incept:    SRE, Tue Feb 15 10:10:40 2005 [St. Louis]
 *
 * Purpose:   Converts an old format secondary structure string <kh>
 *            to shorthand WUSS format <ss>. <ss> must be allocated at least
 *            as large as <kh>. <ss> can be identical to <kh>, in which
 *            case the conversion is done in-place.
 *
 * Note:      Character > is converted to <  (left base of base pairs)
 *            Character < is converted to >  (right base of base pairs)
 *            A space is converted to .      (just in case)      
 *
 * Returns:   <eslOK> on success.
 */
int
esl_kh2wuss(char *kh, char *ss)
{
  while (*kh != '\0')
    {
      if      (*kh == '>') *ss = '<';
      else if (*kh == '<') *ss = '>';
      else if (*kh == ' ') *ss = '.';
      else *ss = *kh;
      kh++;
      ss++;
    }
  *ss = '\0';
  return eslOK;
}


/* Function:  esl_wuss_nopseudo()
 * Incept:    SRE, Tue Feb 15 11:02:43 2005 [St. Louis]
 *
 * Purpose:   Given a WUSS format annotation string <ss1>,
 *            removes all pseudoknot annotation to create a new 
 *            WUSS string <ss2> that contains only a "canonical"
 *            (nonpseudoknotted) structure. <ss2> must be allocated to
 *            be at least as large as <ss1>. <ss1> and <ss2>
 *            may be the same, in which case the conversion is
 *            done in place. Pseudoknot annotation in <ss1> is
 *            simply replaced by <.> in <ss2>; the resulting
 *            <ss2> WUSS string is therefore in valid simplified format,
 *            but may not be valid full format WUSS.
 *
 * Returns:   <eslOK>.
 */
int
esl_wuss_nopseudo(char *ss1, char *ss2)
{
  while (*ss1 != '\0') 
    {
      if (isalpha(*ss1)) *ss2 = '.';
      else *ss2 = *ss1;
      ss1++;
      ss2++;
    }
  *ss2 = '\0';
  return eslOK;
}


#ifdef eslWUSS_TESTDRIVE
/* gcc -g -Wall -o test -I. -DeslWUSS_TESTDRIVE wuss.c stack.c easel.c
 * ./test
 */
#include <stdlib.h>

#include <easel/easel.h>
#include <easel/stack.h>
#include <easel/wuss.h>

int
main(int argc, char **argv)
{
  char ss[] = "<<<<...AAAA>>>>aaaa...";
  int  len;
  int *ct1;
  int *ct2;
  int  i;
  int  nbp;

  len = strlen(ss);
  ESL_MALLOC(ct1, sizeof(int) * (len+1));
  ESL_MALLOC(ct2, sizeof(int) * (len+1));

  esl_wuss2ct(ss, len, ct1);
  nbp = 0;
  for (i = 1; i <= len; i++)
    if (ct1[i] > i) nbp++;
  if (nbp != 8) abort();

  esl_wuss2kh(ss, ss);
  esl_kh2wuss(ss, ss);
  esl_wuss2ct(ss, len, ct2);
  
  for (i = 1; i <= len; i++)
    if (ct1[i] != ct2[i]) abort();

  esl_wuss_nopseudo(ss, ss);
  esl_wuss2ct(ss, len, ct1);
  nbp = 0;
  for (i = 1; i <= len; i++)
    if (ct1[i] > i) nbp++;
  if (nbp != 4) abort();

  return 0;
}
#endif /*eslWUSS_TESTDRIVE*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
