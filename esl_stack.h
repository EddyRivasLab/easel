/* Pushdown stacks for integers, pointers, and characters.
 *
 * nstack - SRE 1 March 2000. [Seattle]
 * mstack - SRE, Fri Oct 10 10:18:16 2003 [St. Louis]
 * cstack - SRE, Mon Oct 13 12:57:56 2003 [St. Louis]
 * Incorp into easel - SRE, Sun Dec 26 07:39:02 2004 [Zaragoza]
 * SVN $Id$
 */
#ifndef ESL_STACK_INCLUDED
#define ESL_STACK_INCLUDED

#define ESL_STACK_INITALLOC 128	/* initial allocation; realloc by doubling  */

typedef struct esl_stack_s {
  int   *idata;			/* integer data stack                       */
  void **pdata;			/* pointer data stack                       */
  char  *cdata;			/* character data stack                     */

  int  n;			/* current (topmost) elem in data           */
  int  nalloc;			/* # of elems allocated right now           */
} ESL_STACK;

extern ESL_STACK *esl_stack_ICreate(void);
extern ESL_STACK *esl_stack_CCreate(void);
extern ESL_STACK *esl_stack_PCreate(void);

extern int        esl_stack_Reuse(ESL_STACK *s);
extern void       esl_stack_Destroy(ESL_STACK *s);

extern int esl_stack_IPush(ESL_STACK *ns, int x);
extern int esl_stack_CPush(ESL_STACK *cs, char c);
extern int esl_stack_PPush(ESL_STACK *ps, void *p);

extern int esl_stack_IPop(ESL_STACK *ns, int   *ret_x);
extern int esl_stack_CPop(ESL_STACK *cs, char  *ret_c);
extern int esl_stack_PPop(ESL_STACK *ps, void **ret_p);

extern int esl_stack_ObjectCount(ESL_STACK *s);

extern char *esl_stack_Convert2String(ESL_STACK *cs);
extern int   esl_stack_DiscardTopN(ESL_STACK *s, int n);

#endif /*ESL_STACK_INCLUDED*/
