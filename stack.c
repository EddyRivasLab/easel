/* stack.c
 * SRE 1 March 2000 [Seattle]
 * Incorp into Easel SRE, Sun Dec 26 07:42:12 2004 [Zaragoza]
 *
 * Implementations of pushdown stacks for integers, pointers, and characters.
 *
 * Stacks are kept as growable arrays. A stack's memory is
 * grown when necessary by adding some block size. The 
 * initial allocation and the block size are set to 100
 * by default. The block size can be changed by the caller.
 * 
 *****************************************************************
 *
 * Basic API (using integer stack as an example):
 * 
 *   say I want to push the numbers 42, 7, and 3 onto a stack,
 *   then pop them off and print them: 
 *   
 *   #include <easel/easel.h>
 *   #include <easel/stack.h>
 *
 *   ESL_STACK *ns;
 *   int       x;
 *
 *   ns = esl_stack_ICreate();
 *   esl_stack_IPush(ns, 42);
 *   esl_stack_IPush(ns, 7);
 *   esl_stack_IPush(ns, 3);
 *   while (esl_stack_IPop(ns, &x) != ESL_EOD) 
 *      printf("%d\n", x);
 *   esl_stack_Destroy(ns);   
 * 
 * Diagnostics:
 *   Create() functions return NULL on an allocation failure.
 *   Push()   functions throw  ESL_EMEM on an allocation failure.
 *   Pop()    functions return ESL_EOD when the stack is empty.
 *
 * Other functions:
 *   esl_stack_ObjectCount(ns)    :  returns # of objects in the stack
 *   esl_stack_Convert2String(ns) :  converts a char stack to a string, destroying the stack.
 *
 * SVN $Id$
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <easel/easel.h>
#include <easel/stack.h>


/* Function:  esl_stack_ICreate()
 * Incept:    SRE, Sun Dec 26 09:11:50 2004 [Zaragoza]
 *
 * Purpose:   Creates an integer stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_ICreate(void)
{
  ESL_STACK *ns;
  
  if ((ns = ESL_MALLOC(sizeof(ESL_STACK))) == NULL) return NULL;

  ns->nalloc   = ESL_STACK_INITALLOC;
  ns->pdata    = NULL;
  ns->cdata    = NULL;
  if ((ns->idata = ESL_MALLOC(sizeof(int) * ns->nalloc)) == NULL) { free(ns); return NULL; }
  ns->n        = 0;
  return ns;
}

/* Function:  esl_stack_CCreate()
 * Incept:    SRE, Sun Dec 26 09:15:35 2004 [Zaragoza]
 *
 * Purpose:   Creates a character stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_CCreate(void)
{
  ESL_STACK *cs;
  
  if ((cs = ESL_MALLOC(sizeof(ESL_STACK))) == NULL) return NULL;
  cs->nalloc   = ESL_STACK_INITALLOC;
  cs->idata    = NULL;
  cs->pdata    = NULL;
  if ((cs->cdata = ESL_MALLOC(sizeof(char) * cs->nalloc)) == NULL) { free(cs); return NULL; }
  cs->n        = 0;
  return cs;
}

/* Function:  esl_stack_PCreate()
 * Incept:    SRE, Sun Dec 26 09:16:07 2004 [Zaragoza]
 *
 * Purpose:   Creates a pointer stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_PCreate(void)
{
  ESL_STACK *ps;
  
  if ((ps = ESL_MALLOC(sizeof(ESL_STACK))) == NULL) return NULL;

  ps->nalloc   = ESL_STACK_INITALLOC;
  ps->idata    = NULL;
  ps->cdata    = NULL;
  if ((ps->pdata = ESL_MALLOC(sizeof(void *) * ps->nalloc)) == NULL) { free(ps); return NULL; }
  ps->n        = 0;
  return ps;
}

/* Function:  esl_stack_Reuse()
 * Incept:    SRE, Tue Dec 28 04:21:36 2004 [Zaragoza]
 *
 * Purpose:   Empties stack <s> so it can be reused without
 *            creating a new one. The stack <s>
 *            can be of any data type; it retains its original
 *            type.
 *
 * Returns:   <ESL_OK>
 */
int
esl_stack_Reuse(ESL_STACK *s)
{
  s->n = 0;	/* it's that simple in this implementation */
  return ESL_OK;
}

/* Function:  esl_stack_Destroy()
 * Incept:    SRE, Sun Dec 26 09:16:24 2004 [Zaragoza]
 *
 * Purpose:   Destroys a created stack <s>, of any data type.
 */
void
esl_stack_Destroy(ESL_STACK *s)
{
  if (s->idata != NULL) free(s->idata);
  if (s->cdata != NULL) free(s->cdata);
  if (s->pdata != NULL) free(s->pdata);
  free(s);
}

/* Function:  esl_stack_IPush()
 * Incept:    SRE, Sun Dec 26 09:17:17 2004 [Zaragoza]
 *
 * Purpose:   Push an integer <x> onto an integer stack <ns>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EMEM> on reallocation failure.
 */
int
esl_stack_IPush(ESL_STACK *ns, int x)
{
  int *ptr;

  if (ns->n == ns->nalloc) {
    ns->nalloc += ns->nalloc;	/* reallocate by doubling */
    ptr = ESL_REALLOC(ns->idata, sizeof(int) * ns->nalloc);
    if (ptr == NULL) return ESL_EMEM; else ns->idata = ptr;
  }
  ns->idata[ns->n] = x;
  ns->n++;
  return ESL_OK;
}

/* Function:  esl_stack_CPush()
 * Incept:    SRE, Sun Dec 26 09:18:24 2004 [Zaragoza]
 *
 * Purpose:   Push a character <c> onto a character stack <cs>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EMEM> on reallocation failure.
 */
int
esl_stack_CPush(ESL_STACK *cs, char c)
{
  char *ptr;

  if (cs->n == cs->nalloc) {
    cs->nalloc += cs->nalloc;	/* reallocate by doubling */
    ptr = ESL_REALLOC(cs->cdata, sizeof(char) * cs->nalloc);
    if (ptr == NULL) return ESL_EMEM; else cs->cdata = ptr;
  }
  cs->cdata[cs->n] = c;
  cs->n++;
  return ESL_OK;
}

/* Function:  esl_stack_PPush()
 * Incept:    SRE, Sun Dec 26 09:18:49 2004 [Zaragoza]
 *
 * Purpose:   Push a pointer <p> onto a pointer stack <ps>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EMEM> on reallocation failure.
 */
int
esl_stack_PPush(ESL_STACK *ps, void *p)
{
  void *ptr;

  if (ps->n == ps->nalloc) {
    ps->nalloc += ps->nalloc;	/* reallocate by doubling */
    ptr = ESL_REALLOC(ps->pdata, sizeof(void *) * ps->nalloc);
    if (ptr == NULL) return ESL_EMEM; else ps->pdata = ptr;
  }
  ps->pdata[ps->n] = p;
  ps->n++;
  return ESL_OK;
}

/* Function:  esl_stack_IPop()
 * Incept:    SRE, Sun Dec 26 09:19:12 2004 [Zaragoza]
 *
 * Purpose:   Pops an integer off the integer stack <ns>, and returns
 *            it through <ret_x>.
 *
 * Returns:   <ESL_OK> on success. <ESL_EOD> if stack is empty.
 */
int
esl_stack_IPop(ESL_STACK *ns, int *ret_x)
{
  if (ns->n == 0) {*ret_x = 0; return ESL_EOD;}
  ns->n--;
  *ret_x = ns->idata[ns->n];
  return ESL_OK;
}

/* Function:  esl_stack_CPop()
 * Incept:    SRE, Sun Dec 26 09:21:27 2004 [Zaragoza]
 *
 * Purpose:   Pops a character off the character stack <cs>, and returns
 *            it through <ret_c>.
 *
 * Returns:   <ESL_OK> on success. <ESL_EOD> if stack is empty.
 */
int
esl_stack_CPop(ESL_STACK *cs, char *ret_c)
{
  if (cs->n == 0) {*ret_c = 0; return ESL_EOD;}
  cs->n--;
  *ret_c = cs->cdata[cs->n];
  return ESL_OK;
}

/* Function:  esl_stack_PPop()
 * Incept:    SRE, Sun Dec 26 09:21:56 2004 [Zaragoza]
 *
 * Purpose:   Pops a pointer off the pointer stack <ps>, and returns
 *            it through <ret_p>.
 *
 * Returns:   <ESL_OK> on success. <ESL_EOD> if stack is empty.
 */
int
esl_stack_PPop(ESL_STACK *ps, void **ret_p)
{
  if (ps->n == 0) {*ret_p = 0; return ESL_EOD;}
  ps->n--;
  *ret_p = ps->pdata[ps->n];
  return ESL_OK;
}

/* Function:  esl_stack_ObjectCount()
 * Incept:    SRE, Sun Dec 26 09:22:41 2004 [Zaragoza]
 *
 * Purpose:   Returns the number of data objects stored in the
 *            stack <s>. The stack may be of any datatype.
 */
int 
esl_stack_ObjectCount(ESL_STACK *s)
{
  return s->n;
}


/* Function:  esl_stack_Convert2String()
 * Incept:    SRE, Sun Dec 26 09:23:36 2004 [Zaragoza]
 *
 * Purpose:   Converts a character stack <cs> to a NUL-terminated
 *            string, and returns a pointer to the string. The
 *            characters in the string are in the same order they
 *            were pushed onto the stack.  The stack is destroyed by
 *            this operation, as if <esl_stack_Destroy()> had been
 *            called on it. The caller becomes responsible for
 *            free'ing the returned string.
 *
 * Returns:   Pointer to the string; caller must <free()> this.
 *
 * Throws:    NULL if a reallocation fails.
 */
char *
esl_stack_Convert2String(ESL_STACK *cs)
{
  char *s;

  if (esl_stack_CPush(cs, '\0') != ESL_OK)
    { free(cs->cdata); free(cs); return NULL; } /* nul-terminate the data or self-destruct */
  s = cs->cdata;		           /* data is already just a string - just return ptr to it */
  free(cs);			           /* free the stack around it. */
  return s;
}

/* Function:  esl_stack_DiscardTopN()
 * Incept:    SRE, Tue Dec 28 04:33:06 2004 [St. Louis]
 *
 * Purpose:   Throw away the top <n> elements on stack <s>.
 *            Equivalent to <n> calls to a <Pop()> function.
 *            If <n> equals or exceeds the number of elements 
 *            currently in the stack, the stack is emptied
 *            as if <esl_stack_Reuse()> had been called.
 *
 * Returns:   <ESL_OK> on success.
 */
int
esl_stack_DiscardTopN(ESL_STACK *s, int n)
{
  if (n <= s->n) s->n -= n;
  else           s->n = 0;
  return ESL_OK;
}


/*****************************************************************
 * Test driver and API example for the pushdown stack module.
 * To compile:
 *    gcc -g -Wall -I. -DESL_STACK_TESTDRIVE -o test stack.c easel.c
 * To run:
 *    ./test
 * Returns 0 (success) w/ no output, or returns 1 and says why.
 *****************************************************************/
#ifdef ESL_STACK_TESTDRIVE
int 
main(void)
{
  ESL_STACK *s;
  int       *obj;
  int        x;
  char       c;
  char      *str;
  int        n1, n2;
  int        i;

  /* Exercise of integer stacks.
   * 
   * Put 257 integers on the stack and pop them off;
   * do this twice, once with a "while pop" loop, and once
   * with a "while stack not empty" loop.
   *
   * With the default ALLOCSIZE of 128, putting 257 objects on the
   * stack forces two reallocations.
   * 
   */
  /* creation: */
  if ((s = esl_stack_ICreate()) == NULL) 
    { fprintf(stderr, "memory allocation failed\n"); return 1; }

  /* pushing data: */
  n1 = 257;
  for (i = 0; i < n1; i++)
    if (esl_stack_IPush(s, i) != ESL_OK) ESL_ERROR(ESL_EMEM, "push failed");

  /* popping data: */
  n2 = 0;
  while (esl_stack_IPop(s, &x) != ESL_EOD) n2++; 
  if (n1 != n2)
    { fprintf(stderr, "Put %d integers on; got %d off\n", n1, n2); return EXIT_FAILURE;}

  /* same again, but use ObjectCount instead of EOD for popping
   */
  n1 = 257;
  for (i = 0; i < n1; i++)
    if (esl_stack_IPush(s, i) != ESL_OK) ESL_ERROR(ESL_EMEM, "push failed");
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_IPop(s, &x) != ESL_OK) { fprintf(stderr, "pop failed\n"); return 1; }
    n2++; 
  }
  if (n1 != n2)
    { fprintf(stderr, "Put %d objects on; got %d off\n", n1, n2); return 1; }

  /* free the stack.
   */
  esl_stack_Destroy(s);



  /* Exercises of the ptr stack functions/API; same as above,
   * but with pointers.
   */
  if ((s = esl_stack_PCreate()) == NULL) 
    { fprintf(stderr, "memory allocation failed\n"); return 1; }

  n1 = 257;
  for (i = 0; i < n1; i++)
    {
      if ((obj = malloc(sizeof(int) * 64)) == NULL) 
	{ fprintf(stderr, "memory allocation failed\n"); return 1; }

      if (esl_stack_PPush(s, obj) != ESL_OK) 
	{ fprintf(stderr, "esl_stack_PPush failed\n"); return 1; }
    }

  n2 = 0;
  while (esl_stack_PPop(s, (void **) &obj) != ESL_EOD) {
    free(obj); 
    n2++; 
  }
  if (n1 != n2)
    { fprintf(stderr, "Put %d objects on; but popped %d off\n", n1, n2); return 1; }


  n1 = 257;
  for (i = 0; i < n1; i++)
    {
      if ((obj = malloc(sizeof(int) * 64)) == NULL)
	{ fprintf(stderr, "memory allocation failed\n"); return 1; }

      if (esl_stack_PPush(s, obj) != ESL_OK) 
	{ fprintf(stderr, "esl_stack_PPush failed\n"); return 1; }
    }
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_PPop(s, (void **) &obj) == ESL_EOD) 
      { fprintf(stderr, "pop failed\n"); return 1; }
    free(obj); 
    n2++; 
  }
  if (n1 != n2)
    { fprintf(stderr, "Put %d objects on; got %d off\n", n1, n2); return 1; }
  esl_stack_Destroy(s);



  /* Exercise of char stack functions/API;
   * same as above, also including esl_stack_ToString().
   */
  if ((s = esl_stack_CCreate()) == NULL) 
    { fprintf(stderr, "memory allocation failed\n"); return 1; }

  n1 = 257;
  for (i = 0; i < n1; i++)
    if (esl_stack_CPush(s, 'X') != ESL_OK) ESL_ERROR(ESL_EMEM, "push failed");
  n2 = 0;
  while (esl_stack_CPop(s, &c) != ESL_EOD) {
    if (c != 'X') {
      fprintf(stderr, "Put X's on; got a %c off\n", c); 
      return 1;
    }
    n2++; 
  }
  if (n1 != n2)
    { fprintf(stderr, "Put %d characters on; got %d off\n", n1, n2); return 1; }

  for (i = 0; i < n1; i++)
    if (esl_stack_CPush(s, 'X') != ESL_OK) ESL_ERROR(ESL_EMEM, "push failed");
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_CPop(s, &c) != ESL_OK) { fprintf(stderr, "pop failed\n"); return 1; }
    n2++; 
  }
  if (n1 != n2)
    { fprintf(stderr, "Put %d characters on; got %d off\n", n1, n2);  return 1; }

  n1 = 257;
  for (i = 0; i < n1; i++)
    if (esl_stack_CPush(s, 'X') != ESL_OK) ESL_ERROR(ESL_EMEM, "push failed");

  /* Note: the Convert2String() call destroys the stack!
   */
  str = esl_stack_Convert2String(s);
  if ((n2 = strlen(str)) != n1) {
    fprintf(stderr, "thought I'd get %d chars in that string, but got %d\n", n1, n2);
    return 1;
  }
  free(str);

  return EXIT_SUCCESS;
}
#endif /*ESL_STACK_TESTDRIVE*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
