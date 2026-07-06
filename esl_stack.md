# esl_stack - pushdown stacks for integers, chars, and pointers

The `stack` module implements pushdown stacks for storing integers,
characters, or arbitrary pointers (objects).

The module uses a convention of prepending `I`, `C`, `P` to
`Create()`, `Push()`, and `Pop()` function names, to indicate the
stack's datatype as integer, character, or pointer, respectively. For
example, `esl_stack_PCreate()` creates a stack for pointer storage.
(This is also the naming convention in the `vectorops` module.) Stacks
can be thought of as typed objects, with the type defined by the
`Create()` call. Types may not be mixed for any particular created
stack.


## example

The `eslSTACK_EXAMPLE` example at the end of `esl_stack.c`
uses an integer stack, pushes 42, 7, and 3 on, then pops them
off.

The `Create()` functions create a stack for a particular purpose.
`esl_stack_ICreate()` creates a stack for integers,
`esl_stack_CCreate()` creates a stack for characters, and
`esl_stack_PCreate()` creates a stack for pointers. They throw NULL if
an allocation fails. All three stack types are free'd by a call to
`esl_stack_Destroy()`. All three types can also be reused without
reallocation or recreation by `esl_stack_Reuse()`. A `Reuse()`'d stack
retains its original datatype.

The `Push()` functions push one datum onto the stack, of the
appropriate type. They throw `eslEMEM` if the stack needs to
reallocate internally but fails.

The `Pop()` functions pop one datum off the stack, returning it
through a passed pointer. They return `eslOK` on success, and `eslEOD`
if the stack is empty.

`esl_stack_ObjectCount()` returns the number of objects stored in the
stack.

A special function, `esl_stack_Convert2String()`, operates only on
character stacks. It converts the stack structure to a NUL-terminated
string, with characters in the same order they were pushed. The stack
is destroyed by this operation, leaving a `char *` behind.


## allocation strategy

Stacks are initially allocated for a certain number of objects,
defined by a compile-time constant `ESL_STACK_INITALLOC` in
`esl_stack.h`. The default is 128. Whenever a stack needs to grow, it
reallocates by doubling its current allocation.
