/* Vectorized routines for PowerPC, using Altivec.
 * 
 */
#ifdef HAVE_VMX
#ifndef ESL_VMX_INCLUDED
#define ESL_VMX_INCLUDED

#include "easel.h"

#include <stdio.h>
#include <altivec.h>


extern vector float esl_vmx_logf(vector float x);
extern vector float esl_vmx_expf(vector float x);
extern void         esl_vmx_dump_vecfloat(FILE *fp, vector float v);

#endif /*ESL_VMX_INCLUDED*/
#endif /*HAVE_VMX*/
