#ifndef eslBITFIELD_INCLUDED
#define eslBITFIELD_INCLUDED
#include "esl_config.h"

#include "easel.h"

typedef uint64_t ESL_BITFIELD;

static inline void
esl_bitfield_Set(ESL_BITFIELD *b, int i)
{
  b[i/64] |= (1ull << (i%64));
}

static inline void
esl_bitfield_Clear(ESL_BITFIELD *b, int i)
{
  b[i/64] &= ~(1ull << (i%64));
}

static inline void
esl_bitfield_Toggle(ESL_BITFIELD *b, int i)
{
  b[i/64] ^= (1ull << (i%64));
}

static inline int
esl_bitfield_IsSet(ESL_BITFIELD *b, int i)
{
  return ((b[i/64] & (1ull << (i%64))) ? TRUE : FALSE);
}


extern ESL_BITFIELD *esl_bitfield_Create(int nb);
extern void          esl_bitfield_Destroy(ESL_BITFIELD *b);

#endif //eslBITFIELD_INCLUDED



