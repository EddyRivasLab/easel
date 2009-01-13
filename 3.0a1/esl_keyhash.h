/* Storing keys in hash tables, similar to Perl's associative arrays.
 * 
 * SRE, Sun Jan 30 08:55:17 2005;  from squid's gki.h, 1999.
 * SVN $Id$
 */
#ifndef eslKEYHASH_INCLUDED
#define eslKEYHASH_INCLUDED

#include <stdio.h>		/* for FILE */

/* ESL_KEYHASH:
 *    a dynamically resized hash structure; 
 *    contains a hash table and associated data
 */
typedef struct {
  int      *hashtable;          /* hashtable[0..hashsize-1] = index of first elem, or -1 */
  uint32_t  hashsize;	        /* size of the hash table                                */

  int      *key_offset;		/* key [idx=0..nkeys-1] starts at smem + key_offset[idx] */
  int      *nxt;		/* nxt [idx=0..nkeys-1], next "pointers" in hash table   */
  int       nkeys;		/* number of keys stored                                 */
  int       kalloc;		/* number of keys allocated for                          */

  char *smem;	 	        /* Array of memory for storing key strings (w/ \0's)     */
  int   salloc;			/* current allocated size of <key_mem>                   */
  int   sn; 			/* current used size of key strings, inclusive \0's      */
} ESL_KEYHASH;

extern ESL_KEYHASH *esl_keyhash_Create(void);
extern ESL_KEYHASH *esl_keyhash_Clone(const ESL_KEYHASH *kh);
extern char *       esl_keyhash_Get(const ESL_KEYHASH *kh, int idx);
extern int          esl_keyhash_GetNumber(const ESL_KEYHASH *kh);
extern void         esl_keyhash_Destroy(ESL_KEYHASH *kh);
extern void         esl_keyhash_Dump(FILE *fp, const ESL_KEYHASH *kh);

extern int  esl_key_Store (      ESL_KEYHASH *kh, const char *key, int *ret_index);
extern int  esl_key_Lookup(const ESL_KEYHASH *kh, const char *key, int *ret_index);


#endif /* eslKEYHASH_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
