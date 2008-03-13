/* Storing keys in hash tables, similar to Perl's associative arrays.
 * 
 * SRE, Sun Jan 30 08:55:17 2005;  from squid's gki.h, 1999.
 * SVN $Id$
 */
#ifndef eslKEYHASH_INCLUDED
#define eslKEYHASH_INCLUDED

#include <stdio.h>		/* for FILE */

/* esl_key_elem:
 *    key, array index pairs are kept in linked list structures.
 */
struct esl_key_elem {
  char *key;		
  int   idx;
  struct esl_key_elem *nxt;
};

/* ESL_KEYHASH:
 *    a dynamically resized hash structure; 
 *    contains a hash table and associated data
 */
typedef struct {
  struct esl_key_elem **table;	  /* hash table; heads of linked lists, [0..nhash-1]   */
  int primelevel;		  /* 0..eslKEY_NPRIMES-1; which hash table size we are */
  int nhash;			  /* redundant with key_primes[primelevel]             */
  int nkeys;			  /* number of keys stored in the <table>              */

  struct esl_key_elem *table_mem; /* Pool of memory for storing key_elem's             */
  int                  talloc;	  /* current allocated # in <table_mem>                */

  char *key_mem;		  /* Pool of memory for storing key strings (w/ \0's)  */
  int   kalloc;			  /* current allocated size of <key_mem>               */
  int   kn;			  /* current used size of key strings, inclusive \0's  */
} ESL_KEYHASH;

extern ESL_KEYHASH *esl_keyhash_Create(void);
extern ESL_KEYHASH *esl_keyhash_Clone(const ESL_KEYHASH *kh);
extern void         esl_keyhash_Destroy(ESL_KEYHASH *kh);
extern void         esl_keyhash_Dump(FILE *fp, const ESL_KEYHASH *kh);

extern int  esl_key_Store (      ESL_KEYHASH *kh, const char *key, int *ret_index);
extern int  esl_key_Lookup(const ESL_KEYHASH *kh, const char *key, int *ret_index);


#endif /* eslKEYHASH_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
