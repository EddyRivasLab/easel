/* keyhash.c [gki]
 * "Generic key index" module: emulation of Perl hashes.
 * Maps keys (ASCII char strings) to array index. Dynamically
 * resizes the hash table. 
 * 
 * SVN $Id$
 * SRE, Sun Jan 30 09:14:21 2005
 * From squid's gki.c, 1999.
 *****************************************************************
 * 
 * Limitations:
 *     - hash table can only grow; no provision for deleting keys
 *       or downsizing the hash table.
 *     - Maximum hash table size set at 100003. Performance 
 *       will degrade for key sets much larger than this.
 *     - Assumes that integers are 32 bits (or greater). 
 * 
 *****************************************************************
 * 
 * API for storing/reading stuff: 
 * moral equivalent of Perl's $foo{$key} = whatever, $bar{$key} = whatever:
 *       #include <easel/easel.h>
 *       #include <easel/keyhash.h>
 *     
 *       ESL_KEYHASH  *hash;
 *       int   idx;
 *       char *key;
 *       
 *       hash = esl_keyhash_Create();
 * (Storing:) 
 *       (foreach key) {
 *          esl_gki_Store(hash, key, &idx);       
 *          (reallocate foo, bar as needed)
 *          foo[idx] = whatever;
 *          bar[idx] = whatever;
 *       }     
 * (Reading:)
 *       (foreach key) {
 *          idx = esl_gki_Lookup(hash, key);
 *          if (idx == -1) {no_such_key; }
 *          (do something with) foo[idx];
 *          (do something with) bar[idx];
 *       }   
 *       esl_keyhash_Destroy();
 *       
 *****************************************************************
 *
 * Timings on wrasse for 45402 keys in /usr/dict/words using
 * Tests/test_gki: 
 *      250 msec store      (6 usec/store)
 *      140 msec retrieve   (3 usec/retrieve)
 * and using the 13408 names of Pfam's GP120.full alignment:
 *       70 msec store      (5 usec/store)
 *       50 msec retrieve   (4 usec/retrieve)     
 * 
 * CVS $Id: gki.c 878 2003-04-14 16:00:17Z eddy $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <easel/easel.h>
#include <easel/keyhash.h>

/* gki_primes[] defines the ascending order of hash table sizes
 * that we use in upsizing the hash table dynamically.
 *
 * Best hash table sizes are prime numbers (see Knuth vol 3, Sorting
 * and Searching). Useful site for testing primes:
 *   http://www.idbsu.edu/people/jbrennan/algebra/numbers/sieve.html
 *
 * Because of the way gki_hashvalue works, the largest number
 * must be < INT_MAX / 128 / 128 : 131072 on a 32 bit machine.
 */
static int gki_primes[]  = { 101, 1009, 10007, 100003 };
#define eslGKI_NPRIMES      4
#define eslGKI_ALPHABETSIZE 128

static ESL_KEYHASH *gki_alloc(int primelevel);
static int          gki_hashvalue(ESL_KEYHASH *hash, char *key);
static int          gki_upsize(ESL_KEYHASH *old);

/* Function:  esl_keyhash_Create()
 * Incept:    SRE, Sun Jan 30 09:17:20 2005 [St. Louis]
 *
 * Purpose:   Create a hash table for key indexing.  
 *           
 * Note:      A wrapper around a level 0 gki_alloc().
 *
 * Returns:   a new <hash>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_KEYHASH *
esl_keyhash_Create(void)
{
  ESL_KEYHASH *hash;
  hash = gki_alloc(0);
  return hash;
}

/* Function:  esl_keyhash_Destroy()
 * Incept:    SRE, Sun Jan 30 09:19:19 2005 [St. Louis]
 *
 * Purpose:   Destroys a hash table <hash>.
 *
 * Returns:   (void)
 */
void
esl_keyhash_Destroy(ESL_KEYHASH *hash)
{
  struct esl_gki_elem *ptr;
  int i;

  if (hash == NULL) return;	/* tolerate a NULL */

  for (i = 0; i < hash->nhash; i++)
    while (hash->table[i] != NULL)
      {
	ptr = hash->table[i]->nxt;
			/* NULL keys can occur after we've gki_upsize'd */
	if (hash->table[i]->key != NULL) free(hash->table[i]->key);
	free(hash->table[i]);
	hash->table[i] = ptr;
      }
  free(hash->table);
  free(hash);
}

/* Function:  esl_keyhash_Dump()
 * Incept:    SRE, Sun Jan 30 09:42:22 2005 [St. Louis]
 *
 * Purpose:   (Mainly for debugging purposes.) Dump 
 *            some information about the hash table <h>
 *            to the stream <fp>, which might be stderr
 *            or stdout.
 */
void
esl_keyhash_Dump(FILE *fp, ESL_KEYHASH *h)
{
  struct esl_gki_elem *ptr;
  int i;
  int nkeys;
  int nempty  = 0;
  int maxkeys = -1;
  int minkeys = INT_MAX;

  for (i = 0; i < h->nhash; i++)
    {
      nkeys = 0;
      for (ptr = h->table[i]; ptr != NULL; ptr = ptr->nxt)
	nkeys++;

      if (nkeys == 0)      nempty++;
      if (nkeys > maxkeys) maxkeys = nkeys;
      if (nkeys < minkeys) minkeys = nkeys;
    }

  fprintf(fp, "Total keys:        %d\n", h->nkeys);
  fprintf(fp, "Hash table size:   %d\n", h->nhash);
  fprintf(fp, "Average occupancy: %.1f\n", (float) h->nkeys /(float) h->nhash);
  fprintf(fp, "Unoccupied slots:  %d\n", nempty);
  fprintf(fp, "Most in one slot:  %d\n", maxkeys);
  fprintf(fp, "Least in one slot: %d\n", minkeys);
}



/* Function: esl_gki_Store()
 * Incept:   SRE, Sun Jan 30 09:21:13 2005 [St. Louis]
 *
 * Purpose:  Store a string <key> in the key index hash table <h>.
 *           Associate it with a unique "key index", counting from
 *           0. (It's this index that lets us map the hashed keys to
 *           integer-indexed C arrays, clumsily emulating Perl's
 *           hashes.) Returns this index through <ret_index>.
 *
 *           Does *not* check to see if the key's already in the
 *           table, so it's possible to store multiple copies of a key
 *           with different indices; this is probably not what you
 *           want. If you're not sure the key is unique, check the
 *           table first with <esl_keyhash_GetIndex()>.
 *
 * Returns:  <eslOK> on success, and sets <ret_index>.
 *
 * Throws:   <eslEMEM> on allocation failure, and <ret_index> is -1.
 */
int
esl_gki_Store(ESL_KEYHASH *h, char *key, int *ret_index)
{
  int val;
  struct esl_gki_elem *new;
  int status;

  if (ret_index != NULL) *ret_index = -1;
  val = gki_hashvalue(h, key);
  
  /* Create the new element; don't change the hash until this
   * allocation succeeds.
   * 
   * We could optimize these mallocs by keeping a pool
   * of memory, rather than malloc'ing every individual key.
   */
  ESL_MALLOC(new, sizeof(struct esl_gki_elem));
  if ((status = esl_strdup(key, -1, &(new->key))) != eslOK)
    { free(new); return status; }
  new->idx = h->nkeys;
  
  /* Insert the new element at hash->table[val], at the head
   * of the linked list.
   */
  new->nxt      = h->table[val];
  h->table[val] = new;
  h->nkeys++;

  /* Time to upsize? If we're 3x saturated, expand the hash table.
   */
  if (h->nkeys > 3*h->nhash && h->primelevel < eslGKI_NPRIMES-1)
    gki_upsize(h);

  if (ret_index != NULL) *ret_index = h->nkeys-1; 
  return eslOK;
}

/* Function:  esl_gki_Lookup()
 * Incept:    SRE, Sun Jan 30 09:38:53 2005 [St. Louis]
 *
 * Purpose:   Look up a <key> in the hash table <h> and 
 *            returns its array index (0..nkeys-1), or -1
 *            if <key> isn't found.
 */
int
esl_gki_Lookup(ESL_KEYHASH *h, char *key)
{
  struct esl_gki_elem *ptr;
  int val;
  
  val = gki_hashvalue(h, key);
  for (ptr = h->table[val]; ptr != NULL; ptr = ptr->nxt)
    if (strcmp(key, ptr->key) == 0) return ptr->idx;
  return -1;
}




/* gki_alloc():
 * SRE, Sun Jan 30 09:45:47 2005 [St. Louis]
 * 
 * Allocate a hash table structure with the
 * size given by primelevel.
 *
 * Args:     primelevel - level 0..GKI_NPRIMES-1, specifying
 *                        the size of the table; see gki_primes[]
 *                        array.
 *
 * Returns:  An allocated hash table structure; or NULL on failure.
 */
static ESL_KEYHASH *
gki_alloc(int primelevel)
{
  ESL_KEYHASH *hash;
  int  i;

  if (primelevel < 0 || primelevel >= eslGKI_NPRIMES) 
    ESL_ERROR_NULL(eslEINCONCEIVABLE, "bad primelevel in gki_alloc()");
  
  if ((hash = malloc(sizeof(ESL_KEYHASH))) == NULL)
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  hash->primelevel = primelevel;
  hash->nhash      = gki_primes[hash->primelevel];
  hash->table      = malloc(sizeof(struct esl_gki_elem) * hash->nhash);
  if (hash->table == NULL) 
    { free(hash); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  for (i = 0; i < hash->nhash; i++)
    hash->table[i] = NULL;
  hash->nkeys = 0;
  return hash;
}  


/* gki_hashvalue()
 * SRE, Sun Jan 30 09:50:45 2005 [St. Louis]
 *
 * Calculate the hash value for a key. Usually we expect a one-word
 * key, but the function will hash any ASCII string effectively. The
 * hash function is a simple one (see p. 233 of Sedgewick, Algorithms
 * in C).  Slightly optimized: does two characters at a time before
 * doing the modulo; this gives us a significant speedup.
 *
 * Args:     hash - the gki structure (we need to know the hash table size)
 *           key  - a string to calculate the hash value for;
 *                  this must be 7-bit ASCII, we assume all chars are 0..127.
 *
 * Returns:  a hash value, in the range 0..hash->nhash-1.
 *           hash table is unmodified.
 */
static int
gki_hashvalue(ESL_KEYHASH *hash, char *key)
{
  int val = 0;

  for (; *key != '\0'; key++)
    {
      val = eslGKI_ALPHABETSIZE*val + *key; 
      if (*(++key) == '\0') { val = val % hash->nhash; break; }
      val = (eslGKI_ALPHABETSIZE*val + *key) % hash->nhash;
    }
  return val;
}

/* gki_upsize()
 * SRE, Sun Jan 30 09:50:39 2005 [St. Louis]
 *
 * Grow the hash table to the next available size.
 *
 * Args:     old - the GKI hash table to reallocate.
 *
 * Returns:  <eslOK> on success (the hash table is changed);
 *           <eslEMEM> on failure; the table is already at its maximum size,
 *              or an allocation failed. The hash table is returned unchanged.
 */
static int
gki_upsize(ESL_KEYHASH *old)
{
  ESL_KEYHASH *new;
  int       i;
  struct esl_gki_elem *optr;
  struct esl_gki_elem *nptr;
  int       val;

  if (old->primelevel >= eslGKI_NPRIMES-1) return eslEMEM;
  new = gki_alloc(old->primelevel+1);

  /* Read the old, store in the new, while *not changing*
   * any key indices. Because of the way the lists are
   * treated as LIFO stacks, all the lists are reversed 
   * in the new structure.
   */
  for (i = 0; i < old->nhash; i++)
    {
      optr = old->table[i];
      while (optr != NULL)
	{
	  val = gki_hashvalue(new, optr->key);

	  nptr = new->table[val];
	  new->table[val]      = optr;
	  optr                 = optr->nxt;
	  new->table[val]->nxt = nptr;
	}
    }
  free(old->table);

  /* Now swap within the interior of the structures, so the old
   * structure is updated to the new structure.
   * (nkeys is identical, so we don't need to swap that element.)
   */
  old->primelevel = new->primelevel;
  old->nhash      = new->nhash;
  old->table      = new->table;
  free(new);
  return eslOK;
}

/******************************************************************************
 * Example and test driver
 *****************************************************************************/

#ifdef eslKEYHASH_EXAMPLE
/* gcc -g -Wall -o example -I. -DeslKEYHASH_EXAMPLE keyhash.c easel.c 
 * time ./example /usr/share/dict/words /usr/share/dict/words
 */
#include <stdio.h>

#include <easel/easel.h>
#include <easel/keyhash.h>

int
main(int argc, char **argv)
{
  FILE        *fp;
  char         buf[256];
  char        *s, *tok;
  ESL_KEYHASH *h;
  int          idx;
  int          nstored, nsearched, nshared;

  h = esl_keyhash_Create();

  /* Read/store keys from file 1.
   */
  if ((fp = fopen(argv[1], "r")) == NULL)
    { fprintf(stderr, "couldn't open %s\n", argv[1]); exit(1); }
  nstored = 0;
  while (fgets(buf, 256, fp) != NULL)
    {
      s = buf;
      esl_strtok(&s, " \t\r\n", &tok, NULL);
      esl_gki_Store(h, tok, &idx);
      nstored++;
    }
  fclose(fp);
  printf("Stored %d keys.\n", nstored);

  /* Look up keys from file 2.
   */
  if ((fp = fopen(argv[2], "r")) == NULL)
    { fprintf(stderr, "couldn't open %s\n", argv[2]); exit(1); }
  nsearched = nshared = 0;
  while (fgets(buf, 256, fp) != NULL)
    {
      s = buf;
      esl_strtok(&s, " \t\r\n", &tok, NULL);

      idx = esl_gki_Lookup(h, tok);
      if (idx != -1) nshared++;
      nsearched++;
    }
  fclose(fp);
  printf("Looked up %d keys.\n", nsearched);
  printf("In common: %d keys.\n", nshared);

  esl_keyhash_Destroy(h);
  return 0;
}
#endif /*eslKEYHASH_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
