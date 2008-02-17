/* Partial emulation of Perl hashes (associative arrays),
 * mapping keys (ASCII char strings) to array indices.
 * 
 * Contents:
 *    1. The <ESL_KEYHASH> object.
 *    2. Storing and retrieving keys.
 *    3. Internal functions.        
 *    4. Benchmark driver.
 *    5. Unit tests.
 *    6. Test driver.
 *    7. Example.
 *    8. Copyright and license information.
 * 
 * SRE, Sun Jan 30 09:14:21 2005; from squid's gki.c, 1999.
 * SVN $Id$
 *
 *****************************************************************
 * Limitations:
 *     - hash table can only grow; no provision for deleting keys
 *       or downsizing the hash table.
 *     - Maximum hash table size set at 100003. Performance 
 *       will degrade for key sets much larger than this.
 *     - Assumes that integers are 32 bits (or greater). 
 * 
 *****************************************************************
 * 
 * API for storing/reading keys (strings) and associating
 * them with integer indices in an array:
 * moral equivalent of Perl's $foo{$key} = whatever, $bar{$key} = whatever:
 *
 *       #include <easel.h>
 *       #include <esl_keyhash.h>
 *     
 *       ESL_KEYHASH  *hash;
 *       int   idx;
 *       char *key;
 *       
 *       hash = esl_keyhash_Create();
 * (Storing:) 
 *       (foreach key) {
 *          esl_key_Store(hash, key, &idx);       
 *          (reallocate foo, bar as needed)
 *          foo[idx] = whatever;
 *          bar[idx] = whatever;
 *       }     
 * (Reading:)
 *       (foreach key) {
 *          idx = esl_key_Lookup(hash, key);
 *          if (idx == -1) {no_such_key; }
 *          (do something with) foo[idx];
 *          (do something with) bar[idx];
 *       }   
 *       esl_keyhash_Destroy();
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"
#include "esl_keyhash.h"

/* key_primes[] defines the ascending order of hash table sizes
 * that we use in upsizing the hash table dynamically.
 *
 * Best hash table sizes are prime numbers (see Knuth vol 3, Sorting
 * and Searching). Useful site for testing primes:
 *   http://www.idbsu.edu/people/jbrennan/algebra/numbers/sieve.html
 *
 * Because of the way key_hashvalue works, the largest number
 * must be < INT_MAX / 128 / 128 : 131072 on a 32 bit machine.
 */
static int key_primes[]  = { 101, 1009, 10007, 100003 };
#define eslKEY_NPRIMES      4
#define eslKEY_ALPHABETSIZE 128

static ESL_KEYHASH *key_alloc(int primelevel);
static int          key_hashvalue(ESL_KEYHASH *hash, char *key);
static int          key_upsize(ESL_KEYHASH *old);


/*****************************************************************
 * 1. The <ESL_KEYHASH> object
 *****************************************************************/ 

/* Function:  esl_keyhash_Create()
 * Synopsis:  Allocates a new keyhash.
 * Incept:    SRE, Sun Jan 30 09:17:20 2005 [St. Louis]
 *
 * Purpose:   Create a new hash table for key indexing, and returns
 *            a pointer to it.
 *           
 * Note:      A wrapper around a level 0 key_alloc().
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_KEYHASH *
esl_keyhash_Create(void)
{
  ESL_KEYHASH *kh = key_alloc(0);
  return kh;
}

/* Function:  esl_keyhash_Clone()
 * Synopsis:  Duplicates a keyhash.
 * Incept:    SRE, Fri Feb 15 18:57:50 2008 [Janelia]
 *
 * Purpose:   Allocates and duplicates a keyhash <kh>. Returns a
 *            pointer to the duplicate.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_KEYHASH *
esl_keyhash_Clone(ESL_KEYHASH *kh)
{
  ESL_KEYHASH *nw;		
  int          h;
  struct esl_key_elem *newprv;
  struct esl_key_elem *newelem;
  struct esl_key_elem *oldelem;
  int          status;

  if ((nw = key_alloc(kh->primelevel)) == NULL) return NULL;

  for (h = 0; h < kh->nhash; h++)
    {
      for (newprv = NULL, oldelem = kh->table[h]; oldelem != NULL; oldelem = oldelem->nxt)
	{
	  ESL_ALLOC(newelem, sizeof(struct esl_key_elem)); /* we really ought to malloc more efficiently */
	  if (esl_strdup(oldelem->key, -1, &(newelem->key)) != eslOK) goto ERROR;  /* again, inefficient */
	  newelem->idx = oldelem->idx;
	  newelem->nxt = NULL;

	  if (newprv == NULL) nw->table[h] = newelem;
	  else                newprv->nxt  = newelem;
	  newprv = newelem;
	}
    }
  nw->nkeys = kh->nkeys;
  return nw;
  
 ERROR:
  return NULL;
}


/* Function:  esl_keyhash_Destroy()
 * Synopsis:  Frees a keyhash.
 * Incept:    SRE, Sun Jan 30 09:19:19 2005 [St. Louis]
 *
 * Purpose:   Destroys <kh>.
 *
 * Returns:   (void)
 */
void
esl_keyhash_Destroy(ESL_KEYHASH *kh)
{
  struct esl_key_elem *ptr;
  int i;

  if (kh == NULL) return;	/* tolerate a NULL */

  for (i = 0; i < kh->nhash; i++)
    while (kh->table[i] != NULL)
      {
	ptr = kh->table[i]->nxt;
			/* NULL keys can occur after we've key_upsize'd */
	if (kh->table[i]->key != NULL) free(kh->table[i]->key);
	free(kh->table[i]);
	kh->table[i] = ptr;
      }
  free(kh->table);
  free(kh);
}

/* Function:  esl_keyhash_Dump()
 * Synopsis:  Dumps debugging information about a keyhash.
 * Incept:    SRE, Sun Jan 30 09:42:22 2005 [St. Louis]
 *
 * Purpose:   (Mainly for debugging purposes.) Dump 
 *            some information about the hash table <kh>
 *            to the stream <fp>, which might be stderr
 *            or stdout.
 */
void
esl_keyhash_Dump(FILE *fp, ESL_KEYHASH *kh)
{
  struct esl_key_elem *ptr;
  int i;
  int nkeys;
  int nempty  = 0;
  int maxkeys = -1;
  int minkeys = INT_MAX;

  for (i = 0; i < kh->nhash; i++)
    {
      nkeys = 0;
      for (ptr = kh->table[i]; ptr != NULL; ptr = ptr->nxt)
	nkeys++;

      if (nkeys == 0)      nempty++;
      if (nkeys > maxkeys) maxkeys = nkeys;
      if (nkeys < minkeys) minkeys = nkeys;
    }

  fprintf(fp, "Total keys:        %d\n", kh->nkeys);
  fprintf(fp, "Hash table size:   %d\n", kh->nhash);
  fprintf(fp, "Average occupancy: %.1f\n", (float) kh->nkeys /(float) kh->nhash);
  fprintf(fp, "Unoccupied slots:  %d\n", nempty);
  fprintf(fp, "Most in one slot:  %d\n", maxkeys);
  fprintf(fp, "Least in one slot: %d\n", minkeys);
}
/*--------------- end, <ESL_KEYHASH> object ---------------------*/




/*****************************************************************
 * 2. Storing and retrieving keys
 *****************************************************************/ 

/* Function: esl_key_Store()
 * Synopsis: Store a key, get a key index for it.
 * Incept:   SRE, Sun Jan 30 09:21:13 2005 [St. Louis]
 *
 * Purpose:  Store a string <key> in the key index hash table <kh>.
 *           Associate it with a unique "key index", counting from
 *           0. (It's this index that lets us map the hashed keys to
 *           integer-indexed C arrays, clumsily emulating Perl's
 *           hashes.) Returns this index through <ret_index>.
 *
 * Returns:  <eslOK> on success; stores <key>; <ret_index> is set to
 *           the next higher index value.
 *           Returns <eslEDUP> if <key> was already stored in the table;
 *           <ret_index> is set to the existing index for <key>.
 *
 * Throws:   <eslEMEM> on allocation failure, and sets <ret_index> to -1.
 */
int
esl_key_Store(ESL_KEYHASH *kh, char *key, int *ret_index)
{
  int val;
  struct esl_key_elem *new = NULL;
  int status;

  if (ret_index != NULL) *ret_index = -1;
  val = key_hashvalue(kh, key);
  
  /* Was this key already stored?
   */
  for (new = kh->table[val]; new != NULL; new = new->nxt)
    if (strcmp(key, new->key) == 0) { *ret_index = new->idx; return eslEDUP;}

  /* If not, create the new element; don't change the hash until this
   * allocation succeeds.
   * 
   * We could optimize these allocations by keeping a pool
   * of memory, rather than allocating every individual key.
   */
  ESL_ALLOC(new, sizeof(struct esl_key_elem));
  new->key = NULL;
  if ((status = esl_strdup(key, -1, &(new->key))) != eslOK) goto ERROR;
  new->idx = kh->nkeys;
  
  /* Insert the new element at kh->table[val], at the head
   * of the linked list.
   */
  new->nxt       = kh->table[val];
  kh->table[val] = new;
  kh->nkeys++;

  /* Time to upsize? If we're 3x saturated, expand the hash table.
   */
  if (kh->nkeys > 3*kh->nhash && kh->primelevel < eslKEY_NPRIMES-1)
    if ((status = key_upsize(kh)) != eslOK) goto ERROR;

  if (ret_index != NULL) *ret_index = kh->nkeys-1; 
  return eslOK;

 ERROR:
  if (new != NULL) {
    if (new->key != NULL) free(new->key);
    free(new);
  }
  if (ret_index != NULL) *ret_index = -1;
  return status;
}

/* Function:  esl_key_Lookup()
 * Incept:    SRE, Sun Jan 30 09:38:53 2005 [St. Louis]
 *
 * Purpose:   Look up a <key> in the hash table <h> and 
 *            returns its array index (0..nkeys-1), or -1
 *            if <key> isn't found.
 */
int
esl_key_Lookup(ESL_KEYHASH *h, char *key)
{
  struct esl_key_elem *ptr;
  int val;
  
  val = key_hashvalue(h, key);
  for (ptr = h->table[val]; ptr != NULL; ptr = ptr->nxt)
    if (strcmp(key, ptr->key) == 0) return ptr->idx;
  return -1;
}
/*---------- end, API for storing/retrieving keys ---------------*/




/*****************************************************************
 * 3. Internal functions
 *****************************************************************/ 

/* key_alloc():
 * SRE, Sun Jan 30 09:45:47 2005 [St. Louis]
 * 
 * Allocate a hash table structure with the
 * size given by primelevel.
 *
 * Args:     primelevel - level 0..KEY_NPRIMES-1, specifying
 *                        the size of the table; see key_primes[]
 *                        array.
 *
 * Returns:  An allocated hash table structure; or NULL on failure.
 */
static ESL_KEYHASH *
key_alloc(int primelevel)
{
  int status;
  ESL_KEYHASH *hash = NULL;
  int  i;

  if (primelevel < 0 || primelevel >= eslKEY_NPRIMES) 
    ESL_XEXCEPTION(eslEINCONCEIVABLE, "bad primelevel in key_alloc()");
  
  ESL_ALLOC(hash, sizeof(ESL_KEYHASH));
  hash->table = NULL;

  hash->primelevel = primelevel;
  hash->nhash      = key_primes[hash->primelevel];
  ESL_ALLOC(hash->table, sizeof(struct esl_key_elem *) * hash->nhash);
  for (i = 0; i < hash->nhash; i++)
    hash->table[i] = NULL;
  hash->nkeys = 0;
  return hash;

 ERROR:
  if (hash != NULL) {
    if (hash->table != NULL) free(hash->table);
    free(hash);
  }
  return NULL;
}  


/* key_hashvalue()
 * SRE, Sun Jan 30 09:50:45 2005 [St. Louis]
 *
 * Calculate the hash value for a key. Usually we expect a one-word
 * key, but the function will hash any ASCII string effectively. The
 * hash function is a simple one (see p. 233 of Sedgewick, Algorithms
 * in C).  Slightly optimized: does two characters at a time before
 * doing the modulo; this gives us a significant speedup.
 *
 * Since we expect primarily alphabetic strings, we could probably
 * find a better hashfunction than this. 
 * 
 * Args:     hash - the key structure (we need to know the hash table size)
 *           key  - a string to calculate the hash value for;
 *                  this must be 7-bit ASCII, we assume all chars are 0..127.
 *
 * Returns:  a hash value, in the range 0..hash->nhash-1.
 *           hash table is unmodified.
 */
static int
key_hashvalue(ESL_KEYHASH *hash, char *key)
{
  int val = 0;

  for (; *key != '\0'; key++)
    {
      val = eslKEY_ALPHABETSIZE*val + *key; 
      if (*(++key) == '\0') { val = val % hash->nhash; break; }
      val = (eslKEY_ALPHABETSIZE*val + *key) % hash->nhash;
    }
  return val;
}

/* key_upsize()
 * SRE, Sun Jan 30 09:50:39 2005 [St. Louis]
 *
 * Grow the hash table to the next available size.
 *
 * Args:     old - the KEY hash table to reallocate.
 *
 * Returns:  <eslOK> on success. 'Success' includes the case
 *           where the hash table is already at its maximum size,
 *           and cannot be upsized any more.
 *           
 * Throws:   <eslEMEM> on allocation failure, and
 *           the hash table is left in its initial state.
 */
static int
key_upsize(ESL_KEYHASH *old)
{
  ESL_KEYHASH *new;
  int       i;
  struct esl_key_elem *optr;
  struct esl_key_elem *nptr;
  int       val;

  if (old->primelevel >= eslKEY_NPRIMES-1) return eslOK; /* quasi-success */

  new = key_alloc(old->primelevel+1);
  if (new == NULL) return eslEMEM; /* percolation */

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
	  val = key_hashvalue(new, optr->key);

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
/*--------------- end, internal functions -----------------*/


/*****************************************************************
 * 4. Benchmark driver
 *****************************************************************/
#ifdef eslKEYHASH_BENCHMARK
/* gcc -g -O2 -o keyhash_benchmark -I. -L. -DeslKEYHASH_BENCHMARK esl_keyhash.c -leasel -lm
 * time ./keyhash_benchmark /usr/share/dict/words /usr/share/dict/words
 */
#include <stdio.h>
#include "easel.h"
#include "esl_keyhash.h"

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
      esl_key_Store(h, tok, &idx);
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

      idx = esl_key_Lookup(h, tok);
      if (idx != -1) nshared++;
      nsearched++;
    }
  fclose(fp);
  printf("Looked up %d keys.\n", nsearched);
  printf("In common: %d keys.\n", nshared);

  esl_keyhash_Destroy(h);
  return 0;
}
#endif /*eslKEYHASH_BENCHMARK*/
/*------------------- end, benchmark driver ---------------------*/


/*****************************************************************
 * 5. Unit tests
 *****************************************************************/


/*---------------------- end, unit tests ------------------------*/

/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef eslKEYHASH_TESTDRIVE
/* gcc -g -Wall -o test -I. -DeslKEYHASH_TESTDRIVE keyhash.c easel.c 
 * ./test
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <easel.h>
#include <esl_keyhash.h>

#define NSTORE  1200
#define NLOOKUP 1200
#define KEYLEN  2

int
main(int argc, char **argv)
{
  ESL_KEYHASH *h;
  char keys[NSTORE+NLOOKUP][KEYLEN+1]; 
  int  i,j,nk,k42;
  int  nmissed;
  int  status;

  /* Generate 2400 random k=2 keys. 26^2 = 676 possible.
   * We'll store the first 1200 and search on the remaining
   * 1200. We're ~1.775x saturated; expect Poisson P(0) = 17% miss
   * rate on the searches, so we ought to exercise hits and
   * misses on the lookups.
   */
  srand(31);
  for (i = 0; i < NSTORE+NLOOKUP; i++)
    {
      for (j = 0; j < KEYLEN; j++)
	keys[i][j] = 'a' + (rand() % 26); /* yeah, low-order bits; so sue me */
      keys[i][j] = '\0';
    }
  /* spike a known one in (XX.. at key 42).
   */
  for (j = 0; j < KEYLEN; j++)
    keys[42][j] = 'X';

  h = esl_keyhash_Create();
  nk = 0;
  for (i = 0; i < NSTORE; i++)
    {
      status = esl_key_Store(h, keys[i], &j);
      if      (status == eslOK)   { assert(j==nk); nk++; }
      else if (status == eslEDUP) { assert(j<nk); }
      else esl_fatal("store failed.");

      if (i == 42) { k42 = j;}	/* remember where key 42 went */
    }
  nmissed = 0;
  for (i = NSTORE; i < NSTORE+NLOOKUP; i++)
    {
      j = esl_key_Lookup(h, keys[i]);
      if (j == -1) nmissed++;
    }
  j = esl_key_Lookup(h, keys[42]);
  assert(j==k42);

  /* 
  printf("missed %d/%d (%.1f%%)\n", nmissed, NLOOKUP, 
	 100. * (float) nmissed / (float) NLOOKUP);
  esl_keyhash_Dump(stdout, h);
  */

  esl_keyhash_Destroy(h);
  exit (0);
}
#endif /*eslKEYHASH_TESTDRIVE*/

/*--------------------- end, test driver ------------------------*/



/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef eslKEYHASH_EXAMPLE
/* gcc -g -Wall -o example -I. -DeslKEYHASH_EXAMPLE keyhash.c easel.c 
 * time ./example /usr/share/dict/words /usr/share/dict/words
 */
#include <stdio.h>
#include <easel.h>
#include <esl_keyhash.h>

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
      esl_key_Store(h, tok, &idx);
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

      idx = esl_key_Lookup(h, tok);
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
/*----------------------- end, example --------------------------*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
