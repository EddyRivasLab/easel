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
 *       #include "easel.h"
 *       #include "esl_keyhash.h"
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
 *          if (esl_key_Lookup(hash, key, &idx) != eslOK) { no_such_key; }
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

static ESL_KEYHASH *keyhash_create(int primelevel, int init_table_alloc, int init_key_alloc);
static int          key_hashvalue(const ESL_KEYHASH *hash, const char *key);
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
  return keyhash_create(0,     /* primelevel 0 (hashtable has 101 entries)          */
			256,   /* initial alloc for 256 key_elem structures         */
			2048); /* initial alloc for keys totalling up to 2048 chars */
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
esl_keyhash_Clone(const ESL_KEYHASH *kh)
{
  ESL_KEYHASH *nw;		
  int          h;
  struct esl_key_elem *newprv;
  struct esl_key_elem *newelem;
  struct esl_key_elem *oldelem;
  int          status;

  if ((nw = keyhash_create(kh->primelevel, kh->talloc, kh->kalloc)) == NULL) goto ERROR;

  nw->nkeys = 0;
  for (h = 0; h < kh->nhash; h++)
    {
      for (newprv = NULL, oldelem = kh->table[h]; oldelem != NULL; oldelem = oldelem->nxt)
	{
	  newelem      = &(kh->table_mem[nw->nkeys]);   /* get a new structure from the pool */
	  newelem->key = nw->key_mem + nw->kn;          /* get memory pointer for key        */
	  strcpy(newelem->key, oldelem->key);           /* copy old key to new               */
	  newelem->idx = oldelem->idx;
	  newelem->nxt = NULL;
	  nw->nkeys++;
	  nw->kn      += strlen(newelem->key)+1;	/* bookkeeping in key string pool    */

	  if (newprv == NULL) nw->table[h] = newelem;   /* attach new elem to hash table     */
	  else                newprv->nxt  = newelem;
	  newprv = newelem;
	}
    }
  if (nw->nkeys != kh->nkeys) ESL_XEXCEPTION(eslEINCONCEIVABLE, "oops, that can't happen");
  return nw;
  
 ERROR:
  esl_keyhash_Destroy(nw);
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
  if (kh == NULL) return;	
  if (kh->table     != NULL) free(kh->table);
  if (kh->table_mem != NULL) free(kh->table_mem);
  if (kh->key_mem   != NULL) free(kh->key_mem);
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
esl_keyhash_Dump(FILE *fp, const ESL_KEYHASH *kh)
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

  fprintf(fp, "Total keys:             %d\n", kh->nkeys);
  fprintf(fp, "Hash table size:        %d\n", kh->nhash);
  fprintf(fp, "Average occupancy:      %.1f\n", (float) kh->nkeys /(float) kh->nhash);
  fprintf(fp, "Unoccupied slots:       %d\n", nempty);
  fprintf(fp, "Most in one slot:       %d\n", maxkeys);
  fprintf(fp, "Least in one slot:      %d\n", minkeys);
  fprintf(fp, "Keys allocated for:     %d\n", kh->talloc);
  fprintf(fp, "Key string space alloc: %d\n", kh->kalloc);
  fprintf(fp, "Key string space used:  %d\n", kh->kn);
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
esl_key_Store(ESL_KEYHASH *kh, const char *key, int *ret_index)
{
  struct esl_key_elem *new = NULL;
  int val = key_hashvalue(kh, key); /* hash the key */
  int idx = -1;			    /* this'll be the key array index we return */
  int n;
  void *p;  
  int   h;
  int status;

  /* Was this key already stored?  */
  for (new = kh->table[val]; new != NULL; new = new->nxt)
    if (strcmp(key, new->key) == 0) { *ret_index = new->idx; return eslEDUP; }

  /* Get a pointer to the new element, reallocating element memory if needed */
  if (kh->nkeys == kh->talloc) { 
    struct esl_key_elem *oldp = kh->table_mem;
    ESL_RALLOC(kh->table_mem, p, sizeof(struct esl_key_elem)*kh->talloc*2);
    kh->talloc *=2;
    /* By reallocating, you've just invalidated all pointers to key_elems. Fix them. */
    for (h = 0; h < kh->nkeys; h++) 
      if (kh->table_mem[h].nxt != NULL)
	kh->table_mem[h].nxt += kh->table_mem - oldp;
    for (h = 0; h < kh->nhash; h++)
      if (kh->table[h]         != NULL)
	kh->table[h]         += kh->table_mem - oldp;
  }
  new = &(kh->table_mem[kh->nkeys]);

  /* Get a pointer to space for the new key, reallocating key memory if needed */
  n = strlen(key);
  if (kh->kn + n + 1 > kh->kalloc) {
    do {
      char *oldp = kh->key_mem;
      ESL_RALLOC(kh->key_mem, p, sizeof(char) * kh->kalloc * 2);
      kh->kalloc *= 2;
      /* By reallocating, you've just invalidated all pointers to keys. Fix them. */
      for (h = 0; h < kh->nkeys; h++) kh->table_mem[h].key += kh->key_mem - oldp;
    } while (kh->kn + n + 1 > kh->kalloc);
  }
  new->key = kh->key_mem + kh->kn;

  /* Copy the key, assign an index */
  idx      = kh->nkeys;
  strcpy(new->key, key);
  kh->kn  += n+1;
  new->idx = idx;

  /* Insert the new element at head of kh->table[val] */
  new->nxt       = kh->table[val];
  kh->table[val] = new;
  kh->nkeys++;

  /* Time to upsize? If we're 3x saturated, expand the hash table */
  if (kh->nkeys > 3*kh->nhash && kh->primelevel < eslKEY_NPRIMES-1)
    if ((status = key_upsize(kh)) != eslOK) goto ERROR;

  if (ret_index != NULL) *ret_index = idx;
  return eslOK;

 ERROR:
  if (ret_index != NULL) *ret_index = -1;
  return status;
}

/* Function:  esl_key_Lookup()
 * Synopsis:  Lookup a key's array index.
 * Incept:    SRE, Sun Jan 30 09:38:53 2005 [St. Louis]
 *
 * Purpose:   Look up a <key> in the hash table <kh>.
 *            If <key> is found, return <eslOK>, and optionally set <*opt_index>
 *            to its array index (0..nkeys-1).
 *            If <key> is not found, return <eslENOTFOUND>, and
 *            optionally set <*opt_index> to -1.
 */
int
esl_key_Lookup(const ESL_KEYHASH *kh, const char *key, int *opt_index)
{
  struct esl_key_elem *ptr;
  int val  = key_hashvalue(kh, key);
  
  for (ptr = kh->table[val]; ptr != NULL; ptr = ptr->nxt)
    if (strcmp(key, ptr->key) == 0) 
      { 
	if (opt_index != NULL) *opt_index = ptr->idx;
	return eslOK; 
      }

  if (opt_index != NULL) *opt_index = -1;
  return eslENOTFOUND;
}
/*---------- end, API for storing/retrieving keys ---------------*/




/*****************************************************************
 * 3. Internal functions
 *****************************************************************/ 

/* keyhash_create()
 * SRE, Sun Jan 30 09:45:47 2005 [St. Louis]
 * 
 * The real creation function, which takes arguments for memory sizes.
 * This is abstracted to a static function because it's used by both
 * Create() and Clone() but slightly differently.
 *
 * Args:     primelevel - level 0..KEY_NPRIMES-1, specifying
 *                        the size of the table; see key_primes[]
 *                        array.
 *           init_table_alloc - initial allocation for key_elem structures.
 *           init_key_alloc   - initial allocation for key strings.
 *
 * Returns:  An allocated hash table structure; or NULL on failure.
 */
ESL_KEYHASH *
keyhash_create(int primelevel, int init_table_alloc, int init_key_alloc)
{
  ESL_KEYHASH *kh = NULL;
  int  i;
  int  status;

  if (primelevel < 0 || primelevel >= eslKEY_NPRIMES) 
    ESL_XEXCEPTION(eslEINCONCEIVABLE, "bad primelevel");
  
  ESL_ALLOC(kh, sizeof(ESL_KEYHASH));
  kh->table     = NULL;
  kh->table_mem = NULL;
  kh->key_mem   = NULL;

  kh->primelevel = primelevel;
  kh->nhash      = key_primes[kh->primelevel];

  ESL_ALLOC(kh->table, sizeof(struct esl_key_elem *) * kh->nhash);
  for (i = 0; i < kh->nhash; i++)
    kh->table[i] = NULL;
  kh->nkeys = 0;

  ESL_ALLOC(kh->table_mem, sizeof(struct esl_key_elem) * init_table_alloc);
  ESL_ALLOC(kh->key_mem,   sizeof(char)                * init_key_alloc);
  kh->talloc = init_table_alloc;
  kh->kalloc = init_key_alloc;
  kh->kn     = 0;
  return kh;

 ERROR:
  esl_keyhash_Destroy(kh);
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
key_hashvalue(const ESL_KEYHASH *kh, const char *key)
{
  int val = 0;

  for (; *key != '\0'; key++)
    {
      val = eslKEY_ALPHABETSIZE*val + *key; 
      if (*(++key) == '\0') { val = val % kh->nhash; break; }
      val = (eslKEY_ALPHABETSIZE*val + *key) % kh->nhash;
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
  struct esl_key_elem **newtbl = NULL;
  struct esl_key_elem *optr;
  struct esl_key_elem *nptr;
  int       old_nhash;
  int       h;
  int       val;
  int       status;

  if (old->primelevel >= eslKEY_NPRIMES-1) return eslOK; /* quasi-success (can't grow any more) */

  /* The catch here is that when you upsize the table, all the hash functions
   * change; so you have to go through all the keys, recompute their hash functions,
   * and store them again in the new table.
   */

  /* Allocate a new, larger hash table. (Don't change <kh> until this succeeds) */
  old_nhash     = old->nhash;
  ESL_ALLOC(newtbl, sizeof(struct esl_key_elem *) * (key_primes[old->primelevel+1]));
  old->primelevel++;
  old->nhash    = key_primes[old->primelevel];

  for (h = 0; h < old->nhash; h++)
    newtbl[h] = NULL;

  /* Traverse the old, store in the new, while *not changing* any key
   * indices. Because of the way the lists are treated as LIFO stacks,
   * all the lists are reversed in the new structure.
   */
  for (h = 0; h < old->nkeys; h++)
    {
      optr = &(old->table_mem[h]);
      val  = key_hashvalue(old, optr->key);

      nptr             = newtbl[val]; /* pick up the right head node in the new table  */
      newtbl[val]      = optr;        /* attach the key structure as the new head node */
      newtbl[val]->nxt = nptr;        /* and reattach the rest of the list to it       */
    }

  /* Now swap the tables, free'ing the old. */
  free(old->table);
  old->table      = newtbl;
  return eslOK;

 ERROR:
  if (newtbl != NULL) free(newtbl);
  return eslEMEM;
}
/*--------------- end, internal functions -----------------*/


/*****************************************************************
 * 4. Benchmark driver
 *****************************************************************/
#ifdef eslKEYHASH_BENCHMARK
/* 
   gcc -g -O2 -o keyhash_benchmark -I. -L. -DeslKEYHASH_BENCHMARK esl_keyhash.c -leasel -lm
   time ./keyhash_benchmark /usr/share/dict/words /usr/share/dict/words
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
  ESL_KEYHASH *kh;
  int          idx;
  int          nstored, nsearched, nshared;

  kh = esl_keyhash_Create();

  /* Read/store keys from file 1.
   */
  if ((fp = fopen(argv[1], "r")) == NULL)
    { fprintf(stderr, "couldn't open %s\n", argv[1]); exit(1); }
  nstored = 0;
  while (fgets(buf, 256, fp) != NULL)
    {
      s = buf;
      esl_strtok(&s, " \t\r\n", &tok, NULL);
      esl_key_Store(kh, tok, &idx);
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

      if (esl_key_Lookup(kh, tok, &idx) == eslOK) nshared++;
      nsearched++;
    }
  fclose(fp);
  printf("Looked up %d keys.\n", nsearched);
  printf("In common: %d keys.\n", nshared);

  esl_keyhash_Destroy(kh);
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
#include "easel.h"
#include "esl_keyhash.h"

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
      if (esl_key_Lookup(h, keys[i], &j) != eslOK) nmissed++;
    }
  esl_key_Lookup(h, keys[42], &j);
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

      if (esl_key_Lookup(h, tok, &idx) == eslOK) nshared++;
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
