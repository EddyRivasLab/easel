/* Huffman codes for digitized alphabets.
 * 
 * Contents:
 *   1. Internal functions and structures; components of creating huffman codes
 *   2. The ESL_HUFFMAN object
 *   3. Huffman encoding
 *   4. Huffman decoding
 *   5. Debugging, development
 *   6. Example driver
 *   7. Copyright and license information.
 */

#include "easel.h"
#include "esl_quicksort.h"
#include "esl_huffman.h"


/*****************************************************************
 * 1. Internal functions and structures
 *****************************************************************/

struct hufftree_s {
  float val;    // Sum of frequencies of all leaves under this node
  int   depth;  // Depth of node
  int   left;   // index of left child in array of tree nodes (0..N-2; 0 is the root)
  int   right;  //  "" for right child
};

/* sort_floats_decreasing()
 * Sorting function for esl_quicksort(), putting 
 * symbol frequencies in decreasing order.
 */
static int
sort_floats_decreasing(void *data, int e1, int e2)
{
  float *fq = (float *) data;
  if (fq[e1] > fq[e2]) return -1;
  if (fq[e1] < fq[e2]) return 1;
  return 0;
}

/* sort_canonical()
 * Sorting function for esl_quicksort(), putting symbols into
 * canonical Huffman order: primarily by ascending code length,
 * secondarily by ascending symbol code.
 */
static int
sort_canonical(void *data, int e1, int e2)
{
  ESL_HUFFMAN *hm = (ESL_HUFFMAN *) data;
  int          L1 = hm->len[e1];
  int          L2 = hm->len[e2];
  
  if      (L2 == 0) return -1;   // len=0 means symbol isn't encoded at all, doesn't occur
  else if (L1 == 0) return 1;
  else if (L1 < L2) return -1;
  else if (L1 > L2) return 1;
  else if (e1 < e2) return -1;
  else if (e1 > e2) return 1;
  else              return 0;
}

/* Build the Huffman tree, joining nodes/leaves of smallest frequency.
 * This takes advantage of having the fq[] array sorted, and the fact
 * that the internal node values also come out sorted... i.e. we don't
 * have to re-sort, we can always find the smallest leaves/nodes by 
 * looking at the last ones.
 * 
 * Input: 
 *   hm->sorted_at[] lists symbol indices from largest to smallest freq.
 *   hm->Ku          is the number of syms w/ nonzero freq; tree has Ku-1 nodes
 *   htree           blank, allocated for at least Ku-1 nodes
 *   
 * Output:
 *   htree's left, right, val fields are filled.
 *   
 * Returns:
 *   <eslOK> on success.  
 */
static int
huffman_tree(ESL_HUFFMAN *hm, struct hufftree_s *htree, float *fq)
{
  int r = hm->Ku-1;   // r = smallest leaf symbol that hasn't been included in tree yet; r+1 = # of leaves left
  int k = hm->Ku-2;   // k = smallest internal node not used as a child yet; k-j = # nodes not used as child yet
  int j;

  for (j = hm->Ku-2; j >= 0; j--)  // j = index of next node we add; we add one per iteration
    {       
      /* Should we join two leaves?
       *   If we have no internal nodes yet (because we're just starting),
       *   or the two smallest frequencies are <= the smallest unjoined node's value
       */
      if ( (j == hm->Ku-2) ||  (r >= 1 && fq[hm->sorted_at[r]] <= htree[k].val))
	{
	  htree[j].right = -hm->sorted_at[r];    // leaves are signified by negative indices in tree
	  htree[j].left  = -hm->sorted_at[r-1];
	  htree[j].val   = fq[hm->sorted_at[r]] + fq[hm->sorted_at[r-1]];
	  r -= 2;
	}

      /* Or should we join two nodes?
       *  If we have no leaves left, 
       *  or (we do have two nodes) and both are smaller than smallest unjoined leaf's value
       */
       else if (r == -1  || (k-j >= 2 && htree[k-1].val < fq[hm->sorted_at[r]]))
	 {
	   htree[j].right = k;
	   htree[j].left  = k-1;
	   htree[j].val   = htree[k].val + htree[k-1].val;
	   k -= 2;
	 }
      
      /* Otherwise, we join smallest node and smallest leaf. */
       else 
	 {
	   htree[j].right = -hm->sorted_at[r];
	   htree[j].left  = k;
	   htree[j].val   = fq[hm->sorted_at[r]] + htree[k].val;
	   r--;
	   k--;
	 }
    }
  return eslOK;
}


/* Calculate code lengths, equal to the depth of each node. 
 * Traverse the tree, calculating depth of each node, starting with
 * depth 0 for root 0. We don't need a stack for this traversal,
 * tree is already indexed in traversal order.
 * 
 * Input:
 *   hm->Ku          is the number of syms w/ nonzero freqs; tree has Ku-1 nodes.
 *   htree[0..Ku-2]  is the constructed Huffman tree, with right/left/val set.
 *   htree[].len     has been initialized to 0 for all symbols 0..K
 *
 * Output:
 *   htree's depth field is set.
 *   hm->len is set for all encoded symbols (left at 0 for unused symbols)
 *   hm->Lmax is set
 *   
 * Return: 
 *   <eslOK> on success
 *   <eslERANGE> if max code length > eslHUFFMAN_MAXCODE and won't fit in uint32_t   
 */
static int
huffman_codelengths(ESL_HUFFMAN *hm, struct hufftree_s *htree, float *fq)
{
  int i;

  htree[0].depth = 0;
  for (i = 0; i < hm->Ku-1; i++)
    {
      if (htree[i].right <= 0) hm->len[-htree[i].right]    = htree[i].depth + 1;
      else                     htree[htree[i].right].depth = htree[i].depth + 1;
      
      if (htree[i].left <= 0)  hm->len[-htree[i].left]     = htree[i].depth + 1;
      else                     htree[htree[i].left].depth  = htree[i].depth + 1;
    }

  hm->Lmax = 0;
  for (i = 0; i < hm->K; i++)
    hm->Lmax = ESL_MAX(hm->len[i], hm->Lmax);

  return (hm->Lmax > eslHUFFMAN_MAXCODE ? eslERANGE : eslOK);
}


/* huffman_canonize()
 * Given code lengths, now we calculate the canonical Huffman encoding.
 * 
 * Input:
 *   hm->len[]  code lengths are set for all K (0 for unused symbols)
 *   hm->code[] have been initialized to 0 for all K
 *   
 * Output:  
 *   hm->code[] have been set for all used symbols.
 *   hm->D      number of different code lengths is set
 *   
 * Returns:
 *  <eslOK> on success.  
 */
static int
huffman_canonize(ESL_HUFFMAN *hm)
{
  int i,r;

  /* Sort symbols according to 1) code length; 2) order in digital alphabet (i.e. symbol code itself)
   * Reuse/reset <sorted_at>.
   * You can't just sort the encoded Ku; you have to sort all K, because
   * quicksort expects indices to be contiguous (0..K-1).
   */
  esl_quicksort(hm, hm->K, sort_canonical, hm->sorted_at);

  /* Assign codes. (All K have been initialized to zero already.) */
  for (r = 1; r < hm->Ku; r++)
    hm->code[hm->sorted_at[r]] =
      (hm->code[hm->sorted_at[r-1]] + 1) << (hm->len[hm->sorted_at[r]] - hm->len[hm->sorted_at[r-1]]);


  /* Set D, the number of different code lengths */
  hm->D = 1;
  for (r = 1; r < hm->Ku; r++)
    if (hm->len[hm->sorted_at[r]] > hm->len[hm->sorted_at[r-1]]) hm->D++;

  return eslOK;
}


/* huffman_decoding_table()
 * Given a canonical Huffman code; build the table that lets us
 * efficiently decode it.
 * 
 * Input:
 *   hm->K         is set: total # of symbols (inclusive of unused ones)
 *   hm->Ku        is set: total # of encoded/used symbols
 *   hm->code      is set: canonical Huffman codes for symbols 0..K-1
 *   hm->len       is set: code lengths for symbols 0..K-1
 *   hm->sorted_at is set: canonical Huffman sort order 
 *   hm->Lmax      is set: maximum code length
 *   hm->D         is set: # of different code lengths
 *
 *   hm->dt_len    is allocated for hm->D, but otherwise uninitialized
 *   hm->dt_lcode  is allocated for hm->D, but otherwise uninitialized
 *   hm->dt_rank   is allocated for hm->D, but otherwise uninitialized
 *   
 * Output:  
 *   hm->dt_len    is set: lengths of each used code length 0..D-1
 *   hm->dt_lcode  is set: left-flushed first code for each code length [d]
 *   hm->dt_rank   is set: rank r for 1st code for each used code length [d]
 */
static int 
huffman_decoding_table(ESL_HUFFMAN *hm)
{
  int r;
  int D = 0;

  hm->dt_len[0]   = hm->len[hm->sorted_at[0]];
  hm->dt_lcode[0] = hm->code[hm->sorted_at[0]] << (eslHUFFMAN_MAXCODE - hm->len[hm->sorted_at[0]]);
  hm->dt_rank[0]  = 0;
  for (r = 1; r < hm->Ku; r++)
    if (hm->len[hm->sorted_at[r]] > hm->len[hm->sorted_at[r-1]]) 
      {
	D++;
	hm->dt_len[D]   = hm->len[hm->sorted_at[r]];
	hm->dt_lcode[D] = hm->code[hm->sorted_at[r]] << (eslHUFFMAN_MAXCODE - hm->len[hm->sorted_at[r]]);
	hm->dt_rank[D]  = r;
      }
  ESL_DASSERT1(( hm->D == D ));
  return eslOK;
}


void
dump_uint32(FILE *fp, uint32_t v, int L)
{
  uint32_t mask;
  int      i;
  
  for (mask = 1 << (L-1), i = L; i >= 1; i--, mask = mask >> 1)
    putc( ((v & mask) ? '1' : '0'), fp);
}





/*****************************************************************
 * 2. The ESL_HUFFMAN object
 *****************************************************************/

/* Function:  esl_huffman_Create()
 * Synopsis:  Create a new Huffman code.
 * Incept:    SRE, Thu Nov 12 11:08:09 2015
 *
 * Purpose:   Create a canonical Huffman code for observed symbol 
 *            frequencies <fq[0..K]> for <K> possible symbols.
 *            
 *            If you're encoding an Easel digital alphabet, 
 *            <K = abc->Kp>, inclusive of ambiguity codes, gaps,
 *            missing data, rare digital codes.
 *
 *            If you're encoding 7-bit ASCII text, K=128 and the
 *            symbols codes are ASCII codes.
 *            
 *            If you're encoding MTF-encoded ASCII text, K=128
 *            and the "symbol" codes are 0..127 offsets in 
 *            the move-to-front encoding.
 *
 *            Unobserved symbols (with <fq[] = 0>) will not be encoded;
 *            they get a code length of 0, and a code of 0.
 */
ESL_HUFFMAN *
esl_huffman_Create(float *fq, int K)
{
  ESL_HUFFMAN       *hm    = NULL;
  struct hufftree_s *htree = NULL;  // only need tree temporarily, during code construction.
  int                i,j,k,r;
  int                status;

  ESL_ALLOC(hm, sizeof(ESL_HUFFMAN));
  hm->len       = NULL;
  hm->code      = NULL;
  hm->sorted_at = NULL;
  hm->dt_len    = NULL;
  hm->dt_lcode  = NULL;
  hm->dt_rank   = NULL;

  hm->K         = K;
  hm->Ku        = 0;
  hm->D         = 0;
  hm->Lmax      = 0;

  ESL_ALLOC(hm->len,       sizeof(int)               * hm->K);
  ESL_ALLOC(hm->code,      sizeof(uint32_t)          * hm->K);
  ESL_ALLOC(hm->sorted_at, sizeof(int)               * hm->K);

  for (i = 0; i < hm->K; i++) hm->len[i]  = 0;
  for (i = 0; i < hm->K; i++) hm->code[i] = 0;
  
  /* Sort the symbol frequencies, largest to smallest */
  esl_quicksort(fq, hm->K, sort_floats_decreasing, hm->sorted_at);
  
  /* Figure out how many are nonzero */
  for (r = hm->K-1; r >= 0; r--)
    if (fq[hm->sorted_at[r]] > 0.) break;
  hm->Ku = r+1;

  ESL_ALLOC(htree,         sizeof(struct hufftree_s) * (hm->Ku-1));
  if ( (status = huffman_tree       (hm, htree, fq)) != eslOK) goto ERROR;
  if ( (status = huffman_codelengths(hm, htree, fq)) != eslOK) goto ERROR; // can fail eslERANGE on maxlen > 32
  if ( (status = huffman_canonize   (hm))            != eslOK) goto ERROR;


  ESL_ALLOC(hm->dt_len,   sizeof(int)      * hm->D);
  ESL_ALLOC(hm->dt_lcode, sizeof(uint32_t) * hm->D);
  ESL_ALLOC(hm->dt_rank,  sizeof(int)      * hm->D);
  if ( (status = huffman_decoding_table(hm))         != eslOK) goto ERROR;

  free(htree);
  return hm;

 ERROR:
  if (hm)    esl_huffman_Destroy(hm);
  if (htree) free(htree);
  return NULL;
}



/* Function:  esl_huffman_Destroy()
 * Synopsis:  Free an <ESL_HUFFMAN> code.
 * Incept:    SRE, Thu Nov 12 11:07:39 2015
 */
void
esl_huffman_Destroy(ESL_HUFFMAN *hm)
{
  if (hm) {
    if (hm->len)       free(hm->len);
    if (hm->code)      free(hm->code);
    if (hm->sorted_at) free(hm->sorted_at);
    if (hm->dt_len)    free(hm->dt_len);
    if (hm->dt_lcode)  free(hm->dt_lcode);
    if (hm->dt_rank)   free(hm->dt_rank);
    free(hm);
  }
}

  

/*****************************************************************
 * 3. Encoding
 *****************************************************************/

/* huffman_pack()
 * <X[i]> is the current uint32_t unit in the encoded buffer <X>. It
 * has <a> bits in it, maximally left-shifted. 32-a bits are
 * available.
 * 
 * <code> is the next Huffman code to pack into the buffer, of length <L>.
 * 
 * If L < 32-a, then we just pack it into X[i]. Else, we pack
 * what we can into X[i], and leave the remainder in X[i+1].
 * 
 * We update <i> and <a> for <X> accordingly... so we pass them by
 * reference in <ip> and <ap>.
 */
static void
huffman_pack(uint32_t *X, int *ip, int *ap, uint32_t code, int L)
{
  int w = 32 - (*ap+L);
  
  if (w > 0)      // code can pack into X[i]'s available space.
    {
      X[*ip] = X[*ip] | (code << w);
      *ap += L;
    }
  else if (w < 0) // code packs partly in X[i], remainder in X[i+1].
    {
      X[*ip] = X[*ip] | (code >> (-w));
      (*ip)++;
      X[*ip] = code >> (32+w);
      (*ap) = -w;
    }
  else           // code packs exactly.
    {
      X[*ip] = X[*ip] | (code << w);
      *ap =  0;
      *ip += 1;
    }
}

int
esl_huffman_Encode(ESL_HUFFMAN *hm, uint8_t *T, int n, uint32_t **ret_X, int *ret_nX, int *ret_nb)
{
  int       xalloc = 4096;
  uint32_t *X      = malloc(sizeof(uint32_t) * xalloc);
  int       pos      = 0;
  int       nb;
  int       i;

  X[0] = 0;
  nb     = 0;
  for (i = 0; i < n; i++)
    {
      huffman_pack(X, &pos, &nb, hm->code[T[i]], hm->len[T[i]]);
      
      if (pos+1 == xalloc) {
	xalloc *= 2;
	X       = realloc(X, sizeof(uint32_t) * xalloc);
      }
    }
  *ret_X  = X;
  *ret_nX = (pos+1);
  *ret_nb = 32*pos + nb;
  return eslOK;
}


/*****************************************************************
 * 4. Decoding
 *****************************************************************/

/* huffman_unpack()
 * *vp  : ptr to v; v = next 32 bits
 * *X   : encoded input
 *  n   : length of input (in uint32_t)
 * *ip  : current position in <X>
 * *ap  : number of bits left in X[*ip]
 * 
 * If we have to buffer X (say, if we're reading it from
 * a long input) we'll have to redesign. Right now we assume
 * it's just an array.
 */
static void
huffman_unpack(ESL_HUFFMAN *hm, uint32_t *vp, uint32_t *X, int n, int *ip, int *ap, uint8_t *ret_x, int *ret_L)
{
  int  L,D;
  int  idx;

  //printf("Unpacking...\n");
  //dump_uint32(stdout, *vp, 32);               fputc('\n', stdout);

  for (D = 0; D < hm->D-1; D++)
    if ((*vp) < hm->dt_lcode[D+1]) break;
  L = hm->dt_len[D];
  /* L is now the next code's length (prefix of v) */
  //printf("code length = %d\n", L);

  /* Take advantage of lexicographic sort/numerical order of canonical code, within each L */
  idx     = hm->dt_rank[D] +  ( ((*vp) - hm->dt_lcode[D]) >> (eslHUFFMAN_MAXCODE-L) );
  
  /* Now refill v, as much as we can, from bits in X[i] and X[i+1], and update i, a */
  *vp  = (*vp << L);                    // Take L bits from *vp by leftshifting it.
  //dump_uint32(stdout, *vp, 32); fputc('\n', stdout);

  if (*ip < n) {                        // Take either L or all *ap bits from X[i], if it exists.
    *vp   |= (X[*ip] >> (32-L));
    //dump_uint32(stdout, *vp, 32); fputc('\n', stdout);
    X[*ip] =  X[*ip] << L;
    *ap -= L;
    if (*ap == 0)                       // If we exactly finished off X[i]:
      {
	(*ip)++;
	*ap = 32;
      }
    else if (*ap < 0)                   // If we finished off X[i] but still need some bits
      { 
	(*ip)++;                        //   then go on to X[i+1].
	if (*ip < n)                    //   If it exists...
	  {                             //       (...no, I don't like all these branches either...)
	    *ap += 32;                  //     then we're going to leave it w/ <*ap> bits
	    *vp |= (X[*ip] >> *ap);     //     after taking the bits we need to fill v
	    X[*ip] = X[*ip] << (32 - *ap);  //     and leftshifting X[i+1] by that many.
	  }
	else 
	  {
	    *ap = 0;                      //   If X[i+1] doesn't exist, leave *ip = n and *ap = 0; out of data in X (though not necessarily in v)
	  }
      }
  }

  *ret_x  = hm->sorted_at[idx];
  *ret_L  = L;
}

int
esl_huffman_Decode(ESL_HUFFMAN *hm, uint32_t *X, int n, int nb, uint8_t **ret_T, int *ret_nT)
{
  uint8_t *T   = NULL;
  int      allocT;
  uint32_t v   = X[0];
  int      i   = 1;
  int      a   = (n > 1 ? 32 : 0);
  int      pos = 0;
  int      L;
  
  allocT = 4096;
  T      = malloc(sizeof(uint8_t) * allocT);

  while (nb > 0)
    {
      huffman_unpack(hm, &v, X, n, &i, &a, &(T[pos]), &L);
      nb -= L;

      if (++pos == allocT) {
	allocT *= 2;
	T = realloc(T, sizeof(uint8_t) * allocT);
      }
    }
  
  *ret_T  = T;
  *ret_nT = pos;
  return eslOK;
}




/*****************************************************************
 * x. Debugging, development
 *****************************************************************/
    
int
esl_huffman_Dump(FILE *fp, ESL_HUFFMAN *hm)
{
  int r,x;
  int d,L;

  for (r = 0; r < hm->Ku; r++)
    {
      x = hm->sorted_at[r];
      fprintf(fp, "%3d %3d ", x, hm->len[x]);
      dump_uint32(fp, hm->code[x], hm->len[x]);
      fprintf(fp, "\n");
    }
  fputc('\n', fp);


  if (hm->dt_len)
    for (d = 0; d < hm->D; d++)
      {
	L = hm->dt_len[d];
	fprintf(fp, "L=%2d  r=%3d (%3d) ", L, hm->dt_rank[d], hm->sorted_at[hm->dt_rank[d]]);
	dump_uint32(fp, hm->dt_lcode[d], eslHUFFMAN_MAXCODE);
	fputc('\n', fp);
      }

  return eslOK;
}



/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef eslHUFFMAN_EXAMPLE

#include "easel.h"
#include "esl_huffman.h"

#include <stdio.h>
#include <string.h>

/* Given an open <fp> for reading;
 * input text from it and return it as a single string.
 * Optionally return the number of characters in <opt_n>.
 * Convert all \n to spaces.
 */
static uint8_t *
read_text(FILE *fp, int *opt_n)
{
  int   maxlinelen = 4096;
  char *text       = malloc(sizeof(char) * maxlinelen);
  int   n          = 0;
  char *p;

  while (fgets(text+n, maxlinelen-1, fp) != NULL)
    {
      for (p = text+n; *p != '\0'; p++) 
	if (*p == '\n') *p = ' ';
      n   += strlen(text+n);
      text = realloc(text, sizeof(char) * (n+maxlinelen));
    }

  if (opt_n) *opt_n = n;
  return (uint8_t *) text;
}

int
main(int argc, char **argv)
{
  FILE     *fp = fopen(argv[1], "r");
  int       n;
  uint8_t  *T  = read_text(fp, &n);
  uint32_t *X  = NULL;
  float     fq[128];
  int       c,i;
  int       nX, nb;
  ESL_HUFFMAN *hm = NULL;
  uint8_t     *newT = NULL;
  int          nT;

  for (c = 0; c < 128; c++) fq[c]     = 0.;
  for (i = 0; i < n;   i++) fq[T[i]] += 1.;
  
  hm = esl_huffman_Create(fq, 128);
  esl_huffman_Dump(stdout, hm);

  esl_huffman_Encode(hm, T, n, &X, &nX, &nb);

  printf("Original:   %d bytes\n", n);
  printf("Compressed: %d bytes (%d bits)\n", nX*4, nb);

  for (i = 0; i < 30; i++) {
    dump_uint32(stdout, X[i], 32);
    fputc('\n', stdout);
  }

  esl_huffman_Decode(hm, X, nX, nb, &newT, &nT);
  for (i = 0; i < 30; i++)
    fputc(newT[i], stdout);
  fputc('\n', stdout);

  esl_huffman_Destroy(hm);
  free(T);
  fclose(fp);
  return 0;
}
#endif /*eslHUFFMAN_EXAMPLE*/
  
