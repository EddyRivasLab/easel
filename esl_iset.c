/* Find subset of vertices such that no pair is adjacent (independent set)
 * Or find a pair of disjoint subsets X and Y such that no pair 
 * one in X and one in Y are adjacent (bipartite independent pair)
 * (A pair of vertices are adjacent if their corresponding sequences are
 * >t% indentical)
 *
 * Contents:
 *     1. Array tools: print and shuffle
 *     2. Functions for validating independent sets and bipartite independent pairs
 *     3. Random splitting algorithm 
 *     4. Cobalt splitting algorithms
 *     5. Blue splitting algorithms
 */
#include "esl_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_cluster.h"
#include "esl_iset.h"


 /*****************************************************************
  * 1. Array tools: print and shuffle
 *****************************************************************/

/* Function: shuffle_array()
 * Synopsis: Randomly permutes the elements in an array
 * Incept:   SNP, Oct 16, 2020
 *
 * Purpose:  To place the elements of an array in a random order
 *
 * Notes:    Implements Fisher-Yates algorithm
 *
 * Args:     a - array
 *           n - shuffle first n elements of array
 *           r - source of randomness
 *
 * Returns:   <eslOK> on success; 
 *
 * Throws:   error codes of esl_rnd_Roll or ESL_SWAP if applicable
 */


static int
shuffle_array(ESL_RANDOMNESS *r, int a[], int n)
{
    int w;
    while (n > 1) {
    w = esl_rnd_Roll(r, n);	/* shuffling algorithm: swap last elem with w, decrement n. */
    ESL_SWAP(a[w], a[n-1], int);
    n--;
  }
  return eslOK;
}

/* Function: print_array()
 * Synopsis: print the elements in an array
 * Incept:   SNP, Oct 16, 2020
 *
 * Purpose:  Print out array; useful for debugging. Not required for iset algorithms.
 *
 * Args:     a - array
 *           n - print first n elements of array
 */

static void print_array(int array[], int n)
{
   int i;

   for(i = 0; i < n; i++)
      printf("%d ", array[i]);

    printf("\n");
}


/*****************************************************************
 * 2. Functions for validating independent sets and bipartite independent pairs
      (For debugging and unit tests only)
*****************************************************************/

/* Function: check_iset()
 * Synopsis: Verify that a subset of vertices is an independent set
 * Incept:   SNP, Oct 16 2020
 *
 * Purpose:  Given a subset of vertices, verify whether they form an 
 *           independent set 
 *
 * Args:     base        - pointer to array of n fixed-size vertices in graph
 *           n           - number of vertices
 *           size        - size of each vertex element
 *           linkfunc    - pointer to caller's function for defining linked pairs (edges)
 *           param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *           assignments - array of 0/1s; 1 indicates a vertex is in the subset, 0 indicates
 *                         vertex not in the subset
 *
 * Returns:   <eslOK> if the subset is an independent set 
 *
 * Throws:   esl_fatal error if not an independent set
 */

static int check_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments)
{
   int i,j;
   int status;
   int do_link;

   for (i = 0; i < n; i++){
     for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==1) {
            //printf("evaluating pair %d , %d \t", i,j);
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	       if (do_link){
              esl_fatal("FAILED iset test on pair %d , %d \n", i,j);
           }
        }
     }
   }
   return eslOK;

   ERROR:
   return status;

}

/* Function: check_bi_iset()
 * Synopsis: Verify that a pair of subsets of vertices form a bipartite independent
 *           pair
 *
 * Incept:   SNP, Oct 16 2020
 *
 * Purpose:  Given a pair of disjoint subsets of vertices, verify that the pair is a 
 *           bipartite independent pair 
 *
 * Args:     base        - pointer to array of n fixed-size vertices in graph
 *           n           - number of vertices
 *           size        - size of each vertex element
 *           linkfunc    - pointer to caller's function for defining linked pairs (edges)
 *           param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *           assignments - array of 0/1/2s; 1 indicates the vertex is in one subset, 2 indicates a
 *                         the vertex in in the other subset, 0 indicates the vertex is in neither
 *
 * Returns:   <eslOK> if the pair forms a bipartite independent pair 
 *
 * Throws:   esl_fatal error if not a bipartite independent set pair
 */

static int check_bi_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments)
{

   int i,j;
   int status;
   int do_link;

   for (i = 0; i < n; i++){
     for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==2) {
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	         if (do_link){
              esl_fatal("FAILED bi_iset test on pair %d , %d \n", i,j);
           }
        }
          else if (assignments[i]==2 && assignments[j]==1) {
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	         if (do_link){
              esl_fatal("FAILED bi_iset test on pair %d , %d \n", i,j);
           }
        }


     }

   }

   return eslOK;

   ERROR:
   return status;

}

/*****************************************************************
 * 3. Random splitting algorithm
*****************************************************************/

/* Function:  esl_bi_iset_Random()
 * Synopsis:  Random bipartite independent pair algorithm
 * Incept:    SNP,  Oct  16 2020
 *
 * Purpose:   Produces a bipartite independent pair by randomly selecting
 *            vertices for one group and placing all eligible vertices 
 *            in the other group.
 *
 *            For each vertex v:
 *                With probability t_prob, place vertex v in set 1 
 *            For each vertex w not in set 1:
 *                If w is not adjacent to any vertex in set 1, place w in set 2
 *
 *            Two vertices are adjacent if their corresponding sequences are >t% identical
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            assignments - RETURN: assignments to sets (caller provides n*sizeof(int) space)
 *            r           - source of randomness
 *            t_prob      - probability of set 1
 *
 * Returns:   <eslOK> on success; the <assignments[0..nseq-1]> array contains
 *            set indices: 
 *            0 - vertex not in bipartite independent pair
 *            1 - vertex in one set of bipartite independent pair
 *            2 - vertex in other set of bipartite independent pair
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */

int
esl_bi_iset_Random(void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *assignments, ESL_RANDOMNESS *r, double t_prob)
{

  int i,j;			/* indices  */
  int do_link;
  int status;
 
  for (i=0; i<n; i++){
      if (esl_random(r)< t_prob) assignments[i]=1;
      else assignments[i]=2;
  }
    
  for (i=0; i<n; i++){
      if (assignments[i]==2){
          /*check if adjacent to anyone on side 1*/
          for (j = 0; j<n; j++){
              if (assignments[j]==1){
                  if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
                  if (do_link){ /* is adjacent */
                      /* adjacent, break out of loop  */
                      assignments[i]=0;
                      break;
                  }
              }
          }
      }
  }
                  
  //check_bi_iset( base, n, size, linkfunc, param, assignments);
  return eslOK;

  ERROR:
    return status;

}



/*****************************************************************
 * 4. Cobalt splitting algorithms
*****************************************************************/

/* Function:  esl_iset_Cobalt()
 * Synopsis:  Greedy algorithm for independent set, with a random order
 * Incept:    SNP,  Oct  16 2020
 *
 * Purpose:   Produces an independent set.
 *
 *            U= empty set
 *            For each vertex v:
 *                If v is not adjacent to any vertex in U, add v to U
 *            return U
 *
 *            Two vertices are adjacent if their corresponding sequences are >t% identical
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            assignments - RETURN: assignments to sets (caller provides n*sizeof(int) space)
 *            r           - source of randomness
 *
 * Returns:   <eslOK> on success; the <assignments[0..nseq-1]> array contains
 *            set indices: 
 *            0 - vertex not in independent set
 *            1 - vertex is in independent set
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */

int
esl_iset_Cobalt(void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *workspace, int *assignments, ESL_RANDOMNESS *r)
{
  int *a = NULL;		/* array of vertices      */
  int nb, *b = NULL; 		/* array that will contain vertices added to iset; nb is number of vertices added so far*/
  int *c = NULL;		/*  assignments to groups; 1 if in ISET 0 if not */
  int i,j,v;			/* indices and vertices  */
  int do_link;
  int status;
  int adj;

  a = workspace;
  b = workspace + n;
  c = assignments; 

  for (i=0; i<n; i++){
    c[i]=0;
  }

  for (v = 0; v < n; v++) a[v] = n-v-1; /* initialize by putting all vertices into an array*/
  nb = 0;

  /* shuffle  vertices randomly */
  shuffle_array(r, a, n);

  for (j = 0; j<n; j++){

    v = a[j];	/* to decide whether v goes in iset */

    /*check if adjacent to any vertex in b*/
    adj=FALSE;
    for (int i = n; i < n+nb; i++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
          /* adjacent, break out of loop  */
          adj=TRUE;
          break;
	      }
    }

    /* if exited loop early, v is adjacent to a vertex in b, v will not go in iset*/
    if (adj) c[v]=0;

    /*if ran through loop without exiting early, v is not adjacent to any vertex in b, v will go in iset*/
    else{
      c[v]=1;
      b[nb]= v;
      nb++;
    }
  }
   
   // check_iset( base, n, size, linkfunc, param, assignments);
  return eslOK;

  ERROR:
    return status;

}

/* Function:  esl_bi_iset_Cobalt()
 * Synopsis:  Greedy algorithm for bipartite independent set, with a random order
 * Incept:    SNP,  Oct  16 2020
 *
 * Purpose:   Produces an bipartite independent pair
 *
 *            S= empty set
 *            T= empty set
 *            for each vertex v:
 *                With probability 1/2:
 *                    if v is not adjacent to any vertex in S, add v to S
 *                    else, if v is not adjacent to any vertex in T, add v to T
 *                Alternately (with probability 1/2):
 *                    if v is not adjacent to any vertex in T, add v to T
 *                    else, if v is not adjacent to any vertex in S, add v to S
 *            return S,T
 *
 *            Two vertices are adjacent if their corresponding sequences are >t% identical
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            assignments - RETURN: assignments to sets (caller provides n*sizeof(int) space)
 *            r           - source of randomness
 *
 * Returns:   <eslOK> on success; the <assignments[0..nseq-1]> array contains
 *            set indices: 
 *            0 - vertex not in bipartite independent pair
 *            1 - vertex is in set S of bipartite independent pair
 *            2 - vertex is in set T of bipartite independent pair
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */


int
esl_bi_iset_Cobalt(void *base, size_t n, size_t size,
int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
int *workspace, int *assignments, int *ret_larger, ESL_RANDOMNESS *r)
{
  int *a = NULL;    /* array of vertices */
  int nb1, *b1 = NULL;    /* array of vertices to be added to side 1; number of vertices added so far */
  int nb2, *b2 = NULL;    /* array of vertices to be added to side 2; number of vertices added so far */
  int  *c = NULL;    /* assignments; *1 if in b1, 2 if in b2, 0 if in neither */
  int v,i,j;      /* indices and vertices */
  int do_link;
  int status;
  int larger;
  int adj1, adj2;

  a = workspace;
  b1 = workspace + n;
  b2 = workspace + 2*n;
  c = assignments; 

  for (i=0; i<n; i++){
    c[i]=0;
  }

  for (v = 0; v < n; v++) a[v] = n-v-1; /* initialize by putting all vertices into an array*/
  nb1 = 0;
  nb2 = 0;

  /* shuffle  vertices randomly */
  shuffle_array(r, a, n);


  for (j = 0; j<n; j++){

    v = a[j]; /* to decide whether v goes in side1 or side2 or neither */
    double roll = esl_random(r);  /* uniform 0.0 <= x < 1.0 */

    if (roll<=0.5){
      
      /* first try to put v in b2 */
      /*check if v is adjacent to any vertex in b1*/

      adj1=FALSE;
      for (i = n; i < n+nb1; i++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
          /* adjacent, break out of loop  */
          adj1=TRUE;
          break;
        }
      }

      /* if exited loop early, v is adjacent to a vertex in b1, v will not go in b2; try putting v in b1*/
      if (adj1) {

        adj2=FALSE;
        for (int i = 2*n; i < 2*n+nb2; i++){
          if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
          if (do_link){ /* is adjacent */
            /* adjacent, break out of loop  */
             adj2=TRUE;
             break;
          }
        }

        /*  v is adjacent to something in b2, v will not go in b1*/
        if (adj2){
          c[v]=0;
        }
        /*  v is not adjacent to something in b2, v will go in b1*/
        else{
          c[v]=1;
          b1[nb1]= v;
          nb1++;
        }
      }

      /*if ran through loop without exiting early, v is not adjacent to any vertex in b1, v will go in b2*/
      else{
        c[v]=2;
        b2[nb2]= v;
        nb2++;
      }

    }

    /*roll is > 0.5 */
    else{
    /* first try to put v in b1 */
    /*check if v is adjacent to any vertex in b2*/
      adj2=FALSE;
      for (int i = 2*n; i < 2*n+nb2; i++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
        /* adjacent, break out of loop  */
         adj2=TRUE;
         break;
        }
      }

      /* if exited loop early, v is adjacent to a vertex in b2, v will not go in b1; try putting v in b2*/
      if (adj2) {
        adj1=FALSE;
        for (int i = n; i < n+nb1; i++){
          if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
          if (do_link){ /* is adjacent */
            /* adjacent, break out of loop  */
             adj1=TRUE;
             break;
          }
        }

        /*  v is adjacent to something in b1, v will not go in b2*/
        if (adj1){
          c[v]=0;
        }
        /*  v is not adjacent to something in b1, v will go in b2*/
        else{
          c[v]=2;
          b2[nb2]= v;
          nb2++;
        }
      }

      /*if ran through loop without exiting early, v is not adjacent to any vertex in b2, v will go in b1*/
      else{
        c[v]=1;
        b1[nb1]= v;
        nb1++;
      }
    }
  }

  if (nb1>= nb2) larger=1;
  else larger=2;

  *ret_larger=larger;
  
  //check_bi_iset( base, n, size, linkfunc, param, assignments);

  return eslOK;

  ERROR:
  return status;
}


/*****************************************************************
 * 4. Blue splitting algorithms
*****************************************************************/

/* Function:  esl_iset_Blue()
 * Synopsis:  Algorithm for independent set via a multi-round election process
 * Incept:    SNP,  Oct  16 2020
 *
 * Purpose:   Produces an independent set.
 *
 *            U= empty set
 *            L= all vertices
 *            while L is non-empty:
 *                Place vertices of L in a random order v_1,...v_k
 *                Assign each vertex in L a value ~ unif[0,1]
 *                for i=1 to k:
 *                    if label of v_i < label of w for all neighbors w of v_i in L:
 *                        Add v_i to U
 *                        Remove all neighbors of v_i from L
 *            return U
 *
 *            Two vertices are adjacent if their corresponding sequences are >t% identical
 *
 * Notes:     Pseudocode above is given for intuition. Here we implement the following pseudocode
 *            which produces the same result as the above pseudocode with fewer edge queries. 
 *
 *            U= empty set
 *            status_d = dictionary with keys=vertices; all values initially 0
 *                    // keeps track of current status: status_d[i] = -1 if i in iset, -3 if i removed from graph, 
 *                    // >=0 if i still eligible, value is next position of to_add that needs to be checked
 *            k= number of vertices still in graph // will represent number of eligible vertices (i.e. the number
 *                                                 // of vertices for which status_d[v] is non-negative)
 *
 *            while k>0:
 *                dec_o= array with eligible vertices placed in a random order // order in which we will make 
 *                                                                             // decisions about vertices
 *                label_o= array with eligible vertices placed in a random order // instead of labeling the vertices with
 *                                                                               // random values, place them in a random 
 *                                                                               // order representing lowest to highest label
 *                to_add = empty array // array of vertices to be added to iset in this round
 *                lta= 0 //length of to_add array
 *                
 *                Iterate through the vertices v according to dec_o:
 *                    if status_d[v] <0, continue // vertex already removed, nothing to do here
 *
 *                    //first check if v is adjacent to a vertex in to_add                     
 *                    for i=status_d[v] to lta:
 *                        if to_add[i] is adjacent to v:
 *                            status_d[v]=-3 // remove v's eligibility since v is adjacent to vertex in iset
 *                            break
 *                    
 *                    status_d[v]=lta // next time checking for adjacencies between v and vertices in to_add, start at lta
 *                  
 *                    //now check if v is adjacent to a vertex with a lower label
 *                    found_self=False //will become true after v is reached during iteration through label_o
 *                    adj=False //will become true if v is adjcent to a vertex in label_o that is eligible
 *                    j=0// current vertex in label order being evaluated
 *                   
 *                    while !found_self && !adj:
 *                        w=label_o[j]
 *                        if w==v, found_self=TRUE break
 *                        if status_d[w]>=0: //nothing to do in other cases; if status_d[w], w no longer eligible
 *                                           //we already know that v is not adjacent to any vertex in to_add
 *                             if w and v are adjacent:
 *                                // need to check if w is actually eligible (or whether w is ineligible because an adjacency in to_add)
 *                                w_there=True
 *                                for l=status_d[w] to lta:
 *                                    if w and to_add[l]:
 *                                        status_d[w]=-3 // w is not eligible
 *                                        w_there=False
 *                                        break
 *                                if w_there==True:
 *                                    status_d[w]=lta // next time checking for adjacencies between w and vertices 
 *                                                    // in to_add, start at lta
 *                                    adj=True // v is adjacent to a vertex (w) with a lower label
 *                                    break
 *                        j++
 *
 *                    if found_self==True: // v is not adjacent to a vertex with a lower label, can add v to iset!
 *                        to_add[lta]=v
 *                        lta++
 *                        status_d[v]=-1
 *
 *                // Remove eligibility of vertices adjacent to vertices in to_add (This is necesary when a vertex y
 *                // in to_add is adjacent to an eligible vertex x whose fate was decided before y)
 *                for i=0 to k:
 *                    v=dec_o[i] //check whether v is adjacent to any vertex in to_add
 *                    for j=status_d[v] to lta:
 *                        if to_add[j] and v are adjacent:
 *                            status_d[v]=-3 //remove eligibility of v    
 *            
 *                Add all vertices in to_add to U 
 *                Reset status_d[v]=0 for all eligible vertices (i.e. vertices with status_d[v]>=0)
 *                k=number of eligible vertices remaining
 *                Clear dec_o, label_o
 *            
 *            return U
 *                     
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            assignments - RETURN: assignments to sets (caller provides n*sizeof(int) space)
 *            r           - source of randomness
 *
 * Returns:   <eslOK> on success; the <assignments[0..nseq-1]> array contains
 *            set indices: 
 *            0 - vertex not in independent set
 *            1 - vertex is in independent set
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */

int
esl_iset_Blue(void *base, size_t n, size_t size,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
        int *workspace, int *assignments, ESL_RANDOMNESS *r)
{
  int k=n;
  int i;
  int status;
  int lta=0; /*length of to_add*/
  int *dec_o = NULL;  /* vertices in order in which decisions about them will be made */
  int *label_o = NULL;  /* vertices in order of smallest to largest label */
  int *status_d = NULL;   /* keeps track of current status of each vertex like a dictionary; status_d[i] is the status of vertex i */
  /* -1 in iset, -3 removed from graph, >=0 still in graph- value is (most recent index of to_add checked against) + 1 */
  int *to_add= NULL;  /* vertices to add to independent set */


  dec_o = workspace;
  label_o = workspace + n;
  status_d= workspace +2*n;
  to_add= workspace +3*n;

  for (i=0; i<n; i++){
    dec_o[i]=i;
    label_o[i]=i;
    status_d[i]=0;
    assignments[i]=0;
  }


  while (k>0){

    if ((status=i_select(base, n, size, k, linkfunc, param, dec_o, label_o ,status_d, to_add, &lta))!= eslOK) goto ERROR;
    i_update_workspace(dec_o, label_o ,status_d, to_add, assignments, n, &k, &lta, r);

  }

  //check_iset( base, n, size, linkfunc, param, assignments);
  ERROR:
    return status;

  return eslOK;

}


/* Helper function for esl_iset_Blue() that fills to_add */

static int
i_select(void *base, size_t n, size_t size, int k,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
         int *dec_o, int *label_o, int *status_d, int *to_add, int *ret_lta)
{
  int v, w;   /* vertices  */
  int lta=0; /* length of to_add */
  int do_link;
  int status;
  int adj; /* keeps track if adjacency is found*/
  int found_self=FALSE; /* keeps track of whether have found self in label_o*/
  int w_there=FALSE; /* keeps track of whether a vertex is still in graph*/
  int i,j,l; /*indices for for loops*/

  for(i=0; i<k; i++){

    v=dec_o[i]; /* decide fate of this vertex v*/

    /* if vertex v has already been removed, nothing to decide so skip this iteration */
    if (status_d[v]<0) continue;

    /* check if adjacent to any vertex in to_add*/
    adj=FALSE;
    for (j=status_d[v]; j< lta; j++){
      /* if v is adjacent to to_add[j], remove v from graph and break out of loop */
      if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + to_add[j]*size, param, &do_link)) != eslOK) goto ERROR;
      if (do_link){ /* is adjacent */
        status_d[v]= -3;
        adj=TRUE;
        break;
      }

    }

    /* if not adjacent to any vertex in to_add, check if v is not adjacent to any vertex with a lower label*/
    if (!adj){

      status_d[v]= lta; /* next time check adjacencies between v and to_add, start at index lta */
      adj=FALSE; /* becomes true when v is determined to be adjacent to a vertex with a lower label that is still in the graph*/
      found_self=FALSE; /* becomes true when v is reached in label_o */

      /* iterate through label_o until find v or find that v is adjacent to a vertex with a lower label that is still in graph*/
      j=0;
      while (!found_self && !adj){
        w=label_o[j]; /*w is a vertex with a lower label than v's label*/

        /* check if w is v, if so break*/
        if (w==v){
          found_self=TRUE;
          break;
        }

        if (status_d[w] >=0 ){

          /* check whether w and v are adjacent */
          if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
          if (do_link){ /* is adjacent */

            /* check whether w should really be in the graph*/
            w_there=TRUE;
            for(l=status_d[w]; l< lta; l++){
              if ((status = (*linkfunc)( (char *) base + w*size, (char *) base + to_add[l]*size, param, &do_link)) != eslOK) goto ERROR;
              if (do_link){
                /* remove w from graph*/
                status_d[w]=-3;
                w_there=FALSE;
                break;
              }
            }

            /* w is in the graph and so v does not get added to iset (since it is adjacent to w, which is a vertex in the graph with a lower label)*/
            if (w_there){
              status_d[w]= lta; /* next time check adjacencies between w and to_add, start at index lta */
              adj=TRUE; /* v is adjacent to w, which has a lower label */
              break;
            }
          }
        }
        j++;
      }

      /* if v is not adjacent to any vertex with a lower label, v should be added to iset */
      if (found_self){
        to_add[lta]=v;
        lta++;
        status_d[v]=-1;
      }
    }
  }

  /* check if vertices are adjacent to vertices in to_add that were added after them, if so remove vertex*/
  for(i=0; i<k; i++){
   
    v=dec_o[i]; /* vertex to check v*/
   
    if (status_d[v]>=0){
      for (j=status_d[v]; j< lta; j++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + to_add[j]*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){
          /* remove v from graph*/
          status_d[v]=-3;
          break;
        }
      }
    }
  }

  *ret_lta=lta;
  return eslOK;

 ERROR:
   return status;
}


/* Helper function for esl_iset_Blue() that resets dec_o, label_o, and status_d */

static void
i_update_workspace(int *dec_o, int *label_o, int *status_d, int *to_add, int *assignments, size_t n, int *k, int *lta, ESL_RANDOMNESS *r){
  int i;
  int d=0;


  /* add all vertices in to_add to iset and clear to_add*/
  for (i=0; i<*lta; i++){
    assignments[to_add[i]]=1;
    to_add[i]=0;
  }

  /* clear decison order */
  for (i=0; i<*k; i++){
    dec_o[i]=0;
  }

  /*put all vertices left in graph (i.e. status order is >=0 and is in label order) into decision order */
  for (i=0; i<*k; i++){
    if (status_d[label_o[i]]>=0){
      dec_o[d]=label_o[i];
      d++;
      status_d[label_o[i]]=0;
    }
    if (status_d[label_o[i]]==-3){
      assignments[label_o[i]]=0;
    }
    label_o[i]=0;
  }

  /* copy decision order to label order */
  for (i=0; i<d; i++){
    label_o[i]=dec_o[i];
  }

  /*shuffle label_o and dec_o */
  shuffle_array(r, dec_o, d);
  shuffle_array(r, label_o, d);


  *k=d;
  *lta=0;


}


/* Function:  esl_bi_iset_Blue()
 * Synopsis:  Algorithm for bipartite independent pair via a multi-round election process
 * Incept:    SNP,  Oct  16 2020
 *
 * Purpose:   Produces an bipartite independent pair.
 *
 *            S,T= empty set
 *            L_S, L_T= all vertices // represent eligibility for S and T sets 
 *            while L_T or L_S is non-empty:
 *                
 *                C_S, C_T = empty set // represents S-candidates and T-candidates for the round
 *                For each vertex that is in L_S, but not in L_T, add the vertex to C_S
 *                For each vertex that is in L_T, but not in L_S, add the vertex to C_T
 *                For each vertex in both L_S and L_T, randomly assign to C_S (exclusive) or C_T
 *
 *                
 *                Place vertices of C_S in a random order v_1,...v_k
 *                Assign each vertex in L a value ~ unif[0,1]
 *                For i=1 to k:
 *                    If label of v_i < label of w for all neighbors w of v_i in both L_T and C_T:
 *                        Add v_i to S
 *                        Remove all neighbors of v_i from L_T
 *                        Remove v_i from L_T and L_S
 *                
 *                For all vertices w in both L_T and C_T:
 *                    Add w to T, remove w from L_T, remove w and all of its neighbors from L_S
 *                
 *            Return S, T
 *
 *            Two vertices are adjacent if their corresponding sequences are >t% identical
 *
 * Notes:     Pseudocode above is given for intuition. Here we implement the following pseudocode
 *            which produces the same result as the above pseudocode with fewer edge queries. 
 *            
 *            S,T=empty set
 *
 *            elig= dictionary with key=vertices; all values initially 3 // keeps track of eligibility of the vertices
 *                  // 0 removed from graph (in one side of iset or disqualified because no longer eligibile for either side)
 *                  // 1 eligibile for 1 only, 2 eligible for 2 only, 3 eligible for both 1 and 2 
 *            
 *            while some vertices are still eligible: 
                d=0 // number of 1-side candidates
 *              l=0 // number of 2-side candidate
 *              dec_o= empty array of length number of total vertices // will store 1-candidates on left side, and 2-candidates on right side
 *              // assign candidacy to vertices
 *              for each vertex:
 *                  if v is 1-eligible and not 2-eligible:
 *                      // make v 1-candidate by placing v in next open spot of dec_o on left side
 *                      dec_o[d]=v
 *                      d++
 *                  if v is 2-eligible and not 1-eligible:
 *                      // make v a 2-candidate by placing v in next open spot of dec_o on right side
 *                      dec_o=[n-1-l]=v
 *                      l++
 *                  if v is 1-eligible and 2-eligible:
 *                      with prob 1/2 make v 1-candidate by placing v in the next open spot of dec_o on left side (dec_o[d]=v, d++)
 *                      alternately make v a 2-candidate by placing v in the next open spot of dec_o on right side (dec_o=[n-1-l]=v, l++)
 *            
 *              label_o= dec_o+n-l // to avoid indexing from the right, we break off the dec_o array into an array called label_o
 *              status_d=dictionary with keys=vertices; 
 *                        //if v is a 1-candidate, v is before position status_d[i] in label order (there is never a need to compare
 *                        //the labels of two 1-candidates, so it is fine if multiple 1-candidates are in the same position of the label order) 
 *                        //if v is a 2-candidate, value of status_d[v] is next position of to_add that needs to be checked 
 *              for all 2-candidates, initialize status_d[v]=0
 *              for all 1-candidates, status_d[v]= random integer in [0, l+1]
 *
 *              shuffle 1-candidates within dec_0 (i.e. apply random permutation to first d elements of dec_o)
 *              shuffle 2-candidates within label_0 (i.e. apply random permutation to the l elements of label_o)
 *
 *             // elect vertices to 1-side
 *             lta1=0, lta2=0 // length of to_add for each of the sides
 *             to_add=empty array // left side will store vertices added to 1-side, right side will store vertices added to 2-side
 *             Iterate through the vertices v according to dec_o (only the 1-candidates):
 *                  should_add=True
 *                  for j=0 to status_d[v]: //check if v is adjacent to a vertex that is a 2-candidate with a lower label
 *                      w=label_o[j]
 *                      update_2_elig(j) // see helper function below, updates 2-elig of label_o[j]
 *                      if elig[w]==2 or 3 // then w is 2-eligible
 *                          if w and v are adjacent:
 *                              status_d[v]=j // keep track that v has been compared to first j vertices in label order (new meaning)
 *                              should_add=False
 *                              break
 *                  if should_add:
 *                      to_add[lta1]=v; lta++; elig[v]=0
 *
 *              // add vertices to 2-side
 *              for j=0 to l:
 *                  w=label_o[j] // 2-candidate vertex to decide 
 *                  update_2_elig(j) // see helper function below, updates 2-elig of label_o[j]
 *               
 *                  if elig[w]==2 or 3 // then w is 2-eligible so add w to 2-side
 *                      to_add[n-1-lta2]=w; lta2++; elig[w]=0
 *                      
 *                      //remove eligibility of 1-side candidates that are adjacent to w
 *                      for i=0 to d:
 *                          v=dec_o[d] // decide 2-eligibility of v
 *                          if elig[v]==1 or 3:  
 *                              // since v was not added to 1-side, status_d[v] has new meaning, v is not adjacent to first j-1 vertices                       
 *                              // of label_o and is adjacent to label_o[j]
 *                              if status_d[v]==j, elig[v]=elig[v]-1 // w is adj to v, w is on 2-side, so must remove 1-eligibility of v
 *                              else if status_d[v]< j:
 *                                  if v and w are adjacent, elig[v]=elig[v]-1 // w is adj to v, w is on 2-side, so must remove 1-eligibility of v
 *
 *              // remove 1-elig of 2-candidates that are adjacent ot a vertex in 2-side of to-add
 *              for j=0 to l:
 *                  w=label_o[j]
 *                  if elig[w]==1 or 3:
 *                      for k=n-1 to n-1-lta2: //iterate through vertices added to 2-side
 *                          u=to_add[k]
 *                          if u and w are adjacent, elig[w]=elig[w]-1, break // w is adajcent to u (which is in 2-side), so w is no longer 1-eligible
 *
 *              // remove 2-elig of 1-candidates that are adjacent to a vertex in 1-side of to_add
 *              for i=0 to d:
 *                  v=dec_o[i]
 *                  if elig[v]==2 or 3: //iterate through vertices added to 1-side
 *                      for k=0 to lta1:
 *                          u=to_add[k]
 *                            if u and v are adjacent, elig[v]=elig[v]-2 // v is adjacent to u (which is in 1-side) so must remove 2-elig of v
 *            
 *              add vertices on left side of to_add (positions 0 to lta1) to S
 *              add vertices on right side of to add (position n-lta2 to n-1) to T
 *        
 *          return S,T
 *
 *          PSEUDOCODE for update_2_elig
 *          
 *          update_2_elig(j):
 *              w=label_o[j]
 *              if elig[w]==2 or 3:
 *                  //check 1-side of to_add for adjacencies with w
 *                  for i=status_d[w] to lta:
 *                       v=to_add[i]
 *                       // if v has a higher label that j, v and u were already compared and determined to be non-adjacent before v was added 
 *                      if status_d[v]<=j:
 *                          if v and w are adjacent, elig[w]=elig[w]-2, status_d[w]=lta1, break // w is not 2 eligible because it is adjacent 
 *                                                                                              // to v, which is in 1-side of to_add
 *
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            assignments - RETURN: assignments to sets (caller provides n*sizeof(int) space)
 *            r           - source of randomness
 *
 * Returns:   <eslOK> on success; the <assignments[0..nseq-1]> array contains
 *            set indices: 
 *            0 - vertex not in bipartite independent pair
 *            1 - vertex is in set S of bipartite independent pair
 *            2 - vertex is in set T of bipartite independent pair
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */

int
esl_bi_iset_Blue(void *base, size_t n, size_t size,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
        int *workspace, int *assignments, int *ret_larger, ESL_RANDOMNESS *r)
{
  int i;
  int status;
  int larger;
  int nb1=0; /* number of vertices selected for 1-side*/
  int nb2=0; /* number of vertices selected for 2-side*/
  int d=0; /* number of 1-side candidates */
  int l=0; /* number of 2-side candidate*/
  int lta1=0, lta2=0; /*length of to_add*/
  int *dec_o = NULL;  /* vertices in order in which decisions about them will be made */
  int *label_o = NULL;  /* vertices in order of smallest to largest label */
  int *status_d = NULL;   /* keeps track of current status of each vertex like a dictionary; status_d[i] is the status of vertex i */
  /* if v is a 1-candidate v is before status_d[i] in label order; if v is a 2-candidate status_d keeps track of most recent member of to_add compared to + 1 */
  int *to_add= NULL;  /* vertices to add to independent set */
  int *elig= NULL;  /* dictionary to keep track of eligibility of the vertices */
  /* 0 removed from graph (in one side of iset or disqualified because no longer eligibile for either side), 1 eligibile for 1 only, 2 eligible for 2 only, 3 eligible for both 1 and 2 */


  label_o=workspace;
  dec_o = workspace;
  status_d= workspace +n;
  to_add= workspace +2*n;
  elig= workspace +3*n;

  /* initialize assignments to avoid funny business; assignments should not have 0 or 1*/
  /* initialize to_add */
  for (i=0; i<n; i++){
    assignments[i]=0;
    to_add[i]=-1; /*should never try to add vertex -1*/

  }

  
  for (i=0; i<n; i++){
    /* all vertices initially eligible for both sides */
    elig[i]=3;
  }

  bi_update_workspace_blue(dec_o, label_o, status_d, to_add, elig, assignments, (int) n, &d, &l, &lta1, &lta2, &nb1, &nb2, r);
  label_o=dec_o+n-l;   
  
  while (l+d>0){
      
    if ((status=bi_select_blue(base, (int) n, size, linkfunc, param, dec_o, label_o ,status_d, to_add, elig, d, l, &lta1, &lta2))!= eslOK) goto ERROR;
      
    bi_update_workspace_blue(dec_o, label_o, status_d, to_add, elig, assignments, (int) n, &d, &l, &lta1,&lta2, &nb1, &nb2, r);
    label_o=dec_o+n-l;  
      
  }


  if (nb1>= nb2) larger=1;
  else larger=2;

  *ret_larger=larger;

  return eslOK;
  
  //check_bi_iset( base, n, size, linkfunc, param, assignments);
  ERROR:
    return status;


}

/* Helper function for esl_bi_iset_Blue() */

static void
bi_update_workspace_blue(int *dec_o, int *label_o, int *status_d, int *to_add, int *elig, int *assignments, int n, int *d, int *l, int *lta1, int *lta2, int *nb1, int *nb2, ESL_RANDOMNESS *r){
  
  *d=0;
  *l=0;
  int i;

  /* add all vertices on left side of to_add to side 1 and clear to_add*/
  for (i=0; i<*lta1; i++){
    assignments[to_add[i]]=1;
    (*nb1)++;
    to_add[i]=-1;
  }

  /* add all vertices on right side of to_add to side 2 and clear to_add*/
  for (i=n-1; i>= n-*lta2; i--){
    assignments[to_add[i]]=2;
    (*nb2)++;
    to_add[i]=-1;
  }


  for (i=0; i<n; i++){
    if (elig[i]==3) {
      /* randomly assign to a side*/
      if (esl_random(r)< .5){
      /* vertex is a 1-candidate, put into left side of order */
        dec_o[*d]=i;
        (*d)++;
      } 
      else{
       /* vertex is a 2-candidate, put into right side of order */
        dec_o[n-1-*l]=i;
        (*l)++;
        status_d[i]=0;
      }
    }

    else if (elig[i]==1){
        /* vertex is a 1-candidate, put into left side of order */
        dec_o[*d]=i;
        (*d)++;
    }

    else if (elig[i]==2){
       /* vertex is a 2-candidate, put into right side of order */
        dec_o[n-1-*l]=i;
        (*l)++;
        status_d[i]=0;

    }
  }

  /* right side of dec_o is the label order*/
  label_o=dec_o+n-*l;

  /*place 1-side candidates in label order*/
  for (i=0; i<*d; i++){
    /* vertex v is before this position in label order */
    /* choose random value between 0 and l inclusive */
    status_d[dec_o[i]]=esl_rnd_Roll(r, (*l)+1);
  }

  shuffle_array(r, dec_o, *d);
  shuffle_array(r, label_o, *l);


  *lta1=0;
  *lta2=0;
  
}



/* Helper function for esl_bi_iset_Blue() */

static int
update_2_elig(int j, void *base, int n, size_t size,
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,int *label_o, int *status_d, int *to_add, int *elig, const int lta1)
{
  int w,v,i;
  int do_link;
  int status;
 
  
  w=label_o[j];
  
  /* not 2-eligible, nothing to do */
  if (elig[w]==2 || elig[w]==3) {

    /* check 1- side of to_add for adjacencies with w */
    for (i=status_d[w]; i< lta1; i++){
      v=to_add[i];
      /* if v has a higher label than j, v and u were already compared and determined to be non-adjacent before v was added to_add */
      if (status_d[v]<=j){
      if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
            elig[w]=elig[w]-2;
            status_d[w]=lta1;
            break;
        }
      }
    }
  }


  return eslOK;

  ERROR:
    return status;

}

static int
bi_select_blue(void *base, int n, size_t size, 
        int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
         int *dec_o, int *label_o, int *status_d, int *to_add, int *elig, const int d, const int l, int *ret_lta1, int *ret_lta2)
{

  int v, w, u;     /* vertices  */
  int lta1=0, lta2=0; /* length of to_add */
  int do_link;
  int status;
  int i,j,k; /*indices for for loops*/
  int should_add;

  /* select 1-candidates for 1-side */
  /* iterate over 1-candidates, all of which are in dec_o*/
  for(i=0; i<d; i++){

    v=dec_o[i]; /* decide fate of this vertex v*/
    /* iterate over 2-candidate vertices that have a smaller label than v */
    should_add=TRUE;
    for (j=0; j< status_d[v]; j++){
      
      update_2_elig(j, base, n, size, linkfunc, param, label_o, status_d, to_add, elig, lta1); /*update eligibility of w*/
      w=label_o[j];

      /* if w is still 2-eligible and is adjacent to v, v should not be added*/
      if (elig[w]==2 || elig[w]==3){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
            status_d[v]= j; /*keep track that v got up to j in label order*/
            should_add=FALSE;
            break;
        }
      }
    }

    if (should_add){
      to_add[lta1]=v;
      lta1++;
      elig[v]=0;
    }

  }

  /* select 2-candidates for 2-side */
  /* iterate over 2-candidates, all of which are in label_o*/
  for (j=0; j<l; j++){
    
    update_2_elig(j, base, n, size, linkfunc, param, label_o, status_d, to_add, elig, lta1);
    w=label_o[j]; /* decide whether w goes into 2-side*/
    if (elig[w]==2 || elig[w]==3){
      /* add to 2-side*/
      to_add[n-1-lta2]=w;
      lta2++;
      elig[w]=0;
      
      /* remove 1-side candidates that are adjacent to w */
      for (i=0; i<d; i++){
        v=dec_o[i];
        if (elig[v]==0) continue;
        if (elig[v]==1 || elig[v]==3){
          /* since v was not added, status_d[v] represents last position in label order that v was compared to and v is adjacent to vertex at that position */
          if (status_d[v]==j) elig[v]=elig[v]-1; /* v is adjacent to w, which was just added to 2-side, so remove 1-elig of v */
          /* v only checked adjacencies up to label_o[status_d[v]] and so v and w=label_o[j] have never been compared */
          else if (status_d[v] < j){
            if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
            if (do_link) elig[v]=elig[v]-1; /* v is adjacent to w, which was just added to 2-side, so remove 1-elig of v */
          }
        }
      }
    }
  }

  /* remove 1-elig of 2-candidates that are adjacent to a vertex in 2-side of to_add */
  for (j=0; j<l; j++){
    w=label_o[j];
    if (elig[w]==1 || elig[w]==3){
      for (k=n-1; k> n-1-lta2; k--){
        u=to_add[k];
        if ((status = (*linkfunc)( (char *) base + u*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ 
          elig[w]=elig[w]-1; /* w is adjacent to u, which is in 2-side, so remove 1-elig of w */
          break;
        }
      }
    }
  }

  /* remove 2-elig of 1-candidates that are adjacent to a vertex in 1-side of to_add */
  for (i=0; i<d; i++){
    v=dec_o[i];
    if (elig[v]==2 || elig[v]==3){
      for (k=0; k< lta1; k++){
        u=to_add[k];
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + u*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ 
          elig[v]=elig[v]-2; /* v is adjacent to u, which is in 1-side, so remove 2-elig of v */
          break;
        }
      }
    }
  }

  *ret_lta1=lta1;
  *ret_lta2=lta2;

 return eslOK;

 ERROR:
   //printf("in error\n");
   return status;
}

/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef eslISET_TESTDRIVE
/* gcc -g -Wall -o msacluster_utest -I. -L. -DeslMSACLUSTER_TESTDRIVE esl_msacluster.c -leasel -lm
 */
#include "esl_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msafile.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for iset module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  ESL_MSA        *msa     = esl_msa_CreateFromString("\
# STOCKHOLM 1.0\n\
\n\
seq0  AAAAAAAAAA\n\
seq1  AAAAAAAAAA\n\
seq2  AAAAAAAAAC\n\
seq3  AAAAAAAADD\n\
seq4  AAAAAAAEEE\n\
seq5  AAAAAAFFFF\n\
seq6  AAAAAGGGGG\n\
seq7  AAAAHHHHHH\n\
seq8  AAAIIIIIII\n\
seq9  AAKKKKKKKK\n\
seq10 ALLLLLLLLL\n\
seq11 MMMMMMMMMM\n\
//",   eslMSAFILE_STOCKHOLM);



  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /* eslISET_TESTDRIVE*/

