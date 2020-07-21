/* Generalized single linkage clustering.
 *
 * Contents:
 *     1. Single linkage clustering, generalized
 *     2. Unit tests
 *     3. Test driver
 *     4. Example
 */
#include "esl_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_cluster.h"
#include "esl_iset.h"


/*****************************************************************
 * 1. Single linkage clustering, generalized
 *****************************************************************/

/* Function:  esl_cluster_SingleLinkage()
 * Synopsis:  Generalized single linkage clustering.
 * Incept:    SRE, Mon Jan  7 08:35:10 2008 [Janelia]
 *
 * Purpose:   Given a set of vertices, cluster them by single-linkage
 *            clustering.
 *
 *            The data describing each vertex is provided in an array
 *            starting at <base>, consisting of <n> vertices. Each
 *            vertex can be of any type (structure, scalar, pointer)
 *            so long as each vertex element is of fixed size <n>
 *            bytes.
 *
 *            A pointer to the clustering function is provided in
 *            <(*linkfunc)()>, and a pointer to any necessary
 *            parameters for that function (for example, any
 *            thresholds) is provided in <param>.
 *
 *            The <int (*linkfunc)()> must be written by the
 *            caller. It takes arguments <(void *v1, void *v2, void
 *            *param, int *ret_link)>: pointers to two vertices to
 *            test for linkage and a pointer to any necessary
 *            parameters, and it passes the answer <TRUE> (1) or
 *            <FALSE> (0) back in <*ret_link>. The <(*linkfunc)()>
 *            returns <eslOK> (0) on success, and a nonzero error code
 *            on failure (see <easel.h> for a list of Easel's error
 *            codes).
 *
 *            The caller provides an allocated <workspace> with space
 *            for at least <2n> integers. (Allocation in the caller
 *            allows the caller to reuse memory and save
 *            allocation/free cycles, if it has many rounds of
 *            clustering to do.)
 *
 *            The caller also provides allocated space in
 *            <assignments> for <n> integers which, upon successful
 *            return, contains assignments of the <0..n-1> vertices to
 *            <0..C-1> clusters. That is, if <assignments[42] = 1>,
 *            that means vertex 42 is assigned to cluster 1.  The
 *            total number of clusters is returned in <ret_C>.
 *
 *            The algorithm runs in $O(N)$ memory; importantly, it
 *            does not require a $O(N^2)$ adjacency matrix. Worst case
 *            time complexity is $O(N^2)$ (multiplied by any
 *            additional complexity in the <(*linkfunc()> itself), but
 *            the worst case (no links at all; <C=n> clusters) should
 *            be unusual. More typically, time scales as about $N \log
 *            N$. Best case is $N$, for a completely connected graph
 *            in which all vertices group into one cluster. (More
 *            precisely, best case complexity arises when vertex 0 is
 *            connected to all other <n-1> vertices.)
 *
 * Notes:    I don't know if this algorithm is published. I
 *           haven't seen it in graph theory books, but that might
 *           be because it's so obvious that nobody's bothered.
 *
 *           In brief, we're going to do a breadth-first search of the
 *           graph, and we're going to calculate links on the fly
 *           rather than precalculating them into a standard adjacency
 *           matrix.
 *
 *           While working, we keep two stacks of maximum length N:
 *                a : list of vertices that are still unconnected.
 *                b : list of vertices that we've connected to
 *                    in our current breadth level, but we haven't
 *                    yet tested for other connections to a.
 *           The current length (number of elements in) a and b are
 *           kept in na, nb.
 *
 *           We store our results in an array of length N:
 *                c : assigns each vertex to a component. for example
 *                    c[4] = 1 means that vertex 4 is in component 1.
 *                    nc is the number of components. Components
 *                    are numbered from 0 to nc-1. We return c and nc
 *                    to our caller.
 *
 *           The algorithm is:
 *
 *           Initialisation:
 *                a  <-- all the vertices
 *                na <-- N
 *                b  <-- empty set
 *                nb <-- 0
 *                nc <-- 0
 *
 *           Then:
 *                while (a is not empty)
 *                  pop a vertex off a, push onto b
 *                  while (b is not empty)
 *                    pop vertex v off b
 *                    assign c[v] = nc
 *                    for each vertex w in a:
 *                       compare v,w. If w is linked to v, remove w
 *                       from a, push onto b.
 *                  nc++
 *           q.e.d.
 *
 * Args:      base        - pointer to array of n fixed-size vertices to be clustered.
 *            n           - number of vertices
 *            size        - size of each vertex element
 *            linkfunc    - pointer to caller's function for defining linked pairs
 *            param       - pointer to any data that needs to be provided to <(*linkfunc)>
 *            workspace   - caller provides at least 2n*sizeof(int) of workspace
 *            assignments - RETURN: assignments to clusters (caller provides n*sizeof(int) space)
 *            ret_C       - RETURN: number of clusters
 *
 * Returns:   <eslOK> on success; <assignments[0..n-1]> contains cluster assigments
 *            <0..C-1> for each vertex, and <*ret_C> contains the number of clusters
 *            <C>
 *
 * Throws:    status codes from the caller's <(*linkfunc)> on failure; in this case,
 *            the contents of <*assignments> is undefined, and <*ret_C> is 0.
 */



 /*****************************************************************
  * Array tools: print and shuffle
 *****************************************************************/

int
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



void print_array(int array[], int n)
{
    int i;

   for(i = 0; i < n; i++)
      printf("%d ", array[i]);

    printf("\n");
}


/*****************************************************************
 * Functions to check accuracy of sets found (for debugging)
*****************************************************************/

static int check_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments)
{
   //printf("check iset called\n");
   //printf("param is %lf\n", *(double *) param);
   int i,j;
   int status;
   int do_link;
   int passed=1;
   //print_array(assignments, n);
   for (i = 0; i < n; i++){
     for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==1) {
            //printf("evaluating pair %d , %d \t", i,j);
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	       if (do_link){
              printf("FAILED iset test on pair %d , %d \n", i,j);
              passed=0;
           }
           //else printf("pair ok\n");

        }


     }

   }
    if (passed==1) printf("PASSED iset test!\n");

     ERROR:
   return status;



}


static int check_bi_iset( void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *assignments)
{
   //printf("check iset called\n");
   //printf("param is %lf\n", *(double *) param);
   int i,j;
   int status;
   int do_link;
   int passed=1;
   //print_array(assignments, n);
   for (i = 0; i < n; i++){
     for (j= i+1; j<n; j++){
        if (assignments[i]==1 && assignments[j]==2) {
            //printf("evaluating pair %d , %d \t", i,j);
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	       if (do_link){
              printf("FAILED bi_iset test on pair %d , %d \t", i,j);
              passed=0;
           }
        }
          else if (assignments[i]==2 && assignments[j]==1) {
            //printf("evaluating pair %d , %d \t", i,j);
           if ((status = (*linkfunc)( (char *) base + j*size, (char *) base + i*size, param, &do_link)) != eslOK) goto ERROR;
	       if (do_link){
              printf("FAILED bi_iset test on pair %d , %d \n", i,j);
              passed=0;
           }
           //else printf("pair ok\n");

        }


     }

   }
    if (passed==1) printf("PASSED bi_iset test!\n");

     ERROR:
   return status;

}


/*****************************************************************
 * Cobalt algorithm
*****************************************************************/

int
esl_iset_Cobalt(void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *workspace, int *assignments, ESL_RANDOMNESS *r)
{
  int *a = NULL;		/* stack of available vertices (still unconnected)       */
  int nb, *b = NULL; 		/* stack of connected but unextended vertices            */
  int *c = NULL;		/*  assignments to groups */
  int v;			/* indices of vertices                                   */
  int do_link;
  int status;


  a = workspace;
  b = workspace + n;
  c = assignments; /*1 if in ISET 0 if not*/


   //printf("starting cobalt iset construction\n");

   for (v = 0; v < n; v++) a[v] = n-v-1; /* initialize by pushing all vertices into an array*/
   nb = 0;

   // print_array(a, n);
    /* shuffle  vertices randomly */

   shuffle_array(r, a, n);

  // print_array(a, n);

    int adj;

   for (int j = 0; j<n; j++){
       //print_array(b,n);
       v = a[j];	/* to decide whether v goes in iset */

     /*check if adjacent to any vertex in b*/
       adj=FALSE;
       for (int i = n; i < n+nb; i++){
           if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
	       if (do_link) /* is adjacent */
		 {
              /* adjacent, break out of loop  */
               adj=TRUE;
               break;
		 }
       }

       /* if exited loop early, v is adjacent to a vertex in b, v will not go in iset*/
        if (adj)
            c[v]=0;


    /*if ran through loop without exiting early, v is not adjacent to any vertex in b, v will go in iset*/
        else{
            c[v]=1;
            b[nb]= v;
            nb++;
        }


     }
    //esl_iset_Cobalt((void *) msa->aseq, (size_t) msa->nseq, sizeof(char *),
				      // msacluster_clinkage, (void *) &maxid,
				       //workspace, assignment, r);
    //printf("param is %lf \n", *(double *) param);

    //check_iset( base, n, size, linkfunc, param, assignments);
    return eslOK;

 ERROR:
   return status;
}



int
esl_bi_iset_Cobalt(void *base, size_t n, size_t size,
int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
int *workspace, int *assignments, int *ret_larger, ESL_RANDOMNESS *r)
{
  int *a = NULL;    /* stack of available vertices (still unconnected)       */
  int nb1, *b1 = NULL;    /* stack of connected but unextended vertices            */
  int nb2, *b2 = NULL;    /* stack of connected but unextended vertices            */
  int  *c = NULL;    /* array of results: # clusters, assignments to clusters */
  int v;      /* indices of vertices                                   */
  int do_link;
  int status;
  int larger;


  a = workspace;
  b1 = workspace + n;
  b2 = workspace + 2*n;
  c = assignments; /*1 if in b1, 2 if in b2, 0 if in neither*/




  //printf("called esl_bi_iset_Cobalt \n");

  for (v = 0; v < n; v++) a[v] = n-v-1; /* initialize by pushing all vertices into an array*/
  nb1 = 0;
  nb2 = 0;

  /* shuffle  vertices randomly */

  shuffle_array(r, a, n);


  int adj1, adj2;
    //printf("this is a\n")
  //print_array(a,n);

  for (int j = 0; j<n; j++){
    // print_array(b1,n);
   // print_array(b2,n);

    v = a[j]; /* to decide whether v goes in iset */
    //printf("current vertex %d\n", v);
    double roll = esl_random(r);  /* uniform 0.0 <= x < 1.0 */

    if (roll<=0.5){
      /* first try to put v in b2 */
      /*check if v is adjacent to any vertex in b1*/
      //printf("trying to add to b2\n");

      adj1=FALSE;
      for (int i = n; i < n+nb1; i++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + a[i]*size, param, &do_link)) != eslOK) goto ERROR;
        if (do_link){ /* is adjacent */
          /* adjacent, break out of loop  */
          adj1=TRUE;
          break;
        }
      }


      /* if exited loop early, v is adjacent to a vertex in b1, v will not go in b2; try putting v in b1*/
      if (adj1) {
        //printf("v adjacent to vertex in b1 \n");

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
 * Cyan algorithm
*****************************************************************/


/* workspace should have size 4*n*/


static int
i_select(void *base, size_t n, size_t size, int k,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			   int *dec_o, int *label_o, int *status_d, int *to_add, int *ret_lta)
{

  int v, w;			/* vertices  */
  int lta=0; /* length of to_add */
  int do_link;
  int status;
  int adj; /* keeps track if adjacency is found*/
  int found_self=FALSE; /* keeps track of whether have found self in label_o*/
  int w_there=FALSE; /* keeps track of whether have found self in label_o*/
  int i,j,l; /*indices for for loops*/


  //printf("i select called \n");
  for(i=0; i<k; i++){
    v=dec_o[i]; /* decide fate of this vertex v*/

    /* if vertex v has already been removed, nothing to decide so skip this iteration */
    if (status_d[v]<0) continue;

    /* check if adjacent to any vertex in to_add*/
    adj=FALSE;
    for (j=0; j< lta; j++){
      /* if v is adjacent to to_add[j], remove v from graph and break out of loop*/
      //printf("checking adjacencies with vertices in to_add\n");
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

        /*only compare v and w if w is still in graph */
        if (status_d[w] >=0 ){

          /* check whether w and v are adjacent */
          //printf("checking w and v adjacency\n");

          if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + w*size, param, &do_link)) != eslOK) goto ERROR;
          if (do_link){ /* is adjacent */

            /* check whether w should really be in the graph*/
            w_there=TRUE;
            for(l=status_d[w]; l< lta; l++){
              //printf("checking adjacencies to see if w still in graph\n");

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
              adj=TRUE;
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
  //printf("this is k\t%d\n", k);
  for(i=0; i<k; i++){
    //printf("checking adjacencies with vertices added later\n");
    v=dec_o[i]; /* vertex to check v*/
    //printf("this is v\t%d\n", v);
    //printf("this is status_d[v]\t%d\n", status_d[v]);
    if (status_d[v]>=0){
      for (j=status_d[v]; j< lta; j++){
        if ((status = (*linkfunc)( (char *) base + v*size, (char *) base + to_add[j]*size, param, &do_link)) != eslOK){
      //    printf("sending to error\n");
          goto ERROR;
        }
        if (do_link){
          /* remove v from graph*/
          status_d[v]=-3;
          break;
        }
      }
    }
  }

  *ret_lta=lta;

  //printf("end of i_select\n");



 return eslOK;

 ERROR:
   //printf("in error\n");
   return status;
}



static void
i_update_workspace(int *dec_o, int *label_o, int *status_d, int *to_add, int *assignments, size_t n, int *k, int *lta, ESL_RANDOMNESS *r){
  int i;
  int d=0;

  //printf("i update workspace called \n");

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



int
esl_iset_Cyan(void *base, size_t n, size_t size,
			  int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
			  int *workspace, int *assignments, ESL_RANDOMNESS *r)
{
  int k=n;
  int i;
  int status;
  int lta=0; /*length of to_add*/
  int *dec_o = NULL;	/* vertices in order in which decisions about them will be made */
  int *label_o = NULL; 	/* vertices in order of smallest to largest label */
  int *status_d = NULL;		/* keeps track of current status of each vertex like a dictionary; status_d[i] is the status of vertex i */
  /* -1 in iset, -3 removed from graph, >=0 still in graph- value is (most recent index of to_add checked against) + 1 */
  int *to_add= NULL;  /* vertices to add to independent set */


  dec_o = workspace;
  label_o = workspace + n;
  status_d= workspace +2*n;
  to_add= workspace +3*n;

  //printf("called esl_iset_cyan\n");
  //printf("k:\t%d\n",k);

  for (i=0; i<n; i++){
    dec_o[i]=i;
    label_o[i]=i;
    status_d[i]=0;
  }


  while (k>0){
    //printf("top of while loop\n");
    //print_array(assignments,n);
    if ((status=i_select(base, n, size, k, linkfunc, param, dec_o, label_o ,status_d, to_add, &lta))!= eslOK) goto ERROR;
    //printf("i select done... calling update workspace\n");
    i_update_workspace(dec_o, label_o ,status_d, to_add, assignments, n, &k, &lta, r);
    //printf("bottom of while loop\n");
    //printf("k:\t%d\n",k);
    //print_array(assignments,n);
  }

  //printf("finishing esl_iset_cyan\n");
  //print_array(assignments,n);
  check_iset( base, n, size, linkfunc, param, assignments);
  ERROR:
    return status;

 return eslOK;
}



/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef eslCLUSTER_TESTDRIVE
#include <math.h>

static int
test_linkage_definition(const void *v1, const void *v2, const void *param, int *ret_link)
{
  double a         = *((double *) v1); /* you have to cast a void ptr before you can dereference it */
  double b         = *((double *) v2);
  double threshold = *((double *) param);

  *ret_link =  ((fabs(a-b) <= threshold) ? TRUE : FALSE);
  return eslOK;
}

static void
utest_singlelinkage(double *testdata, int n, double threshold, int *correct_assignment, int correct_C)
{
  int   *workspace;
  int   *assignment;
  int    C;
  int    v;

  if ((workspace  = malloc(sizeof(int) * n * 2)) == NULL) esl_fatal("allocation failed");
  if ((assignment = malloc(sizeof(int) * n))     == NULL) esl_fatal("allocation failed");

  if (esl_cluster_SingleLinkage(testdata, n, sizeof(double),
				test_linkage_definition, &threshold,
				workspace, assignment, &C) != eslOK) esl_fatal("single linkage clustering failed");

  if (C != correct_C) esl_fatal("expected %d clusters, but got %d\n", correct_C, C);
  for (v = 0; v < n; v++)
    if (correct_assignment[v] != assignment[v])
      esl_fatal("expected vertex %d to be in cluster %d, but it's in %d\n", v, correct_assignment[v], assignment[v]);

  free(workspace);
  free(assignment);
}
#endif /* eslCLUSTER_TESTDRIVE */




/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef eslCLUSTER_TESTDRIVE
/* gcc -g -Wall -o test -I. -L. -DeslCLUSTER_TESTDRIVE esl_cluster.c -leasel -lm
 */
#include "esl_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_cluster.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for cluster module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  double vertex[]      = { 1.0, 2.0, 4.0, 5.0, 7.0, 8.0 };
  int    na1 = 3, a1[] = { 0,   0,   1,   1,   2,   2   };     /* correct answer when threshold = 1.5 */
  int    na2 = 6, a2[] = { 0,   1,   2,   3,   4,   5   };     /* correct answer when threshold < 1.0 */
  int    na3 = 1, a3[] = { 0,   0,   0,   0,   0,   0   };     /* correct answer when threshold > 2.0 */
  int    n         = sizeof(vertex) / sizeof(double);

  utest_singlelinkage(vertex, n, 1.5, a1, na1);
  utest_singlelinkage(vertex, n, 0.5, a2, na2);
  utest_singlelinkage(vertex, n, 2.5, a3, na3);

  esl_getopts_Destroy(go);
  return 0;
}
#endif /* eslCLUSTER_TESTDRIVE*/





/*****************************************************************
 * 4. Example
 *****************************************************************/
#ifdef eslCLUSTER_EXAMPLE
/*::cexcerpt::cluster_example::begin::*/
/* gcc -g -Wall -o example -I. -L. -DeslCLUSTER_EXAMPLE esl_cluster.c easel.c -lm  */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_cluster.h"

static int
my_linkage_definition(const void *v1, const void *v2, const void *param, int *ret_link)
{
  double a         = *((double *) v1); /* you have to cast a void ptr before you can dereference it */
  double b         = *((double *) v2);
  double threshold = *((double *) param);

  *ret_link =  ((fabs(a-b) <= threshold) ? TRUE : FALSE);
  return eslOK;
}

int
main(int argc, char **argv)
{
  double vertex[]  = { 1.0, 2.0, 4.0, 5.0, 7.0, 8.0 };
  int    n         = sizeof(vertex) / sizeof(double);
  double threshold = 1.5;
  int   *workspace;
  int   *assignment;
  int    C;
  int    v;

  workspace  = malloc(sizeof(int) * n * 2);
  assignment = malloc(sizeof(int) * n);

  esl_cluster_SingleLinkage(vertex, n, sizeof(double),
			    my_linkage_definition, &threshold,
			    workspace, assignment, &C);

  printf("There are %d clusters.\n", C);
  for (v = 0; v < n; v++) printf("vertex %d is in cluster %d\n", v, assignment[v]);

  free(workspace);
  free(assignment);
  return 0;
}
/*::cexcerpt::cluster_example::end::*/
#endif /*eslCLUSTER_EXAMPLE*/
