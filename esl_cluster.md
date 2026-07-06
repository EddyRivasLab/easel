# esl_cluster - single linkage clustering

The `cluster` module implements a generalized, efficient discrete
single linkage clustering algorithm.

The clustering algorithm in `esl_cluster_SingleLinkage()` tests for
links on the fly, thus avoiding construction of an $O(N^2)$ adjacency
matrix. This gives an algorithm of $O(N)$ memory, $O(N^2)$ time
worst-case complexity for $N$ vertices. Average time typically scales
much better than this, as efficiently as $O(N)$ time in the case of a
densely connected graph that forms a single cluster.

In order to work on generalized vertices of any data type, the
implementation uses an interface akin to that of the C `qsort()`
utility: the caller provides a void pointer to an untyped array of
vertices, the number of vertices, and the size of each vertex data
element, and a function that can take untyped pointers to two vertices
and compute whether they are linked or not.

The thing to pay most attention to in the example is the mechanism of
dealing with vertices via generic untyped pointers; in particular, the
way the caller-provided linkage-determining function takes its `void
*` arguments and immediately casts them back to data types that the
caller wants to use in computing whether the two vertices are linked.

In the example, the linkage function needs only one parameter from the
caller, so a pointer to `threshold` itself is passed into the API. If
your linkage function needs more parameters, you would define a
structure that bundles them together, then pass a pointer to that
structure into `esl_cluster_SingleLinkage()`.
