# esl_vectorops

The `vectorops` module contains routines for simple operations on
vectors of various types, of length <n>. 

Different functions allow an operation to be performed in vectors
containing elements of different scalar types (`double`, `float`,
`int`, `int64_t`, `char`). The appropriate routine is prefixed with
`D`, `F`, `I`, `L`, or `C`. For example, `esl_vec_DSet()` is the Set
routine for a vector of doubles; `esl_vec_LSet()` is for 64-bit
("Long") integers.

Floating point `*Sum()` functions use the Kahan compensated summation
algorithm to reduce numerical error [[Kahan,
1965]](https://doi.org/10.1145/363707.363723).
	
Floating point `*Compare()` functions take an additional argument
`tol`, and use `esl_*Compare(v[i], w[i], tol)` to compare each element
with **relative** tolerance `tol`. We are phasing `esl_*Compare()`
functions out. They are to be replaced by `esl_*CompareNew()` which
use `atol` and `rtol` absolute and relative tolerances.

## notes to future 

* We could provide SIMD accelerated versions and runtime dispatchers. 

* The vector length <n> is always an `int`, so vectors can be up to a
  couple billion ($2^{31}-1$) long. If we ever need longer vectors
  (!), we'll write a `esl_vec64` module. That is, the contents (type)
  of the vector is a different issue from its length.  Don't be
  confused to see the `L` versions of the functions using `int64_t
  *vec` and `int n`, and using `int i` to index it; this is
  deliberate.



