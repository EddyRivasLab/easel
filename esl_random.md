# esl_random - pseudorandom numbers and sampling

The `random` module contains routines for generating uniformly
distributed pseudorandom numbers and sampling random deviates from
distributions. The heart of the module is the `esl_random()`
pseudorandom number generator.

The `esl_random()` random number generator is portable, reentrant, and
threadsafe. It gives reproducible results on all platforms.

The default `esl_random()` generator implements the Mersenne Twister
algorithm MT19937 [Matsumoto98]. MT19937 has strong properties,
including a period of $2^{19937}-1$ and equidistribution over $2^{32}$
values. The default Mersenne Twister should be suitable for all but a
few speed-critical applications.

Alternatively, a simple and classic linear congruential generator
(LCG) can be chosen [Knu-81a]. The LCG is much faster to initialize
(about 20x) and somewhat faster to generate samples (about 25%), while
still generating pseudorandom numbers suitable for most
applications. Because of its initialization speed, the LCG is
advantageous when a small number of reasonably random samples is needed
in a speed-critical application. However, it has a relatively short
period ($2^{32}$), making it unsuitable for large simulations.

A new generator can either be seeded with a number that you provide, or
using an arbitrary (quasirandom) seed. If you seed it yourself, you can
guarantee reproducibility: two generators seeded with the same seed
will give exactly the same sequence, even on different hardware
platforms and operating systems. If you let it seed itself arbitrarily,
you will get different sequences. Multiple instances of running the
same program will get quite different arbitrary seeds even if you start
them at the same time.

Bit streams from `esl_random()`'s two generators were tested against a
National Institute of Standards and Technology statistical benchmark
for random number generators [NIST08]. The default Mersenne Twister
passes the benchmark suite as expected [Matsumoto98]. The fast LCG
passes most but not all of the NIST tests.

`esl_random()` returns a double-precision floating point sample on the
interval $0.0 \leq x < 1$.

The module implements one object, `ESL_RANDOMNESS`, which contains
state information for the random number generator. This makes random
number generation reentrant and threadsafe. You can have more than one
active generator and they will not interfere with each other. The
object is meant to be opaque; you should not need to use its contents.

## example

The `eslRANDOM_EXAMPLE` example at the end of `esl_random.c`
initializes the random number generator with a seed provided on the
command line, then samples 10 random numbers using `esl_random()`.

When a `ESL_RANDOMNESS` object is created with
`esl_randomness_Create()`, it needs to be given a _seed_, an integer
$\geq 0$, which specifies the initial state of the generator. After a
generator is seeded, it is typically never seeded again. A series of
`esl_random()` calls generates a pseudorandom number sequence from that
starting point. If you create two `ESL_RANDOMNESS` objects seeded
identically, they are guaranteed to generate the same random number
sequence on all platforms. This makes it possible to reproduce
stochastic simulations. Thus, if you run the example multiple times,
you get the same ten numbers, because the generator is always seeded
with 42.

Often one wants different runs to generate different random number
sequences, which creates a chicken and the egg problem: how can we
select a pseudorandom seed for the pseudorandom number generator?
Calling `esl_randomness_Create(0)` (i.e., a seed argument of 0) causes
Easel to select an arbitrary seed. The arbitrary seed is constructed by
a combination of the current wall clock time (in seconds), the elapsed
cpu time since starting the program (in milliseconds or microseconds),
and (if available) the process id. Two different `ESL_RANDOMNESS`
objects created this way are expected to always produce different
pseudorandom number sequences.

Specifically, the arbitrary seed is constructed by a bitwise mixing
function that combines input from `time()`, `clock()`, and
`getpid()`. On some platforms, `getpid()` is not available, and an
arbitrary constant is used instead; on those platforms, arbitrary
seeds are a little less arbitrary, but are still quite randomly
distributed. It is improbable to get two generators with the same
arbitrary seed; to try, you would have to start two generators in the
same process at the same time.

## References

* **[Knu-81a]** DE Knuth, _The Art of Computer Programming: Seminumerical Algorithms_, 2nd edition, Addison-Wesley (1981).
* **[Matsumoto98]** M Matsumoto, T Nishimura, "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator," ACM Transactions on Modeling and Computer Simulation 8:3-30 (1998).
* **[NIST08]** A Rukhin, J Soto, J Nechvatal, M Smid, E Barker, S Leigh, M Levenson, M Vangel, D Banks, A Heckert, J Dray, S Vo, LE Bassham III, ["A Statistical Test Suite for Random and Pseudorandom Number Generators for Cryptographic Applications"](http://csrc.nist.gov/publications/nistpubs/800-22-rev1/SP800-22rev1.pdf), National Institute of Standards and Technology (2008). Special Publication 800-22, Revision 1.
