## esl_alloc : portable aligned memory allocation

The `esl_alloc` module provides for portable aligned allocation. This
is generally only needed for SIMD vector code.

### rationale

Yes, the C99 standard states that malloc() is _suitably aligned so
that it may be assigned to a pointer to any type of object_. But SIMD
vector types are not part of the C99 standard, vector types may be
wider than any C99 object type, and vector memory should be aligned.

Most, if not all systems we work on will provide posix_memalign().
But I don't trust it to be there; plus we have a policy of making it
as easy as possible to port to non-POSIX platforms when we can.

We use POSIX's posix_memalign(), C11's aligned_alloc, or Intel's
_mm_malloc() (in that preference order), if available on the system.
Easel configure.ac tests for them, and sets HAVE_POSIX_MEMALIGN,
HAVE_ALIGNED_ALLOC, and/or HAVE__MM_MALLOC as appropriate. If none are
available, we fall back to a handrolled implementation.

### the quasi-portable fallback implementation

The fallback implementation does unspeakable things, things that are
technically undefined-behaviour in C99, but which happen to work
everywhere I know of. Specifically, given a pointer to a malloc()
allocation, we cast that pointer to an integer (of type `uintptr_t`),
mask off its low-order bits to achieve alignment, and store those
low-order bits in a byte preceding our allocation:

```
  pp = pointer that malloc() gives us
  v
  ...X[------------------------]
     ^^
     | \
     \  p = aligned pointer we give to the caller
      \
        one byte storing r-1, where r=(p-pp), total alignment shift
        r is 1..256, so we can store r-1 in an unsigned byte.
```
       
When we free, we use the alignment shift to reconstruct what p was, so
we can call free() with the pointer that malloc() originally gave us.

#### alignment is limited to <= 256 bytes

We only use one byte for the shift r, because we don't anticipate
needing to align on more than a 256 byte boundary. Currently the
largest vectors are AVX-512's 64-byte vectors, and Intel is
projecting an AVX-1024 with 128-byte vectors. We can revisit if
needed.

#### there is an overallocation cost

In the best case, malloc() gives us an allocation that's off by
exactly one byte from a properly aligned location; r=1 and we store
0 in the byte. In the worst case, malloc() gives us a properly
aligned allocation, in which case our extra byte looks pretty
stupid, r=V, and we store V-1 in the byte.

Because the worst case behavior means we overallocate by V bytes,
for a pointer that was already properly aligned, the fallback
implementation is potentially wasteful, and to minimize the
wastage, you should minimize allocation calls where possible. For
example, it'd be better to do 2D arrays by setting pointers into a
single allocation, for example.

It would be desirable to know when the system malloc() already
returns a suitably aligned pointer, and we could then just call
malloc() directly - but I don't know a reliable way to test for that.

#### may cause unnecessary unit test failure

Currently, the unit tests deliberately compile and test the fallback
implementation, even if `esl_alloc_aligned()` is using a system call
like `posix_memalign()`. Thus it may happen that `esl_alloc_aligned()`
is working fine, but the unit test fails because `esl_alloc_aligned()`
doesn't work on some system (perhaps because of the unspeakable things
it does).

