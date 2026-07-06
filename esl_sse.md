# esl_sse - SIMD minivectors on Intel/AMD

The `sse` module provides a few vectorized functions that use
the Intel/AMD SSE (Streaming SIMD Extensions) Intrinsics: most
importantly, vectorized `logf()` and `expf()` routines.

The `sse` module is only available on platforms that support
SSE, SSE2, and SSE4.1 instructions. This includes all modern Intel and
AMD processors (since Intel Penryn in 2007, and AMD Bulldozer, Jaguar,
and Piledriver processors since 2011), but not PowerPC processors. By
default, the Easel configure script enables SSE if it is available on
the compilation machine.

The `eslSSE_EXAMPLE` example at the end of `esl_sse.c` calculates
`logf()` and `expf()` on an SSE `__m128` vector
containing four floats. It also shows a useful `union` idiom for
accessing four floats either as an SSE vector or as individual floats.
