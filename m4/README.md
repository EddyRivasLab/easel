
# autoconf macros

Autoconf macros used by the Easel configure script.

## build features (compiler, make)

__ax_check_gnu_make.m4__: 


__ax_compiler_vendor.m4__:

__ax_check_compile_flag.m4__:

__ax_gcc_func_attribute.m4__:


## vector instruction set support

| macro            | shell variables set                    |
|------------------|----------------------------------------|
| __esl_sse.m4__   | `esl_have_sse`, `esl_sse_cflags`       |
| __esl_avx.m4__   | `esl_have_avx`, `esl_avx_cflags`       |
| __esl_avx512.m4__| `esl_have_avx512`, `esl_avx512_cflags` |
| __esl_neon.m4__  | `esl_have_neon`,   `esl_neon_cflags`   |



## parallelization support

__ax_mpi.m4__:

__ax_pthread.m4__:











