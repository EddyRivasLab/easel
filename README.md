# Easel - a C library for biological sequence analysis

[![](https://travis-ci.org/EddyRivasLab/easel.svg?branch=develop)](https://travis-ci.org/EddyRivasLab/easel)
![](http://img.shields.io/badge/license-BSD-brightgreen.svg)

[Easel](http://bioeasel.org) is an ANSI C code library developed by
the [Eddy/Rivas laboratory](http://eddylab.org) at Harvard. Easel
supports our work on computational analysis of biological sequences
using probabilistic models. Easel is used by
[HMMER](http://hmmer.org), the profile hidden Markov model software
that underlies several protein and DNA sequence family databases such
as [Pfam](http://pfam.xfam.org), and by
[Infernal](http://eddylab.org/infernal), the profile stochastic
context-free grammar software that underlies the
[Rfam](http://rfam.xfam.org) RNA family database. Easel aims to make
similar applications more robust and easier to develop, by providing a
set of reusable, documented, and well-tested functions.

Easel is not (yet) released on its own. It is part of the HMMER and
Infernal releases.

To participate in Easel development, visit us at
[github](https://github.com/EddyRivasLab/easel).


## Building Easel source code from github:

```bash
    % git clone https://github.com/EddyRivasLab/easel
    % cd easel
    % autoconf
    % ./configure
    % make
    % make check
```

This procedure gives you our `master` branch - the most recent Easel
release that was packaged with a release of HMMER, Infernal, or other
lab software. If you want a feature that we've put in more recently,
you want our `develop` branch. After the `git clone` and `cd easel`,
do:

```bash
   % git checkout develop
```

and proceed to `autoconf`.

## Developers

The main development happens on the `develop` branch.
The build system is configured in the files `Makefile.in` and `configure.ac`.
The `autoconf` command will generate the `configure` script which is then used
to generate the `Makefile` file responsible for actually building the source
code, testing, and housekeeping of generated files:

```bash
   % autoconf
   % ./configure
```

The `./configure` command will also generate the header `esl_config.h` from
the template `esl_config.h.in`.
This header is meant to be included in every source code file to fine-tune
the compilation according to the platform being used.
The software is built by entering the usual `make` command.

**Hint:** Use `make -jN` to parallelize the compilation process, where `N` can
be ~1.5x the number of cores you have.

### Testing

After the binaries are built via `make`, one can run `make check` for testing
the software.
The test consisting in running a list of binaries called _exercises_ declared
in `testsuite/testsuite.sqc`.

During the development, it is common that the programmer wants to quickly run
a small subset of all the tests, or even to run specific tests as needed.
The command `make dev` can be used to only generate the binaries used by the
`exercises`, without actually running then.
The programmer can then simply run one of the generated `*_utest` binaries 
that are pertinent to his source code modifications.

To avoid recompiling tests that are not relevant, the programmer can override some
of the variables used in `Makefile`. For example:

```bash
   % make dev ALL_UTESTS=esl_mem_utest ALL_EXAMPLES= ALL_BENCHMARKS= ALL_EXPERIMENTS= PROGS=
```

The above command will drastically decrese the number of source files needed
to compile in order to run the `esl_mem_utest` test.
