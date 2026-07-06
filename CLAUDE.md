Easel is a library used by other Eddy lab software packages including
HMMER and Infernal. 

Our coding style is described in `Codestyle.md`.

The code for the Easel library is in `esl_<foo>.[ch]` files.  Many
modules have an additional documentation file `esl_<foo>.md`,
providing some overall description of what the module is doing,
especially when there are issues that are best explained with
illustrations or well-formatted equations.

Easel includes a suite of small tools, run by the `easel` command,
that we call "miniapps": `easel reformat` and `easel seqstat`, for
example. These are in the `miniapps/` directory. Each miniapp has a
`cmd_<foo>.c` source code file and a `easel-<foo>.man.in` man page.
`miniapps/easel.c` is the `easel` program.

The tests run by `make check` are in the `testsuite/` directory.
Each module file `esl_<foo>.c` has code for compiling a "test driver"
`esl_<foo>_utest` with unit tests. The `testsuite/` directory
contains additional integrated test scripts, primarily in Python, 
including one for each miniapp (`easel-<foo>-itest.py`). `make check`
runs these with a custom harness named `sqc` (see `devkit/sqc`),
using an input list of tests in `testsuite/testsuite.sqc`.

Easel has been under continuous development for almost 40 years. We
aim to keep the code compact, well-organized, clear, and
well-documented so we can come back to it at any time, understand it,
and improve it. 
