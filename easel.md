# easel - the foundation

The easel (esl) module implements a small set of functionality shared
by all the modules: notably, the error-handling system, which is
described in [the Easel coding style](documentation/codestyle.md), as
well as some other capabilities, described below.

## secure temporary files

A program may need to write and read temporary files.  Many of the
methods for creating temporary files, even using standard library
calls, are known to create exploitable security holes
[Wheeler03, ChenDeanWagner04].

Easel provides a secure and portable POSIX procedure for obtaining an
open temporary file handle, `esl_tmpfile()`. This replaces the
ANSI C `tmpfile()` function, which is said to be insecurely
implemented on some platforms.  Because closing and reopening a
temporary file can create an exploitable race condition under certain
circumstances, `esl_tmpfile()` does not return the name of the
invisible file it creates, only an open `FILE *` handle to
it. The tmpfile is not persistent, meaning that it automatically
vanishes when the `FILE *` handle is closed. The tmpfile is
created in the usual system world-writable temporary directory, as
indicated by `TMPDIR` or `TMP` environment variables, or
`/tmp` if neither environment variable is defined.

Still, it is sometimes useful, even necessary, to close and reopen a
temporary file. For example, Easel's own test suites generate a
variety of input files for testing input parsers.  Easel also provides
the `esl_tmpfile_named()` procedure for creating a persistent
tmpfile, which returns both an open `FILE *` handle and the
name of the file. Because the tmpfile name is known, the file may be
closed and reopened.  `esl_tmpfile_named()` creates its files
relative to the current working directory, not in `TMPDIR`, in
order to reduce the chances of creating the file in a shared directory
where a race condition might be exploited. Nonetheless, secure use of
`esl_tmpfile_named()` requires that you must only reopen a
tmpfile for reading only, not for writing, and moreover, you must not
trust the contents.  (It may be possible for an attacker to replace
the tmpfile with a symlink to another file.)

Example code at the end of `easel.c` shows how to create tmpfiles.

### input maps

An _input map_ is for converting input ASCII symbols to
internal encodings. It is a many-to-one mapping of the 128 7-bit ASCII
symbol codes (0..127) onto new ASCII symbol codes. It is defined as
an `unsigned char inmap[128]` or a `unsigned char *`
allocated for 128 entries.

Input maps are used in two contexts: for filtering ASCII text input
into internal text strings, and for converting ASCII input or internal
ASCII strings into internal digitized sequences (an `alphabet`
object contains an input map that it uses for digitization).

The rationale for input maps is the following. The ASCII strings that
represent biosequence data require frequent massaging. An input file
might have sequence data mixed up with numerical coordinates and
punctuation for human readability. We might want to distinguish
characters that represent residues (that should be input) from
characters for coordinates and punctuation (that should be ignored)
from characters that aren't supposed to be present at all (that should
trigger an error or warning). Also, in representing a sequence string
internally, we might want to map the symbols in an input string onto a
smaller internal alphabet. For example, we might want to be
case-insensitive (allow both T and t to represent thymine), or we
might want to allow an input T to mean U in a program that deals with
RNA sequence analysis, so that input files can either contain RNA or
DNA sequence data.  Easel reuses the input map concept in routines
involved in reading and representing input character sequences, for
example in the `alphabet`, `sqio`, and `msa`
modules.


## References

* **[ChenDeanWagner04]** H Chen, D Dean, D Wagner, "Model Checking One
  Million Lines of C Code," _Network and Distributed System Security
  Symposium_, pp. 171-185 (2004).

* **[Wheeler03]** DA Wheeler, [_Secure Programming for Linux and Unix
  HOWTO_](http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/index.html)
  (2003).
