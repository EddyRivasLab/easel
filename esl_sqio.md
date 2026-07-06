# esl_sqio - unaligned sequence file input/output

The `sqio` module handles input from unaligned sequence data
files, such as FASTA files.

`sqio` can automatically recognize and parse several different
file formats, including FASTA, UniProt, Genbank, DDBJ, and EMBL.
Additionally, it can read individual unaligned sequences from multiple
alignment files in several different formats, including Stockholm,
Clustal, aligned FASTA, PSI-BLAST, and Phylip.

Sequences can be read from normal files, directly from the
`stdin` pipe, or from `gzip`-compressed files.

Sequence files can be automatically looked for in a list of one or
more database directories, specified by an environment variable (such
as `HMMERDB`).

The module uses an `ESL_SQFILE` object which works much like an
ANSI C `FILE`, maintaining information for an open sequence file
while it's being read.


## examples

### reading sequences from a file

The `eslSQIO_EXAMPLE` example at the end of `esl_sqio.c` 
opens a file, reads sequences from it one at a time, then
closes the file.

A FASTA file named `seqfile` is opened for reading by calling
`esl_sqfile_Open(filename, format, env, &sqfp)`, which creates a new
`ESL_SQFILE` and returns it through the `sqfp` pointer. If the
`format` is passed as `eslSQFILE_UNKNOWN`, then the format of the file
is autodetected. Here, we bypass autodetection by asserting that the
file is in FASTA format by passing a `eslSQFILE_FASTA` code.  The
optional `env` argument is described below too; here, we're passing
`NULL` and not using it.

Several things can go wrong in trying to open a sequence file that are
beyond the control of Easel or your application, so it's important
that you check the return code. `esl_sqfile_Open()` returns
`eslENOTFOUND` if the file can't be opened, and
`eslEFORMAT` if the file is empty, or if autodetection can't
determine its format.

Additionally, internal errors might be thrown, which you should check
for if you installed a nonfatal error handler.

The file is then read one sequence at a time by calling
`esl_sq_Read(sqfp, sq)`. This function returns `eslOK`
if it read a new sequence, and leaves that sequence in the `sq`
object that the caller provided. When there is no more data in the
file, `esl_sq_Read()` returns `eslEOF`.

If at any point the file does not appear to be in the proper format,
`esl_sq_Read()` returns `eslEFORMAT`. The application
must check for this. The API can provide a little information about
what went wrong and where. `sqfp->filename` is the name of the
file that we were parsing (not necessarily the same as
`seqfile`; `sqfp->filename` can be a full pathname if we
used an `env` argument to look for `seqfile` in installed
database directories). The function `esl_sqfile_GetErrorBuf()`
should be called to get a pointer to the generated error message. The
buffer is a brief explanatory message that gets filled in when a
`eslEFORMAT` error occurs.

Unlike in the MSA module, you don't get access to the current line
text; some of sqio's parsers use fast block-based (`fread()`) input
instead of line-based input.

We can reuse the same `ESL_SQ` object for subsequent sequences
by calling `esl_sq_Reuse()` on it when we're done with the
previous sequence. If we wanted to load a set of sequences, we'd
`_Create()` an array of `ESL_SQ` objects.

Finally, to clean up properly, a `ESL_SQ` that was created is
destroyed with `esl_sq_Destroy(sq)`, and a `ESL_SQFILE`
is closed with `esl_sqfile_Close()`.

## digital sequence input mode

Most Easel-based programs manipulate sequences in Easel's digital
sequence format, using `alphabet`, as opposed to manipulating
them as plaintext. The `sqio` reader can be used either in text
mode or digital mode. In text mode, you get the `sq->seq` field
of the `ESL_SQ`; in digital mode, you get `sq->dsq`.

To use digital mode, both the `ESL_SQFILE` reader and the
`ESL_SQ` sequence object must be set to digital mode. The
reader, because it has an input map that maps plaintext input
characters to internal residue codes (plaintext or digital), or
errors. The sequence object, because it needs to have either its
`seq` or `dsq` field allocated. Both also carry a copy of
the pointer to the alphabet.

The `eslSQIO_EXAMPLE2` example in `esl_sqio.c` shows the standard
idiom for opening files in digital mode, autoguessing their format and
alphabet by default, and allowing format and alphabet to be specified
on the command line.


## accepted formats

Accepted unaligned sequence file formats (and their Easel format
codes) include:

| name | code | description |
|------|------|-------------|
| fasta | `eslSQFILE_FASTA` | FASTA format |
| embl | `eslSQFILE_EMBL` | EMBL DNA database format |
| genbank | `eslSQFILE_GENBANK` | GenBank DNA database format |
| ddbj | `eslSQFILE_DDBJ` | DDBJ DNA database format |
| uniprot | `eslSQFILE_UNIPROT` | UniProt protein database format |

The above names, case-insensitive, are what a user uses to specify a
format on a command line: i.e. `--informat fasta` or
`--informat FASTA`.

The codes are what you use as a developer to specify a format to an
Easel function call.

Additionally, there is a code `eslSQFILE_UNKNOWN`. It
tells `esl_sqfile_Open()` to perform format
autodetection.

There are some other formats as well, which we don't advertise
because they're less well supported.

## reading from stdin and compressed files

There are two special cases for input files.

The module can read sequence input from a stdin pipe. If the
`seqfile` argument is "-", `esl_sqfile_Open()` "opens"
standard input (really, it just associates `stdin`, which is
always open, with the `ESL_SQFILE`).

The module can read compressed sequence files. If the `seqfile`
argument to `esl_sqfile_Open()` ends in `.gz`, the file
is assumed to be compressed with `gzip`; instead of opening it
normally, `esl_sqfile_Open()` opens it as a pipe from
`gzip -dc`. Your system must support pipes to use
this. 

Obviously, the user must also have `gzip` installed and in his PATH.
Also, the OS must support the `popen()` system call (POSIX.2 compliant
operating systems do). The `configure` script automatically checks
this at compile-time and defines `HAVE_POPEN` appropriately.

For both special cases, the catch is that you can't use format
autodetection; you must provide a valid known format code when you
read from stdin or from a compressed file. Pipes are not rewindable,
and format autodetection currently relies on a two-pass algorithm: it
reads partway into the file to determine the format, then rewinds to
start parsing for real.

The `msafile` module is more advanced. Its parsers are based on the
newer `buffer` module which provides rewindable input buffers even for
stdin and pipes.

## adding a sequence parser

New parsers for new formats can be plugged into `sqio` without
any API changes. Existing Easel-based programs don't need code changes
to use the new parser.

A new parser will need a format type, a structure for parser specific
data, API functions and a hook into the `sqfile_open` function.

The list of formats are defined in `esl_sqio.h`. A new
`#define` gets added to the existing formats.

A data structure for parser specific data will need to be added
to `ESL_SQDATA`. This structure is a union of the
different parser specific data structures.

```c
typedef union {
  ESL_SQASCII_DATA ascii;
  ESL_SQNCBI_DATA  ncbi;
} ESL_SQDATA;
```

Finally, a set of parser specific function pointers need to be
defined. The functions in `esl_sqio.c` in turn call these
function pointers. The `esl_sqfile_Open` function initializes
the function pointers to NULL, so if they are not set, an exception
will occur when the function is called. At a minimum, the function
should be defined to return an `eslEUNIMPLEMENTED`. Below is a
map the the function pointers to their respective function.

| Function pointer | `sqio` function |
|------------------|-----------------|
| position | `esl_sqfile_Position` |
| close | `esl_sqfile_Close` |
| set_digital | `esl_sqfile_SetDigital` |
| guess_alphabet | `esl_sqfile_GuessAlphabet` |
| read | `esl_sqio_Read` |
| read_info | `esl_sqio_ReadInfo` |
| read_seq | `esl_sqio_ReadSequence` |
| read_window | `esl_sqio_ReadWindow` |
| echo | `esl_sqio_Echo` |
| read_block | `esl_sqio_ReadBlock` |
| open_ssi | `esl_sqfile_OpenSSI` |
| pos_by_key | `esl_sqfile_PositionByKey` |
| pos_by_number | `esl_sqfile_PositionByNumber` |
| fetch | `esl_sqio_Fetch` |
| fetch_info | `esl_sqio_FetchInfo` |
| fetch_subseq | `esl_sqio_FetchSubseq` |
| is_rewindable | `esl_sqfile_IsRewindable` |
| get_error | `esl_sqfile_GetErrorBuf` |

A hook needs to be added to the function `sqfile_open`.
This hook will try to open the specified file. If successful,
the `ESL_SQFILE` structure should be filled in with function
pointers and the parser specific data and the open hook
return `eslOK`. If the sequence files were not found for
the specific parser, an `eslENOTFOUND` is returned and
the next parser tries to open the file. Below is an example
of code that tries to open an NCBI BLAST database if not
successful, then the ASCII sequence parsers try to open the
file.

```c
    if (format == eslSQFILE_NCBI && status == eslENOTFOUND)
      status = esl_sqncbi_Open(sqfp->filename, sqfp->format, sqfp);
    if (status == eslENOTFOUND)
      status = esl_sqascii_Open(sqfp->filename, sqfp->format, sqfp);
```
