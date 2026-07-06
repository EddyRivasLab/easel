# esl_msa - Multiple sequence alignments and i/o

The `msa` module reads and writes multiple sequence alignment
files.

The module uses two objects. An `ESL_MSA` holds a multiple
sequence alignment. An `ESL_MSAFILE` is an alignment file,
opened for input. No object is needed for output of an alignment file;
a normal C `FILE` stream is used for output.

MSAs can be handled in "text mode" or "digital mode", and
converted back and forth between modes.

Large MSA database files like Pfam or Rfam can be indexed with SSI,
allowing fast random access. The `esl_msafile_Open()` and
`esl_msafile_OpenDigital()` functions automatically open an
accompanying SSI index, if it is present.

Some things about the use of the API in the example are worth noting:

1. The format of the alignment file can either be automatically
   detected, or set by the caller when the file is opened.
   Autodetection is invoked when the caller passes a format code
   (here, `fmt`) of
   `eslMSAFILE_UNKNOWN`. Autodetection is a "best effort"
   guess, but it is not 100% reliable - especially if the input
   file isn't an alignment file at all. So autodetection is a
   convenient default, but the caller will probably want to provide
   a way for the user to specify the input file format and override
   autodetection, just in case.

2. Errors can occur either in opening or reading the file that you
   must check for. This error checking could be as simple as making
   sure that `esl_msafile_Open()` and
   `esl_msa_Read()` returned `eslOK`, but the example
   shows how to catch all the normal errors returned by these
   calls, and how to format some reasonably informative error
   messages for the user. For example, when parsing the file fails
   and `esl_msa_Read()` returns an `eslEFORMAT`
   error, information about the problem is stored in `afp`:
   the caller can use `afp->linenumber`, `afp->buf`,
   and `afp->errbuf` to get the line number in the file that
   the error occurred, the text that was on that line, and a short
   error message about what was wrong with it, respectively.

3. To output (write) an alignment, open a normal C `FILE`
   stream, write the alignment(s) with `esl_msa_Write()`,
   and close the stream with C's `fclose()`. Here, the
   example is regurgitating the alignments it reads to
   `stdout`.


## accessing alignment data

The information in the `ESL_MSA` object is meant to be accessed
directly, so you need to know what it contains. This object is defined
and documented in `esl_msa.h`. It contains various information,
as follows:

### important/mandatory information

The following information is always available in an MSA (except
digital-mode alignments, which replace `aseq[][]` with
`ax[][]`, as described later):

```c
  char  **aseq;       /* alignment itself, [0..nseq-1][0..alen-1], \0-terminated */
  char  **sqname;     /* sequence names [0..nseq-1][], \0-terminated             */
  double *wgt;        /* sequence weights [0..nseq-1], default 1.0               */
  int64_t alen;       /* length of alignment (columns); or (if growable) -1      */
  int     nseq;       /* number of seqs in alignment; or (if growable) blocksize */
  int     flags;      /* flags for what info has been set                        */
```

The alignment contains `nseq` sequences, each of which contains
`alen` characters.

`aseq[i]` is the i'th aligned sequence, numbered
`0..nseq-1`.

`aseq[i][j]` is the j'th character in aligned sequence i,
numbered `0..alen-1`.

`sqname[i]` is the name of the i'th sequence.

`wgt[i]` is a non-negative real-valued weight for sequence
i. This defaults to 1.0 if the alignment file did not provide weight
data. You can determine whether weight data was parsed by checking
`(flags & eslMSA_HASWGTS)`.

### optional information

The following information is optional. It is usually only provided by
annotated Stockholm alignments (for instance, Pfam and Rfam database
alignments):

```c
  char  *name;      /* name of alignment, or NULL                                           */
  char  *desc;      /* description of alignment, or NULL                                    */
  char  *acc;       /* accession of alignment, or NULL                                      */
  char  *au;        /* "author" information, or NULL                                        */
  char  *ss_cons;   /* consensus sec structure, or NULL;  [0..alen-1], even in digital mode */
  char  *sa_cons;   /* consensus surface access, or NULL; [0..alen-1], even in digital mode */
  char  *pp_cons;   /* consensus posterior prob, or NULL; [0..alen-1], even in digital mode */
  char  *rf;        /* reference coord system, or NULL;   [0..alen-1], even in digital mode */
  char  *mm;        /* model mask, or NULL;   [0..alen-1], even in digital mode             */
  char **sqacc;     /* accession numbers for sequences i                                    */
  char **sqdesc;    /* description lines for sequences i                                    */
  char **ss;        /* per-seq secondary structures, or NULL    (string, \0-term)           */
  char **sa;        /* per-seq surface accessibilities, or NULL (string, \0-term)           */
  char **pp;        /* posterior prob per residue, or NULL.     (string, \0-term)           */
  float  cutoff[eslMSA_NCUTS];  /* NC/TC/GA cutoffs propagated to Pfam/Rfam                 */
  int    cutset[eslMSA_NCUTS];  /* TRUE if a cutoff is set; else FALSE                      */
```

These should be self-explanatory; but for more information, see the
Stockholm format documentation. Each of these fields corresponds to
Stockholm markup.

These pointers will be NULL for any optional annotation that was not
present in the alignment file. This is true at any level; for
instance, `ss` will be NULL if no secondary structures are
available for any sequence, and `ss[i]` will be NULL if some
secondary structures are available, but not for sequence i.

The `cutoff` array contains Pfam/Rfam curated trusted, gathering
and noise score cutoffs. They are indexed as follows:

```c
#define eslMSA_TC1     0
#define eslMSA_TC2     1
#define eslMSA_GA1     2
#define eslMSA_GA2     3
#define eslMSA_NC1     4
#define eslMSA_NC2     5
#define eslMSA_NCUTS   6
```

### unparsed information

The MSA object also stores additional "unparsed" information from
Stockholm files; that is, tags that are present but not recognized by
the MSA module. This information is stored so that it may be
regurgitated if the application needs to faithfully output the entire
alignment file, even the bits that it didn't understand. If you need
to access unparsed Stockholm tags, see the comments in
`esl_msa.h`.

### off-by-one issues in indexing alignment columns

With one exception, all arrays over alignment columns are normal C
string arrays, indexed `0..alen-1`. This includes optional
information such as `msa->rf[]` (the reference annotation line)
and `msa->cs[]` (the consensus structure annotation line).

The exception is a digitized sequence alignment, `msa->ax[][]`
(see below), where columns are indexed 1..alen and sentinel bytes at
positions 0 and alen+1, following Easel's convention for digitized
sequences.

Thus, when your code is manipulating a digitized alignment and using
optional information like the reference annotation line or the
consensus structure line, you must be careful of the off-by-one
difference in how the two types of data are indexed.


## accepted formats

Currently, the MSA module only parses Stockholm format.

Stockholm format and other alignment formats are documented in a later
chapter.


## digital versus text representation

The multiple alignment can be stored either in text or digital mode.

A text-mode MSA stores ASCII text symbols in a 2D array `char ** aseq[0..nseq-1][0..alen-1]`. These strings are stored exactly as
they appeared in the original file; they aren't converted to upper or
lower case, for example.

A digital-mode MSA is digitized in the Easel internal alphabet. This
enables more consistent, robust, and speedy handling of the sequence
data.

Text mode is the default behavior. An `ESL_MSA` is in digital
mode if its `eslMSA_DIGITAL` flag is up (`msa->flags & eslMSA_DIGITAL` is `TRUE`). When the alignment data are in
digital mode, they are stored internally as a 2D digital sequence
array `ESL_DSQ ** ax[0..nseq-1][1..alen]`, and the `aseq`
field is `NULL`.

To use a digital internal representation, it is most efficient to read
directly as digital data, using a `esl_msafile_OpenDigital()`
call in place of `esl_msafile_Open()`. You can also change the
mode of an MSA from text to digital using
`esl_msa_Digitize()`, and digital to text using
`esl_msa_Textize()`.

Suppose you want to open an alignment file and read its alignments in
digital mode, but you don't know whether the file contains DNA or
protein alignments. You can't use `esl_msafile_OpenDigital()`
unless you have an alphabet; but you can't see the alphabet until
you've read an alignment. Easel provides
`esl_msafile_GuessAlphabet()` to peek at the first alignment
and guess its alphabet [^1], and `esl_msafile_SetDigital()` to set an
already-open `ESL_MSAFILE` so that all subsequent alignments
are read in digital mode.

[^1]: Because the stream that alignments are
being read from may be non-rewindable, the implementation of
`esl_msafile_GuessAlphabet()` reads and caches the first
alignment.


## reading from stdin or gzip-compressed files

The module can read compressed alignment files. If the
`filename` passed to `esl_msafile_Open()` ends in
`.gz`, the file is assumed to be compressed with gzip. Instead
of opening it normally, `esl_msafile_Open()` opens it as a pipe
through `gzip -dc`. Obviously this only works on a POSIX
system -- pipes have to work, specifically the `popen()` system
call -- and `gzip` must be installed and in the PATH.

The module can also read from a standard input pipe. If the
`filename` passed to `esl_msafile_Open()` is `-`,
the alignment is read from `STDIN` rather than from a file.

Because of the way format autodetection works, you cannot use it when
reading from a pipe or compressed file. The application must know the
appropriate format and pass that code it calls
`esl_msafile_Open()`.
