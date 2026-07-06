# esl_ssi - indexing flatfiles for record retrieval

The `ssi` module is for creating and using "SSI"
(sequence/subsequence index) files. SSI indexes flatfile databases by
names and/or accessions, enabling fast sequence record retrieval.

An SSI index is a binary file that stores sequence names or accessions
as _keys_, associating them with information about the sequence
record, including its location (file name and disk offset), so that it
can be looked up rapidly. It differentiates between _primary
keys_ and _secondary keys_ (aliases). There is one and only one
primary key per record. There can be more than one secondary key
(alias) per sequence. Both primary and secondary keys must be unique
identifiers (no two records have the same key). For example, a program
for sequence retrieval might create an SSI index with accessions as
primary keys and names as secondary keys (or the other way around).

Records can also be retrieved by number from the list of primary keys.
This may be useful for distributed data-parallel applications, which
can use SSI to rapidly position individual processes at different
record ranges in a flatfile database.

A single SSI file can index a sequence database that consists of more
than one individual sequence file. For example, the GenBank database
is distributed as a large number of flatfiles; one SSI file can index
them all.

When records are consistently formatted, SSI indices can allow a
specific subsequence in a sequence record to be identified
rapidly. This is useful when the sequence records are very large, such
as whole assembled genomes or chromosomes.

Although SSI indices are designed with sequence databases in mind, SSI
can also be used to index records in other types of flatfile
databases. For example, HMMER uses SSI to index HMM databases, and the
`msa` module can use SSI to index Stockholm format multiple
alignment databases like Pfam and Rfam.

SSI can handle large amounts of data. It is capable of indexing tens
of thousands of files and hundreds of trillions of sequence records.
The lengths of primary keys, secondary keys, or filenames are
effectively unlimited, and individual sequence records may be many
trillions of residues long, orders of magnitude larger than the
largest complete chromosomes. SSI indexing is effectively limited only
by the size of the SSI index file itself (up to 8 exabytes, on 64-bit
filesystems).

Binary SSI indices are portable between different machines. The sole
exception is that SSI indices built for 64-bit filesystems might not
be readable on a fully 32-bit filesystem.

A `ESL_SSI` object is used for reading an index, and a
`ESL_NEWSSI` object is used for creating one. There is also a
set of functions for portable binary file i/o.


## examples

### creating an SSI index

The `eslSSI_EXAMPLE` example at the end of `esl_ssi.c` 
creates an SSI index
for a FASTA sequence file, in which sequence records start with a line
like:

```c
 >SEQ_NAME  Rest of the line is a free-text description.
```

* A new index is created (`esl_newssi_Open()`). You can
  name an SSI index anything you want, but for an index of a
  single file, Easel generally defaults to assuming an
  `.ssi` suffix is appended to the filename. That's what
  the example does here.

* Each file to be indexed is added to the index by a call to
  `esl_newssi_AddFile()`. This returns a _file handle_
  (`fh`) that you will need when you add primary keys. In
  this example, there is only one file and only one file handle.

* You need to determine the disk offset at the exact beginning of
  each sequence record. You retrieve your current position in the
  file using an `ftello()` call.

* You add each primary key to the index with a
  `esl_newssi_AddKey()` call. You provide the handle of the
  file that key is in, and the offset to the start of this key's
  sequence record.

* The `esl_fgets()` function (part of the `easel`
  foundation module) is a way of reading text files line by line,
  no matter how long each line might be: `esl_fgets()`
  reallocates its buffer as needed.

* `esl_newssi_Write()` saves the index to an open file.

* Finally, the index structure is freed by
  `esl_newssi_Close()`.



### using an SSI index

The `eslSSI_EXAMPLE2` example in `esl_ssi.c` retrieves a FASTA
sequence record by its name, using an SSI index.

* `esl_ssi_Open()` opens the SSI index file.

* `esl_ssi_FindName()` looks up the record by its name.
  Primary keys are checked first, then secondary keys. If it is
  found, `fh` contains a file handle (what file it's in),
  and `offset` contains the position of the desired record
  in that file.

* The file handle `fh` is looked up in the file index with
  `esl_ssi_FileInfo()`, and the name of the file and a
  format code are returned. The format code is useful if you need
  to hand the filename off to different kinds of file parsers,
  depending on what file type it is. (SSI can index files in
  heterogeneous formats.)

* After that, you use the retrieved information however you need,
  independent of the SSI index. The example emphasizes this, by
  freeing the SSI index immediately with `esl_ssi_Close()`
  after it knows `fafile` and `offset`. The example
  opens the file, positions the disk with `fseeko()`, and
  reads a sequence record out of it one line at a time, until it
  reaches EOF or the start of the next sequence record.

## SSI file format

There are four sections to the SSI file:

* **Header** — Contains a magic number indicating SSI version number, followed by
  information about the number and sizes of items in the index.

* **Files** — Contains one or more _file records_, one per sequence file that's
  indexed. These contain information about the individual files.

* **Primary keys** — Contains one or more _primary key records_, one per primary key.

* **Secondary keys** — Contains one or more _secondary key records_, one per secondary key.

All numeric quantities are stored as fixed-width unsigned integers in
network (bigendian) order, for crossplatform portability of the index
files, using types `uint16_t`, `uint32_t`, and
`uint64_t`.[^2] Values may need to be cast to
signed quantities elsewhere in Easel, so only half of their dynamic
range is valid (e.g. 0..32,767 for values of type `uint16_t`;
0..2,146,483,647 (2 billion) for `uint32_t`; and 0..9.22e18 (9
million trillion) for `uint64_t`).

[^2]: These types are available on C99-compliant systems. On other
systems, Easel automatically defines appropriate substitutes at
configuration time.

File offsets (type `off_t`) are assumed to be either 32-bit or
64-bit signed integers. Easel uses 64-bit offsets if at all possible
on your system. The size of `off_t` for the system that created
the SSI file is stored in the SSI header, for portability to other
systems that try to read the SSI file.

### header section

The header section contains:

| Variable | Description | Bytes | Type |
|----------|-------------|-------|------|
| `magic` | SSI version magic number. | 4 | `uint32_t` |
| `flags` | Optional behavior flags (see below) | 4 | `uint32_t` |
| `offsz` | `off_t` size on system that made index | 4 | `uint32_t` |
| `nfiles` | Number of files in file section. | 2 | `uint16_t` |
| `nprimary` | Number of primary keys indexed. | 8 | `uint64_t` |
| `nsecondary` | Number of secondary keys indexed. | 8 | `uint64_t` |
| `flen` | Length of filenames (incl. '`\0`') | 4 | `uint32_t` |
| `plen` | Length of primary key names (incl. '`\0`') | 4 | `uint32_t` |
| `slen` | Length of sec. key names (incl. '`\0`') | 4 | `uint32_t` |
| `frecsize` | # of bytes in a file record | 4 | `uint32_t` |
| `precsize` | # of bytes in a primary key record | 4 | `uint32_t` |
| `srecsize` | # of bytes in a sec. key record | 4 | `uint32_t` |
| `foffset` | disk offset, start of file records | † | `off_t` |
| `poffset` | disk offset, start of primary key recs | † | `off_t` |
| `soffset` | disk offset, start of sec. key records | † | `off_t` |

The `flags` field is currently unused. It is stored for possible
future use, for any optional behaviors that may need to be
implemented.

The reason to explicitly record various record sizes
(`frecsize`, `precsize`, `srecsize`) and index file
positions (`foffset`, `poffset`, `soffset`) is to
allow for future extensions. Using explicit offsets means we can add
more fields in future versions of SSI without breaking older SSI
parsers. The format is meant to be both forwards- and
backwards-compatible.

### file section

The file section consists of `nfiles` file records. Each record
is `frecsize` bytes long, and contains:

| Variable | Description | Bytes | Type |
|----------|-------------|-------|------|
| `filename` | Name of file (possibly including full path) | `flen` | `char *` |
| `format` | Format code for file | 4 | `uint32_t` |
| `flags` | Optional behavior flags | 4 | `uint32_t` |
| `bpl` | Bytes per sequence data line | 4 | `uint32_t` |
| `rpl` | Residues per sequence data line | 4 | `uint32_t` |

When a SSI file is written, `frecsize` is equal to the sum of
the sizes above. When a SSI file is read by a parser, it is possible
that `frecsize` is larger than the parser expects, if the parser
is expecting an older version of the SSI format, because additional
fields might be present. The parser will only try to parse data up to
the `frecsize` it expected to see, but still knows the (possibly
larger) `frecsize` that is operative in this SSI file, for
purposes of skipping around in the index file.

An SSI index might reside in the same directory as the data file(s) it
indexes, so `filename` might be relative to the location of the
SSI index. Alternatively, `filename` might be a full path. These
semantics are not enforced by the `ssi` module. Rather, this is
an issue for an SSI-enabled application to define for
itself. SSI-enabled applications would typically include program(s)
for creating indices and program(s) for using them. Different
applications might employ different conventions for where the indices
are expected to be, relative to the sequence files, so long as that
convention is consistently applied by both index creator and index
user.

Similarly, the `ssi` module does not specify the meaning of the
`format` code. An SSI-enabled application may use this field to
associate any useful format code (or indeed, any other number) with
each indexed file. A typical use, though, would be sequence file
format codes like `eslSQFILE_FASTA` or
`eslMSAFILE_STOCKHOLM` from the `sqio` or `msa`
modules.

Only one possible optional behavior flag is currently defined:

| Flag | Value | Note |
|------|-------|------|
| `eslSSI_FASTSUBSEQ` | $1 \ll 0$ | Fast subseq retrieval is possible for this file. |

When `eslSSI_FASTSUBSEQ` is set, `bpl` and `rpl`
are nonzero. These can then be used to calculate the offset of
subsequence positions in the data file. This optional behavior is
described in detail a bit later.

### primary key section

The primary key section consists of `nprimary` records. Each
record is `precsize` bytes long, and contains:

| Variable | Description | Bytes | Type |
|----------|-------------|-------|------|
| `key` | Key name (seq name, identifier, accession) | `plen` | `char *` |
| `fnum` | File number (0..nfiles-1) | 2 | `uint16_t` |
| `r_off` | Offset to start of record | ‡ | `off_t` |
| `d_off` | Offset to start of sequence data | ‡ | `off_t` |
| `len` | Length of data (e.g. seq length, residues) | 8 | `uint64_t` |

The two offsets are sequence file offsets that may be either 8 or 4
bytes (indicated by ‡ above). They are usually 64-bit (8 byte)
signed integers. If an SSI index is created on an older system that
only allows 32-bit offsets (and hence cannot have files >2 GB), they
are 32-bit (4-byte) signed integers.

`r_off` (the _record offset_) indicates the position of
the first byte of the record.

`d_off` (the _data offset_) is optional. It indicated the
position of the first byte of the data in the record (the sequence
itself, for example), after any header information. If
`eslSSI_FASTSUBSEQ` is set on this key's file, `d_off`
can be used to calculate a disk position close to (and possibly
exactly at) the start of any subsequence.

`len` is the length of the data record. It is optional, because
some kinds of records that SSI might be used to index may not have a
meaningful length. The units of length are application-defined (i.e.
defined by whatever creates the SSI index for a particular file); but
for sequences, `len` is almost certainly in residues. In
subsequence retrieval, a `len` in residues is necessary for
bounds checking.

### secondary key section

The secondary key section consists of `nsecondary` records. Each
record is `srecsize` bytes long, and contains:

| Variable | Description | Bytes | Type |
|----------|-------------|-------|------|
| `key` | Key name (seq name, identifier, accession) | `slen` | `char *` |
| `pkey` | Primary key | `plen` | `char *` |

That is, secondary keys are simply associated with primary keys as
_aliases_. There can be many secondary keys for a given record.
However, all keys (primary and secondary) must be unique: no key can
occur more than once in the index.

## fast subsequence retrieval

In some files (notably whole chromosomal DNA sequences) the size of
each sequence is large. It may be slow (even prohibitively slow) to
extract a desired subsequence, even if an SSI index says how to find
the sequence record quickly, if you had to read the entire sequence
into memory just to extract the right part of it.

SSI uses a simple but effective technique to find subsequences.
Provided that he sequence data file is consistently formatted so that
each line in each record (except the last one) is of the same length,
in both bytes and residues, we can determine a disk offset of the
start of any subsequence by arithmetic. Easel refers to such a file
as "well-formatted". For example, a simple well-formatted FASTA file
with 50 residues per line might have 51 bytes on every sequence line
(counting the '`\0`') but for the last line in each record
(`bpl`=51, `rpl`=50). Position $i$ in a sequence $1..L$
will be on line $l = (i-1)/\mathrm{rpl}$, and line $l$ starts at
disk offset $l * \mathrm{bpl}$ relative to the start of the
sequence data.

If there are no nonsequence characters in the data line except the
terminal '`\0`' (which is true iff `bpl` = `rpl`+1
and 1 residue = 1 byte), we can precisely identify the disk position
of any residue $i$ (_single residue resolution_):

$$
\mathrm{relative\ offset\ of\ residue\ } i =
\left((i-1)/\mathrm{rpl}\right)*\mathrm{bpl} + (i-1) \bmod \mathrm{rpl}
$$

Even for sequence data lines with extra characters (e.g. spaces,
coordinates, whatever), we can still identify the start of the text
line that residue $i$ is on (_line resolution_). A parser can be
positioned at the beginning of the appropriate line $l$, which starts
at residue $(l*\mathrm{rpl}) + 1$, and it can start reading from
there (e.g. the line that $i$ is on) rather than the beginning of the
whole sequence record.

When creating an index, your application is responsible for
determining if `bpl` and `rpl` are consistent throughout a
file. If so, you may call `esl_newssi_SetSubseq()` on that
file's handle to set `bpl`, `rpl`, and the
`eslSSI_FASTSUBSEQ` flag. Then, when using that index, you can
use the `esl_ssi_FindSubseq()` call to retrieve not only the
filehandle `fh` and record offset `r_off` for a key; you
also provide a desired start position `requested_start` for the
subsequence you want to retrieve, and the routine gives you back a
data offset `d_off`, which corresponds to a actual starting
position `actual_start` that is also returned. For single
residue resolution, `actual_start` is `requested_start`,
and the data offset `d_off` will position you right at the
residue you want; you position the file with `fseeko()` and
start reading your subsequence immediately. When we can only achieve
line resolution, `actual_start` is $\leq$
`requested_start`; you position the disk to the start of the
appropriate line with `fseeko()`, start reading, and skip zero
or more residues to reach your `requested_start`. Your
application should be prepared to deal with line resolution; it should
not assume that `requested_start` and `actual_start` are
identical.

Data is always read "left to right". To read a reverse complemented
strand in DNA files, you must read your subsequence in forward
orientation first, and reverse complement it later.
