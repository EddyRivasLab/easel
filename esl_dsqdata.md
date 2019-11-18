# esl_dsqdata: efficient sequence database format

The `dsqdata` module implements a binary sequence data format. It
accelerates sequence data input in four ways, compared to Easel
flatfile parsers in `sqio`:

* __Asynchronous input.__
   Disk and CPU resources are used concurrently, using POSIX threads.
   A "loader" thread does essentially nothing but read chunks of
   data.  An "unpacker" thread does CPU work to prepare loaded
   sequence data chunks for consumption. If it takes time $R$ to read
   and $P$ to process the data, instead of overall time $R+P$, with
   asynchronous input we only need time $\mathrm{max}(R,P)$.

* __Predigitization.__ 
  Sequence data in the `dsqdata` format are already encoded in
  Easel digital sequence format.  User-oriented error checking is done
  up front when the `dsqdata` file is created.
                                                                     
* __Bit packing.__ 
  Disk read time is typically rate-limiting in HMMER and other
  Easel-based programs, so minimizing data volume is critical.
  Sequence data are packed bitwise in 32-bit packets to reduce volume
  by a factor of 1.5 (protein) to 3.75 (nucleic). A packet contains
  six 5-bit residues (protein or degenerate nucleic) or fifteen 2-bit
  residues (canonical nucleic) and two control bits.
                                                                     
* __Separate metadata.__
  Sequence data and metadata (name, accession, description, taxonomy
  identifier) are stored separately in `.dsqs` and `.dsqm`
  files. This streamlines unpacking, because these data are handled
  differently.  It also allows a deferred metadata read: sequences may
  be identified simply by index number during an initial processing
  sweep, and metadata can be loaded later by random access for a small
  number of targets of interest.

The following table lists the functions in the `dsqdata` API.

| Function                       | Synopsis                                                     |
|--------------------------------|--------------------------------------------------------------|
| `esl_dsqdata_Open()`           | Open a digital sequence database for reading                 |
| `esl_dsqdata_Read()`           | Read next chunk of sequence data.                            |
| `esl_dsqdata_Recycle()`        | Give a chunk back to the reader.                             |
| `esl_dsqdata_Close()`          | Close a dsqdata reader.                                      |
| `esl_dsqdata_Write()`          | Create a dsqdata database                                    |


## dsqdata format's four files 

The format of a database `mydb` consists of four files:

| File        | Contents  | Description                                                  |
|-------------|-----------|--------------------------------------------------------------|
| `mydb`      | Stub      | Human-readable information about the data                    |
| `mydb.dsqi` | Index     | Disk offsets for each seq in metadata and sequence files     | 
| `mydb.dsqm` | Metadata  | Name, accession, description, and taxonomy ids               |
| `mydb.dsqs` | Sequence  | Sequences (digitized, packed)                                |

The database is specified on command lines by the name of the stub
file (`mydb`), without any suffix. For example,

    % myprogram mydb

says to open `mydb`. The `esl_dsqdata_Open()` call then opens all four
files.


## definition of dsqdata file formats

### the stub file

An example stub file:

```
Easel dsqdata v1 x4019752601

Original file:   refprot.fa
Original format: FASTA
Type:            amino
Sequences:       11432138
Residues:        4358716588
```

The first line is the only line in the stub file that's parsed by the
reader. Its text format matches `/Easel dsqdata v(\d+) x(\d+)/`.  The first
field is a version number for the format, $\geq$ 1. It is currently unused, but
in the future we might need it to parse different versions of the format, if we
need to update it. The second field is a 32-bit unsigned integer tag in the
range 0..$2^{32}-1$. Each of the four files carries the same randomly generated
tag. The tag is used to ascertain that the four files belong together in the
same database, as opposed to one or more of them being inadvertently clobbered
somehow by the user.

After the first line, the rest of the stub file is ignored by the Easel reader,
and can contain anything -- even your own notes, if you want to add any. The
text here is the useful information that the Easel writer writes by default.

### the .dsqi index file

The purpose of the .dsqi index file is to (quickly) give us data offsets (in
bytes) we need to randomly access the metadata and sequence data for sequence
number 0..nseq-1 in the .dsqm and .dsqs files, or for any specified range of
sequences.

The header of the binary index file consists of:

| name         | type       | description                                  |
|--------------|------------|----------------------------------------------|
| magic        | `uint32_t` | magic number (version, byte order)           |
| uniquetag    | `uint32_t` | random integer tag (0..$2^{32}-1$)           |
| alphatype    | `uint32_t` | alphabet type code (1,2,3 = RNA, DNA, amino) |
| flags        | `uint32_t` | Currently 0. Reserved for future flags       |
| max_namelen  | `uint32_t` | Maximum seq name length in metadata          |
| max_acclen   | `uint32_t` | Maximum accession length in metadata         |
| max_desclen  | `uint32_t` | Maximum description length in metadata       |
| max_seqlen   | `uint64_t` | Maximum sequence length                      |
| nseq         | `uint64_t` | Total number of sequences in database        |
| nres         | `uint64_t` | Total number of residues in database         |

The **magic** is used to check that the file is indeed a dsqdata
format file, and to detect byte order swapping. Valid values for the
magic version/byteorder number are:

| value      | derivation          | description              |
|------------|---------------------|--------------------------|
| 0xc4d3d1b1 | "dsq1" + 0x80808080 | dsqdata version 1 format |
| 0xb1d1d3c4 | above, byteswapped  | above, byteswapped       |

The **uniquetag** matches the tag seen in the other files.

The dsqdata packet format is only defined for biological sequence alphabets.
Valid integer values for the **alphatype** code come from a subset of the codes
used in `esl_alphabet.h`:

| value | `esl_alphabet.h` | description |
|-------|------------------|-------------|
| 1     | `eslRNA`         | RNA         |
| 2     | `eslDNA`         | DNA         |
| 3     | `eslAMINO`       | protein     |

The unused **flags** field gives us some flexibility for future
versions of the format.

The maximum lengths of the names, accessions, and descriptions in the
metadata file might someday be useful (in making allocations, for
example) but they are currently unused by the Easel reader.

Likewise, the maximum sequence length, total number of sequences, and
total number of residues in the database may someday be useful (for
making decisions about how to partition a parallel search, for
example), but they are currently unused too.

After the header, the remainder of the file consists of `nseq`
records of type `ESL_DSQDATA_RECORD` (defined in
`esl_dsqdata.h`): 

| element        | type      | description                                                             |
|----------------|-----------|-------------------------------------------------------------------------|
| `metadata_end` | `int64_t` | Position of terminal `\0` of metadata for seq i in .dsqm file, in bytes |
| `psq_end`      | `int64_t` | Position of final packet for sequence i in .dsqs file, in packets       |

Storing _end_ positions instead of _start_ positions allows
us to determine lengths, without needing an n+1'th sentinel record,
albeit at the cost of special casing what happens for the first
sequence i=0. For example:

```
    len[i]   = (i == 0 ? r[i].end + 1 : r[i].end - r[i-1].end)
    start[i] = (i == 0 ? 0            : r[i-1].end + 1)
```

This is equivalent to treating `r[-1].end = -1`. Some of the Easel reader's code
tracks a `last_end` variable for the end of the previous metadata or packed
sequence field i-1, which is initialized to -1. This -1 boundary condition is
why we use _signed_ `int64_t` types.

Packet sequence endpoints are stored in units of unsigned 32-bit binary
_packets_, not in bytes. To convert to a disk offset or a length in bytes you
multiply by 4 (`sizeof(uint32_t)`).

Keeping the size of the dsqdata files as small as possible is critical
because the reading speed is limited by the raw size of the
data. Therefore we don't store separate positions for the different
metadata fields (name/accession/description/taxonomy id); only one
position for all the metadata associated with sequence i. The reader
reads all of it in one chunk, and parses it for the stored `\0`
sentinels.

For the same reason, we don't store any information about _unpacked_ sequence
lengths, only the bare minimum of information that the dsqdata loader and
unpacker need to locate, load, and unpack the packed data for any given sequence
i. The unpacker determines the unpacked sequence length when it unpacks the
data.


### the .dsqm metadata file

The metadata file starts with two header fields, the same two that the
index file starts with:

| name      | type       | description                        |
|-----------|------------|------------------------------------|
| magic     | `uint32_t` | magic number (version, byte order) |
| uniquetag | `uint32_t` | random integer tag (0..$2^32-1$)   |

After the header, the remainder of the file consists of the following
data for each sequence i = 0..nseq-1:

| field        | type                        | description                                         |
|--------------|-----------------------------|-----------------------------------------------------|
| name         | `char` array ending in `\0` | sequence name (1 word, no whitespace); mandatory    |
| accession    | `char` array ending in `\0` | sequence accession (1 word, no whitespace); or "\0" |
| description  | `char` array ending in `\0` | sequence description line; or "\0"                  |
| taxonomy id  | `int32_t`                   | NCBI taxonomy identifier; or -1                     |


The name, accession, and description are variable length strings. The name and
accession are single "words" with no whitespace (`\S+`). The description is one
line, may contain spaces, but may not contain any newlines. All sequences must
have a name, so `strlen(name) > 0`. The accession and description are optional;
if they are not present, these are 0-length strings ("\0").

The taxonomy identifier is an integer in NCBI's taxonomy. Valid taxonomy
identifiers are $\geq 1$.  This field is optional; use a value of -1 to indicate
unset.

(I cannot find any documentation at NCBI on the maximum range of the taxid, nor
can I find a clear statement of whether 0 is valid or not. 0 is currently unused
in the NCBI taxonomy.  1 indicates the top level. That makes it look like it's
safe to treat 0 as "unset" but it seems even safer to go with -1 and a signed
integer. Unless NCBI ends up having more than two billion species. Currently
there are about 1.8 million.)

These names, types, and semantics match the corresponding fields in an `ESL_SQ`.

### the .dsqs sequence file

The sequence file also starts with the same two header fields that the
index and metadata files started with:

| name      | type       | description                        |
|-----------|------------|------------------------------------|
| magic     | `uint32_t` | magic number (version, byte order) |
| uniquetag | `uint32_t` | random integer tag (0..$2^32-1$)   |

After the header, the remainder of the file consists of the packed
sequences, with one packet array for each sequence i = 0..nseq-1. Each
packet array ends with a specially marked sentinel packet. The packet
format is described next.

#### packet format

Each packet is an unsigned 32 bit integer.  The two leading (most
significant) bits are control bits. Bit 31 signals EOD (end of data):
the last packet in a packed sequence. Bit 30 signals the packet
format: 1 for 5-bit, 0 for 2-bit.  The remaining bits are the packed
residue codes:

```
      [31] [30] [29..25]  [24..20]  [19..15]  [14..10]  [ 9..5 ]  [ 4..0 ]
       ^    ^   |------------  6 5-bit packed residues ------------------|
       |    |   []  []  []  []  []  []  []  []  []  []  []  []  []  []  []
       |    |   |----------- or 15 2-bit packed residues ----------------|
       |    |    
       |    "packtype" bit 30 = 0 if packet is 2-bit packed; 1 if 5-bit packed
       "sentinel" bit 31 = 1 if last packet in packed sequence; else 0
       
       (packet & (1 << 31)) tests for end of sequence
       (packet & (1 << 30)) tests for 5-bit packing vs. 2-bit
       ((packet >> shift) && 31) decodes 5-bit, for shift=25..0 in steps of 5
       ((packet >> shift) && 3)  decodes 2-bit, for shift=28..0 in steps of 2
```

Packets without the sentinel bit set are full. They unpack to 15 or 6
residues.
 
5-bit EOD packets may be partial: they unpack to 0..6 residues. The
remaining residue codes are set to 0x1f (11111), indicating EOD within
the packet. The only case in which a partial EOD packet encodes 0
residues is a zero-length sequence: there has to be at least one EOD
packet.

2-bit EOD packets must be full, because there is no way to signal EOD
locally within a 2-bit packet. It can't use 0x03 (11), because that
encodes U/T. Generally, therefore, the last packet(s) of a nucleic
acid sequence must be 5-bit encoded, solely to be able to use sentinel
residues in a partial packet, unless the end happens to come flush at
the end of a 2-bit packet. (If we ever needed to pack an alphabet of 2
or 3 residues, we could use 0x03 as a sentinel.  This seems unlikely
to ever happen, so I'm simply not going to include any code to read
EOD 2-bit partial packets.)

A protein sequence of length L packs into exactly P $= \mathrm{max}(1,
(L+5)/6)$ 5-bit packets. A DNA sequence packs into P $\leq \mathrm{max}(1,
(L+14)/15)$ mixed 2- and 5-bit packets. P $\geq 1$ because even a
zero-length sequence ($L=0$) requires an EOD packet.

A packed sequence consists of an integer number of packets, P, ending
with an EOD packet.
 
A packed amino acid sequence unpacks to $\leq$ 6P residues. All its
packets are 5-bit encoded.
 
A packed nucleic acid sequence unpacks to $\leq$ 15P residues.  The
packets are a mix of 2-bit and 5-bit. Degenerate residues must be
5-bit packed, and the EOD packet usually is too. A 5-bit packet does
not have to contain degenerate residues, because it might have been
necessary to get "in frame" to pack a downstream degenerate
residue. For example, the sequence ACGTACGTNNA... must be packed as
[ACGTAC][CGTNNA]... to get the N's packed correctly.
 

