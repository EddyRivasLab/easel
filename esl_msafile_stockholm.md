# Stockholm format for multiple sequence alignments

Stockholm format was developed by the Pfam Consortium to support
extensible markup and annotation of multiple sequence alignments.

A principal advantage of Stockholm is the ability to add annotation
and metadata to the file, to sequences, to alignment columns, or to
individual residues.

## A minimal Stockholm file

```
# STOCKHOLM 1.0

seq1  ACDEF...GHIKL
seq2  ACDEF...GHIKL
seq3  ...EFMNRGHIKL

seq1  MNPQTVWY
seq2  MNPQTVWY
seq3  MNPQT...
//
```

The first line in the file must be `# STOCKHOLM x.y`, where `x.y` is a
major/minor version number for the format specification. This line
allows a parser to identify the file format. There is currently only
one version of Stockholm format, `1.0`.

In the alignment, each line contains a name followed by the aligned
sequence. Neither the name nor the aligned sequence may contain
whitespace characters. Stockholm does not enforce any other character
conventions on the name or the aligned sequence. Gaps are indicated in
an aligned sequence by a dash, underscore, or period.

If the alignment is too long to fit on one line, the alignment may be
split into multiple blocks, with blocks separated by blank lines. The
number of sequences, their order, and their names must be the same in
every block. Within a given block, each (sub)sequence (and any
associated `#=GR` and `#=GC` markup, see below) is of equal length,
called the _block length_. Block lengths may differ from block to
block; the block length must be at least one residue, and there is no
maximum.

Stockholm files in a single-block format (no interleaving) can be
specifically referred to as "Pfam" format. That is, if you specify
`--informat stockholm`, the MSA can be single-block or multi-block; if
`--informat pfam`, it must be single-block.

Any line starting with a `#` is a comment, and is ignored.

Other blank lines in the file are ignored.

All other annotation is added using a tag/value comment style. The
tag/value format is inherently extensible, and readily made
backwards-compatible; unrecognized tags will simply be ignored. Extra
annotation includes consensus and individual RNA or protein secondary
structure, sequence weights, a reference coordinate system for the
columns, and database source information including name, accession
number, and coordinates (for subsequences extracted from a longer
source sequence). See below for details.

It is usually easy to convert other alignment formats into a least
common denominator Stockholm format. For instance, SELEX, GCG's MSF
format, and the output of the CLUSTALW multiple alignment program are
all closely related interleaved formats.

## Syntax of Stockholm markup

There are four types of Stockholm markup annotation, for per-file,
per-sequence, per-column, and per-residue annotation:

* **`#=GF <tag> <s>`**

  Per-file annotation. `<s>` is a free format text line of annotation
  type `<tag>`. For example, `#=GF DATE April 1, 2000`. Can occur
  anywhere in the file, but usually all the `#=GF` markups occur in a
  header.

* **`#=GS <seqname> <tag> <s>`**

  Per-sequence annotation. `<s>` is a free format text line of
  annotation type `<tag>` associated with the sequence named
  `<seqname>`. For example, `#=GS seq1 SPECIES_SOURCE Caenorhabditis
  elegans`. Can occur anywhere in the file, but in single-block formats
  (e.g. the Pfam distribution) will typically follow on the line after
  the sequence itself, and in multi-block formats (e.g. HMMER output),
  will typically occur in the header preceding the alignment but
  following the `#=GF` annotation.

* **`#=GC <tag> <..s..>`**

  Per-column annotation. `<..s..>` is an aligned text line of
  annotation type `<tag>`. `#=GC` lines are associated with a sequence
  alignment block; `<..s..>` is aligned to the residues in the
  alignment block, and has the same length as the rest of the block.
  Typically `#=GC` lines are placed at the end of each block.

* **`#=GR <seqname> <tag> <..s..>`**

  Per-residue annotation. `<..s..>` is an aligned text line of
  annotation type `<tag>`, associated with the sequence named
  `<seqname>`. `#=GR` lines are associated with one sequence in a
  sequence alignment block; `<..s..>` is aligned to the residues in
  that sequence, and has the same length as the rest of the block.
  Typically `#=GR` lines are placed immediately following the aligned
  sequence they annotate.

## Semantics of Stockholm markup

Any Stockholm parser will accept syntactically correct files, but is
not obligated to do anything with the markup lines. It is up to the
application whether it will attempt to interpret the meaning (the
semantics) of the markup in a useful way. At the two extremes are the
Belvu alignment viewer and the HMMER profile hidden Markov model
software package.

Belvu simply reads Stockholm markup and displays it, without trying to
interpret it at all. The tag types (`#=GF`, etc.) are sufficient to
tell Belvu how to display the markup: whether it is attached to the
whole file, sequences, columns, or residues.

HMMER (and other Easel-based software, including Infernal) uses
Stockholm markup to pick up a variety of information from the Pfam
multiple alignment database. The Pfam (and Rfam) consortia therefore
agree on additional syntax for certain tag types, so HMMER and
Infernal can parse some markups for useful information. This
additional syntax is imposed by software and databases, not by
Stockholm format per se. You can think of Stockholm as akin to XML,
and what my software reads as akin to an XML DTD, if you're into that
sort of structured data format lingo.

The Stockholm markup tags that are parsed semantically by Easel are as
follows:

### Recognized #=GF annotations

* **`ID  <s>`**

  Identifier. `<s>` is a name for the alignment; e.g. "rrm". One word.
  Unique in file.

* **`AC  <s>`**

  Accession. `<s>` is a unique accession number for the alignment;
  e.g. "PF00001". Used by the Pfam database, for instance. Often an
  alphabetical prefix indicating the database (e.g. "PF") followed by a
  unique numerical accession. One word. Unique in file.

* **`DE  <s>`**

  Description. `<s>` is a free format line giving a description of the
  alignment; e.g. "RNA recognition motif proteins". One line. Unique in
  file.

* **`AU  <s>`**

  Author. `<s>` is a free format line listing the authors responsible
  for an alignment; e.g. "Bateman A". One line. Unique in file.

* **`GA  <f> <f>`**

  Gathering thresholds. Two real numbers giving HMMER bit score
  per-sequence and per-domain cutoffs used in gathering the members of
  Pfam full alignments. See Pfam and HMMER documentation for more
  detail.

* **`NC  <f> <f>`**

  Noise cutoffs. Two real numbers giving HMMER bit score per-sequence
  and per-domain cutoffs, set according to the highest scores seen for
  unrelated sequences when gathering members of Pfam full alignments.
  See Pfam and HMMER documentation for more detail.

* **`TC  <f> <f>`**

  Trusted cutoffs. Two real numbers giving HMMER bit score per-sequence
  and per-domain cutoffs, set according to the lowest scores seen for
  true homologous sequences that were above the GA gathering
  thresholds, when gathering members of Pfam full alignments. See Pfam
  and HMMER documentation for more detail.

### Recognized #=GS annotations

* **`WT  <f>`**

  Weight. `<f>` is a nonnegative real number giving the relative weight
  for a sequence, usually used to compensate for biased representation
  by downweighting similar sequences. Usually the weights average 1.0
  (e.g. the weights sum to the number of sequences in the alignment)
  but this is not required. Either every sequence must have a weight
  annotated, or none of them can.

* **`AC  <s>`**

  Accession. `<s>` is a database accession number for this sequence.
  (Contrast to `#=GF AC` markup, which gives an accession for the whole
  alignment.) One word.

* **`DE  <s>`**

  Description. `<s>` is one line giving a description for this
  sequence. (Contrast to `#=GF DE` markup, which gives a description
  for the whole alignment.)

### Recognized #=GC annotations

* **`RF`**

  Reference line. Any character is accepted as a markup for a column.
  The intent is to allow labeling the columns with some sort of mark.

* **`SS_cons`**

  Secondary structure consensus. For protein alignments, DSSP codes or
  gaps are accepted as markup: `[HGIEBTSCX.-_]`, where H is alpha
  helix, G is 3/10-helix, I is p-helix, E is extended strand, B is a
  residue in an isolated b-bridge, T is a turn, S is a bend, C is a
  random coil or loop, and X is unknown (for instance, a residue that
  was not resolved in a crystal structure). For RNA alignments, the
  annotation is in WUSS format. Minimally, the symbols `<` and `>`
  indicate a base pair, `.` indicate single-stranded positions, and RNA
  pseudoknots are represented by alphabetic characters, with upper case
  letters representing the 5' side of the helix and lower case letters
  representing the 3' side. Note that this limits the annotation to a
  maximum of 26 pseudoknots per sequence.

* **`SA_cons`**

  Surface accessibility consensus. 0-9, gap symbols, or X are accepted
  as markup. 0 means $<10\%$ accessible residue surface area, 1 means
  $<20\%$, 9 means $<100\%$, etc. X means unknown structure.

### Recognized #=GR annotations

* **`SS`**

  Secondary structure consensus. See `#=GC SS_cons` above.

* **`SA`**

  Surface accessibility consensus. See `#=GC SA_cons` above.
