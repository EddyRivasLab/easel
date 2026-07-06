# esl_sq - single biological sequences

The `sq` module provides Easel's object for single biological
sequences: an `ESL_SQ`.

Sequence objects invariably become complicated, even though their
designer intends them to be simple. There's many things we want to do
with a sequence, and useful features naturally accrete over time. If a
library isn't careful to balance creeping featuritis against having an
easy way to start using the object in simple applications, then the
sequence object - possibly the most fundamental object of a
biosequence library - can become a barrier to anyone else actually
using the library. All those useful features won't matter much if you
can't figure out how to turn your sequence data into an object, or get
it back out. Easel expects you to have your own preferred way of
dealing with sequence data that's not necessarily Easel's way, so it
provides simple ways to create sequence objects from elemental (C
string) data, and simple ways to get elemental C strings back out.
This lets you minimize your exposure to Easel's more complicated
capabilities if you like.

The most basic use of an `ESL_SQ` object is to hold one
complete sequence, simply as a plain C string. A sequence may also
have a name, an accession, and a description line. This is called a
_text mode_ sequence. In text mode, Easel doesn't know
whether the sequence is DNA, RNA, protein, or something else; it's
just an ASCII character string. This limits some of Easel's more
powerful abilities, such as the ability to check the sequence for
errors, or to automatically deal with degenerate residue codes; but
it's a simple mode that's easy to start using.

Alternatively, a sequence may be in _digital mode_. In digital
mode, sequences are predigested and encoded into Easel's internal
format, which makes many sequence routines more robust, efficient, and
powerful.

In addition to storing a complete sequence, an `ESL_SQ` is
designed to be used in three other situations:

* to hold a _subsequence_ of a larger source sequence. The
  object maintains source and coordinate information necessary for
  crossreferencing the subsequence's coordinate system to the original
  source coordinate system.

* to hold a _window_ of a larger source sequence. This is
  like a subsequence, but is more specifically intended for reading a
  sequence from a file in overlapping windows. This avoids having to
  suck an entire chromosome (for example) into memory at any one
  time. The stored subsequence is composed of two segments, a
  _previous context_ that gets saved from the previous window,
  and a _new window_ of fresh residues. The size of both the
  context and the window are configurable at the time each new window
  is read.

* to hold only _information_ about a sequence, such as its
  name, its length, and its position in a file, excluding the sequence
  (and optional secondary structure annotation) itself. This is handy
  for example when indexing a sequence file, when we'd rather not read
  any (possibly prohibitively large) sequence into memory until after
  we've mapped out how big it is.

To keep all this straight, the object contains a bunch of internal
bookkeeping data.

Sequence objects are growable and reusable, for efficiency in memory
allocation. If you're going to go through many different sequences
sequentially, you would typically just allocate a single
`ESL_SQ` object and `esl_sq_Reuse()` it for each new
sequence, rather than creating and destroying a lot of objects.

A sequence object can also store an optional secondary structure
annotation line for the sequence, one character per residue.

An interface to the `ESL_MSA` multiple alignment object
provides the ability to extract single unaligned sequences from a
multiple alignment.

You would often use the `sq` module in conjunction with
`sqio`, which provides the ability to read and write
`ESL_SQ` objects from and to files.

## examples

### getting data in and out of an `ESL_SQ`

The easiest way to create a new `ESL_SQ` object is with the
`esl_sq_CreateFrom()` function, which just takes character
strings for a sequence and its name (and also, optionally, an
accession, description, and/or secondary structure annotation string).

You can also build up (and/or change and manipulate) the contents of
an `ESL_SQ` object by accessing the name, accession,
description, sequence, and structure annotation line more directly.

The `eslSQ_EXAMPLE` example at the end of `esl_sq.c` illustrates both
approaches.

A few things to notice about that example:

* Every sequence has a name and a sequence. If we didn't want to
  add the optional accession, description, or structure annotation
  line, we'd pass `NULL` for those arguments to
  `esl_sq_CreateFrom()`.

* An RNA secondary structure annotation line is shown here as part
  of the example, but it's really sort of a more advanced
  feature. It's good to know it's there (see the `wuss` module
  for more information about how Easel annotates RNA structure) but
  you can ignore it if you're getting started.

* The `esl_sq_Set*` functions use the same syntax as C's
  `*printf()` family, which gives you a flexible way to create
  new sequence names, accessions, and descriptions automatically.

* The sequence in `sq->seq` is just a C string. (Here it's a
  copy of the `testseq` string.) That has a couple of
  implications. One is that it's a verbatim copy of what you provided;
  Easel doesn't know (or care) whether it's DNA or protein sequence,
  upper case or lower case, or if it contains illegal non-sequence
  characters. With a text mode sequence, that's _your_ problem!
  For more robustness and defined biosequence alphabets, read on below
  about digital mode sequences. The second implication is that, as a C
  string, the `n` residues are indexed `0..sq->n-1`, not
  `1..sq->n`.

* If you're going to directly copy a sequence of length `n`
  into a `sq->seq` field, note the `esl_sq_GrowTo()`
  call, which makes sure the sequence object is allocated with enough
  space for `n` residues; and don't forget to set `sq->n`.

* The structure annotation `sq->ss` is also a C string,
  indexed identically to `sq->seq`, but it's optional, and isn't
  allocated by default; `esl_sq_GrowTo()` calls will only
  reallocate for the structure annotation string after it's been
  allocated at least once. Hence the `esl_strdup` call in the
  example, which duplicates (allocates and copies) the annotation into
  `sq->ss`.

To get simple character strings back out of an `ESL_SQ` object,
you're encouraged to peek inside the object. (Yeah, I know, object
oriented design says that there should be methods for this,
independent of the object's implementation; but I balance that against
simplicity, and here, simplicity wins.) The object is defined and
documented in `esl_sq.h`. 

Ignore the `dsq` field for now; we're about to get to it, when
we talk about digital mode sequences.

The `ESL_SQ` object itself doesn't particularly care about the
contents of the annotation text fields, so long as they're C strings, and so
long as `n` is the length of the `seq` (and optional
`ss`, if it's non-`NULL`) strings. However, sequence file
formats do impose some expectations on the annotation strings, and it
would be a Good Idea to adhere to them:

* **name** — A sequence name is almost always expected to be
  a single "word" (no whitespace), like `SNRPA_HUMAN`.

* **acc** — An accession is also usually expected to be a
  single "word" with no whitespace, like `P09012`. Database
  accessions only make sense if you know what database they're for, so
  when sequences might be from different databases, you'll sometimes
  see accessions prefixed with a code indicating the source database,
  as in something like `UniProt:P09012`. Again, Easel itself
  isn't enforcing the format of this string, so your application is
  free to create its own accession/version/database format as needed.

* **desc** — A description line is something like `U1
  small nuclear ribonucleoprotein A (U1 snRNP protein A) (U1A protein)
  (U1-A).`; a one-line summary of what the sequence is. You can expect
  the description line to show up in the tabular output of sequence
  analysis applications, so ideally you want it to be short and sweet
  (so it fits on one line with a name, accession, score, coords, and
  other information from an analysis app). You also don't want the
  description line to end in a newline (`\n`) character, or the
  description line will introduce unexpected line breaks in these
  tabular output files.

You can reach into a `ESL_SQ` and copy or modify any of these
strings, but don't try to overwrite them with a larger string unless
You Know What You're Doing. Their memory allocations are managed by
the `ESL_SQ` object. Instead, use the appropriate
`esl_sq_Set*` function to overwrite an annotation field.

The `sq` module isn't much use by itself; it's a building block
for several other modules. For example, one of the main things you'll
want to do with sequences is to read them from a file. For examples
and documentation of sequence input, see the `sqio` module.

### using a digital `ESL_SQ`

What follows might make more sense if you've read about the
`alphabet` module first. `alphabet`'s documentation
explains how Easel uses an internal digital biosequence "alphabet",
where residues are encoded as small integers, suitable for direct use
as array indices. 

The `eslSQ_EXAMPLE2` example at the end of `esl_sq.c` creates and
accesses a digital mode sequence.

Things to notice about that example code:

* An `ESL_SQ` object has a `sq->seq` if it's in text
  mode, and `sq->dsq` if its in digital mode. These two fields are
  mutually exclusive; one of them is `NULL`.

* If you looked at the contents of `sq->dsq` in either of
  the objects, you'd see that each residue is encoded as a value
  `0..3`, representing (for an RNA alphabet) the residues
  `ACGU`.

* That representation is defined by the digital RNA alphabet
  `abc`, which was the first thing we created.

* In digital mode, both the sequence residues and the optional
  secondary structure characters are indexed `1..n`.

* To make the digital sequence in the first sequence object, we
  created a digital sequence `dsq` by encoding the
  `testseq` using `esl_abc_CreateDsq()`; this
  function allocated new memory for `dsq`, so we have to
  free it. An `ESL_DSQ *` is just a special character array;
  it's not a full-fledged Easel object, and so there's no
  conventional `Create()`,`Destroy()` function pair.

* In the second sequence object, we used
  `esl_abc_Digitize()` to encode the `testseq` directly
  into space that the `sq2` object already had allocated, saving
  us the temporary allocation of another `dsq`, because we
  created it in digital mode (`esl_sq_CreateDigital()`) and
  made it big enough to hold `n` digital residues with
  `esl_sq_GrowTo()`. Notice that `esl_sq_GrowTo()` is
  smart enough to know whether to grow the digital or the text mode
  sequence field.

* By convention, when using digital sequences, we usually keep
  track of (and pass as arguments) both a digital sequence `dsq`
  and its length `n`, and we also need to have the digital
  alphabet itself `abc` available to know what the `dsq`
  means; with text mode sequences, we usually just pass the string
  pointer. Thus the `esl_sq_CreateDigitalFrom()` function
  takes `abc`, `dsq`, and `n` as arguments, whereas
  the counterpart text mode `esl_sq_CreateDigitalFrom()` only
  took a C string `seq`. This is solely a convention - digital
  sequences begin and end with a special sentinel character, so we
  could always count the length of a `dsq` if we had to (using
  `esl_abc_dsqlen()`, for example), much as we can use ANSI
  C's `strlen()` to count the number of chars in a C string up
  to the sentinel `\0` `NUL` character at the end.

* To get the structure annotation to be indexed `1..n` for
  consistency with the `dsq`, even though the annotation string
  is still just an ASCII string, it's offset by one, and the leading
  character is set by convention to a `\0`. Therefore to access
  the whole structure string (for printing, for instance), you want to
  access `sq->ss+1`. This is a hack, but it's a simple one, so
  long as you don't forget about the convention.

* Because the original sequence has been encoded, you may not get
  the original sequence back out when you decode the digital values as
  alphabet symbols. `abc->sym[sq2->dsq[3]]`, for example, takes
  the third digital residue and looks it up in the alphabet's symbol
  table, returning the canonical character it's
  representing. Upper/lower case distinctions are lost, for example;
  digital alphabet symbol tables are uniformly upper case. And this
  example shows another example, where the input `testseq`
  contains T's, but since the digital alphabet was declared as RNA,
  the symbol table represents those residues as U's when you access
  them.

* In that respect, a more careful example should have checked the
  return status of the `esl_abc_CreateDsq()` and
  `esl_abc_Digitize()` calls. These have a normal failure
  mode, when the input text sequence contains one or more ASCII
  characters that are unrecognized and therefore invalid in the
  digital alphabet. If this had happened, these functions would have
  returned `eslEINVAL` instead of `eslOK`. We can get away
  without checking, however, because the functions just replace any
  invalid character with an "any" character (representing `N`
  for DNA or RNA, `X` for protein).

For more information about how digital sequence alphabets work, see
the `alphabet` module.
