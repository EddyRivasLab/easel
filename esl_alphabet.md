# esl_alphabet - digitized seqs

The `alphabet` module contains routines for digitizing alphabetic
biosequences.

It is convenient to represent nucleotides and amino acids as array
indices 0..3 or 0..19, respectively, for efficiency and other
reasons. It is also convenient to index the residues in a biosequence
in 1..L coordinates instead of the C language's 0..L-1 array
representation, partly for human readability, and also because some
codes (dynamic programming alignment algorithms, for example) have
boundary conditions where initializing a boundary at coordinate 0 is
convenient.

Real biosequences do not consist of just four or twenty different
canonical symbols, though. The `alphabet` module provides
mechanisms for dealing with several other biosequence coding issues:

* Degenerate residue symbols representing uncertainties,
  including both IUPAC/IUBMB standard one-letter nomenclature
  and nonstandard extensions such as the use of `J` to
  mean isoleucine or leucine (`I|L`) in protein sequences
  determined by mass spec;

* Symbols for "any residue" (N in DNA/RNA; X in amino acid
  sequences) and "not a residue" (a "translated" stop codon '*'
  in protein sequence;

* Standard and nonstandard symbols for unusual residues, such as
  selenocysteine (`U`) and pyrrolysine (`O`) in
  proteins;

* _Ad hoc_ symbols representing modified residues, such as
  the slew of one-letter codes used to annotate
  posttranscriptionally modified nucleotides in the Sprinzl tRNA
  database [Sprinzl98];

* Case-insensitivity of input sequences, for instance allowing
  both `a` and `A` to mean alanine in amino acid
  sequences;

* Tolerating common malpractices in the field, like the use of
  `X` instead of `N` as a degeneracy code in nucleic
  acid sequence;

* The semantic difference between a gap symbol and a missing
  data symbol in sequence alignments.

The `alphabet` module provides standard defaults for protein,
RNA, and DNA alphabets which follow both community standards and
IUPAC/IUBMB nomenclature for representing sequence residues in
one-letter ASCII characters. Additionally, the design of the
`alphabet` module is flexible enough to allow an application to
customize its own alphabet, to deal with these issues almost any way
it chooses.

Easel maintains alphabet information in an `ESL_ALPHABET`
structure. 

A digitized sequence `dsq` is an `ESL_DSQ *` (`unsigned char *`) array
of length `L+2`, where `dsq[1..L]` are digitized residues, and
`dsq[0]` and `dsq[L+1]` are sentinel bytes (of value `eslSENTINEL`,
127).

The `eslALPHABET_EXAMPLE` example in `esl_alphabet.c` creates a DNA
alphabet and digitizes a short sequence.


## terminology

A _symbol_ is a 7-bit ASCII input character, representing a
residue, gap, nonresidue, or degeneracy. A _code_ is the
digital internal representation of the symbol as an `unsigned char`
in the range $0..127$, suitable for use as an array index. The
`alphabet` module translates input symbols into internal
digital codes.

We distinguish between an input alphabet, an internal alphabet, and a
canonical alphabet. The _input alphabet_ consists of all the
symbols that Easel allows in an input biosequence character
string. The _internal alphabet_ is the standardized one-letter
alphabet that Easel deals with. The _canonical alphabet_ is the
fundamental set of 4 nucleotides or 20 amino acids.

Easel deals with all of the complications of sequence encoding using
two concepts, equivalency and degeneracy. _Equivalency_
defines how the input alphabet maps to the internal
alphabet. _Degeneracy_ defines how the internal alphabet maps
to the canonical alphabet.

Equivalent residues are symbols that are accepted in an input sequence
character string and silently translated into an appropriate internal
code. Characters in the input alphabet are mapped many-to-one to the
internal alphabet using an _input map_. One use of equivalency
is to map both lower and upper case input to the same internal
symbol. Another use is to allow several different input characters to
mean a gap, 'any' symbol, or 'nonresidue' symbol. Another use is to
silently accept and "fix" nonstandard but common input "errors",
such as tolerating the use of X to mean N in nucleic acid sequences.

Degenerate residues are codes in the internal alphabet that are mapped
one-to-many onto canonical residue codes, using a _degeneracy
map_. In addition to mapping the degeneracy codes onto the canonical
alphabet, the degeneracy mechanism is also used to deal with unusual
and modified residues. Selenocysteine, for instance, is represented by
default as a `U`, but treated as a degenerate code for `C`
(cysteine). The rationale for this will be described in more detail
below.

## the internal alphabet

Easel's internal alphabet is a string (`a->sym`) of length
`Kp`, which contains:

* the `K` symbols of the canonical alphabet;
* a standard gap symbol;
* (optionally) any other degenerate, unusual, or modified residue codes;
* a mandatory "any" symbol (a completely degenerate residue);
* a mandatory "not-a-residue" symbol;
* a mandatory "missing data" symbol.

Residues `0..K-1` must be the canonical alphabet. Residue
`K` must be the gap symbol. Residues `K+1..Kp-4` must be
the degenerate and modified residue symbols (there can be zero of
these). Residue `Kp-3` must be the completely degenerate symbol
(such as `X` for protein sequence or `N` for nucleic acid
sequence); all alphabets must have such a symbol. Residue `Kp-2`
must be the not-a-residue symbol. Residue `Kp-1` must be the
missing data symbol. Because the 'any' symbol, 'not-a-residue'
symbol, and the two kinds of gap symbols are mandatory in any
alphabet, `Kp` $\geq$ `K+4`. Aside from these
constraints, symbols may occur in any order.

The digital code used for each residue is then the index of a residue
in this string, `0..Kp-1`. The only other value that can appear
in a digitized sequence is `eslSENTINEL` (127), the sentinel
byte in positions `0` and `L+1` of a digitized sequence of
length `L`.

The rationale for the ordering is the following. Most applications
will define residue scores in vectors and matrices that are smaller
than the full range of the internal alphabet; for instance, it's
common to only have `K` scores for the canonical residues. As
much as possible, we want array indices to be the same whether we're
accessing the full internal alphabet or a smaller score vector or
matrix. So: we expect many applications to have score vectors or
matrices that only contain the `K` canonical residues, so the
canonical residues go first. We expect some applications to treat
gaps as an extra symbol, and provide `K+1` position-specific
scores or a `K+1` $\times$ `K+1` score matrix, so the gap
character is next. We expect a few applications to optimize degeneracy
scoring by precalculating them in `Kp-2` vectors or $Kp-2 \times
Kp-2$ matrices, so the degeneracies go next (the gap character at $K$
might then go unused in the score vectors and matrices, but that's a
minor inefficiency). The 'any' symbol should be at a predictable
position in the degeneracy list, so it's arbitrarily placed at the
end, in position `Kp-3`. The most robust applications will also
handle the not-a-residue symbol (they may see translated stop codons),
so it's next. Finally, the missing data symbol is expected to always
require special handling when it occurs, rather than appearing in a
score vector or matrix, so it's put last.

## the standard alphabets: DNA, RNA, protein

The three standard internal alphabets are:

| Type | `sym` | equivs | gaps | `K` | `Kp` |
|------|-------|--------|------|-----|------|
| `eslRNA` | `ACGU-RYMKSWHBVDN*~` | T=U; X=N | `-_.` | 4 | 18 |
| `eslDNA` | `ACGT-RYMKSWHBVDN*~` | U=T; X=N | `-_.` | 4 | 18 |
| `eslAMINO` | `ACDEFGHIKLMNPQRSTVWY-BJZOUX*~` |  | `-_.` | 20 | 29 |

The `sym` string contains all the symbols that can be handled
internally, and all the residues that can be represented when a
digitized sequence is converted back to text. An application might
still convert some characters for its own purposes before displaying
an alphabetic string; for instance, to use different gap symbols for
insertions versus deletions, or to use upper/lower case conventions to
represent match/insert positions in a profile HMM alignment.

The standard DNA and RNA alphabets follow published IUBMB
recommendations ("Nomenclature for incompletely specified bases in
nucleic acid" [IUBMB85]), with an addition of X as an
equivalence for N (acquiescing to the _de facto_ BLAST filter
standard of using X's to mask residues), and equivalencing T to U in
RNA sequences (and vice versa in DNA).

The one-letter code for amino acids follows section 3AA-21 of the
IUPAC recommendations [IUPAC84]. The code is augmented by U for
selenocysteine, as recommended in 1999 by the JCBN/NC-IUBMB Newsletter
(<http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html>).
It is also augmented by O for pyrrolysine and J for a
leucine/isoleucine ambiguity (from a mass spectrometry experiment),
following usage in the RESID database
(<http://www.ebi.ac.uk/RESID/>).

## degenerate residues

The symbols from `K+1..Kp-4` in the internal alphabet are all
treated as degenerate residues.

When creating a custom alphabet, each degenerate symbol is initialized
by calling `esl_alphabet_SetDegeneracy(alphabet, c, syms)` to
assign degenerate alphabetic symbol `c` to the alphabetic
symbols in the string `syms`. For example,
`esl_alphabet_SetDegeneracy(a, 'R', "AG")` assigns R
(purine) to mean A or G. For the standard biosequence alphabets, this
is done automatically to define the proper degeneracy codes.

For amino acid alphabets, the default code is:

```c
  esl_alphabet_SetDegeneracy(a, 'B', "ND");
  esl_alphabet_SetDegeneracy(a, 'J', "IL");
  esl_alphabet_SetDegeneracy(a, 'Z', "QE");
```

For RNA alphabets, the default code is:

```c
  esl_alphabet_SetDegeneracy(a, 'R', "AG");
  esl_alphabet_SetDegeneracy(a, 'Y', "CU");
  esl_alphabet_SetDegeneracy(a, 'M', "AC");
  esl_alphabet_SetDegeneracy(a, 'K', "GU");
  esl_alphabet_SetDegeneracy(a, 'S', "CG");
  esl_alphabet_SetDegeneracy(a, 'W', "AU");
  esl_alphabet_SetDegeneracy(a, 'H', "ACU");
  esl_alphabet_SetDegeneracy(a, 'B', "CGU");
  esl_alphabet_SetDegeneracy(a, 'V', "ACG");
  esl_alphabet_SetDegeneracy(a, 'D', "AGU");  
```

For DNA alphabets, the calls are is the same as for RNA code, but with
`T` in place of `U`.

### implementation: the degeneracy map

The alphabet's degeneracy map is implemented in an array
`a->degen[0..Kp-1][0..K-1]` of 1/0 (TRUE/FALSE) flags.
`a->degen[x][y] == TRUE` indicates that the residue set $D(x)$
for degeneracy code `x` contains base residue `y`.
`a->ndegen[x]` contains the cardinality $|D(x)|$, how many base
residues are represented by degeneracy code `x`.

For the two kinds of gap symbols, the degeneracy map is empty; all
flags are FALSE and the cardinality is 0.

Because character `Kp-3` in the internal alphabet is
automatically assumed to be an "any" character (such as 'N' for DNA
or RNA, 'X' for protein), `a->degen[Kp-3][i] = 1` for all
$i=0..K-1$, and `a->ndegen[Kp-3] = K`.

The storage of the degeneracy map is a little wasteful. We really only
need rows `a->degen[K+1..Kp-3]`, but optimizing this would
create some index translation hassles, and it doesn't seem worth it.

## equivalent residues

The concept of equivalent residues allows an input symbol to be mapped
to a different internal symbol. One use of equivalence is to map both
lower and upper case input to the same internal representation.
Another use is to allow several different input characters to mean a
gap. Another use is to silently accept and "fix" nonstandard but
common input "errors", such as the use of T instead of U in RNA
sequences (or vice versa in DNA), or the use of X instead of N as an
ambiguity code in nucleic acid sequences.

The call `esl_alphabet_SetEquiv(a, 'U', 'T')`, for example,
makes an alphabet interpret `U` as a `T` (encoding both as
`3`, in the case of the standard DNA and RNA alphabets).

All three standard alphabets accept `_` or `.` symbols
as equivalences for the standard gap symbol `-`. An application
can define additional gap characters, such as `,`, by calling
`esl_alphabet_SetSynonym(a, ',', '-')` on one of the standard
alphabets to define additional equivalences (that is, you don't have
to create a custom alphabet to add new equivalences).

`esl_alphabet_SetCaseInsensitive()` maps both upper case and
lower case input alphabetic characters map to their equivalent in the
internal alphabet in a case-insensitive manner. This function works
only on residues that have already been declared to be part of the
alphabet, so when defining a custom alphabet, it must be called after
all individual equivalences have been defined. The standard alphabets
are always set to be case insensitive.

### implementation of equivalent residues: the input map

Internally, an **input map**, `a->inmap[0..127]`, specifies
how an input ASCII 7-bit text symbol is converted to digital
code. `a->inmap['T'] = 3` in the standard DNA alphabet, for
example, and the call `esl_alphabet_SetSynonym(a, 'U', 'T')`
sets `a->inmap['U'] = a->inmap['T']`.

The elements in input maps are of type `unsigned char`. Legal
values are 0..127 (values that can be cast to the `unsigned char`
codes in a digitized sequence) and two additional flags with
negative values, `eslILLEGAL_CHAR` (255) and
`eslIGNORED_CHAR` (254).

## unusual or modified residues

In addition to the canonical 4 or 20 residues and their ambiguity
codes, there are many unusual and/or modified residues. For instance,
there are many posttranscriptional or posttranslational modifications
on residues in RNAs and proteins. Some databases try to capture this
information in a single-letter alphabetic code, such as the Sprinzl
transfer RNA database [Sprinzl98].

Additionally, and perhaps more importantly, proteins are known to
contain at least two additional genetically encoded amino acids,
selenocysteine and pyrrolysine. Selenocysteine is represented by a
`U` according to the IUPAC standard, and pyrrolysine is
represented by a `O` in the RESID database at EBI.

Unusual one-letter residue codes pose a tradeoff issue for sequence
analysis applications. On the one hand, an application should
recognize symbols for unusual or modified residues, and be able to
represent them both internally and in any sequence output. For
example, no application should read an input selenocysteine residue
(`U`) and output it as a cysteine (`C`) -- this changes
the original sequence and causes data corruption.[^1] On the other
hand, most sequence analysis applications would
not want to take the trouble to define a canonical alphabet larger
than the usual 4 or 20 residues, and then have to parameterize that
alphabet, just to be able to handle a few rare residues. (Pyrrolysine,
for example, has only been found in a handful of proteins in a few
Archaea.) It is useful to be able to deal with probability parameters
and scores only for the canonical alphabet. However (on yet another
hand!) in some cases one _would_ want to write a specialized
application that parameterizes unusual residues as part of its
canonical alphabet -- for instance, an application for analyzing
posttranscriptional tRNA modifications, for example.

[^1]: However, at least one the main public protein databases
(UniProt) has already chosen to replace all selenocysteines with
`C` and all pyrrolysines with `K`, for fear of breaking legacy
sequence analysis software. So, this data corruption is already a fact
of life.

Therefore, Easel must not force an input selenocysteine or pyrrolysine
(or any other unusual residue) to be recoded as an arbitrary symbol
(such as cysteine or lysine). That is, unusual symbols cannot be
treated as equivalences, but must be allowed to be part of the
internal alphabet. However, Easel _can_ allow unusual symbols to
be treated as noncanonical, and _score_ them as some other
arbitrary residue, as a reasonable approximation. Thus for most
purposes, unusual symbols are best handled as a special kind of
degeneracy, with a one-to-one degeneracy map from the unusual symbol
to the "closest" canonical residue.

Therefore, the default amino acid alphabet accepts selenocysteine
(`U`) and pyrrolysine symbols (`O`) and represents them in
the internal alphabet, and maps them as "degeneracies" onto cysteine
(`C`) and lysine (`K`), respectively.

When that behavior is not suitable, an application can also define any
custom alphabet it chooses, as described below.

## creating a custom alphabet

The `eslALPHABET_EXAMPLE2` example in `esl_alphabet.c` creates a
customized 22-letter amino acid alphabet that includes the `U`
selenocysteine code and the `O` pyrrolysine code.


## scoring degenerate residues

To score a degenerate residue code $x$, Easel provides two strategies.
One set of functions assigns an average score:

$$
  S(x) =  \frac{\sum_{y \in D(x)}  S(y) } { |D(x)| },
$$

where $D(x)$ is the set of residues $y$ represented by degeneracy code
$x$ (for example, $D(\mathrm{R}) = \{ \mathrm{A,G} \}$), $| D(x) |$ is the
number of residues that the degeneracy code includes, and $S(y)$ is
the score of a base residue $y$. Because scores $S(y)$ are commonly
kept as integers, floats, or doubles, depending on the application,
three functions are provided that differ only in the storage type of
the scores: `esl_abc_IAvgScore(a,x,sc)`,
`esl_abc_FAvgScore(a,x,sc)`, and
`esl_abc_DAvgScore(a,x,sc)` calculate and return the average
score of residue `x` in alphabet `a` given a base score
vector `sc[0]..sc[K-1]` for integers, floats, and doubles,
respectively.

A second set of functions assigns an expected score, weighted by an
expected occurrence probability $p_y$ of the residues $y$ (often the
random background frequencies):

$$
  S(x) =  \frac{\sum_{y \in D(x)}  p_y S(y) } { \sum_{y \in D(x)} p_y },
$$

These three functions are `esl_abc_IExpectScore(a,x,sc,p)`,
`esl_abc_FExpectScore(a,x,sc,p)`, and
`esl_abc_DExpectScore(a,x,sc,p)`. For the integer and float
versions, the probability vector is in floats; for the double version,
the probability vector is in doubles.

For efficiency reasons, an application might choose to preculate
scores for all possible degenerate codes it might see. HMMER, for
example, turns probability vectors of length `K` into score
residues of length `Kp`.

An application might also choose to score residues on-the-fly, using
score vectors of length `K`. Each input residue `x` would
then have to be tested to see if it is degenerate, before scoring it
appropriately. `esl_abc_IsBasic(a, x)` returns `TRUE`
if `x` is in the basic set of `K` residues in alphabet
`a`, and `FALSE` otherwise. Similarly,
`esl_abc_IsGap(a,x)` tests whether $x$ is a gap, and
`esl_abc_IsDegenerate(a,x)` tests whether $x$ is a degenerate
residue.

## References

* **[IUBMB85]** A Cornish-Bowden, "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences," NAR 13:3021-3030 (1985).
* **[IUPAC84]** IUPAC-IUB Joint Commission on Biochemical Nomenclature (JCBN), "Nomenclature and Symbolism for Amino Acids and Peptides: Recommendations 1983," Biochem. J. 219:345-373 (1984).
* **[Sprinzl98]** M Sprinzl, C Horn, M Brown, A Ioudovitch, S Steinberg, ["Compilation of tRNA Sequences and Sequences of tRNA Genes"](https://pubmed.ncbi.nlm.nih.gov/9399820/), NAR 26:148-153 (1998).
