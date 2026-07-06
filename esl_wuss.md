# esl_wuss: WUSS notation for RNA secondary structure

Easel supports RNA secondary structure annotation using a linear string
representation called "WUSS notation" (Washington University Secondary
Structure notation), as originally developed for Infernal and the Rfam
database.

WUSS notation is a superset of the common dot-bracket notation for RNA
secondary structures, where open- and close-bracket symbols (or
parentheses) are used to annotate base pairing partners: for example,
`((((...))))`  indicates a four-base stem with a three-base loop.  For
anything much larger than a simple stem-loop, dot-bracket notation is
difficult for humans to look at.  WUSS notation makes it somewhat
easier to interpret the annotation for larger structures -- such as in
structural alignments output by Infernal.

The following example shows the key elements of WUSS notation, giving
the WUSS notation string for a small structure composed of two
hairpins sticking out of a multiloop closed by a stem:

```
  ::((((,<<<___>>>,,,<<-<<____>>-->>,))-))
  AACGGAACCAACAUGGAUUCAUGCUUCGGCCCUGGUCGCG
```

## Full (output) WUSS notation

Symbols used by WUSS notation in _output_ structure
annotation strings are as follows:

* **Base pairs**

  Base pairs are annotated by nested matching pairs of symbols `<>`,
  `()`, `[]`, or `{}`. The different symbols indicate the depth of
  the helix in the RNA structure: `<>` are used for simple
  terminal stems; `()` are used for "internal" helices enclosing a
  multifurcation of all terminal stems; `[]` are used for internal
  helices enclosing a multifurcation that includes at least one
  annotated `()` stem already; and `{}` are used for all internal
  helices enclosing deeper multifurcations.

* **Hairpin loops**

  Hairpin loop residues are indicated by underscores, `_`. Simple stem
  loops stand out as, e.g. `<<<<____>>>>`.

* **Bulge, interior loops**

  Bulge and interior loop residues are indicated by dashes, `-`.

* **Multifurcation loops**

  Multifurcation loop residues are indicated by commas, `,`. The
  mnemonic is "stem 1, stem 2", e.g. `<<<___>>>,,<<<___>>>`.

* **External residues**

  Unstructured single stranded residues completely outside the
  structure (unenclosed by any base pairs) are annotated by colons,
  `:`.

* **Insertions**

  Insertions relative to a known structure are indicated by periods,
  `.`. Regions where local structural alignment was invoked, leaving
  regions of both target and query sequence unaligned, are indicated by
  tildes, `~`. These symbols only appear in alignments of a known
  (query) structure annotation to a target sequence of unknown
  structure.

* **Pseudoknots**

  WUSS notation allows pseudoknots to be annotated as pairs of upper
  case/lower case letters: for example, `<<<<_AAAA____>>>>aaaa`
  annotates a simple pseudoknot; additional pseudoknotted stems could
  be annotated by `Bb`, `Cc`, etc.

  This is not a general notation, since it is limited to 26
  pseudoknotted stems, but it has sufficed so far.


An example of WUSS notation for a complicated structure (_E. coli_
RNase P) is shown below. The structure is from Jim Brown's RNase P
database [Brown99]. The WUSS notation annotates the _E. coli_ RNase P
sequence; note that the P4 and P6 pseudoknots are annotated as A's and
B's.

```
           {{{{{{{{{{{{{{{{{{,<<<<<<<<<<<<<-<<<<<____>>>>>>>>>->>>>>>>>
         1 GAAGCUGACCAGACAGUCGCCGCUUCGUCGUCGUCCUCUUCGGGGGAGACGGGCGGAGGG 60

           >,,,,AAA-AAAAA[[[[---BBBB-[[[[[<<<<<_____>>>>><<<<____>>>->(
        61 GAGGAAAGUCCGGGCUCCAUAGGGCAGGGUGCCAGGUAACGCCUGGGGGGGAAACCCACG 120

           (---(((((,,,,,,,,,,,,<<<<<--<<<<<<<<____>>>>>->>>>>>-->>,,,,
       121 ACCAGUGCAACAGAGAGCAAACCGCCGAUGGCCCGCGCAAGCGGGAUCAGGUAAGGGUGA 180

           ,,,<<<<<<_______>>>>>><<<<<<<<<____>>>->>>>>->,,)))--))))]]]
       181 AAGGGUGCGGUAAGAGCGCACCGCGCGGCUGGUAACAGUCCGUGGCACGGUAAACUCCAC 240

           ]]]]]],,,<<<<------<<<<<<----<<<<<_bbbb>>>>>>>>>>>----->>>>,
       241 CCGGAGCAAGGCCAAAUAGGGGUUCAUAAGGUACGGCCCGUACUGAACCCGGGUAGGCUG 300

           ,,,,,<<<<<<<<____>>>>>>>>,,,,,,,,,,}}}}}}}----------aaaaaaaa
       301 CUUGAGCCAGUGAGCGAUUGCUGGCCUAGAUGAAUGACUGUCCACGACAGAACCCGGCUU 360

           -}-}}}}}}}}}}::::
       361 AUCGGUCAGUUUCACCU 377
```

A second example illustrates the use of the local alignment annotation
symbols. It shows a local structural alignment of _B. subtilis_ RNase P
to _E. coli_ RNase P. The _E. coli_ query structure (from Jim Brown's
RNase P database [Brown99]) is aligned to the _B. subtilis_ sequence.
The local structural alignment is in four pieces; three regions of the
structure (102, 37, and 64 nt long) are skipped over, and one
additional stem is treated as a 40 nt insertion. The output below is
from Infernal:

```
hit 0   :      4    399    52.56 bits
           {{{{{{{{{{{{{{{{{{,<<<<<<<<<<<<<-<<<<<____>>>>>>>>>->>>>>>>>
         1 ggAGuggGgcaGgCaguCGCugcuucggccuuGuucaguuaacugaaaaggAccgaagga 60
           +: :::G::C:GG:A:UCGCU+C::::            U+            ::::G+A
         4 CUUAACGUUCGGGUAAUCGCUGCAGAUC-----------UUG----------AAUCUGUA 42

           >,,,,,,,,,,,,,[[[.[--------[[[[[~~~~~~~((---(((((,,,,~~~~~~)
        61 GAGGAAAGUCCGGGCUC.CACAGGGCAgGGUG*[ 29]*GGAAAGUGCCACAG*[96]*G 229
           GAGGAAAGUCC  GCUC C  A GG   :G G       :GAAAGUGCCACAG      G
        43 GAGGAAAGUCCAUGCUCgC--ACGGUGCUGAG*[102]*UGAAAGUGCCACAG*[37]*G 226

           ))--))))]]]]]].]]],,,~~~~~~,,,,,,,,,,}}}}}}}--..............
       230 GUAAACCCCACCcG.GAGCAA*[77]*CuAGAUGAAUGacuGcCCA.............. 344
           GUAAACC:C C: G GAG AA       UAGAU++AUGA:U:CC
       227 GUAAACCCCUCGAGcGAGAAA*[64]*GUAGAUAGAUGAUUGCC--gccugaguacgagg 342

           ..........................-----------------}-}}}}}}}}}}::::
       345 ..........................CGACAGAACCCGGCUUAuagcCccaCUccucuu 377
                                       ACA AAC  GGCUUA:AG::C::: :+ C
       343 ugaugagccguuugcaguacgaugga--ACAAAACAUGGCUUACAGAACGUUAGACCAC 399
```

## Shorthand (input) WUSS notation

While WUSS notation makes it easier to visually interpret Infernal
_output_ structural annotation, it would be painful to require people
to _input_ all structures in full WUSS notation. Therefore when
software like Infernal reads input secondary structure annotation, it
also allows simpler rules:

* **Base pairs**

  Any matching nested pair of `()`, `<>`, `[]`, `{}` symbols indicates
  a base pair; the exact choice of symbol has no meaning, so long as
  the left and right partners match up. Similarly, pseudoknotted pairs
  can also be annotated by matching nested pairs of any alphabet
  character, such as `Aa`, `Bb`, etc.

* **Single stranded residues**

  All other symbols `_-,:.~` indicate single stranded residues. The
  choice of symbol has no special meaning. Annotated pseudoknots
  (nested matched pairs of upper/lower case alphabetic characters) are
  also interpreted as single stranded residues in Infernal input.

Thus, for instance, `<<<<....>>>>` and `((((____))))` and
`<(<(._._)>)>` all indicate a four base stem with a four base loop (the
last example is legal but weird).

Remember that the key property of canonical (nonpseudoknotted) RNA
secondary structure is that the pairs are _nested_. `((<<....))>>` is
not a legal annotation string: the pair symbols don't match up
properly.

Because many other RNA secondary structure analysis programs use a
simple bracket notation for annotating structure, the ability to input
the simple format makes it easier to use data generated by other RNA
software packages. Conversely, converting output WUSS notation to
simple bracket notation is a matter of a simple Perl or sed script,
substituting the symbols appropriately.

## References

* **[Brown99]** JW Brown, ["The Ribonuclease P
  Database"](https://pubmed.ncbi.nlm.nih.gov/9847214/), _Nucleic Acids
  Research_ 27:314 (1999).

* **[RivasEddy99]** E Rivas, SR Eddy, ["A Dynamic Programming Algorithm
  for RNA Structure Prediction Including
  Pseudoknots"](https://doi.org/10.1006/jmbi.1998.2436), _Journal of
  Molecular Biology_ 285:2053-2068 (1999).
