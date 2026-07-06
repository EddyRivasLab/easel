# esl_randomseq - sampling random sequences

The `randomseq` module contains routines for generating,
shuffling, and randomizing sequences.

"Generating" means sampling a sequence from a given statistical
distribution. "Shuffling" means taking a given sequence and
randomizing it in some way that preserves at least the exact residue
composition of that sequence, and possibly other higher-order
statistics. "Randomizing" means sampling a sequence from a
statistical distribution estimated from a given sequence.

Routines like this are useful in testing sequence analysis algorithms;
for example, calculating expected score statistics under simple null
models.

When shuffling sequences, it is desirable to sample uniformly among
possible permutations. Many shuffling algorithms (even some published
ones [Fitch83]) are not correct in this sense. Easel's
monoresidue sequence shuffling implements the so-called Fisher/Yates
algorithm (Knuth's "Algorithm P") [Durstenfeld64, Knu-81a]. A
nontrivial additional power of the module is the ability to shuffle
sequences while preserving exact diresidue composition; the
`esl_rsq_CShuffleDP()` and `esl_rsq_XShuffleDP()`
routines implement the Altschul/Erickson method
[AltschulErickson86].[^1]

A more efficient method than Altschul/Erickson is known
[KandelWinkler96, Coward99] but it has not yet been implemented for
Easel.

The base routines work on any alphabetic text string. Augmentation
with the `alphabet` module adds routines for shuffling
digitized sequences.


## References

* **[AltschulErickson86]** SF Altschul, BW Erickson, "Optimal Sequence Alignment Using Affine Gap Costs," Bull. Math. Biol. 48:603-616 (1986).
* **[Coward99]** E Coward, ["Shufflet: Shuffling Sequences While Conserving the k-let Counts"](https://pubmed.ncbi.nlm.nih.gov/10745997/), Bioinformatics 15:1058-1059 (1999).
* **[Durstenfeld64]** R Durstenfeld, "Algorithm 235: Random Permutation," CACM 7:420 (1964).
* **[Fitch83]** WM Fitch, ["Random Sequences"](https://pubmed.ncbi.nlm.nih.gov/6842586/), JMB 163:171-176 (1983).
* **[KandelWinkler96]** D Kandel, Y Matias, R Unger, P Winkler, "Shuffling Biological Sequences," Discrete Appl. Math. 71:171-185 (1996).
* **[Knu-81a]** DE Knuth, _The Art of Computer Programming: Seminumerical Algorithms_, vol 2, second edition, Addison-Wesley (1981).
