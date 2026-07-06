# esl_distance - distance calculations

The `distance` module implements routines for inferring
mutational distances between pairs of aligned sequences, and for
constructing distance matrices from multiple sequence alignments.

## definition of pairwise identity and pairwise difference

Given a pairwise sequence alignment of length $L$, between two
sequences of $n_1$ and $n_2$ residues ($n_1 \leq L$, $n_2 \leq L$),
where the $L$ aligned symbol pairs are classified and counted as
$c_{\mathrm{ident}}$ identities, $c_{\mathrm{mismat}}$ mismatches, and
$c_{\mathrm{indel}}$ pairs that have a gap symbol in either or both
sequences ($c_{\mathrm{ident}} + c_{\mathrm{mismat}} + c_{\mathrm{indel}} =
L$), _pairwise sequence identity_ is defined as:

$$
   \mathrm{pid} = \frac{c_{\mathrm{ident}}}{\mathrm{MIN}(n_1, n_2)},
$$

and _pairwise sequence difference_ is defined as

$$
   \mathrm{diff} = 1 - \mathrm{pid} = \frac{\mathrm{MIN}(n_1,n_2) - c_{\mathrm{ident}}}{\mathrm{MIN}(n_1, n_2)}.
$$

Both pid and diff range from 0 to 1.

In the unusual case where $\mathrm{MIN}(n_1,n_2)=0$ -- that is, one of
the aligned sequences consists entirely of gaps -- the percent
identity $0/0$ is defined as 0. The calculation is robust against
length 0 sequences, which do arise in real applications. (Not just in
bad input, either. For example, this arises when dealing with subsets
of the columns of a multiple alignment.)

There are many ways that pairwise identity might be calculated,
because there are a variety of choices for the denominator. In Easel,
identity calculations are used primarily in _ad hoc_ sequence
weight calculations for multiple sequence alignments, as part of
profile HMM or profile SCFG construction. Multiple alignments will
often contain short sequence fragments. We want to deal robustly with
cases where two short fragments may have little overlap, or none at
all. The most obvious calculation of pairwise identity,
$c_{\mathrm{ident}} / c_{\mathrm{ident}} + c_{\mathrm{mismat}}$, is not
robust, because alignments with few aligned residues (either because
they are highly gappy, or they are partially overlapping fragments)
may receive artifactually high identities. Other definitions,
$c_{\mathrm{ident}} / L$ or $c_{\mathrm{ident}} / \mathrm{MEAN}(n_1, n_2)$
or $c_{\mathrm{ident}} / \mathrm{MAX}(n_1, n_2)$ are also not robust,
sharing the disadvantage that good alignments of fragments to longer
sequences would be scored as artifactually low identities.


## generalized Jukes-Cantor distances

The Jukes-Cantor model of DNA sequence evolution assumes that all
substitutions occur at the same rate $\alpha$
[JukesCantor69]. It is a reversible, multiplicative evolutionary
model. It implies equiprobable stationary probabilities. The
_Jukes/Cantor distance_ is the maximum likelihood estimate of
the number of substitutions per site that have occurred between the
two sequences, correcting for multiple substitutions that may have
occurred the same site. Given an ungapped pairwise alignment of length
$L$ consisting of $c_{\mathrm{ident}}$ identities and
$c_{\mathrm{mismat}}$ mismatches (observed substitutions)
($c_{\mathrm{ident}} + c_{\mathrm{mismat}} = L$, the fractional observed
difference $D$ is defined as

$$
  D = \frac{c_{\mathrm{mismat}}}{c_{\mathrm{ident}} + c_{\mathrm{mismat}}},
$$

and the Jukes-Cantor distance $d$ is defined in terms of $D$ as:

$$
  d = -\frac{3}{4} \log \left( 1 - \frac{4}{3} D \right)
$$

The Jukes/Cantor model does not allow insertions or deletions. When
calculating "Jukes/Cantor distances" for gapped alignments, gap
symbols are simply ignored, and the same calculations above are
applied.

The Jukes-Cantor model readily generalizes from the four-letter DNA
alphabet to any alphabet size $K$, using the same definition of
observed fractional difference $D$. A _generalized Jukes-Cantor
distance_ is:

$$
  d = -\frac{K-1}{K} \log \left( 1 - \frac{K}{K-1} D \right).
$$

The large-sample variance of this estimate $d$ is:

$$
   \sigma^2 = e^\frac{2Kd}{K-1} \frac{D(1-D)}{L'}
$$

where $L'$ is the total number of columns counted, exclusive of gaps,
$L' = c_{\mathrm{ident}} + c_{\mathrm{mismat}}$.

If the observed $D \geq \frac{K-1}{K}$, the maximum likelihood
Jukes-Cantor distance is infinity, as is the variance. In this case,
both $d$ and $V$ are returned as `HUGE_VAL`.


## References

* **[JukesCantor69]** TH Jukes, CR Cantor, "Evolution of Protein
  Molecules," in HN Munro, ed., _Mammalian Protein Metabolism_,
  pp. 21-132, Academic Press (1969).
