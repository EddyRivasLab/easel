# esl_scorematrix - Pairwise alignment scores

The `scorematrix` module implements standard residue pairwise
alignment scoring systems (such as BLOSUM matrices), and their
probabilistic interpretation [Altschul91, YuAltschul03].

The module implements one object,
`ESL_SCOREMATRIX`, which holds a score matrix with integer
scores.

The `eslSCOREMATRIX_EXAMPLE` example at the end of `esl_scorematrix.c`
reads a matrix file from disk, reverse engineers it [YuAltschul03] to
obtain its target probabilities, background frequencies, and lambda,
then prints information about it.



## References

* **[Altschul91]** SF Altschul, ["Amino Acid Substitution Matrices from
  an Information Theoretic Perspective"](https://pubmed.ncbi.nlm.nih.gov/2051488/),
  _JMB_ 219:555-565 (1991).

* **[YuAltschul03]** YK Yu, JC Wootton, SF Altschul, ["The Compositional
  Adjustment of Amino Acid Substitution Matrices"](https://pubmed.ncbi.nlm.nih.gov/14663142/),
  _PNAS_ 100:15688-15693 (2003).
