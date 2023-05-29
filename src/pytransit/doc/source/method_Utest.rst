
Mann-Whitney U-test (utest)
===========================

This is a method for comparing datasets from a TnSeq library evaluated in
two different conditions, analogous to resampling.
This is a *rank-based* test on whether the level of insertions in a
gene or chromosomal region are significantly higher or lower in one
condition than the other.  Effectively, the insertion counts at the TA
sites in the region are pooled and sorted.  Then the combined ranks of the counts
in region A are compared to those in region B, and p-value is calculated
that reflects whether there is a significant difference in the ranks.
The advantage of this method is that it is less sensitive to outliers
(a unusually high insertion count at just a single TA site).
A reference for this method is `(Santa Maria et al., 2014)
<https://www.ncbi.nlm.nih.gov/pubmed/25104751>`__.


|

.. rst-class:: transit_sectionend
----
