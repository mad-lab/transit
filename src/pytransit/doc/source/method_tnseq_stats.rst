.. _tnseq_stats:

tnseq_stats
========

This command is useful for calculating some metrics on input wig (or
combined_wig files) for assessement of data quality.

Similar information can be generated through the GUI via
the menu options 'View->Quality Control' (select samples in main window first).

Usage:
------

::

  > usage: python3 src/transit.py tnseq_stats <file.wig>+ [-o <output_file>]
           python3 src/transit.py tnseq_stats -c <combined_wig> [-o <output_file>]


It generates a table (tab-separated text file that can be opened in Excel) with the following statistics in it:

=============  ==============================================  =============================================================================================================
Column Header  Column Definition                                 Comments
=============  ==============================================  =============================================================================================================
dataset        Name of sample (wig file)
density        Fraction of sites with insertions.              "Well-saturated" Himar1 datasets have >30% saturation. Below this, statistical methods may have trouble.
mean_ct        Average read-count over all TA sites.
NZMean         Average read-count, excluding empty sites.       A value between 30-200 is usually good for Himar1 datasets. Too high or too low can indicate problems.
NZMedian       Median read-count, excluding empty sites.        As read-counts can often have spikes, median serves as a good robust estimate.
max_ct         Largest read-count at any TA site                Useful to determine whether there are outliers/spikes, which may indicate sequencing issues.
total_cts      Sum of total read-counts in the sample.          Indicates how much sequencing material was obtained. Typically >1M reads is desired for Himar1 datasets.
skewness       3rd-order moment of read-count distribution.     Sharp peak? Large skew may indicate issues with a dataset. Typically a skew < 50 is desired. May be higher when
kurtosis       4th-order moment of read-counts distribution     Lop-sided peak?
PTI            Pickand Tail Index                               Another statistical measure that indicates presence of individual TA sites with high outlier counts. PTI>1.0 is not good.
=============  ==============================================  =============================================================================================================



Here is an example:

::

  > python3 src/transit.py tnseq_stats -c src/pytransit/data/cholesterol_glycerol_combined.dat
  dataset density mean_ct NZmean  NZmedian        max_ct  total_cts       skewness        kurtosis        pickands_tail_index
  src/pytransit/data/cholesterol_H37Rv_rep1.wig   0.439   139.6   317.6   147     125355.5        10414005        54.8    4237.7  0.973
  src/pytransit/data/cholesterol_H37Rv_rep2.wig   0.439   171.4   390.5   148     704662.8        12786637        105.8   14216.2 1.529
  src/pytransit/data/cholesterol_H37Rv_rep3.wig   0.359   173.8   484.2   171     292294.8        12968502        42.2    2328.0  1.584
  src/pytransit/data/glycerol_H37Rv_rep1.wig      0.419   123.3   294.5   160     8813.3  9195672 4.0     33.0    0.184
  src/pytransit/data/glycerol_H37Rv_rep2.wig      0.516   123.8   240.1   127     8542.5  9235984 4.0     33.5    0.152


In this example, you can see the 5 samples have saturations in the
range of 35.9-51.6% (which is decent).  The NZMeans are in the range
123-139, but this is *post-normalization*.  (TTR normalization had
already been applied to this combined_wig file, so the means can be
expected to be scaled to around 100.)  If you want to see the NZmeans
for the *raw data*, re-generate the combined_wig file using '-n
nonorm', to skip the automatic normalization step.  These samples also
exhibit *skewness* that is on the high side (33.0-105.8).  This is
probably related to the fact that some individual TA sites have very
high counts.  For example, rep2 of cholesterol has a max count of
704662 at a single TA site, representing over 5% of the 12.78M total
insertion counts.  TTR is supposed to be robust by ignoring the top 5%
of most abundant sites during normalization, but still the rest of the
distribution of counts could be skewed.  This sample also has a high
*Pickands' tail index* of 1.53 (which is above 1.0), also suggesting
skew.  While, we currently don't have recommendations for hard cutoffs
to use for identifying bad samples (e.g. that might need to be
re-sequenced), I would say that skew>30 and/or PTI>1.0 are signs that
a sample might be noisy or lower-quality.  See :ref:`Quality Control
<transit_quality_control>` for more discussion about assessing quality
of TnSeq datasets.  Nonetheless, doing resampling on this data still
yielded insights into many genes required for cholesterol metabolism
in *M. tuberculosis* `(Griffin et al, 2009)
<https://pubmed.ncbi.nlm.nih.gov/21980284/>`_. 
