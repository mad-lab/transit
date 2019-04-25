
.. _`analysis_methods`:

Analysis Methods
================


TRANSIT has analysis methods capable of analyzing **Himar1** and **Tn5** datasets.
Below is a description of some of the methods.

|

.. _gumbel:

Gumbel
------

The Gumbel can be used to determine which genes are essential in a
single condition. It does a gene-by-gene analysis of the insertions at
TA sites with each gene, makes a call based on the longest consecutive
sequence of TA sites without insertion in the genes, calculates the
probability of this using a Bayesian model.

.. NOTE::
   Intended only for **Himar1** datasets.

How does it work?
~~~~~~~~~~~~~~~~~

| For a formal description of how this method works, see our paper [DeJesus2013]_:

|  DeJesus, M.A., Zhang, Y.J., Sassettti, C.M., Rubin, E.J.,
  Sacchettini, J.C., and Ioerger, T.R. (2013).
| `Bayesian analysis of gene essentiality based on sequencing of transposon insertion libraries. <http://www.ncbi.nlm.nih.gov/pubmed/23361328>`_ *Bioinformatics*, 29(6):695-703.

Example
~~~~~~~

::

 python transit.py gumbel <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -b <integer>    :=  Number of Burn-in samples. Default -b 500
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -t <integer>    :=  Trims all but every t-th value. Default: -t 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0




Parameters
~~~~~~~~~~

-  **Samples:** Gumbel uses Metropolis-Hastings (MH) to generate samples
   of posterior distributions. The default setting is to run the
   simulation for 10,000 iterations. This is usually enough to assure
   convergence of the sampler and to provide accurate estimates of
   posterior probabilities. Less iterations may work, but at the risk of
   lower accuracy.

-  **Burn-In:** Because the MH sampler many not have stabilized in the
   first few iterations, a "burn-in" period is defined. Samples obtained
   in this "burn-in" period are discarded, and do not count towards
   estimates.

-  **Trim:** The MH sampler produces Markov samples that are correlated.
   This parameter dictates how many samples must be attempted for every
   sampled obtained. Increasing this parameter will decrease the
   auto-correlation, at the cost of dramatically increasing the
   run-time. For most situations, this parameter should be left at the
   default of "1".

-  **Minimum Read:** The minimum read count that is considered a true
   read. Because the Gumbel method depends on determining gaps of TA
   sites lacking insertions, it may be susceptible to spurious reads
   (e.g. errors). The default value of 1 will consider all reads as true
   reads. A value of 2, for example, will ignore read counts of 1.

-  **Replicates:** Determines how to deal with replicates by averaging
   the read-counts or summing read counts across datasets. This should
   not have an affect for the Gumbel method, aside from potentially
   affecting spurious reads.

|

Outputs and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~

The Gumbel method generates a tab-separated output file at the location
chosen by the user. This file will automatically be loaded into the
Results Files section of the GUI, allowing you to display it as a table.
Alternatively, the file can be opened in a spreadsheet software like
Excel as a tab-separated file. The columns of the output file are
defined as follows:

+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Column Header   | Column Definition                                                                                                             |
+=================+===============================================================================================================================+
| ORF             | Gene ID.                                                                                                                      |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Name            | Name of the gene.                                                                                                             |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Description     | Gene description.                                                                                                             |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| k               | Number of Transposon Insertions Observed within the ORF.                                                                      |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| n               | Total Number of TA dinucleotides within the ORF.                                                                              |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| r               | Length of the Maximum Run of Non-Insertions observed.                                                                         |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| s               | Span of nucleotides for the Maximum Run of Non-Insertions.                                                                    |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| zbar            | Posterior Probability of Essentiality.                                                                                        |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+
| Call            | Essentiality call for the gene. Depends on FDR corrected thresholds. E=Essential U=Uncertain, NE=Non-Essential, S=too short   |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------+

|
|  Note: Technically, Bayesian models are used to calculate posterior
  probabilities, not p-values (which is a concept associated with the
  frequentist framework). However, we have implemented a method for
  computing the approximate false-discovery rate (FDR) that serves a
  similar purpose. This determines a threshold for significance on the
  posterior probabilities that is corrected for multiple tests. The
  actual thresholds used are reported in the headers of the output file
  (and are near 1 for essentials and near 0 for non-essentials). There
  can be many genes that score between the two thresholds (t1 < zbar <
  t2). This reflects intrinsic uncertainty associated with either low
  read counts, sparse insertion density, or small genes. If the
  insertion\_density is too low (< ~30%), the method may not work as
  well, and might indicate an unusually large number of Uncertain or
  Essential genes.

|

Run-time
~~~~~~~~

The Gumbel method takes on the order of 10 minutes for 10,000 samples.
Run-time is linearly proportional to the 'samples' parameter, or length
of MH sampling trajectory. Other notes: Gumbel can be run on multiple
replicates; replicate datasets will be automatically merged.

|

.. rst-class:: transit_sectionend
----

griffin
-------

This is an earlier version of the Gumbel method that
identifies essential genes based on how unlikely 'gaps'
(or consecutive runs of TA sites with 0 insertions) are,
given the overall level of saturation.
It is a frequentist (non-Bayesian) model that uses
the Gumbel Extreme-Value Distribution as a likelihood function.
This is the analysis used in our paper on
`cholesterol catabolism (Griffin et al., 2011)
<http://www.ncbi.nlm.nih.gov/pubmed/21980284>`_.
All things considered, you are probably better off using the
hierarchical-Bayesian Gumbel model above, which does a better job
of estimating internal parameters.

|


.. rst-class:: transit_sectionend
----

.. _`tn5gaps`:

Tn5Gaps
--------

The Tn5Gaps method can be used to determine which genes are essential
in a single condition for **Tn5** datasets. It does an analysis of the
insertions at each site within the genome, makes a call for a given
gene based on the length of the most heavily overlapping run of sites
without insertions (gaps), calculates the probability of this using a
the Gumbel distribution.

.. NOTE::
   Intended only for **Tn5** datasets.



|

How does it work?
~~~~~~~~~~~~~~~~~

This method is loosely is based on the original gumbel analysis
method described in this paper:

Griffin, J.E., Gawronski, J.D., DeJesus, M.A., Ioerger, T.R., Akerley, B.J., Sassetti, C.M. (2011).
`High-resolution phenotypic profiling defines genes essential for mycobacterial survival and cholesterol catabolism. <http://www.ncbi.nlm.nih.gov/pubmed/21980284>`_  *PLoS Pathogens*, 7(9):e1002251.


The Tn5Gaps method modifies the original method in order to work on
Tn5 datasets, which have significantly lower saturation of insertion sites
than Himar1 datasets. The main difference comes from the fact that
the runs of non-insertion (or "gaps") are analyzed throughout the whole
genome, including non-coding regions, instead of within single genes.
In doing so, the expected maximum run length is calculated and a
p-value can be derived for every run. A gene is then classified by
using the p-value of the run with the largest number of nucleotides
overlapping with the gene.

This method was tested on a Salmonella Tn5 dataset presented in this paper:

Langridge GC, Phan MD, Turner DJ, Perkins TT, Parts L, Haase J,
Charles I, Maskell DJ, Peters SE, Dougan G, Wain J, Parkhill J, Turner
AK. (2009). `Simultaneous assay of every Salmonella Typhi gene using one million
transposon mutants. <http://www.ncbi.nlm.nih.gov/pubmed/19826075>`_ *Genome Res.* , 19(12):2308-16.

This data was downloaded from SRA (located `herei <http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP000051>`_) , and used to make
wig files (`base <http://orca1.tamu.edu/essentiality/transit/data/salmonella_base.wig>`_ and `bile <http://orca1.tamu.edu/essentiality/transit/data/salmonella_bile.wig>`_) and the following 4 baseline datasets
were merged to make a wig file: (IL2_2122_1,3,6,8). Our analysis
produced 415 genes with adjusted p-values less than 0.05, indicating
essentiality, and the analysis from the above paper produced 356
essential genes. Of these 356 essential genes, 344 overlap with the
output of our analysis.

|

Example
~~~~~~~
::

    python PATH/src/transit.py tn5gaps genomes/Salmonella-Ty2.prot_table data/salmonella_2122_rep1.wig,data/salmonella_2122_rep2.wig results/test_console_tn5gaps.dat -m 2 -r Sum

        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum

Parameters
~~~~~~~~~~


+ **Minimum Read:** The minimum read count that is considered a true read. Because the Gumbel method depends on determining gaps of TA sites lacking insertions, it may be suceptible to spurious reads (e.g. errors). The default value of 1 will consider all reads as true reads. A value of 2, for example, will ignore read counts of 1.


+ **Replicates:** Determines how to deal with replicates by averaging the read-counts or suming read counts accross datasets. This should not have an affect for the Gumbel method, aside from potentially affecting spurious reads.



|

Outputs and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~

The Tn5Gaps method generates a tab-separated output file at the
location chosen by the user. This file will automatically be loaded
into the Results Files section of the GUI, allowing you to display it
as a table. Alternatively, the file can be opened in a spreadsheet
software like Excel as a tab-separated file. The columns of the output
file are defined as follows:


+-----------------+--------------------------------------------------------------------------------------------------+
| Column Header   | Column Definition                                                                                |
+=================+==================================================================================================+
| ORF             | Gene ID.                                                                                         |
+-----------------+--------------------------------------------------------------------------------------------------+
| Name            | Name of the gene.                                                                                |
+-----------------+--------------------------------------------------------------------------------------------------+
| Desc            | Gene description.                                                                                |
+-----------------+--------------------------------------------------------------------------------------------------+
| k               | Number of Transposon Insertions Observed within the ORF.                                         |
+-----------------+--------------------------------------------------------------------------------------------------+
| n               | Total Number of TA dinucleotides within the ORF.                                                 |
+-----------------+--------------------------------------------------------------------------------------------------+
| r               | Length of the Maximum Run of Non-Insertions observed.                                            |
+-----------------+--------------------------------------------------------------------------------------------------+
| ovr             | The number of nucleotides in the overlap with the longest run partially covering the gene.       |
+-----------------+--------------------------------------------------------------------------------------------------+
| lenovr          | The length of the above run with the largest overlap with the gene.                              |
+-----------------+--------------------------------------------------------------------------------------------------+
| pval            | P-value calculated by the permutation test.                                                      |
+-----------------+--------------------------------------------------------------------------------------------------+
| padj            | Adjusted p-value controlling for the FDR (Benjamini-Hochberg).                                   |
+-----------------+--------------------------------------------------------------------------------------------------+
| call            | Essentiality call for the gene. Depends on FDR corrected thresholds. Essential or Non-Essential. |
+-----------------+--------------------------------------------------------------------------------------------------+

|

Run-time
~~~~~~~~
The Tn5Gaps method takes on the order of 10 minutes.
Other notes: Tn5Gaps can be run on multiple replicates; replicate
datasets will be automatically merged.

|


.. rst-class:: transit_sectionend
----

.. _HMM:

HMM
---

The HMM method can be used to determine the essentiality of the entire genome, as opposed to gene-level analysis of the other methods. It is capable of identifying regions that have unusually high or unusually low read counts (i.e. growth advantage or growth defect regions), in addition to the more common categories of essential and non-essential.

.. NOTE::
   Intended only for **Himar1** datasets.

|

How does it work?
~~~~~~~~~~~~~~~~~

| For a formal description of how this method works, see our paper [DeJesus2013HMM]_:
|
|  DeJesus, M.A., Ioerger, T.R. `A Hidden Markov Model for identifying essential and growth-defect regions in bacterial genomes from transposon insertion sequencing data. <http://www.ncbi.nlm.nih.gov/pubmed/24103077>`_ *BMC Bioinformatics.* 2013. 14:303

|


Example
~~~~~~~

::


  python transit.py hmm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file>
        Optional Arguments:
            -r <string>     :=  How to handle replicates. Sum, Mean. Default: -r Mean
            -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Off.
            -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
            -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0


Parameters
~~~~~~~~~~

The HMM method automatically estimates the necessary statistical
parameters from the datasets. You can change how the method handles
replicate datasets:

-  **Replicates:** Determines how the HMM deals with replicate datasets
   by either averaging the read-counts or summing read counts across
   datasets. For regular datasets (i.e. mean-read count > 100) the
   recommended setting is to average read-counts together. For sparse
   datasets, it summing read-counts may produce more accurate results.

|

Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

| The HMM method outputs two files. The first file provides the most
  likely assignment of states for all the TA sites in the genome. Sites
  can belong to one of the following states: "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage). In
  addition, the output includes the probability of the particular site
  belonging to the given state. The columns of this file are defined as
  follows:

+------------+-----------------------------------------------------------------------------------------------------+
| Column #   | Column Definition                                                                                   |
+============+=====================================================================================================+
| 1          | Coordinate of TA site                                                                               |
+------------+-----------------------------------------------------------------------------------------------------+
| 2          | Observed Read Counts                                                                                |
+------------+-----------------------------------------------------------------------------------------------------+
| 3          | Probability for ES state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 4          | Probability for GD state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 5          | Probability for NE state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 6          | Probability for GA state                                                                            |
+------------+-----------------------------------------------------------------------------------------------------+
| 7          | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+------------+-----------------------------------------------------------------------------------------------------+
| 8          | Gene(s) that share(s) the TA site.                                                                  |
+------------+-----------------------------------------------------------------------------------------------------+

|
|  The second file provides a gene-level classification for all the
  genes in the genome. Genes are classified as "E" (Essential), "GD"
  (Growth-Defect), "NE" (Non-Essential), or "GA" (Growth-Advantage)
  depending on the number of sites within the gene that belong to those
  states.

+-------------------+-----------------------------------------------------------------------------------------------------+
| Column Header     | Column Definition                                                                                   |
+===================+=====================================================================================================+
| Orf               | Gene ID                                                                                             |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Name              | Gene Name                                                                                           |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Desc              | Gene Description                                                                                    |
+-------------------+-----------------------------------------------------------------------------------------------------+
| N                 | Number of TA sites                                                                                  |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n0                | Number of sites labeled ES (Essential)                                                              |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n1                | Number of sites labeled GD (Growth-Defect)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n2                | Number of sites labeled NE (Non-Essential)                                                          |
+-------------------+-----------------------------------------------------------------------------------------------------+
| n3                | Number of sites labeled GA (Growth-Advantage)                                                       |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Avg. Insertions   | Mean insertion rate within the gene                                                                 |
+-------------------+-----------------------------------------------------------------------------------------------------+
| Avg. Reads        | Mean read count within the gene                                                                     |
+-------------------+-----------------------------------------------------------------------------------------------------+
| State Call        | State Classification (ES = Essential, GD = Growth Defect, NE = Non-Essential, GA = Growth-Defect)   |
+-------------------+-----------------------------------------------------------------------------------------------------+

|
|  Note: Libraries that are too sparse (e.g. < 30%) or which contain
  very low read-counts may be problematic for the HMM method, causing it
  to label too many Growth-Defect genes.

|

Run-time
~~~~~~~~

| The HMM method takes less than 10 minutes to complete. The parameters
  of the method should not affect the running-time.

|

.. rst-class:: transit_sectionend
----

.. _resampling:

Resampling
----------

The resampling method is a comparative analysis the allows that can be
used to determine conditional essentiality of genes. It is based on a
permutation test, and is capable of determining read-counts that are
significantly different across conditions.

See :ref:`Pathway Enrichment Analysis <GSEA>` for post-processing the hits to
determine if the hits are associated with a particular functional catogory
of genes or known biological pathway.


.. NOTE::
   Can be used for both **Himar1** and **Tn5** datasets


|

How does it work?
~~~~~~~~~~~~~~~~~

This technique has yet to be formally published in the context of
differential essentiality analysis. Briefly, the read-counts at each
genes are determined for each replicate of each condition. The total
read-counts in condition A is subtracted from the total read counts at
condition B, to obtain an observed difference in read counts. The TA
sites are then permuted for a given number of "samples". For each one of
these permutations, the difference is read-counts is determined. This
forms a null distribution, from which a p-value is calculated for the
original, observed difference in read-counts.

|


Usage
~~~~~


::

  python transit.py resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive resampling. Default: Turned Off.
        -ez             :=  Exclude rows with zero accross conditions. Default: Turned off
                            (i.e. include rows with zeros).
        -pc             :=  Pseudocounts to be added at each site.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias.
                            Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        --ctrl_lib      :=  String of letters representing library of control files in order
                            e.g. 'AABB'. Default empty. Letters used must also be used in --exp_lib
                            If non-empty, resampling will limit permutations to within-libraries.

        --exp_lib       :=  String of letters representing library of experimental files in order
                            e.g. 'ABAB'. Default empty. Letters used must also be used in --ctrl_lib
                            If non-empty, resampling will limit permutations to within-libraries.


Parameters
~~~~~~~~~~

The resampling method is non-parametric, and therefore does not require
any parameters governing the distributions or the model. The following
parameters are available for the method:

-  **Samples:** The number of samples (permutations) to perform. The
   larger the number of samples, the more resolution the p-values
   calculated will have, at the expense of longer computation time. The
   resampling method runs on 10,000 samples by default.

-  **Output Histograms:**\ Determines whether to output .png images of
   the histograms obtained from resampling the difference in
   read-counts.

-  **Adaptive Resampling:** An optional "adaptive" version of resampling
   which accelerates the calculation by terminating early for genes
   which are likely not significant. This dramatically speeds up the
   computation at the cost of less accurate estimates for those genes
   that terminate early (i.e. deemed not significant). This option is
   OFF by default.

-  **Include Zeros:** Select to include  sites that are zero. This is the
   preferred behavior, however, unselecting this (thus ignoring sites that)
   are zero accross all dataset (i.e. completely empty), is useful for
   decreasing running time (specially for large datasets like Tn5).

-  **Normalization Method:** Determines which normalization method to
   use when comparing datasets. Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences. See the :ref:`Normalization <normalization>` section for a description
   of normalization method available in TRANSIT.

-  **--ctrl_lib, --exp_lib** These are for doing resampling with datasets from multiple libraries, see below.

|


Doing resampling with a combined_wig file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Resampling can also now take a combined_wig_ file as input (containing insertion counts
for multiple sample), along with a samples_metadata_ file
that describes the samples. This mode is indicated with a '-c' flag.
If you want to compare more than two conditions, see :ref:`ZINB <zinb>`.


::

  usage: 

  python transit.py resampling -c <combined_wig> <samples_metadata> <control_condition_name> <experimental_condition_name> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

  example:

  python transit.py resampling -c antibiotic_combined_wig.txt antibiotic_samples_metadata.txt Untreated Isoniazid H37Rv.prot_table results.txt -a


Doing resampling with datasets from different libraries.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In most cases, comparisons are done among samples (replicates) from
the same library evaluated in two different conditions.  But if the
samples themselves come from different libraries, then this could
introduce extra variability, the way resampling is normally done.  To
compensate for this, if you specify which libraries each dataset comes
from, the permutations will be restricted to permuting counts only
among samples within each library.  Statistical significance is still
determined from all the data in the end (by comparing the obversed
difference of means between the two conditions to a null distribution).
Of course, this method makes most sense when you have at least 1 replicate
from each library in each condition.

|


Doing resampling between different strains.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most common case is that resampling is done among replicates all
from the same Tn library, and hence all the datasets (fastq files) are
mapped to the same refence genome.  Occasionally, it is useful to
compare TnSeq datasets between two different strains, such as a
reference strain and a clinical isolate from a different lineage.
Suppose for simplicity that you want to compare one replicate from
strain A (e.g. H37Rv) and one replicate from strain B (e.g. CDC1551).
Resampling was not originally designed to handle this case.  The
problem is that the TA sites in the .wig files with insertion counts
might have different coordinates (because of shifts due to indels
between the genomes).  Furthermore, a given gene might not even have
the same number of TA sites in the two strains (due to SNPs).  A
simplistic solution is to just map both datasets to the same genome
sequence (say H37Rv, for example).  Then a resampling comparison could
be run as usual, because the TA sites would all be on the same
coordinate system. This is not ideal, however, because some reads of
strain B might not map properly to genome A due to SNPs or indels
between the genomes.  In fact, in more divergent organisms with higher
genetic diversity, this can cause entire regions to look artificially
essential, because reads fail to map in genes with a large number of
SNPs, resulting in the apparent absence of transposon insertions.

A better approach is to map each library to the custom genome sequence
of its own strain (using TPP).  It turns out the resampling can still
be applied (since it is fundamentally a test on the difference of the
*mean* insertion count in each gene).  The key to making this work,
aside from mapping each library to its own genome sequence, is that
you need an annotation (prot_table) for the second strain that has
been "adapted" from the first strain.  This is because,
to do a comparison between conditions for a gene, Transit needs to be
able to determine which TA sites fall in that gene for each strain.
This can be achieved by producing a "modified" prot_table, where the
START and END coordinates of each ORF in strain B have been adjusted
according to an alignment between genome A and genome B.  You can use
this web app: `Prot_table Adjustment Tool
<http://saclab.tamu.edu/cgi-bin/iutils/app.cgi/>`__, to create a
modifed prot_table, given the prot_table for one strain and the fasta
files for both genomes (which will be aligned).  In other words, the
app allows you to create 'B.prot_table' from 'A.prot_table' (and 'A.fna'
and 'B.fna').

Once you have created B.prot_table, all you need to do is provide
*both* prot_tables to resampling (either through the GUI, or on the
command-line), as a comma-separated list.  For example:

::

  > python transit.py resampling Rv_1_H37Rv.wig,Rv_2_H37Rv.wig 632_1_632WGS.wig,632_2_632WGS.wig H37Rv.prot_table,632WGS.prot_table resampling_output.txt -a

In this example, 2 replicates from H37Rv (which had been mapped to
H37Rv.fna by TPP) were compared to 2 replicates from strain 632 (which
had been mapped to 632WGS.fna, the custom genome seq for strain 632).
The important point is that **two annotations** are given in the 3rd
arg on the command-line: **H37Rv.prot_table,632WGS.prot_table**.  The
assumption is that the ORF boundaries for H37Rv will be used to find
TA sites in Rv_1_H37Rv.wig and Rv_2_H37Rv.wig, and the ORF boundaries
in 632WGS.prot_table (which had been adapted from H37Rv.prot_table
using the web app above) will be used to find TA sites in the
corrsponding regions in 632_1_632WGS.wig and 632_2_632WGS.wig.


Note that, in contrast to handling datasets from different libraries
disucssed above, in this case, the assumption is that all replicates
in condition A will be from one library (and one strain), and all
replicates in condition B will be from another library (another strain).


|

Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

The resampling method outputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+-----------------+-----------------------------------------------------------------+
| Column Header   | Column Definition                                               |
+=================+=================================================================+
| Orf             | Gene ID.                                                        |
+-----------------+-----------------------------------------------------------------+
| Name            | Name of the gene.                                               |
+-----------------+-----------------------------------------------------------------+
| Description     | Gene description.                                               |
+-----------------+-----------------------------------------------------------------+
| N               | Number of TA sites in the gene.                                 |
+-----------------+-----------------------------------------------------------------+
| TAs Hit         | Number of TA sites with at least one insertion.                 |
+-----------------+-----------------------------------------------------------------+
| Sum Rd 1        | Sum of read counts in condition 1.                              |
+-----------------+-----------------------------------------------------------------+
| Sum Rd 2        | Sum of read counts in condition 2.                              |
+-----------------+-----------------------------------------------------------------+
| Delta Rd        | Difference in the sum of read counts.                           |
+-----------------+-----------------------------------------------------------------+
| p-value         | P-value calculated by the permutation test.                     |
+-----------------+-----------------------------------------------------------------+
| p-adj.          | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+-----------------+-----------------------------------------------------------------+

|

Run-time
~~~~~~~~

A typical run of the resampling method with 10,000 samples will take
around 45 minutes (with the histogram option ON). Using the *adaptive
resampling* option (-a), the run-time is reduced to around 10 minutes.

|

.. rst-class:: transit_sectionend
----

Mann-Whitney U-test (utest)
---------------------------

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

.. _genetic-interactions:

Genetic Interactions
--------------------

The genetic interactions (GI) method is a comparative analysis used
used to determine genetic interactions. It is a Bayesian method
that estimates the distribution of log fold-changes (logFC) in two
strain backgrounds under different conditions, and identifies significantly
large changes in enrichment (delta-logFC) to identify those genes
that imply a genetic interaction.


.. NOTE::
   Can be used for both **Himar1** and **Tn5** datasets


|

How does it work?
~~~~~~~~~~~~~~~~~

| For a formal description of how this method works, see our paper [DeJesus20170NAR]_:
|
|  DeJesus, M.A., Nambi, S., Smith, C.M., Baker, R.E., Sassetti, C.M., Ioerger, T.R. `Statistical analysis of genetic interactions in Tn-Seq data. <https://www.ncbi.nlm.nih.gov/pubmed/28334803>`_ *Nucleic Acids Research.* 2017. 45(11):e93. doi: 10.1093/nar/gkx128.

|


Example
~~~~~~~

::

  python transit.py GI <comma-separated .wig control files condition A> <comma-separated .wig control files condition B> <comma-separated .wig experimental files condition A> <comma-separated .wig experimental files condition B> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        --rope <float>  :=  Region of Practical Equivalence. Area around 0 (i.e. 0 +/- ROPE) that is NOT of interest. Can be thought of similar to the area of the null-hypothesis. Default: --rope 0.5
        -n <string>     :=  Normalization method. Default: -n TTR
        -iz             :=  Include rows with zero accross conditions.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0

You can think of 'control' and 'experimental' samples as 'untreated' vs. 'treated'.

Parameters
~~~~~~~~~~

The resampling method is non-parametric, and therefore does not require
any parameters governing the distributions or the model. The following
parameters are available for the method:



-  **Samples:** The number of samples (permutations) to perform. The
   larger the number of samples, the more resolution the p-values
   calculated will have, at the expense of longer computation time. The
   resampling method runs on 10,000 samples by default.


-  **ROPE:** Region of Practical Equivalence. This region defines an area
   around 0.0 that represents differences in the log fold-change that are
   practically equivalent to zero. This aids in ignoring spurious changes
   in the logFC that would otherwise be identified under a strict
   null-hypothesis of no difference.

-  **Include Zeros:** Select to include  sites that are zero. This is the
   preferred behavior, however, unselecting this (thus ignoring sites that)
   are zero accross all dataset (i.e. completely empty), is useful for
   decreasing running time (specially for large datasets like Tn5).

-  **Normalization Method:** Determines which normalization method to
   use when comparing datasets. Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences. See the :ref:`Normalization <normalization>` section for a description
   of normalization method available in TRANSIT.





Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

The GI method outputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+-----------------------------------------+----------------------------------------------------+
| Column Header                           | Column Definition                                  |
+=========================================+====================================================+
| Orf                                     | Gene ID.                                           |
+-----------------------------------------+----------------------------------------------------+
| Name                                    | Name of the gene.                                  |
+-----------------------------------------+----------------------------------------------------+
| Number of TA Sites                      | Number of TA sites in the gene.                    |
+-----------------------------------------+----------------------------------------------------+
| Mean count (Strain A Condition 1)       | Mean read count in strain A, condition 1           |
+-----------------------------------------+----------------------------------------------------+
| Mean count (Strain A Condition 2)       | Mean read count in strain A, condition 2           |
+-----------------------------------------+----------------------------------------------------+
| Mean count (Strain B Condition 1)       | Mean read count in strain B, condition 1           |
+-----------------------------------------+----------------------------------------------------+
| Mean count (Strain B Condition 2)       | Mean read count in strain B, condition 2           |
+-----------------------------------------+----------------------------------------------------+
| Mean logFC (Strain A)                   | The log2 fold-change in read-count for strain A    |
+-----------------------------------------+----------------------------------------------------+
| Mean logFC (Strain B)                   | The log2 fold-change in read-count for strain B    |
+-----------------------------------------+----------------------------------------------------+
| Mean delta logFC                        | The difference in log2 fold-change between B and A |
+-----------------------------------------+----------------------------------------------------+
| Lower Bound delta logFC                 | Lower bound of the difference (delta logFC)        |
+-----------------------------------------+----------------------------------------------------+
| Upper Bound delta logFC                 | Upper bound of the difference (delta logFC)        |
+-----------------------------------------+----------------------------------------------------+
| Prob. of delta-logFC being within ROPE  | Portion of the delta-logFC within ROPE             |
+-----------------------------------------+----------------------------------------------------+
| Adjusted Probability (BFDR)             | Posterior probability adjusted for comparisons     |
+-----------------------------------------+----------------------------------------------------+
| Is HDI outside ROPE?                    | True/False whether the delta-logFC overlaps ROPE   |
+-----------------------------------------+----------------------------------------------------+
| Type of Interaction                     | Final classification.                              |
+-----------------------------------------+----------------------------------------------------+

|


.. rst-class:: transit_sectionend
----

.. rst-class:: transit_clionly

.. _anova:

ANOVA
-----

The Anova (Analysis of variance) method is used to determine which genes
exhibit statistically significant variability of insertion counts across multiple conditions.
Unlike other methods which take a comma-separated list of wig files as input,
the method takes a *combined_wig* file (which combined multiple datasets in one file)
and a *samples_metadata* file (which describes which samples/replicates belong
to which experimental conditions).

|

How does it work?
~~~~~~~~~~~~~~~~~

The method performs the `One-way anova test <https://en.wikipedia.org/wiki/Analysis_of_variance?oldformat=true#The_F-test>`_ for each gene across conditions.
It takes into account variability of normalized transposon insertion counts among TA sites
and among replicates,
to determine if the differences among the mean counts for each condition are significant.


Example
~~~~~~~

::

  python transit.py anova <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        --ignore-conditions <cond1,cond2> :=  Comma separated list of conditions to ignore, for the analysis. Default --ignore-conditions Unknown
        --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.
        -iN <float>     :=  Ignore TAs occuring within given percentage of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring within given percentage of the C terminus. Default: -iC 0.0

The output file generated by ZINB identifies which genes exhibit statistically
significant variability in counts across conditions (see Output and Diagnostics below).

Note: the combined_wig input file can be generated from multiple wig
files through the Transit GUI
(File->Export->Selected_Datasets->Combined_wig), or via the 'export'
command on the command-line (see combined_wig_).

Format of the samples metadata file: a tab-separated file (which you can edit in Excel)
with 3 columns: Id, Condition, and Filename (it must have these headers).  You can include
other columns of info, but do not include additional rows.  Individual rows can be
commented out by prefixing them with a '#'.  Here is an example of a samples metadata file:
The filenames should match what is shown in the header of the combined_wig (including pathnames, if present).

::

  ID      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig



Parameters
~~~~~~~~~~

The following parameters are available for the method:

-  **Ignore Conditions:** Ignores the given set of conditions from the Anova test.

-  **Normalization Method:** Determines which normalization method to
   use when comparing datasets. Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences. See the :ref:`Normalization <normalization>` section for a description
   of normalization method available in TRANSIT.


Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

The anova method outputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+-----------------+-----------------------------------------------------------------+
| Column Header   | Column Definition                                               |
+=================+=================================================================+
| Orf             | Gene ID.                                                        |
+-----------------+-----------------------------------------------------------------+
| Name            | Name of the gene.                                               |
+-----------------+-----------------------------------------------------------------+
| TAs             | Number of TA sites in Gene                                      |
+-----------------+-----------------------------------------------------------------+
| <Condition Mean>| Mean readcounts for the Gene, by condition                      |
+-----------------+-----------------------------------------------------------------+
| p-value         | P-value calculated by the Anova test.                           |
+-----------------+-----------------------------------------------------------------+
| p-adj           | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+-----------------+-----------------------------------------------------------------+
| status          | Debug information (If any)                                      |
+-----------------+-----------------------------------------------------------------+

|

Run-time
~~~~~~~~

A typical run of the anova method takes less than 1 minute for a combined wig file with 6 conditions, 3 replicates per condition.

|


.. rst-class:: transit_sectionend
----

.. rst-class:: transit_clionly

.. _zinb:

ZINB
----

The ZINB (Zero-Inflated Negative Binomial) method is used to determine
which genes exhibit statistically significant variability in either
the magnitude of insertion counts or local saturation, across multiple
conditions.  Like :ref:`ANOVA <anova>`, the ZINB method takes a
*combined_wig* file (which combines multiple datasets in one file) and
a *samples_metadata* file (which describes which samples/replicates
belong to which experimental conditions).

ZINB can be applied to two or more conditions at a time.  Thus it
subsumes :ref:`resampling <resampling>`.  Our testing suggests that
ZINB typically identifies 10-20% more varying genes than resampling
(and vastly out-performs ANOVA for detecting significant variability
across conditions).  Furthermore, because of how ZINB treats magnitude
of read counts separately from local saturation in a gene, it
occasionally identifies genes with variability not detectable by
resampling analysis.

Note: ZINB analysis requires R (statistical analysis software)
to be installed on your system.  See :ref:`Installation Instructions <install-zinb>`.

|

How does it work?
~~~~~~~~~~~~~~~~~

| For a formal description of how this method works, see our `paper on bioRxiv <https://www.biorxiv.org/content/10.1101/590281v1>`_.



Example
~~~~~~~

::

  python transit.py zinb <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        --ignore-conditions <cond1,cond2> :=  Comma separated list of conditions to ignore, for the analysis.
        --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.
        -iN <float>     :=  Ignore TAs occuring within given percentage of the N terminus. Default: -iN 5.0
        -iC <float>     :=  Ignore TAs occuring within given percentage of the C terminus. Default: -iC 5.0
        --covars <covar1,covar2>     :=  Comma-separated list of covariates to include, for the analysis.



.. _combined_wig:

Combined wig files
~~~~~~~~~~~~~~~~~~

Transit now supports a new file format called 'combined_wig' which basically
combines multiple wig files into one file (with multiple columns).  This is
used for some of the new analysis methods for larger collections of datasets, like :ref:`Anova <anova>`, :ref:`ZINB <zinb>`.
Combined_wig files can created through the Transit GUI
(File->Export->Selected_Datasets->Combined_wig), or via the command line.
You can specify the normalization method you want to use with a flag.
TTR is the default, but other relevant normalization options would be 'nonorm'
(i.e. preserve raw counts) and 'betageom' (this corrects for skew, but is slow).


::

  > python src/transit.py export combined_wig --help

  usage: python src/transit.py export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file>

  > python ../transit/src/transit.py export combined_wig Rv_1_H37RvRef.wig,Rv_2_H37RvRef.wig,Rv_3_H37RvRef.wig H37Rv.prot_table clinicals_combined_TTR.wig -n TTR



.. _samples_metadata:

Samples Metadata File
~~~~~~~~~~~~~~~~~~~~~

Format of the *samples_metadata* file: a tab-separated file (which you
can edit in Excel) with 3 columns: Id, Condition, and Filename (it
must have these headers).  You can include other columns of info, but
do not include additional rows.  Individual rows can be commented out
by prefixing them with a '#'.  Here is an example of a samples
metadata file: The filenames should match what is shown in the header
of the combined_wig (including pathnames, if present).

::

  ID      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig

Parameters
~~~~~~~~~~

The following parameters are available for the method:

-  **Ignore Conditions:** Ignores the given set of conditions from the ZINB test.
-  **Include Conditions:** Includes the given set of conditions from the ZINB test. Conditions not in this list are ignored.
-  **Normalization Method:** Determines which normalization method to
   use when comparing datasets. Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences. See the :ref:`Normalization <normalization>` section for a description
   of normalization method available in TRANSIT.
- **Covariates:** If additional covariates distinguishing the samples are available, such as library, timepoint, or genotype, they may be incorporated in the test.

Covariates and Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

While ZINB is focus on identifying variability of insertion counts across conditions,
the linear model also allows you to take other variables into account.
There are two types of auxilliary variables: *covariates* and *interactions*.
Covariates are attributes of the individuals samples that could have a systematic
effect on the insertion counts which we would like to account for and subsequently ignore
(like nuissance variables).
Examples include things like batch or library.  These can be provided as extra
columns in the samples metadata file.

Interactions are extra variables for which we want to test their effect on the 
main variable (or condition).  For example, suppose we collect TnSeq data at several
different timepoints (e.g. length of incubation or infection).  If we just test
time as the condition, we will be identifying genes that vary over time (if timepoints
are numeric, think of the model as fitting a 'slope' to the counts).
But suppose we have data both a wild-type and knock-out strain.  Then we might be
interested in genes for which the time-dependent behavior *differs* between the two
strains (think: different 'slopes'). In such a case, we would say strain and time interact.


If covariates distinguishing the samples are available,
such as batch or library, they may be
incorporated in the ZINB model by using the `covars` flag and samples
metadata file. For example, consider the following samples metadata
file, with a column describing the batch information of each
replicate.

::

  ID      Condition    Filename                                     Batch
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig        B1
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig        B2
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig     B1
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig     B2
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig     B2

This information can be included to eliminate variability due to batch by using
the **\\-\\-covars** flag.

::

 python transit.py zinb combined.wig samples.metadata prot.table output.file --covars Batch


Similarly, one or more interaction may be included in the model.
These are specified by the user with the **\\-\\-interactions** flag,
followed by the name of a column in the samples metadata to test as the interaction
with the condition.  If there are multiple interactions, they may be given as a comma-separated list.

To give an example,
consider an experiment where the condition represents
a treatment (e.g. with values 'treated' and 'control'), and we have another column
called Strain (with values 'wild-type' and 'mutant').
If we want to test whether the effect of the treatment (versus control)
differs depending on the strain, we could do this:

::

 python transit.py zinb combined.wig samples.metadata prot.table output.file --interactions Strain
 
In this case, the condition is implicitly assumed to be the column in the samples metadata file
labeled 'Condition'.  If you want to specify a different column to use as the primary condition to 
test (for example, if Treatment were a distinct column), you can use the **\\-\\-condition** flag:

::

 python transit.py zinb combined.wig samples.metadata prot.table output.file --condition Treatment --interactions Strain
 


The difference between how covariates and interactions are handeled in the model 
is discussed below in the section on Statistical Significance.  

Categorical vs Numeric Covariates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, covariates are intended to be treated as categorical
variables, like 'batch' or 'library' or 'medium'.
In other cases, a covariate might be a numeric value, such as
'time' or 'concentration', in which the ordering of values is
relevant.  The ZINB implementation tries to guess the type of each covariate.
If they are strings, they are treated as discrete factors (each with
their own distinct parameter).  If the given covariate can
be parsed as numbers, the model interprets them as real values.  In this
case, the covariate is treated as a linear factor (regressor), and is
incorporated in the model as a single coefficient, capturing the slope or
trend in the insertion counts as the covariate value increases.


Statistical Significance - What the P-values Mean in the ZINB Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Formally, the P-value is from a likelihood ratio test (LRT) between a
condition-dependent ZINB model (:math:`m_1`) and a
condition-independent (null) ZINB model (:math:`m_0`).

.. math::

  2 \ ln \frac{L(m_1)}{L(m_0)} \sim \chi^2_{df}

where L(.) is the ZINB likelihood function, and :math:`\chi^2_{df}` is
the chi-squared distribution with degrees of freedom (df) equal to
difference in the number of parameters bewteen the two models.  The p-value is
calculated based on this distribution.

In a simple case where variability across a set of conditions X is being tested,
you can think of the model approximately as:


.. math::

  m_1: ln \ \mu = \alpha_0+\vec\alpha X

where :math:`\mu` is an estimate of the mean (non-zero) insertion
count in a gene (a parameter in the likelihood function for ZINB), 
:math:`\alpha_0` is a constant (the mean across all
conditions), and :math:`\vec\alpha` is a vector of coefficients
representing the deviation of the mean count in each condition.
(There is a corresponding equation for estimating the saturation as a
function of condition.)

To evaluate whether the variability across conditions is significant, we
compare to a null model, where the counts are estimated by the global mean only
(dropping the condition variable X).

.. math::

  m_0: ln \ \mu = \alpha_0

When a covariate C is available, it is incorporated in both models (additively),
to account for the effect of the covariate in :math:`m_1`. Coefficients in :math:`\vec\beta`
represent systematic effects on the mean count due to the covariate, and effectively
get subtracted out of the condition coefficients, but :math:`\vec\beta` is also
included in the null model :math:`m_0`, since we want to discount the effect of C on the
likelihood and focus on evaluting the effect of X.


.. math::

  m_1: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta C

  m_0: ln \ \mu = \alpha_0 + \vec\beta C


When an interaction I is being tested, it is incorporated *multiplicatively* in
the main model :math:`m_1` and *additively* in the null model :math:`m_0`:

.. math::

  m_1: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta I + \vec\gamma X*I

  m_0: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta I

The meaning of this is that the coefficients :math:`\vec\alpha` and
:math:`\vec\beta` capture the additive effects of how the mean
insertion count in a gene depends on the condition variable and the
interaction variable, respectively, and the X*I term captures
additional (non-additive) deviations (which is the traditional way
interactions are handled in generalized linear models, GLMs).  Thus,
if there were no interaction, one would expect the mean in datasets
representing the *combination* of X and I to be predicted by the
offsets for each independently.  To the extend that this is not the
case, we say that X and I interaction, and the coefficients
:math:`\gamma` for X*I capture these deviations (non-additive
effects).

For example, think of condition X as Strain (e.g. wild-type vs mutant),
and interaction I as Treatment (e.g. treated versus control).
Then the main model would look like this:

.. math::

  m_1: ln \ \mu = \alpha_0 + \alpha_1 WT  + \alpha_2 mutant + \beta_1 control + \beta_2 treated + \gamma mutant * treated

and this would be compared to the following null model (without the interaction term):

.. math::

  m_0: ln \ \mu = \alpha_0 + \alpha_1 WT  + \alpha_2 mutant + \beta_1 control + \beta_2 treated




Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

The ZINB method outputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+---------------------+-----------------------------------------------------------------+
| Column Header       | Column Definition                                               |
+=====================+=================================================================+
| Orf                 | Gene ID.                                                        |
+---------------------+-----------------------------------------------------------------+
| Name                | Name of the gene.                                               |
+---------------------+-----------------------------------------------------------------+
| TAs                 | Number of TA sites in Gene                                      |
+---------------------+-----------------------------------------------------------------+
| <Condition Mean>    | Mean read-counts for the gene, by condition                     |
+---------------------+-----------------------------------------------------------------+
| <Condition NZMean>  | Non-zero Mean read-counts for the gene, by condition            |
+---------------------+-----------------------------------------------------------------+
| p-value             | P-value calculated by the ZINB test.                            |
+---------------------+-----------------------------------------------------------------+
| p-adj               | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+---------------------+-----------------------------------------------------------------+
| status              | Diagnositic information (explanation for genes not analyzed)    |
+---------------------+-----------------------------------------------------------------+

|

Run-time
~~~~~~~~

A typical run of the ZINB method takes ~5 minutes to analze a combined wig
file with 6 conditions, 3 replicates per condition. It will, of
course, run more slowly if you have many more conditions.

|


.. rst-class:: transit_sectionend
----

.. _normalization:

Normalization
-------------


Proper normalization is important as it ensures that other sources of
variability are not mistakenly treated as real differences in
datasets. TRANSIT provides various normalization methods, which are
briefly described below:

- **TTR:**
    Trimmed Total Reads (TTR), normalized by the total
    read-counts (like totreads), but trims top and bottom 5% of
    read-counts. **This is the recommended normalization method for most cases**
    as it has the beneffit of normalizing for difference in
    saturation in the context of resampling.

- **nzmean:**
    Normalizes datasets to have the same mean over the
    non-zero sites.

- **totreads:**
    Normalizes datasets by total read-counts, and scales
    them to have the same mean over all counts.

- **zinfnb:**
    Fits a zero-inflated negative binomial model, and then
    divides read-counts by the mean. The zero-inflated negative
    binomial model will treat some empty sites as belonging to the
    "true" negative binomial distribution responsible for read-counts
    while treating the others as "essential" (and thus not influencing
    its parameters).

- **quantile:**
    Normalizes datasets using the quantile normalization
    method described by `Bolstad et al.
    (2003) <http://www.ncbi.nlm.nih.gov/pubmed/12538238>`_. In this
    normalization procedure, datasets are sorted, an empirical
    distribution is estimated as the mean across the sorted datasets
    at each site, and then the original (unsorted) datasets are
    assigned values from the empirical distribution based on their
    quantiles.

- **betageom:**
    Normalizes the datasets to fit an "ideal" Geometric
    distribution with a variable probability parameter *p*. Specially
    useful for datasets that contain a large skew. See :ref:`BGC` .

- **nonorm:**
    No normalization is performed.

Command-line
~~~~~~~~~~~~

In addition to choosing normalization for various analyses in the GUI,
you can also call Transit to normalize wig files from the command-line,
as shown in this example:

Example
~~~~~~~

::

  > python src/transit.py normalize --help

  usage: python src/transit.py normalize <input.wig> <output.wig> [-n TTR|betageom]
     or: python src/transit.py normalize -c <combined_wig> <output> [-n TTR|betageom]

  > python src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_TTR.wig -n TTR

  > python src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_BG.wig -n betageom

The normalize command now also works on combined_wig_ files too.
If the input file is a combined_wig file, indicate it with a '-c' flag.



.. rst-class:: transit_sectionend
----

.. _GSEA:

.. rst-class:: transit_clionly
Pathway Enrichment Analysis
---------------------------

How does it work?
~~~~~~~~~~~~~~~~~
Pathway Enrichment Analysis provides a method to
identify enrichment of functionally-related genes among those that are
conditionally essential (i.e.
significantly more or less essential between two conditions).
The analysis is typically applied as post-processing step to the hits identified
by a comparative analysis, such as resampling.
Four analytical method are provided:
a Hypergeometric approach, GSEA by `Subramanian et al (2005) <https://www.ncbi.nlm.nih.gov/pubmed/16199517>`_, and GSEA-Z, GSEA-Chi proposed by `Irizarry et al. (2009) <https://www.ncbi.nlm.nih.gov/pubmed/20048385>`_.
For the Hypergeometic test, genes in the resampling output file with adjusted p-value < 0.05 are taken as hits,
and evaluated for overlap with functional categories of genes.
The GSEA methods use the whole list of genes, ranked in order of statistical significance
(without requiring a cutoff), to calculated enrichment.

Two systems of categories are provided for *M. tuberculosis* (but you can add your own):
the Sanger functional categories of genes determined by Stewart Cole in the
original annotation of the H37Rv genome (1998, with updates), and
also GO terms (Gene Ontology).  They are in the src/pytransit/data/ directory.

For now, pathway enrichment analysis is only implemented as a command-line function,
and is not available in the Transit GUI.

GSEA is very slow, which is why we included
two faster methods that use approximations (GSEA-Z and GSEA-Chi).


Example
~~~~~~~

::

    python ../../transit.py pathway_enrichment <resampling files> <annotation file> <output file> [-p] [-S] [-M GSEA|HYPE|GSEA-Z|GSEA-CHI]

|


Parameters
~~~~~~~~~~
- **Resampling File**
    The resampling file is the one obtained after using the resampling method in Transit. (It is a tab separated file with 11 columns.) GSEA method makes usage of the last column (adjusted P-value)
- **Annotation File**
    This file is a tab-separated text file (which you could open in Excel).  The file format has 3 columns: the id of the pathway (e.g. "GO:0090502"), a description (e.g. "RNA phosphodiester bond hydrolysis"), and a space-separated list of genes (ORF ids, like Rv1339).  Examples for Sanger categories and GO terms for H37Rv are in src/pytransit/data/.
- **Output File**
    This parameter is used to specify the output file name and the path where it will be created.
- **p**
   This parameter is optional, and can be used to adjust the stringency of the GSEA calculation. GSEA method calculates a weighted Kolgomorov-Smirnov statistics. The default value is 1, when p=0, the enrichment score is the Kolgomorov-Smirnov statistics.
- **S**
   Number of samples for generating the null distribution.  In order to estimate the significance, the enrichment score is compared to a null distribution computed with S randomly assigned phenotypes. The default S value is 1000.  Changing S affects the speed of the calculation linearly.
- **M**
    Methodology to be used. (HYPE, or hypergeometric test, is the default)
  -**GSEA**
    This method was proposed by Subramanian in:

    Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).  `Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles <http://www.pnas.org/content/102/43/15545.short>`_ . Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
  -**HYPE**
    This is a traditional Hypergeometric approach, which test whether the
proportion of genes with a given function/pathway among the list of hits from resampling
is significantly higher or lower than expected, relative to the background distribution (proportion of the category/pathway among all the ORFS in the genome)
  -**GSEA-Z**
    A faster approximation of GSEA.  This method is the one proposed by Irizarry, R. A et al in :

    Irizarry, R. A., Wang, C., Zhou, Y., & Speed, T. P. (2009). `Gene set enrichment analysis made simple  <http://journals.sagepub.com/doi/abs/10.1177/0962280209351908>`_. Statistical methods in medical research, 18(6), 565-575.
  -**GSEA-CHI**
    A faster approximation of GSEA.  This method is the one proposed by Irizarry, R. A et al in :

    Irizarry, R. A., Wang, C., Zhou, Y., & Speed, T. P. (2009). `Gene set enrichment analysis made simple  <http://journals.sagepub.com/doi/abs/10.1177/0962280209351908>`_. Statistical methods in medical research, 18(6), 565-575.

Run-time
~~~~~~~~
The GSEA method, proposed by Subramanian, might take some hours to calculate the p-value.
The other methods are much faster.


Outputs and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~
The output file is a tab separated file and according to the method, the file has specific columns.

-*Output File*
   -*GSEA*
      - ID descr: ID of the pathway, functional category, or GO term, with its description. This information comes from the annotation file

      - Total genes: The number of genes in the pathway

      - score: Enrichment Score

      - P-Value: Statistical Significance

      - P-Adjust : FDR Correction

      - list of genes in the category, ranked by their log-fold-change (from resampling). Rank values near 0 or near N (total number of genes in genome) are the most interesting (representing pathways with high enrichment)

   -*HYPE*

      - ID descr: ID of the pathway, functional category, or GO term, with its description. This information comes from the annotation file

      - Total genes, The number of genes in the pathway. The genes considered hits are those with p-value < 0.05

      - Total in the intersection, The number of genes that are in the pathway and the whole genome

      - P-Value, the statistical significance

      - P-Adjust, FDR correction of the p-value (using Benjamini-Hochberg)

      - genes in intersection of the pathway and list of hits (conditional essentials)

   -*GSEA-Z, GSEA-CHI*

      - ID descr: ID of the pathway, functional category, or GO term, with its description. This information comes from the annotation file

      - Total genes, The number of genes in the pathway

      - Score, Either Z or Chi, according to the one selected to calculate

      - P-Value, the statistical significance

      - P-Adjust, FDR correction of the p-value

      - Rank of Genes, list of genes in the pathway and their position in the resampling file (LFC order).



.. rst-class:: transit_sectionend
----

.. _tnseq_stats:

tnseq_stats command
~~~~~~~~~~~~~~~~~~~

You can generate the same table to statistics as on the Quality Control panel in the GUI
from the command-line using the 'tnseq_stats' command.  Here is an example:

::

  > python src/transit.py tnseq_stats --help

  usage: python src/transit.py tnseq_stats <file.wig>+ [-o <output_file>]
         python src/transit.py tnseq_stats -c <combined_wig> [-o <output_file>]

  > python src/transit.py tnseq_stats -c src/pytransit/data/cholesterol_glycerol_combined.dat

  dataset density mean_ct NZmean  NZmedian        max_ct  total_cts       skewness        kurtosis
  src/pytransit/data/cholesterol_H37Rv_rep1.wig   0.44    139.6   317.6   147     125355.5        10414005.0      54.8    4237.7
  src/pytransit/data/cholesterol_H37Rv_rep2.wig   0.44    171.4   390.5   148     704662.8        12786637.9      105.8   14216.2
  src/pytransit/data/cholesterol_H37Rv_rep3.wig   0.36    173.8   484.2   171     292294.8        12968502.500000002      42.2    2328.0
  src/pytransit/data/glycerol_H37Rv_rep1.wig      0.42    123.3   294.5   160     8813.3  9195672.4       4.0     33.0
  src/pytransit/data/glycerol_H37Rv_rep2.wig      0.52    123.8   240.1   127     8542.5  9235984.2       4.0     33.5






