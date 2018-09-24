

Analysis Methods
================


TRANSIT has analysis methods capable of analyzing **Himar1** and **Tn5** datasets.
Below is a description of some of the methods.

|

Gumbel
------

The Gumbel can be used to determine which genes are essential in a
single condition. It does a gene-by-gene analysis of the insertions at
TA sites with each gene, makes a call based on the longest consecutive
sequence of TA sites without insertion in the genes, calculates the
probability of this using a Bayesian model.

.. NOTE::
   Intended only for **Himar1** datasets.


|

How does it work?
~~~~~~~~~~~~~~~~~

| For a formal description of how this method works, see our paper [DeJesus2013]_:

|

|  DeJesus, M.A., Zhang, Y.J., Sassettti, C.M., Rubin, E.J.,
  Sacchettini, J.C., and Ioerger, T.R. (2013).
| `Bayesian analysis of gene essentiality based on sequencing of transposon insertion libraries. <http://www.ncbi.nlm.nih.gov/pubmed/23361328>`_ *Bioinformatics*, 29(6):695-703.

|

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

This method was tested on a salmonella Tn5 dataset presented in this
paper:

Langridge GC1, Phan MD, Turner DJ, Perkins TT, Parts L, Haase J,
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

--------------

|

.. _resampling:

Re-sampling
-----------

The re-sampling method is a comparative analysis the allows that can be
used to determine conditional essentiality of genes. It is based on a
permutation test, and is capable of determining read-counts that are
significantly different across conditions.


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

Parameters
~~~~~~~~~~

The resampling method is non-parametric, and therefore does not require
any parameters governing the distributions or the model. The following
parameters are available for the method:

-  **Samples:** The number of samples (permutations) to perform. The
   larger the number of samples, the more resolution the p-values
   calculated will have, at the expense of longer computation time. The
   re-sampling method runs on 10,000 samples by default.

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

|

Output and Diagnostics
~~~~~~~~~~~~~~~~~~~~~~

The re-sampling method outputs a tab-delimited file with results for each
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

A typical run of the re-sampling method with 10,000 samples will take
around 45 minutes (with the histogram option ON). Using the adaptive
resampling option, the run-time is reduced to around 10 minutes.



--------------




|

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

Parameters
~~~~~~~~~~

The resampling method is non-parametric, and therefore does not require
any parameters governing the distributions or the model. The following
parameters are available for the method:



-  **Samples:** The number of samples (permutations) to perform. The
   larger the number of samples, the more resolution the p-values
   calculated will have, at the expense of longer computation time. The
   re-sampling method runs on 10,000 samples by default.


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


--------------

|

.. _anova:

ANOVA (command line only)
--------------

The Anova (Analysis of variance) method, is used to determine conditionally essential genes, across multiple conditions. See :ref:`Example <anova_example>` for usage.

--------------

.. NOTE::
   Can be used for both **Himar1** and **Tn5** datasets

|

How does it work?
~~~~~~~~~~~~~~~~~

The method performs the `One-way anova test <https://en.wikipedia.org/wiki/Analysis_of_variance?oldformat=true#The_F-test>`_ for each gene, across conditions.

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

The re-sampling method outputs a tab-delimited file with results for each
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
| p-adj.          | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+-----------------+-----------------------------------------------------------------+

|

Run-time
~~~~~~~~

Takes in the order of 30 seconds for a combined wig file with 10 conditions, with 3 replicates per condition.

|

.. _normalization:

Normalization
-------------


Proper normalization is important as it ensures that other sources of variability are not mistakenly treated
as real differences in datasets. TRANSIT provides various normalization methods, which are briefly described below:

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
    useful for datasets that contain a large skew.

- **nonorm:**
    No normalization is performed.

-----------------------------------



.. _GSEA:

Gene Set Enrichment Analysis (GSEA)
-----------------------------------
How does it work?
~~~~~~~~~~~~~~~~~
Gene Set Enrichment Analysis provides a score to measure a set of functionally related genes are related to a pathway. The programmed analysis here provides 4 methodologies to perform GSEA, GSEA by Subramanian et al, a Hypergeometric approach, and GSEA-Z, GSEA-Chi proposed by Irizarry et al.

Parameters
~~~~~~~~~~
- **Resampling File**
    The resampling file is the one obtained after using the resampling method in Transit. It is a tab separated file with 11 columns. GSEA method makes usage of the first two, Genes and Name, and column 10, P-value.
- **Annotation File**
    This file is a text, the first column represents the pathway to check, followed by a description joined by a dash (-), and the list of genes, separated by a space, in that specific pathway. For example, GO:0090502-RNA phosphodiester bond hydrolysis,Rv0627 Rv1102c Rv1339 Rv1955 Rv2063A Rv2228c Rv2407 Rv2549c Rv2752c Rv2801c Rv2816c Rv2902c Rv2925c Rv3796 Rv3923c
- **Output File**
    This parameter is used to set the output file name and the path where it will be created.
- **p**
   This parameter is optional. GSEA method calculates a weighted Kolgomorov-Smirnov statistics. The default value is 1, when p=0, the enrichment score is the Kolgomorov-Smirnov statistics.
- **S**
   In order to estimate the significance, the enrichment score is compared to a null distribution computed with S randomly assigned phenotypes. The default S value is 1000.
- **M**
    This analysis can be performed using different methodologies
  -**GSEA**
    This method was proposed by Subramanian in:
    
    Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).  `Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles <http://www.pnas.org/content/102/43/15545.short>`_ . Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
  -**HYPE**
    This is a traditional Hypergeometric approach.
  -**GSEA-Z**
    This method is the one proposed by Irizarry, R. A et al in :

    Irizarry, R. A., Wang, C., Zhou, Y., & Speed, T. P. (2009). `Gene set enrichment analysis made simple  <http://journals.sagepub.com/doi/abs/10.1177/0962280209351908>`_. Statistical methods in medical research, 18(6), 565-575.
  -**GSEA-CHI**
    This method is the one proposed by Irizarry, R. A et al in :
        
    Irizarry, R. A., Wang, C., Zhou, Y., & Speed, T. P. (2009). `Gene set enrichment analysis made simple  <http://journals.sagepub.com/doi/abs/10.1177/0962280209351908>`_. Statistical methods in medical research, 18(6), 565-575.

Run-time
~~~~~~~~
GSEA method, proposed by Subramanian, might take some hours to calculate the p-value.


Outputs and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~
The output file is a tab separated file and according to the method, the file has specific columns.

- *Output File*
  -*GSEA*
      -ID descr: ID of the pathway, SAGER Category, GO Term, with its description. This information comes from the annotation file
      -Total genes: The number of genes in the pathway
      -score: Enrichment Score
      -pval: Statistical Significance
      -padj : FDR Correction
      -rank of genes: According to Subramanian, "The goal of GSEA is to determine whether the member of S are randomly distributed throughout L or primarily find at the top or bottom "
  -*HYPE*
      -ID descr, ID of the pathway, SAGER Category, GO Term, with its description. This information comes from the annotation file
      -Total genes, The number of genes in the pathway      
      -Total in the intersection, The number of genes that are in the pathway and the whole genome
      -pval, the statistical significance
      -padj, FDR correction of the p-value
      -genes in intersection
  -*GSEA-Z, GSEA-CHI*
         -ID descr, ID of the pathway, SAGER Category, GO Term, with its description. This information comes from the annotation file
         -Total genes, The number of genes in the pathway      
         -Score, Either Z or Chi, according to the one selected to calculate
         -P-Value, the statistical significance
         -P-Adjust, FDR correction of the p-value

-----------------------------------
