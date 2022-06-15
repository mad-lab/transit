.. _gumbel:

======
Gumbel
======

The Gumbel can be used to determine which genes are essential in a
single condition. It does a gene-by-gene analysis of the insertions at
TA sites with each gene, makes a call based on the longest consecutive
sequence of TA sites without insertion in the genes, calculates the
probability of this using a Bayesian model.

.. NOTE::
   Intended only for **Himar1** datasets.

How does it work?
-----------------

|

This method for identifying essential genes is based on analyzing
'gaps', or consecutive sequences of TA sites lacking insertions.
The statistical significance of the length of a gap is determined
using the Gumbel distribution, which is a form of an Extreme-Value distribution.

For a formal description of how this method works, see our paper [DeJesus2013]_:

|
  DeJesus, M.A., Zhang, Y.J., Sassettti, C.M., Rubin, E.J.,
  Sacchettini, J.C., and Ioerger, T.R. (2013).  `Bayesian analysis of
  gene essentiality based on sequencing of transposon insertion
  libraries. <http://www.ncbi.nlm.nih.gov/pubmed/23361328>`_
  *Bioinformatics*, 29(6):695-703.

|
**Update (2021) - Binomial**

Since the Gumbel method depends on the overall
saturation (percent of TA sites with insertions), it can sometimes
call a lot of smaller genes 'Uncertain'.  This might reduce the total
number of essentials detected.  For example, if the saturation is
~30%, no genes with fewer than ~10 TA sites might confidently be
labeled as Essential.  In particular, there are often genes with no
insertions that look like they should obviously be called essential,
and yet, they are too short to be confidently called essential in
low-saturation datasets by the conservative Gumbel model.

To compensate for this, we have added a simple **Binomial** model for
detecting small genes totally lacking insertions (described in `(Choudhery et al, 2021)
<https://journals.asm.org/doi/full/10.1128/mSystems.00876-21>`_).  In
the output file, they are labled as 'EB' to distinguish them from
essentials (E) based on the Gumbel model.  The EB genes supplement
genes marked E by Gumbel, and the combination of the 2 groups (E and
EB) should be considered as 'essentials'.  The number of genes in
each category is reported in a table in the header of the output file, like this:

::

 (for an M. tuberculosis TnSeq dataset with 60% saturation:)
 #Summary of Essentiality Calls:
 #  E  =  551 (essential based on Gumbel)
 #  EB =   94 (essential based on Binomial)
 #  NE = 2829 (non-essential)
 #  U  =  262 (uncertain)
 #  S  =  254 (too short)

 (for an M. tuberculosis TnSeq dataset with 40% saturation:)
 #  E  =  194 (essential based on Gumbel)
 #  EB =  315 (essential based on Binomial)
 #  NE = 2441 (non-essential)
 #  U  =  774 (uncertain)
 #  S  =  266 (too short)



Usage
-----

::

  > python3 transit.py gumbel <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -b <integer>    :=  Number of Burn-in samples. Default -b 500
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -t <integer>    :=  Trims all but every t-th value. Default: -t 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring at given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given percentage (as integer) of the C terminus. Default: -iC 0




Parameters
----------

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
-----------------------

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
--------

The Gumbel method takes on the order of 10 minutes for 10,000 samples.
Run-time is linearly proportional to the 'samples' parameter, or length
of MH sampling trajectory. Other notes: Gumbel can be run on multiple
replicates; replicate datasets will be automatically merged.

|

.. rst-class:: transit_sectionend
----
