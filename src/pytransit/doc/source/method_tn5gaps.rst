.. _`tn5gaps`:

Tn5Gaps
=======

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
-----------------

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

This data was downloaded from SRA (located `here <http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP000051>`_) , and used to make
wig files (`baseline <http://orca1.tamu.edu/essentiality/transit/data/salmonella_baseline.wig>`_ and `bile <http://orca1.tamu.edu/essentiality/transit/data/salmonella_bile.wig>`_) and the following 4 baseline datasets
were merged to make a wig file: (IL2_2122_1,3,6,8). Our analysis
produced 415 genes with adjusted p-values less than 0.05, indicating
essentiality, and the analysis from the above paper produced 356
essential genes. Of these 356 essential genes, 344 overlap with the
output of our analysis.

|

Usage
-----
::

    python3 ../../../transit.py tn5gaps <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0


Parameters
----------


+ **Minimum Read:** The minimum read count that is considered a true read. Because the Gumbel method depends on determining gaps of TA sites lacking insertions, it may be suceptible to spurious reads (e.g. errors). The default value of 1 will consider all reads as true reads. A value of 2, for example, will ignore read counts of 1.


+ **Replicates:** Determines how to deal with replicates by averaging the read-counts or suming read counts accross datasets. This should not have an affect for the Gumbel method, aside from potentially affecting spurious reads.

+ **-iN:** Trimming of insertions in N-terminus (given as percentage of ORF length, e.g. "5" for 5%; default=0)

+ **-iC:** Trimming of insertions in C-terminus (given as percentage of ORF length, e.g. "5" for 5%; default=0)

Example
-------
::

    python3 PATH/src/transit.py tn5gaps salmonella_baseline.wig Salmonella-Ty2.prot_table salmonella_baseline_tn5gaps_trimmed.dat -m 2 -r Sum -iN 5 -iC 5


These input and output files can be downloaded from the **Example Data** section on the `Transit home page <http://saclab.tamu.edu/essentiality/transit/index.html>`_ .

|

Outputs and diagnostics
-----------------------

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
--------
The Tn5Gaps method takes on the order of 10 minutes.
Other notes: Tn5Gaps can be run on multiple replicates; replicate
datasets will be automatically merged.

|


.. rst-class:: transit_sectionend
----
