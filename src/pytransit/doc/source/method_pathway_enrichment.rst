
.. _GSEA:

.. rst-class:: transit_clionly
Pathway Enrichment Analysis
===========================

Pathway Enrichment Analysis provides a method to
identify enrichment of functionally-related genes among those that are
conditionally essential (i.e.
significantly more or less essential between two conditions).
The analysis is typically applied as post-processing step to the hits identified
by a comparative analysis, such as *resampling*.
Several analytical method are provided:
Fisher's exact test (FET, hypergeometric distribution), GSEA (Gene Set Enrichment Analysis)
by `Subramanian et al (2005) <https://www.ncbi.nlm.nih.gov/pubmed/16199517>`_,
and `Ontologizer <https://www.ncbi.nlm.nih.gov/pubmed/17848398>`_.
For Fisher's exact test,
genes in the resampling output file with adjusted p-value < 0.05 are taken as hits,
and evaluated for overlap with functional categories of genes.
The GSEA methods use the whole list of genes, ranked in order of statistical significance
(without requiring a cutoff), to calculate enrichment.

Three systems of categories are provided for (but you can add your own):
the Sanger functional categories of genes determined in the
original annotation of the H37Rv genome (`Cole et al, 1998 <https://www.ncbi.nlm.nih.gov/pubmed/9634230>`_,
with subsequent updates),
COG categories (`Clusters of Orthologous Genes <https://www.ncbi.nlm.nih.gov/pubmed/25428365>`_) and
also GO terms (Gene Ontology).  The supporting files for *M. tuberculosis*
H37Rv are in the src/pytransit/data/ directory.

For other organisms, it might be possible to download COG categories from
`http://www.ncbi.nlm.nih.gov/COG/ <http://www.ncbi.nlm.nih.gov/COG/>`_
and GO terms from `http://www.geneontology.org <http://www.geneontology.org>`_
or `http://patricbrc.org <http://patricbrc.org>`_.
If these files can be obtained for your organism, they will have to be converted into
the *associations* file format described below. (The *pathways* files for COG categories and GO terms
in the Transit data directory should still work, because they just encode pathways names for all terms/ids.)

At present, pathway enrichment analysis is only implemented as a command-line function,
and is not available in the Transit GUI.


Usage
-----

::

  > python3 ../../transit.py pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FET|GSEA|GO>] [-PC <int>] [-ranking SLPV|LFC] [-p <float>] [-Nperm <int>] [-Pval_col <int>] [-Qval_col <int>]  [-LFC_col <int>]

  Optional parameters:
     -M FET|GSEA|ONT:     method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Analysis (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)
     -Pval_col <int>    : indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2, i.e. second-to-last column)
     -Qval_col <int>    : indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)
 for GSEA...
     -ranking SLPV|LFC  : SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling 
     -LFC_col <int>     : indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)
     -p <float>         : exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)
     -Nperm <int>       : number of permutations to simulate for null distribution to determine p-value (default=10000)
 for FET...
     -PC <int>          :  pseudo-counts to use in calculating p-value based on hypergeometric distribution (default=2)

|


Parameters
----------
- **Resampling File**
    The resampling file is the one obtained after using the resampling method in Transit. (It is a tab separated file with 11 columns.) GSEA method makes usage of the last column (adjusted P-value)
- **Associations File**
   This is a tab-separated text file with 2 columns: pathway id, and pathway name. If a gene is in multiple pathways, the associated ids should be listed on separate lines.  It is OK if there are no associations listed for some genes.  Important: if pathways are hierarchical, you should expand this file to explicitly include associations of each gene with all parent nodes. Files with GO term associations will have to be pre-processed this way too.

::

  Example: H37Rv_sanger_roles.dat

  Rv3823c	II.C.4
  Rv3823c	II.C
  Rv3823c	II
  Rv0337c	I.D.2
  Rv0337c	I.D
  Rv0337c	I
  ...

- **Pathways File**
   This is a tab-separated text file with 2 columns: pathway id, and pathway name.

::

  Example: sanger_roles.dat

  I	Small-molecule metabolism
  I.A	Degradation
  I.A.1	Carbon compounds
  I.A.2	Amino acids and amines
  I.A.3	Fatty acids
  I.A.4	Phosphorous compounds
  ...

- **\-\-Pval_col <int>**: indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2, i.e. second-to-last column)

- **\-\-Qval_col <int>**: indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)


- **-M [FET|GSEA|ONT]**
    Methodology to be used. FET is used by default (even without specifying -M).

  **-M FET**
    This implements Fisher's Exact Test (hypergeometric distribution) to determine a p-value for each pathway, based on the proportion of pathway member observed in list of hits (conditionally essential gene by resampling, padj<0.05) compared to the background proportion in the overall genome, and p-values are adjusted post-hoc by the Benjamini-Hochberg procedure to limit the FDR to 5%.

    In the output file, an "enrichment score" is reported, which is the ratio of the observed number of pathway members among the hits to the expected number.  Pseudocounts of 2 are included in the calculation to reduce the bias toward small pathways with only a few genes; this can be adjusted with the -PC flag (below).

    FET can be used with GO terms.

    Additional flags for FET:

    - **-PC <int>**: Pseudocounts used in calculating the enrichment score and p-value by hypergeometic distribution. Default: PC=2.

  **-M GSEA**
    Gene Set Enrichment Analysis. GSEA assess the significance of a pathway by looking at how the members fall in the ranking of all genes.  The genes are first ranked by significance from resampling.  Specifically, they are sorted by signed-log-p-value, SLPV=sign(LFC)*(log(pval)), which puts them in order so that the most significant genes with negative LFC are at the top, the most significant with positive LFC are at the bottom, and insignificant genes fall in the middle.  Roughly, GSEA computes the mean rank of pathway members, and evaluates significance based on a simulated a null distribution.  p-values are again adjusted at the end by BH.

    `Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).  `ene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences, 102(43), 15545-15550. <http://www.pnas.org/content/102/43/15545.short>`_

    GSEA can be used with GO terms.

    Additional flags for GSEA:

    - **-ranking SLPV|LFC**: method used to rank all genes; SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling

    - **-p <float>**: exponent to use in calculating enrichment score; recommend trying '-p 0' (default) or '-p 1' (as used in Subramaniam et al, 2005)

    - **-Nperm <int>**: number of permutations to simulate for null distribution to determine p-value (default=10000)

    - **\-\-LFC_col <int>**: indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)


  **-M ONT**
    Ontologizer is a specialized method for GO terms that takes parent-child relationships into account among nodes in the GO hierarchy.  This can enhance the specificity of pathways detected as significant.  (The problem is that there are many GO terms in the hierarchy covering similar or identical sets of genes, and often, if one node is significantly enriched, then several of its ancestors will be too, which obscures the results with redundant hits; Ontologizer reduces the significance of nodes if their probability distribution among hits can be explained by their parents.) Hierarhical relationships among GO terms are encoded in an OBO file, which is included in the src/pytransit/data/ directory.

    `Grossmann S, Bauer S, Robinson PN, Vingron M. Improved detection of overrepresentation of Gene-Ontology annotations with parent child analysis. Bioinformatics. 2007 Nov 15;23(22):3024-31. <https://www.ncbi.nlm.nih.gov/pubmed/17848398>`_

  For the ONT method in pathway_enrichment, the enrichment for a given
  GO term can be expressed (in a simplified way, leaving out the
  pseudocounts) as:

::

  enrichment = log (  (b/q) / (m/p)  )
|

  where:

*    b is the number of genes with this GO term in the subset of hits (e.g. conditional essentials from resampling, with qval<0.05)
*    q is the number of genes in the subset of hits with a parent of this GO term
*    m is the total number of genes with this GO term in the genome
*    p is the number of genes in the genome with a parent of this GO term

  So enrichment is the log of the ratio of 2 ratios:

  1. the relative abundance of genes with this GO term compared to those with a parent GO term   among the hits
  2. the relative abundance of genes with this GO term compared to those with a parent GO term   in the whole genome


Auxilliary Pathway Files in Transit Data Directory
--------------------------------------------------

::

These files for pathway analysis are distributed in the Transit data directory
(e.g. transit/src/pytransit/data/).

Note: The "Sanger" roles are custom pathway associations for
*M. tuberculosis* defined in the original Nature paper on
the H37Rv genome sequence `(Cole et al., 1998)
<https://www.nature.com/articles/31159>`_ (Table 1).  They are more specific
that COG categories, but less specific than GO terms.  For other
organisms, one should be able to find GO terms (e.g. on PATRIC,
Uniprot, or geneontology.org) and COG roles (from
https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/, `(Galerpin et al, 2021)
<https://academic.oup.com/nar/article/49/D1/D274/5964069>`_ ).

Pathway association files for *M. smegmatis* mc2 155 are also provided in the table below.


+----------+----------+--------------------+--------------------------------------+------------------------------------+
| system   | num roles| applicable methods | associations of genes with roles     | pathway definitions/role names     |
+==========+==========+====================+======================================+====================================+
| COG      | 25       | FET*, GSEA         | H37Rv_COG_roles.dat;                 | COG_roles.dat                      |
|          |          |                    | smeg_COG_roles.dat                   |                                    |
+----------+----------+--------------------+--------------------------------------+------------------------------------+
| Sanger   | 153      | FET*, GSEA*        | H37Rv_sanger_roles.dat               | sanger_roles.dat                   |
+----------+----------+--------------------+--------------------------------------+------------------------------------+
| GO       | 2545     | ONT*               | H37Rv_GO_terms.txt;                  | gene_ontology.1_2.3-11-18.obo      |
|          |          |                    | smeg_GO_terms.txt                    |                                    |
+----------+----------+--------------------+--------------------------------------+------------------------------------+
|          |          | FET, GSEA          | H37Rv_GO_terms.txt;                  | GO_term_names.dat                  |
|          |          |                    | smeg_GO_terms.txt                    |                                    |
+----------+----------+--------------------+--------------------------------------+------------------------------------+

'\*' means *recommended* combination of method with system of functional categories


Current Recommendations
-----------------------

Here are the recommended combinations of pathway methods to use for different systems of functional categories:

 * For COG, use '-M FET'
 * For Sanger roles, try both FET and GSEA
 * For GO terms, use 'M -ONT'


Examples
--------

::

    # uses Fisher's exact test by default (with PC=2 as pseudocounts)
    > transit pathway_enrichment resampling_glyc_chol.txt $DATA/H37Rv_sanger_roles.dat $DATA/sanger_roles.dat pathways_glyc_chol_Sanger.txt

    # can do this with GO terms too
    > transit pathway_enrichment resampling_glyc_chol.txt $DATA/H37Rv_GO_terms.txt $DATA/GO_term_names.dat pathways_glyc_chol_GO.txt

    # with COG categories
    > transit pathway_enrichment resampling_glyc_chol.txt $DATA/H37Rv_COG_roles.dat $DATA/COG_roles.dat pathways_glyc_chol_COG.txt

    # can also do GSEA method (on any system of functional categories)
    > transit pathway_enrichment resampling_glyc_chol.txt $DATA/H37Rv_sanger_roles.dat $DATA/sanger_roles.dat pathways_Sanger_GSEA.txt -M GSEA

    # Ontologizer is a specialized method for GO terms
    > transit pathway_enrichment resampling_glyc_chol.txt $DATA/H37Rv_GO_terms.txt $DATA/GO_term_names.dat pathways_Ontologizer.txt -M ONT

The $DATA environment variable in these examples refers to the Transit data directory, e.g. src/pytransit/data/.


.. rst-class:: transit_sectionend
------
