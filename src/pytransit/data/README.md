Transit Data directory
======================

Files for Pathway Enrichment analysis
-------------------------------------

These files can be used for pathway_enrichment analysis.
Files are included for *M. tuberculosis* H37Rv and *M. smegmatis* mc2 155.

```
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
```
'\*' means *recommended* combination of method with system of functional categories


These file are to be used as inputs for the transit pathway_enrichment 
command (see the [documentation](https://transit.readthedocs.io/en/latest/transit_methods.html#pathway-enrichment-analysis)).



Test files for resampling
-------------------------

These wig files are from the 
[Griffin et al (2011) paper in PLOS Pathogens](http://www.ncbi.nlm.nih.gov/pubmed/21980284).
on M. tuberculosis genes required for growth on glycerol.

 * glycerol_H37Rv_rep1.wig, glycerol_H37Rv_rep2.wig, glycerol_H37Rv_rep3.wig
 * cholesterol_H37Rv_rep1.wig, cholesterol_H37Rv_rep2.wig, cholesterol_H37Rv_rep3.wig
