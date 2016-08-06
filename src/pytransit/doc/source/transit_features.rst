

Features
========


TRANSIT has several useful features to help inspect the quality of datasets as
and export them to different formats.

|

Quality Control
---------------

As you add datasets to the control or experimental sections, TRANSIT
automatically provides some metrics like density, average, read-counts and
max read-count to give you an idea of how the quality of the dataset.

However, TRANSIT provides more in-depth statistics in the Quality Control
window. To use this feature, add the annotation file for your organism
(in .prot_table or GFF3 format). Next, add and highlight/select the desired
read-count datasets in .wig format. Finally, click on View -> Quality Control.
This will open up a new window containing a table of metrics for the datasets
as well as figures corresponding to whatever dataset is currently highlighted.

.. image:: _images/transit_quality_control_window.png
   :width: 600
   :align: center


QC Metrics Table
~~~~~~~~~~~~~~~~

The Quality Control window contains a table of the datasets and metrics, similar
to the one in the main TRANSIT interface. This table has an extended set of
metrics to provide a better picture of the quality of the datasets:



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




QC Figures
~~~~~~~~~~

Useful figures at the top of the Quality Control window contain
