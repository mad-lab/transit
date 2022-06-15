.. transit documentation master file, created by
   sphinx-quickstart on Wed May  4 13:37:43 2016.

Welcome to TRANSIT's documentation!
===================================

.. image:: https://img.shields.io/github/tag/mad-lab/transit.svg
    :target: https://github.com/mad-lab/transit
    :alt: GitHub last tag

Transit is python-based software for analyzing TnSeq data
(sequencing data from transposon mutant libraries)
to determine essentiality of bacterial genes under different conditions.

This page contains the documentation for TRANSIT. Below are a few
quick links to some of the most important sections of the
documentation, followed by a brief overview of TRANSIT's features.

Quick Links
~~~~~~~~~~~
.. _quick-link:


* :ref:`install-link`
* :ref:`manual-link`
* :ref:`tutorial-link`
* :ref:`tpp-link`
* :ref:`code-link`
* `PDF manual with overview of analysis methods in Transit <https://orca1.tamu.edu/essentiality/transit/transit-manual.pdf>`_

Features
~~~~~~~~
TRANSIT offers a variety of features including:
 
*   More than **8 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

*   Ability to analyze **himar1 or tn5 transposons** datasets.

*   **TrackView** to help visualize read-counts accross the genome.

*   Can **export datasets** into a variety of formats, including **IGV**.

*   Includes a **variety of normalization methods**.

*   **Quality Control** diagnostics, to idenfity poor quality datasets.

*   Ability to install as a **python package**, to import and use in your own personal scripts.   



.. 
 Mailing List
 ~~~~~~~~~~~~
 
 You can join our mailing list to get announcements of new versions, discuss any bugs, or request features! Just head over to the following site and enter your email address:
 
 
  + `https://groups.google.com/forum/#!forum/tnseq-transit/join <https://groups.google.com/forum/#!forum/tnseq-transit/join>`_
 
..


.. _manual-link:

.. toctree::
   :maxdepth: 3
   :caption: TRANSIT MANUAL

   transit_overview
   transit_install
   transit_running
   transit_features

.. toctree::
   :maxdepth: 3
   :caption: ANALYSIS METHODS

   method_gumbel
   method_griffin
   method_tn5gaps
   method_HMM
   method_resampling
   method_Utest
   method_GI
   method_anova
   method_zinb
   method_normalization
   method_pathway_enrichment
   method_tnseq_stats
   method_corrplot
   method_heatmap
   method_ttnfitness

.. _tutorial-link:

.. toctree::
   :maxdepth: 3
   :caption: TRANSIT Tutorials

   transit_essentiality_tutorial
   transit_genome_tutorial
   transit_comparative_tutorial
   transit_interactions_tutorial
   transit_normalization_tutorial
   transit_export_tutorial
   transit_console_cheatsheet

.. _tpp-link:

.. toctree::
   :maxdepth: 3
   :caption: TPP Manual

   tpp.rst

.. _code-link:

.. toctree::
   :maxdepth: 3
   :caption: Code Documentation

   transit

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Support
~~~~~~~

For any questions or comments, please contact Dr. Thomas Ioerger, ioerger@cs.tamu.edu.


