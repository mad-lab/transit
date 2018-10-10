


# TRANSIT 2.3.0


[![Build Status](https://travis-ci.org/mad-lab/transit.svg?branch=master)](https://travis-ci.org/mad-lab/transit)   [![Documentation Status](https://readthedocs.org/projects/transit/badge/?version=latest)](http://transit.readthedocs.io/en/latest/?badge=latest) 


**Version 2.3.0 changes (Oct, 2018)**
- Added calculation of Pathway Enrichment as post-processing for resampling, to determine if conditionally essential genes over-represent a particular functional category or pathway (such as for GO terms)
- Added ANOVA analysis for identifying genes with significant variability of counts across multiple conditions
- Updated Documentation - especially for "Quality Control/TnSeq Statistics"; also added more command-line examples under "Analysis Methods"
- Fixed bugs (including TrackView in the GUI)
- Upgraded dependencies, including wxPython 4.0 (required)


**Version 2.2.0 changes (June, 2018)**
- Added analysis method for Genetic Interactions.
- Added Mann-Whitney U-test for comparative analysis.
- Made TRANSIT compatible with wxPython 4.0 (Phoenix).
- Datasets now automatically selected when they are added to TRANSIT.
- TRANSIT window now starts maximized.
- Updated documentation.
- Fixed bug with plots of finished results files.
- Fixed bug in packaging of TPP, causing problem with console mode in new setuptools.
- Other misc. bugs fixes


**Version 2.1.2 changes (May, 2018)**
- Improved resampling on comparisons with unbalanced number of replicates.
- Improved speed of TPP.
- Added Barseq functionality to TPP.
- Fixed bug with creating resampling histograms on GUI mode.
- Miscellaneous code improvements.


**Version 2.1.1 changes (July, 2017)**
- Miscellaneous bug fixes


**Version 2.1.0 changes (June, 2017)**
- Added tooltips next to most parameters to explain their functionality.
- Added Quality Control window, with choice for normalization method.
- Added more normalization options to the HMM method.
- Added LOESS correction functionality back to TRANSIT
- Added ability to scale Track View based on mean-count of the window.
- Added ability to scale individual tracks in Track View.
- Added ability to add tracks of features to Track View.
- Better status messages for TrackView
- TPP can now accept empty primer prefix (in case reads have been trimmed).
- TPP can now process reads obtained using Mme1 enzyme and protocol.
- TPP can now pass flags to BWA.
- Lots of bug fixes.


**Version 2.0.2 changes (August, 2016)**
- Added support for for custom primers in TPP.
- Added support for annotations in GFF3 format.
- Ability to specify pseudocounts in resampling.
- Misc. Bug fixes
- **New [mailing list](https://groups.google.com/forum/#!forum/tnseq-transit/join)**


**New in Version 2.0+**
 - Support for Tn5 datasets.
 - New analysis methods.
 - New way to export normalized datasets.



Welcome! This is the distribution for the TRANSIT and TPP tools developed by the Ioerger Lab.

TRANSIT is a tool for the analysis of Tn-Seq data. It provides an easy to use graphical interface and access to three different analysis methods that allow the user to determine essentiality in a single condition as well as between conditions.


## Features
TRANSIT offers a variety of features including:
    
-   More than **8 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

-   Ability to analyze datasets from libraries constructed using  **himar1 or tn5 transposons**.

-   **TrackView** to help visualize read-counts accross the genome.

-   Can **export datasets** into a variety of formats, including **IGV**.

-   Includes a **variety of normalization methods**.

-   **Quality Control** diagnostics, to idenfity poor quality datasets.

-   Ability to install as a **python package**, to import and use in your own personal scripts.





## Mailing List

You can join our mailing list to get announcements of new versions, discuss any bugs, or request features! Just head over to the following site and enter your email address:

https://groups.google.com/forum/#!forum/tnseq-transit/join




## Instructions

For full instructions on how to install and run TRANSIT (and the optional pre-processor, TPP), please see the documentation included in this distribution ("doc" folder) or visit the following web page:


http://saclab.tamu.edu/essentiality/transit/transit.html


## Datasets

The TRANSIT distribution comes with some example .wig files in the data/ directory, as well as an example annotation file (.prot\_table format) in the genomes/ directory. Additional genomes may be found on the following website:

http://saclab.tamu.edu/essentiality/transit/genomes/


## Copyright Information

Source code for TRANSIT and TPP are available open source under the terms of the GNU General Public License (Version 3.0) as published by the Free Software Foundation. For more information on this license, please see the included LICENSE.md file or visit their website at:

http://www.gnu.org/licenses/gpl.html
