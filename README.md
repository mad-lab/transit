# TRANSIT

[![Version](https://img.shields.io/github/tag/mad-lab/transit.svg)](https://github.com/mad-lab/transit)   [![Build Status](https://travis-ci.org/mad-lab/transit.svg?branch=master)](https://travis-ci.org/mad-lab/transit)   [![Documentation Status](https://readthedocs.org/projects/transit/badge/?version=latest)](http://transit.readthedocs.io/en/latest/?badge=latest)   [![Downloads](https://pepy.tech/badge/tnseq-transit)](https://pepy.tech/project/tnseq-transit)

=======

**NOTE: TRANSIT v3.0+ now requires python3.6+. If you want to use TRANSIT with python2, use version < 3.0.**

Welcome! This is the distribution for the TRANSIT and TPP tools developed by the [Ioerger Lab](http://orca2.tamu.edu/tom/iLab.html) at Texas A&M University.

TRANSIT is a tool for processing and statistical analysis of Tn-Seq data.
It provides an easy to use graphical interface and access to three different analysis methods that allow the user to determine essentiality in a single condition as well as between conditions.

TRANSIT Home page: http://saclab.tamu.edu/essentiality/transit/index.html

TRANSIT Documentation: https://transit.readthedocs.io/en/latest/transit_overview.html

[Changelog](https://github.com/mad-lab/transit/blob/master/CHANGELOG.md)


## Features
TRANSIT offers a variety of features including:

-   More than **10 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

-   Ability to analyze datasets from libraries constructed using  **himar1 or tn5 transposons**.

-   **TrackView** to help visualize read-counts across the genome.

-   Can **export datasets** into a variety of formats, including **IGV**.

-   Includes a **variety of normalization methods**.

-   **Quality Control** diagnostics, to idenfity poor quality datasets.

-   Ability to install as a **python package**, to import and use in your own personal scripts.



## Support

For any questions or comments, please contact Dr. Thomas Ioerger, ioerger@cs.tamu.edu.




## Instructions

For full instructions on how to install and run TRANSIT (and the optional pre-processor, TPP), please see the documentation included in this distribution ("src/pytransit/doc" folder) or visit the following web page:


https://transit.readthedocs.io/en/latest/


## Datasets

The TRANSIT distribution comes with some example .wig files in the data/ directory, as well as an example annotation file (.prot\_table format) in the genomes/ directory. Additional genomes may be found on the following website:

http://saclab.tamu.edu/essentiality/transit/genomes/


## Copyright Information

Source code for TRANSIT and TPP are available open source under the terms of the GNU General Public License (Version 3.0) as published by the Free Software Foundation. For more information on this license, please see the included LICENSE.md file or visit their website at:

http://www.gnu.org/licenses/gpl.html
