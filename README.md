


# TRANSIT 2.1.0


[![Build Status](https://travis-ci.org/mad-lab/transit.svg?branch=master)](https://travis-ci.org/mad-lab/transit)   [![Documentation Status](https://readthedocs.org/projects/transit/badge/?version=latest)](http://transit.readthedocs.io/en/latest/?badge=latest) 



**Version 2.1.0 changes (June, 20017)**
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
