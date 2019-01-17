# TRANSIT 2.4.1

[![Build Status](https://travis-ci.org/mad-lab/transit.svg?branch=master)](https://travis-ci.org/mad-lab/transit)   [![Documentation Status](https://readthedocs.org/projects/transit/badge/?version=latest)](http://transit.readthedocs.io/en/latest/?badge=latest) 

**TRANSIT 2.4.1 by Robert Jenquin and William Matern**
- Added the ability to accept multiple replicons in the form of either multiline reference genomes or multiple reference genome files.
- Added `-bwa-alg` argument, allowing the user to specify `mem` or `aln` to use `bwa mem` or `bwa aln` algorithms
- Now requires `-replicon-id` argument to specify names for the replicons if multiple reference genomes given (respective order to order appearing in reference genome(s)
- Code cleanup: closing dangling file handles
- Bug fix: if adapter is at exact end of R1, it is now properly handled
- Bug fix: trimmed\_reads now counted properly 
- Added support for specifying `-window-size` argument
- **Sample usage:**
```
python2 src/tpp.py -himar1 -bwa /usr/bin/bwa -bwa-alg aln -ref MAC109_genome.fa -replicon-id CP029332 CP029333 CP029334 -reads1 ../HJKK5BCX2_ATGCTG_1.fastq -reads2 ../HJKK5BCX2_ATGCTG_2.fastq -primer AACCTGTTA -mismatches 2 -window-size 6 -output tpp_output/avium
```
- Explanation of arguments
  - `-himar1` specifies that the Himar1 transposon was used in the transposon mutagenesis procedure. Tn5 is also supported (`-tn5`)
  - `-bwa` specifies the path to the `bwa` executable
  - `-bwa-alg` specifies either `mem` or `aln` algorithms for `bwa` to use. `aln` is widely considered obselete to `mem` for reads of length > 70bp. `aln` is default.
  - `-ref` specifies the reference genome(s) in FASTA format to which reads will be mapped. If more than one, they can be specified in either multiple FASTAs, or as a multilined FASTA (or a combination of both).
  - `-replicon-id [contig1 contig2 ...]` specifies the names of the contigs in the genome(s). These are used as filename suffixes for output files (ie \*\_contig1.wig, \*\_contig2.wig, etc). The order of the contigs is assumed to be the same as they appear in the reference genome(s) (as given with `-ref`). Specifying this option is only required if there is more than one contig. Note: While you can technically use any contig name at this step, if you wish to use `wig_gb_to_csv.py` to organize the data you should use the contig names as they appear in the Genbank file (as specified by `wig_gb_to_csv.py -g`).
  - `-reads1` specifies the file containing the raw reads (untrimmed) for read1 in FASTQ or FASTA format
  - `-reads2` specifies the file containing the raw reads (untrimmed) for read2 in FASTQ or FASTA format
  - `-primer` specifies a nucleotide sequence at the end of the transposon, is used to separate transposon DNA from genomic DNA in read 1.
  - `-window-size` specifies how many positions to look for `-primer` within read 1. It should be set to at least the difference between the maximum and minumum expected positions of the first base of genomic DNA in read 1 (and larger if you want to allow for insertions/deletions). For the Long et al 2015 protocol (using a pool of 4 shifting prefixes) the window-size should be at least 6. Default value is 6.
  - `-mismatches` specifies the number of mismatches to allow when searching for the transposon in read 1 (ie number of mismatches to `-primer`).
  - `-output` specifies the filename prefix to be applied to output files. Can include directories, allowing custom paths to be specified.

# TRANSIT 2.3.4


[![Build Status](https://travis-ci.org/mad-lab/transit.svg?branch=master)](https://travis-ci.org/mad-lab/transit)   [![Documentation Status](https://readthedocs.org/projects/transit/badge/?version=latest)](http://transit.readthedocs.io/en/latest/?badge=latest) 


Welcome! This is the distribution for the TRANSIT and TPP tools developed by the Ioerger Lab at Texas A&M University.

TRANSIT is a tool for processing and statistical analysis of Tn-Seq data. 
It provides an easy to use graphical interface and access to three different analysis methods that allow the user to determine essentiality in a single condition as well as between conditions.

TRANSIT Home page: http://saclab.tamu.edu/essentiality/transit/index.html

TRANSIT Documentation: https://transit.readthedocs.io/en/latest/transit_overview.html

[Changelog](https://github.com/mad-lab/transit/blob/master/CHANGELOG.md)


## Features
TRANSIT offers a variety of features including:
    
-   More than **8 analysis methods**, including methods for determining **conditional essentiality** as well as **genetic interactions**.

-   Ability to analyze datasets from libraries constructed using  **himar1 or tn5 transposons**.

-   **TrackView** to help visualize read-counts accross the genome.

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
