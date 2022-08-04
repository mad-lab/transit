# Change log
All notable changes to this project will be documented in this file.


## Version 3.2.6 2021-08-03
#### TRANSIT:

Major changes:
 - added a parameter 'alpha' to ANOVA to make the F-test less sensitive to genes with low counts, cutting down on 'irrelevant' genes with significant variability
 - updated the online documentation to describe this

Minor changes:
 - fixed a (recently-introduced) bug that was causing the GUI to crash when running resampling
 - updated 'export combined_wig' to include ref genome and column headers


## Version 3.2.5 2021-06-15
#### TRANSIT:

Minor changes:
 - update dependencies for pillow and sklearn
 - refactor documentation (replace transit_methods.rst with separate .rst files)
 - added rpy2 warning (if not installed) for corrplot and heatmap


## Version 3.2.4 2021-06-05
#### TRANSIT:

Major changes:
 - added 'ttnfitness' analysis method (to categorize growth-defect genes in single (reference) conditions, and compute TTN-fitness ratio to quantify the magnitude of growth defect based on comparison of observed insertion counts to expected counts at each TA site (based on surrounding nucleotides)
 - added winzorization (-winz flag) to resampling, ANOVA, and ZINB (to help mitigate effects due to sites with outlier counts)
 - fixed bug in ANOVA that assumed files in combined_wig and metadata were listed in same order (now they don't have to be)

Minor changes:
 - switched back to original implementation of mmfind() 
 - added pathway assocation files for M. smegmatis to data dir 
 - added --Pval_col and --Qval_col to pathway_enrichment.py 
 - added --prot_table flag to zinb.py
 - updated header info in output files for HMM, resampling, and ZINB
 - updated explanation of -signif in documentation for Genetic Interactions
 - cleaned up documentation

	
## Version 3.2.3 2021-10-16
#### TRANSIT:
  - added Binomial essentials (EB) to Gumbel analysis (supplementing genes categorized as E), to help with low-saturation datasets
  - modified ANOVA and ZINB so that --include-conditions and --exclude-conditions refer to original Conditions column in samples metadata file (instead of whatever is specified by --conditions) 
	
#### TPP:
  - improved metrics reported in *.tn_stats by TPP, to help diagnose why reads don't map
	
## Version 3.2.2 2021-09-08
#### TRANSIT:
 - fixed bug in converting gff_to_prot_table
 - fixed bug in tn5gaps (fixes some false negative calls)
 - fixed some bugs in pathway_enrichment (GSEA calculations)
 - fixed links to Salmonella Tn5 data in docs
 - fixed problem with margins in heatmap.py that was causing R to fail 
 - added --ref to anova.py and zinb.py (for computing LFCs relative to designated reference condition)
 - added --low_mean_filter for heatmap.py (for excluding genes with low counts, even if they are significant by ANOVA or ZINB)
 - add dependency on pypubsub<4.0
	
## Version 3.2.1 2020-12-22
#### TRANSIT:
 - maintenance release
  - fixed a bug in the GUI caused by changes in wxPython 4.1.0
  - added GO terms for M. smegmatis in the data directory for doing pathway analysis

## Version 3.2.0 2020-10-26
#### TRANSIT:
 - improvements to pathway_enrichment analysis
  - added '--ranking' flag for GSEA to sort genes based on LFC or SLPV
  - implemented Ontologizer method (-M ONT), which works better for GO terms
  - updated auxilliary files in transit data directory for different systems of functional categories (COG, Sanger, GO)
 - added '-signif' flag to GI (Genetic Interaction analysis) (options: HDI, prob, BFDR, FWER)
  - updated description of methods for determining significant interactions in documentation
 - various improvements to other methods, including corrplot and heatmap

## Version 3.1.0 2020-03-08
#### TRANSIT:
 - added 'corrplot' and 'heatmap' commands
 - pathway_enrichment: 
  - completely re-done so it is faster and simpler
  - now implements Fisher's exact test and GSEA
  - can be used with COG categories and GO terms
  - switch to 2-column format for associations files
 - resampling: 
  - changed semantics of pseudocounts from "fake sites" (-pc, dropped) to calculation of log-fold-changes (-PC, new)
 - anova: 
  - put columns for condition means in correct order
  - added columns for log-fold-changes for each condition to output
 - zinb: 
  - improved handling of --include-conditions and --ignore-conditions
  - now prints out a summary of how many samples are in each condition (including cross-product with covars and interactions)
 - make pseudocounts flag (-PC) work uniformly for resampling, anova, and zinb

	
## Version 3.0.2 2019-12-21
#### TRANSIT:
 - Mostly cosmetic fixes
 - Updated some command-line and GUI messages
 - Updated documentation (especially for GI and resampling)
 - Removed "warning: high stderr" from gene status in ZINB
 - Added LFCs in ZINB output
 - updated 'convert gff_to_prot_table' so it works with gff3 files downloaded from NCBI

## Version 3.0.1 2019-08-01
#### TRANSIT:
 - Add check for python3 (TRANSIT 3+ requires python3.6+)
 - Minor fixes in GI and Pathway Enrichment

## Version 3.0.0 2019-07-18
#### TRANSIT:
 - TRANSIT now supports python 3. (To use with python 2, use releases < 3.0.0)
 - Improved speed of GSEA in Pathway Enrichment analysis.

## Version 2.5.2 2019-05-16
#### TRANSIT:
 - Made some improvements in command-line version of 'tn5gaps'
 - Added flags for trimming insertions in N- and C-termini of genes for tn5gaps (-iN and -iC) 

## Version 2.5.1 2019-04-25
#### TRANSIT:
  - Add support for [handling interactions in ZINB](https://transit.readthedocs.io/en/latest/transit_methods.html#covariates-and-interactions)
  - Fix selection bug for gff3 in GUI

## Version 2.5.0 2019-03-28
#### TRANSIT:
  - Added analysis method for Zero-Inflated Negative Binomial ([ZINB](https://transit.readthedocs.io/en/latest/transit_methods.html#zinb))
  - Fix LOESS flag bug in resampling 2.4.2
  - Resampling supports combined_wig files
  - Change ordering of metadata and annotation file in ANOVA cmd

## Version 2.4.2 2019-03-15
#### TPP:
 - updated docs for TPP; expanded discussion of protocols, including Mme1 
 - for Mme1, change min read length from 20bp to 15bp (for genomic part of read1)
 - replaced '-himar1' and 'tn5' flags with '-protocol [sassetti|tn5|mme1]'
 - added 'auto' for -replicon-ids
 - added 'pre-trimmed' as option for transposon in TPP GUI (prefix="")
#### TRANSIT:
 - [resampling can now be done between TnSeq libraries from different strains](https://transit.readthedocs.io/en/latest/transit_methods.html#re-sampling)
 - add documentation for 'griffin' and Mann-Whitney 'utest' analysis methods

	
## Version 2.4.1 2019-03-04
#### TPP:
 - allow the primer sequence to be the empty string (i.e. -primer "" on command-line; for pre-trimmed reads)
 - do not throw an error if header ids in read1 and read2 fastq files happen to match identically
 - minor bug fixes:
 - fixed problem of order of data in tn_stats table when there are multiple contigs but only single-ended reads
 - fixed name of flag from "replicon-id" to "replicon-ids"
 - prevent div-by-zero error in cases where no reads map
	

## Version 2.4.0 2019-02-28
#### TPP:
 - **can now handle genomes with multiple contigs** (thanks to modifications by Robert Jenquin and William Matern); it creates multiple .wig files as output
 - BWA: switched from using 'aln' to 'mem' by default
 - added flags to set the nucleotide window for searching for start of primer sequence (-primer-window-start)
 - fixed bug in counting misprimed reads, and reads mapped to both R1 and R2
 - added some fields to TPP GUI, and made it more consistent about saving/reading parameters in the tpp.cfg config file
#### Transit:	
 - fixed bug in handling '-minreads' flag in Gumbel analysis
 - updated support for converting .gff files to .prot_table format (in GUI and on command line)
 - added a status field to ANOVA output
 - TrackView scales all plots simultaneously by default
- updated documentation

## Pull Request 18 by Robert Jenquin and William Matern (Jan, 2019)

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

## Version 2.3.4 2019-01-14
- TRANSIT:	
  - Minor bug fixes related to flags in Resampling and HMM

## Version 2.3.3 2018-12-06
- TRANSIT:	
  - Minor bug fixes related to flags in HMM

## Version 2.3.2 2018-11-09
- TRANSIT:	
  - Minor bug fixes related to changing parameters in TPP GUI

## Version 2.3.1 2018-10-19
- TRANSIT:	
  - Removed dependence on PyPubSub (can run Transit in command-line mode without it, but needed for GUI)

	
## Version 2.3.0 2018-10-10
- TRANSIT:	
  - Added calculation of Pathway Enrichment as post-processing for resampling, to determine if conditionally essential genes over-represent a particular functional category or pathway (such as for GO terms)
- Added ANOVA analysis for identifying genes with significant variability of counts across multiple conditions
  - Updated Documentation - especially for "Quality Control/TnSeq Statistics"; also added more command-line examples under "Analysis Methods"
  - Fixed bugs (including TrackView in the GUI)
  - Upgraded dependencies, including wxPython 4.0 (required)
	


## Version 2.2.0 - 2018-06-04
- TRANSIT:
    - Added analysis method for Genetic Interactions.
    - Added Mann-Whitney U-test for comparative analysis.
    - Made TRANSIT compatible with wxPython 4.0 (Phoenix).
    - Datasets now automatically selected when they are added to TRANSIT.
    - Fixed bug in packaging of TPP, causing problem with console mode in new setuptools.
    - Miscellaneous bugs fixes 

	

## Version 2.1.0 - 2017-06-23
- TRANSIT:
    - Added tooltips next to most parameters to explain their functionality.
    - Added Quality Control window, with choice for normalization method.
    - Added more normalization options to the HMM method.
    - Added LOESS correction functionality back to TRANSIT
    - Added ability to scale Track View based on mean-count of the window.
    - Added ability to scale individual tracks in Track View.
    - Added ability to add tracks of features to Track View.
    - New documentation on normalization.

- TPP:
    - TPP can now accept empty primer prefix (in case reads have been trimmed).
    - TPP can now process reads obtained using Mme1 enzyme and protocol.
    - TPP can now pass flags to BWA.



## Version 2.0.2 - 2016-08-19
    
- TRANSIT:
    - Now accepts GFF3 formatted annotations.
    - Added ability to specify pseudocounts for resampling.
    - Added extra columns to resampling output.
    - Fixed bug with some log2FC calculations.
    - Export to combined wig format now asks for normalization BEFORE file name.
    - Fixed bug preventing Quality Control window from opening.
    - Miscellanous bug fixes.
    - Updates to Documentation

- TPP:
    - Now accepts custom primer sequences.
    - Reporting additional diagnostic statistics for reads mapping to phiMycoMarT7, and Illumina adapters.
    - Miscellaneous bug fixes.
   


## Version 2.0.1 - 2016-07-05

-TRANSIT:
    - Fixed crash in TPP.
    - Misc changes for outputs.


## Version 2.0.0 - 2016-06-16

- TRANSIT: 
    - Added new method for datasets created with Tn5 transposons.
    - Added label indicating intended transposons for the methods.
    - Added textbox with short description of the chosen method.
    - Changed methods choices to be in menu (on top).
    - Changed the file display window.
    - Added Help menu with link to online documentation.
    - Added new logo.
    - Added option to export (normalized) datasets to IGV or combined wig format.
    - Can now select multiple .wig files at the same time (Ctrl + select).
    - Lots of changes under the hood. 


## Version 1.4.5 - 2016-01-10
- TRANSIT:
    - Added Binomial analysis method as an option to TRANSIT.
    - Added DE-HMM analysis method as an option to TRANSIT.
 

## Version 1.4.3.1 - 2016-01-02
- TRANSIT:
    - Fixed bug causing TRANSIT not to open on some Windows systems.


## Version 1.4.3 - 2015-12-04
- TRANSIT:
    - Precision of resampling p-values in output file now increases with sample size
    - Added preliminary Quality Control functionality. Select some datasets and click View -> Quality Control
    - In resampling, changed logFC to divide by number of replicates
    - Changed plotting of results files to be more versitile
    - Fixed bug causing HMM_sites output not to be added to list of files
    - Fixed bug causing LOESS correction not to work in HMM


## Version 1.4.2 - 2015-07-29
- TRANSIT:
    - Added Total Trimmed Reads normaliztion (TTR) as the default option. This is the recommended normalization method at this point.
    - Added BetaGeomtric Correction (betageom) as a normalization option. This is recommended for datasets that are very skewed.
    - Fixed bug that caused transit to create histograms when not desired.
    - Added a pseudo-count when calculating log-FC to genes without reads.
    - Increased size of result windows so that all columns are immediately visible.




## Version 1.4.1 - 2015-06-5
- TRANSIT:
    - TRANSIT now accepts read-counts in floating-point precision, not just integers.
    - Made transit work with most recent versions of matplotlib.


## Version 1.4.0 - 2015-05-27
- TRANSIT:
    - Added option to correct for genomic position bias (using LOESS)
    - Added more options for normalization, including zero-inflated negative binomial and quantile normalization.


- TPP:
    - Eliminated soft-clipped reads.
    - Modified template_counts() to be much more memory efficient (does not need gigabytes of RAM any more to process large datasets)
    - Added ability to process Tn5 datasets





## Version 1.3.0 - 2015-03-31
- TRANSIT:
    - Fixed threading issue for volcano plot.
    - Improved format and quality of the output messages.
    - Fixed direction of log-fold change in volcano plots.
    - Added log-fold change column to resampling output file.
    - Made adaptive resampling work better with custom sample sizes.


- TPP:
    - Fixed genomic portion for single ends.
    - Added usage help as part of command line arguments.



## Version 1.2.33 - 2015-03-06

- TRANSIT:
    - Fixed issue with histograms create using adaptive resampling.



## Version 1.2.32 - 2015-03-05

- TRANSIT:
    - Put .pyc files in in new src/ directory.
    - Fixed error that sometimes occurred when plotting volcano plots.
    - Made TRANSIT default to the current working directory when opening file dialogs.

- TPP:
    - TPP can now process files with single-end reads.
    - TPP can now process *.fasta and compressed files with "*.fastq.gz" extension



## Version 1.2.7 - 2015-02-25

- TRANSIT:
    - Fixed error that occured when displaying graphs after running an analysis.
    - Updated datasets included in the data/ directory.

- TPP:
    - Removed the requirement for wxPython when running TPP on command-line mode.



## Version 1.1 - 2015-02-20

- TRANSIT:
    - Fixed error in HMM results file table, which was not correctly showing breakdown of genes.
    - Made TRANSIT work from the command-line, without displaying GUI. See documentation for arguments/flags.
    - Added ability to convert annotation files between several formats (.prot_table, ptt.table, gff3).

- TPP:
    - User can supply reads in either FastA or FastQ format.
    - Added an option to specify number of mismatches (default=1) when looking
       for sequence patterns such as the transposon prefix in read 1.
    - Added command-line arguments so TPP can be run in batch mode without the GUI.
    - Number of mapped reads for R1 and R2 independently is also now reported.
    - Modified how barcodes are extracted from read 2.  It now looks for specific
       sequence patterns, even if they are shifted.  This should greatly increase
       the number of mapped reads (esp. the genomic part of R2) for certain datasets.
    - Properly handle short fragments, ie. for reads where the insert size is shorter
       than the read length.  In such cases, the adapter from other end appears
       at the end of read 1, and this suffix is now stripped off so these reads 
       will map too.



## Version 1.0  - 2015-02-10

- First limited-release version of TRANSIT
- Released to close collaborators first and presented in teleconference to get feedback.


