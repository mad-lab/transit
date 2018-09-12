# Change log
All notable changes to this project will be documented in this file.


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


