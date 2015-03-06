# Change log
All notable changes to this project will be documented in this file.


## Version 1.2.32 - 2015-03-05

-TRANSIT:
    - Put .pyc files in in new src/ directory.
    - Fixed error that sometimes occurred when plotting volcano plots.
    - Made TRANSIT default to the current working directory when opening file dialogs.

- TPP:
    - TPP can now process files with single-end reads.
    - TPP can now process compressed files with "*.fasta.gz" extension



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


