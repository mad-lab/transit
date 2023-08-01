CRISRi-DR - Chemical-Genetic Interaction Analysis of CRISPRi libraries based on a Dose-Response model
============================================================



Overview 
--------------------
This software is designed to analyze CRISPRi libraries from CGI experiments and identify significant CGIs ie genes that affect sensitivity to the drug when depleted. 
[REF: TBA]


Workflow
--------
FLOWCHART?
Big picture: start with fastq files, extract barcode counts, manually create a metadata file, run the model. the output file lists genes and their statistacal parameters and significance. option to visual specific genes at the sgRNA level.

**genes with sig interactions are those with qval<0.05 and |z|>2 on the slope coeffcient**




_Preprocessing: Fastq to Count Files_

This is a longer process, taking a few minutes. However, the number of reads processed is printed to the console to indicate progress
> Usage : python3 ../src/transit.py CGI extract_counts <fastq_file> <ids_file> > <counts_file>

* ids_file : list of sgRNAs used in the experiment, where each row is one sgRNA id
  * for H37Rv experiments, the ids file is available in transit/src/pytransit/data/IDs.H37Rv.CRISPRi.lib.txt


_Step 1: Combine Individual Counts File to a Combined Counts File_

Pass in n individual counts files seperated by spaces and the name of the combined counts file of your choice
All the counts files should have sgRNA ids as their first column, and can have any number of columns. 
The comma seperated headers must be in order of the columns in the counts files entered. 

> Usage : python3 ../src/transit.py CGI combined_counts <comma seperated headers> <counts_file_1> <counts_file_2>  ... <counts_file_n> > <combined_counts_file>

_Step 2: Extract Fractional Abundances_

This step is to turn the barcodes extracted into relative normalized abundances. These are normalized within samples, relative to the abundances in the ATC-induced 0-concentration file (the <no_depletion_abundances_file>), essentially fractions. This is a relatively quick process, taking less than a minute.

> Usage : python3 ../src/transit.py CGI extract_abund <combined_counts_file> <metadata_file> <reference_condition> <extrapolated_LFCs_file> <no_depletion_abundances_file> <drug> <days>  >  <frac_abund_file>

* metadata_file (USER created):
  * The columns expected in this file: column_name,drug,conc_xMIC,days_predepletion
  * See ShuquiCGI_metadata.txt for an example of this type file
  * You do not need equal number of replicates for all concentrations
  * see [Shuqi REF] for explanation of days_predepletion
* reference_condition: the condition to calculate relative abundances from as specificed in the 'drug' columns of the metadata file; typically an ATC-induced, no drug concentration
* extrapolated_LFCs_file: A file that contains metadata for each sgRNA in the combined counts file
  * The last column must be extrapolated LFCs calculated through a passaging experiment. This column will be used as a measurement of sgRNA strength in the CRISPRi-DR model 
* no_depletion_abundances_file: A two column file of sgRNAs and their abundances in -ATC-induced (no ATC) with 0 drug concentration 
* drug : Name of the drug in the "drug" column of the metadata file passed in to be fit in the model
* days: Sampled from predepletion day as listed in the "days_predepletion" column of the metadata file to be used in the analysis


_Step 3: Run the CRISPRi-DR model_

This is a relatively quick process, taking less than a minute. This step fits the CRISPRi-DR model (statistical analysis of concentration dependence for each gene) to each gene in the file and prints each output to the <CRISPRi-DR_results>
In this spreadhseet, siginificant interacting genes are those with adjusted P-val (Q-val) < 0.05 and |Z slope|>2

> Usage : python3 ../src/transit.py CGI run_model <abund_file>  >  <CRISPRi-DR_results>



Example
-------

Note that the first step requires some data files.
* ShiquiCGI_metadata.txt - describes the samples
* Bosch21_TableS2.txt - contains betaE estimates for each sgRNA
* Bosch21_TableS2_extended.txt - contains extrapolated LFCs for each sgRNA
* no_depletion_abundances.txt - pre-calculated abundance for -ATC (no induction of target depletion)

Preliminary step: download raw counts from github
  > git clone https://github.com/rock-lab/CGI_nature_micro_2022

  I suggest linking this in the local (CGI) directory as follows 
    (but I can go anywhere, and you provide the path to the data/counts/ dir on the command line for extract_abund)

  > ln -s CGI_nature_micro_2022/data data

usage: 
  python3 ../src/transit.py CGI extract_abund <metadata_file> <data_dir> <extrapolated_LFCs_file> <no_drug_file> <no_depletion_abundances_file> <drug> <days>  >  <output_file>
  python3 ../src/transit.py CGI run_model <abund_file>  >  <logsigmodfit_file>
  python3 ../src/transit.py CGI post_process <logsigmoidfit_file>  >  <results_file>

example of pipeline:

> python3 ../src/transit.py CGI extract_abund ShiquiCGI_metadata.txt data/counts/ Bosch21_TableS2_extended.txt counts_1972_DMSO_D5.txt no_depletion_abundances.txt RIF 5 > frac_abund.RIF_D5.txt

  gathers relevant samples (at all available concs, and DMSO representing 0xMIC)
  calculates fractional abundances (normalized within samples, relative to no-depletion abundances, essentially fractions)


> python3 ../src/transit.py CGI run_model frac_abund.RIF_D5.txt > logsigmoidfit.RIF_D5.txt

  runs linear regressions for a log-sigmoid model (in R)

> python3 ../src/transit.py CGI post_process logsigmoidfit.RIF_D5.txt > CGI_results.RIF_D5.txt

  outputs statistical analysis of concentration dependence for each gene
  can open as spreadsheet in Excel
  genes that interact significantly with drug are those with adjusted P-val (Q-val) < 0.05



python3 ../transit/src/transit.py CGI create_combined_counts DMSO_1,DMSO_2,DMSO_3,VAN_0_0625_1,VAN_0_0625_2,VAN_0_0625_3,VAN_0_125_1,VAN_0_125_2,VAN_0_125_3,VAN_0_25_1,VAN_0_25_2,VAN_0_25_3 CGI_nature_micro_2022/data/counts/counts_1952_DMSO_D10.txt CGI_nature_micro_2022/data/counts/counts_1953_VAN_0_0625X_D10.txt CGI_nature_micro_2022/data/counts/counts_1954_VAN_0_125X_D10.txt CGI_nature_micro_2022/data/counts/counts_1955_VAN_0_25X_D10.txt > combined_VAN_D10.txt

python3 ../transit/src/transit.py CGI extract_abund combined_VAN_D10.txt VAN_D10_metadata.txt DMSO ../transit/CGI/Bosch21_TableS2_extended.txt ../transit/CGI/no_depletion_abundances.txt VAN 10  >  VAN_D10_frac_abund.txt

python3 ../transit/src/transit.py CGI run_model VAN_D10_frac_abund.txt > logsig.txt