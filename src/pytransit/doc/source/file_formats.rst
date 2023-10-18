.. _input_files:

File Formats for Transit
=======================


Wig Files
---------

TPP processes raw sequencing files (.fastq) to extract transposon insertion counts 
at TA sites.
These are stored as ".wig" files.
Wig files primarily have 2 columns of data: a) coordinates of TA sites, 
and b) the number of insertions observed at those sites.
By convention, there are also 2 header lines.
The first header line can be arbitrary info, such as where the data came from.
The second has to contain "variableStep chrom=" followed by the name of the
reference genome used for mapping the reads (e.g. H37Rv.fna).
Here is an example (transit/src/pytransit/data/glycerol_H37Rv_rep2.wig):

::
  
  # from Griffin et al, (2011).  PLOS Pathogens, e1002251.
  variableStep chrom=H37Rv
  60 0
  72 0
  102 0
  ...
  1522 124
  1552 0
  1635 201
  1779 0
  1782 0
  1788 0
  1847 0
  1858 342
  1921 21
  ...
  4411440 99
  4411479 0
  4411508 0
  4411526 0

  (full file has ~75,000 lines)

|

.. _annotation_files:

Prot_Tables (Genome Annotations)
------------------

The annotation of a genome contains information about genes, such as
coordinates, strand, locus tag, gene name, and functional description.
Transit uses a custom format for annotations called "prot_table"s,
e.g. H37Rv.prot_table.  Prot_tables are **tab-separated text files**
containing the gene information in 9 specific columns:

**Prot_table file format:**

1. gene function description
2. start coordinate
3. end coordinate
4. strand
5. length of protein product (in amino acids)
6. don't care
7. don't care
8. gene name (like "dnaA")
9. ORF id (like Rv0001)

Here is an example (transit/src/pytransit/genomes/H37Rv.prot_table):

::

  chromosomal replication initiation protein 	1	1524	+	507	15607143	885041	dnaA	Rv0001
  DNA polymerase III subunit beta 	2052	3260	+	402	15607144	887092	dnaN	Rv0002
  recombination protein F 	3280	4437	+	385	15607145	887089	recF	Rv0003
  hypothetical protein Rv0004 	4434	4997	+	187	15607146	887088	-	Rv0004
  DNA gyrase subunit B 	5123	7267	+	714	15607147	887081	gyrB	Rv0005
  DNA gyrase subunit A 	7302	9818	+	838	15607148	887105	gyrA	Rv0006
  ... 

  (full file has ~4000 lines)
  

.. NOTE::

  *It is crucial* that the annotation file (.prot_table) used for
  analyses in Transit corresponds to exactly the same genome sequence
  (.fasta or .fna) that was used to generate the .wig files with TPP,
  because it is used to determine which TA sites are contained in which
  genes (by coordinates). For example, H37Rv.fna is paired with
  H37Rv.prot_table, both derived from GenBank sequence NC_000962.3.


In many cases, users might often obtain annotations for their genome
in **.gff (or .gff3)** file format, such as downloaded from NCBI.  .gff
files contains essentially the same information about genes.  However,
there is a bit more flexibility in the .gff file format (especially in
the tags used in the right-most column), and the information about
genes is not always encoded in a uniform way, making it difficult to
use arbitrary .gff filess for analyses in Transit.  
Therefore, there is a
simple procedure in Transit to convert a .gff file to .prot_table
format (via GUI or command-line).  This
step only has to be done once, and then the .prot_table can be used
for all subsequent analyses in Transit.
(The routine specifically looks for the 'locus_tag', 'gene', and 'product'
tags in info field of CDS records.)

::

  > python3 transit.py convert gff_to_prot_table <input.gff_file> <output.prot_table>



.. _combined_wig:

Combined_wig Files
------------------

Transit now supports a new file format called 'combined_wig' which basically
combines multiple wig files into one file (with multiple columns).  This is
especially useful for managing larger collections of datasets.

Combined_wig files can created through the Transit GUI
(in the menu: **File -> Export -> Selected Datasets -> to_Combined_Wig**), or via the command line:

::

  Usage:

  python3 transit.py export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file> [-n <norm>]

  Example:

  > cd src/pytransit/data/
  > python3 ../../transit.py export combined_wig glycerol_H37Rv_rep1.wig,glycerol_H37Rv_rep2.wig,cholesterol_H37Rv_rep1.wig,cholesterol_H37Rv_rep2.wig,cholesterol_H37Rv_rep3.wig ../genomes/H37Rv.prot_table glyc_chol.combined_wig.txt


You can specify the normalization method you want to use with a flag.
TTR is the default, but other relevant normalization options would be 'nonorm'
(i.e. preserve raw counts) and 'betageom' (this corrects for skew, but is slow).

The format of a combined_wig is simply a multi-column (tab-separated) file with
the first column being the coordinates of TA sites, followed by 
N columns of counts (for N samples), possibly with a final column indicating
the gene annotation information.
A combined_wig file can have header lines, prefixed by '#'.

Importantly, a combined_wig file must include sample identifiers
(filenames) that are prefixed with **"#File: "**.  These header lines
are automatically included by the 'export combined_wig' Transit
command the creates a combined_wig from multiple .wig files.  These "#File:"
header lines should be given in the same order as the sample columns,
and the filenames are used as identifiers to cross-reference
information about each sample in the metadata file (see below).

Here is the output file from the example above (glyc_chol.combined_wig.txt):

::

 #Converted to CombinedWig with TRANSIT.
 #normalization method: TTR
 #Normalization Factors: 2.286197257382532 1.4371695843583265 2.6365107688705036 3.2329912278198036 2.865270203473558
 #RefGenome: H37Rv
 #File: glycerol_H37Rv_rep1.wig
 #File: glycerol_H37Rv_rep2.wig
 #File: cholesterol_H37Rv_rep1.wig
 #File: cholesterol_H37Rv_rep2.wig
 #File: cholesterol_H37Rv_rep3.wig
 #TA_coord	glycerol_H37Rv_rep1.wig	glycerol_H37Rv_rep2.wig	cholesterol_H37Rv_rep1.wig	cholesterol_H37Rv_rep2.wig	cholesterol_H37Rv_rep3.wig
 60	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 72	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 ...
 1345	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 1423	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 1522	0.0	178.2	0.0	0.0	0.0	Rv0001 (dnaA)
 1552	61.7	0.0	0.0	29.1	0.0	
 1635	32.0	288.9	142.4	22.6	0.0	
 ...

 (full file has ~75,000 lines)

(DnaA is essential, which is why there are no insertion counts in the first few TA sites)

Note that the ORF id and gene names are appended for TA sites in CDS regions, for convenience.
If you open a combined_wig file in Excel (with tab as separator character), you will
see the last line of the header (starting with "TAcoord...") provides the input wig filenames as column headers.

The '#RefGenome' is extracted from information in the individual .wig files,
which records the reference genome sequence that was used with TPP; 
the coordinates of TA sites are determined from the genome sequence.
The **reference genome must be the same for all wigs** being combined.


|


.. _metadata_files:

Samples Metadata File
---------------------

The **metadata** file describes the sample ids, filenames,
and conditions they represent (e.g. different media, growth
conditions, knockout strains, animal passaging, etc., or whatever
treatments and controls your TnSeq experiment involves) for all the
samples in a combined_wig file.  

Format of the *samples_metadata* file: a tab-separated file (which you
can edit in Excel) with 3 columns: Id, Condition, and Filename (it
must have these headers). The Condition column should be as specific
as possible, indicating **how to group replicates**.
You can include other columns of info, but
do not include additional rows.  Individual rows can be commented out
by prefixing them with a '#'.  Here is an example of a samples
metadata file: The filenames should match what is shown in the header
of the combined_wig (including pathnames, if present).

Note: the Condition column should have a unique label for each
distinct condition (the same label shared only among replicates).  If
there are attributes that distinguish the conditions (such as strain,
treatment, etc), they could be included as **additional columns**
(e.g. covariates, like Carbon_source or Drug or Days or Strain...).

::

  Id      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig




|

.. rst-class:: transit_sectionend
----
