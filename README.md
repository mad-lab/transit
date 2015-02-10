This archive contains the files necessary to run the TRANSIT analysis
tool. The directory hierarchy is explained below.


For information and installation instructions, see
the manual by opening doc/transit.html in your browser.

For genomes and annotation files to use with TRANSIT, see
 http://saclab.tamu.edu/essentiality/transit/genomes/genomes.html


File directory:


* /  - Contains the TRANSIT python files

    - transit.py - main transit code.
    - transit_tools.py - global tools for reading files/data. Should be used by all methods.
    - gumbelMH.py - Gumbel method code.
    - MH_tools.py - tools used by gumbel method



* doc/ - Contains documentation (html format) for TRANSIT and TPP tools.
    - transit.html

    - TPP.html

* data/

    - Contains 5 input files (.wig format) from Griffin et al. 2011
      paper (http://www.ncbi.nlm.nih.gov/pubmed/21980284) on H37Rv
      grown on glycerol versus cholesterol.

    - Also contains represenative output files for the different
      analyses offered by TRANSIT.

* genomes/

    - Contains .fna and .prot_table files for 4 Mycobacterial organisms.