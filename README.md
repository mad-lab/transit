# TRANSIT


This repository contains the code and files necessary to run TRANIST and TPP. For more information see:

http://saclab.tamu.edu/essentiality/transit/transit.html

http://saclab.tamu.edu/tom/TPP.html


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