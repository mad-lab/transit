


Running TRANSIT
===============


TRANSIT can be run in either GUI mode or in console mode. GUI Mode
will be run if TRANSIT is not given any command-line arguments. If any
arguments are given, TRANSIT will run in console-mode.

The exact commands will vary depending on the method of
installation. Details are given below:

|

GUI Mode
--------

In general, if you installed TRANSIT as a python package (e.g. using
*pip install tnseq-transit*), then the proper way to run TRANSIT in
GUI mode is simply to type the following into a console:

::

    transit


.. NOTE::
    In windows, you will likely have to navigate to C:\\Python2.7\\Scripts to be able to recognize the transit.exe file.



If, however, you installed transit by downloading and extracting the source-code archive, you can run TRANSIT in GUI mode by typing in the command line:

::

    python PATH/src/transit.py

where PATH is the path to the TRANSIT installation directory. You might be able to double-click on icon for transit.py, if your OS associates .py files with python and automatically runs them.


.. NOTE::
    Note, because TRANSIT has a graphical user interface, if you are trying to run TRANSIT in GUI mode across a network, for example by running on a unix server but displaying on a desktop machine, you will probably need to use 'ssh -Y' and a local X11 client (like Xming or Cygwin/X on PCs). This will allow the GUI component to be properly displayed accross the network connection.


|

Command line Mode
-----------------
TRANSIT can also be run purely the command line, without a GUI interface. This is convenient if you want to run many analyses in batch, as you can write a script that automatically runs several analyses in parallel or in sequence

If you installed TRANSIT as a python package, you can get a list of possible arguments by typing:


::

    transit -h


Or if you installed it by downloading and extracting an archive with the source code:

::

    python PATH/src/transit.py -h



In most cases TRANSIT expects the user to specify which analysis method they wish to run as their first argument. The user will need to type the short-name of the analysis method desired, e.g. "gumbel", "hmm", or "resampling". By choosing a method, and adding the "-h" flag, you will get a list of all the necessary parameters and optional flags for the chosen method.


If you installed TRANSIT as a python package, you can achieve this by typing:


::

    transit gumbel -h


Or if you installed it by downloading and extracting an archive with the source code:

::

    python PATH/src/transit.py gumbel -h


|

See example usages of supported methods in :ref:`Analysis Methods <analysis_methods>` section.

|

Prot_tables (Annotations)
-------------------------

Most of the methods in Transit use a custom format for genome annotations called a '.prot_table'.
It is a simple tab-separated text file with specific columns, as originally defined for genomes
in Genbank many years ago.

The required columns are:

1. gene function description
2. start coordinate
3. end coordinate
4. strand
5. length of protein product (in amino acids)
6. don't care
7. don't care
8. gene name (like "dnaA")
9. ORF id (like Rv0001)

It is crucial to use the same .prot_table corresponding to the genome sequence that was
used to generate the wig file (count insertions) by TPP.  This is because the
coordinates of TA sites in the wig file and the coordinates of ORF boundaries
must use the same coordinate system (which can be thrown out of register by indels).

Suppose you have a .prot_table for genome A, and you want to map reads to 
another genome B which is closely related, but for which you do not have an annotation.
You can use the following web-app ( `Prot_table Adjustment Tool <http://saclab.tamu.edu/cgi-bin/iutils/app.cgi>`_ ) 
to convert the annotation for A to B
by adjusting all the coordinates of ORFs from A to B according to a genome alignment.
For example, you could use this to map known ORFs in H37Rv to sequences of other strains, like HN878 or CDC1551.
(Even though they have their own annotations, it might be helpful to use the genes as defined in H37Rv)

While some Transit methods can also work with .gff (or .gff3) files,
the flexibility of the .gff format makes it difficult to anticipate all possible encoding schemes.
Therefore, to simplify things, we recommend you convert your .gff file to .prot_table format
once at the beginning and then use that for all work with Transit,
which can be done through the GUI (under 'Convert' in menu), or on the command-line as follows:


::

  > python transit.py convert gff_to_prot_table <.gff> <.prot_table>

|


.. _tn5-main-overview:

Tn5 Datasets
------------

Transit can now process and analyze Tn5 datasets  This is a different transposon than Himar1.
The major difference is Tn5 can insert at any site in the genome, and is not restricted
to TA dinucleotides (and saturation is typically much lower).  This affects
the statistical analyses (which were originally designed for Himar1 and can't directly
be applied to Tn5). Therefore, :ref:`Resampling <resampling>` was extended to handle Tn5 for comparative analysis, and
:ref:`Tn5Gaps <tn5gaps>` is a new statistical model for identifying essential genes in single Tn5 datasets.
Amplification of Tn5 libraries
uses different primers, and this affects the pre-processing by TPP.  But TPP has
be modified to recognize the primer sequence for the most widely
used protocol for Tn5.  Furthermore, TPP now has an option for users to define their
own primer sequences, if they use a different sample prep protocol.



