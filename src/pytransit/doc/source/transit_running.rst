


Running TRANSIT
===============


TRANSIT can be run in either GUI mode or in console mode. GUI Mode will be run if TRANSIT is not given any command-line argumentsl. If any arguments are given, TRANSIT will run in console-mode.

The exact commands will vary depending on the method of installation. Details are given below:

|

GUI Mode
--------

In general, if you installed TRANSIT as a python package (e.g. using *pip install tnseq-transit*), then the proper way to run TRANSIT in GUI mode is simply to type the following into a console:

::

    transit


.. NOTE::
    In windows, you will likely have to navigate to C:\Python2.7\Scripts to be able to recognize the transit.exe file.


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

Gumbel
~~~~~~

To run the Gumbel analysis from the command line, type "python PATH/src/transit.py gumbel" followed by the following arguments:


+----------------+----------------+----------------+----------------+----------------+
| Argument       | Type           | Description    | Default        | Example        |
+================+================+================+================+================+
| annotation     | Required       | Path to        |                | genomes/H37Rv. |
|                |                | annotation     |                | prot\_table    |
|                |                | file in        |                |                |
|                |                | .prot\_table   |                |                |
|                |                | format         |                |                |
+----------------+----------------+----------------+----------------+----------------+
| control\_files | Required       | Comma-separate |                | data/glycerol\ |
|                |                | d              |                | _reads\_rep1.w |
|                |                | list of paths  |                | ig,data/glycer |
|                |                | to the \*.wig  |                | ol\_reads\_rep |
|                |                | replicate      |                | 2.wig          |
|                |                | datasets       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| output\_file   | Required       | Name of the    |                | results/gumbel |
|                |                | output file    |                | \_glycerol.dat |
|                |                | with the       |                |                |
|                |                | results.       |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -s SAMPLES     | Optional       | Number of      | 10000          | -s 20000       |
|                |                | samples to     |                |                |
|                |                | take.          |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -m MINREAD     | Optional       | Smallest       | 1              | -m 2           |
|                |                | read-count     |                |                |
|                |                | considered to  |                |                |
|                |                | be an          |                |                |
|                |                | insertion.     |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -b BURNIN      | Optional       | Burn in        | 500            | -b 100         |
|                |                | period, Skips  |                |                |
|                |                | this number of |                |                |
|                |                | samples before |                |                |
|                |                | getting        |                |                |
|                |                | estimates. See |                |                |
|                |                | documentation. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -t TRIM        | Optional       | Number of      | 1              | -t 2           |
|                |                | samples to     |                |                |
|                |                | trim. See      |                |                |
|                |                | documentation. |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -r REP         | Optional       | How to handle  | Sum            | -r Mean        |
|                |                | replicates     |                |                |
|                |                | read-counts:   |                |                |
|                |                | 'Sum' or       |                |                |
|                |                | 'Mean'.        |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iN IGNOREN    | Optional       | Ignore TAs     | 5              | -iN 0          |
|                |                | occuring at X% |                |                |
|                |                | of the N       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+
| -iC IGNOREC    | Optional       | Ignore TAs     | 5              | -iC 10         |
|                |                | occuring at X% |                |                |
|                |                | of the C       |                |                |
|                |                | terminus.      |                |                |
+----------------+----------------+----------------+----------------+----------------+



::

    python PATH/src/transit.py gumbel genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig test_console_gumbel.dat -s 20000 -b 1000




|

Tn5 Gaps
~~~~~~~~

To run the Tn5 Gaps analysis from the command line, type "python
PATH/src/transit.py tn5gaps" followed by the following arguments:

Argument Type Description Default Example annotation Required Path to
annotation file in .prot_table format genomes/Salmonella-
Ty2.prot_table control_files Required Comma-separated list of paths to
the \*.wig replicate datasets
data/salmonella_2122_rep1.wig,data/salmonella_2122_rep2.wig
output_file Required Name of the output file with the results.
results/test_console_tn5gaps.dat -m MINREAD Optional Smallest read-
count considered to be an insertion. 1 -m 2 -r REP Optional How to
handle replicates read-counts: 'Sum' or 'Mean'. Sum -r Sum

Example Tn5 Gaps command:

::

    python PATH/src/transit.py tn5gaps genomes/Salmonella-Ty2.prot_table data/salmonella_2122_rep1.wig,data/salmonella_2122_rep2.wig results/test_console_tn5gaps.dat -m 2 -r Sum





Example HMM command:

::

    python PATH/src/transit.py hmm genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig test_console_hmm.dat -r Sum


| 

Resampling
~~~~~~~~~~

To run the Resampling analysis from the command line, type "python
PATH/src/transit.py resampling" followed by the following arguments:

+------------------+----------------+--------------------------+----------------+----------------+
| Argument         | Type           | Description              | Default        | Example        |
+==================+================+==========================+================+================+
| annotation       | Required       | Path to                  |                | genomes/H37Rv. |
|                  |                | annotation               |                | prot\_table    |
|                  |                | file in                  |                |                |
|                  |                | .prot\_table             |                |                |
|                  |                | format                   |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| control\_files   | Required       | Comma-separate           |                | data/glycerol\ |
|                  |                | d                        |                | _reads\_rep1.w |
|                  |                | list of paths            |                | ig,data/glycer |
|                  |                | to the \*.wig            |                | ol\_reads\_rep |
|                  |                | replicate                |                | 2.wig          |
|                  |                | datasets for             |                |                |
|                  |                | the control              |                |                |
|                  |                | condition                |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| exp\_files       | Required       | Comma-separate           |                | data/cholester |
|                  |                | d                        |                | ol\_reads\_rep |
|                  |                | list of paths            |                | 1.wig,data/cho |
|                  |                | to the \*.wig            |                | lesterol\_read |
|                  |                | replicate                |                | s\_rep2.wig    |
|                  |                | datasets for             |                |                |
|                  |                | the                      |                |                |
|                  |                | experimental             |                |                |
|                  |                | condition                |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| output\_file     | Required       | Name of the              |                | results/gumbel |
|                  |                | output file              |                | \_glycerol.dat |
|                  |                | with the                 |                |                |
|                  |                | results.                 |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -s SAMPLES       | Optional       | Number of                | 10000          | -s 5000        |
|                  |                | permutations             |                |                |
|                  |                | performed.               |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -h               | Optional       | Creates                  | Not set        | -H             |
|                  |                | histograms of            |                |                |
|                  |                | the                      |                |                |
|                  |                | permutations             |                |                |
|                  |                | for all genes.           |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -a               | Optional       | Performs                 | Not set        | -a             |
|                  |                | adaptive                 |                |                |
|                  |                | appoximation             |                |                |
|                  |                | to resampling.           |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -N               | Optional       | Select which             | nzmean         | -N nzmean      |
|                  |                | normalizing              |                |                |
|                  |                | procedure to             |                |                |
|                  |                | use. Can                 |                |                |
|                  |                | choose between           |                |                |
|                  |                | 'TTR',                   |                |                |
|                  |                | 'nzmean',                |                |                |
|                  |                | 'totreads',              |                |                |
|                  |                | 'zinfnb',                |                |                |
|                  |                | 'betageom',              |                |                |
|                  |                | and 'nonorm'.            |                |                |
|                  |                | See the                  |                |                |
|                  |                | parameters               |                |                |
|                  |                | section for              |                |                |
|                  |                | the                      |                |                |
|                  |                | `Re-sampling             |                |                |
|                  |                | method <http:/           |                |                |
|                  |                | /saclab.tamu.e           |                |                |
|                  |                | du/essentialit           |                |                |
|                  |                | y/transit/tran           |                |                |
|                  |                | sit.html#resam           |                |                |
|                  |                | pling>`__                |                |                |
|                  |                | for a                    |                |                |
|                  |                | description of           |                |                |
|                  |                | these                    |                |                |
|                  |                | normalization            |                |                |
|                  |                | options.                 |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -iN IGNOREN      | Optional       | Ignore TAs               | 5              | -iN 0          |
|                  |                | occuring at X%           |                |                |
|                  |                | of the N                 |                |                |
|                  |                | terminus.                |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| -iC IGNOREC      | Optional       | Ignore TAs               | 5              | -iC 10         |
|                  |                | occuring at X%           |                |                |
|                  |                | of the C                 |                |                |
|                  |                | terminus.                |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| --ctrl_lib <str> | Optional       | String of letters (ABC)  | Empty          | -ctrl_lib AB   |
|                  |                | representing the         |                |                |
|                  |                | libraries of the control |                |                |
|                  |                | datasets being used.     |                |                |
|                  |                | If set, permutatations   |                |                | 
|                  |                | will be within libraries |                |                |
+------------------+----------------+--------------------------+----------------+----------------+
| --exp_lib <str>  | Optional       | String of letters (ABC)  | Empty          | -exp_lib ABB   |
|                  |                | representing the         |                |                |
|                  |                | libraries of the experi- |                |                |
|                  |                | mental datasets used.    |                |                |  
|                  |                | If set, permutatations   |                |                | 
|                  |                | will be within libraries |                |                |
+------------------+----------------+--------------------------+----------------+----------------+





Example Resampling command:

::

    python PATH/src/transit.py resampling genomes/H37Rv.prot_table data/glycerol_reads_rep1.wig,data/glycerol_reads_rep2.wig data/cholesterol_reads_rep1.wig,data/cholesterol_reads_rep2.wig,data/cholesterol_reads_rep3.wig test_console_resampling.dat -h -s 10000 -n nzmean


