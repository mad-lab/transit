


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

Use as a Python Package
-----------------------------------------------------


TRANSIT can be (optionally) installed as a python package. This can simplify the installation process as it will automatically install most of the requirements. In addition, it will allow users to use some of transit functions in their own scripts if they desire. Below is a brief example of importing transit functions into python. In this example, pair of .wig files are parsed into their read-counts (data) and genomic positions (position), and then normalization factors are calculated. See the documentation of the package for further examples:

::

        >>> import pytransit.norm_tools as norm_tools
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print(data)
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> factors = norm_tools.TTR_factors(data)
        >>> print(factors)
        array([[ 1.        ],
               [ 0.62862886]])


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



