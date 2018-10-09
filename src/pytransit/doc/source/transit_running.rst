


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

See example usages of supported methods in :ref:`Analysis Methods <analysis_methods>` section.
