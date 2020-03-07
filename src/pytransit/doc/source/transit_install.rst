

.. _install-link:

Installation
============
TRANSIT can be downloaded from the public GitHub server,
`http://github.com/mad-lab/transit <http://github.com/mad-lab/transit>`_. It is released under a GPL
License. An archive with the lastest version fo the source code can be downloaded at the following link:


`Source code.zip <https://github.com/mad-lab/transit/archive/master.zip>`_



If you know how to utilize git, you can clone the git respository as follows:

::


    git clone https://github.com/mad-lab/transit/


TRANSIT is python-based You must have python installed (installed by
default on most systems). In addition, TRANSIT relies on some python
packages/libraries/modules that you might need to install (see `Requirements`_).

If you encounter problems, please :ref:`contact us <developers>` or head to the :ref:`Troubleshooting section <transit-troubleshoot>`.


|

Requirements
------------

TRANSIT runs on both python2.7 and python3. But the dependencies vary slightly.

Python 2.7:
-----------

The following libraries/modules are required to run TRANSIT:

+ `Python 2.7 <http://www.python.org>`_
+ `Numpy <http://www.numpy.org/>`_ (tested on 1.15.0)
+ `Statsmodels <https://pypi.org/project/statsmodels/>`_ (tested on 0.9.0)
+ `Scipy <http://www.scipy.org/>`_ (tested on 1.1)
+ `matplotlib <http://matplotlib.org/users/installing.html>`_ (tested on 2.2)
+ `Pillow 5.0 <https://github.com/python-pillow/Pillow>`_
+ `wxpython 4+ <http://www.wxpython.org/>`_
+ `PyPubSub 3.3 <https://pypi.org/project/PyPubSub/>`_ (Version 4.0 does not support python2 `Github Issue <https://github.com/schollii/pypubsub/issues/9>`_)

All of these dependencies can be installed using the following command.

::

   pip install numpy scipy pillow "pypubsub<4.0" "matplotlib<3.0" statsmodels wxPython

Pip and Python are usually preinstalled in most modern operating systems.

|

Python 3:
-----------

The following libraries/modules are required to run TRANSIT:

+ `Python 3+ <http://www.python.org>`_
+ `Numpy <http://www.numpy.org/>`_ (tested on 1.16.0)
+ `Statsmodels <https://pypi.org/project/statsmodels/>`_ (tested on 0.9.0)
+ `Scipy <http://www.scipy.org/>`_ (tested on 1.2)
+ `matplotlib <http://matplotlib.org/users/installing.html>`_ (tested on 3.0)
+ `Pillow 6.0 <https://github.com/python-pillow/Pillow>`_
+ `wxpython 4+ <http://www.wxpython.org/>`_
+ `PyPubSub 4+ <https://pypi.org/project/PyPubSub/>`_ (tested on 4.0.3)

All of these dependencies can be installed using the following command.

::

   pip3 install numpy scipy pillow pypubsub matplotlib statsmodels wxPython

Pip and Python are usually preinstalled in most modern operating systems.

|
.. _install-zinb:

Additional Requirements: R (statistical analysis package)
~~~~~~~~~~~~~~~~~~~~~~~~~~ 

R is called by Transit for certain commands, such as :ref:`ZINB <zinb>`, corrplot, and heatmap.
As of now, installing R is optional, and requires these additional steps...

Additional Installation Requirements for R:

 - install `R <https://www.r-project.org/>`_ (tested on v3.5.2)
 - R packages: **MASS, pscl, corrplot, gplots** (run "install.packages(MASS)" etc. in R console)
 - Python packages (for python3): rpy2 (v>=3.0) (run "pip3 install rpy2" on command line) 
 - Python packages (for python2.7): rpy2 (v<2.9.0) (run "pip install 'rpy2<2.9.0' " on command line)



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



You can use pip to install the TRANSIT package.

::

    sudo pip install tnseq-transit

This will automatically download and install TRANSIT as a package, and all remaining required python packages. Once TRANSIT is installed as a package, it can be executed as


.. NOTE::
   If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.

.. NOTE::
   The Transit package *does not* install wxPython. For graphical interface usage, this has to be done by the user. See :ref:`install wxPython <install-wxpython>`

|

Optional: Install BWA to use with TPP pre-processor
---------------------------------------------------

If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_.




.. _bwa-unix:

Linux & OSX Instructions
~~~~~~~~~~~~~~~~~~~~~~~~

Download the source files:


 + `http://sourceforge.net/projects/bio-bwa/files/ <http://sourceforge.net/projects/bio-bwa/files/>`_


Extract the files:

::


    tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::


    cd bwa-0.7.12
    make


.. _bwa-win:

Windows Instructions
~~~~~~~~~~~~~~~~~~~~

For Windows, we provide a windows executable (.exe) for Windows 64 bit:

  + `bwa-0.7.12_windows.zip <http://saclab.tamu.edu/essentiality/transit/bwa-0.7.12_windows.zip>`_



The 32-bit version of Windows is not recommended as it is limited in the amount of system memory that can be used.


|

.. _transit-upgrade:

Upgrading
---------

The process of upgrading transit will depend on how you installed transit initially.


Method 1: Upgrading package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If you installed TRANSIT as a package, then to upgrade, simply use pip to install tnseq-transit again, but this time include the '--upgrade' flag. For example:


::

    sudo pip install tnseq-transit --upgrade

This will automatically download and install the latest version of TRANSIT, as well as upgrade any of its requirements if necessary for compatability.


Method 2: Upgrading source installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed TRANSIT by downloading the raw source, then you can upgrade TRANSIT simply by replacing the old source code with the latest version. You can obtain a .zip archive with the latest version of the source through the following link:

https://github.com/mad-lab/transit/archive/master.zip

Simply exctract the code, and replace your existing files or delete the directory with the old source doe and use the newest version.

|

.. NOTE::
   If an an older version of wxPython is already installed (< 4.0), you may have to remove it and install version 4.0+.

|

.. _install-wxpython:

Installing wxPython
-------------------

wxPython 4+ can be installed using pip

::

   pip install wxPython

If the above command fails and you already have wxPython < 4.0 installed, you may have to manually remove it.
See https://stackoverflow.com/questions/50688630/cannot-uninstall-wxpython-3-0-2-0-macos for details.

.. _transit-troubleshoot:

Troubleshooting
---------------

1. No window appears when running in GUI mode.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This problem is likely due to running OSX and previously unsuported versions of matplotlib.
Please upgrade matplotlib to the latest version using:

::

    pip install 'matplotlib' --upgrade

|

2. pip: SystemError: Cannot compile 'Python.h'.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have the development libraries for python. You can fix this by installing the python-dev packages:


::

    sudo apt-get install python-dev


|

3. pip: "The following required packages can not be built: freetype,png," etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually matplotlib). Installing the following linux dependencies should fix this:

::

    sudo apt-get install libpng-dev libjpeg8-dev libfreetype6-dev


|

4. pip: "No lapack/blas resources found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually numpy/scipy). Installing the following linux dependencies should fix this:


::

    sudo apt-get install libblas-dev liblapack-dev libatlas-base-dev gfortran


|

5. "resources.ContextualVersionConflict (six 1.5.2)..."
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs some of the python modules are out of date. You can use pip to upgrade them as follows:


::

    sudo pip install six --upgrade
