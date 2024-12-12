

.. _install-link:

Installation
============
TRANSIT can be downloaded from the public GitHub server,
`http://github.com/ioerger/transit <http://github.com/ioerger/transit>`_. It is released under a GPL
License. An archive with the lastest version fo the source code can be downloaded at the following link:


`Source code.zip <https://github.com/ioerger/transit/archive/master.zip>`_



If you know how to utilize git, you can clone the git respository as follows:

::


    git clone https://github.com/ioerger/transit/


TRANSIT is python-based. You must have python installed (installed by
default on most systems). In addition, TRANSIT relies on some python
packages/libraries/modules that you might need to install (see `Requirements`_).

If you encounter problems, please :ref:`contact us <developers>` or head to the :ref:`Troubleshooting section <transit-troubleshoot>`.


|


Installing Transit (and its Dependencies) using pip
---------------------------------------------------

You can use pip to install the TRANSIT package.

::

    sudo pip install transit1

This will automatically download and install TRANSIT as a package, and all remaining required python packages. Once TRANSIT is installed as a package, it can be executed at the command-line as 'transit', which should pop-up the GUI.


.. NOTE::
   In Dec 2024, we switched the package name on PyPi from 'tnseq-transit' to 'transit1'.


.. NOTE::
   If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.


.. NOTE::
   The Transit package *does not* install R, which would have to be done manually.  R is not required for all of Transit, just certain methods (like ZINB).



Requirements
------------

Starting with release v3.0, TRANSIT now requires python3. 

TRANSIT runs on both python2.7 and python3. But the dependencies vary slightly.

To use TRANSIT with python2, use a TRANSIT release prior to 3.0 (e.g. v2.5.2)

.. NOTE::
  **Note about Python3.12**: The update to Python3.12 made changes in the setuptools package that 
  disrupted the installation of Transit as well as other python packages (including some Transit
  dependencies, like wxPython).
  
  We are working on fixing these problems caused by 3.12.  In the meantime,
  if you have difficulty installing Transit with Python3.12, we recommend trying it with Python3.11, with
  which Transit should work fine.

|


Python 3:
~~~~~~~~~

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

Python 2.7:
~~~~~~~~~~~

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

Install BWA to use with TPP pre-processor (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _install-wxpython:

Installing wxPython
~~~~~~~~~~~~~~~~

wxPython 4+ can be installed using pip

::

   pip install wxPython

If the above command fails and you already have wxPython < 4.0 installed, you may have to manually remove it.
See https://stackoverflow.com/questions/50688630/cannot-uninstall-wxpython-3-0-2-0-macos for details.

.. NOTE::

  Installing *wxPython* can be a bit finicky.  It might require installing the
  development version of GTK first.  There are at least two versions currently, 
  *gtk2* and *gtk3*.
  Transit should work with both, although there can be small differences in the 
  visual look of the GUI.  To get *wxPython* to install, you might try doing this:

    > sudo apt-get install libgtk-2-dev

    or

    > sudo apt-get install libgtk-3-dev

  depending on which version of *libgtk* you have installed.

|

Windows Instructions
~~~~~~~~~~~~~~~~~~~~

For Windows, we provide a windows executable (.exe) for Windows 64 bit:

  + `bwa-0.7.12_windows.zip <http://saclab.tamu.edu/essentiality/transit/bwa-0.7.12_windows.zip>`_



The 32-bit version of Windows is not recommended as it is limited in the amount of system memory that can be used.


|

.. _install-zinb:

Installing R (statistical analysis package) (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~ 

R is called by Transit for certain commands, such as :ref:`ZINB <zinb>`, corrplot, and heatmap.
As of now, installing R is optional, and requires these additional steps...

Additional Installation Requirements for R:

 - install `R <https://www.r-project.org/>`_ (tested on v3.5.2)
 - R packages: **MASS, pscl, corrplot, gplots** (run "install.packages(MASS)" etc. in R console)
 - Python packages (for python3): rpy2 (v>=3.0) (run "pip3 install rpy2" on command line) 
 - Python packages (for python2.7): rpy2 (v<2.9.0) (run "pip install 'rpy2<2.9.0' " on command line)


|


.. _transit-upgrade:

Upgrading
---------

The process of upgrading transit will depend on how you installed transit initially.


Method 1: Upgrading package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If you installed TRANSIT as a package, then to upgrade, simply use pip to install 'transit1' again, but this time include the '--upgrade' flag. For example:


::

    sudo pip install transit1 --upgrade

This will automatically download and install the latest version of TRANSIT, as well as upgrade any of its requirements if necessary for compatability.


.. NOTE::
   In Dec 2024, we switched the package name on PyPi from 'tnseq-transit' to 'transit1'. Hopefully, users who previously installed transit using 'tnseq-transit' should also be able to upgrade, and it will automatically upgrade to transit1.  You might have to use this command: "sudo pip install tnseq-transit --upgrade' for your package name.



Method 2: Upgrading source installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed TRANSIT by downloading the raw source, then you can upgrade TRANSIT simply by replacing the old source code with the latest version. You can obtain a .zip archive with the latest version of the source through the following link:

https://github.com/ioerger/transit/archive/master.zip

Simply exctract the code, and replace your existing files or delete the directory with the old source doe and use the newest version.

|

.. NOTE::
   If an an older version of wxPython is already installed (< 4.0), you may have to remove it and install version 4.0+.

|

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
