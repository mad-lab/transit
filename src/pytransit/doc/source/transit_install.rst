

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
The following libraries/modules are required to run TRANSIT:


+ `Python 2.7 <http://www.python.org>`_
+ `Numpy <http://www.numpy.org/>`_ (tested on 1.13.0)
+ `Scipy <http://www.scipy.org/>`_ (tested on 0.19.1)
+ `matplotlib <http://matplotlib.org/users/installing.html>`_ (tested on 2.0.2)
+ `wxpython 2.8.0+ <http://www.wxpython.org/>`_ (for Mac OSX, use the **cocoa** version of wxPython; If using El Capitan, please see :ref:`OSX El Capitan notice <osx_el_capitan>` for special instructions)
+ `PIL (Python Imaging Library) <http://www.pythonware.com/products/pil/>`_ or Pillow.


Generally, these requirements are install using the appropriate
methods for your operating system, i.e. apt-get or yum for unix
machines, pip or easy_install for OSX, or binary installers on
Windows. Below more detailed instructions are provided.

|



Use as a Python Package
-----------------------------------------------------


TRANSIT can be (optionally) installed as a python package. This can simplify the installation process as it will automatically install most of the requirements. In addition, it will allow users to use some of transit functions in their own scripts if they desire. Below is a brief example of importing transit functions into python. In this example, pair of .wig files are parsed into their read-counts (data) and genomic positions (position), and then normalization factors are calculated. See the documentation of the package for further examples:

    :Example:
        >>> import pytransit.norm_tools as norm_tools
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> factors = norm_tools.TTR_factors(data)
        >>> print factors
        array([[ 1.        ],
               [ 0.62862886]])

    .. seealso:: :class:`transit`






Detailed Instructions: Linux
----------------------------


Method 1: Install as a Python Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Most of the requirements are available in default package sources in
most Linux distributions. The following commands will install python, pip, wxPython, and several other dependencies needed by the python modules:


::

    #Ubuntu:
    sudo apt-get install python python-dev python-pip pkg-config python-wxgtk3.0 libpng-dev libjpeg8-dev libfreetype6-dev

    #Fedora:
    sudo yum install python python-dev python-pip pkg-config python-wxgtk3.0 libpng-dev libjpeg8-dev libfreetype6-dev


Finally you can use pip to install the TRANSIT package:

::

    sudo pip install tnseq-transit

This will automatically download and install TRANSIT as a package, and all remaining required python packages. Once TRANSIT is installed as a package, it can be executed as


.. NOTE::
   If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.


|


Method 2: Install Source Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Most of the requirements are available in default package sources in
most Linux distributions. The following commands will install python,
numpy, scipy, matplotlib on the Ubuntu or Fedora Linux distributions:

::


    #Ubuntu:
    sudo apt-get install python python-numpy python-scipy python-matplotlib python-wxgtk3.0

    #Fedora:
    sudo yum install python numpy scipy python-matplotlib python-wxgtk3.0


The final requirement left to install is Pillow. First you need
install pip which simplifies the process of installing certain python
modules like Pillow:


::


    #Ubuntu:
    sudo apt-get install pip

    #Fedora:
    sudo yum install pip


Next, using pip you must have a clean installation of Pillow, and the
desired libraries. You can achieve this through the following
commands:

::


    #Ubuntu:
    pip uninstall pillow
    pip uninstall Pillow
    sudo apt-get install libjpeg-dev zlib1g-dev
    pip install -I Pillow

    #Fedora:
    pip uninstall pillow
    pip uninstall Pillow
    sudo yum install install libjpeg-dev zlib1g-dev
    pip install -I Pillow




After all the requirements have been installed, you can download an archive with the lastest version of the TRANSIT source code through the following link:


`Source code.zip <https://github.com/mad-lab/transit/archive/master.zip>`_


Or, if you prefer to utilize git, you can clone the git respository as follows:

::

    git clone https://github.com/mad-lab/transit/



.. NOTE::
       If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.



|

Detailed Instructions: OSX
--------------------------





Method 1: Install as a Python Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


First, download and install the latest Python 2.7.x installation file from the official python website:



 + `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


Next make sure you have pip installed. Pip can be installed through easy_install, which should come with OSX:

::


    sudo easy_install pip


Download and install the OSX binary of wxpython (cocoa version) for python 2.7:

 + `http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg <http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg>`_

.. _osx_el_capitan:

.. NOTE::
   If you are running OSX El Capitan or later, you will need to use a repackaged version of the
   wxpython installer. You can `download a repackaged version from our servers <http://orca1.tamu.edu/essentiality/transit/wxPython3.0-osx-cocoa-py2.7_mad_elcapitan.pkg>`_ or you can follow `these detailed instructions to repackage the installer <http://davixx.fr/blog/2016/01/25/wxpython-on-os-x-el-capitan/>`_ if you prefer.




Finally you can use pip to install the TRANSIT package:


::

    sudo pip install tnseq-transit

This will automatically download and install TRANSIT and all remaining requirements.


|


.. NOTE::
   If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.

|


Method 2: Install Source Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download and install the latest Python 2.7.x installation file from the official python website:


 + `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


Next make sure you have pip installed. Pip can be installed through easy_install, which should come with OSX:

::


    sudo easy_install pip


Next install numpy, scipy, and matplotlib and pillow using pip:

::


    sudo pip install numpy
    sudo pip install scipy
    sudo pip install matplotlib
    sudo pip install pillow


Download and install the OSX binary of wxpython (cocoa version) for python 2.7:


 + `http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg <http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg>`_

.. NOTE::
   If you are running OSX El Capitan or later, you will need to use a repackaged version of the
   wxpython installer. You can `download a repackaged version from our servers <http://orca1.tamu.edu/essentiality/transit/wxPython3.0-osx-cocoa-py2.7_mad_elcapitan.pkg>`_ or you can follow `these detailed instructions to repackage the installer <http://davixx.fr/blog/2016/01/25/wxpython-on-os-x-el-capitan/>`_ if you prefer.



After all the requirements have been installed, you can download an archive with the lastest version of the TRANSIT source code through the following link:


`Source code.zip <https://github.com/mad-lab/transit/archive/master.zip>`_


Or, if you prefer to utilize git, you can clone the git respository as follows:

::

    git clone https://github.com/mad-lab/transit/



.. NOTE::
     If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-unix>`.



|

Detailed Instructions: Windows
------------------------------


Method 1: Install as a Python Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download and install the latest Python 2.7.x installation file
from the official python website:


 + `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


Next, you will need to install pip. If you are using python 2.7.9+
then pip will come pre-installed and included in the default script
directory (i.e. C:\Python27\Scripts ). If you are using python 2.7.8
or older, you will need to manually install pip by downloading and
running the `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`_ script:


::


    python.exe get-pip.py


Make sure that "wheel" is installed. This is necessary to allow you to
install .whl (wheel) files:

::

    pip.exe install wheel

Next install the transit package using pip:

::

    pip.exe install tnseq-transit




To use transit in GUI mode you will need to install wxPython versions 3.0 or earlier. We have provided .whl files which you can download and install below. (Note: Make sure to
choose the files that match your Windows version i.e. 32/64 bit)


  + `wxPython-3.0.2.0-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win32.whl>`_


  + `wxPython_common-3.0.2.0-py2-none-any.whl <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_


Finally, install the files using pip:

::


    pip.exe install wxPython-3.0.2.0-cp27-none-win_amd64.whl
    pip.exe install wxPython_common-3.0.2.0-py2-none-any.whl


making sure to replace the name with the file you downloaded (i.e. 32bit vs 64 bit)


.. NOTE::
    If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-win>`.




Method 2: Install Source Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download and install the latest Python 2.7.x installation file
from the official python website:


 + `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


Next, you will need to install pip. If you are using python 2.7.9+
then pip will come pre-installed and included in the default script
directory (i.e. C:\Python27\Scripts ). If you are using python 2.7.8
or older, you will need to manually install pip by downloading and
running the `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`_ script:


::


    python.exe get-pip.py


Make sure that "wheel" is installed. This is necessary to allow you to
install .whl (wheel) files:

::


    pip.exe install wheel


Download the .whl files for all the requirements (Note: Make sure to
choose the files that match your Windows version i.e. 32/64 bit)

  + `numpy-1.9.2+mkl-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/numpy-1.9.2+mkl-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/numpy-1.9.2+mkl-cp27-none-win32.whl>`_


  + `scipy-0.15.1-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/scipy-0.15.1-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/scipy-0.15.1-cp27-none-win32.whl>`_


  + `matplotlib-1.4.3-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/matplotlib-1.4.3-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/matplotlib-1.4.3-cp27-none-win32.whl>`_


  + `Pillow-2.8.2-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/Pillow-2.8.2-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/Pillow-2.8.2-cp27-none-win32.whl>`_


  + `wxPython-3.0.2.0-cp27-none-win_amd64.whl <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win_amd64.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython-3.0.2.0-cp27-none-win32.whl>`_


  + `wxPython_common-3.0.2.0-py2-none-any.whl <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_ or `[32 bit] <http://saclab.tamu.edu/essentiality/transit/wxPython_common-3.0.2.0-py2-none-any.whl>`_






Source: These files were obtained from the `Unofficial Windows Binaries for Python Extension Packages by Christoph Gohlke, Laboratory for Fluorescence Dynamics, University of California, Irvine. <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_


Finally, install the files using pip:

::


    pip.exe install numpy-1.9.2+mkl-cp27-none-win_amd64.whl
    pip.exe install scipy-0.15.1-cp27-none-win_amd64.whl
    pip.exe install matplotlib-1.4.3-cp27-none-win_amd64.whl
    pip.exe install Pillow-2.8.1-cp27-none-win_amd64.whl
    pip.exe install wxPython-3.0.2.0-cp27-none-win_amd64.whl
    pip.exe install wxPython_common-3.0.2.0-py2-none-any.whl



After all the requirements have been installed, you can download an archive with the lastest version of the TRANSIT source code through the following link:


`Source code.zip <https://github.com/mad-lab/transit/archive/master.zip>`_


Or, if you prefer to utilize git, you can clone the git respository as follows:

::

    git clone https://github.com/mad-lab/transit/




.. NOTE::
       If you will be using the pre-processor, TPP, you will also need to install :ref:`install BWA <bwa-win>`.


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

|
----

|

|

.. _transit-upgrade:

Upgrading
---------

The process of upgrading transit will depend on how you installed transit initially.



Method 1: Upgrading package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If you installed TRANSIT as a package, then to upgrade, simply use pip to install tnseq-transit again, but this time include the '--upgrade' flag. For example:


::

    sudo pip install tnseq-transit

This will automatically download and install the latest version of TRANSIT, as well as upgrade any of its requirements if necessary for compatability.




Method 2: Upgrading source installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed TRANSIT by downloading the raw source, then you can upgrade TRANSIT simply by replacing the old source code with the latest version. You can obtain a .zip archive with the latest version of the source through the following link:

https://github.com/mad-lab/transit/archive/master.zip

Simply exctract the code, and replace your existing files or delete the directory with the old source doe and use the newest version.


|

|
----


|

|



.. _transit-troubleshoot:

Troubleshooting
---------------

|

1. Gtk-ERROR \*\*: GTK+ 2.x symbols detected
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This error can occur if you have GTK2 already installed and then install wxPython version 3.0+. To fix this, please try installing version 2.8 of wxPython or install a new version of GTK3. More information on this error to come. 


|

2. wxPython & OSX: "The Installer could not install the software because there was no software found to install."
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are running OSX El Capitan or later, you will need to use a repackaged
version of the wxpython installer as OSX El Capitan has removed support for older packaging methods still used by wxPython. You can `download a repackaged version of wxPython
from our servers <http://orca1.tamu.edu/essentiality/transit/wxPython3.0-osx-cocoa-py2.7_mad_elcapitan.pkg>`_ or you can follow `these detailed instructions to repackage the installer <http://davixx.fr/blog/2016/01/25/wxpython-on-os-x-el-capitan/>`_ if you prefer.


|

3. No window appears when running in GUI mode.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This problem is likely due to running OSX and previously unsuported versions of matplotlib.
Please upgrade matplotlib to the latest version using:

::

    pip install 'matplotlib' --upgrade


|

4. Unable to locate package python-wxgtk3.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your version of Linux might not have the repository address that includes python-wxgtk3.0. You can attempt to install version 2.8 instead:

::

    sudo apt-get install python-wxgtk2.8



or you can add the repository that includes version 3.0 and install it:

::

    # Add repo for 14.04
    sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu utopic main restricted universe"

    #Update repo information
    sudo apt-get update

    #Install wxPython 3.0
    sudo apt-get install python-wxgtk3.0

    #Remove repo to prevent version conflicts
    sudo add-apt-repository --remove "deb http://archive.ubuntu.com/ubuntu utopic main restricted universe"


|

5. pip: SystemError: Cannot compile 'Python.h'.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have the development libraries for python. You can fix this by installing the python-dev packages:


::

    sudo apt-get install python-dev


|

6. pip: "The following required packages can not be built: freetype,png," etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually matplotlib). Installing the following linux dependencies should fix this:

::

    sudo apt-get install libpng-dev libjpeg8-dev libfreetype6-dev


|

7. pip: "No lapack/blas resources found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs when you do not have some dependencies that are necessary to build some of the python modules TRANSIT requires (usually numpy/scipy). Installing the following linux dependencies should fix this:


::

    sudo apt-get install libblas-dev liblapack-dev libatlas-base-dev gfortran


|

8. "resources.ContextualVersionConflict (six 1.5.2)..."
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This occurs some of the python modules are out of date. You can use pip to upgrade them as follows:


::

    sudo pip install six --upgrade
