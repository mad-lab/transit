


Installation
------------
TRANSIT can be downloaded from the public GitHub server,
`http://github.com/mad-lab/transit <http://github.com/mad-lab/transit>`_. It is released under a GPL
License. It can be downloaded with git as follows:

::

    
    
    git clone https://github.com/mad-lab/transit/
    

TRANSIT is python-based You must have python installed (installed by
default on most systems). In addition, TRANSIT relies on some python
packages/libraries/modules that you might need to install. Below are
the list of requirements:


|

Requirements
~~~~~~~~~~~~
The following libraries/modules are required to run TRANSIT:


+ `Python 2.7 <http://www.python.org>`_
+ `Numpy 1.6.1+ <http://www.numpy.org/>`_
+ `Scipy 0.14.0+ <http://www.scipy.org/>`_
+ `matplotlib 1.1.1+ <http://matplotlib.org/users/installing.html>`_
+ `wxpython 2.8.0+ <http://www.wxpython.org/>`_ (for Mac OSX, use the **cocoa** version of wxPython)
+ `PIL (Python Imaging Library) <http://www.pythonware.com/products/pil/>`_ or Pillow.


Generally, these requirements are install using the appropriate
methods for your operating system, i.e. apt-get or yum for unix
machines, pip or easy_install for OSX, or binary installers on
Windows. Below more detailed instructions are provided.

|

Detailed Instructions: Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_. Download the source files:

::

    
    `http://sourceforge.net/projects/bio-bwa/files/`_


Extract the files:

::

    
    tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::

    
    cd bwa-0.7.12
    make



|

Detailed Instructions: OSX
~~~~~~~~~~~~~~~~~~~~~~~~~~
First, download and install the latest Python 2.7.x installation file from the official python website:


    
    `http://www.python.org/downloads/ <http://www.python.org/downloads/>`_


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

::

    
    `http://downloads.sourceforge.net/wxpython/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg`_

Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_ . Download the source files:

::

    
    `http://sourceforge.net/projects/bio-bwa/files/`_


Extract the files:

::

    
    tar -xvjf bwa-0.7.12.tar.bz2


Go to the directory with the extracted source-code, and run make to create the executable files:

::

    
    cd bwa-0.7.12
    make




|

Detailed Instructions: Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, download and install the latest Python 2.7.x installation file
from the official python website:


::

    
    `http://www.python.org/downloads/`_


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


Optional: If you will be using the pre-processor, TPP, you will also need to install `BWA <http://bio-bwa.sourceforge.net/>`_. We provide a windows executable (.exe) for Windows 64 bit:

`bwa-0.7.12_windows.zip <http://saclab.tamu.edu/essentiality/transit/bwa-0.7.12_windows.zip>`_






