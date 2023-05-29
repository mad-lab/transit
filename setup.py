"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Command
# To use a consistent encoding
from codecs import open
from os import path
from shutil import rmtree

import os

here = path.abspath(path.dirname(__file__))


# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Get current version
import sys
sys.path.insert(1, "src/")
import pytransit
version =  pytransit.__version__[1:]

class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def yes_or_no(self, question):
        while True:
            reply = str(raw_input(question +' (y/n): ')).lower().strip()
            if reply == 'y':
                return True
            return False

    def run(self):
        try:
            self.status('Removing previous builds...')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        if not self.yes_or_no("Have you done the following? \n" +
                    "- Updated README/Documentation?\n"
                    "- Are in the master branch, and have you merged version branch into master?\n"
                    "- Have you run the tests ('pytest tests/')?\n"
                    "- Have you updated CHANGELOG?\n"
                    "- Have you updated Transit Essentiality page?\n"
                    "- Updated version in src/pytransit/__init__.py (used to set git tag)?\n"
                    "- Is version v{0} correct".format(version)):
            self.status("Exiting...")
            sys.exit()

        self.status('Building Source and Wheel (universal) distribution...')
        os.system('{0} setup.py sdist bdist_wheel'.format(sys.executable))

        if self.yes_or_no("Add tag and push to public github? tag:v{0}".format(version)):
            self.status('Adding and pushing git tags to origin and public...')
            os.system('git tag v{0}'.format(version))
            os.system('git push origin --tags')
            os.system('git push https://github.com/mad-lab/transit master')
            os.system('git push https://github.com/mad-lab/transit --tags')
        else:
            self.status("Exiting...")
            sys.exit()

        if self.yes_or_no("Proceed with publish to PyPI? version: v{0}, tag:v{0}".format(version)):
            self.status('Uploading the package to PyPI via Twine...')
            os.system('twine upload dist/*')
        else:
            self.status("Exiting...")
            sys.exit()

        sys.exit()

setup(
    name='tnseq-transit',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version,

    description='TRANSIT is a tool for the analysis of Tn-Seq data. It provides an easy to use graphical interface and access to three different analysis methods that allow the user to determine essentiality in a single condition as well as between conditions.',
    long_description=long_description,
    long_description_content_type='text/markdown',

    # The project's main homepage.
    url='https://github.com/mad-lab/transit',
    download_url='https://github.com/mad-lab/transit',

    # Author details
    author='Michael A. DeJesus',
    author_email='mad@cs.tamu.edu',

    # Choose your license
    license='GNU GPL',
    python_requires='>=3.6',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers

    classifiers=[
        #'Development Status :: 3 - Alpha',
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    # What does your project relate to?
    keywords=['tnseq', 'analysis', 'biology', 'genome'],

    #package_dir = {'tnseq-transit': 'src/pytransit'},

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = find_packages('src', exclude=['contrib', 'tests']),
    #packages = ['pytransit'],
    package_dir = {'pytransit': 'src/pytransit',  'pytpp': 'src/pytpp'},
    include_package_data=True,
    #py_modules = ['tpp'],

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # 'pypubsub<4.0' and 'wxPython' are needed for GUI only, but go ahead and install them
    # the reason for restriction on pypubsub is that version>=4.0 does not work with python2 - I can probably get rid of this restriction, since everybody must be using python3 by now
    install_requires=['setuptools', 'numpy~=1.16', 'scipy~=1.2', 'matplotlib~=3.0', 'pillow', 'scikit-learn', 'statsmodels~=0.9', 'pypubsub<4.0', 'wxPython'],

    #dependency_links = [
    #	"git+https://github.com/wxWidgets/wxPython.git#egg=wxPython"
	#],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    #extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'pytransit': ['pytransit/data/*', 'pytransit/doc/*.*', 'pytransit/doc/images/*', 'pytransit/genomes/*']
    },

    #scripts=['src/tpp.py', 'src/transit.py'],

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('transitdata', ['package_data.dat'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'transit=pytransit.__main__:run_main',
            'tpp=pytpp.__main__:run_main',
        ],
    },
    cmdclass={
        'upload': UploadCommand,
    },
)

