
.. _cheat-link:

Console Mode Cheat-Sheet
========================


TRANSIT has the capability of running in Console mode, without 
depending on libraries for GUI elements. More hands-on users
can utilize transit in this manner to quickly run multiple
jobs in parallel. Below is brief 

|


Analysis Methods
++++++++++++++++
TRANSIT has the capacity of determining essentiality within a
single condition, or between conditions to determine
conditional essentiality. 



Single Condition Essentiality 
-----------------------------

Analysis methods in a single condition require at least 4
positional arguments followed by optional flags.


.. code-block:: python

   python transit.py <method> <wig-files> <annotation> <output>


+----------------------+-------------------------------------------------------------------------+
| Positional Arguments | Definition                                                              |
+======================+=========================================================================+
| <method>             | Short name of the desired analysis method  e.g. gumbel, resampling, hmm |
+----------------------+-------------------------------------------------------------------------+
| <wig-files>          | Comma-separated list of paths read-count datasets in .wig format        |
+----------------------+-------------------------------------------------------------------------+
| <annotation>         | Path to the annotation in .prot_table or .GFF3 format.                  |
+----------------------+-------------------------------------------------------------------------+
| <output>             | Desired path and name of the output file                                |
+----------------------+-------------------------------------------------------------------------+

Example
~~~~~~~

.. code-block:: python

   python transit.py gumbel glycerol_H37Rv_rep1.wig,glycerol_H37Rv_rep2.wig H37Rv.prot_table glycerol_TTR.txt -r Sum -s 10000



Conditional Essentiality
------------------------

Analysis methods between two conditions require at least 5
positional arguments followed by optional flags.



+----------------------+----------------------------------------------------------------------------------------------+
| Positional Arguments | Definition                                                                                   |
+======================+==============================================================================================+
| <method>             | Short name of the desired analysis method  e.g. gumbel, resampling, hmm                      |
+----------------------+----------------------------------------------------------------------------------------------+
| <control-files>      | Comma-separated list of paths read-count files in .wig format for the control datasets       |
+----------------------+----------------------------------------------------------------------------------------------+
| <experimental-files> | Comma-separated list of paths read-count files in .wif format for the experimental datasets  |
+----------------------+----------------------------------------------------------------------------------------------+
| <annotation>         | Path to the annotation in .prot_table or .GFF3 format.                                       |
+----------------------+----------------------------------------------------------------------------------------------+
| <output>             | Desired path and name of the output file                                                     |
+----------------------+----------------------------------------------------------------------------------------------+

|


Example
~~~~~~~

.. code-block:: python

   python transit.py resampling glycerol_H37Rv_rep1.wig,glycerol_H37Rv_rep2.wig cholesterol_H37Rv_rep1.wig,cholesterol_H37Rv_rep2.wig H37Rv.prot_table glycerol_TTR.txt -n TTR -s 10000




Normalizing datasets
++++++++++++++++++++

TRANSIT also allows users to normalize datasets and export them afterwards. To normalize datasets, 3 positional arguments followed by optional flags.




Positional Arguments
--------------------


+----------------------+-------------------------------------------------------------------------+
| Positional Arguments | Definition                                                              |
+======================+=========================================================================+
| <wig-files>          | Comma-separated list of paths read-count datasets in .wig format        |
+----------------------+-------------------------------------------------------------------------+
| <annotation>         | Path to the annotation in .prot\_table or .GFF3 format.                 |
+----------------------+-------------------------------------------------------------------------+
| <output>             | Desired path and name of the output file                                |
+----------------------+-------------------------------------------------------------------------+



Optional Arguments
------------------

+----------------------+-----------------------------------------------------+
| Argument             | Definition                                          |
+======================+=====================================================+
| -n <String>          | Short name of the normalization method, e.g. -n TTR |
+----------------------+-----------------------------------------------------+



Example
~~~~~~~

.. code-block:: python

   python transit.py norm glycerol_H37Rv_rep1.wig,glycerol_H37Rv_rep2.wig H37Rv.prot_table glycerol_TTR.txt -n TTR







