


Tutorial: Essentiality Analysis of the Entire Genome
====================================================

To illustrate how TRANSIT works, we are going to go through a tutorial
where we analyze datasets of H37Rv M. tuberculosis grown on glycerol
and cholesterol.

|

Run TRANSIT
-----------
Navigate to the directory containing the TRANSIT files, and run
TRANSIT:

::

    
    python PATH/src/transit.py




Adding the annotation file
--------------------------
Before we can analyze datasets, we need to add an annotation file for
the organism corresponding to the desired datasets. Click on the file
dialog button, on the top of the TRANSIT window (see image below), and
browse and select the appropriate annotation file. Note: Annotation
files must be in ".prot_table" format, described above.


.. image:: _images/transit_tutorial_annotation.png
   :width: 600
   :align: center


Adding the control datasets
---------------------------
We want to analyze datasets grown in glycerol to those grown in
cholesterol. We are choosing the datasets grown in glycerol as the
"Control" datasets. To add these, we click on the control sample file
dialog (see image below), and select the desired datasets (one by
one). In this example, we have two replicates:


.. image:: _images/transit_tutorial_control_datasets.png
   :width: 600
   :align: center



As we add the datasets they will appear in the table in the Control
samples section. This table will provide the following statistics
about the datasets that have been loaded so far: Total Number of
Reads, Density, Mean Read Count and Maximum Count. These statistics
can be used as general diagnostics of the datasets.


Visualizing read counts
-----------------------
TRANSIT allows us to visualize the read-counts of the datasets we have
already loaded. To do this, we must select the desired datasets
("Control+Click") and then click on "View -> Track View" in the menu
bar at the top of the TRANSIT window. Only those selected datasets
will be displayed:


.. image:: _images/transit_dataset_track_view.png
   :width: 600
   :align: center



This will open a window that allows that shows a visual representation
of the read counts at the TA sites throughout the genome. The scale of
the read counts can be set by changing the value of the "Max Read"
textbox on the right. We can browse around the genome by clicking on
the left and right arrowm, or search for a specific gene with the
search text box.


This window also allows us to save a .png image of the canvas for
future reference if desired (i.e. Save Img button).


Scatter plot
------------
We can also view a scatter plot of read counts of two selected
datasets. To achieve this we select two datasets (using "Control +
Clicck") and then clicking on "View -> Scatter Plot" in the menu bar
at the top of the TRANSIT window.


.. image:: _images/transit_dataset_scatter_graph.png
   :width: 600
   :align: center


A new window will pop-up, show a scatter plot of both of the selected
datasets. This window contains controls to zoom in and out (magnifying
glass), allowing us to focus in on a specific area. This is
particularly useful when large outliers may throw off the scale of the
scatter plot.



Essentiality analysis with the HMM method
-----------------------------------------
An alternative method for determining essentiality is the HMM method.
This method differs from the Gumbel method in that is capable of
assessing the essentiality of the entire genome, and is not limited to
a gene-level analysis (See above for discussions of the pros and cons
of each method). To run the HMM method we select it from the list of
(Himar1) methods on the Analysis at the top. This automatically displays 
the available options for the HMM methid. Because the HMM method estimates
parameters by examining the datasets, there is no need to set parameters
for the model. One important option provided is how to deal with replicate
datasets. Because the glycerol replicates had a mean read-count
between 53-85, we decide to sum read-counts together by selecting
"Sum" from the drop-down option.


.. image:: _images/transit_hmm_options.png
   :width: 600
   :align: center


Finally we click on the "Run HMM" button, and wait for the method to
finish. Once the analysis finishes, two new files will be created and
automatically added to the list of files in the Results Files section.
One file contains the output of states for each TA site in the genome.
The other file contains the analysis for each gene. We can display
each of the files be selecting them (individually) and clicking on the
"Display Table" button (one at a time).
Like for the Gumbel method, a break down of the states is provided at
the top of the table. In the case of glycerol, the HMM analysis
classifies 16.3% of the genome as belonging to the "Essential" state,
5.4% belonging to the Growth-Defect state, 77.1% to the Non-Essential
state, and 1.2% to the Growth Advantage state. This break down can be
used as a diagnostic, to see if the results match our expectations.
For example, in datasets with very low read-counts, or very low
density, the percentage of Growth-Defect states may be higher (e.g. >
10% ), which could indicate a problem.


The HMM sites file contains the state assignments for all the TA sites
in the genome. This file is particularly useful to browse for browsing
the different types of regions in the genome. We can use this file to
see how regions have different impacts on the growth-advantage (or
disadvantage) of the organism. For example, the PDIM locus, which is
required for virulance in vivo, results in a Growth-Advantage for the
organism when disrupted. We can see this in the HMM Sites file by
scrolling down to this region (Rv2930-Rv2939) and noticing the large
read-counts at these sites, and the how they are labeled "GA".






