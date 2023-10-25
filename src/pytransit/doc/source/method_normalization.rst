
.. _normalization:

Normalization
=============


Proper normalization is important as it ensures that other sources of
variability are not mistakenly treated as real differences in
datasets. TRANSIT provides various normalization methods, which are
briefly described below:

- **TTR:**
    Trimmed Total Reads (TTR), normalized by the total
    read-counts (like totreads), but trims top and bottom 5% of
    read-counts. **This is the recommended normalization method for most cases**,
    as it has the benefit of compensating for differences in
    saturation (which is especially important for resampling).

- **nzmean:**
    Normalizes datasets to have the same mean over the
    non-zero sites.

- **totreads:**
    Normalizes datasets by total read-counts, and scales
    them to have the same mean over all counts.

- **zinfnb:**
    Fits a zero-inflated negative binomial model, and then
    divides read-counts by the mean. The zero-inflated negative
    binomial model will treat some empty sites as belonging to the
    "true" negative binomial distribution responsible for read-counts
    while treating the others as "essential" (and thus not influencing
    its parameters).

- **quantile:**
    Normalizes datasets using the quantile normalization
    method described by `Bolstad et al.
    (2003) <http://www.ncbi.nlm.nih.gov/pubmed/12538238>`_. In this
    normalization procedure, datasets are sorted, an empirical
    distribution is estimated as the mean across the sorted datasets
    at each site, and then the original (unsorted) datasets are
    assigned values from the empirical distribution based on their
    quantiles.  This actually doesn't work well on TnSeq data if a large fraction
    of TA sites have counts of 0 (ties).

- **betageom:**
    Normalizes the datasets to fit an "ideal" Geometric
    distribution with a variable probability parameter *p*. Specially
    useful for datasets that contain a large skew. 

- **nonorm:**
    No normalization is performed.



.. _BGC:

Beta-Geometric Correction (BGC)
~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a "bad" or poorly-behaving or "skewed" dataset (e.g. with mostly low
counts, dominated by a few high counts), right now the only remedy you
can try is applying the **Beta-Geometric correction (BGC)**, which is a
non-linear adjustment to the insertion counts in a wig file to make
them more like an ideal Geometric distribution 
(`DeJesus & Ioerger, 2016 <https://www.ncbi.nlm.nih.gov/pubmed/26932272>`_). 

Note, BGC is the only non-linear normalization available in Transit.  All the other
normalization methods, like TTR, are linear adjustments (scaling) of counts,
and so they can't correct for skewing.

BGC normalization helps 'de-skew' a dataset by bring the highest insertion
counts down into the range of a few thousands.  This sometimes improves
statistical analyses like resampling, etc by reducing noise.

In the GUI, when you are looking at wig files loaded, you can change
the normalization (e.g. from TTR to betageom) using the drop-down.  Be
aware that the Beta-Geometric normalization is *computationally-intensive* and
**might take few minutes per sample**.

If it looks like it might help (i.e. if the QQ-plot fits the diagonal better using BG
normalization),
you can create BG-corrected versions of individual wig files by
exporting them using the *normalize* command
on the command-line with '-n betageom' to specify normalization.



Command-line
------------

In addition to choosing normalization for various analyses in the GUI,
you can also call Transit to normalize wig files from the command-line,
as shown in this example:

Example
-------

::

  > python3 src/transit.py normalize --help

  usage: python3 src/transit.py normalize <input.wig> <output.wig> [-n TTR|betageom]
     or: python3 src/transit.py normalize -c <combined_wig> <output> [-n TTR|betageom]

  > python3 src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_TTR.wig -n TTR

  > python3 src/transit.py normalize Rv_1_H37RvRef.wig Rv_1_H37RvRef_BG.wig -n betageom

The normalize command now also works on :ref:`combined_wig <combined_wig>` files too.
If the input file is a combined_wig file, indicate it with a '-c' flag.



.. rst-class:: transit_sectionend
----
