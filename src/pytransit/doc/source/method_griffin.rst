.. _`griffin`:

=======
Griffin
=======

This is an earlier version of the Gumbel method that
identifies essential genes based on how unlikely 'gaps'
(or consecutive runs of TA sites with 0 insertions) are,
given the overall level of saturation.
It is a frequentist (non-Bayesian) model that uses
the Gumbel Extreme-Value Distribution as a likelihood function.
This is the analysis used in our paper on
`cholesterol catabolism (Griffin et al., 2011)
<http://www.ncbi.nlm.nih.gov/pubmed/21980284>`_.
All things considered, you are probably better off using the
hierarchical-Bayesian Gumbel model above, which does a better job
of estimating internal parameters.

|


.. rst-class:: transit_sectionend
----
