
.. rst-class:: transit_clionly

.. _zinb:

ZINB
====

The ZINB (Zero-Inflated Negative Binomial) method is used to determine
which genes exhibit *statistically significant variability across multiple
conditions*, in either the magnitude of insertion counts or local saturation, agnostically (in any one condition compared to the others).
Like :ref:`ANOVA <anova>`, the ZINB method takes a
*combined_wig* file (which combines multiple datasets in one file) and
a *samples_metadata* file (which describes which samples/replicates
belong to which experimental conditions).

ZINB can be applied to two or more conditions at a time.  Thus it
subsumes :ref:`resampling <resampling>`.  Our testing suggests that
ZINB typically identifies 10-20% more varying genes than resampling
(and vastly out-performs ANOVA for detecting significant variability
across conditions).  Furthermore, because of how ZINB treats magnitude
of read counts separately from local saturation in a gene, it
occasionally identifies genes with variability not detectable by
resampling analysis.

Note: ZINB analysis requires R (statistical analysis software)
to be installed on your system, along with the 'pscl' R package.
See :ref:`Installation Instructions <install-zinb>`.

|

How does it work?
-----------------

| For a formal description of how this method works, see our paper [ZINB]_:
|
|  Subramaniyam S, DeJesus MA, Zaveri A, Smith CM, Baker RE, Ehrt S, Schnappinger D, Sassetti CM, Ioerger TR. (2019).  `Statistical analysis of variability in TnSeq data across conditions using Zero-Inflated Negative Binomial regression. <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3156-z>`_, *BMC Bioinformatics*. 2019 Nov 21;20(1):603. doi: 10.1186/s12859-019-3156-z.



Example
-------

::

  > python3 transit.py zinb <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
        Optional Arguments:
        -n <string>                 :=  Normalization method. Default: -n TTR
        --prot_table <filename>     := for appending annotations of genes
        --condition                 :=  columnname (in samples_metadata) to use as the Condition. Default: "Condition"
        --covars <covar1,...>       :=  Comma separated list of covariates (in metadata file) to include, for the analysis.
        --interactions <covar1,...> :=  Comma separated list of covariates to include, that interact with the condition for the analysis.
        --exclude-conditions <cond1,...> :=  Comma separated list of conditions to ignore for the analysis. Default: None
        --include-conditions <cond1,...> :=  Comma separated list of conditions to include for the analysis. Default: All
        --ref <cond>                := which condition(s) to use as a reference for calculating LFCs (comma-separated if more than one) (by default, LFCs for each condition are computed relative to the grandmean across all condintions)
        -iN <float>                 :=  Ignore TAs occuring within given percentage of the N terminus. Default: -iN 5
        -iC <float>                 :=  Ignore TAs occuring within given percentage of the C terminus. Default: -iC 5
        -PC <N>                     :=  Pseudocounts used in calculating LFCs in output file. Default: -PC 5
        -winz                       := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
        -v                          := verbose, print out the model coefficients for each gene.
        --gene <Orf id or Gene name>:= Run method for one gene and print model output.


.. _combined_wig:

Combined wig files
------------------

Transit now supports a new file format called 'combined_wig' which basically
combines multiple wig files into one file (with multiple columns).  This is
used for some of the new analysis methods for larger collections of datasets, like :ref:`Anova <anova>`, :ref:`ZINB <zinb>`.
Combined_wig files can created through the Transit GUI
(File->Export->Selected_Datasets->Combined_wig), or via the command line.
You can specify the normalization method you want to use with a flag.
TTR is the default, but other relevant normalization options would be 'nonorm'
(i.e. preserve raw counts) and 'betageom' (this corrects for skew, but is slow).


::

  > python3 src/transit.py export combined_wig --help

  usage: python3 src/transit.py export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file>

  > python3 ../transit/src/transit.py export combined_wig Rv_1_H37RvRef.wig,Rv_2_H37RvRef.wig,Rv_3_H37RvRef.wig H37Rv.prot_table clinicals_combined_TTR.wig -n TTR



.. _samples_metadata:

Samples Metadata File
---------------------

Format of the *samples_metadata* file: a tab-separated file (which you
can edit in Excel) with 3 columns: Id, Condition, and Filename (it
must have these headers).  You can include other columns of info, but
do not include additional rows.  Individual rows can be commented out
by prefixing them with a '#'.  Here is an example of a samples
metadata file: The filenames should match what is shown in the header
of the combined_wig (including pathnames, if present).

Note: the Condition column should have a unique label for each distinct condition (the same label shared only among replicates).
If there are attributes that distinguish the conditions (such as strain, treatment, etc), they could be included as additional columns (e.g. covariates).

::

  ID      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig

Parameters
----------

The following parameters are available for the ZINB method:

-  **\-\-prot_table <filename>:** for appending annotations of genes
-  **\-\-include-conditions <cond1,...>:** Includes the given set of conditions from the ZINB test. Conditions not in this list are ignored. Note: this is useful for specifying the order in which the columns are listed in the output file.
-  **\-\-exclude-conditions <cond1,...>:** Ignores the given set of conditions from the ZINB test.
-  **\-\-ref <cond>:** which condition to use as a reference when computing LFCs in the output file. By default, LFCs for each condition are computed relative to the grandmean across all condintions.
-  **-n <TTR|betageom|nonorm...>:** Determines which Normalization method to 
   use when comparing datasets (Default: -n TTR). Proper normalization is important as it
   ensures that other sources of variability are not mistakenly treated
   as real differences. See the :ref:`Normalization <normalization>` section for a description
   of normalization method available in TRANSIT.
-  **\-\-condition <covar>:** column name (in samples_metadata) to use as the primary Condition being evaluated (to test for significant variability of insertions among conditions). Default: "Condition" (column name in metadata)
-  **\-\-covars <covar1,...>:** Comma separated list of covariates (columns in metadata file) to include, for the analysis.  If additional covariates distinguishing the samples are available, such as library, timepoint, or genotype, they may be "factored out" of the test of the primary condition. (variation due to covars is accounted for in the model, but not considered in evaluating the effect on variability due to the primary condition)
-  **\-\-interactions <covar1,...>:** Comma separated list of covariates (cols in metadata) to include, that interact with the condition for the analysis. (variation due to these variables *is* included in testing the effect of the main condition)
-  **-PC <N>:** Pseudocounts used in calculating LFCs in output file. (Default: -PC 5)
-  **-winz**: `winsorize <https://en.wikipedia.org/wiki/Winsorizing>`_ insertion counts for each gene in each condition. 
   Replace max count in each gene with 2nd highest.  This can help mitigate effect of outliers.

Covariates and Interactions
---------------------------

While ZINB is focus on identifying variability of insertion counts across conditions,
the linear model also allows you to take other variables into account.
There are two types of auxilliary variables: *covariates* and *interactions*. These can be provided as extra columns in the samples metadata file.
Covariates are attributes of the individual samples that could have a systematic
effect on the insertion counts which we would like to account for and subsequently ignore
(like nuissance variables). Examples include things like batch or library.

Interactions are extra variables for which we want to test their effect on the
main variable (or condition).  For example, suppose we collect TnSeq data at several
different timepoints (e.g. length of incubation or infection).  If we just test
time as the condition, we will be identifying genes that vary over time (if timepoints
are numeric, think of the model as fitting a 'slope' to the counts).
But suppose we have data for both a wild-type and knock-out strain.  Then we might be
interested in genes for which the time-dependent behavior *differs* between the two
strains (think: different 'slopes'). In such a case, we would say strain and time interact.


If covariates distinguishing the samples are available,
such as batch or library, they may be
incorporated in the ZINB model by using the **\-\-covars** flag and samples
metadata file. For example, consider the following samples metadata
file, with a column describing the batch information of each
replicate.

::

  ID      Condition    Filename                                     Batch
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig        B1
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig        B2
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig     B1
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig     B2
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig     B2

This information can be included to eliminate variability due to batch by using
the **\-\-covars** flag.

::

 python3 transit.py zinb combined.wig samples.metadata prot.table output.file --covars Batch


Similarly, an interaction variable may be included in the model.
This is specified by the user with the **\-\-interactions** flag,
followed by the name of a column in the samples metadata to test as the interaction
with the condition. If there are multiple interactions, they may be given as a comma-separated list.

To give an example,
consider an experiment where the condition represents
a treatment (e.g. with values 'treated' and 'control'), and we have another column
called Strain (with values 'wild-type' and 'mutant').
If we want to test whether the effect of the treatment (versus control)
differs depending on the strain, we could do this:

::

 python3 transit.py zinb combined.wig samples.metadata prot.table output.file --interactions Strain

In this case, the condition is implicitly assumed to be the column in the samples metadata file
labeled 'Condition'.  If you want to specify a different column to use as the primary condition to
test (for example, if Treatment were a distinct column), you can use the **\-\-condition** flag:

::

 python3 transit.py zinb combined.wig samples.metadata prot.table output.file --condition Treatment --interactions Strain



The difference between how covariates and interactions are handeled in the model
is discussed below in the section on Statistical Significance.

Categorical vs Numeric Covariates
---------------------------------

In some cases, covariates are intended to be treated as categorical
variables, like 'batch' or 'library' or 'medium'.
In other cases, a covariate might be a numeric value, such as
'time' or 'concentration', in which the ordering of values is
relevant.  The ZINB implementation tries to guess the type of each covariate.
If they are strings, they are treated as discrete factors (each with
their own distinct parameter).  If the given covariate can
be parsed as numbers, the model interprets them as real values.  In this
case, the covariate is treated as a linear factor (regressor), and is
incorporated in the model as a single coefficient, capturing the slope or
trend in the insertion counts as the covariate value increases.


Statistical Significance - What the P-values Mean in the ZINB Output
--------------------------------------------------------------------

Formally, the P-value is from a likelihood ratio test (LRT) between a
condition-dependent ZINB model (:math:`m_1`) and a
condition-independent (null) ZINB model (:math:`m_0`).

.. math::

  2 \ ln \frac{L(m_1)}{L(m_0)} \sim \chi^2_{df}

where L(.) is the ZINB likelihood function, and :math:`\chi^2_{df}` is
the chi-squared distribution with degrees of freedom (df) equal to
difference in the number of parameters bewteen the two models.  The p-value is
calculated based on this distribution.

In a simple case where variability across a set of conditions X is being tested,
you can think of the model approximately as:


.. math::

  m_1: ln \ \mu = \alpha_0+\vec\alpha X

where :math:`\mu` is an estimate of the mean (non-zero) insertion
count in a gene (a parameter in the likelihood function for ZINB),
:math:`\alpha_0` is a constant (the mean across all
conditions), and :math:`\vec\alpha` is a vector of coefficients
representing the deviation of the mean count in each condition.
(There is a corresponding equation for estimating the saturation as a
function of condition.)

To evaluate whether the variability across conditions is significant, we
compare to a null model, where the counts are estimated by the global mean only
(dropping the condition variable X).

.. math::

  m_0: ln \ \mu = \alpha_0

When a covariate C is available, it is incorporated in both models (additively),
to account for the effect of the covariate in :math:`m_1`. Coefficients in :math:`\vec\beta`
represent systematic effects on the mean count due to the covariate, and effectively
get subtracted out of the condition coefficients, but :math:`\vec\beta` is also
included in the null model :math:`m_0`, since we want to discount the effect of C on the
likelihood and focus on evaluting the effect of X.


.. math::

  m_1: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta C

  m_0: ln \ \mu = \alpha_0 + \vec\beta C


When an interaction I is being tested, it is incorporated *multiplicatively* in
the main model :math:`m_1` and *additively* in the null model :math:`m_0`:

.. math::

  m_1: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta I + \vec\gamma X*I

  m_0: ln \ \mu = \alpha_0 + \vec\alpha X + \vec\beta I

The meaning of this is that the coefficients :math:`\vec\alpha` and
:math:`\vec\beta` capture the additive effects of how the mean
insertion count in a gene depends on the condition variable and the
interaction variable, respectively, and the X*I term captures
additional (non-additive) deviations (which is the traditional way
interactions are handled in generalized linear models, GLMs).  Thus,
if there were no interaction, one would expect the mean in datasets
representing the *combination* of X and I to be predicted by the
offsets for each independently.  To the extend that this is not the
case, we say that X and I interaction, and the coefficients
:math:`\gamma` for X*I capture these deviations (non-additive
effects).

For example, think of condition X as Strain (e.g. wild-type vs mutant),
and interaction I as Treatment (e.g. treated versus control).
Then the main model would look like this:

.. math::

  m_1: ln \ \mu = \alpha_0 + \alpha_1 WT  + \alpha_2 mutant + \beta_1 control + \beta_2 treated + \gamma mutant * treated

and this would be compared to the following null model (without the interaction term):

.. math::

  m_0: ln \ \mu = \alpha_0 + \alpha_1 WT  + \alpha_2 mutant + \beta_1 control + \beta_2 treated




Output and Diagnostics
----------------------

The ZINB method outputs a tab-delimited file with results for each
gene in the genome. P-values are adjusted for multiple comparisons using
the Benjamini-Hochberg procedure (called "q-values" or "p-adj."). A
typical threshold for conditional essentiality on is q-value < 0.05.

+---------------------+-----------------------------------------------------------------+
| Column Header       | Column Definition                                               |
+=====================+=================================================================+
| Orf                 | Gene ID.                                                        |
+---------------------+-----------------------------------------------------------------+
| Name                | Name of the gene.                                               |
+---------------------+-----------------------------------------------------------------+
| TAs                 | Number of TA sites in Gene                                      |
+---------------------+-----------------------------------------------------------------+
| Means...            | Mean read-counts for each condition                             |
+---------------------+-----------------------------------------------------------------+
| LFCs...             | Log-fold-change (base 2) of mean insertion count relative to    |
|                     | mean across all conditions. Pseudo-counts of 5 are added.       |
|                     | If only 2 conditions, LFC is based on ratio of second to first. |
+---------------------+-----------------------------------------------------------------+
| NZmeans...          | Mean read-counts at non-zero zites for each condition           |
+---------------------+-----------------------------------------------------------------+
| NZpercs...          | Saturation (percentage of non-zero sites) for each condition    |
+---------------------+-----------------------------------------------------------------+
| p-value             | P-value calculated by the ZINB test.                            |
+---------------------+-----------------------------------------------------------------+
| p-adj               | Adjusted p-value controlling for the FDR (Benjamini-Hochberg)   |
+---------------------+-----------------------------------------------------------------+
| status              | Diagnostic information (explanation for genes not analyzed)     |
+---------------------+-----------------------------------------------------------------+


**LFCs** (log-fold-changes):
For each condition, the LFC is calculated as the log-base-2 of the
ratio of mean insertion count in that condition **relative to the
mean of means across all the conditions** (by default).
However, you can change this by desginating a specific reference condition using the flag **\-\-ref**.
(If there are multiple reference conditions, they may be given as a comma separated list.)
(If you are using interactions, it is more complicated to specify a reference condition by name because they have to include the interactions, e.g. as shown in the column headers in the output file.)
Pseudocount are incorporated to reduce the impact of noise on LFCs, based on the formula below.
The pseudocounts can be adjusted using the -PC flag.
Changing the pseudocounts (via -PC) can reduce the artifactual appearance of genes with
high-magnitude LFCs but that have small overall counts (which are susceptible to noise).
Changing the pseudocounts will not affect the analysis of statistical significance and hence number of varying genes, however.

::

  LFC = log2((mean_insertions_in_condition + PC)/(mean_of_means_across_all_conditions + PC))

|

Run-time
--------

A typical run of the ZINB method takes ~5 minutes to analze a combined wig
file with 6 conditions, 3 replicates per condition. It will, of
course, run more slowly if you have many more conditions.

|


.. rst-class:: transit_sectionend
----
