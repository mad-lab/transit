Example Scripts (for Developers)
===============

|

Here are two example python scripts for developers showing how to use the *pytransit* package.

Example 1: src/pytransit_resampling.py
~~~~~~~~~~~~~~~~~~~~~

This is a stand-alone script that uses the pytransit library to do
resampling between two groups of wig files, analogous to running
"python3 transit.py resampling...".  It is a simplified version of
src/pytransit/analysis/resampling.py, which has more options.  

To run this, cd into your transit/src/ directory.  It prints output to
stdout. **Be sure to put transit/src/ in your PYTHON_PATH** 
(e.g. using export in bash, or setenv in csh).

::

  usage: python3 pytransit_resampling.py <comma_separated_wigs_condA> <comma_separated_wigs_condB> <prot_table>

  > cd transit/src/
  > python3 pytransit_resampling.py pytransit/data/glycerol_H37Rv_rep1.wig,pytransit/data/glycerol_H37Rv_rep2.wig pytransit/data/cholesterol_H37Rv_rep1.wig,pytransit/data/cholesterol_H37Rv_rep2.wig,pytransit/data/cholesterol_H37Rv_rep3.wig pytransit/genomes/H37Rv.prot_table 

The important functions illustrated in this script are:

 * **transit_tools.get_validated_data()**
 * **norm_tools.normalize_data()**
 * **stat_tools.resampling()**

.. code-block:: python

  # transit/src/pytransit_resampling.py
  import sys,numpy
  import pytransit.transit_tools as transit_tools
  import pytransit.tnseq_tools as tnseq_tools
  import pytransit.stat_tools as stat_tools
  import pytransit.norm_tools as norm_tools
  from statsmodels.stats.multitest import fdrcorrection
  
  if len(sys.argv)!=4:
    print("usage: python3 pytransit_resampling.py <comma_separated_wigs_condA> <comma_separated_wigs_condB> <prot_table>")
    sys.exit(0)
  
  wigsA = sys.argv[1].split(',')
  wigsB = sys.argv[2].split(',')
  nA = len(wigsA)
  prot_table = sys.argv[3]
  
  (counts,sites) = transit_tools.get_validated_data(wigsA+wigsB) # data is a DxN numpy array, D=num datasets, N=num TA sites
  print(counts.shape)

  print("normalizing data with TTR")
  (data, factors) = norm_tools.normalize_data(counts, "TTR")

  # this divides the raw data (counts at all TA sites) into a smaller array of counts for each gene based on its coordinates
 genes = tnseq_tools.Genes(None,prot_table,data=data, position=sites)

  results,pvals = [],[]
  for gene in genes:
    ii = gene.reads # coordinates of TA sites in gene
    if gene.n==0: continue # skip genes with 0 TA sites
    data1 = gene.reads[:nA] # counts at TA sites in gene for wigs for condition A
    data2 = gene.reads[nA:]
    stats = stat_tools.resampling(data1, data2, adaptive=True) # use defaults for other params
    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) = stats #     unpack
    vals = [gene.orf,gene.name,gene.n,round(mean1,1),round(mean2,1),round(log2FC,3),pval_2tail]
    results.append(vals)

  # do multiple-tests correction
  pvals = [x[-1] for x in results]
  qvals = fdrcorrection(pvals)[1]
  
  print('\t'.join("orf gene numTA meanA meanB LFC Pval Qval".split()))
  for vals,qval in zip(results,qvals):
    print('\t'.join([str(x) for x in vals+[qval]]))
  print("%s/%s hits (qval<0.05)" % (sum([q<0.05 for q in qvals]),len(results)))
  



Example 2: src/allpairs_resampling.py
~~~~~~~~~~~~~~~~~~~~~

While the example above shows how to read-in and process individual
wig files, this examples shows how to work with :ref:`combined_wig and
sample metadata files <combined_wig>`.  It does resampling between 
each pair of conditions, and prints out a matrix of hits
(statistically-signficant conditionally-essential genes).

To run this, cd into your transit/src/ directory.  It prints output to
stdout. **Be sure to put transit/src/ in your PYTHON_PATH** 
(e.g. using export in bash, or setenv in csh).

The important parts illustrated in this example script are: 

 * the **tnseq_tools.read_combined_wig()** function
 * how to read the metadata file
 * selecting counts from the data matrix for the samples associated with each condition

::

  import sys,numpy
  import pytransit.transit_tools as transit_tools
  import pytransit.tnseq_tools as tnseq_tools
  import pytransit.stat_tools as stat_tools
  import pytransit.norm_tools as norm_tools
  from statsmodels.stats.multitest import fdrcorrection
  
  # this is a stand-alone script that uses the pytransit library to do resampling between all pairs of conditions in a combined_wig file
  # metadata file indicates which replicates belong to which conditions
  # put transit/src/ in your PYTHON_PATH (e.g. using export in bash, or setenv in csh)
  # prints output to stdout
  
  if len(sys.argv)!=4:
    print("usage: python3 allpairs_resampling.py <combined_wig_file> <metadata_file> <prot_table>")
    sys.exit(0)
  
  combined_wig_file = sys.argv[1]
  metadata_file = sys.argv[2]
  prot_table = sys.argv[3]
  
  #################################
  
  print("reading data")
  (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(combined_wig_file)
  print("data.shape =",data.shape)
  
  print("normalizing using TTR")
  (data, factors) = norm_tools.normalize_data(data, "TTR")
  
  # there is a tnseq_tools.read_samples_metadata(), but it is kind of complicated
  # so just read through tab-separated file and collect sample Filenames associated with Conditions
  # metadata files have at least 3 columns: Id, Filename, Condition
  # the columns in the combined_wig are cross-referenced by Filename (from the original wigs)
  
  Conditions,Samples = [],{}
  CondCol,FnameCol = -1,-1
  for line in open(metadata_file):
    w = line.rstrip().split('\t') # tab-separated
    if CondCol==-1: CondCol,FnameCol = w.index("Condition"),w.index("Filename"); continue # error if headers not found
    cond,fname = w[CondCol],w[FnameCol]
    if cond not in Conditions: Conditions.append(cond)
    if cond not in Samples: Samples[cond] = []
    Samples[cond].append(fname)
  
  print("\nConditions\t: Samples")
  print("-------------------------")
  for i,cond in enumerate(Conditions):
    print("%s:%-8s" % (i+1,cond),"\t: ",", ".join(Samples[cond]))
   
  #################################
  
  genes = tnseq_tools.Genes(None,prot_table,data=data, position=sites) # divides counts at all TAs sites into groups by orf
  
  print()
  print("running resampling on each pair of conditions...")
  print("reporting number of conditionally essential genes (Qval<0.05)")
  for i,condA in enumerate(Conditions):
    for j,condB in enumerate(Conditions):
      if i<j:
        pvals = []
        for gene in genes:
          if gene.n==0: continue # skip genes with 0 TA sites
          idxA = [filenamesInCombWig.index(s) for s in Samples[condA]]
          idxB = [filenamesInCombWig.index(s) for s in Samples[condB]]
          data1 = gene.reads[idxA] # counts at TA sites in gene for wigs for condition A
          data2 = gene.reads[idxB]
          stats = stat_tools.resampling(data1, data2, adaptive=True) # use defaults for other params
          (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) = stats # unpack
          pvals.append(pval_2tail)
  
        qvals = fdrcorrection(pvals)[1]
        numhits = sum([q<0.05 for q in qvals])
        vals = ["%s:%-8s" % (i+1,condA),"%s:%-8s" % (j+1,condB),numhits]
        print("\t".join([str(x) for x in vals]))
  

Here is the output of this script for data from growth
of M. tuberculosis H37Rv on media containing iron supplied by
different vehicles (e.g. mycobactin, carboxymycobactin, hemin,
hemoglobin...), which requires genes in different pathways for uptake
`(Zhang et al., 2020) <https://pubmed.ncbi.nlm.nih.gov/32069330/>`_.
The raw data (wig files, with insertion counts at TA sites) have been
combined into a **combined_wig file** and a **metatdata file** that
describes which samples belong to which condition.  These files 
(*iron_combined_wig4.txt* and *iron_samples_metadata.txt*) can be found in
the transit data directory, transit/src/pytransit/data/.  You can
also compare this to the heatmap shown on the page for :ref:`corrplot
<corrplot>`.

::

  > cd transit/src/
  > python3 allpairs_resampling.py pytransit/data/iron_combined_wig4.txt pytransit/data/iron_samples_metadata.txt pytransit/genomes/H37Rv.prot_table
  reading data
  data.shape = (14, 74605)
  normalizing using TTR
  
  Conditions	: Samples
  -------------------------
  1:HighFeMBT 	:  HighFeMBT.wig, HighFeMBT2.wig, HighFeMBT3.wig
  2:LowFeMBT 	:  LowFeMBT.wig, LowFeMBT2.wig
  3:FeCMBT   	:  FeCMBT.wig, FeCMBT2b.wig
  4:Hemin    	:  Hemin.wig, Hemin2b.wig, Hemin3b.wig
  5:Hemoglobin 	:  Hemoglobin.wig, Hemoglobin2b.wig
  6:HeminMBT 	:  HeminMBT.wig, HeminMBT2.wig
  
  running resampling on each pair of conditions...
  reporting number of conditionally essential genes (Qval<0.05)
  1:HighFeMBT	2:LowFeMBT	38
  1:HighFeMBT	3:FeCMBT  	10
  1:HighFeMBT	4:Hemin   	136
  1:HighFeMBT	5:Hemoglobin	134
  1:HighFeMBT	6:HeminMBT	30
  2:LowFeMBT	3:FeCMBT  	35
  2:LowFeMBT	4:Hemin   	70
  2:LowFeMBT	5:Hemoglobin	71
  2:LowFeMBT	6:HeminMBT	5
  3:FeCMBT  	4:Hemin   	104
  3:FeCMBT  	5:Hemoglobin	106
  3:FeCMBT  	6:HeminMBT	44
  4:Hemin   	5:Hemoglobin	4
  4:Hemin   	6:HeminMBT	77
  5:Hemoglobin	6:HeminMBT	89


.. NOTE::

  Note, in allpairs_resampling.py, the FDR correction 
  to adjust P-values for testing all genes in parallel 
  is applied within each pairwise comparison.
  This correction is then repeated independently for each pair
  of conditions analyzed.  Formally, it would be more rigorous to apply
  the FDR correction one time at the end, to adjust P-values over all pairs and
  over all genes (which would be G*N*(N-1)/2 probabilities, where G in
  the number of genes in the genome, and N is the number of conditions).
  This would likely further reduce the number of significant
  conditionally-essential genes.  But for simplicity, this example
  script does not do that, because it is just designed to illustrate
  using the *pytransit* package.

