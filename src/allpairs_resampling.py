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

