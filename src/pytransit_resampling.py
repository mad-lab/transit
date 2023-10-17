import sys,numpy
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.stat_tools as stat_tools
import pytransit.norm_tools as norm_tools
from statsmodels.stats.multitest import fdrcorrection

# this is a stand-alone script that uses the pytransit library to do resampling between 2 groups of wig files
# analogous to running "python3 transit.py resampling..."; see also src/pytransit/analysis/resampling.py
# put transit/src/ in your PYTHON_PATH (e.g. using export in bash, or setenv in csh)
# prints output to stdout

if len(sys.argv)!=4:
  print("usage: python3 pytransit_resampling.py <comma_separated_wigs_condA> <comma_separated_wigs_condB> <prot_table>")
  sys.exit(0)

wigsA = sys.argv[1].split(',')
wigsB = sys.argv[2].split(',')
nA = len(wigsA)
prot_table = sys.argv[3]

(counts,sites) = transit_tools.get_validated_data(wigsA+wigsB) # data is a DxN numpy array, D=num datasets, N=num TA sites
print("counts.shape =",counts.shape)

print("normalizing data with TTR")
(data, factors) = norm_tools.normalize_data(counts, "TTR")

# this divides the raw data (counts at all TA sites) into a smaller array of counts for each gene based on its coordinates
genes = tnseq_tools.Genes(None,prot_table,data=data, position=sites)

results,pvals = [],[]
for gene in genes:
  if gene.n==0: continue # skip genes with 0 TA sites
  data1 = gene.reads[:nA] # counts at TA sites in gene for wigs for condition A
  data2 = gene.reads[nA:]
  stats = stat_tools.resampling(data1, data2, adaptive=True) # use defaults for other params
  (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) = stats # unpack
  vals = [gene.orf,gene.name,gene.n,round(mean1,1),round(mean2,1),round(log2FC,3),pval_2tail]
  results.append(vals)

# do FDR correction
pvals = [x[-1] for x in results]
qvals = fdrcorrection(pvals)[1]

print('\t'.join("orf gene numTA meanA meanB LFC Pval Qval".split()))
for vals,qval in zip(results,qvals):
  print('\t'.join([str(x) for x in vals+[qval]]))
print("%s/%s hits (qval<0.05)" % (sum([q<0.05 for q in qvals]),len(results)))
