import sys,random,numpy
import pytransit.tnseq_tools
from statsmodels.stats.multitest import fdrcorrection

def read_wig(fname):
  coords,counts = [],[]
  for line in open(fname):
    if line[0] not in "0123456789": continue
    w = line.rstrip().split()
    coord,cnt = int(w[0]),int(w[1])
    coords.append(coord)
    counts.append(cnt)
  return coords,counts

# remove all runs of zeros of length >= W

def remove_essential_regions(wig,W):
  runs = []
  i,n = 0,len(wig)
  while i<n:
    if wig[i]>0: i += 1
    else:
      j = i
      while j<n and wig[j]==0: j += 1
      runs.append((i,j))
      i = j
  counts = []
  for k,(i,j) in enumerate(runs):
    if j-i<W: counts += wig[i:j]
    if k<len(runs)-1:
      next = runs[k+1][0]
      counts += wig[j:next]
  return counts

def get_counts(coords,counts,gene):
  cnts = []
  start,end = gene['start'],gene['end']
  for i,co in enumerate(coords):
    if co>=start and co<=end: cnts.append(counts[i])
  return cnts

def sample_counts(counts,size,times):
  samples = []
  for i in range(times):
    samples.append(random.sample(counts,size)) # without replacement
  return samples

###################################

if len(sys.argv)<3:
  print "usage: python fitness_defect.py <comma_separated_list_of_wig_files> <prot_table>"
  sys.exit(0)

print "# command: python",
for x in sys.argv: print x,
print

coords,counts = read_wig(sys.argv[1])
genes = pytransit.tnseq_tools.read_genes(sys.argv[2])

noness = remove_essential_regions(counts,5) # run length
print '# sites: %s, noness: %s' % (len(counts),len(noness))
noness_arr = numpy.array((noness))
noness_NZvals = noness_arr[noness_arr>0]
print '# noness: zeros=%s, NZmean=%0.1f' % (noness_NZvals.size,numpy.mean(noness_NZvals))

cache = {}

results = []
for gene in genes:
  sys.stderr.write("%s %s\n" % (gene['rv'],gene['gene']))
  cnts = get_counts(coords,counts,gene)
  n = len(cnts)
  if n==0: continue
  nonzeros = [cnts[x] for x in numpy.nonzero(cnts)[0]]
  zeros,NZmean = n-len(nonzeros),numpy.mean(nonzeros) if len(nonzeros)>0 else 0
  sat = len(nonzeros)/float(n)
  tot = sum(cnts)
  mn = tot/float(n)

  # determine p-value by comparing to sum of counts for random draws of sites from noness
  N,alpha = 10000,0.05
  if n in cache: sample = cache[n]
  else: 
    sample = sample_counts(noness,n,N)
    cache[n] = sample
  samplesums = [sum(lst) for lst in sample]
  meansum = numpy.mean(samplesums)
  PC = 1 # pseudo-counts
  rel = (sum(cnts)+PC)/float(meansum+PC)
  LFC = numpy.log2(rel)

  lesser = len(list(filter(lambda x: x<=tot,samplesums)))
  greater = len(list(filter(lambda x: x>=tot,samplesums)))
  pval = min(lesser,greater)/float(N)

  vals = [gene[x] for x in "rv gene start end strand".split()]
  vals += [len(cnts),"%s" % tot,"%0.1f" % mn,"%0.3f" % sat,len(nonzeros),"%0.1f" % NZmean]
  vals += [int(meansum),"%0.3f" % rel,"%0.3f" % LFC,pval]

  results.append(vals)
  #if gene['rv']=='Rv0020c': break

pvals = [x[-1] for x in results]
qvals = list(fdrcorrection(pvals)[1])
results = [x+["%0.6f" % y] for x,y in zip(results,qvals)]

print '\t'.join("ORF gene start end strand TAs sum mean sat NZsites NZmean expec_sum FC LFC pval qval".split())
for vals in results:
  print '\t'.join([str(x) for x in vals])
