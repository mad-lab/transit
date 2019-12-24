import sys

try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0

if hasWx:
    import wx.xrc
    from wx.lib.buttons import GenBitmapTextButton
    from pubsub import pub
    import wx.adv

import os
import time
import math
import random
import numpy
import scipy.stats
from scipy.stats import norm
import datetime
import operator

from pytransit.analysis import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools

import io
from scipy.stats import hypergeom
import copy
from statsmodels.stats import multitest
# from datetime import datetime
import math

############# Description ##################

short_name = "pathway_enrichment"
long_name = "pathway_enrichment"
short_desc = "Pathway enrichment analysis"
long_desc = "Pathway enrichment analysis"
transposons = [] ##What's this for?
columns = ["[ID][descr]","Total genes","score","pval","padj","rank of genes"]

############# Analysis Method ##############

class PathwayAnalysis(base.TransitAnalysis):
  def __init__(self):    
    base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, PathwayMethod, PathwayGUI, [PathwayFile])


################## FILE ###################

class PathwayFile(base.TransitFile):

  def __init__(self):
    base.TransitFile.__init__(self, "#Example", columns)

  def getHeader(self, path):
    text = """This is file contains analysis of pathways enriched among sigificant genes from resampling."""
    return text


################# GUI ##################

# right now, there is no GUI interface for this analysis

class PathwayGUI(base.AnalysisGUI):

  def __init__(self):
    base.AnalysisGUI.__init__(self)

########## METHOD #######################

class PathwayMethod(base.AnalysisMethod):

  def __init__(self,resamplingFile,associationsFile,pathwaysFile,outputFile,method,PC=0,N=10000):
    base.AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, open(outputFile,"w"), None) # no annotation file
    self.resamplingFile = resamplingFile
    self.associationsFile = associationsFile
    self.pathwaysFile = pathwaysFile
    self.outputFile = outputFile 
    # self.output is the opened file, which will be set in base class
    self.method = method
    self.PC = PC # for FISHER
    self.N = N # for GSEA
    self.useZscores = False # for GSEA

  @classmethod
  def fromGUI(self, wxobj):
    pass

  @classmethod
  def fromargs(self, rawargs): 
    (args, kwargs) = transit_tools.cleanargs(rawargs)
    resamplingFile = args[0]    
    associations = args[1]
    pathways = args[2]
    output = args[3]
    method = kwargs.get("M", "FISHER")
    N = int(kwargs.get("N", "10000"))
    PC = int(kwargs.get("PC","2"))
    return self(resamplingFile,associations,pathways,output,method,PC=PC,N=N)

  @classmethod
  def usage_string(self):
    return """python %s pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FISHER|GSEA|GO>] """ % (sys.argv[0])

  def Run(self):
    self.transit_message("Starting Pathway Enrichment Method")
    start_time = time.time()

    # self.output in base class should be open by now
    self.print("# command: "+' '.join(sys.argv))
    self.print("# date: "+str(datetime.datetime.now()))

    if self.method=="FISHER": self.fisher_exact_test()
    elif self.method =="GSEA": self.GSEA3()
    else:
      method = "Not a valid method"
      self.progress_update("Not a valid method", 100)

  ############### GSEA ######################

  def makeindex(self,lst):
    index = {}
    for i in range(len(lst)): index[lst[i]] = i
    return index
    
  # based on GSEA paper (Subramanian et al, 2005, PNAS)
  # assume that Bindex is a dictionary that maps all genes into ranks
  # Bscores is a sorted list of values (e.g. correlations coeffs, LFCs, Z-scores)
    
  def enrichment_score(self,A,Bindex,Bscores,p=0):
    n2 = int(len(Bindex.keys())/2)
    ranks = [Bindex.get(x,n2) for x in A] # default to middle if not found
    ranks.sort()
    scores = [Bscores[i] for i in ranks]
    powers = [math.pow(abs(x),p) for x in scores]
    NR = sum(powers)
    Nmiss = len(Bscores)-len(ranks)
    best = -1
    powersum = 0
    for i in range(len(ranks)):
      powersum += powers[i]    
      Phit = powersum/float(NR)
      Pmiss = (ranks[i]-i)/float(Nmiss)
      es = abs(Phit-Pmiss) # looking for max deviation
      if es>best: best = es
    return best
    
  def mean_rank(self,A,Bindex): 
    n2 = len(Bindex.keys())/2
    return round(numpy.mean([Bindex.get(x,n2) for x in A]),1)

  # during initialization, self.resamplingFile etc have been set, and self.output has been opened    
    
  def GSEA3(self):
    N = self.N
    p = 1

    data,hits,headers = self.read_resampling_file(self.resamplingFile) # hits are not used in GSEA3
    headers = headers[-1].rstrip().split('\t')
    associations = self.read_associations(self.associationsFile)
    ontology = self.read_pathways(self.pathwaysFile)
    genenames = {}
    for gene in data: genenames[gene[0]] = gene[1]
    n2 = int(len(data)/2)
    terms = list(ontology.keys())
    terms2orfs = associations

    self.print("# method=GSEA, N=%d" % N)
    self.print("# total genes: %s, mean rank: %s" % (len(data),n2))

    pairs = [] # pair are: rv and LFC (or Zscore)
    if self.useZscores: 
      self.transit_message("using Z-scores to rank genes")
      if "Z-score" not in headers: self.transit_error("error: to rank genes by Z-score, you have to have used the -Z flag in resampling"); sys.exit(0)
      for w in data: pairs.append((w[0],float(w[headers.index("Z-score")]))) 
    else: 
      self.transit_message("using LFCs to rank genes")
      for w in data: pairs.append((w[0],float(w[headers.index("log2FC")]))) 
    pairs.sort(key=lambda x: x[1])
    sorted_rvs = [x[0] for x in pairs]
    sorted_scores = [x[1] for x in pairs]
    MainIndex = self.makeindex(sorted_rvs)

    Nperm = N
    permutations = []
    for i in range(Nperm): 
      lst = [x for x in sorted_rvs]
      random.shuffle(lst)
      permutations.append(lst)
    indexes = [self.makeindex(perm) for perm in permutations]
    
    results,Total = [],len(terms)

    for i,term in enumerate(terms):
      sys.stdout.flush()
      rvs = terms2orfs.get(term,[])
      if len(rvs)<=1: continue
      ranks = [MainIndex.get(x,n2) for x in rvs] # use n2 (avg rank) if gene not found
      es = self.enrichment_score(rvs,MainIndex,sorted_scores) # always positive, even if negative deviation, since I take abs
      mr = self.mean_rank(rvs,MainIndex)
      higher = 0
      for n,perm in enumerate(indexes):
        if self.enrichment_score(rvs,perm,sorted_scores)>es: higher += 1
        if n>100 and higher>10: break # adaptive
      pval = higher/float(n)
      vals = ['#',term,len(rvs),mr,es,pval,ontology.get(term,"?")]
      #sys.stderr.write(' '.join([str(x) for x in vals])+'\n')
      pctg=(100.0*i)/Total
      text = "Running Pathway Enrichment Method... %5.1f%%" % (pctg)
      self.progress_update(text, i)      
      results.append((term,mr,es,pval))
    
    results.sort(key=lambda x: x[-1])
    
    pvals = [x[-1] for x in results]
    rej,qvals = multitest.fdrcorrection(pvals)
    results = [tuple(list(res)+[q]) for res,q in zip(results,qvals)]

    self.output.write('\t'.join("#pathway description num_genes mean_rank GSEA_score pval qval genes".split())+'\n')
    for term,mr,es,pval,qval in results:
      rvs = terms2orfs[term]
      rvinfo = [(x,genenames.get(x,"?"),MainIndex.get(x,n2)) for x in rvs]
      rvinfo.sort(key=lambda x: x[2])
      rvs = ["%s/%s (%s)" % x for x in rvinfo]
      rvs = ' '.join(rvs)
      vals = [term,ontology.get(term,"?"),len(terms2orfs[term]),"%0.1f" % mr]+["%0.6f" % x for x in [es,pval,qval]]+[rvs]
      self.output.write('\t'.join([str(x) for x in vals])+'\n')
    self.output.close()

  #################################

  def read_resampling_file(self,filename):
    genes,hits,headers = [],[],[]
    for line in open(filename):
      if line[0]=='#': headers.append(line); continue
      w = line.rstrip().split('\t')
      genes.append(w)
      qval = float(w[-1])
      if qval<0.05: hits.append(w[0])
    return genes,hits,headers

  # assume these are listed as pairs (tab-sep)
  # return bidirectional hash (genes->[terms], terms->[genes]; each can be one-to-many, hence lists)
  def read_associations(self,filename):
    associations = {}
    for line in open(filename):
      if line[0]=='#': continue
      w = line.rstrip().split('\t')
      # store mappings in both directions
      for (a,b) in [(w[0],w[1]),(w[1],w[0])]:
        if a not in associations: associations[a] = []
        associations[a].append(b)
    return associations

  def read_pathways(self,filename):
    pathways = {}
    for line in open(filename):
      if line[0]=='#': continue
      w = line.rstrip().split('\t')
      pathways[w[0]] = w[1]    
    return pathways

  # HYPERGEOMETRIC 
  # scipy.stats.hypergeom.sf() is survival function (1-cdf), so only enriched genes will be significant
  # M = all genes
  # n = category members overall
  # N = sample size (resampling hits)
  # k = number of hits in category (intersection)

  def hypergeometric(self,k,M,n,N):
    return hypergeom.sf(k,M,n,N)

  def print(self,msg): self.output.write(msg+"\n")

  def fisher_exact_test(self):
    genes,hits,headers = self.read_resampling_file(self.resamplingFile)
    associations = self.read_associations(self.associationsFile)
    pathways = self.read_pathways(self.pathwaysFile)

    # how many genes are there, and how many have associations?
    # how many genes in associations are not listed in resampling file?
    # do all associations have a definition in pathways?
    # how many pathways have >1 gene? (out of total?) what is max?

    genes_with_associations = 0
    for gene in genes: 
      orf = gene[0]
      if orf in associations: genes_with_associations += 1
    self.print("# method=FISHER, PC=%s" % self.PC)
    self.print("# genes with associations=%s out of %s total" % (genes_with_associations,len(genes)))
    self.print("# significant genes (qval<0.05): %s" % (len(hits)))

    terms = list(pathways.keys())
    terms.sort()
    term_counts = [len(associations.get(term,[])) for term in terms]
    goodterms = []
    for term,cnt in zip(terms,term_counts):
      if cnt>1: goodterms.append(term)
    self.print("# %s out of %s pathways have >=1 gene; max has %s" % (len(goodterms),len(terms),term_counts[term_counts.index(max(term_counts))]))

    results = []
    for term in goodterms:
      n = len(associations[term]) # number of pathway members overall
      M = len(genes) # total genes
      N = len(hits) # number of resampling hits
      intersection = list(filter(lambda x: x in associations[term],hits))
      k = len(intersection)
      # add pseudo-counts
      PC = self.PC
      k_PC = int(k+PC)
      n_PC = n+int(M*PC/float(N)) # add same proportion to overall, round it
      expected = round((N*n/float(M)),2)
      enrichment = round((k+PC)/(expected+PC),3)
      pval = self.hypergeometric(k_PC,M,n_PC,N)
      results.append([term,M,n,N,k,expected,k_PC,n_PC,enrichment,pval])

    pvals = [x[-1] for x in results]
    rej,qvals = multitest.fdrcorrection(pvals)
    results = [x+[y] for x,y in zip(results,qvals)]

    genenames = {}
    for gene in genes: genenames[gene[0]] = gene[1]

    header = "#pathway total_genes(M) genes_in_path(n) significant_genes(N) signif_genes_in_path(k) expected k+PC n_adj_by_PC enrichement pval qval description genes"
    self.print('\t'.join(header.split()))

    results.sort(key=lambda x: x[-2]) # pvals
    for res in results:
      vals = res
      term = res[0]
      vals.append(pathways[term])
      intersection = list(filter(lambda x: x in associations[term],hits))      
      intersection = ["%s/%s" % (x,genenames[x]) for x in intersection]
      vals.append(' '.join(intersection))
      self.print('\t'.join([str(x) for x in vals]))

    self.transit_message("Adding File: %s" % (self.outputFile))
    self.add_file(filetype="Pathway Enrichment")
    self.finish()
    self.transit_message("Finished Pathway Enrichment Method") 

if __name__ == "__main__":

  app = PathwayMethod.fromargs(sys.argv[1:])
  app.Run()



