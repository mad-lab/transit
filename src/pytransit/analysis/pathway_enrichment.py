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

  def __init__(self,resamplingFile,associationsFile,pathwaysFile,outputFile,method,PC=0,N=10000,p=1):
    base.AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, open(outputFile,"w"), None) # no annotation file
    self.resamplingFile = resamplingFile
    self.associationsFile = associationsFile
    self.pathwaysFile = pathwaysFile
    self.outputFile = outputFile 
    # self.output is the opened file, which will be set in base class
    self.method = method
    self.PC = PC # for FISHER
    self.N = N # for GSEA
    self.p = p # for GSEA

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
    N = int(kwargs.get("N", "10000")) # for GSEA?
    p = int(kwargs.get("p","1")) # for GSEA?
    PC = int(kwargs.get("PC","2"))
    return self(resamplingFile,associations,pathways,output,method,PC=PC,N=N,p=p)

  @classmethod
  def usage_string(self):
    return """python3 %s pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FISHER|GSEA|GO>] [-PC <int>]""" % (sys.argv[0])

  def Run(self):
    self.transit_message("Starting Pathway Enrichment Method")
    start_time = time.time()

    # self.output in base class should be open by now
    self.write("# command: "+' '.join(sys.argv))
    self.write("# date: "+str(datetime.datetime.now()))

    if self.method=="FISHER": self.fisher_exact_test()
    elif self.method =="GSEA": self.GSEA()
    else:
      method = "Not a valid method"
      self.progress_update("Not a valid method", 100)

  ############### GSEA ######################

  def makeindex(self,lst):
    index = {}
    for i in range(len(lst)): index[lst[i]] = i
    return index
    
  # based on GSEA paper (Subramanian et al, 2005, PNAS)
  # xxx assume that Bindex is a dictionary that maps all genes into ranks, and Bscores maps all genes to SLPV
  # ranks and scores are hashes from genes into ranks and SLPV
  # when p=0, ES(S) reduces to the standard K-S statistic; p=1 is used in PNAS paper
    
  def enrichment_score(self,A,ranks,scores,p=0):
    n = len(ranks); n2 = int(n/2.0)
    Aranks = [ranks.get(x,n2) for x in A] # default to middle if not found
    Ascores = [scores.get(x,0) for x in A] # default to 0 if not found
    pairs = list(zip(Aranks,Ascores))
    pairs.sort() # sort A by ranks
    Aranks,Ascores = [x[0] for x in pairs],[x[1] for x in pairs]
    powers = [math.pow(abs(x),p) for x in Ascores]
    NR = sum(powers)
    if NR==0: return 0 # special case
    Nmiss = n-len(A) # totalGenes-hits
    powersum,best = 0,-1
    for i in range(len(powers)):
      powersum += powers[i]    
      Phit = powersum/float(NR)
      Pmiss = (Aranks[i]-i)/float(Nmiss)
      es = abs(Phit-Pmiss) # looking for max deviation
      if es>best: best = es
    return best
    
  def mean_rank(self,A,orfs2ranks): 
    n2 = len(orfs2ranks.keys())/2
    return round(numpy.mean([orfs2ranks.get(x,n2) for x in A]),1)

  # during initialization, self.resamplingFile etc have been set, and self.output has been opened    
    
  def GSEA(self):
    data,hits,headers = self.read_resampling_file(self.resamplingFile) # hits are not used in GSEA()
    headers = headers[-1].rstrip().split('\t')
    associations = self.read_associations(self.associationsFile)
    ontology = self.read_pathways(self.pathwaysFile)
    genenames = {}
    for gene in data: genenames[gene[0]] = gene[1]
    n2 = int(len(data)/2)
    terms = list(ontology.keys())
    terms2orfs = associations
    allgenes = [x[0] for x in data]

    self.write("# method=GSEA, using SLPV to rank genes, Nperm=%d" % self.N)
    self.write("# total genes: %s, mean rank: %s" % (len(data),n2))

    # rank by SLPV=sign(LFC)*log10(pval)
    # note: genes with lowest p-val AND negative LFC have highest scores (like positive correlation)
    # there could be lots of ties with pval=0 or 1, but that's OK
    LFC_col = headers.index("log2FC")
    Pval_col = headers.index("p-value")
    pairs = [] # pair are: rv and score (SLPV)
    for w in data:
      orf,LFC,Pval = w[0],float(w[LFC_col]),float(w[Pval_col])
      SLPV = (-1 if LFC<0 else 1)*math.log(Pval+0.000001,10)
      pairs.append((orf,SLPV))
    pairs.sort(key=lambda x: x[1],reverse=True) # emulate ranking genes with *higher* correlation at top
    orfs2rank,orfs2score = {},{}
    for i,(orf,score) in enumerate(pairs): 
      orfs2score[orf] = score
      orfs2rank[orf] = i

    Nperm = self.N
    results,Total = [],len(terms)
    for i,term in enumerate(terms):
      sys.stdout.flush()
      orfs = terms2orfs.get(term,[])
      if len(orfs)<=1: continue
      mr = self.mean_rank(orfs,orfs2rank)
      es = self.enrichment_score(orfs,orfs2rank,orfs2score,p=self.p) # always positive, even if negative deviation, since I take abs
      higher = 0
      for n in range(Nperm):
        perm = random.sample(allgenes,len(orfs)) # compare to ES for random sets of genes of same size
        if self.enrichment_score(perm,orfs2rank,orfs2score,p=self.p)>es: higher += 1
        if n>100 and higher>10: break # adaptive
      pval = higher/float(n)
      vals = ['#',term,len(orfs),mr,es,pval,ontology.get(term,"?")]
      #sys.stderr.write(' '.join([str(x) for x in vals])+'\n')
      pctg=(100.0*i)/Total
      text = "Running Pathway Enrichment Method... %5.1f%%" % (pctg)
      self.progress_update(text, i)      
      results.append((term,mr,es,pval))
    
    results.sort(key=lambda x: x[1])
    pvals = [x[-1] for x in results]
    rej,qvals = multitest.fdrcorrection(pvals)
    results = [tuple(list(res)+[q]) for res,q in zip(results,qvals)]

    n2 = int(len(data)/2)
    up,down = 0,0
    for term,mr,es,pval,qval in results:
      if qval<0.05:
        if mr<n2: up += 1
        else: down += 1

    self.write("# significant pathways enriched for conditionally ESSENTIAL genes: %s (qval<0.05, mean_rank<%s) (includes genes that are MORE required in condition B than A)" % (up,n2))
    for term,mr,es,pval,qval in results:
      if qval<0.05 and mr<n2: self.write("#   %s %s (mean_rank=%s)" % (term,ontology.get(term,"?"),mr))
    self.write("# significant pathways enriched for conditionally NON-ESSENTIAL genes: %s (qval<0.05, mean_rank>%s) (includes genes that are LESS required in condition B than A)" % (down,n2))
    for term,mr,es,pval,qval in results:
      if qval<0.05 and mr>n2: self.write("#   %s %s (mean_rank=%s)" % (term,ontology.get(term,"?"),mr))
    self.write("# pathways sorted by mean_rank")

    self.output.write('\t'.join("#pathway description num_genes mean_rank GSEA_score pval qval genes".split())+'\n')
    for term,mr,es,pval,qval in results:
      rvs = terms2orfs[term]
      rvinfo = [(x,genenames.get(x,"?"),orfs2rank.get(x,n2)) for x in rvs]
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

  def write(self,msg): self.output.write(msg+"\n")

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
    self.write("# method=FISHER, PC=%s" % self.PC)
    self.write("# genes with associations=%s out of %s total" % (genes_with_associations,len(genes)))
    self.write("# significant genes (qval<0.05): %s" % (len(hits)))

    terms = list(pathways.keys())
    terms.sort()
    term_counts = [len(associations.get(term,[])) for term in terms]
    goodterms = []
    for term,cnt in zip(terms,term_counts):
      if cnt>1: goodterms.append(term)
    self.write("# %s out of %s pathways have >=1 gene; max has %s" % (len(goodterms),len(terms),term_counts[term_counts.index(max(term_counts))]))

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
    self.write('\t'.join(header.split()))

    results.sort(key=lambda x: x[-2]) # pvals
    for res in results:
      vals = res
      term = res[0]
      vals.append(pathways[term])
      intersection = list(filter(lambda x: x in associations[term],hits))      
      intersection = ["%s/%s" % (x,genenames[x]) for x in intersection]
      vals.append(' '.join(intersection))
      self.write('\t'.join([str(x) for x in vals]))

    self.transit_message("Adding File: %s" % (self.outputFile))
    self.add_file(filetype="Pathway Enrichment")
    self.finish()
    self.transit_message("Finished Pathway Enrichment Method") 

####################################################

if __name__ == "__main__":

  app = PathwayMethod.fromargs(sys.argv[1:])
  app.Run()



