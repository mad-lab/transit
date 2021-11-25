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

  def __init__(self,resamplingFile,associationsFile,pathwaysFile,outputFile,method,PC=0,Nperm=10000,p=0,ranking="SLPV",Pval_col=-2,Qval_col=-1,LFC_col=6): # default cols are for resampling files
    base.AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, open(outputFile,"w"), None) # no annotation file
    self.resamplingFile = resamplingFile
    self.associationsFile = associationsFile
    self.pathwaysFile = pathwaysFile
    self.outputFile = outputFile 
    # self.output is the opened file, which will be set in base class
    self.method = method
    self.Pval_col = Pval_col
    self.Qval_col = Qval_col
    self.LFC_col = LFC_col
    self.PC = PC # for FET
    self.Nperm = Nperm # for GSEA
    self.p = p # for GSEA
    self.ranking = ranking # for GSEA

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
    method = kwargs.get("M", "FET")
    Pval_col = int(kwargs.get("Pval_col","-2")) # default cols are for resampling files
    Qval_col = int(kwargs.get("Qval_col","-1"))
    LFC_col = int(kwargs.get("LFC_col","6"))
    PC = int(kwargs.get("PC","2")) # for FET
    Nperm = int(kwargs.get("Nperm", "10000")) # for GSEA
    p = float(kwargs.get("p","0")) # for GSEA
    ranking = kwargs.get("ranking","SLPV") # for GSEA

    if method not in "FET GSEA ONT".split(): 
      print("error: method %s not recognized" % method)
      print(self.usage_string()); 
      sys.exit(0)

    return self(resamplingFile,associations,pathways,output,method,PC=PC,Nperm=Nperm,p=p,ranking=ranking,Pval_col=Pval_col,Qval_col=Qval_col,LFC_col=LFC_col)

  @classmethod
  def usage_string(self):
    return """python3 %s pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FET|GSEA|GO>] [-PC <int>] [-ranking SLPV|LFC] [-p <float>] [-Nperm <int>] [-Pval_col <int>] [-Qval_col <int>]  [-LFC_col <int>]

Optional parameters:
 -M FET|GSEA|ONT:     method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Analysis (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)
 -Pval_col <int>    : indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2)
 -Qval_col <int>    : indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)
 for GSEA...
 -ranking SLPV|LFC  : SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling 
 -LFC_col <int>     : indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)
 -p <float>         : exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)
 -Nperm <int>       : number of permutations to simulate for null distribution to determine p-value (default=10000)
 for FET...
 -PC <int>          :  pseudo-counts to use in calculating p-value based on hypergeometric distribution (default=2)
""" % (sys.argv[0])

  def Run(self):
    self.transit_message("Starting Pathway Enrichment Method")
    start_time = time.time()

    # self.output in base class should be open by now
    self.write("# command: python3 "+' '.join(sys.argv))
    now = str(datetime.datetime.now())
    now = now[:now.rfind('.')]
    self.write("# date: "+now)
    self.write("# pathway method: "+self.method)

    if self.method=="FET": self.fisher_exact_test()
    elif self.method =="GSEA": self.GSEA()
    elif self.method =="ONT": self.Ontologizer()
    else:
      method = "Not a valid method"
      self.progress_update("Not a valid method", 100)

  def write(self,msg): self.output.write(msg+"\n")

  def read_resampling_file(self,filename):
    genes,hits,headers = [],[],[]
    for line in open(filename):
      if line[0]=='#': headers.append(line); continue
      w = line.rstrip().split('\t')
      genes.append(w)
      qval = float(w[self.Qval_col])
      if qval<0.05: hits.append(w[0])
    return genes,hits,headers

  # assume these are listed as pairs (tab-sep)
  # return bidirectional hash (genes->[terms], terms->[genes]; each can be one-to-many, hence lists)
  # filter could be a subset of genes we want to focus on (throw out the rest)

  def read_associations(self,filename,filter=None):
    associations = {}
    for line in open(filename):
      if line[0]=='#': continue
      w = line.rstrip().split('\t')
      if filter!=None and w[0] not in filter: continue # skip genes in association file that are not relevant (i.e. not in resampling file)
      # store mappings in both directions
      for (a,b) in [(w[0],w[1]),(w[1],w[0])]:
        if a not in associations: associations[a] = []
        if b not in associations[a]: associations[a].append(b) # ignore duplicates
    return associations

  def read_pathways(self,filename):
    pathways = {}
    for line in open(filename):
      if line[0]=='#': continue
      w = line.rstrip().split('\t')
      pathways[w[0]] = w[1]    
    return pathways

  ############### GSEA ######################

  def makeindex(self,lst):
    index = {}
    for i in range(len(lst)): index[lst[i]] = i
    return index
    
  # based on GSEA paper (Subramanian et al, 2005, PNAS)
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
    orfs_in_resampling_file = [w[0] for w in data]
    headers = headers[-1].rstrip().split('\t') # last line prefixed by '#'
    associations = self.read_associations(self.associationsFile,filter=orfs_in_resampling_file) # bidirectional map; includes term->genelist and gene->termlist
    # filter: project associations (of orfs to pathways) onto only those orfs appearing in the resampling file

    ontology = self.read_pathways(self.pathwaysFile)
    genenames = {}
    for gene in data: genenames[gene[0]] = gene[1]
    n2 = int(len(data)/2)
    terms = list(ontology.keys())
    terms2orfs = associations
    allgenes = [x[0] for x in data]

    self.write("# method=GSEA, Nperm=%d, p=%d" % (self.Nperm,self.p))
    self.write("# ranking genes by %s" % self.ranking)
    self.write("# total genes: %s, mean rank: %s" % (len(data),n2))

    # rank by SLPV=sign(LFC)*log10(pval)
    # note: genes with lowest p-val AND negative LFC have highest scores (like positive correlation)
    # there could be lots of ties with pval=0 or 1, but so be it (there are probably fewer such ties than with Qvals) (randomize order below) 
    pairs = [] # pair are: rv and score (SLPV)
    for w in data:
      orf = w[0]
      if self.ranking=="SLPV": 
        #Pval_col = headers.index("p-value")
        Pval = float(w[self.Pval_col])
        LFC = float(w[self.LFC_col])
        SLPV = (-1 if LFC<0 else 1)*math.log(Pval+0.000001,10)
        pairs.append((orf,SLPV))
      elif self.ranking=="LFC": 
        #LFC_col = headers.index("log2FC")
        LFC = float(w[self.LFC_col])
        pairs.append((orf,LFC))

    # pre-randomize ORFs, to avoid genome-position bias in case of ties in pvals (e.g. 1.0)
    indexes = range(len(pairs))
    indexes = numpy.random.permutation(indexes).tolist() 
    pairs = [pairs[i] for i in indexes]

    pairs.sort(key=lambda x: x[1],reverse=True) # emulate ranking genes with *higher* correlation at top
    orfs2rank,orfs2score = {},{}
    for i,(orf,score) in enumerate(pairs): 
      orfs2score[orf] = score
      orfs2rank[orf] = i

    Nperm = self.Nperm
    results,Total = [],len(terms)
    for i,term in enumerate(terms):
      sys.stdout.flush()
      orfs = terms2orfs.get(term,[])
      num_genes_in_pathway = len(orfs)
      if num_genes_in_pathway<2: continue # skip pathways with less than 2 genes
      mr = self.mean_rank(orfs,orfs2rank)
      es = self.enrichment_score(orfs,orfs2rank,orfs2score,p=self.p) # always positive, even if negative deviation, since I take abs
      higher = 0
      for n in range(Nperm):
        perm = random.sample(allgenes,num_genes_in_pathway) # compare to enrichment score for random sets of genes of same size
        e2 = self.enrichment_score(perm,orfs2rank,orfs2score,p=self.p)
        if e2>es: higher += 1
        if n>100 and higher>10: break # adaptive: can stop after seeing 10 events (permutations with higher ES)
      pval = higher/float(n)
      vals = ['#',term,num_genes_in_pathway,mr,es,pval,ontology.get(term,"?")]
      #sys.stderr.write(' '.join([str(x) for x in vals])+'\n')
      pctg=(100.0*i)/Total
      text = "Running Pathway Enrichment Method... %5.1f%%" % (pctg)
      self.progress_update(text, i)      
      results.append((term,mr,es,pval))
    
    results.sort(key=lambda x: x[1]) # sort on mean rank
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

  ########## Fisher Exact Test ###############

  # HYPERGEOMETRIC 
  # scipy.stats.hypergeom.sf() is survival function (1-cdf), so only enriched genes will be significant
  # M = all genes
  # n = category members overall
  # N = sample size (resampling hits)
  # k = number of hits in category (intersection)

  def hypergeometric(self,k,M,n,N):
    return hypergeom.sf(k,M,n,N)

  def fisher_exact_test(self):
    genes,hits,headers = self.read_resampling_file(self.resamplingFile) # use self.Qval_col to determine hits
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
    self.write("# method=FET, PC=%s" % self.PC)
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

  ########## Ontologizer ###############

  # this method is restricted to GO terms

  # this implements the union method of:
  #  Grossman et al (2007). Improved Detection of overrepresentation of Gene-Ontology 
  #  annotation with parent-child analysis. Bioinformatics, 23(22):3024-3031.

  def Ontologizer(self):

    def warning(s): sys.stderr.write("%s\n" % s) # use self.warning? prepend method name?
  
    # returns True if ch is a descendant of GO (as GO terms, like "GO:0006810")
    
    def descendant_of(ch,GO):
      if ch==GO: return True
      for node in parents.get(ch,[]):
        if descendant_of(node,GO)==True: return True
      return False
    
    # visited is a hash
    # assume B is above A
    # do upward DFS following parent pointers till hit NULL
    # return all nodes on all paths? include A and B?
    
    def get_ancestors(A,visited):
      if A in visited: return visited
      visited[A] = 1
      for parent in parents.get(A,[]): get_ancestors(parent,visited)
      return visited
    
    def depth(go):
      if go not in parents: return 0
      return 1+min([depth(x) for x in parents[go]])
    
    ontology,parents = {},{}
  
    GOannot = self.associationsFile
    OBOfile = self.pathwaysFile

    for line in open(OBOfile):
      if line[:3]=="id:": id = line[4:-1]
      if line[:5]=="is_a:":
        parent = line.split()[1]
        if id not in parents: parents[id] = []
        parents[id].append(parent)
      if len(line)<2: id = None
      if line[:5]=="name:": ontology[id] = line[6:-1]
  
    rv2gos,go2rvs = {},{}
    MINTERMS,MAXTERMS = 2,300
  
    for line in open(GOannot):
      w = line.rstrip().split('\t')
      rv,go = w[0],w[1]
      if rv not in rv2gos: rv2gos[rv] = []
      if go not in rv2gos[rv]: rv2gos[rv].append(go)
      if go not in go2rvs: go2rvs[go] = []
      if rv not in go2rvs[go]: go2rvs[go].append(rv)
      # expand to all parents...
      BP,CC,MF = ["GO:0008150","GO:0005575","GO:0003674"]
      for g in get_ancestors(go,{}):
        if g not in rv2gos[rv]: 
          rv2gos[rv].append(g)
          if g not in go2rvs: go2rvs[g] = []
          go2rvs[g].append(rv)
  
    warning("GO terms with at least one ORF: %s" % len(go2rvs.keys())) # what about between MIN and MAX?
    for go in go2rvs.keys():
      if go not in ontology: warning("not found: %s" % go) # also indicate which gene?

    # could use class method, but would have to adapt it: 
    # genes,hits,headers = self.read_resampling_file(self.resamplingFile)
  
    genes,pvals = {},[]
    allorfs,studyset = [],[]
    for line in open(self.resamplingFile):
      if line[0]=="#": continue
      w = line.rstrip().split('\t')
      genes[w[0]] = w
      allorfs.append(w[0])
      #pval,qval = float(w[-2]),float(w[-1])
      pval,qval = float(w[self.Pval_col]),float(w[self.Qval_col])
      if qval<0.05: studyset.append(w[0])
      pvals.append((w[0],pval))
    pvals.sort(key=lambda x: x[1])
    ranks = {}
    for i,(rv,pval) in enumerate(pvals): ranks[rv] = i+1
    self.write("# number of resampling hits (qval<0.05): %s" % len(studyset))
  
    counts = []
    n,a = len(allorfs),len(studyset)
    for go in go2rvs.keys():
      if go not in ontology: continue
      m = len(go2rvs[go]) # orfs in the GO terms in the genome overall
      if m>=MINTERMS and m<=MAXTERMS:
        b = len(list(filter(lambda x: x in go2rvs[go],studyset))) # orfs with GO term in studyset
        if b==0: continue
        # calc p=len(rvs for parents of go), q=len(subset of parent rvs in studyset)
        P = set()
        if go not in parents: continue
        npar = len(parents[go])
        for par in parents[go]: P = P.union(go2rvs[par])
        p,q = len(P),len(P.intersection(studyset))
        enrich = ((b+1)/float(q+1))/((m+1)/float(p+1)) # enrichment
        # parents: p overall, q in studyset; GO in studyset: b out of q (subset of studyset labeled with a parent GO term)
        # assume m orfs have GO term inside parent (same as overall)         
        # subtract 1 from b possibly because we are using 1-cdf? this is needed to get same numeric results as 'phypergeometric' in Ontologizer
        pval = 1.0-scipy.stats.hypergeom.cdf(b-1,p,m,q) 
        counts.append([go,n,m,p,q,a,b,npar,enrich,pval])
    
    pvals = [x[-1] for x in counts]
    qvals = multitest.fdrcorrection(pvals)[1]
  
    counts = [x+[y] for x,y in zip(counts,qvals)]
    counts.sort(key=lambda x: x[-1])

    self.write("# number of GO terms that are significantly enriched (qval<0.05): %s" % len(list(filter(lambda x: x[-1]<0.05,counts))))
  
    # 1-sided pvals, only report enriched terms, not depleted
  
    self.write('\t'.join("GO_term description total_orfs orfs_in_GO orfs_in_parents hits_in_parents hits GO_in_hits num_parent_nodes enrichment pval qval genes_in_intersection_of_hits_and_GO_term".split())) # assume self.outputFile has already been opened
    for (go,n,m,p,q,a,b,npar,enrich,pval,qval) in counts:
      hits = filter(lambda x: x in go2rvs[go],studyset)
      hits = [(x,genes[x][1],ranks[x]) for x in hits]
      hits.sort(key=lambda x: x[2]) # sort on ranks
      hits = ["%s/%s(%s)" % (rv,gene,rank) for (rv,gene,rank) in hits]
      hits = ','.join(hits)    # add gene names
      vals = [go,ontology[go],n,m,a,b,p,q,npar,round(enrich,3),round(pval,6),round(qval,6),hits]
      self.write('\t'.join([str(x) for x in vals]))

    self.transit_message("Adding File: %s" % (self.outputFile))
    self.add_file(filetype="Pathway Enrichment")
    self.finish()
    self.transit_message("Finished Pathway Enrichment Method") 


####################################################

if __name__ == "__main__":

  app = PathwayMethod.fromargs(sys.argv[1:])
  app.Run()



