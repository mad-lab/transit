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
short_desc = "Gene set enrichment analysis"
long_desc = "Gene set enrichment analysis"
transposons = [] ##What's this for?
columns = ["[ID][descr]","Total genes","score","pval","padj","rank of genes"]

############# Analysis Method ##############

class GSEAAnalysis(base.TransitAnalysis):
  def __init__(self):    
    base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, GSEAMethod, GSEAGUI, [GSEAFile])


################## FILE ###################

class GSEAFile(base.TransitFile):

  def __init__(self):
    base.TransitFile.__init__(self, "#Example", columns)

  def getHeader(self, path):
    text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
    return text


################# GUI ##################

class GSEAGUI(base.AnalysisGUI):

  def __init__(self):
    base.AnalysisGUI.__init__(self)

########## METHOD #######################

class GSEAMethod(base.SingleConditionMethod):
  """   
  Example
 
  """
  def __init__(self,
        resamplingFile,        
        geneSetFile,
        output_file,p,N,M):

    base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, resamplingFile, geneSetFile, output_file, p,N,M)
    self.resamplingFile=resamplingFile
    self.geneSetFile = geneSetFile
    self.output = output_file
    self.p=p
    self.N=N
    self.M=M

####################################################################################
  ###############FILES################################
  def loadD(self,fileName):
    file = open(fileName,"r")
    dict=[]
    ORFNameDict={}    
    for line in file:      
      if not(line.startswith("#")):        
        line=line.strip().split("\t") #this will return three elements        
        ORFNameDict[line[0]]=line[1]
        dict.append([line[0],float(line[10])])  
    file.close()
    return dict,ORFNameDict

    # file = open(fileName,"r")
    # line = file.readline()
    # while line.startswith("#"):
    #   line = file.readline()    
    # dict=[]
    # line = line.split("\t")  
    # dict.append([line[0],float(line[10])])
    # ORFNameDict={line[0]:line[1]}
    # for line in file:    
    #   line = line.strip().split("\t")
    #   dict.append([line[0],float(line[10])])
    #   ORFNameDict[line[0]]=line[1]
    # return dict,ORFNameDict

  def getM(self,protTable):
    file = open(protTable,"r")
    return len(file.read().splitlines())

  def loadGoTermsGoTermAsKey(self,fileName):
    dict={}
    descr={}
    file = open(fileName,"r")
    for f in file:
      if not(f.startswith("#")):        
        line = f.strip().split("\t") # It will return the id , description and list of ORFS        
        if len(line)!=3:           
          self.output.write("Format Error in"+fileName+"\n")
        if not(line[0] in dict):
          dict[line[0]]=[]          
        dict[line[0]]+=line[2].split()
        descr[line[0]] = line[1]
    return dict,descr


  def saveInterestingPaths(self,fileName,m,n):
    inputF = open(fileName,"r")
    o = fileName.split(".")
    outputF = open(o[0]+"_"+str(m)+"_"+str(n)+"."+o[1],"w")
    for line in inputF:
      cad = line
      line = line.split(",")[1].split()
      lenLine = len(line)    
      if lenLine>m and lenLine<n:
        outputF.write(cad)
    outputF.close()
    

  ################PREPROCESSING###################
  #Keeping only the value with p-value < 0.05
  #I receives a Dictionary with the gene as key
  #DEPRECATED!! D only has two columns
  def filteredPValue(self,D):  
    return {d:D[d] for d in D if D[d][4]<0.05}

  ################STATISTICAL ANALYSIS############
  #M = the whole Genome
  #n = hits
  #N = sample size
  def hipergeometric(self,k,M,n,N):
    return hypergeom.sf(k,M,n,N)


  # HYPERGEOMETRIC

  def hyperGeometricTest(self,dict_k,M,dict_n,N):
    result={}
    for key in dict_n:
      n = len(dict_n[key])
      I = set(dict_k) & set(dict_n[key])
      k = len(I)
      result[key] = {"parameters":str(n)+"\t"+str(k),"p-value": hypergeom.sf(k,M,n,N), "Intersection":I}
    self.padjustForHypergeom(result)
    return result


  def padjustForHypergeom(self,result):
    keys = list(result.keys())
    pvalues = [result[k]["p-value"] for k in keys]
    b,adj=multitest.fdrcorrection(pvalues, alpha=0.05, method='indep')
    for i in range(len(keys)):
      result[keys[i]]["padjust"]=adj[i]



  def saveHyperGeometricTest(self,results,ORFNameDict,DESCR):
    #GoTermsDescription = loadGoTermsDescriptions("GO_terms_used_in_H37Rv.csv")
    # f = open(fileNameOut,"w")
    for k in results:
      cad = " ".join([x+"/"+ORFNameDict[x] for x in results[k]["Intersection"]])
      self.output.write(k+"\t"+DESCR[k]+"\t"+results[k]["parameters"]+"\t"+str(results[k]["p-value"])+"\t"+str(results[k]["padjust"])+"\t"+cad+"\n")

  # HYPERGEOMETRIC END

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
    
    
  def GSEA3(self,resampling_file,GO_associations_file,output_file,p=1,N=100):

    GOterms,ontology,GOrvs = [],{},{}
    for line in open(GO_associations_file):
      w = line.rstrip().split('\t')
      ontology[w[0]] = w[1]
      rvs  = w[2].split()
      GOterms.append(w[0])
      GOrvs[w[0]] = rvs

    data,genenames,headers = [],{},[]
    for line in open(resampling_file):
      if line[0]=='#': headers = line.rstrip().split('\t'); continue # last header contains col names
      w = line.rstrip().split('\t')
      data.append(w)
      genenames[w[0]] = w[1]
    n2 = int(len(data)/2)

    pairs = [] # pair are: rv and (LFC or Zscore)
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
    
    results,Total = [],len(GOterms)
    
    for i,go in enumerate(GOterms):
      sys.stdout.flush()
      rvs = GOrvs[go]
      ranks = [MainIndex.get(x,2000) for x in rvs]
      es = self.enrichment_score(rvs,MainIndex,sorted_scores) # always positive, even if negative deviation, since I take abs
      mr = self.mean_rank(rvs,MainIndex)
      higher = 0
      for n,perm in enumerate(indexes):
        if self.enrichment_score(rvs,perm,sorted_scores)>es: higher += 1
        if n>100 and higher>10: break # adaptive
      pval = higher/float(n)
      vals = ['#',go,len(rvs),mr,es,pval,ontology.get(go,"?")]
      #sys.stderr.write(' '.join([str(x) for x in vals])+'\n')
      pctg=(100.0*i)/Total
      text = "Running Pathway Enrichment Method... %5.1f%%" % (pctg)
      self.progress_update(text, i)      
      results.append((go,mr,es,pval))
    
    results.sort(key=lambda x: x[-1])
    
    pvals = [x[-1] for x in results]
    rej,qvals = multitest.fdrcorrection(pvals)
    results = [tuple(list(res)+[q]) for res,q in zip(results,qvals)]

    self.output.write('\t'.join("GO_term description num_genes mean_rank enrichment_score pval qval".split())+'\n')
    for go,mr,es,pval,qval in results:
      rvs = GOrvs[go]
      rvinfo = [(x,genenames.get(x,"?"),MainIndex.get(x,n2)) for x in rvs]
      rvinfo.sort(key=lambda x: x[2])
      rvs = ["%s/%s (%s)" % x for x in rvinfo]
      rvs = ' '.join(rvs)
      vals = [go,ontology.get(go,"?"),len(GOrvs[go]),"%0.1f" % mr]+["%0.6f" % x for x in [es,pval,qval]]+[rvs]
      self.output.write('\t'.join([str(x) for x in vals])+'\n')
    self.output.close()


  #################Z Test#################################
  # D is the whole set with p values <=0.05
  # S is the Sanger Categories with their genes.
  
  def Ztest(self,D,S):    
    auxD = [d[0] for d in D]    
    DictD={d[0]:d[1] for d in D}
    ES={}
    rank={}    
    for key in S: #GoTerms or Category in S      
      r=[DictD[k] for k in S[key] if k in DictD]
      lenr=len(r)
      if(lenr>1):        
        meanR = numpy.mean(r)
        es = math.sqrt(lenr)*meanR
        ES[key]= [lenr,es,norm.sf(es)]        
        rank[key]={i: auxD.index(i) for i in S[key] if i in auxD}
    self.padjustTest(ES)    
    return ES, rank

  def ChiSquareTest(self,D,S):    
    auxD = [d[0] for d in D]
    DictD={d[0]:d[1] for d in D}
    ES={}
    rank={}

    for key in S: #GoTerms or Category in S
      r=[DictD[k] for k in S[key] if k in DictD]
      lenr=len(r)
      if(lenr>19):
        meanR = numpy.mean(r)      
        lenSkey = lenr-1
        es=(sum([(ri-meanR)**2 for ri in r]) - (lenSkey)) / (2*lenSkey)
        ES[key]= [lenr,es, norm.sf(es)]
        rank[key]={i: auxD.index(i) for i in S[key] if i in auxD}
    self.padjustTest(ES)
    return ES,rank


  def printTest(self,test,ORFNameDict,rank,DESCR):    
    for key in test:
      cad = " ".join([str(rank[key][g])+":"+g+"/"+ORFNameDict[g] for g in rank[key]])
      self.output.write(key+"\t"+DESCR[key]+"\t"+str(test[key][0])+"\t"+str(test[key][1])+"\t"+str(test[key][2])+"\t"+str(test[key][3])+"\t"+cad+"\n")
    self.output.close()  

  def padjustTest(self,test):
    keys = list(test.keys())
    pvals = [test[k][2] for k in keys]
    b,adj=multitest.fdrcorrection(pvals, alpha=0.05, method='indep')
    for i in range(len(keys)):
      test[keys[i]]+=[adj[i]]

  def t_ishEstimators(self,D):
    lenD = len(D)
    for i in range(lenD):
      if D[i][1] == 0:
        D[i][1]=-4.0
      elif D[i][1] == 1:
        D[i][1]=3.0
      else:
        D[i][1]=norm.ppf(D[i][1])
####################################################################################
  @classmethod
  def fromGUI(self, wxobj):
    pass

  @classmethod
  def fromargs(self, rawargs): 
    (args, kwargs) = transit_tools.cleanargs(rawargs)

    p=int(kwargs.get("p", 1))
    N=int(kwargs.get("S", 1000))
    M = kwargs.get("M", "GSEA")
    if "Z" in kwargs: self.useZscores = True
    else: self.useZscores = False
    resamplingFile = args[0]    
    geneSetFile = args[1]
    outpath = args[2]
    output = open(outpath, "w") #The outputfile is opened here!!!

    
    return self(resamplingFile,
        geneSetFile,
        output,
        p,N,M)

  def Run(self):
    self.transit_message("Starting Pathway Enrichment Method")
    start_time = time.time()

    if self.M =="GSEA":
      self.GSEA3(self.resamplingFile,self.geneSetFile,self.output,p=0,N=10000) # what if self.N is set by kwargs?
      sys.exit(0)

    elif self.M =="HYPE":
      method = "HYPERGEOMETRIC"
      D,ORFNameDict = self.loadD(self.resamplingFile)
      GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)
      DE = [d[0] for d in D if d[1]<=0.05]
      S = len(DE) #Sample Size in hyperGeometric  
      W = len(D) #Whole Genes      
      # k is the number of genes in DE that are in GoTermsWithRv
      results=self.hyperGeometricTest(DE,W,GoTermsWithRV,S)

    elif self.M=="GSEA-Z":
      method = "GSEA-Z"
      D,ORFNameDict = self.loadD(self.resamplingFile)
      GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)

      self.t_ishEstimators(D)
      results,rank = self.Ztest(D,GoTermsWithRV)

    elif self.M=="GSEA-CHI":
      method = "GSEA-CHI"
      D,ORFNameDict = self.loadD(self.resamplingFile)
      GoTermsWithRV,DESCR = self.loadGoTermsGoTermAsKey(self.geneSetFile)
      self.t_ishEstimators(D)
      results,rank= self.ChiSquareTest(D,GoTermsWithRV)

    else:
      method = "Not a valid option"
      self.progress_update("Not a valid option", 100)

    self.output.write("#Pathway Enrichment\n")
    if self.wxobj:
      members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
      memberstr = ""
      for m in members:
        memberstr += "%s = %s, " % (m, getattr(self, m))
      self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata).encode('utf-8'), self.annotation_path.encode('utf-8'), self.output.name.encode('utf-8')))
    else:
      self.output.write("#Console: python %s\n" % " ".join(sys.argv))

    self.output.write("#Data: %s\n" % self.resamplingFile) 
    self.output.write("#Annotation path: %s\n" % self.annotation_path.encode('utf-8')) 
    self.output.write("#Time: %s\n" % (time.time() - start_time))
    self.output.write("#Methodology: %s\n"%method)
    if self.M =="GSEA":
      columns = ["[ID][descr]","Total genes","score","pval","padj","rank of genes"]
      self.output.write("#%s\n" % "\t".join(columns))
      self.saveExit(gseaVal,PathDict,rank,ORFNameDict,DESCR)
    elif self.M =="HYPE":
      columns=["cat id", "descr","Total genes","Total in intersection","pval","padj","genes in intersection"]
      self.output.write("#%s\n" % "\t".join(columns))
      self.saveHyperGeometricTest(results,ORFNameDict,DESCR)
    elif self.M =="GSEA-Z":
      columns=["#ID-Description","Total Genes","Score","P-Value","P-Adjust","genes"]
      self.output.write("#%s\n" % "\t".join(columns))
      self.printTest(results,ORFNameDict,rank,DESCR)
    elif self.M =="GSEA-CHI":      
      columns=["#ID-Description","Total Genes","Score","P-Value","P-Adjust","genes"]
      self.output.write("#%s\n" % "\t".join(columns))
      self.printTest(results,ORFNameDict,rank,DESCR)
    self.transit_message("") # Printing empty line to flush stdout 
    self.transit_message("Adding File: %s" % (self.output.name))
    self.add_file(filetype="Pathway Enrichment")
    self.finish()
    self.transit_message("Finished Pathway Enrichment Method") 


  @classmethod
  def usage_string(self):
    return """python %s pathway_enrichment <resampling_file> <pathway_associations_file> <output_file> [-p <float> -Z -S <int> -M <GSEA|HYPE|Z|CHI>] """ % (sys.argv[0])

if __name__ == "__main__":

  (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

  print("ARGS:", args)
  print("KWARGS:", kwargs)

  G = GSEAMethod.fromargs(sys.argv[1:])

  print(G)
  G.console_message("Printing the member variables:")   
  G.print_members()

  print("")
  print("Running:")

  G.Run()



