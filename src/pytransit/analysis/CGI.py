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
import datetime

from pytransit.analysis import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools

hasR = False
try:
    import rpy2.robjects
    hasR = True
except Exception as e:
    hasR = False

if hasR:
    from rpy2.robjects import r, DataFrame, globalenv, IntVector, FloatVector, StrVector, packages as rpackages

############# Description ##################

short_name = "CGI"
long_name = "Chemical Genetic Analysis"
short_desc = "CGI Analysis of CRISPRi libraries"
long_desc = "CGI Analysis of CRISPRi libraries"
transposons = []

columns = ["Position","Reads","Genes"] # ???

############# Analysis Method ##############

class CGI(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, CGI_Method, CGI_GUI, []) 

################## FILE ###################

# there is no output file that could be loaded into the GUI

#class HeatmapFile(base.TransitFile):
#
#    def __init__(self):
#        base.TransitFile.__init__(self, "#CombinedWig", columns) 
#
#    def getHeader(self, path):
#        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
#        return text

################# GUI ##################

# right now, CGI is just intended for the command-line; TRI

class CGI_GUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################

class CGI_Method(base.SingleConditionMethod):
    def __init__(self):
                ctrldata=None # initializers for superclass
                annotation_path=""
                output_file=""
                replicates="Sum"
                normalization="nonorm" 
                LOESS=False
                ignoreCodon=True
                NTerminus=0.0
                CTerminus=0.0
                wxobj=None
                # this initialization seems pointless for CGI, but must do this for base class...
                base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

    @classmethod
    def usage_string(self):
        return """usage (3 sub-commands):
  python3 ../src/transit.py CGI extract_abund <metadata_file> <data_dir> <betaE_file> <no_drug_file> <no_depletion_abundances_file> <drug> <days>  >  <output_file>
  python3 ../src/transit.py CGI run_model <abund_file>  >  <logsigmodfit_file>
  python3 ../src/transit.py CGI post_process <logsigmoidfit_file>  >  <results_file>
note: redirect output from stdout to output files as shown above"""


    @classmethod
    def fromargs(self, rawargs): 
        if not hasR:
            print("Error: R and rpy2 (~= 3.0) required to run heatmap.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)
        if len(args)<1: print(self.usage_string())
        self.cmd = args[0]
        self.args = args[1:]
        self.kwargs = kwargs

        return self()

    def Run(self):
        cmd,args,kwargs = self.cmd,self.args,self.kwargs

        if cmd=="extract_abund":
          if len(args)<6: print(self.usage_string())

          metadata_file = args[0]
          data_dir = args[1]
          betaE_file = args[2]
          no_drug_file = args[3]
          no_dep_abund = args[4]
          drug = args[5]
          days = args[6]

          self.extract_abund(metadata_file,data_dir,betaE_file,no_drug_file,no_dep_abund,drug,days)
  
        else: print(self.usage_string())

    ######################################

    # alternatively, I could pass in a list of sample ids to load(including DMSO), 
    #   or provide a row selection string (like "days=5 and drug=RIF")

    # PC=pseudocounts for fractional abundances 
    #   (important, since it sets a lower bound on vals, below which is assumed to be noise)
    #   how will this differ between libraries or sequencing runs?

    def extract_abund(self,metadata_file,data_dir,betaE_file,no_drug_file,no_dep_abund,drug,days,PC=1e-8):
      print("in extract_abund: betaE_file=%s, days=%s" % (betaE_file,days))

      #################
      # read in all the files with supporting data
  
      betaE = {}
      skip = 1
      for line in open(betaE_file): # Bosch21_TableS2.txt
        if skip>0: skip -= 1; continue
        w = line.rstrip().split('\t')
        if len(w[9])<2: continue # some betaE entries are empty
        id,b = w[0],float(w[9])
        id = id[:id.rfind("_")] # strip off v4PAMscore...
        betaE[id] = b
      
      metadata = Spreadsheet(metadata_file) # ShiquiCGI_metadata.txt
      
      files = []
      concs = []
      for i in range(metadata.nrows):
        if metadata.get(i,"drug")==drug and metadata.get(i,"days_predep")==days:
          files.append(metadata.get(i,"filename"))
          concs.append(float(metadata.get(i,"conc_xMIC")))
      
      # also need to get DMSO to represent 0xMIC
      # I should make the DMSO file an input arg too...
      
      #if days=="1": files.append("counts_1962_DMSO_D1.txt"); concs.append(0) # what about different pools?
      #if days=="5": files.append("counts_1972_DMSO_D5.txt"); concs.append(0)
      #if days=="10": files.append("counts_1952_DMSO_D10.txt"); concs.append(0)
      files.append(no_drug_file); concs.append(0)
      
      pairs = sorted(zip(concs,files))
      concs,files = [x[0] for x in pairs],[x[1] for x in pairs]
      
      for i in range(len(files)):
        sys.stderr.write("%s\t%s\n" % (concs[i],files[i]))
      
      no_dep = {}
      IDs = []
      for line in open(no_dep_abund): # no_depletion_abundances.txt
        w = line.rstrip().split('\t')
        id,abund = w[0],float(w[1])
        id = id[:id.rfind("_")]
        no_dep[id] = abund
        IDs.append(id)
      
      #################
      # process the counts, calculate fractional normalized abundances, and generate output
      
      # assume there are 3 replicates in each file
      
      Abund = []
      Concs = [] # one for each replicate, include those for DMSO (0)
      for file,conc in zip(files,concs):
        counts = {} # rows of 3, indexed by id
        skip = 1
        for line in open("%s/%s" % (data_dir,file)):
          if skip>0: skip -= 1; continue
          if "Negative" in line or "Empty" in line: continue
          w = line.rstrip().split('\t')
          id = w[0]
          id = id[:id.rfind("_")]
          vals = [int(x) for x in w[1:]]
          counts[id] = vals # rows of 3
        for i in range(3): # assume there are 3 replicates in each file
          cnts = [x[i] for x in counts.values()] # a column
          tot = sum(cnts)
          abund = []
          for id in IDs: abund.append(counts[id][i]/float(tot)) # a column for fracs, parallel to IDs
          Abund.append(abund)
          Concs.append(conc)
      
      # assumes that sgRNA ids embed orf and gene ids, e.g. "RVBD00067:rpoB" 
      
      print('\t'.join("orf gene id abund betaE".split()+[str(x) for x in Concs]))
      a,b = 0,0
      for i,id in enumerate(IDs):
        #vals = [id]+["%0.8f" % x[i] for x in Abund]
        if id not in betaE: a += 1; continue; #sys.stderr.write("%s not in betaE\n" % id); continue
        if id not in no_dep: b += 1; continue # sys.stderr.write("%s not in no_dep\n" % id); continue
        temp = id[:id.find("_")].split(":") # assumes that sgRNA ids embed orf and gene ids, e.g. "RVBD00067:rpoB" 
        orf,gene = temp[0],temp[1]
        vals = [orf,gene,id,"%0.8f" % (no_dep[id]),str(betaE[id])]+["%0.6f" % ((x[i]+PC)/(no_dep[id]+PC)) for x in Abund]
        print('\t'.join(vals))
      
      sys.stderr.write("warning: betaE values not found for %s gRNAs\n" % a)
      sys.stderr.write("warning: no_dep values not found for %s gRNAs\n" % b)

  #####################################################

  # derived from logsigmoidfit.R
  # see heatmap.py for example of how to put data in a pandas.DataFrame and call an R function like make_heatmapFunc()

  # def run_model(self): 

################################

class Spreadsheet:
  def __init__(self,filename): # ,key="Id"
    self.keys,self.data,self.headers = [],[],None
    self.rowhash,self.colhash = {},{}
    for line in open(filename,'rU'): # also handle mac files with x0d vs x0a
      if len(line)==0 or line[0]=='#': continue
      w = line.rstrip().split('\t') 
      if self.headers==None: 
        self.headers = w
        for i,h in enumerate(self.headers): 
          if h in self.colhash: print("error: '%s' appears multiple times in first row of '%s' - column headers must be unique" % (h,filename)); sys.exit(-1)
          self.colhash[h] = i
        continue
      self.data.append(w)
      key = w[0]
      if key=="": continue
      if key in self.rowhash: print("error: '%s' appears multiple times in first column of '%s' - keys must be unique" % (key,filename)); sys.exit(-1)
      self.rowhash[key] = len(self.keys)
      self.keys.append(key) # check if unique?
    self.ncols = max([len(row) for row in self.data])
    self.nrows = len(self.keys)
  def getrow(self,r):
    if r in self.rowhash: r = self.rowhash[r] 
    # check that it is an integer (in range)?
    return self.data[r]
  def getcol(self,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    return [r[c] for r in self.data] # what if not all same length?
  def get(self,r,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    if r in self.rowhash: r = self.rowhash[r]
    return self.data[r][c]

################################

if __name__ == "__main__":

    G = CGI.fromargs(sys.argv[1:])
    G.Run()


