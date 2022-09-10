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

short_name = "heatmap"
long_name = "Heatmap"
short_desc = "Heatmap among Conditions"
long_desc = "Heatmap among Conditions"
transposons = ["himar1", "tn5"]

columns = ["Position","Reads","Genes"] # ???

############# Analysis Method ##############

class Heatmap(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, HeatmapMethod, HeatmapGUI, []) 

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

# right now, tnseq_stats is just intended for the command-line; TRI

class HeatmapGUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################

# should Heatmap be a SingleConditionMethod? args like normalization are irrelevant

class HeatmapMethod(base.SingleConditionMethod):
    """   
    Norm
 
    """
    def __init__(self,gene_means,outfile): 
                ctrldata=None # initializers for superclass
                annotation_path=""
                output_file=outfile
                replicates="Sum"
                normalization="nonorm" 
                LOESS=False
                ignoreCodon=True
                NTerminus=0.0
                CTerminus=0.0
                wxobj=None
                base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)


    @classmethod
    def fromargs(self, rawargs): 
        if not hasR:
            print("Error: R and rpy2 (~= 3.0) required to run heatmap.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)
        if len(rawargs)<3: print(self.usage_string()); sys.exit(-1)
        self.filetype = None
        if kwargs.get("anova",False): self.filetype = "anova"
        elif kwargs.get("zinb",False): self.filetype = "zinb"
        else: print(self.usage_string()); sys.exit(-1)
        self.infile = args[0]
        self.outfile = args[1]
        self.qval = float(kwargs.get("qval",0.05))
        self.topk = int(kwargs.get("topk",-1))
        self.low_mean_filter = int(kwargs.get("low_mean_filter",5)) # filter out genes with grandmean<5 by default
        return self(self.infile,outfile=self.outfile)

    def Run(self):

        if self.filetype!="anova" and self.filetype!="zinb":
          print("filetype not recognized: %s" % self.filetype); sys.exit(-1)
        
        headers = None
        data,hits = [],[]
        n = -1 # number of conditions

        for line in open(self.infile):
          w = line.rstrip().split('\t')
          if line[0]=='#' or ('pval' in line and 'padj' in line): # check for 'pval' for backwards compatibility
            headers = w; continue # keep last comment line as headers
          # assume first non-comment line is header
          if n==-1: 
            # ANOVA header line has names of conditions, organized as 3+2*n+3 (2 groups (means, LFCs) X n conditions)
            # ZINB header line has names of conditions, organized as 3+4*n+3 (4 groups X n conditions)
            if self.filetype=="anova": n = int((len(w)-8)/2) 
            elif self.filetype=="zinb": n = int((len(headers)-6)/4) 
            headers = headers[3:3+n]
            headers = [x.replace("Mean_","") for x in headers]
          else:
            means = [float(x) for x in w[3:3+n]] # take just the columns of means
            lfcs = [float(x) for x in w[3+n:3+n+n]] # take just the columns of LFCs
            qval = float(w[-2])
            data.append((w,means,lfcs,qval))

        data.sort(key=lambda x: x[-1])
        hits,LFCs = [],[]
        for k,(w,means,lfcs,qval) in enumerate(data):
          if (self.topk==-1 and qval<self.qval) or (self.topk!=-1 and k<self.topk): 
            mm = round(numpy.mean(means),1)
            if mm<self.low_mean_filter: print("excluding %s/%s, mean(means)=%s" % (w[0],w[1],mm))
            else: hits.append(w); LFCs.append(lfcs)

        print("heatmap based on %s genes" % len(hits))
        genenames = ["%s/%s" % (w[0],w[1]) for w in hits]
        hash = {}
        headers = [h.replace("Mean_","") for h in headers]
        for i,col in enumerate(headers): hash[col] = FloatVector([x[i] for x in LFCs])
        df = DataFrame(hash)
        heatmapFunc = self.make_heatmapFunc()
        heatmapFunc(df,StrVector(genenames),self.outfile)

    def make_heatmapFunc(self):
      r('''
make_heatmap = function(lfcs,genenames,outfilename) { 
rownames(lfcs) = genenames
suppressMessages(require(gplots))
colors <- colorRampPalette(c("red", "white", "blue"))(n = 200)

C = length(colnames(lfcs))
R = length(rownames(lfcs))
W = 300+C*30
H = 300+R*15

png(outfilename,width=W,height=H)
#defaults are lwid=lhei=c(1.5,4)
#heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),lwid=c(2,6),lhei=c(0.1,2),trace="none",cexCol=1.4,cexRow=1.4,key=T) # make sure white=0
#heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # make sure white=0 # setting margins was causing failures, so remove it 8/22/21
heatmap.2(as.matrix(lfcs),col=colors,margin=c(12,12),trace="none",cexCol=1.2,cexRow=1.2,key=T) # actually, margins was OK, so the problem must have been with lhei and lwid
dev.off()
}
      ''')
      return globalenv['make_heatmap']

    @classmethod
    def usage_string(self):
        return "usage: python3 %s heatmap <anova_or_zinb_output> <heatmap.png> -anova|-zinb [-topk <int>] [-qval <float>] [-low_mean_filter <int>]\n note: genes are selected based on qval<0.05 by default" % sys.argv[0]


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    G = Norm.fromargs(sys.argv[1:])
    G.Run()


