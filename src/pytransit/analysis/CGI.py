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
  python3 ../src/transit.py CGI extract_abund <metadata_file> <data_dir> <betaE_file> <no_depletion_abundances_file> <drug> <days>  >  <output_file>
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
          if len(args)<5: print(self.usage_string())
          metadata_file = args[0]
          data_dir = args[1]
          betaE_file = args[2]
          no_dep_abund = args[3]
          drug = args[4]
          days = args[5]
          self.extract_abund(metadata_file,data_dir,betaE_file,no_dep_abund,drug,days)
  
        else: print(self.usage_string())

    def extract_abund(self,metadata_file,data_dir,betaE_file,no_dep_abund,drug,days):
      print("in extract_abund: betaE_file=%s, days=%s" % (betaE_file,days))


    # see heatmap.py for example of how to put data in a pandas.DataFrame and call an R function like make_heatmapFunc()


if __name__ == "__main__":

    G = CGI.fromargs(sys.argv[1:])
    G.Run()


