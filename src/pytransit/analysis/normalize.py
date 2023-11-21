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


############# Description ##################

short_name = "normalize"
long_name = "Normalize"
short_desc = "Normalization method"
long_desc = "Method for normalizing datasets."
transposons = ["himar1", "tn5"]
columns = ["Position","Reads","Genes"]


############# Analysis Method ##############

class Normalize(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, NormalizeMethod, NormalizeGUI, [NormalizeFile])


################## FILE ###################

# maybe this is only needed if it is going to get loaded into the output-files panel - TRI

class NormalizeFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#CombinedWig", columns)

    def getHeader(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################

class NormalizeGUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################


class NormalizeMethod(base.SingleConditionMethod):
    """
    Norm

    """
    def __init__(self,infile,outfile,normalization):
                ctrldata=[infile]
                annotation_path=""
                output_file=outfile
                replicates="Sum"
                normalization=normalization
                LOESS=False
                ignoreCodon=True
                NTerminus=0.0
                CTerminus=0.0
                wxobj=None
                base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)


    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        isCombinedWig = 'c' in kwargs
        if (not isCombinedWig and len(args) < 2) or (isCombinedWig and len(args) < 1):
            raise base.InvalidArgumentException("Must provide all necessary arguments")
        if isCombinedWig:
            self.infile = kwargs.get("c") # only 1 input wig file
            self.outfile = args[0] # if no arg give, could print to screen
        else:
            self.infile = args[0] # only 1 input wig file
            self.outfile = args[1] # if no arg give, could print to screen
        self.normalization = kwargs.get("n", "TTR") # should check if it is a legal method name
        self.combined_wig = isCombinedWig

        return self(self.infile,self.outfile,self.normalization)

    def Run(self):

        self.transit_message("Starting Normalization")
        start_time = time.time()

        infile = self.infile
        outputPath = self.outfile # output file exists, should I require -overwrite flag?

        # determine ref genome from first; assume they are all the same; assume wigs have 2 header lines
        line2 = "variableStep chrom=" # unknown
        for line in open(infile):
          if line.startswith("variableStep"): line2 = line.rstrip(); break

        if self.combined_wig==True: (sites,data,files) = tnseq_tools.read_combined_wig(self.ctrldata[0])
        else: (data, sites) = tnseq_tools.get_data(self.ctrldata)
        (data,factors) = norm_tools.normalize_data(data,self.normalization)

        print("writing",outputPath)
        file = open(outputPath,"w")
        file.write("# %s normalization of %s\n" % (self.normalization,infile))
        if self.combined_wig==True:
          for f in files: file.write("#File: %s\n" % f)
          for i in range(len(sites)): file.write('\t'.join([str(sites[i])]+["%0.1f" % x for x in list(data[...,i])])+"\n")
        else:
          file.write(line2+"\n")
          for j in range(len(sites)):
            file.write("%s %s\n" % (sites[j],int(data[0,j])))
        file.close()

        self.finish()
        self.transit_message("Finished Normalization")

    @classmethod
    def usage_string(self):
        return """
python3 %s normalize <input.wig> <output.wig> [-n TTR|betageom]
---
OR
---
python3 %s normalize -c <input combined_wig> <output.wig> [-n TTR|betageom]

        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        """ % (sys.argv[0], sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    G = Norm.fromargs(sys.argv[1:])
    G.Run()


