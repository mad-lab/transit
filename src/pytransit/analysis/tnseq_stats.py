import sys

try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0
    print "EXCEPTION:", str(e)

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

import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools


############# Description ##################

short_name = "tnseq_stats"
long_name = "TnSeq Statistics"
short_desc = "Statistical Metrics for TnSeq datasets"
long_desc = "Statistical Metrics for TnSeq datasets"
transposons = ["himar1", "tn5"]
columns = ["Position","Reads","Genes"]


############# Analysis Method ##############

class TnseqStats(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, TnseqStatsMethod, TnseqStatsGUI, [TnseqStatsFile])


################## FILE ###################

# I am not sure what this means in the context of tnseq_stats(); is it referring to the output file? TRI

class TnseqStatsFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#CombinedWig", columns) 

    def getHeader(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################

# right now, tnseq_stats is just intended for the command-line; TRI

class TnseqStatsGUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################


class TnseqStatsMethod(base.SingleConditionMethod):
    """   
    Norm
 
    """
    def __init__(self,wigs): # list of wig files
                ctrldata=wigs # initializers for superclass
                annotation_path=""
                output_file=""
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
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if len(args) < 2:
            raise base.InvalidArgumentException("Must provide all necessary arguments")

        self.wigs = args
        return self(self.wigs)

    def Run(self):

        self.transit_message("Starting TnseqStats")
        start_time = time.time()

        (data, sites) = tnseq_tools.get_data(self.wigs)

        # should write this to an output file?
        print "dataset\tdensity\tmean_ct\tNZmean\tNZmedian\tmax_ct\ttotal_cts\tskewness\tkurtosis"
        for i in range(data.shape[0]):
          density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis = tnseq_tools.get_data_stats(data[i,:])
          nzmedianrd = int(nzmedianrd) if numpy.isnan(nzmedianrd)==False else 0
          vals = [self.wigs[i], "%0.2f" % density, "%0.1f" % meanrd, "%0.1f" % nzmeanrd, "%d" % nzmedianrd, maxrd, totalrd, "%0.1f" % skew, "%0.1f" % kurtosis]
          print '\t'.join([str(x) for x in vals])

        self.finish()
        self.transit_message("Finished TnseqStats") 

    @classmethod
    def usage_string(self):
        return """usage: python %s tnseq_stats <file1.wig> [<file2.wig>...]
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    G = Norm.fromargs(sys.argv[1:])
    G.Run()


