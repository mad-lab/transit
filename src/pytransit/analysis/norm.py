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

short_name = "norm"
long_name = "Normalization"
short_desc = "Normalization method"
long_desc = "Method for normalizing datasets and outputting into CombinedWig file."
transposons = ["himar1", "tn5"]
columns = ["Position","Reads","Genes"]


############# Analysis Method ##############

class NormAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, NormMethod, NormGUI, [NormFile])


################## FILE ###################

class NormFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#CombinedWig", columns)

    def getHeader(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################

class NormGUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################


class NormMethod(base.SingleConditionMethod):
    """   
    Norm
 
    """
    def __init__(self,
                ctrldata,
                annotation_path,
                output_file,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)


    @classmethod
    def fromargs(self, rawargs): 
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if len(args) < 3:
            raise base.InvalidArgumentException("Must provide all necessary arguments")
            

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = "Sum"
        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0

        return self(ctrldata,
                annotationPath,
                output_file,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting Normalization")
        start_time = time.time()
        outputPath = self.output.name
        # Normalize Data
        transit_tools.convertToCombinedWig(self.ctrldata, self.annotation_path, outputPath, normchoice=self.normalization)       

        self.finish()
        self.transit_message("Finished Normalization") 

    @classmethod
    def usage_string(self):
        return """
python3 %s norm <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    G = Norm.fromargs(sys.argv[1:])
    G.Run()


