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

from pytransit.export import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools


############# Description ##################

short_name = "combined_wig"
long_name = "Method to export datasets in 'Combined Wig' format."
description = "A method to export and normalized datasets in 'Combined Wig' format."
label = "to Combined Wig"
transposons = ["himar1", "tn5"]

############# Analysis Method ##############

class CombinedWigExport(base.TransitExport):
    def __init__(self):
        base.TransitExport.__init__(self, short_name, long_name, description, label, transposons, CombinedWigMethod, CombinedWigGUI,)


################# GUI ##################

class CombinedWigGUI(base.ExportGUI):

    def __init__(self):
        base.ExportGUI.__init__(self)

########## METHOD #######################

class CombinedWigMethod(base.SingleConditionMethod):
    """
    CombinedWig

    """
    def __init__(self,
                ctrldata,
                annotation_path,
                output_file,
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, description, label, ctrldata, annotation_path, output_file, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)




    @classmethod
    def fromGUI(self, wxobj):
        """ """

        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        #Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata, transposons):
            return None

        # Choose normalization method
        normalization = wxobj.chooseNormalization()


        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0


        #Get output path
        defaultFileName = "combined_wig_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(ctrldata,
                annotationPath,
                output_file,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (len(args) != 3): # wigs prot_table output
            print("Error: Incorrect number of args. See usage")
            print(self.usage_string())
            sys.exit(0)

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0

        return self(ctrldata,
                annotationPath,
                output_file,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting Combined Wig Export")
        start_time = time.time()

        #Get orf data
        self.transit_message("Getting Data")
        (fulldata, position) = tnseq_tools.get_data(self.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(fulldata, self.normalization,
            self.ctrldata, self.annotation_path)
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.annotation_path)
        rv2info = transit_tools.get_gene_info(self.annotation_path)

        self.transit_message("Normalizing")
        self.output.write("#Converted to CombinedWig with TRANSIT.\n")
        self.output.write("#normalization method: %s\n" % self.normalization)
        if self.normalization != "nonorm":
            if type(factors[0]) == type(0.0):
                self.output.write("#Normalization Factors: %s\n" % "\t".join(["%s" % f for f in factors.flatten()]))
            else:
                self.output.write("#Normalization Factors: %s\n" % " ".join([",".join(["%s" % bx for bx in b]) for b in factors]))

        # get ref genome from first wig file (we could check that they are all the same)
        # by this point, we know all the wig files exist and have same number of TA sites
        # this assumes 2nd line of each wig file is "variableStep chrom=XXX"; if not, set ref to "unknown"
        self.ref = "unknown"
        wig = self.ctrldata[0]
        with open(wig) as file:
          line = file.readline()
          line = file.readline()
          if line.startswith("variableStep"): 
            w = line.rstrip().split("=")
            if len(w)>=2: self.ref = w[1]

        (K,N) = fulldata.shape
        self.output.write("#RefGenome: %s\n" % self.ref)
        for f in self.ctrldata:
            self.output.write("#File: %s\n" % f)
        self.output.write("#TA_coord\t%s\n" % ('\t'.join(self.ctrldata)))

        for i,pos in enumerate(position):
            #self.output.write("%d\t%s\t%s\n" % (position[i], "\t".join(["%1.1f" % c for c in fulldata[:,i]]),",".join(["%s (%s)" % (orf,rv2info.get(orf,["-"])[0]) for orf in hash.get(position[i], [])])   ))
            if self.normalization!='nonorm': vals = "\t".join(["%1.1f" % c for c in fulldata[:,i]])
            else: vals = "\t".join(["%d" % c for c in fulldata[:,i]]) # no decimals if raw counts
            self.output.write("%d\t%s\t%s\n" % (position[i],vals,",".join(["%s (%s)" % (orf,rv2info.get(orf,["-"])[0]) for orf in hash.get(position[i], [])])   ))
            # Update progress
            text = "Running Export Method... %5.1f%%" % (100.0*i/N)
            if i%1000==0: self.progress_update(text, i)
        self.output.close()



        self.transit_message("") # Printing empty line to flush stdout
        self.finish()
        self.transit_message("Finished Export")

#

    @classmethod
    def usage_string(self):
        return """python %s export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file> [-n normalization_method]\ndefault normalization_method=TTR""" % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    print("ARGS:", args)
    print("KWARGS:", kwargs)

    G = Example.fromargs(sys.argv[1:])

    print(G)
    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()


