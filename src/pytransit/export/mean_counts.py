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

short_name = "mean_counts"
long_name = "Method to export datasets in 'Mean Gene Counts' format."
description = "A method to export and normalized datasets in 'Mean Gene Counts' format."
label = "to Mean Gene Counts"
transposons = ["himar1", "tn5"]

############# Analysis Method ##############

class MeanCountsExport(base.TransitExport):
    def __init__(self):
        base.TransitExport.__init__(self, short_name, long_name, description, label, transposons, MeanCountsMethod, MeanCountsGUI,)


################# GUI ##################

class MeanCountsGUI(base.ExportGUI):

    def __init__(self):
        base.ExportGUI.__init__(self)

########## METHOD #######################

class MeanCountsMethod(base.SingleConditionMethod):
    """   
    Mean Gene Counts
 
    """
    def __init__(self,
                combined_wig,
                ctrldata,
                annotation_path,
                output_file,
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, description, label, ctrldata, annotation_path, output_file, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.combined_wig = combined_wig # boolean, interprete ctrldata as combined_wig file or comma-separated list of wig files?


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
        defaultFileName = "gene_mean_counts_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")

        # could add a checkbox for combined_wig
        combined_wig = False

        return self(combined_wig,ctrldata,
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
        print("ARGS="+str(args))
        print("KWARGS="+str(kwargs))

        combined_wig = kwargs.get("c",False)
        ctrldata = args[0].split(",") 
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        normalization = kwargs.get("n", "TTR")
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0

        return self(combined_wig,ctrldata,
                annotationPath,
                output_file,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting Gene Mean Counts Export")
        start_time = time.time()
        
        #Get orf data
        self.transit_message("Getting Data")
        if self.combined_wig: (position, fulldata, datasets) = tnseq_tools.read_combined_wig(self.ctrldata[0])
        else: (fulldata, position) = tnseq_tools.get_data(self.ctrldata)
        (fulldata, factors) = norm_tools.normalize_data(fulldata, self.normalization, self.ctrldata, self.annotation_path)
        position = position.astype(int)

        hash = transit_tools.get_pos_hash(self.annotation_path)
        rv2info = transit_tools.get_gene_info(self.annotation_path)

        self.transit_message("Normalizing")
        self.output.write("#Summarized to Mean Gene Counts with TRANSIT.\n")
        if self.normalization != "nonorm":
            self.output.write("#Reads normalized using '%s'\n" % self.normalization)
            if type(factors[0]) == type(0.0):
                self.output.write("#Normalization Factors: %s\n" % "\t".join(["%s" % f for f in factors.flatten()]))
            else:
                self.output.write("#Normalization Factors: %s\n" % " ".join([",".join(["%s" % bx for bx in b]) for b in factors]))


        self.output.write("#Files:\n")
        names = datasets if self.combined_wig else self.ctrldata 
        for f in names:
            self.output.write("#%s\n" % f)


        K,Nsites = fulldata.shape
        # Get Gene objects
        if self.combined_wig: G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, norm=self.normalization,data=fulldata,position=position)
        else: G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, norm=self.normalization)
        N = len(G)
        self.progress_range(N)
        if self.combined_wig: dataset_header = '\t'.join(datasets)
        else: dataset_header = "\t".join([transit_tools.fetch_name(D) for D in self.ctrldata])
        self.output.write("#Orf\tName\tNumber of TA sites\t%s\n" % dataset_header)
        for i,gene in enumerate(G):
            if gene.n > 0:
                data_str = "\t".join(["%1.2f" % (M) for M in numpy.mean(gene.reads, 1)])
            else:
                data_str = "\t".join(["%1.2f" % (Z) for Z in numpy.zeros(K)])
            self.output.write("%s\t%s\t%s\t%s\n" % (gene.orf, gene.name, gene.n, data_str))

            # Update progress
            text = "Running Export Method... %5.1f%%" % (100.0*i/N)
            self.progress_update(text, i)
        self.output.close()



        self.transit_message("") # Printing empty line to flush stdout 
        self.finish()
        self.transit_message("Finished Export") 

#

    @classmethod
    def usage_string(self):
        return """python %s export mean_counts <comma-separated .wig files>|<combined_wig> <annotation .prot_table> <output file> [-c]\n note: append -c if inputing a combined_wig file\n""" % (sys.argv[0])


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


