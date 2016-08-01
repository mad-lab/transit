import sys

try:
    import wx
    hasWx = True
    #Check if wx is the newest 3.0+ version:
    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False
except Exception as e:
    hasWx = False
    newWx = False

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

#method_name = "griffin"


############# GUI ELEMENTS ##################

short_name = "griffin"
long_name = "Basic frequentist analysis of essentiality using gaps."
description = "Analysis of gaps used in Griffin et al. 2011"
transposons = ["himar1"]
columns = ["Orf","Name","Desc","k","n","r","s","t","Expected Run","p-value", "p-adjusted"]



############# Analysis Method ##############

class GriffinAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, description, transposons, GriffinMethod, GriffinGUI, [GriffinFile])


################## FILE ###################

class GriffinFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#Griffin", columns)

    def getHeader(self, path):
        ess=0; unc=0; non=0; short=0
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if float(tmp[-1]) < 0.05:
                ess+=1
            else:
                non+=1

        text = """Results:
    Essentials: %s
    Non-Essential: %s
            """ % (ess,non)
        return text


################## GUI ###################

class GriffinGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        griffinPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        griffinSection = wx.BoxSizer( wx.VERTICAL )

        griffinLabel = wx.StaticText( griffinPanel, wx.ID_ANY, u"griffin Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        griffinLabel.Wrap( -1 )
        griffinSection.Add( griffinLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        griffinSizer1 = wx.BoxSizer( wx.HORIZONTAL )

        griffinSection.Add( griffinSizer1, 1, wx.EXPAND, 5 )

        griffinButton = wx.Button( griffinPanel, wx.ID_ANY, u"Run griffin", wx.DefaultPosition, wx.DefaultSize, 0 )
        griffinSection.Add( griffinButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        griffinPanel.SetSizer( griffinSection )
        griffinPanel.Layout()
        griffinSection.Fit( griffinPanel )

        #Connect events
        griffinButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )

        self.panel = griffinPanel





########## CLASS #######################

class GriffinMethod(base.SingleConditionMethod):
    """   
    griffin
 
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

        base.SingleConditionMethod.__init__(self, short_name, long_name, description, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)


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
        if not transit_tools.validate_filetypes(ctrldata, transposons):
            return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        replicates = "Sum"
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "griffin_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(ctrldata,
                annotationPath,
                output_file,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs): 
        (args, kwargs) = transit_tools.cleanargs(rawargs)


        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = "Sum"
        normalization = None
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

        self.transit_message("Starting Griffin Method")
        start_time = time.time()
       
 
        #Get orf data
        self.transit_message("Getting Data")
        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)

        N = len(G)
        self.progress_range(N)
        count = 0
        pins = G.global_theta()
        pnon = 1.0 - pins
        results = []
        for gene in G:
            if gene.n == 0:
                results.append([gene, 0.0, 1.000])
            else:
                B = 1.0/math.log(1.0/pnon)
                u = math.log(gene.n*pins, 1.0/pnon)
                exprun = tnseq_tools.ExpectedRuns(gene.n, pnon)
                pval = 1.0 - tnseq_tools.GumbelCDF(gene.r, u, B)
                results.append([gene, exprun, pval])

            text = "Running Griffin Method... %2.0f%%" % (100.0*(count+1)/(N))
            self.progress_update(text, count)
            self.transit_message_inplace(text)
            count+=1


        pval = [row[-1] for row in results]
        padj = stat_tools.BH_fdr_correction(pval)
        for i in range(len(results)):
            results[i].append(padj[i])
        results.sort()
        
        self.output.write("#Griffin\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata), self.annotation_path, self.output.name))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Annotation path: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))
        
        for (gene, exprun, pval, padj) in results:
            self.output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.1f\t%1.5f\t%1.5f\n" % (gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, gene.s, gene.t, exprun, pval, padj))

        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Griffin")
        self.finish()
        self.transit_message("Finished Griffin Method")

    @classmethod
    def usage_string(self):
        return """python %s griffin <comma-separated .wig files> <annotation .prot_table> <output file>""" % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)


    G = GriffinMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


