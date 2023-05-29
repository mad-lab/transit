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

#method_name = "griffin"


############# GUI ELEMENTS ##################

short_name = "griffin"
long_name = "Griffin"
short_desc = "Basic frequentist analysis of essentiality using gaps."
long_desc = "Analysis of gaps used in Griffin et al. 2011"
transposons = ["himar1"]
columns = ["Orf","Name","Desc","k","n","r","s","t","Expected Run","p-value", "p-adjusted"]



############# Analysis Method ##############

class GriffinAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, GriffinMethod, GriffinGUI, [GriffinFile])


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

        griffinLabel = wx.StaticText( griffinPanel, wx.ID_ANY, u"griffin Options", wx.DefaultPosition, (120,-1), 0 )
        griffinLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
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
                minread=1,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, ignoreCodon=ignoreCodon, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.minread = minread


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


        #
        minread = 1

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
                minread,
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

        minread = int(kwargs.get("m", 1))
        replicates = kwargs.get("r", "Sum")
        normalization = None
        LOESS = False
        ignoreCodon = not kwargs.get("sC", False)
        NTerminus = float(kwargs.get("iN", 0.0))
        CTerminus = float(kwargs.get("iC", 0.0))


        return self(ctrldata,
                annotationPath,
                output_file,
                minread,
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

        (data, position) = transit_tools.get_validated_data(self.ctrldata, wxobj=self.wxobj)
        (K,N) = data.shape

        if self.normalization and self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata, self.annotation_path)

        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, minread=1, reps=self.replicates, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data, position=position)





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

            text = "Running Griffin Method... %5.1f%%" % (100.0*(count+1)/(N))
            self.progress_update(text, count)
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
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata).encode('utf-8'), self.annotation_path.encode('utf-8'), self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode('utf-8'))) 
        self.output.write("#Annotation path: %s\n" % self.annotation_path.encode('utf-8')) 
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
        return """python3 %s griffin <comma-separated .wig files> <annotation .prot_table> <output file> [Optional Arguments]

        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -sC             :=  Include stop-codon (default is to ignore).
        -iN <float>     :=  Ignore TAs occuring at given fraction (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given fraction (as integer) of the C terminus. Default: -iC 0
        """ % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)


    G = GriffinMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print("")
    print("Running:")

    G.Run()


