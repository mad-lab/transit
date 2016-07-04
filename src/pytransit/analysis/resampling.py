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
import ntpath
import math
import random
import numpy
import scipy.stats
import datetime

import matplotlib.pyplot as plt

import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "resampling"
long_name = "Resampling test of conditional essentiality between two conditions"
description = """Method for determining conditional essentiality based on resampling (i.e. permutation test). Identifies significant changes in mean read-counts for each gene after normalization."""

transposons = ["himar1", "tn5"]
columns = ["Orf","Name","Desc","Sites","Mean A","Mean B","Delta sum","log2FC","pvalue","adj. pvalue"]

class ResamplingAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, description, transposons, ResamplingMethod, ResamplingGUI, [ResamplingFile])



############# FILE ##################

class ResamplingFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#Resampling", columns)

    def getHeader(self, path):
        DE=0; poslogfc=0; neglogfc=0;
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if float(tmp[-1]) < 0.05:
                DE +=1
                if float(tmp[-3]) > 0:
                    poslogfc+=1
                else:
                    neglogfc+=1

        text = """Results:
    Conditionally - Essentials: %s
        More Essential in Experimental datasets: %s
        Less Essential in Experimental datasets: %s
            """ % (DE, poslogfc, neglogfc)
        return text


    def getMenus(self):
        menus = []
        menus.append(("Display in Track View", self.displayInTrackView))
        menus.append(("Display Histogram", self.displayHistogram))
        return menus

    def displayHistogram(self, displayFrame, event):
            gene = displayFrame.grid.GetCellValue(displayFrame.row, 0)
            filepath = os.path.join(ntpath.dirname(displayFrame.path), transit_tools.fetch_name(displayFrame.path))
            filename = os.path.join(filepath, gene+".png")
            if os.path.exists(filename):
                imgWindow = pytransit.fileDisplay.ImgFrame(None, filename)
                imgWindow.Show()
            else:
                ShowError(MSG="Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
                print "Error Displaying File. Histogram image does not exist."


        

############# GUI ##################

class ResamplingGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        resamplingPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        resamplingSizer = wx.BoxSizer( wx.VERTICAL )

        resamplingLabel = wx.StaticText( resamplingPanel, wx.ID_ANY, u"resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingLabel.Wrap( -1 )
        resamplingSizer.Add( resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        resamplingTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        resamplingTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        resamplingLabelSizer = wx.BoxSizer( wx.VERTICAL )

        # Samples Label
        resamplingSampleLabel = wx.StaticText( resamplingPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSampleLabel.Wrap( -1 )
        resamplingLabelSizer.Add( resamplingSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        # Norm Label
        resamplingNormLabel = wx.StaticText( resamplingPanel, wx.ID_ANY, u"Normalization", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingNormLabel.Wrap( -1 )
        resamplingLabelSizer.Add( resamplingNormLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )


        resamplingTopSizer2.Add( resamplingLabelSizer, 1, wx.EXPAND, 5 )

        resamplingControlSizer = wx.BoxSizer( wx.VERTICAL )

        # Samples Text
        self.wxobj.resamplingSampleText = wx.TextCtrl( resamplingPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingControlSizer.Add( self.wxobj.resamplingSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        # Norm Choices
        resamplingNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        self.wxobj.resamplingNormChoice = wx.Choice( resamplingPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, resamplingNormChoiceChoices, 0 )
        self.wxobj.resamplingNormChoice.SetSelection( 0 )
        resamplingControlSizer.Add( self.wxobj.resamplingNormChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )


        # Adaptive Check
        self.wxobj.resamplingAdaptiveCheckBox = wx.CheckBox(resamplingPanel, label = 'Adaptive Resampling (Faster)')

        # Histogram Check
        self.wxobj.resamplingHistogramCheckBox = wx.CheckBox(resamplingPanel, label = 'Generate Resampling Histograms')

        # Zeros Check
        self.wxobj.resamplingZeroCheckBox = wx.CheckBox(resamplingPanel, label = 'Include sites with all zeros')


        resamplingTopSizer2.Add( resamplingControlSizer, 1, wx.EXPAND, 5 )

        resamplingTopSizer.Add( resamplingTopSizer2, 1, wx.EXPAND, 5 )


        resamplingSizer.Add( resamplingTopSizer, 1, wx.EXPAND, 5 )
        resamplingSizer.Add( self.wxobj.resamplingAdaptiveCheckBox, 0, wx.EXPAND, 5 )
        resamplingSizer.Add( self.wxobj.resamplingHistogramCheckBox, 0, wx.EXPAND, 5 )
        resamplingSizer.Add( self.wxobj.resamplingZeroCheckBox, 0, wx.EXPAND, 5 )

        resamplingButton = wx.Button( resamplingPanel, wx.ID_ANY, u"Run resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSizer.Add( resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

 
        resamplingPanel.SetSizer( resamplingSizer )
        resamplingPanel.Layout()
        resamplingSizer.Fit( resamplingPanel )

        #Connect events
        resamplingButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )

        self.panel = resamplingPanel



########## CLASS #######################

class ResamplingMethod(base.DualConditionMethod):
    """   
    resampling
 
    """
    def __init__(self,
                ctrldata,
                expdata,
                annotation_path,
                output_file,
                normalization="TTR",
                samples=10000,
                adaptive=False,
                doHistogram=False,
                includeZeros=False,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.DualConditionMethod.__init__(self, short_name, long_name, description, ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        self.includeZeros = includeZeros




    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get selected files
        ctrldata = wxobj.ctrlSelected()
        expdata = wxobj.expSelected()
        if not transit_tools.validate_both_datasets(ctrldata, expdata):
            return None

        #Validate transposon types
        if not transit_tools.validate_filetypes(ctrldata+expdata, transposons):
            return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        samples = int(wxobj.resamplingSampleText.GetValue())
        normalization = wxobj.resamplingNormChoice.GetString(wxobj.resamplingNormChoice.GetCurrentSelection())
        replicates="Sum"
        adaptive = wxobj.resamplingAdaptiveCheckBox.GetValue()
        doHistogram = wxobj.resamplingHistogramCheckBox.GetValue()

        includeZeros = wxobj.resamplingZeroCheckBox.GetValue()

        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = False

        #Get output path
        defaultFileName = "resampling_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")


        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs):

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        ctrldata = args[0].split(",")
        expdata = args[1].split(",")
        annotationPath = args[2]
        output_path = args[3]
        output_file = open(output_path, "w")

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False)
        doHistogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")
        includeZeros = kwargs.get("iz", False)
    
        
        LOESS = kwargs.get("l", False)
        ignoreCodon = True
        NTerminus = float(kwargs.get("iN", 0.00))
        CTerminus = float(kwargs.get("iC", 0.00))

        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)



    def Run(self):

        self.transit_message("Starting resampling Method")
        start_time = time.time()
       

        if self.doHistogram:
            histPath = os.path.join(os.path.dirname(self.output.name), transit_tools.fetch_name(self.output.name)+"_histograms")
            if not os.path.isdir(histPath):
                os.makedirs(histPath)
        else:
            histPath = ""
 



        Kctrl = len(self.ctrldata)
        Kexp = len(self.expdata)
        #Get orf data
        self.transit_message("Getting Data")
        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
        G = tnseq_tools.Genes(self.ctrldata+self.expdata, self.annotation_path, norm=self.normalization, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)


        #Resampling
        data = []
        N = len(G)
        count = 0
        self.progress_range(N)
        for gene in G:
            count+=1
            if gene.k == 0 or gene.n == 0:
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) = (0, 0, 0, 0, 1.00, 1.00, 1.00, [])
            else:

                if not self.includeZeros:
                    ii = numpy.sum(gene.reads,0) > 0
                else:
                    ii = numpy.ones(gene.n) == 1
                
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  stat_tools.resampling(gene.reads[:Kctrl,ii].flatten(), gene.reads[Kctrl:,ii].flatten(), S=self.samples, testFunc=stat_tools.F_sum_diff_flat, adaptive=self.adaptive)


            if self.doHistogram:
                if testlist:
                    n, bins, patches = plt.hist(testlist, normed=1, facecolor='c', alpha=0.75, bins=100)
                else:
                    n, bins, patches = plt.hist([0], normed=1, facecolor='c', alpha=0.75, bins=100)
                plt.xlabel('Delta Sum')
                plt.ylabel('Probability')
                plt.title('%s - Histogram of Delta Sum' % gene.orf)
                plt.axvline(test_obs, color='r', linestyle='dashed', linewidth=3)
                plt.grid(True)
                genePath = os.path.join(histPath, gene.orf +".png")
                plt.savefig(genePath)
                plt.clf()



            data.append([gene.orf, gene.name, gene.desc, gene.n, mean1, mean2, test_obs, log2FC, pval_2tail])
            self.progress_update("resampling", count)
            self.transit_message_inplace("Running Resampling Method... %1.1f%%" % (100.0*count/N))


        #
        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Performing Benjamini-Hochberg Correction")
        data.sort() 
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])
       
 
        self.output.write("#Resampling\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: norm=%s, samples=%s, adaptive=%s, histogram=%s, includeZeros=%s, output=%s\n" % (self.normalization, self.samples, self.adaptive, self.doHistogram, self.includeZeros, self.output))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))
        self.output.write("#Control Data: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Experimental Data: %s\n" % (",".join(self.expdata))) 
        self.output.write("#Annotation path: %s\n" % (self.annotation_path))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail) = row
            self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%1.5f\t%1.5f\n" % (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Resampling")
        self.finish()
        self.transit_message("Finished resampling Method") 


    @classmethod
    def usage_string(self):
        return """python %s resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive resampling. Default: Turned Off.
        -iz             :=  Include rows with zero accross conditions.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = ResamplingMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


