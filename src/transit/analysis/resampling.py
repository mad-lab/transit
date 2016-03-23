import sys
import wx
import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

import matplotlib.pyplot as plt

import base
import transit_tools

import tnseq_tools
import norm_tools
import stat_tools



############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.resamplingPanel.Hide()

def Show(wxobj):
    wxobj.resamplingPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run Resampling" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.resamplingPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    #wxobj.resamplingPanel.SetMinSize( wx.Size( 50,1 ) )
    #wxobj.resamplingPanel.SetMaxSize( wx.Size( 250,-1 ) )


    resamplingSizer = wx.BoxSizer( wx.VERTICAL )

    wxobj.resamplingLabel = wx.StaticText( wxobj.resamplingPanel, wx.ID_ANY, u"resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.resamplingLabel.Wrap( -1 )
    resamplingSizer.Add( wxobj.resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    resamplingTopSizer = wx.BoxSizer( wx.HORIZONTAL )

    resamplingTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

    resamplingLabelSizer = wx.BoxSizer( wx.VERTICAL )

    wxobj.resamplingSampleLabel = wx.StaticText( wxobj.resamplingPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.resamplingSampleLabel.Wrap( -1 )
    resamplingLabelSizer.Add( wxobj.resamplingSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.resamplingNormLabel = wx.StaticText( wxobj.resamplingPanel, wx.ID_ANY, u"Normalization", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.resamplingNormLabel.Wrap( -1 )
    resamplingLabelSizer.Add( wxobj.resamplingNormLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    resamplingTopSizer2.Add( resamplingLabelSizer, 1, wx.EXPAND, 5 )

    resamplingControlSizer = wx.BoxSizer( wx.VERTICAL )

    wxobj.resamplingSampleText = wx.TextCtrl( wxobj.resamplingPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
    resamplingControlSizer.Add( wxobj.resamplingSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    resamplingNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
    wxobj.resamplingNormChoice = wx.Choice( wxobj.resamplingPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, resamplingNormChoiceChoices, 0 )
    wxobj.resamplingNormChoice.SetSelection( 0 )
    resamplingControlSizer.Add( wxobj.resamplingNormChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )


    resamplingTopSizer2.Add( resamplingControlSizer, 1, wx.EXPAND, 5 )

    resamplingTopSizer.Add( resamplingTopSizer2, 1, wx.EXPAND, 5 )

    resamplingSizer.Add( resamplingTopSizer, 1, wx.EXPAND, 5 )

    wxobj.resamplingButton = wx.Button( wxobj.resamplingPanel, wx.ID_ANY, u"Run resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
    resamplingSizer.Add( wxobj.resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

 
    wxobj.resamplingPanel.SetSizer( resamplingSizer )
    wxobj.resamplingPanel.Layout()
    resamplingSizer.Fit( wxobj.resamplingPanel )

    #Connect events
    wxobj.resamplingButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.resamplingPanel


def updateProgressBar(wxobj, count):
    wxobj.resamplingProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.resamplingProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.resamplingButton.Enable()





########## CLASS #######################

class Resampling(base.DualConditionMethod):
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
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.DualConditionMethod.__init__(self, "Resampling", "resampling Method", "The permutation test to determing change in read-counts between conditions.", ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        




    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get selected ctrl files
        ctrl_selected = wxobj.ctrlSelected()
        if len(ctrl_selected) ==0:
            wxobj.ShowError("Error: No Control dataset selected.")
            return None

        exp_selected = wxobj.expSelected()
        if len(exp_selected) ==0:
            wxobj.ShowError("Error: No Experimental dataset selected.")
            return None


        #Get Annotation file
        annotationPath = wxobj.annotationFilePicker.GetPath()
        if not annotationPath:
            wxobj.ShowError("Error: No annotation file selected.")
            return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        samples = int(wxobj.resamplingSampleText.GetValue())
        normalization = wxobj.resamplingNormChoice.GetString(wxobj.resamplingNormChoice.GetCurrentSelection())
        replicates="Sum"
        adaptive = False
        doHistogram = False

        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = False

        #Get output path
        defaultFileName = "resampling_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")


        return self(ctrl_selected,
                exp_selected,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs):

        print "RAW:", rawargs
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        print "ARGS:", args
        print "KWARGS:", kwargs

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
        if self.normalization != "none":
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
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  stat_tools.resampling(gene.reads[:Kctrl,:].flatten(), gene.reads[Kctrl:,:].flatten(), S=self.samples, testFunc=stat_tools.F_sum_diff_flat, adaptive=self.adaptive)


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
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata), self.annotation_path, self.output))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Annotation path: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#Orf\tName\tDesc\tSites\tMean A\tMean B\tDelta sum\tlog2FC\tpvalue\tadj. pvalue\n")

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail) = row
            self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%1.5f\t%1.5f\n" % (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
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
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    print "ARGS:", args
    print "KWARGS:", kwargs

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = Resampling.fromargs(sys.argv)

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


