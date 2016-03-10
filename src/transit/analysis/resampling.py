import sys
import wx
import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

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

    resamplingSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.resamplingLabel = wx.StaticText( wxobj.resamplingPanel, wx.ID_ANY, u"resampling Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.resamplingLabel.Wrap( -1 )
    resamplingSection.Add( wxobj.resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    resamplingSizer1 = wx.BoxSizer( wx.HORIZONTAL )
    resamplingSection.Add( resamplingSizer1, 1, wx.EXPAND, 5 )

    wxobj.resamplingButton = wx.Button( wxobj.resamplingPanel, wx.ID_ANY, u"Run resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
    resamplingSection.Add( wxobj.resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.resamplingPanel.SetSizer( resamplingSection )
    wxobj.resamplingPanel.Layout()
    resamplingSection.Fit( wxobj.resamplingPanel )

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
                exp,data,
                annotation_path,
                output_file,
                normalization=None,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.DualConditionMethod.__init__(self, "Resampling", "resampling Method", "The permutation test to determing change in read-counts between conditions.", ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)



    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get selected ctrl files
        ctrl_selected = wxobj.ctrlSelected()
        if len(ctrl_selected) ==0:
            wxobj.ShowError("Error: No Control dataset selected.")
            return None

        exp_selected = wxobj.expSelected()
        if len(wxp_selected) ==0:
            wxobj.ShowError("Error: No Experimental dataset selected.")
            return None


        #Get Annotation file
        annotationPath = wxobj.annotationFilePicker.GetPath()
        if not annotationPath:
            wxobj.ShowError("Error: No annotation file selected.")
            return None


        #Read the parameters from the wxPython widgets
        ctrldata = all_selected
        ignoreCodon = True
        wxobj.resamplingNormChoice.GetString(wxobj.resamplingNormChoice.GetCurrentSelection())
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        normalization = None
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
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromconsole(self): 
        (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

        print "ARGS:", args
        print "KWARGS:", kwargs

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]

        normalization = None
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0

        return self(ctrldata,
                annotationPath,
                outpath,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting resampling Method")
        start_time = time.time()
        
        Kctrl = len(self.ctrldata)
        Kexp = len(self.expdata)
        #Get orf data
        self.transit_message("Getting Data")
        G = tnseq_tools.Genes(self.ctrldata+self.expdata, self.annotation_path, norm=self.normalization, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)

        data = []
        N = len(G)
        count = 0
        self.progress_range(N)
        for gene in G:
            count+=1
            if gene.k == 0 or gene.n == 0:
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail) = (0, 0, 0, 0, 1.00, 1.00, 1.00)
            else:
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail) =  stat_tools.resampling(gene.reads[:Kctrl,:].flatten(), gene.reads[Kctrl:,:].flatten(), S=self.samples)
            #data.append("%s\t%s\t%s\t%s\t%s\t%1.2f\t%1.2f\n" % (gene.orf, gene.name, gene.desc, gene.n, mean1, mean2, test_obs, log2FC, pval_2tail))

            data.append([gene.orf, gene.name, gene.desc, gene.n, mean1, mean2, test_obs, log2FC, pval_2tail])
            self.progress_update("resampling", count)
            self.transit_message_inplace("Running Resampling Method... %1.1f%%" % (100.0*count/N))

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
        self.output.write("#Orf\tName\tDesc\tk\tn\tmean\tnzmean\n")

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail) = row
            self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%1.5f\t%1.5f" % (orf, name, desc, n, mean1, mean2, test_obs, log2FC, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.transit_message("Finished resampling Method") 




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    print "ARGS:", args
    print "KWARGS:", kwargs

    G = resampling.fromconsole()

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


