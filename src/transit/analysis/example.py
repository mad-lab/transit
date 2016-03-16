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

#method_name = "example"


############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.examplePanel.Hide()

def Show(wxobj):
    wxobj.examplePanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run Example" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.examplePanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    #wxobj.examplePanel.SetMinSize( wx.Size( 50,1 ) )
    #wxobj.examplePanel.SetMaxSize( wx.Size( 250,-1 ) )

    exampleSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.exampleLabel = wx.StaticText( wxobj.examplePanel, wx.ID_ANY, u"Example Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.exampleLabel.Wrap( -1 )
    exampleSection.Add( wxobj.exampleLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    exampleSizer1 = wx.BoxSizer( wx.HORIZONTAL )
    exampleSection.Add( exampleSizer1, 1, wx.EXPAND, 5 )

    wxobj.exampleButton = wx.Button( wxobj.examplePanel, wx.ID_ANY, u"Run Example", wx.DefaultPosition, wx.DefaultSize, 0 )
    exampleSection.Add( wxobj.exampleButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.examplePanel.SetSizer( exampleSection )
    wxobj.examplePanel.Layout()
    exampleSection.Fit( wxobj.examplePanel )

    #Connect events
    wxobj.exampleButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.examplePanel


def updateProgressBar(wxobj, count):
    wxobj.exampleProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.exampleProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.exampleButton.Enable()





########## CLASS #######################

class Example(base.SingleConditionMethod):
    """   
    Example
 
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

        base.SingleConditionMethod.__init__(self, "Example", "Example Method", "A basic Example method to show how you could add your own method to TRANSIT.", ctrldata, annotation_path, output_file, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)



    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get selected files
        all_selected = wxobj.ctrlSelected()
        if len(all_selected) ==0:
            wxobj.ShowError("Error: No dataset selected.")
            return None

        #Get Annotation file
        annotationPath = wxobj.annotationFilePicker.GetPath()
        if not annotationPath:
            wxobj.ShowError("Error: No annotation file selected.")
            return None


        #Read the parameters from the wxPython widgets
        ctrldata = all_selected
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(all_selected[0])
        defaultFileName = "example_output.dat"
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

        print "ARGS:", args
        print "KWARGS:", kwargs

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        normalization = None
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

        self.transit_message("Starting Example Method")
        start_time = time.time()
        
        #Get orf data
        self.transit_message("Getting Data")
        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)

        data = []
        N = len(G)
        count = 0
        self.progress_range(N)
        for gene in G:
            count+=1
            mean = numpy.mean(gene.reads)
            if gene.k == 0:
                nzmean = 0.0
            else:
                nzmean = numpy.sum(gene.reads)/float(gene.k)

            data.append("%s\t%s\t%s\t%s\t%s\t%1.2f\t%1.2f\n" % (gene.orf, gene.name, gene.desc, gene.k, gene.n, mean, nzmean))

            
            self.progress_update("gumbel", count)
            self.transit_message_inplace("Running Example Method... %1.1f%%" % (100.0*count/N))
            
        
        self.output.write("#Example\n")
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

        data.sort()
        for line in data:
            self.output.write(line)
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.transit_message("Finished Example Method") 




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    print "ARGS:", args
    print "KWARGS:", kwargs

    G = Example.fromargs(sys.argv)

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


