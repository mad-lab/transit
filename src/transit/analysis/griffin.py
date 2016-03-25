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

#method_name = "griffin"


############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.griffinPanel.Hide()

def Show(wxobj):
    wxobj.griffinPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run griffin" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.griffinPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    #wxobj.griffinPanel.SetMinSize( wx.Size( 50,1 ) )
    #wxobj.griffinPanel.SetMaxSize( wx.Size( 250,-1 ) )

    griffinSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.griffinLabel = wx.StaticText( wxobj.griffinPanel, wx.ID_ANY, u"griffin Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.griffinLabel.Wrap( -1 )
    griffinSection.Add( wxobj.griffinLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    griffinSizer1 = wx.BoxSizer( wx.HORIZONTAL )
    #griffinSizer2 = wx.BoxSizer( wx.HORIZONTAL )
    #griffinLabelSizer = wx.BoxSizer( wx.VERTICAL )
    #griffinControlSizer = wx.BoxSizer( wx.VERTICAL )
    
    
    #wxobj.griffinRepLabel = wx.StaticText( wxobj.griffinPanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
    #wxobj.griffinRepLabel.Wrap(-1)
    #griffinLabelSizer.Add(wxobj.griffinRepLabel, 1, wx.ALL, 5)
    #griffinRepChoiceChoices = [ u"Sum", u"Mean", "TTRMean" ]
    #wxobj.griffinRepChoice = wx.Choice( wxobj.griffinPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, griffinRepChoiceChoices, 0 )
    #wxobj.griffinRepChoice.SetSelection( 2 )
    #griffinControlSizer.Add(wxobj.griffinRepChoice, 0, wx.ALL|wx.EXPAND, 5)
    #griffinSizer2.Add(griffinLabelSizer, 1, wx.EXPAND, 5)
    #griffinSizer2.Add(griffinControlSizer, 1, wx.EXPAND, 5)
    #griffinSizer1.Add(griffinSizer2, 1, wx.EXPAND, 5 )


    griffinSection.Add( griffinSizer1, 1, wx.EXPAND, 5 )

    wxobj.griffinButton = wx.Button( wxobj.griffinPanel, wx.ID_ANY, u"Run griffin", wx.DefaultPosition, wx.DefaultSize, 0 )
    griffinSection.Add( wxobj.griffinButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.griffinPanel.SetSizer( griffinSection )
    wxobj.griffinPanel.Layout()
    griffinSection.Fit( wxobj.griffinPanel )

    #Connect events
    wxobj.griffinButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.griffinPanel


def updateProgressBar(wxobj, count):
    wxobj.griffinProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.griffinProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.griffinButton.Enable()





########## CLASS #######################

class Griffin(base.SingleConditionMethod):
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

        base.SingleConditionMethod.__init__(self, "griffin", "griffin Method", "A basic griffin method to show how you could add your own method to TRANSIT.", ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)



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
        replicates = "Sum"
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(all_selected[0])
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

        print "ARGS:", args
        print "KWARGS:", kwargs

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
        
        self.output.write("#griffin\n")
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
        self.output.write("#Orf\tName\tDesc\tk\tn\tr\ts\tt\tExpected Run\tp-value\n")
        
        for (gene, exprun, pval, padj) in results:
            self.output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.1f\t%1.5f\t%1.5f\n" % (gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, gene.s, gene.t, exprun, pval, padj))

        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.transit_message("Finished Griffin Method")

    @classmethod
    def usage_string(self):
        return """python %s griffin <comma-separated .wig files> <annotation .prot_table> <output file>""" % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    print "ARGS:", args
    print "KWARGS:", kwargs

    G = Griffin.fromargs(sys.argv)

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


