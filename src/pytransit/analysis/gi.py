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
import ntpath
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.analysis import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "gi"
long_name = "Genetic Interactions"
short_desc = "Genetic interactions analysis for change in enrichment"
long_desc = """Method for determining genetic interactions based on changes in enrichment (i.e. delta log fold-change in mean read counts).

NOTE: This method requires 4 groups of datasets. Use the main interface to add datasets for the two strain backgrounds under the first condition. After pressing "Run GI", a subsequent window will allow you to add the datasets under the second condition.
"""

transposons = ["himar1"]
columns = ["Orf","Name","Number of TA Sites","Mean count (Strain A Condition 1)","Mean count (Strain A Condition 2)","Mean count (Strain B Condition 1)","Mean count (Strain B Condition 2)", "Mean logFC (Strain A)", "Mean logFC (Strain B)", "Mean delta logFC","Lower Bound delta logFC","Upper Bound delta logFC", "Prob. of delta-logFC being within ROPE", "Adjusted Probability (BFDR)", "Is HDI outside ROPE?", "Type of Interaction"]


class GIAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, GIMethod, GIGUI, [GIFile])



############# FILE ##################

class GIFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#GI", columns)

    def getHeader(self, path):
        types_to_counts = {}
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[-1] not in types_to_counts: types_to_counts[tmp[-1]] = 0
            types_to_counts[tmp[-1]]+=1


        text = """
Results:
    Aggravating:   %s
    Alleviating:   %s
    Suppressive:   %s
    No Interaction:   %s

""" % (types_to_counts["Aggravating"], types_to_counts["Alleviating"],
        types_to_counts["Suppressive"], types_to_counts["No Interaction"])

        return text


    def getMenus(self):
        menus = []
        menus.append(("Display in Track View", self.displayInTrackView))
        return menus




############# GUI ##################

class GIGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        giPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        giSizer = wx.BoxSizer( wx.VERTICAL )

        giLabel = wx.StaticText( giPanel, wx.ID_ANY, u"GI Options", wx.DefaultPosition, (100,-1), 0 )
        giLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
        giSizer.Add( giLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        giTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        giTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        giLabelSizer = wx.BoxSizer( wx.VERTICAL )

        mainSizer1 = wx.BoxSizer( wx.VERTICAL )


        # Samples
        (giSampleLabel, self.wxobj.giSampleText, sampleSizer) = self.defineTextBox(giPanel, u"Samples:", u"10000", "Number of samples to take when estimating the distributions of means. More samples give more accurate estimates at the cost of computation time.")
        mainSizer1.Add(sampleSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )


        # ROPE
        (giROPELabel, self.wxobj.giROPEText, ROPESizer) = self.defineTextBox(giPanel, u"ROPE:", u"0.5", "Region of Practical Equivalence. Area around 0 (i.e. 0.0 +/- ROPE) that defines changes in enrichment (delta-log2FC) that are NOT of interest. Can be thought of as the area representing a null-hypothesis.")
        mainSizer1.Add(ROPESizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )



        # Norm
        giNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        (giNormLabel, self.wxobj.giNormChoice, normSizer) = self.defineChoiceBox(giPanel, u"Normalization:", giNormChoiceChoices, "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ")
        mainSizer1.Add(normSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )


        giSizer.Add( mainSizer1, 1, wx.EXPAND, 5 )


        # LOESS Check
        (self.wxobj.giLoessCheck, loessCheckSizer) = self.defineCheckBox(giPanel, labelText="Correct for Genome Positional Bias", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
        giSizer.Add( loessCheckSizer, 0, wx.EXPAND, 5 )

        # LOESS Button
        self.wxobj.giLoessPrev = wx.Button( giPanel, wx.ID_ANY, u"Preview LOESS fit", wx.DefaultPosition, wx.DefaultSize, 0 )
        giSizer.Add( self.wxobj.giLoessPrev, 0, wx.ALL|wx.CENTER, 5 )


        # Zeros Check
        (self.wxobj.giZeroCheckBox, zeroSizer) = self.defineCheckBox(giPanel, labelText="Include sites with all zeros", widgetCheck=True, widgetSize=(-1,-1), tooltipText="Includes sites that are empty (zero) accross all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
        giSizer.Add(zeroSizer, 0, wx.EXPAND, 5 )


        giButton = wx.Button( giPanel, wx.ID_ANY, u"Run GI", wx.DefaultPosition, wx.DefaultSize, 0 )
        giSizer.Add( giButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        giPanel.SetSizer( giSizer )
        giPanel.Layout()
        giSizer.Fit( giPanel )

        #Connect events
        giButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )
        self.wxobj.giLoessPrev.Bind(wx.EVT_BUTTON, self.wxobj.LoessPrevFunc)

        self.panel = giPanel





if hasWx:

    class DatasetDialog(wx.Dialog):


        def __init__(self, *args, **kw):

            self.wxobj = args[0]
            self.verbose = self.wxobj.verbose

            from wx.lib.buttons import GenBitmapTextButton
            bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16, 16))

            wx.Dialog.__init__(self, None, title="Dataset Dialog")

            self.ID_DONE = wx.NewId()

            self.SetSize((900, 800))
            self.SetTitle("Please select files for the second condition.")

            mainSizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(mainSizer)

            warningText = """

    The Genetic Interactions method requires a total of four sets of datasets. Typically these are 2 strain backgrounds (e.g. Wildtype and Knockout) each grown under two conditions (e.g. in vitro and in vivo, or rich-media and presence of antibiotic).

    The Control and Experimental datasets added in the main TRANSIT interface are assumed to be the two strain backgrounds grown under the first condition. This interface allows you to add the remaining datasets for the second condition.

            """
            warningStaticBox = wx.StaticText(self, wx.ID_ANY, warningText, (-1,-1), (-1, -1), wx.EXPAND)
            warningStaticBox.Wrap(800)
            mainSizer.Add(warningStaticBox, flag=wx.CENTER, border=5)


            # CONTROL
            ctrlSizerB = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Control Samples - Condition B" ), wx.VERTICAL )

            bSizer2 = wx.BoxSizer( wx.HORIZONTAL )

            self.ctrlRemoveButton = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, (96,-1), 0 )
            bSizer2.Add( self.ctrlRemoveButton, 0, wx.ALL, 5 )

            self.ctrlView = wx.Button( self, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.ctrlView.Hide()

            bSizer2.Add( self.ctrlView, 0, wx.ALL, 5 )

            self.ctrlScatter = wx.Button( self, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.ctrlScatter.Hide()

            bSizer2.Add( self.ctrlScatter, 0, wx.ALL, 5 )

            self.ctrlFilePicker = GenBitmapTextButton(self, 1, bmp, '[Click to add Control Dataset(s)]', size= wx.Size(500, -1))
            bSizer2.Add( self.ctrlFilePicker, 1, wx.ALIGN_CENTER_VERTICAL, 5 )

            ctrlSizerB.Add( bSizer2, 0, wx.EXPAND, 5 )

            self.list_ctrl = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )

            self.list_ctrl.SetMaxSize(wx.Size(940,200))
            ctrlSizerB.Add( self.list_ctrl, 1, wx.ALL|wx.EXPAND, 5 )


            # EXPERIMENTAL

            expSizerB = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Experimental Samples - Condition B" ), wx.VERTICAL )

            bSizer3 = wx.BoxSizer( wx.HORIZONTAL )

            self.expRemoveButton = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, (96,-1), 0 )
            bSizer3.Add( self.expRemoveButton, 0, wx.ALL, 5 )

            self.expView = wx.Button( self, wx.ID_ANY, u"Track View", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.expView.Hide()

            bSizer3.Add( self.expView, 0, wx.ALL, 5 )

            self.expScatter = wx.Button( self, wx.ID_ANY, u"Scatter", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.expScatter.Hide()

            bSizer3.Add( self.expScatter, 0, wx.ALL, 5 )


            self.expFilePicker = GenBitmapTextButton(self, 1, bmp, '[Click to add Experimental Dataset(s)]', size= wx.Size(500, -1))
            bSizer3.Add( self.expFilePicker, 1, wx.ALIGN_CENTER_VERTICAL, 5 )


            expSizerB.Add( bSizer3, 0, wx.EXPAND, 5 )

            self.list_exp = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
            self.list_exp.SetMaxSize(wx.Size(940, 200))
            expSizerB.Add( self.list_exp, 1, wx.ALL|wx.EXPAND, 5 )





            # MAIN



            mainSizer.Add( ctrlSizerB, 1, wx.EXPAND, 5 )
            mainSizer.Add( expSizerB, 1, wx.EXPAND, 5 )


            button_sizer = wx.BoxSizer(wx.HORIZONTAL)
            doneButton = wx.Button(self, self.ID_DONE, label='Done')
            cancelButton = wx.Button(self, wx.ID_CANCEL, label='Cancel')


            button_sizer.Add(doneButton, flag=wx.LEFT, border=5)
            button_sizer.Add(cancelButton, flag=wx.LEFT, border=5)

            self.expFilePicker.Bind( wx.EVT_BUTTON, self.loadExpFileFunc )
            self.ctrlFilePicker.Bind( wx.EVT_BUTTON, self.loadCtrlFileFunc )


            mainSizer.Add(button_sizer,
                flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

            doneButton.Bind(wx.EVT_BUTTON, self.OnClose)
            cancelButton.Bind(wx.EVT_BUTTON, self.OnClose)

            self.ctrlFilePicker.Bind( wx.EVT_BUTTON, self.loadCtrlFileFunc )
            self.expFilePicker.Bind( wx.EVT_BUTTON, self.loadExpFileFunc )
            self.ctrlRemoveButton.Bind( wx.EVT_BUTTON, self.ctrlRemoveFunc )
            self.expRemoveButton.Bind( wx.EVT_BUTTON, self.expRemoveFunc )



            self.index_ctrl = 0
            self.list_ctrl.InsertColumn(0, 'File', width=210)
            self.list_ctrl.InsertColumn(1, 'Total Reads', width=85)
            self.list_ctrl.InsertColumn(2, 'Density', width=85)
            self.list_ctrl.InsertColumn(3, 'Mean Count', width=90)
            self.list_ctrl.InsertColumn(4, 'Max Count', width=85)
            self.list_ctrl.InsertColumn(5, 'Full Path', width=403)

            self.index_exp = 0
            self.list_exp.InsertColumn(0, 'File', width=210)
            self.list_exp.InsertColumn(1, 'Total Reads', width=85)
            self.list_exp.InsertColumn(2, 'Density', width=85)
            self.list_exp.InsertColumn(3, 'Mean Count', width=90)
            self.list_exp.InsertColumn(4, 'Max Count', width=85)
            self.list_exp.InsertColumn(5, 'Full Path',width=403)





        def OnClose(self, event):
            if self.IsModal():
                self.EndModal(event.EventObject.Id)
                self.Close()
            else:
                self.Close()
                self.Destroy()


        def ctrlSelected(self, col=5):
            selected_ctrl = []
            current = -1
            while True:
                next = self.list_ctrl.GetNextSelected(current)
                if next == -1:
                    break
                path = self.list_ctrl.GetItem(next, col).GetText()
                selected_ctrl.append(path)
                current = next
            return selected_ctrl

#

        def loadCtrlFileFunc(self, event):
            try:

                dlg = wx.FileDialog(
                    self, message="Choose a file",
                    defaultDir=self.wxobj.workdir,
                    defaultFile="",
                    wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                    style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR
                    )
                if dlg.ShowModal() == wx.ID_OK:
                    paths = dlg.GetPaths()
                    print("You chose the following Control file(s):")
                    for fullpath in paths:
                        print("\t%s" % fullpath)
                        self.loadCtrlFile(fullpath)
                dlg.Destroy()
            except Exception as e:
                transit_tools.transit_message("Error: %s" % e)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)

#

        def expSelected(self, col=5):
            selected_exp = []
            current = -1
            while True:
                next = self.list_exp.GetNextSelected(current)
                if next == -1:
                    break
                path = self.list_exp.GetItem(next, col).GetText()
                selected_exp.append(path)
                current = next
            return selected_exp


        def loadExpFileFunc(self, event):
            try:

                dlg = wx.FileDialog(
                    self, message="Choose a file",
                    defaultDir=self.wxobj.workdir,
                    defaultFile="",
                    wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                    style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR
                    )
                if dlg.ShowModal() == wx.ID_OK:
                    paths = dlg.GetPaths()
                    print("You chose the following Experimental file(s):")
                    for fullpath in paths:
                        print("\t%s" % fullpath)
                        self.loadExpFile(fullpath)
                dlg.Destroy()
            except Exception as e:
                transit_tools.transit_message("Error: %s" % e)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)

#

        def loadCtrlFile(self, fullpath):
            name = transit_tools.basename(fullpath)
            (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = tnseq_tools.get_wig_stats(fullpath)

            self.list_ctrl.InsertItem(self.index_ctrl, name)
            self.list_ctrl.SetItem(self.index_ctrl, 1, "%1.1f" % (totalrd))
            self.list_ctrl.SetItem(self.index_ctrl, 2, "%2.1f" % (density*100))
            self.list_ctrl.SetItem(self.index_ctrl, 3, "%1.1f" % (meanrd))
            self.list_ctrl.SetItem(self.index_ctrl, 4, "%d" % (maxrd))
            self.list_ctrl.SetItem(self.index_ctrl, 5, "%s" % (fullpath))
            self.list_ctrl.Select(self.index_ctrl)
            self.index_ctrl+=1

#

        def loadExpFile(self, fullpath):
            name = transit_tools.basename(fullpath)
            (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = tnseq_tools.get_wig_stats(fullpath)

            self.list_exp.InsertItem(self.index_exp, name)
            self.list_exp.SetItem(self.index_exp, 1, "%1.1f" % (totalrd))
            self.list_exp.SetItem(self.index_exp, 2, "%2.1f" % (density*100))
            self.list_exp.SetItem(self.index_exp, 3, "%1.1f" % (meanrd))
            self.list_exp.SetItem(self.index_exp, 4, "%d" % (maxrd))
            self.list_exp.SetItem(self.index_exp, 5, "%s" % (fullpath))
            self.list_exp.Select(self.index_exp)
            self.index_exp+=1

#

        def ctrlRemoveFunc(self, event):
            next = self.list_ctrl.GetNextSelected(-1)
            while next != -1:
                if self.verbose:
                    transit_tools.transit_message("Removing control item (%d): %s" % (next, self.list_ctrl.GetItem(next, 0).GetText()))
                self.list_ctrl.DeleteItem(next)
                next = self.list_ctrl.GetNextSelected(-1)
                self.index_ctrl-=1

#

        def expRemoveFunc(self, event):
            next = self.list_exp.GetNextSelected(-1)
            while next != -1:
                if self.verbose:
                    transit_tools.transit_message("Removing experimental item (%d): %s" % (next, self.list_exp.GetItem(next, 0).GetText()))
                self.list_exp.DeleteItem(next)
                next = self.list_exp.GetNextSelected(-1)
                self.index_exp-=1






########## CLASS #######################

class GIMethod(base.QuadConditionMethod):
    """
    GI

    """
    def __init__(self,
                ctrldataA,
                ctrldataB,
                expdataA,
                expdataB,
                annotation_path,
                output_file,
                normalization="TTR",
                samples=10000,
                rope=0.5,
                includeZeros=False,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.QuadConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldataA, ctrldataB, expdataA, expdataB, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.samples = samples
        self.includeZeros = includeZeros
        self.rope = rope
        self.doBFDR = True # TRI
        self.doFWER = True # TRI
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus

    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get selected files
        ctrldataA = wxobj.ctrlSelected()
        expdataA = wxobj.expSelected()
        if not transit_tools.validate_both_datasets(ctrldataA, expdataA):
            return None

        #Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldataA+expdataA, transposons):
            return None



        # Get other datasets:
        dlg = DatasetDialog(wxobj)
        result = dlg.ShowModal()
        if result == dlg.ID_DONE and wxobj:
            status = 1
            ctrldataB = dlg.ctrlSelected()
            expdataB = dlg.expSelected()
            if not transit_tools.validate_both_datasets(ctrldataB, expdataB):
                #dlg.Close()
                #dlg.Destroy()
                return None
            if not transit_tools.validate_transposons_used(ctrldataB+expdataB, transposons):
                dlg.Close()
                dlg.Destroy()
                return None
        else:
            dlg.Close()
            dlg.Destroy()
            return None

        # Close
        dlg.Close()
        dlg.Destroy()

        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        samples = int(wxobj.giSampleText.GetValue())
        rope = float(wxobj.giROPEText.GetValue())
        normalization = wxobj.giNormChoice.GetString(wxobj.giNormChoice.GetCurrentSelection())
        replicates="Sum"

        includeZeros = wxobj.giZeroCheckBox.GetValue()

        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = wxobj.giLoessCheck.GetValue()

        #Get output path
        defaultFileName = "genetic_interactions_output_s%d" % (samples)
        if includeZeros: defaultFileName+= "_iz"
        defaultFileName+=".dat"

        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")


        return self(ctrldataA,
                ctrldataB,
                expdataA,
                expdataB,
                annotationPath,
                output_file,
                normalization,
                samples,
                rope,
                includeZeros,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs):

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        # ctrl-vs-exp = condition 1-vs-2
        # originally, MAD defined order of CL args this way: strA/cond1, strB/cond1, strA/cond2, strB/cond
        #ctrldataA = args[0].split(",")
        #ctrldataB = args[1].split(",")
        #expdataA = args[2].split(",")
        #expdataB = args[3].split(",")

        # TRI changed order of args this way: strA/cond1, strA/cond2, strB/cond1, strB/cond
        ctrldataA = args[0].split(",")
        expdataA = args[1].split(",")
        ctrldataB = args[2].split(",")
        expdataB = args[3].split(",")

        annotationPath = args[4]
        output_path = args[5]
        output_file = open(output_path, "w")

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        rope = float(kwargs.get("-rope", 0.5)) # fixed! changed int to float
        replicates = kwargs.get("r", "Sum")
        includeZeros = kwargs.get("iz", False)


        LOESS = kwargs.get("l", False)
        ignoreCodon = True
        NTerminus = float(kwargs.get("iN", 0.00))
        CTerminus = float(kwargs.get("iC", 0.00))

        return self(ctrldataA,
                ctrldataB,
                expdataA,
                expdataB,
                annotationPath,
                output_file,
                normalization,
                samples,
                rope,
                includeZeros,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)



    def Run(self):


        self.transit_message("Starting Genetic Interactions Method")
        start_time = time.time()
        self.output.write("#GI\n")

        wiglist = self.ctrldataA + self.ctrldataB + self.expdataA + self.expdataB

        Nwig = len(wiglist)
        Na1 = len(self.ctrldataA)
        Nb1 = len(self.ctrldataB)
        Na2 = len(self.expdataA)
        Nb2 = len(self.expdataB)

        # Get data
        self.transit_message("Getting Data")
        (data, position) = transit_tools.get_validated_data(wiglist, wxobj=self.wxobj)

        # Normalize data if specified
        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, wiglist, self.annotation_path)

        # Do LOESS correction if specified
        if self.LOESS:
            self.transit_message("Performing LOESS Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        # Get Gene objects for each condition
        G_A1 = tnseq_tools.Genes([], self.annotation_path, data=data[:Na1], position=position,nterm=self.NTerminus,cterm=self.CTerminus)
        G_B1 = tnseq_tools.Genes([], self.annotation_path, data=data[Na1:(Na1+Nb1)], position=position,nterm=self.NTerminus,cterm=self.CTerminus)
        G_A2 = tnseq_tools.Genes([], self.annotation_path, data=data[(Na1+Nb1):(Na1+Nb1+Na2)], position=position,nterm=self.NTerminus,cterm=self.CTerminus)
        G_B2 = tnseq_tools.Genes([], self.annotation_path, data=data[(Na1+Nb1+Na2):], position=position,nterm=self.NTerminus,cterm=self.CTerminus)

        means_list_a1 = []
        means_list_b1 = []
        means_list_a2 = []
        means_list_b2 = []

        var_list_a1 = []
        var_list_a2 = []
        var_list_b1 = []
        var_list_b2 = []


        # Base priors on empirical observations accross genes.
        for gene in sorted(G_A1):
            if gene.n > 1:
                A1_data = G_A1[gene.orf].reads.flatten()
                B1_data = G_B1[gene.orf].reads.flatten()
                A2_data = G_A2[gene.orf].reads.flatten()
                B2_data = G_B2[gene.orf].reads.flatten()

                means_list_a1.append(numpy.mean(A1_data))
                var_list_a1.append(numpy.var(A1_data))

                means_list_b1.append(numpy.mean(B1_data))
                var_list_b1.append(numpy.var(B1_data))

                means_list_a2.append(numpy.mean(A2_data))
                var_list_a2.append(numpy.var(A2_data))

                means_list_b2.append(numpy.mean(B2_data))
                var_list_b2.append(numpy.var(B2_data))

        # Priors
        mu0_A1 = scipy.stats.trim_mean(means_list_a1, 0.01)
        mu0_B1 = scipy.stats.trim_mean(means_list_b1, 0.01)
        mu0_A2 = scipy.stats.trim_mean(means_list_a2, 0.01)
        mu0_B2 = scipy.stats.trim_mean(means_list_b2, 0.01)

        s20_A1 = scipy.stats.trim_mean(var_list_a1, 0.01)
        s20_B1 = scipy.stats.trim_mean(var_list_b1, 0.01)
        s20_A2 = scipy.stats.trim_mean(var_list_a2, 0.01)
        s20_B2 = scipy.stats.trim_mean(var_list_b2, 0.01)

        k0=1.0
        nu0=1.0
        data = []

        postprob = []
        count = 0
        N = len(G_A1)
        self.progress_range(N)
        # Perform actual analysis
        for gene in G_A1:

            # If there is some data
            if gene.n > 0:
                A1_data = G_A1[gene.orf].reads.flatten()
                B1_data = G_B1[gene.orf].reads.flatten()
                A2_data = G_A2[gene.orf].reads.flatten()
                B2_data = G_B2[gene.orf].reads.flatten()


            #            Time-1   Time-2
            #
            #  Strain-A     A       C
            #
            #  Strain-B     B       D

                try:
                    muA1_post, varA1_post = stat_tools.sample_trunc_norm_post(A1_data, self.samples,
                        mu0_A1, s20_A1, k0, nu0)
                    muB1_post, varB1_post = stat_tools.sample_trunc_norm_post(B1_data, self.samples,
                        mu0_B1, s20_B1, k0, nu0)
                    muA2_post, varA2_post = stat_tools.sample_trunc_norm_post(A2_data, self.samples,
                        mu0_A2, s20_A2, k0, nu0)
                    muB2_post, varB2_post = stat_tools.sample_trunc_norm_post(B2_data, self.samples,
                        mu0_B2, s20_B2, k0, nu0)

                except Exception as e:
                    muA1_post = varA1_post = numpy.ones(self.samples)
                    muB1_post = varB1_post = numpy.ones(self.samples)
                    muA2_post = varA2_post = numpy.ones(self.samples)
                    muB2_post = varB2_post = numpy.ones(self.samples)

                logFC_A_post = numpy.log2(muA2_post/muA1_post)
                logFC_B_post = numpy.log2(muB2_post/muB1_post)
                delta_logFC_post = logFC_B_post - logFC_A_post

                alpha = 0.05

                # Get Bounds of the HDI
                l_logFC_A, u_logFC_A = stat_tools.HDI_from_MCMC(logFC_A_post, 1-alpha)

                l_logFC_B, u_logFC_B = stat_tools.HDI_from_MCMC(logFC_B_post, 1-alpha)

                l_delta_logFC, u_delta_logFC = stat_tools.HDI_from_MCMC(delta_logFC_post, 1-alpha)


                mean_logFC_A = numpy.mean(logFC_A_post)
                mean_logFC_B = numpy.mean(logFC_B_post)
                mean_delta_logFC = numpy.mean(delta_logFC_post)

                # Is HDI significantly different than ROPE?
                not_HDI_overlap_bit = l_delta_logFC > self.rope or u_delta_logFC < -self.rope

                # Probability of posterior overlaping with ROPE
                probROPE = numpy.mean(numpy.logical_and(delta_logFC_post>=0.0-self.rope,  delta_logFC_post<=0.0+self.rope))

            # If there is no data, assume empty defaults
            else:
                A1_data = [0,0]
                B1_data = [0,0]
                A2_data = [0,0]
                B2_data = [0,0]
                muA1_post = varA1_post = numpy.ones(self.samples)
                muB1_post = varB1_post = numpy.ones(self.samples)
                muA2_post = varA2_post = numpy.ones(self.samples)
                muB2_post = varB2_post = numpy.ones(self.samples)
                logFC_A_post = numpy.log2(muA2_post/muA1_post)
                logFC_B_post = numpy.log2(muB2_post/muB1_post)
                delta_logFC_post = logFC_B_post - logFC_A_post

                mean_logFC_A = 0
                mean_logFC_B = 0
                mean_delta_logFC = 0
                l_logFC_A = 0
                u_logFC_A = 0
                l_logFC_B = 0
                u_logFC_B = 0
                l_delta_logFC = 0
                u_delta_logFC = 0
                probROPE = 1.0


            if numpy.isnan(l_logFC_A):
                l_logFC_A = -10
                u_logFC_A = 10
            if numpy.isnan(l_logFC_B):
                l_logFC_B = -10
                u_logFC_B = 10
            if numpy.isnan(l_delta_logFC):
                l_delta_logFC = -10
                u_delta_logFC = 10


            postprob.append(probROPE)
            data.append((gene.orf, gene.name, gene.n, numpy.mean(muA1_post), numpy.mean(muA2_post), numpy.mean(muB1_post), numpy.mean(muB2_post), mean_logFC_A, mean_logFC_B, mean_delta_logFC, l_delta_logFC, u_delta_logFC, probROPE, not_HDI_overlap_bit))


            text = "Running GI Method... %2.0f%%" % (100.0*(count+1)/N)
            self.progress_update(text, count)
            self.transit_message_inplace("Running Export Method... %1.1f%%" % (100.0*count/(N-1)))
            count+=1

        data.sort(key=lambda x: x[-2])

        if self.doBFDR or not self.doFWER:
            postprob = numpy.array(postprob)
            postprob.sort()
            bfdr = numpy.cumsum(postprob)/numpy.arange(1, len(postprob)+1)
            adjusted_prob = bfdr
            adjusted_label = "BFDR"
        elif doFWER:
            fwer = FWER_Bayes(postprob)
            fwer.sort()
            adjusted_prob = fwer
            adjusted_label = "FWER"

        # If not using adjustment for classification, sort correctly
        if not self.doBFDR and not self.doFWER:
            sorted_index = numpy.argsort([d[-1] for d in data])[::-1][:len(data)]
            adjusted_prob = [adjusted_prob[ii] for ii in sorted_index]
            data = [data[ii] for ii in sorted_index]



        # Print(output)
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: norm=%s, samples=%s, includeZeros=%s, output=%s\n" % (self.normalization, self.samples, self.includeZeros, self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))


        self.output.write("#Control Data-A: %s\n" % (",".join(self.ctrldataA).encode('utf-8')))
        self.output.write("#Control Data-B: %s\n" % (",".join(self.ctrldataB).encode('utf-8')))
        self.output.write("#Experimental Data-A: %s\n" % (",".join(self.expdataA).encode('utf-8')))
        self.output.write("#Experimental Data-B: %s\n" % (",".join(self.expdataB).encode('utf-8')))
        self.output.write("#Annotation path: %s\n" % (self.annotation_path.encode('utf-8')))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))


        if self.doBFDR or self.doFWER:
            self.output.write("# Significant interactions are those whose adjusted probability of the delta-logFC falling within ROPE is < 0.05 (Adjusted using %s)\n" % (adjusted_label))
        else:
            self.output.write("# Significant interactions are those genes whose delta-logFC HDI does not overlap the ROPE\n")
        self.output.write("#\n")

        # Write column names
        self.output.write("#ORF\tName\tNumber of TA Sites\tMean count (Strain A Condition 1)\tMean count (Strain A Condition 2)\tMean count (Strain B Condition 1)\tMean count (Strain B Condition 2)\tMean logFC (Strain A)\tMean logFC (Strain B) \tMean delta logFC\tLower Bound delta logFC\tUpper Bound delta logFC\tProb. of delta-logFC being within ROPE\tAdjusted Probability (%s)\tIs HDI outside ROPE?\tType of Interaction\n" % adjusted_label)

        # Write gene results
        for i,row in enumerate(data):
        #1   2    3        4                5              6               7                8            9            10              11             12            13         14
            orf, name, n, mean_muA1_post, mean_muA2_post, mean_muB1_post, mean_muB2_post, mean_logFC_A, mean_logFC_B, mean_delta_logFC, l_delta_logFC, u_delta_logFC, probROPE, not_HDI_overlap_bit = row
            type_of_interaction = "No Interaction"
            if ((self.doBFDR or self.doFWER) and adjusted_prob[i] < 0.05):
                type_of_interaction = self.classify_interaction(mean_delta_logFC, mean_logFC_B, mean_logFC_A)
            elif not (self.doBFDR or self.doFWER) and not_HDI_overlap_bit:
                type_of_interaction = self.classify_interaction(mean_delta_logFC, mean_logFC_B, mean_logFC_A)

            new_row = tuple(list(row[:-1])+[adjusted_prob[i], not_HDI_overlap_bit, type_of_interaction])
            self.output.write("%s\t%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.8f\t%1.8f\t%s\t%s\n" % new_row)


        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="GI")
        self.finish()
        self.transit_message("Finished Genetic Interactions Method")


    @staticmethod
    def classify_interaction(delta_logFC, logFC_KO, logFC_WT):
        if delta_logFC < 0:
            return "Aggravating"
        elif delta_logFC >= 0 and abs(logFC_KO) < abs(logFC_WT):
            return "Alleviating"
        elif delta_logFC >= 0 and abs(logFC_KO) > abs(logFC_WT):
            return "Suppressive"
        else:
            return "N/A"


    @classmethod
    def usage_string(self):
        #return """python3 %s GI <comma-separated .wig control files condition A> <comma-separated .wig control files condition B> <comma-separated .wig experimental files condition A> <comma-separated .wig experimental files condition B> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        return """python3 %s GI <wigs_for_strA_cond1> <wigs_for_strA_cond2> <wigs_for_strB_cond1> <wigs_for_strB_cond2> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

        GI performs a comparison among 4 groups of datasets, strain A and B assessed in conditions 1 and 2 (e.g. control vs treatment).
        It looks for interactions where the response to the treatment (i.e. effect on insertion counts) depends on the strain (output variable: delta_LFC).
        Provide replicates in each group as a comma-separated list of wig files.
        Significant interactions are those with "HDI outside ROPE?"=TRUE, and all genes are sorted by significance using BFDR."

        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        --rope <float>  :=  Region of Practical Equivalence. Area around 0 (i.e. 0 +/- ROPE) that is NOT of interest. Can be thought of similar to the area of the null-hypothesis. Default: --rope 0.5
        -n <string>     :=  Normalization method. Default: -n TTR
        -iz             :=  Include rows with zero accross conditions.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring at given percentage (as integer) of the C terminus. Default: -iC 0
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = GIMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()


