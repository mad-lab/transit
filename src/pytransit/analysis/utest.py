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

import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "utest"
long_name = "Mann-Whitney U-test "
short_desc = "Mann-Whitney U-test of conditional essentiality between two conditions"
long_desc = """Mann-Whitney U-test for determining conditional essentiality. Based on rank order statistics to identify significant changes in mean read-counts between two conditions."""

transposons = ["himar1", "tn5"]
columns = ["Orf","Name","Desc","Sites","Mean Ctrl","Mean Exp","log2FC", "U-Statistic","p-value","Adj. p-value"]

class UTestAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, UTestMethod, UTestGUI, [UTestFile])



############# FILE ##################

class UTestFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#utest", columns)

    def getHeader(self, path):
        DE=0; poslogfc=0; neglogfc=0;
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if float(tmp[-1]) < 0.05:
                DE +=1
                if float(tmp[-4]) > 0:
                    poslogfc+=1
                else:
                    neglogfc+=1

        text = """Results:
    Conditionally - Essentials: %s
        Less Essential in Experimental datasets: %s
        More Essential in Experimental datasets: %s
            """ % (DE, poslogfc, neglogfc)
        return text


    def getMenus(self):
        menus = []
        menus.append(("Display in Track View", self.displayInTrackView))
        return menus




############# GUI ##################

class UTestGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        utestPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        utestSizer = wx.BoxSizer( wx.VERTICAL )

        utestLabel = wx.StaticText( utestPanel, wx.ID_ANY, u"utest Options", wx.DefaultPosition, (120,-1), 0 )
        utestLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
        utestSizer.Add( utestLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        utestTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        utestTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        utestLabelSizer = wx.BoxSizer( wx.VERTICAL )

        mainSizer1 = wx.BoxSizer( wx.VERTICAL )

        #(, , Sizer) = self.defineChoiceBox(utestPanel, u"", u"", "")
        #mainSizer1.Add(Sizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Norm
        utestNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        (utestNormLabel, self.wxobj.utestNormChoice, normSizer) = self.defineChoiceBox(utestPanel, u"Normalization:", utestNormChoiceChoices, "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ")
        mainSizer1.Add(normSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )


        utestSizer.Add( mainSizer1, 1, wx.EXPAND, 5 )


        # LOESS Check
        (self.wxobj.utestLoessCheck, loessCheckSizer) = self.defineCheckBox(utestPanel, labelText="Correct for Genome Positional Bias", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
        utestSizer.Add( loessCheckSizer, 0, wx.EXPAND, 5 )

        # LOESS Button
        self.wxobj.utestLoessPrev = wx.Button( utestPanel, wx.ID_ANY, u"Preview LOESS fit", wx.DefaultPosition, wx.DefaultSize, 0 )
        utestSizer.Add( self.wxobj.utestLoessPrev, 0, wx.ALL|wx.CENTER, 5 )


        # Zeros Check
        (self.wxobj.utestZeroCheckBox, zeroSizer) = self.defineCheckBox(utestPanel, labelText="Include sites with all zeros", widgetCheck=True, widgetSize=(-1,-1), tooltipText="Includes sites that are empty (zero) accross all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
        utestSizer.Add(zeroSizer, 0, wx.EXPAND, 5 )


        utestButton = wx.Button( utestPanel, wx.ID_ANY, u"Run U-test", wx.DefaultPosition, wx.DefaultSize, 0 )
        utestSizer.Add( utestButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        utestPanel.SetSizer( utestSizer )
        utestPanel.Layout()
        utestSizer.Fit( utestPanel )

        #Connect events
        utestButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )
        self.wxobj.utestLoessPrev.Bind(wx.EVT_BUTTON, self.wxobj.LoessPrevFunc)

        self.panel = utestPanel



########## CLASS #######################

class UTestMethod(base.DualConditionMethod):
    """
    U-test

    """
    def __init__(self,
                ctrldata,
                expdata,
                annotation_path,
                output_file,
                normalization="TTR",
                includeZeros=False,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.DualConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

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
        if not transit_tools.validate_transposons_used(ctrldata+expdata, transposons):
            return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        normalization = wxobj.utestNormChoice.GetString(wxobj.utestNormChoice.GetCurrentSelection())
        replicates= None

        includeZeros = wxobj.utestZeroCheckBox.GetValue()

        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = wxobj.utestLoessCheck.GetValue()

        #Get output path
        defaultFileName = "utest_%s_output" % (normalization)
        if includeZeros: defaultFileName+= "_iz"
        defaultFileName+=".dat"

        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")


        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
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
        includeZeros = kwargs.get("iz", False)
        replicates = None


        LOESS = kwargs.get("l", False)
        ignoreCodon = True
        NTerminus = float(kwargs.get("iN", 0.00))
        CTerminus = float(kwargs.get("iC", 0.00))

        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
                includeZeros,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)



    def Run(self):

        self.transit_message("Starting Mann-Whitney U-test Method")
        start_time = time.time()



        Kctrl = len(self.ctrldata)
        Kexp = len(self.expdata)
        #Get orf data
        self.transit_message("Getting Data")
        (data, position) = transit_tools.get_validated_data(self.ctrldata+self.expdata, wxobj=self.wxobj)

        (K,N) = data.shape


        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata+self.expdata, self.annotation_path)

        if self.LOESS:
            self.transit_message("Performing LOESS Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])


        G = tnseq_tools.Genes(self.ctrldata + self.expdata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data, position=position)


        #u-test
        data = []
        N = len(G)
        count = 0
        self.progress_range(N)
        for gene in G:
            count+=1
            if gene.k == 0 or gene.n == 0:
                (test_obs, mean1, mean2, log2FC, u_stat, pval_2tail) = (0, 0, 0, 0, 0.0, 1.00)
            else:

                if not self.includeZeros:
                    ii = numpy.sum(gene.reads,0) > 0
                else:
                    ii = numpy.ones(gene.n) == 1


                data1 = gene.reads[:Kctrl,ii].flatten()
                data2 = gene.reads[Kctrl:,ii].flatten()
                try:
                    u_stat, pval_2tail = scipy.stats.mannwhitneyu(data1, data2,
                        alternative="two-sided")
                except ValueError as e:
                    u_stat, pval_2tail = 0.0, 1.00

                n1 = len(data1)
                n2 = len(data2)

                mean1 = 0
                if n1 > 0:
                    mean1 = numpy.mean(data1)
                mean2 = 0
                if n2 > 0:
                    mean2 = numpy.mean(data2)

                try:
                    # Only adjust log2FC if one of the means is zero
                    if mean1 > 0 and mean2 > 0:
                        log2FC = math.log((mean2)/(mean1),2)
                    else:
                        log2FC = math.log((mean2+1.0)/(mean1+1.0),2)
                except:
                    log2FC = 0.0


            #["Orf","Name","Desc","Sites","Mean Ctrl","Mean Exp","log2FC", "U-Statistic","p-value","Adj. p-value"]


            data.append([gene.orf, gene.name, gene.desc, gene.n, mean1, mean2, log2FC, u_stat, pval_2tail])

            # Update Progress
            text = "Running Mann-Whitney U-test Method... %1.1f%%" % (100.0*count/N)
            self.progress_update(text, count)


        #
        self.transit_message("") # Printing empty line to flush stdout
        self.transit_message("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])


        self.output.write("#utest\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: norm=%s, includeZeros=%s, output=%s\n" % (self.normalization, self.includeZeros, self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))
        self.output.write("#Control Data: %s\n" % (",".join(self.ctrldata).encode('utf-8')))
        self.output.write("#Experimental Data: %s\n" % (",".join(self.expdata).encode('utf-8')))
        self.output.write("#Annotation path: %s\n" % (self.annotation_path.encode('utf-8')))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, log2FC, u_stat, pval_2tail) = row
            self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%1.5f\t%1.5f\n" % (orf, name, desc, n, mean1, mean2, log2FC, u_stat, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="utest")
        self.finish()
        self.transit_message("Finished Mann-Whitney U-test Method")


    @classmethod
    def usage_string(self):
        return """python %s utest <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>     :=  Normalization method. Default: -n TTR
        -iz             :=  Include rows with zero accross conditions.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        """ % (sys.argv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = UTestMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print ""
    print "Running:"

    G.Run()


