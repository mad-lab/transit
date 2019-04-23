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

short_name = "resampling"
long_name = "Resampling (Permutation test)"
short_desc = "Resampling test of conditional essentiality between two conditions"
long_desc = """Method for determining conditional essentiality based on resampling (i.e. permutation test). Identifies significant changes in mean read-counts for each gene after normalization."""

transposons = ["himar1", "tn5"]
columns = ["Orf","Name","Desc","Sites","Mean Ctrl","Mean Exp","log2FC", "Sum Ctrl", "Sum Exp", "Delta Mean","p-value","Adj. p-value"]

class ResamplingAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, ResamplingMethod, ResamplingGUI, [ResamplingFile])



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
        Less Essential in Experimental datasets: %s
        More Essential in Experimental datasets: %s
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
                transit_tools.ShowError(MSG="Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
                print "Error Displaying File. Histogram image does not exist."




############# GUI ##################

class ResamplingGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        resamplingPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        resamplingSizer = wx.BoxSizer( wx.VERTICAL )

        resamplingLabel = wx.StaticText( resamplingPanel, wx.ID_ANY, u"resampling Options", wx.DefaultPosition, (160,-1), 0 )
        resamplingLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
        resamplingSizer.Add( resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        resamplingTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        resamplingTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        resamplingLabelSizer = wx.BoxSizer( wx.VERTICAL )

        mainSizer1 = wx.BoxSizer( wx.VERTICAL )

        #(, , Sizer) = self.defineChoiceBox(resamplingPanel, u"", u"", "")
        #mainSizer1.Add(Sizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Samples
        (resamplingSampleLabel, self.wxobj.resamplingSampleText, sampleSizer) = self.defineTextBox(resamplingPanel, u"Samples:", u"10000", "Number of samples to take when estimating the resampling histogram. More samples give more accurate estimates of the p-values at the cost of computation time.")
        mainSizer1.Add(sampleSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Pseudocount
        (resamplingPseudocountLabel, self.wxobj.resamplingPseudocountText, pseudoSizer) = self.defineTextBox(resamplingPanel, u"Pseudocount:", u"0.0", "Adds pseudo-counts to the each data-point. Useful to dampen the effects of small counts which may lead to deceptively high log-FC.")
        mainSizer1.Add(pseudoSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Norm
        resamplingNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        (resamplingNormLabel, self.wxobj.resamplingNormChoice, normSizer) = self.defineChoiceBox(resamplingPanel, u"Normalization: ", resamplingNormChoiceChoices, "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ")
        mainSizer1.Add(normSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )




        resamplingSizer.Add( mainSizer1, 1, wx.EXPAND, 5 )






        # LOESS Check
        (self.wxobj.resamplingLoessCheck, loessCheckSizer) = self.defineCheckBox(resamplingPanel, labelText="Correct for Genome Positional Bias", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
        resamplingSizer.Add( loessCheckSizer, 0, wx.EXPAND, 5 )

        # LOESS Button
        self.wxobj.resamplingLoessPrev = wx.Button( resamplingPanel, wx.ID_ANY, u"Preview LOESS fit", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSizer.Add( self.wxobj.resamplingLoessPrev, 0, wx.ALL|wx.CENTER, 5 )

        # Adaptive Check
        (self.wxobj.resamplingAdaptiveCheckBox, adaptiveSizer) = self.defineCheckBox(resamplingPanel, labelText="Adaptive Resampling (Faster)", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Dynamically stops permutations early if it is unlikely the ORF will be significant given the results so far. Improves performance, though p-value calculations for genes that are not differentially essential will be less accurate.")
        resamplingSizer.Add( adaptiveSizer, 0, wx.EXPAND, 5 )

        # Histogram Check
        (self.wxobj.resamplingHistogramCheckBox, histSizer) = self.defineCheckBox(resamplingPanel, labelText="Generate Resampling Histograms", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Creates .png images with the resampling histogram for each of the ORFs. Histogram images are created in a folder with the same name as the output file.")
        resamplingSizer.Add(histSizer, 0, wx.EXPAND, 5 )


        # Zeros Check
        (self.wxobj.resamplingZeroCheckBox, zeroSizer) = self.defineCheckBox(resamplingPanel, labelText="Include sites with all zeros", widgetCheck=True, widgetSize=(-1,-1), tooltipText="Includes sites that are empty (zero) accross all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
        resamplingSizer.Add(zeroSizer, 0, wx.EXPAND, 5 )


        resamplingButton = wx.Button( resamplingPanel, wx.ID_ANY, u"Run resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSizer.Add( resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        resamplingPanel.SetSizer( resamplingSizer )
        resamplingPanel.Layout()
        resamplingSizer.Fit( resamplingPanel )

        #Connect events
        resamplingButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )
        self.wxobj.resamplingLoessPrev.Bind(wx.EVT_BUTTON, self.wxobj.LoessPrevFunc)

        self.panel = resamplingPanel

    def GlobalEnable(self):
        self.wxobj.ctrlLibText.Enable()
        self.wxobj.expLibText.Enable()

    def GlobalDisable(self):
        self.wxobj.ctrlLibText.Disable()
        self.wxobj.expLibText.Disable()



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
                pseudocount=0.0,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0,
                ctrl_lib_str="",
                exp_lib_str="",
                wxobj=None, Z = False, diffStrains = False, annotation_path_exp = "", combinedWigParams = None):

        base.DualConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.Z = Z
        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        self.includeZeros = includeZeros
        self.pseudocount = pseudocount
        self.ctrl_lib_str = ctrl_lib_str
        self.exp_lib_str = exp_lib_str
        self.diffStrains = diffStrains
        self.annotation_path_exp = annotation_path_exp if diffStrains else annotation_path
        self.combinedWigParams = combinedWigParams

    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annot_paths = wxobj.annotation.split(",")
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
            diffStrains = True

        if not transit_tools.validate_annotation(annotationPath):
            return None

        if annotationPathExp and not transit_tools.validate_annotation(annotationPathExp):
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
        samples = int(wxobj.resamplingSampleText.GetValue())
        normalization = wxobj.resamplingNormChoice.GetString(wxobj.resamplingNormChoice.GetCurrentSelection())
        replicates="Sum"
        adaptive = wxobj.resamplingAdaptiveCheckBox.GetValue()
        doHistogram = wxobj.resamplingHistogramCheckBox.GetValue()

        includeZeros = wxobj.resamplingZeroCheckBox.GetValue()
        pseudocount = float(wxobj.resamplingPseudocountText.GetValue())
        LOESS = wxobj.resamplingLoessCheck.GetValue()

        # Global Parameters
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        ctrl_lib_str = wxobj.ctrlLibText.GetValue()
        exp_lib_str = wxobj.expLibText.GetValue()


        #Get output path
        defaultFileName = "resampling_output_s%d_pc%1.2f" % (samples, pseudocount)
        if adaptive: defaultFileName+= "_adaptive"
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
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                pseudocount,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus,
                ctrl_lib_str,
                exp_lib_str, wxobj, Z = False, diffStrains = diffStrains, annotation_path_exp = annotationPathExp)

    @classmethod
    def fromargs(self, rawargs):

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        isCombinedWig = True if kwargs.get('c', False) else False
        combinedWigParams = None
        if isCombinedWig:
            if (len(args) != 5):
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string())
                sys.exit(0)
            combinedWigParams = {
                "combined_wig": kwargs.get('c'),
                "samples_metadata": args[0],
                "conditions": [args[1].lower(), args[2].lower()]
            }
            annot_paths = args[3].split(",")
            ctrldata = ""
            expdata = ""
            output_path = args[4]
        else:
            if (len(args) != 4):
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string())
                sys.exit(0)
            ctrldata = args[0].split(",")
            expdata = args[1].split(",")
            annot_paths = args[2].split(",")
            output_path = args[3]
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
            diffStrains = True
        if (diffStrains and isCombinedWig):
            print("Error: Cannot have combined wig and different annotation files.")
            sys.exit(0)

        output_file = open(output_path, "w")

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False)
        doHistogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")
        excludeZeros = kwargs.get("ez", False)
        includeZeros = not excludeZeros
        pseudocount = float(kwargs.get("pc", 0.00))

        Z = True if "Z" in kwargs else False

        LOESS = kwargs.get("l", False)
        ignoreCodon = True

        NTerminus = float(kwargs.get("iN", 0.00))
        CTerminus = float(kwargs.get("iC", 0.00))
        ctrl_lib_str = kwargs.get("-ctrl_lib", "")
        exp_lib_str = kwargs.get("-exp_lib", "")

        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                pseudocount,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus,
                ctrl_lib_str,
                exp_lib_str, Z = Z, diffStrains = diffStrains, annotation_path_exp = annotationPathExp, combinedWigParams = combinedWigParams)

    def preprocess_data(self, position, data):
        (K,N) = data.shape

        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata+self.expdata, self.annotation_path)

        if self.LOESS:
            self.transit_message("Performing LOESS Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        return data

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, None) for f in filenamesInCombWig]

    def filter_wigs_by_conditions(self, data, conditions, included_conditions):
        """
            Filters conditions from wig to ctrl, exp conditions only
            ([[Wigdata]], [ConditionCtrl, ConditionExp]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered = [], []
        if len(included_conditions) != 2:
            self.transit_error("Only 2 conditions expected", included_conditions)
            sys.exit(0)
        for i, c in enumerate(conditions):
          if (c.lower() in included_conditions):
            d_filtered.append(data[i])
            cond_filtered.append(conditions[i])

        return (numpy.array(d_filtered), numpy.array(cond_filtered))

    def Run(self):

        #if not self.wxobj:
        #    # Force matplotlib to use good backend for png.
        #    import matplotlib.pyplot as plt
        #elif "matplotlib.pyplot" not in sys.modules:
        try:
            import matplotlib.pyplot as plt
        except:
            print "Error: cannot do histograms"
            self.doHistogram = False


        self.transit_message("Starting resampling Method")
        start_time = time.time()

        histPath = ""
        if self.doHistogram:
            histPath = os.path.join(os.path.dirname(self.output.name), transit_tools.fetch_name(self.output.name)+"_histograms")
            if not os.path.isdir(histPath):
                os.makedirs(histPath)

        #Get orf data
        self.transit_message("Getting Data")
        if self.diffStrains:
            self.transit_message("Multiple annotation files found")
            self.transit_message("Mapping ctrl data to {0}, exp data to {1}".format(self.annotation_path, self.annotation_path_exp))

        if self.combinedWigParams:
            (position, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combinedWigParams['combined_wig'])
            conditionsByFile, _, _, _ = tnseq_tools.read_samples_metadata(self.combinedWigParams['samples_metadata'])
            conditions = self.wigs_to_conditions(conditionsByFile, filenamesInCombWig)
            data, conditions = self.filter_wigs_by_conditions(data, conditions, self.combinedWigParams['conditions'])
            data_ctrl = numpy.array([d for i, d in enumerate(data) if conditions[i].lower() == self.combinedWigParams['conditions'][0]])
            data_exp = numpy.array([d for i, d in enumerate(data) if conditions[i].lower() == self.combinedWigParams['conditions'][1]])
            position_ctrl, position_exp = position, position
        else:
            (data_ctrl, position_ctrl) = transit_tools.get_validated_data(self.ctrldata, wxobj=self.wxobj)
            (data_exp, position_exp) = transit_tools.get_validated_data(self.expdata, wxobj=self.wxobj)
        (K_ctrl, N_ctrl) = data_ctrl.shape
        (K_exp, N_exp) = data_exp.shape

        if not self.diffStrains and (N_ctrl != N_exp):
            self.transit_error("Error: Ctrl and Exp wig files don't have the same number of sites.")
            self.transit_error("Make sure all .wig files come from the same strain.")
            return
        # (data, position) = transit_tools.get_validated_data(self.ctrldata+self.expdata, wxobj=self.wxobj)

        self.transit_message("Preprocessing Ctrl data...")
        data_ctrl = self.preprocess_data(position_ctrl, data_ctrl)

        self.transit_message("Preprocessing Exp data...")
        data_exp = self.preprocess_data(position_exp, data_exp)

        G_ctrl = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data_ctrl, position=position_ctrl)
        G_exp = tnseq_tools.Genes(self.expdata, self.annotation_path_exp, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data_exp, position=position_exp)

        doLibraryResampling = False
        # If library string not empty
        if self.ctrl_lib_str or self.exp_lib_str:
            letters_ctrl = set(self.ctrl_lib_str)
            letters_exp = set(self.exp_lib_str)

            # Check if using exactly 1 letters; i.e. no different libraries
            if len(letters_ctrl) == 1 and letters_exp==1:
                pass
            # If using more than one letter, then check no differences in set
            else:
                lib_diff = letters_ctrl ^ letters_exp
                # Check that their differences
                if not lib_diff:
                    doLibraryResampling = True
                else:
                    transit_tools.transit_error("Error: Library Strings (Ctrl = %s, Exp = %s) do not use the same letters. Make sure every letter / library is represented in both Control and Experimental Conditions. Proceeding with resampling assuming all datasets belong to the same library." % (self.ctrl_lib_str, self.exp_lib_str))
                    self.ctrl_lib_str = ""
                    self.exp_lib_str = ""

        (data, qval) = self.run_resampling(G_ctrl, G_exp, doLibraryResampling, histPath)
        self.write_output(data, qval, start_time)

        self.finish()
        self.transit_message("Finished resampling Method")

    def write_output(self, data, qval, start_time):

        self.output.write("#Resampling\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: norm=%s, samples=%s, pseudocounts=%1.2f, adaptive=%s, histogram=%s, includeZeros=%s, output=%s\n" % (self.normalization, self.samples, self.pseudocount, self.adaptive, self.doHistogram, self.includeZeros, self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))
        self.output.write("#Control Data: %s\n" % (",".join(self.ctrldata).encode('utf-8')))
        self.output.write("#Experimental Data: %s\n" % (",".join(self.expdata).encode('utf-8')))
        self.output.write("#Annotation path: %s %s\n" % (self.annotation_path.encode('utf-8'), self.annotation_path_exp.encode('utf-8') if self.diffStrains else ''))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        #Z = True # include Z-score column in resampling output?
        global columns # consider redefining columns above (for GUI)
        if self.Z==True: columns = ["Orf","Name","Desc","Sites","Mean Ctrl","Mean Exp","log2FC", "Sum Ctrl", "Sum Exp", "Delta Mean","p-value","Z-score","Adj. p-value"]
        self.output.write("#%s\n" % "\t".join(columns))

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, sum1, sum2, test_obs, log2FC, pval_2tail) = row
            if self.Z==True:
              p = pval_2tail/2 # convert from 2-sided back to 1-sided
              if p==0: p = 1e-5 # or 1 level deeper the num of iterations of resampling, which is 1e-4=1/10000, by default
              if p==1: p = 1-1e-5
              z = scipy.stats.norm.ppf(p)
              if log2FC>0: z *= -1
              self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%0.2f\t%1.5f\n" % (orf, name, desc, n, mean1, mean2, log2FC, sum1, sum2, test_obs, pval_2tail, z, qval[i]))
            else: self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.1f\t%1.2f\t%1.1f\t%1.5f\t%1.5f\n" % (orf, name, desc, n, mean1, mean2, log2FC, sum1, sum2, test_obs, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Resampling")

    def run_resampling(self, G_ctrl, G_exp = None, doLibraryResampling = False, histPath = ""):
        data = []
        N = len(G_ctrl)
        count = 0
        self.progress_range(N)

        for gene in G_ctrl:
            if gene.orf not in G_exp:
                if self.diffStrains:
                    continue
                else:
                    self.transit_error("Error: Gene in ctrl data not present in exp data")
                    self.transit_error("Make sure all .wig files come from the same strain.")
                    return ([], [])

            gene_exp = G_exp[gene.orf]
            count+=1

            if not self.diffStrains and gene.n != gene_exp.n:
                self.transit_error("Error: No. of TA sites in Exp and Ctrl data are different")
                self.transit_error("Make sure all .wig files come from the same strain.")
                return ([], [])

            if (gene.k == 0 and gene_exp.k == 0) or gene.n == 0 or gene_exp.n == 0:
                (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist, data1, data2) = (0, 0, 0, 0, 1.00, 1.00, 1.00, [], [0], [0])
            else:
                if not self.includeZeros:
                    ii_ctrl = numpy.sum(gene.reads,0) > 0
                    ii_exp = numpy.sum(gene_exp.reads,0) > 0
                else:
                    ii_ctrl = numpy.ones(gene.n) == 1
                    ii_exp = numpy.ones(gene_exp.n) == 1

                data1 = gene.reads[:,ii_ctrl].flatten() + self.pseudocount
                data2 = gene_exp.reads[:,ii_exp].flatten() + self.pseudocount

                if doLibraryResampling:
                    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  stat_tools.resampling(data1, data2, S=self.samples, testFunc=stat_tools.F_mean_diff_dict, permFunc=stat_tools.F_shuffle_dict_libraries, adaptive=self.adaptive, lib_str1=self.ctrl_lib_str, lib_str2=self.exp_lib_str)
                else:
                    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  stat_tools.resampling(data1, data2, S=self.samples, testFunc=stat_tools.F_mean_diff_flat, permFunc=stat_tools.F_shuffle_flat, adaptive=self.adaptive, lib_str1=self.ctrl_lib_str, lib_str2=self.exp_lib_str)


            if self.doHistogram:
                import matplotlib.pyplot as plt
                if testlist:
                    n, bins, patches = plt.hist(testlist, density=1, facecolor='c', alpha=0.75, bins=100)
                else:
                    n, bins, patches = plt.hist([0,0], density=1, facecolor='c', alpha=0.75, bins=100)
                plt.xlabel('Delta Mean')
                plt.ylabel('Probability')
                plt.title('%s - Histogram of Delta Mean' % gene.orf)
                plt.axvline(test_obs, color='r', linestyle='dashed', linewidth=3)
                plt.grid(True)
                genePath = os.path.join(histPath, gene.orf +".png")
                if not os.path.exists(histPath):
                    os.makedirs(histPath)
                plt.savefig(genePath)
                plt.clf()


            sum1 = numpy.sum(data1)
            sum2 = numpy.sum(data2)
            data.append([gene.orf, gene.name, gene.desc, gene.n, mean1, mean2, sum1, sum2, test_obs, log2FC, pval_2tail])

            # Update progress
            text = "Running Resampling Method... %5.1f%%" % (100.0*count/N)
            self.progress_update(text, count)


        #
        self.transit_message("") # Printing empty line to flush stdout
        self.transit_message("Performing Benjamini-Hochberg Correction")
        data.sort()
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])

        return (data, qval)

    @classmethod
    def usage_string(self):
        return """
        python %s resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
        ---
        OR
        ---
        python %s resampling -c <combined wig file> <samples_metadata file> <ctrl condition name> <exp condition name> <annotation .prot_table> <output file> [Optional Arguments]
        NB: The ctrl and exp condition names should match Condition names in samples_metadata file.

        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive resampling. Default: Turned Off.
        -ez             :=  Exclude rows with zero accross conditions. Default: Turned off
                            (i.e. include rows with zeros).
        -pc             :=  Pseudocounts to be added at each site.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias.
                            Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring within given percentage of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring within given percentage of the C terminus. Default: -iC 0.0
        --ctrl_lib      :=  String of letters representing library of control files in order
                            e.g. 'AABB'. Default empty. Letters used must also be used in --exp_lib
                            If non-empty, resampling will limit permutations to within-libraries.

        --exp_lib       :=  String of letters representing library of experimental files in order
                            e.g. 'ABAB'. Default empty. Letters used must also be used in --ctrl_lib
                            If non-empty, resampling will limit permutations to within-libraries.

        """ % (sys.argv[0], sys.argv[0])

if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = ResamplingMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print ""
    print "Running:"

    G.Run()


