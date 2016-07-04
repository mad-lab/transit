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

import matplotlib.pyplot as plt

import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "rankproduct"
long_name = "Rank Product test for determining conditional essentiality."
description = "Differential Comparison based on ranks"
transposons = ["himar1", "tn5"]
columns = ["Orf","Name","Desc","Sites","Mean Ctrl","Mean Exp","log2FC","Obs RP","pvalue","adj. pvalue"]



############# Analysis Method ##############

class RankProductAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, description, transposons, RankProductMethod, RankProductGUI, [RankProductFile])


################## FILE ###################

class RankProductFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#RankProduct", columns)


############# GUI ##################

class RankProductGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        rankproductPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    
        rankproductSizer = wx.BoxSizer( wx.VERTICAL )

        rankproductLabel = wx.StaticText( rankproductPanel, wx.ID_ANY, u"rankproduct Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        rankproductLabel.Wrap( -1 )
        rankproductSizer.Add( rankproductLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        rankproductTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        rankproductTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        rankproductLabelSizer = wx.BoxSizer( wx.VERTICAL )

        rankproductSampleLabel = wx.StaticText( rankproductPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
        rankproductSampleLabel.Wrap( -1 )
        rankproductLabelSizer.Add( rankproductSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        rankproductNormLabel = wx.StaticText( rankproductPanel, wx.ID_ANY, u"Normalization", wx.DefaultPosition, wx.DefaultSize, 0 )
        rankproductNormLabel.Wrap( -1 )
        rankproductLabelSizer.Add( rankproductNormLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        rankproductTopSizer2.Add( rankproductLabelSizer, 1, wx.EXPAND, 5 )

        rankproductControlSizer = wx.BoxSizer( wx.VERTICAL )

        self.wxobj.rankproductSampleText = wx.TextCtrl( rankproductPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
        rankproductControlSizer.Add( self.wxobj.rankproductSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        rankproductNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        self.wxobj.rankproductNormChoice = wx.Choice( rankproductPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, rankproductNormChoiceChoices, 0 )
        self.wxobj.rankproductNormChoice.SetSelection( 0 )
        rankproductControlSizer.Add( self.wxobj.rankproductNormChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )


        rankproductTopSizer2.Add( rankproductControlSizer, 1, wx.EXPAND, 5 )

        rankproductTopSizer.Add( rankproductTopSizer2, 1, wx.EXPAND, 5 )

        rankproductSizer.Add( rankproductTopSizer, 1, wx.EXPAND, 5 )

        rankproductButton = wx.Button( rankproductPanel, wx.ID_ANY, u"Run rankproduct", wx.DefaultPosition, wx.DefaultSize, 0 )
        rankproductSizer.Add( rankproductButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

 
        rankproductPanel.SetSizer( rankproductSizer )
        rankproductPanel.Layout()
        rankproductSizer.Fit( rankproductPanel )

        #Connect events
        rankproductButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )

        self.panel = rankproductPanel



########## CLASS #######################

class RankProductMethod(base.DualConditionMethod):
    """   
    rankproduct
 
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

        base.DualConditionMethod.__init__(self, short_name, long_name, description, ctrldata, expdata, annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        




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
        samples = int(wxobj.rankproductSampleText.GetValue())
        normalization = wxobj.rankproductNormChoice.GetString(wxobj.rankproductNormChoice.GetCurrentSelection())
        replicates="Sum"
        adaptive = False
        doHistogram = False

        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        LOESS = False

        #Get output path
        defaultFileName = "rankproduct_output.dat"
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
        samples = int(kwargs.get("s", 100))
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

        self.transit_message("Starting rankproduct Method")
        start_time = time.time()
               


        Kctrl = len(self.ctrldata)
        Kexp = len(self.expdata)
        #Get orf data
        self.transit_message("Getting Data")
        (data, position) = tnseq_tools.get_data(self.ctrldata+self.expdata)
        if self.normalization != "none":
            self.transit_message("Normalizing using: %s" % self.normalization)

            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata+self.expdata, self.annotation_path)           
         

        Gctrl= tnseq_tools.Genes(self.ctrldata + self.expdata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data[:Kctrl,:], position=position)

        Gexp= tnseq_tools.Genes(self.ctrldata + self.expdata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data[Kctrl:,:], position=position)


        Ngenes = len(Gctrl)


        #print numpy.log2(30/10)


        """

        types_list = set([])
        for x in [numpy.mean(g.reads+1.0,1) for g in Gctrl]:
            types_list.add(type(x[0]))

        print [numpy.mean(g.reads+1.0,1) for g in Gctrl]
        print "FLORF", types_list
        """
            
        #print Gctrl[0].reads.shape
        #print numpy.ones((Kctrl,1)).shape

        T = numpy.ones((Kctrl,1))

        meanCtrl = numpy.zeros((Kctrl, Ngenes))
        meanExp = numpy.zeros((Kexp, Ngenes))

        for i in range(Ngenes):
            if numpy.any(Gctrl[i].reads):
                meanCtrl[:,i] = numpy.mean(Gctrl[i].reads,1)
            else:
                meanCtrl[:,i] = numpy.zeros(Kctrl)
            
            if numpy.any(Gexp[i].reads):
                meanExp[:,i] = numpy.mean(Gexp[i].reads,1)
            else:
                meanExp[:,i] = numpy.zeros(Kexp)

            
        #meanCtrl = numpy.array([numpy.mean(numpy.concatenate((g.reads, numpy.ones((Kctrl,1))),1),1) for g in Gctrl])
        #meanExp = numpy.array([numpy.mean(numpy.concatenate((g.reads, numpy.ones((Kexp,1))),1),1) for g in Gexp])

        #meanExp = numpy.nan_to_num(numpy.array([numpy.mean(g.reads+1.0,1) for g in Gexp]))

        #meanCtrl[numpy.isnan(meanCtrl)]

        #print type(meanCtrl[0][0])
        #print type(meanExp[0][0])
        #print type(meanExp[0][0]/meanCtrl[0][0])

        #print numpy.log2(meanExp[0][0]/meanCtrl[0][0])

        """
        logFC2 = numpy.zeros((Ngenes, len(meanCtrl[0])))
        for i,ratio in enumerate(meanExp/meanCtrl):
            #print i, ratio, numpy.log2(ratio)
            logFC2[i] = numpy.log2(ratio)
        """
    
        #print "meanCtrl"
        #print meanCtrl

    
        #print "meanExp"
        #print meanExp


        logFC2 = numpy.log2((meanExp+0.0001)/(meanCtrl+0.0001))
        rank = numpy.array([scipy.stats.rankdata(Lvec)/float(Ngenes) for Lvec in logFC2])
        obsRP = numpy.prod(rank,0)


        #print "logFC2"
        #print logFC2

        #print "rank"
        #print numpy.array([scipy.stats.rankdata(Lvec) for Lvec in logFC2])
        

        #print "obsRP"
        #print obsRP


        permutations = numpy.zeros((self.samples, Ngenes))
        tempranks = scipy.array([numpy.arange(1,Ngenes+1) for rep in range(Kctrl)])
        for s in range(self.samples):
            rankperm = numpy.array([numpy.random.permutation(tr)/float(Ngenes) for tr in tempranks])
            #print numpy.prod(rankperm,0)
            permutations[s] = numpy.prod(rankperm,0)
            

        #print "permutations"
        #print permutations


        rankRP = scipy.stats.rankdata(obsRP)




        #rankproduct
        data = []
        count = 0
        self.progress_range(Ngenes)
        for i,gene in enumerate(Gctrl):
            count+=1
            #print gene, obsRP[i]

            meanctrl = numpy.mean(Gctrl[i].reads)
            meanexp = numpy.mean(Gexp[i].reads)
            log2fc = numpy.log2((meanexp+0.0001)/(meanctrl+0.0001))
            #countbetter = numpy.sum(permutations < obsRP[i])
            countbetter = numpy.sum(permutations < obsRP[i])
            
            pval = countbetter/float(self.samples*Ngenes)
            q_paper = pval/rankRP[i]
 
            data.append([gene.orf, gene.name, gene.desc, gene.n, meanctrl, meanexp, log2fc, obsRP[i], q_paper, pval])
            self.progress_update("rankproduct", count)
            self.transit_message_inplace("Running rankproduct Method... %1.1f%%" % (100.0*count/Ngenes))


        #
        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Performing Benjamini-Hochberg Correction")
        data.sort() 
        qval = stat_tools.BH_fdr_correction([row[-1] for row in data])
       
 
        self.output.write("#RankProduct\n")
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
        self.output.write("#%s\n" % (columns))

        for i,row in enumerate(data):
            (orf, name, desc, n, mean1, mean2, log2FCgene, obsRPgene, q_paper, pval_2tail) = row
            self.output.write("%s\t%s\t%s\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.7f\t%1.7f\t%1.7f\n" % (orf, name, desc, n, mean1, mean2,log2FCgene, obsRPgene, pval_2tail, qval[i]))
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="RankProduct")
        self.finish()
        self.transit_message("Finished rankproduct Method") 


    @classmethod
    def usage_string(self):
        return """python %s rankproduct <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 100
        -n <string>     :=  Normalization method. Default: -n TTR
        -h              :=  Output histogram of the permutations for each gene. Default: Turned Off.
        -a              :=  Perform adaptive rankproduct. Default: Turned Off.
        -l              :=  Perform LOESS Correction; Helps remove possible genomic position bias. Default: Turned Off.
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        """ % (sys.argv[0])





    



if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)

    G = RankProductMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


