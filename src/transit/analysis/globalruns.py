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
import transit
import transit.transit_tools as transit_tools
import transit.tnseq_tools as tnseq_tools
import transit.norm_tools as norm_tools
import transit.stat_tools as stat_tools

#method_name = "example"


############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.globalGumbelPanel.Hide()

def Show(wxobj):
    wxobj.globalGumbelPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run Global Gumbel" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.globalGumbelPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    #wxobj.examplePanel.SetMinSize( wx.Size( 50,1 ) )
    #wxobj.examplePanel.SetMaxSize( wx.Size( 250,-1 ) )

    globalGumbelSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.globalGumbelLabel = wx.StaticText( wxobj.globalGumbelPanel, wx.ID_ANY, u"Global Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.globalGumbelLabel.Wrap( -1 )
    globalGumbelSection.Add( wxobj.globalGumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    globalGumbelSizer1 = wx.BoxSizer( wx.HORIZONTAL )
    #exampleSizer2 = wx.BoxSizer( wx.HORIZONTAL )
    #exampleLabelSizer = wx.BoxSizer( wx.VERTICAL )
    #exampleControlSizer = wx.BoxSizer( wx.VERTICAL )
    
    
    #wxobj.exampleRepLabel = wx.StaticText( wxobj.examplePanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
    #wxobj.exampleRepLabel.Wrap(-1)
    #exampleLabelSizer.Add(wxobj.exampleRepLabel, 1, wx.ALL, 5)
    #exampleRepChoiceChoices = [ u"Sum", u"Mean", "TTRMean" ]
    #wxobj.exampleRepChoice = wx.Choice( wxobj.examplePanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, exampleRepChoiceChoices, 0 )
    #wxobj.exampleRepChoice.SetSelection( 2 )
    #exampleControlSizer.Add(wxobj.exampleRepChoice, 0, wx.ALL|wx.EXPAND, 5)
    #exampleSizer2.Add(exampleLabelSizer, 1, wx.EXPAND, 5)
    #exampleSizer2.Add(exampleControlSizer, 1, wx.EXPAND, 5)
    #exampleSizer1.Add(exampleSizer2, 1, wx.EXPAND, 5 )
    
    # Min read option
    globalGumbelSizer2 = wx.BoxSizer( wx.HORIZONTAL )
    
    globalGumbelSizer2_1 = wx.BoxSizer( wx.VERTICAL )
    wxobj.globalGumbelReadLabel = wx.StaticText( wxobj.globalGumbelPanel, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.globalGumbelReadLabel.Wrap( -1 )
    globalGumbelSizer2_1.Add( wxobj.globalGumbelReadLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
    
    globalGumbelSizer2_2 = wx.BoxSizer( wx.VERTICAL )
    globalGumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
    wxobj.globalGumbelReadChoice = wx.Choice( wxobj.globalGumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, globalGumbelReadChoiceChoices, 0 )
    wxobj.globalGumbelReadChoice.SetSelection( 0 )
    globalGumbelSizer2_2.Add( wxobj.globalGumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    globalGumbelSizer2.Add(globalGumbelSizer2_1, 1, wx.EXPAND, 5)
    globalGumbelSizer2.Add(globalGumbelSizer2_2, 1, wx.EXPAND, 5)
    
    # Replicates option
    globalGumbelSizer3 = wx.BoxSizer( wx.HORIZONTAL )
    
    globalGumbelSizer3_1 = wx.BoxSizer( wx.VERTICAL )
    wxobj.globalGumbelRepLabel = wx.StaticText( wxobj.globalGumbelPanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.globalGumbelRepLabel.Wrap( -1 )
    globalGumbelSizer3_1.Add( wxobj.globalGumbelRepLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
    
    globalGumbelSizer3_2 = wx.BoxSizer( wx.VERTICAL )
    globalGumbelRepChoiceChoices = [ u"Sum", u"Mean" ]
    wxobj.globalGumbelRepChoice = wx.Choice( wxobj.globalGumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, globalGumbelRepChoiceChoices, 0 )
    wxobj.globalGumbelRepChoice.SetSelection( 0 )
    globalGumbelSizer3_2.Add( wxobj.globalGumbelRepChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    globalGumbelSizer3.Add(globalGumbelSizer3_1, 1, wx.EXPAND, 5)
    globalGumbelSizer3.Add(globalGumbelSizer3_2, 1, wx.EXPAND, 5)

    globalGumbelSection.Add( globalGumbelSizer1, 1, wx.EXPAND, 5 )
    globalGumbelSection.Add( globalGumbelSizer2, 1, wx.EXPAND, 5 )
    globalGumbelSection.Add( globalGumbelSizer3, 1, wx.EXPAND, 5 )

    wxobj.globalGumbelButton = wx.Button( wxobj.globalGumbelPanel, wx.ID_ANY, u"Run Global Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
    globalGumbelSection.Add( wxobj.globalGumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.globalGumbelPanel.SetSizer( globalGumbelSection )
    wxobj.globalGumbelPanel.Layout()
    globalGumbelSection.Fit( wxobj.globalGumbelPanel )

    #Connect events
    wxobj.globalGumbelButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.globalGumbelPanel


def updateProgressBar(wxobj, count):
    wxobj.globalGumbelProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.globalGumbelProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.globalGumbelButton.Enable()





########## CLASS #######################

class GlobalGumbel(base.SingleConditionMethod):
    """   
    Example
 
    """
    def __init__(self,
                ctrldata,
                annotation_path,
                output_file,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                minread=1,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, "Global Gumbel", "Global Gumbel Method", "A analysis method using the Gumbel extreme value distribution that uses longest runs over a whole genome instead of individual genes.", ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.minread = minread


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
        minread = int(wxobj.globalGumbelReadChoice.GetString(wxobj.globalGumbelReadChoice.GetCurrentSelection()))
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        replicates="Sum"
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(all_selected[0])
        defaultFileName = "global_gumbel_output.dat"
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
                minread,
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
        minread = 1
        NTerminus = 0.0
        CTerminus = 0.0

        return self(ctrldata,
                annotationPath,
                output_file,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                minread,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting global gumbel method")
        start_time = time.time()
        
        self.transit_message("Getting data")
        genes_obj = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)
        
        # Combine all wigs
        (data,position) = tnseq_tools.get_data(self.ctrldata)
        combined = tnseq_tools.combine_replicates(data, method=self.replicates)
        combined[combined < self.minread] = 0
        counts = combined
        counts[counts > 0] = 1
        num_sites = counts.size
        
        pins = numpy.mean(counts)
        pnon = 1.0 - pins

        # Calculate stats of runs
        exprunmax = tnseq_tools.ExpectedRuns(num_sites, pnon)
        varrun = tnseq_tools.VarR(num_sites, pnon)
        stddevrun = math.sqrt(varrun)
        exp_cutoff = exprunmax + 2*stddevrun

        # Get the runs
        self.transit_message("Getting non-insertion runs in genome")
        run_arr = tnseq_tools.runs_w_info(counts)
        pos_hash = tnseq_tools.get_pos_hash(self.annotation_path)

        # Finally, calculate the results
        self.transit_message("Running global gumbel method")
        results_per_gene = {}
        for gene in genes_obj.genes:
            results_per_gene[gene.orf] = [gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, 0, 0, 1]
        
        N = len(run_arr)
        count = 0
        accum = 0
        self.progress_range(N)
        for run in run_arr: 
            accum += run['length']
            count += 1
            genes = tnseq_tools.get_genes_in_range(pos_hash, run['start'], run['end'])
            for gene_orf in genes:
                gene = genes_obj[gene_orf]
                inter_sz = self.intersect_size([run['start'], run['end']], [gene.start, gene.end]) + 1
                percent_overlap = self.calc_overlap([run['start'], run['end']], [gene.start, gene.end])
                run_len = run['length']
                B = 1.0/math.log(1.0/pnon)
                u = math.log(num_sites*pins, 1.0/pnon)
                pval = 1.0 - tnseq_tools.GumbelCDF(run['length'], u, B)
                
                curr_val = results_per_gene[gene.orf]
                curr_inter_sz = curr_val[6]
                curr_len = curr_val[7]
                if inter_sz > curr_inter_sz:
                    results_per_gene[gene.orf] = [gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, inter_sz, run_len, pval]
            self.progress_update("globalGumbel", count)
            self.transit_message_inplace("Running global gumbel method... %1.1f%%" % (100.0*count/N))
                
        data = list(results_per_gene.values())
        exp_run_len = float(accum)/N
        
        min_sig_len = float('inf')
        sig_genes_count = 0
        pval = [row[-1] for row in data]
        padj = stat_tools.BH_fdr_correction(pval)
        for i in range(len(data)):
            if padj[i] < 0.05:
                sig_genes_count += 1
                min_sig_len = min(min_sig_len, data[i][-2])
            data[i].append(padj[i]);
            data[i].append('Essential' if padj[i] < 0.05 else 'Non-essential');#(data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6], data[i][7], data[i][8], padj[i], 'Essential' if padj[i] < 0.05 else 'Non-essential')
        data.sort(key=lambda l: l[0])
            
        # Output results
        self.output.write("#Global Gumbel\n")
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
        self.output.write("#Essential gene count: %d\n" % (sig_genes_count))
        self.output.write("#Minimum reads: %d\n" % (self.minread))
        self.output.write("#Replicate combination method: %s\n" % (self.replicates))
        self.output.write("#Minimum significant run length: %d\n" % (min_sig_len))
        self.output.write("#Expected run length: %1.5f\n" % (exp_run_len))
        self.output.write("#Expected max run length: %s\n" % (exprunmax))
        self.output.write("#Orf\tName\tDesc\tk\tn\tr\tovr\tlenovr\tpval\tpadj\tcall\n")

        for res in data:
            self.output.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%1.5f\t%1.5f\t%s\n" % (res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10]))
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.transit_message("Finished Global Gumbel Method") 

    @classmethod
    def usage_string(self):
        return """python %s globalruns <comma-separated .wig files> <annotation .prot_table> <output file>""" % (sys.argv[0])



    def intersect_size(self, intv1, intv2):
        first = intv1 if intv1[0] < intv2[0] else intv2
        second = intv1 if first == intv2 else intv2

        if first[1] < second[0]:
            return 0
        
        right_ovr = min(first[1], second[1])
        left_ovr = max(first[0], second[0])

        return right_ovr - left_ovr


    def calc_overlap(self, run_interv, gene_interv):
        first = run_interv if run_interv[0] < gene_interv[0] else gene_interv
        second = run_interv if first == gene_interv else gene_interv
        
        if second[0] > first[0] and second[1] < first[1]:
            return 1.0
        
        intersect = self.intersect_size(run_interv, gene_interv)
        return float(intersect)/(gene_interv[1] - gene_interv[0])




if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    print "ARGS:", args
    print "KWARGS:", kwargs

    G = Example.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


