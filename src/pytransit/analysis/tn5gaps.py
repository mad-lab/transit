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
import math
import random
import numpy
import scipy.stats
import datetime

import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools

#method_name = "example"


############# GUI ELEMENTS ##################

short_name = "tn5gaps"
long_name = "Tn5 Gaps"
short_desc = "Analysis of essentiality on gaps in entire genome (Tn5)."
long_desc = "A analysis method based on the extreme value (Gumbel) distribution that considers longest runs over the whole genome instead of individual genes."
transposons = ["tn5"]
columns = ["Orf","Name","Desc","k","n","r","ovr","lenovr","pval","padj","call"]



############# Analysis Method ##############

class Tn5GapsAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, Tn5GapsMethod, Tn5GapsGUI, [Tn5GapsFile])


################## FILE ###################

class Tn5GapsFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#Tn5 Gaps", columns)

    def getHeader(self, path):
        ess=0; unc=0; non=0; short=0
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[-1] == "Essential": ess+=1
            if tmp[-1] == "Non-essential": non+=1

        text = """Results:
    Essentials: %s
    Non-Essential: %s
            """ % (ess, non)
        return text


################## GUI ###################

class Tn5GapsGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        tn5GapsPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        tn5GapsSection = wx.BoxSizer( wx.VERTICAL )

        tn5GapsLabel = wx.StaticText( tn5GapsPanel, wx.ID_ANY, u"Tn5 Gaps Options", wx.DefaultPosition, (150,-1), 0 )
        tn5GapsLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
        tn5GapsSection.Add( tn5GapsLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        mainSizer1 = wx.BoxSizer( wx.VERTICAL )   
 

        # Min Read
        tn5GapsReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
        (tn5GapsReadLabel, self.wxobj.tn5GapsReadChoice, readSizer) = self.defineChoiceBox(tn5GapsPanel, u"Minimum Read:", tn5GapsReadChoiceChoices, "This is the minimum number of reads to consider a 'true' insertion. Value of 1 will consider all insertions. Larger values allow the method to ignore spurious insertions which might interrupt a run of non-insertions. Noisy datasets or those with many replicates can beneffit from increasing this.")
        mainSizer1.Add(readSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )
   

        # Replicates
        tn5GapsRepChoiceChoices = [ u"Sum", u"Mean" ]
        (tn5GapsRepLabel, self.wxobj.tn5GapsRepChoice, repSizer) = self.defineChoiceBox(tn5GapsPanel, u"Replicates:", tn5GapsRepChoiceChoices, "Determines how to handle replicates, and their read-counts. When using many replicates, summing read-counts may make spurious counts appear to be significantly large and interrupt a run of non-insertions.")
        mainSizer1.Add(repSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )


        tn5GapsSection.Add( mainSizer1, 1, wx.EXPAND, 5 )
 

        tn5GapsButton = wx.Button( tn5GapsPanel, wx.ID_ANY, u"Run Tn5Gaps", wx.DefaultPosition, wx.DefaultSize, 0 )
        tn5GapsSection.Add( tn5GapsButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        tn5GapsPanel.SetSizer( tn5GapsSection )
        tn5GapsPanel.Layout()
        tn5GapsSection.Fit( tn5GapsPanel )
        
        #Connect events
        tn5GapsButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )


        self.panel = tn5GapsPanel




########## CLASS #######################

class Tn5GapsMethod(base.SingleConditionMethod):
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

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.minread = minread



    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        #Validate transposon types
        types = tnseq_tools.get_file_types(ctrldata)
        if 'himar1' in types:
            answer = transit_tools.ShowAskWarning("Warning: One of the selected wig files looks like a Himar1 dataset. This method is designed to work on Tn5 wig files. Proceeding will fill in missing data with zeroes. Click OK to continue.")
            if answer == wx.ID_CANCEL:
                return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        minread = int(wxobj.tn5GapsReadChoice.GetString(wxobj.tn5GapsReadChoice.GetCurrentSelection()))
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        replicates= wxobj.tn5GapsRepChoice.GetString(wxobj.tn5GapsRepChoice.GetCurrentSelection())
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "tn5_gaps_output_m%d_r%s.dat" % (minread, replicates)
        

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

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        replicates = kwargs.get("r", "Sum")
        minread = int(kwargs.get("m", 1))
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
                minread,
                NTerminus,
                CTerminus)




    def Run(self):
        self.transit_message("Starting Tn5 gaps method")
        start_time = time.time()
        
        self.transit_message("Getting data (May take a while)")
        
        # Combine all wigs
        (data,position) = transit_tools.get_validated_data(self.ctrldata, wxobj=self.wxobj)
        combined = tnseq_tools.combine_replicates(data, method=self.replicates)
        combined[combined < self.minread] = 0
        counts = combined
        counts[counts > 0] = 1
        num_sites = counts.size
        
        genes_obj = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data, position=position)
        
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
        pos_hash = transit_tools.get_pos_hash(self.annotation_path)

        # Finally, calculate the results
        self.transit_message("Running Tn5 gaps method")
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
            
            # Update Progress
            text = "Running Tn5Gaps method... %1.1f%%" % (100.0*count/N) 
            self.progress_update(text, count)
                
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
        self.output.write("#Tn5 Gaps\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata).encode('utf-8'), self.annotation_path.encode('utf-8'), self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode('utf-8'))) 
        self.output.write("#Annotation path: %s\n" % self.annotation_path.encode('utf-8')) 
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#Essential gene count: %d\n" % (sig_genes_count))
        self.output.write("#Minimum reads: %d\n" % (self.minread))
        self.output.write("#Replicate combination method: %s\n" % (self.replicates))
        self.output.write("#Minimum significant run length: %d\n" % (min_sig_len))
        self.output.write("#Expected run length: %1.5f\n" % (exp_run_len))
        self.output.write("#Expected max run length: %s\n" % (exprunmax))
        self.output.write("#%s\n" % "\t".join(columns))
        #self.output.write("#Orf\tName\tDesc\tk\tn\tr\tovr\tlenovr\tpval\tpadj\tcall\n")

        for res in data:
            self.output.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%1.5f\t%1.5f\t%s\n" % (res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10]))
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Tn5 Gaps")
        self.finish()
        self.transit_message("Finished Tn5Gaps Method") 


    @classmethod
    def usage_string(self):
        return """python %s resampling <comma-separated .wig control files> <comma-separated .wig experimental files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
    
        Optional Arguments:
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        """ % (sys.argv[0])


        



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

    G = Tn5GapsMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


