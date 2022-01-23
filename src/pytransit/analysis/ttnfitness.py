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
import pandas

from pytransit.analysis import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools

############# Description ##################

short_name = "ttnfitness"
long_name = "TTNFitness"
short_desc = "TTNFitness method that calculates mean read-counts per gene."
long_desc = "A method made to serve as an ttnfitness to implementing other methods."
transposons = ["himar1", "tn5"]
columns = ["Orf","Name","Desc","k","n","mean","nzmean"]
#columns = ["ORF ID","Name","Description","Total # TA Sites","#Sites with insertions","Used in Models","Gene (M0) Coef","Gene (M0) Adj Pval","Gene+TTN (M1) Coef","Gene+TTN (M1) Adj Pval","M0 Fitness Estimation","M1 Fitness Estimation","Mean Actual Count", "TTN-Fitness Assesment"]
############# Analysis Method ##############

class TTNFitnessAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, TTNFitnessMethod, TTNFitnessGUI, [TTNFitnessFile])


################## FILE ###################

class TTNFitnessFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#TTNFitness", columns)

    def getHeader(self, path):
        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
        return text


################# GUI ##################

class TTNFitnessGUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################

class TTNFitnessMethod(base.SingleConditionMethod):
    """
    TTNFitness

    """
    def __init__(self,
                ctrldata,
                annotation_path,
				genome_path,
                output_file,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.genome_path = genome_path


	#Overloading function so it prints out the genome fna file as well!
    def print_members(self):
        #TODO: write docstring
        members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
        for m in members:
             print("%s = %s" % (m, getattr(self, m)))


	#SC neeed to fix so genome file included
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
        if not transit_tools.validate_transposons_used(ctrldata, transposons):
            return None

        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        replicates="Sum"
        normalization = None
        LOESS = False

        #Get output path
        defaultFileName = "ttnfitness_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(ctrldata,
                annotationPath,
				genomePath,
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

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        genomePath = args[2]
        outpath = args[3]
        output_file = open(outpath, "w")

        replicates = "Sum"
        normalization = None
        LOESS = False
        ignoreCodon = True
        NTerminus = 0.0
        CTerminus = 0.0

        return self(ctrldata,
                annotationPath,
				genomePath,
                output_file,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    #read in the fna file as one continous string

    def Run(self):

        self.transit_message("Starting TTNFitness Method")
        start_time = time.time()

        #Get orf data
        self.transit_message("Getting Data")
        (data, position) = transit_tools.get_validated_data(self.ctrldata, wxobj=self.wxobj)
        (K,N) = data.shape

        if self.normalization and self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata, self.annotation_path)

        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, minread=1, reps=self.replicates, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data, position=position)
        N = len(G)

        self.transit_message("Getting Genome")
        genome = ""
        n = 0
        for line in open(self.genome_path):
            if n==0:
                n = 1 # skip first
            else:
                genome += line[:-1]

        self.transit_message("Ceating TA Insertion input to STLM")
        ############################################################
        #Creating the dataset
        orf= []
        name= []
        coords = []
        nucleos = []
        #Get the nucleotides surrounding the TA sites
        genome2 = genome+genome
        all_counts=[]
        for gene in G:
            all_counts.extend(gene.reads[0])
            for pos in gene.position:
                pos -= 1 # 1-based to 0-based indexing of nucleotides
                if pos-20<0: pos += len(genome)
                nucs=genome2[pos-20:pos+22]
                if nucs[20:22]!="TA": sys.stderr.write("warning: site %d is %s instead of TA" % (pos,nucs[20:22]))
                orf.append(gene.orf)
                name.append(gene.name)
                coords.append(pos)
                nucleos.append(nucs)
        #get initial states of the TA Sites
        # compute state labels (ES or NE)
        # for runs of >=R TA sites with cnt=0; label them as "ES", and the rest as "NE"
        # treat ends of genome as connected (circular)

        Nsites = len(all_counts)
        states = ["NE"]*Nsites
        R = 6 # make this adaptive based on saturation?
        MinCount = 2
        i = 0
        while i<Nsites:
            j = i
            while j<Nsites and all_counts[j]<MinCount:
                j += 1
            if j-i>=R:
                for k in range(i,j): states[k] = "ES"
                i = j
            else: i += 1

        #getlocal averages -excludes self
        W = 5
        localmeans = []
        for i in range(Nsites):
            vals = []
            for j in range(-W,W+1):
                if j!=0 and i+j>=0 and i+j<Nsites: # this excludes the site itself
                    if states[i+j]!=states[i]: continue # include only neighboring sites with same state when calculating localmean # diffs2.txt !!!
                    vals.append(float(all_counts[i+j]))
            smoothed = -1 if len(vals)==0 else numpy.mean(vals)
            localmeans.append(smoothed)

        #get LFCs
        LFC_values= []
        PC = 10
        for i in range(len(all_counts)):
            c,m = all_counts[i],localmeans[i]
            lfc = math.log((c+PC)/float(m+PC),2)
            LFC_values.append(lfc)

        TA_sites_df = pandas.DataFrame({"Orf":orf,
        "Name": name,
        "TA site Coord":pos,
        "TA site state": states,
        "TA site Actual Count": all_counts,
        "Local Average": localmeans,
        "LFCs": lfc,
        "Nucleotides":nucleos})
        print(TA_sites_df)
            # Update Progress
            #text = "Running TTNFitness Method... %5.1f%%" % (100.0*count/N)
            #self.progress_update(text, count)
        """
		data = []
        N = len(G)
        count = 0
        self.progress_range(N)
        for gene in G:
            count+=1
            if gene.n == 0:
                mean = 0.0
            else:
                mean = numpy.mean(gene.reads)

            if gene.k == 0:
                nzmean = 0.0
            else:
                nzmean = numpy.sum(gene.reads)/float(gene.k)

            data.append("%s\t%s\t%s\t%s\t%s\t%1.2f\t%1.2f\n" % (gene.orf, gene.name, gene.desc, gene.k, gene.n, mean, nzmean))

        """
        self.output.write("#TTNFitness\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s\n" % (",".join(self.ctrldata).encode('utf-8'), self.annotation_path.encode('utf-8'), self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python3 %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata).encode('utf-8')))
        self.output.write("#Annotation path: %s\n" % self.annotation_path.encode('utf-8'))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))

        """
        data.sort()
        for line in data:
            self.output.write(line)
        self.output.close()
        """

        self.transit_message("") # Printing empty line to flush stdout
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="TTNFitness")
        self.finish()
        self.transit_message("Finished TTNFitness Method")

    @classmethod
    def usage_string(self):
        return """python3 %s ttnfitness <comma-separated .wig files> <annotation .prot_table> <genome .fna> <output file>""" % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    print("ARGS:", args)
    print("KWARGS:", kwargs)

    G = TTNFitnessMethod.fromargs(sys.argv[1:])

    print(G)
    G.console_message("Printing the member variables:")

    G.print_members()

    print("")
    print("Running:")

    G.Run()
