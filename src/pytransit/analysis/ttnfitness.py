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
import itertools
import statsmodels.api as sm
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import statsmodels.stats.multitest
from sklearn.preprocessing import OneHotEncoder

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
#columns = ["ORF ID","Name","Description","Total # TA Sites","#Sites with insertions","Used in Models","Gene (M0) Coef","Gene (M0) Adj Pval","Gene+TTN (M1) Coef","Gene+TTN (M1) Adj Pval","M0 Fitness Estimation","M1 Fitness Estimation","Mean Insertion Count", "TTN-Fitness Assesment"]
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
                gumbelestimations,
                output_file,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.genome_path = genome_path
        self.gumbelestimations = gumbelestimations


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
        gumbelestimations = args[3]
        outpath = args[4]
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
                gumbelestimations,
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

        ### Genome Read -in
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
        self.transit_message("Calculating input to STLM")
        ############################################################
        #Creating the dataset
        orf= []
        name= []
        coords = []
        ttn_vector_list = []
        #Get the nucleotides surrounding the TA sites
        genome2 = genome+genome
        all_counts=[]
        combos=[''.join(p) for p in itertools.product(['A','C','T','G'], repeat=4)]
        gene_obj_dict = {}
        for gene in G:
            gene_obj_dict[gene.orf]= gene
            all_counts.extend(numpy.mean(gene.reads,0)) #mean TA site counts across wig files
            nTAs = len(gene.reads[0])
            for pos in gene.position:
                pos -= 1 # 1-based to 0-based indexing of nucleotides
                if pos-4<0: pos += len(genome)
                nucs=genome2[pos-4:pos+6]
                if nucs[4:6]!="TA":sys.stderr.write("warning: site %d is %s instead of TA" % (pos,nucs[4:6]))
                #convert nucleotides to upstream and downstream TTN
                upseq = nucs[0]+nucs[1]+nucs[2]+nucs[3]
                #reverse complementing downstream
                downseq = ""
                for x in [nucs[9],nucs[8],nucs[7],nucs[6]]:
                    if(str(x)=="A"): downseq+= "T"
                    if(str(x)=="C"): downseq+= "G"
                    if(str(x)=="G"): downseq+= "C"
                    if(str(x)=="T"): downseq+= "A"
                ttn_vector = []
                for c in combos:
                    if upseq==c and downseq==c: ttn_vector.append(int(2)) #up/dwn ttn are same, "bit"=2
                    elif upseq==c or downseq==c:ttn_vector.append(int(1)) #set ttn bit=1
                    else:ttn_vector.append(int(0))
                ttn_vector_list.append(pandas.Series(ttn_vector,index=combos))

                orf.append(gene.orf)
                name.append(gene.name)
                coords.append(pos)

        TA_sites_df = pandas.DataFrame({"Orf":orf, "Name": name, "Coord":pos, "Insertion Count": all_counts})
        TA_sites_df = pandas.concat([TA_sites_df, pandas.DataFrame(ttn_vector_list)],axis=1)
        TA_sites_df=TA_sites_df.sort_values(by=['Coord'],ignore_index=True)


        #get initial states of the TA Sites
        # compute state labels (ES or NE)
        # for runs of >=R TA sites with cnt=0; label them as "ES", and the rest as "NE"
        # treat ends of genome as connected (circular)
        Nsites = len(TA_sites_df["Insertion Count"])
        states = ["NE"]*Nsites
        R = 6 # make this adaptive based on saturation?
        MinCount = 2
        i = 0
        while i<Nsites:
            j = i
            while j<Nsites and TA_sites_df["Insertion Count"].iloc[j]<MinCount:
                j += 1
            if j-i>=R:
                for k in range(i,j): states[k] = "ES"
                i = j
            else: i += 1

        #getlocal averages --excludes self
        W = 5
        localmeans = []
        for i in range(Nsites):
            vals = []
            for j in range(-W,W+1):
                if j!=0 and i+j>=0 and i+j<Nsites: # this excludes the site itself
                    if states[i+j]!=states[i]: continue # include only neighboring sites with same state when calculating localmean # diffs2.txt !!!
                    vals.append(float(TA_sites_df["Insertion Count"].iloc[i+j]))
            smoothed = -1 if len(vals)==0 else numpy.mean(vals)
            localmeans.append(smoothed)

        #get LFCs
        LFC_values= []
        PC = 10
        for i in range(len(TA_sites_df["Insertion Count"])):
            c,m = TA_sites_df["Insertion Count"].iloc[i],localmeans[i]
            lfc = math.log((c+PC)/float(m+PC),2)
            LFC_values.append(lfc)

        TA_sites_df["State"] =states
        TA_sites_df["Local Average"] = localmeans
        TA_sites_df["Actual LFC"] = LFC_values


        self.transit_message("Training the STLM")
        filtered_TA_sites_df= TA_sites_df[TA_sites_df["State"]!="ES"]
        filtered_TA_sites_df.reset_index(inplace=True, drop=True)
        y = filtered_TA_sites_df["Actual LFC"]
        X = filtered_TA_sites_df.drop(["Orf","Name", "Coord","State","Insertion Count","Local Average","Actual LFC"],axis=1)

        #Fit the model and make the TTN predictions
        X = sm.add_constant(X)
        STLM_model = sm.OLS(y,X).fit()

        ####################################################

        self.transit_message("Getting Expected Counts from STLM")
        model_LFC_predictions = STLM_model.predict(X)

        b=0
        LFC_Predictions= [] # to add to the dataframe
        Count_Predictions= []
        for k in range(len(TA_sites_df["Actual LFC"])):
            if TA_sites_df["State"].iloc[k] =="ES":
                LFC_Predictions.append(numpy.nan)
                Count_Predictions.append(numpy.nan)
            else:
                predCount = TA_sites_df["Local Average"].iloc[k]*math.pow(2,model_LFC_predictions[b])
                LFC_Predictions.append(model_LFC_predictions[b])
                Count_Predictions.append(predCount)
                b = b+1

        TA_sites_df["STLM Predicted LFC"] = LFC_Predictions
        TA_sites_df["STLM Predicted Counts"] = Count_Predictions

        self.transit_message("Making Fitness Estimations")
        #Read in Gumbel estimations
        skip_count =0
        gumbel_file = open(self.gumbelestimations,'r')
        for line in gumbel_file.readlines():
                if line.startswith('#'):skip_count = skip_count+1
                else:break
        gumbel_file.close()
        gumbel_df = pandas.read_csv(self.gumbelestimations,sep='\t',skiprows=skip_count, names =["Orf","Name","Desc","k","n","r","s","zbar","Call"],dtype = str)

        saturation = len(TA_sites_df[TA_sites_df["Insertion Count"]>0])/len(TA_sites_df)
        phi = 1.0 - saturation
        significant_n = math.log10(0.05)/math.log10(phi)


        gumbel_gene_calls = {}
        for g in TA_sites_df["Orf"].unique():
            if g=="igr": gene_call=numpy.nan
            else:
                gene_call='U'
                sub_gumbel=gumbel_df[gumbel_df["Orf"]==g]
                if len(sub_gumbel)>0: gene_call = sub_gumbel["Call"].iloc[0]
                #set to ES if greater than n and all 0s
                sub_data = TA_sites_df[TA_sites_df["Orf"]==g]
                if len(sub_data)>significant_n and len(sub_data[sub_data["Insertion Count"]>0])==0: gene_call="EB"
            gumbel_gene_calls[g] = gene_call


        ess_genes = [key for key,value in gumbel_gene_calls.items() if (value=='E') or (value=='EB')]
        filtered_ttn_data = TA_sites_df[~TA_sites_df["Orf"].isin(ess_genes)]
        filtered_ttn_data = filtered_ttn_data.reset_index(drop=True)

        ##########################################################################################
        #Linear Regression
        gene_one_hot_encoded= pandas.get_dummies(filtered_ttn_data["Orf"],prefix='')
        ttn_vectors = filtered_ttn_data.drop(["Coord","Insertion Count","Orf","Name","Local Average","Actual LFC","State", "STLM Predicted LFC","STLM Predicted Counts"],axis=1)

        X0 = pandas.concat([gene_one_hot_encoded],axis=1)
        X0 = sm.add_constant(X0)

        X1 = pandas.concat([gene_one_hot_encoded,ttn_vectors],axis=1)
        X1 = sm.add_constant(X1)
        Y = numpy.log10(filtered_ttn_data["Insertion Count"]+0.5)

        results0 = sm.OLS(Y,X0).fit()
        results1 = sm.OLS(Y,X1).fit()

        filtered_ttn_data["M0 Pred log Count"] = results0.predict(X0)
        filtered_ttn_data["M1 Pred log Count"] = results1.predict(X1)

        predM0Count = []
        predM1Count = []
        for k in range(len(filtered_ttn_data["M0 Pred log Count"])):
            predM0Count.append(math.pow(10,filtered_ttn_data["M0 Pred log Count"].iloc[k])-0.5)
            predM1Count.append(math.pow(10,filtered_ttn_data["M1 Pred log Count"].iloc[k])-0.5)

        filtered_ttn_data["M0 Predicted Count"] = predM0Count
        filtered_ttn_data["M1 Predicted Count"] = predM1Count


        #create Models Summary df
        Models_df = pandas.DataFrame(results0.params[1:],columns=["M0 Coef"])
        Models_df["M0 Pval"] = results0.pvalues[1:]
        Models_df["M0 Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(results0.pvalues[1:],alpha=0.05)[1]
        Models_df["M1 Coef"]= results1.params[1:-256]
        Models_df["M1 Pval"] = results1.pvalues[1:-256]
        Models_df["M1 Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(results1.pvalues[1:-256],alpha=0.05)[1]

        #creating a mask for the adjusted pvals
        Models_df.loc[(Models_df["M1 Coef"]>0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="GA"
        Models_df.loc[(Models_df["M1 Coef"]<0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="GD"
        Models_df.loc[(Models_df["M1 Coef"]==0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="NE"
        Models_df.loc[(Models_df["M1 Adjusted Pval"]>0.05),"Gene+TTN States"]="NE"

        #########################################################################################
        self.transit_message("Write Models Information to CSV")
        #Write Models Information to CSV
        # Columns: ORF ID, ORF Name, ORF Description,M0 Coef, M0 Adj Pval
        gene_dict={} #dictionary to map information per gene
        for g in TA_sites_df["Orf"].unique():
            #ORF Name
            orfName = gene_obj_dict[g].name
            #ORF Description
            orfDescription = gene_obj_dict[g].desc
            #Total TA sites
            numTAsites = len(gene_obj_dict[g].reads[0]) #TRI check this!
            #Sites > 0
            above0TAsites = len([r for r in gene_obj_dict[g].reads[0] if r>0])
            #Insertion Count
            actual_counts = TA_sites_df[TA_sites_df["Orf"]==g]["Insertion Count"]
            mean_actual_counts = numpy.mean(actual_counts)
            #Predicted Count
            if g not in TA_sites_df["Orf"].values or TA_sites_df[TA_sites_df["Orf"]==g]["State"].iloc[0]=="ES":
                M0_ratio=None
                M1_ratio = None
            else:
                M0_ratio = (mean_actual_counts/numpy.mean(filtered_ttn_data[filtered_ttn_data["Orf"]==g]["M0 Predicted Count"]))
                M1_ratio = (mean_actual_counts/numpy.mean(filtered_ttn_data[filtered_ttn_data["Orf"]==g]["M1 Predicted Count"]))
            #M0/M1 info
            if "_"+g in Models_df.index:
                used=True
                M0_coef= Models_df.loc["_"+g,"M0 Coef"]
                M0_adj_pval = Models_df.loc["_"+g,"M0 Adjusted Pval"]
                M1_coef = Models_df.loc["_"+g,"M1 Coef"]
                M1_adj_pval = Models_df.loc["_"+g,"M1 Adjusted Pval"]

            else:
                used=False
                M0_coef= None
                M0_adj_pval = None
                M1_coef = None
                M1_adj_pval = None

            #States
            gumbel_bernoulli_call = gumbel_gene_calls[g]
            if gumbel_bernoulli_call=="E": gene_ttn_call = "ES"
            elif gumbel_bernoulli_call=="EB": gene_ttn_call = "ESB"
            else:
                if "_"+g in Models_df.index: gene_ttn_call = Models_df.loc["_"+g,"Gene+TTN States"]
                else: gene_ttn_call = "Uncertain"
            gene_dict[g] = [g,orfName,orfDescription,numTAsites,above0TAsites,used,M0_coef,M0_adj_pval,M1_coef,M1_adj_pval,M0_ratio,M1_ratio,mean_actual_counts,gene_ttn_call]

        output_df = pandas.DataFrame.from_dict(gene_dict,orient='index')
        output_df.columns=["ORF ID","Name","Description","Total # TA Sites","#Sites with insertions","Used in Models","Gene (M0) Coef","Gene (M0) Adj Pval","Gene+TTN (M1) Coef","Gene+TTN (M1) Adj Pval","M0 Fitness Estimation","M1 Fitness Estimation","Mean Insertion Count", "TTN-Fitness Assesment"]

        assesment_cnt = output_df["TTN-Fitness Assesment"].value_counts()


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
        self.output.write("Saturation of Dataset: %s" % (saturation))
        self.output.write("#Assesment Counts: %s ES, %s ESB, %s GD, %s GA, %s NE \n" % (assesment_cnt["ES"],assesment_cnt["ESB"],assesment_cnt["GD"],assesment_cnt["GA"],assesment_cnt["NE"]))

        output_data = output_df.to_csv(header=True, sep="\t", index=False).split('\n')
        vals = '\n'.join(output_data)
        self.output.write(vals)
        self.output.close()


        self.transit_message("") # Printing empty line to flush stdout
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="TTNFitness")
        self.finish()
        self.transit_message("Finished TTNFitness Method")

    @classmethod
    def usage_string(self):
        return """python3 %s ttnfitness <comma-separated .wig files> <annotation .prot_table> <genome .fna> <gumbel output file> <output file>""" % (sys.argv[0])


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
