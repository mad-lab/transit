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
import tarfile
import time
import math
import random
import numpy
import statistics
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
from statsmodels.iolib.smpickle import load_pickle

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
                #STLM_reg,
                output1_file,
                output2_file,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output1_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.genome_path = genome_path
        self.gumbelestimations = gumbelestimations
        self.output2_file = output2_file
        #self.STLM_reg = STLM_reg

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

        #STLM_pickle_path = args[4]
        #with tarfile.open(STLM_pickle_path+".tar.gz", 'r') as t:
        #    t.extractall(path=".")
        #STLM_reg = sm.load(STLM_pickle_path)
        #os.remove(STLM_pickle_path)

        outpath1 = args[4]
        output1_file = open(outpath1, "w")
        outpath2 = args[5]
        output2_file = open(outpath2, "w")

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
                #STLM_reg,
                output1_file,
		output2_file,
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
        self.transit_message("Processing wig files")
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
        upseq_list = []
        downseq_list=[]
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
                upseq_list.append(upseq)
		#reverse complementing downstream
                downseq = ""
                for x in [nucs[9],nucs[8],nucs[7],nucs[6]]:
                    if(str(x)=="A"): downseq+= "T"
                    if(str(x)=="C"): downseq+= "G"
                    if(str(x)=="G"): downseq+= "C"
                    if(str(x)=="T"): downseq+= "A"
                downseq_list.append(downseq)
                ttn_vector = []
                for c in combos:
                    if upseq==c and downseq==c: ttn_vector.append(int(2)) #up/dwn ttn are same, "bit"=2
                    elif upseq==c or downseq==c:ttn_vector.append(int(1)) #set ttn bit=1
                    else:ttn_vector.append(int(0))
                ttn_vector_list.append(pandas.Series(ttn_vector,index=combos))
                orf.append(gene.orf)
                name.append(gene.name)
                coords.append(pos)

        TA_sites_df = pandas.DataFrame({"Orf":orf, "Name": name, "Coord":coords, "Insertion Count": all_counts,"Upstream TTN":upseq_list,"Downstream TTN":downseq_list})
        TA_sites_df = pandas.concat([TA_sites_df, pandas.DataFrame(ttn_vector_list)],axis=1)
        TA_sites_df=TA_sites_df.sort_values(by=['Coord'],ignore_index=True)
        # get initial states of the TA Sites
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

        ####################################################

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

        self.transit_message("\t + Filtering ES/ESB Genes")
        #function to extract gumbel calls to filter out ES and ESB
        gumbel_bernoulli_gene_calls = {}
        for g in TA_sites_df["Orf"].unique():
            if g=="igr": gene_call=numpy.nan
            else:
                gene_call='U'
                sub_gumbel=gumbel_df[gumbel_df["Orf"]==g]
                if len(sub_gumbel)>0: gene_call = sub_gumbel["Call"].iloc[0]
                #set to ES if greater than n and all 0s
                sub_data = TA_sites_df[TA_sites_df["Orf"]==g]
                if len(sub_data)>significant_n and len(sub_data[sub_data["Insertion Count"]>0])==0: gene_call="EB" #binomial filter
            gumbel_bernoulli_gene_calls[g] = gene_call
        ess_genes = [key for key,value in gumbel_bernoulli_gene_calls.items() if (value=='E') or (value=='EB')]

        self.transit_message("\t + Filtering Short Genes. Labeling as Uncertain")
        #function to call short genes (1 TA site) with no insertions as Uncertain
        uncertain_genes=[]
        for g in TA_sites_df["Orf"].unique():
            sub_data = TA_sites_df[TA_sites_df["Orf"]==g]
            len_of_gene = len(sub_data)
            num_insertions = len(sub_data[sub_data["Insertion Count"]>0])
            saturation = num_insertions/len_of_gene
            if saturation==0 and len_of_gene<=1: uncertain_genes.append(g)

        filtered_ttn_data= TA_sites_df[TA_sites_df["State"]!="ES"]
        filtered_ttn_data=  filtered_ttn_data[filtered_ttn_data["Local Average"]!=-1]
        filtered_ttn_data = filtered_ttn_data[~filtered_ttn_data["Orf"].isin(ess_genes)] #filter out ess genes
        filtered_ttn_data = filtered_ttn_data[~filtered_ttn_data["Orf"].isin(uncertain_genes)] #filter out uncertain genes
        filtered_ttn_data = filtered_ttn_data.reset_index(drop=True)
        ##########################################################################################
        # STLM Predictions
        #self.transit_message("\t + Making TTN based predictions using loaded STLM")
        #X = filtered_ttn_data.drop(["Orf","Name", "Coord","State","Insertion Count","Local Average","Actual LFC","Upseq TTN","Downseq TTN"],axis=1)
        #X = sm.add_constant(X)
        #model_LFC_predictions = self.STLM_reg.predict(X)
        #filtered_ttn_data["STLM Predicted LFC"]=model_LFC_predictions
        #filtered_ttn_data["STLM Predicted Counts"] = filtered_ttn_data["Local Average"].mul(numpy.power(2,filtered_ttn_data["STLM Predicted LFC"]))

        ##########################################################################################
        #Linear Regression
        gene_one_hot_encoded= pandas.get_dummies(filtered_ttn_data["Orf"],prefix='')
        ttn_vectors = filtered_ttn_data.drop(["Coord","Insertion Count","Orf","Name","Local Average","Actual LFC","State","Upstream TTN","Downstream TTN"],axis=1)
        #stlm_predicted_log_counts = numpy.log10(filtered_ttn_data["STLM Predicted Counts"]+0.5)

        Y = numpy.log10(filtered_ttn_data["Insertion Count"]+0.5)

        self.transit_message("\t + Fitting M1")
        X1 = pandas.concat([gene_one_hot_encoded,ttn_vectors],axis=1)
        X1 = sm.add_constant(X1)
        results1 = sm.OLS(Y,X1).fit()
        filtered_ttn_data["M1 Pred log Count"] = results1.predict(X1)
        filtered_ttn_data["M1 Predicted Count"] = numpy.power(10, (filtered_ttn_data["M1 Pred log Count"]-0.5))

        #self.transit_message("\t + Fitting new mod TTN-Fitness")
        #X2 = pandas.concat([gene_one_hot_encoded,stlm_predicted_log_counts],axis=1)
        #X2 = sm.add_constant(X2)
        #results2 = sm.OLS(Y,X2).fit()
        #filtered_ttn_data["mod ttn Pred log Count"] = results2.predict(X2)
        #filtered_ttn_data["mod ttn Predicted Count"] = numpy.power(10, (filtered_ttn_data["mod ttn Pred log Count"]-0.5))

        self.transit_message("\t + Assessing Models")
        #create Models Summary df
        Models_df = pandas.DataFrame(results1.params[1:-256],columns=["M1 Coef"])
        Models_df["M1 Pval"] = results1.pvalues[1:-256]
        Models_df["M1 Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(results1.pvalues[1:-256],alpha=0.05)[1]
        #Models_df["mod ttn Coef"] = results2.params[1:-1]
        #Models_df["mod ttn Pval"] = results2.pvalues[1:-1]
        #Models_df["mod ttn Adjusted Pval"] = statsmodels.stats.multitest.fdrcorrection(results2.pvalues[1:-1],alpha=0.05)[1]

        #creating a mask for the adjusted pvals
        Models_df.loc[(Models_df["M1 Coef"]>0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="GA"
        Models_df.loc[(Models_df["M1 Coef"]<0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="GD"
        Models_df.loc[(Models_df["M1 Coef"]==0) & (Models_df["M1 Adjusted Pval"]<0.05),"Gene+TTN States"]="NE"
        Models_df.loc[(Models_df["M1 Adjusted Pval"]>0.05),"Gene+TTN States"]="NE"

	#mask using mod TTN fitness
        #Models_df.loc[(Models_df["mod ttn Coef"]>0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="GA"
        #Models_df.loc[(Models_df["mod ttn Coef"]<0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="GD"
        #Models_df.loc[(Models_df["mod ttn Coef"]==0) & (Models_df["mod ttn Adjusted Pval"]<0.05),"mod ttn States"]="NE"
        #Models_df.loc[(Models_df["mod ttn Adjusted Pval"]>0.05),"mod ttn States"]="NE"
        #########################################################################################
        self.transit_message("Writing To Output Files")
        #Write Models Information to CSV
        # Columns: ORF ID, ORF Name, ORF Description,M0 Coef, M0 Adj Pval

        gene_dict={} #dictionary to map information per gene
        TA_sites_df["M1 Predicted Count"] = [None]*len(TA_sites_df)
        #TA_sites_df["mod ttn Predicted Count"] = [None]*len(TA_sites_df)
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
            local_saturation = above0TAsites / numTAsites
            #Predicted Count
            if g in filtered_ttn_data["Orf"].values:
                actual_df = filtered_ttn_data[filtered_ttn_data["Orf"]==g]["Insertion Count"]
                coords_orf = filtered_ttn_data[filtered_ttn_data["Orf"]==g]["Coord"].values.tolist()
                for c in coords_orf:
                    TA_sites_df.loc[(TA_sites_df["Coord"]==c),'M1 Predicted Count'] = filtered_ttn_data[filtered_ttn_data["Coord"]==c]["M1 Predicted Count"].iloc[0]
                #TA_sites_df[TA_sites_df["Coord"].isin(coords_orf)]['mod ttn Predicted Count'] = filtered_ttn_data[filtered_ttn_data["Coord"].isin(coords_orf)]["mod ttn Predicted Count"]
            #M1 info
            if "_"+g in Models_df.index:
                M1_coef = Models_df.loc["_"+g,"M1 Coef"]
                M1_adj_pval = Models_df.loc["_"+g,"M1 Adjusted Pval"]
                modified_M1 = math.exp(M1_coef - statistics.median(Models_df["M1 Coef"].values.tolist()))
                #mod_M1_coef = Models_df.loc["_"+g,"mod ttn Coef"]
                #mod_M1_adj_pval = Models_df.loc["_"+g,"mod ttn Adjusted Pval"]
                #mod_modified_M1 = math.exp(mod_M1_coef - statistics.median(Models_df["mod ttn Coef"].values.tolist()))
            else:
                M1_coef = None
                M1_adj_pval = None
                modified_M1 = None
                #mod_M1_coef = None
                #mod_M1_adj_pval = None
                #mod_modified_M1 = None

            #States
            gumbel_bernoulli_call = gumbel_bernoulli_gene_calls[g]
            if gumbel_bernoulli_call=="E":
                gene_ttn_call = "ES"
                #mod_gene_ttn_call = "ES"
            elif gumbel_bernoulli_call=="EB":
                gene_ttn_call = "ESB"
                #mod_gene_ttn_call = "ESB"
            else:
                if "_"+g in Models_df.index:
                     gene_ttn_call = Models_df.loc["_"+g,"Gene+TTN States"]
                     #mod_gene_ttn_call = Models_df.loc["_"+g,"mod ttn States"]
                else:
                    gene_ttn_call = "U" #these genes are in the uncertain genes list
                    #mod_gene_ttn_call = "U"
            TA_sites_df.loc[(TA_sites_df["Orf"]==g), 'TTN-Fitness Assessment'] = gene_ttn_call
            #TA_sites_df.loc[(TA_sites_df["Orf"]==g), 'Mod TTN-Fitness Assessment'] = mod_gene_ttn_call
            gene_dict[g] = [g,orfName,orfDescription,numTAsites,above0TAsites,local_saturation,M1_coef,M1_adj_pval, mean_actual_counts,modified_M1, gene_ttn_call]
        output_df = pandas.DataFrame.from_dict(gene_dict,orient='index')
        output_df.columns=["ORF ID","Name","Description","Total # TA Sites","#Sites with insertions","Gene Saturation","Gene+TTN (M1) Coef","Gene+TTN (M1) Adj Pval","Mean Insertion Count","Fitness Ratio","TTN-Fitness Assessment"]
        assesment_cnt = output_df["TTN-Fitness Assessment"].value_counts()
        #mod_assesment_cnt = output_df["Mod TTN-Fitness Assessment"].value_counts()

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
        self.output.write("#Saturation of Dataset: %s\n" % (saturation))
        self.output.write("#Assesment Counts: %s ES, %s ESB, %s GD, %s GA, %s NE, %s U \n" % (assesment_cnt["ES"],assesment_cnt["ESB"],assesment_cnt["GD"],assesment_cnt["GA"],assesment_cnt["NE"],assesment_cnt["U"]))
        #self.output.write("#Mod Assesment Counts: %s ES, %s ESB, %s GD, %s GA, %s NE, %s U \n" % (mod_assesment_cnt["ES"],mod_assesment_cnt["ESB"],mod_assesment_cnt["GD"],mod_assesment_cnt["GA"],mod_assesment_cnt["NE"],mod_assesment_cnt["U"]))

        TA_sites_df = TA_sites_df[["Coord","Orf","Name","Upstream TTN","Downstream TTN","TTN-Fitness Assessment","Insertion Count","Local Average","M1 Predicted Count"]]

        output2_data = TA_sites_df.to_csv(header=True,sep='\t' ,index=False).split('\n')
        vals = '\n'.join(output2_data)
        self.output2_file.write(vals)
        self.output2_file.close()

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
        return """python3 %s ttnfitness <comma-separated .wig files> <annotation .prot_table> <genome .fna> <gumbel output file> <output1 file> <output2 file>""" % (sys.argv[0])


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
