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

hasR = False
try:
    import rpy2.robjects
    hasR = True
except Exception as e:
    hasR = False

if hasR:
    from rpy2.robjects import r, DataFrame, globalenv, IntVector, FloatVector, StrVector, packages as rpackages

############# Description ##################

short_name = "CGI"
long_name = "Chemical Genetic Analysis"
short_desc = "CGI Analysis of CRISPRi libraries"
long_desc = "CGI Analysis of CRISPRi libraries"
transposons = []

columns = ["Position","Reads","Genes"] # ???

############# Analysis Method ##############

class CGI(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, CGI_Method, CGI_GUI, []) 

################## FILE ###################

# there is no output file that could be loaded into the GUI

#class HeatmapFile(base.TransitFile):
#
#    def __init__(self):
#        base.TransitFile.__init__(self, "#CombinedWig", columns) 
#
#    def getHeader(self, path):
#        text = """This is file contains mean counts for each gene. Nzmean is mean accross non-zero sites."""
#        return text

################# GUI ##################

# right now, CGI is just intended for the command-line; TRI

class CGI_GUI(base.AnalysisGUI):

    def __init__(self):
        base.AnalysisGUI.__init__(self)

########## METHOD #######################

class CGI_Method(base.SingleConditionMethod):
    def __init__(self):
                ctrldata=None # initializers for superclass
                annotation_path=""
                output_file=""
                replicates="Sum"
                normalization="nonorm" 
                LOESS=False
                ignoreCodon=True
                NTerminus=0.0
                CTerminus=0.0
                wxobj=None
                # this initialization seems pointless for CGI, but must do this for base class...
                base.SingleConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

    @classmethod
    def usage_string(self):
        return """usage (6 sub-commands):
    python3 ../src/transit.py CGI extract_counts <fastq file> <ids file> > <counts file>
    python3 ../src/transit.py CGI create_combined_counts <comma seperated headers> <counts file 1> <counts file 2> ... <counts file n> > <combined counts file>
    python3 ../src/transit.py CGI extract_abund <combined counts file> <metadata file> <control condition> <sgRNA efficiency file> <uninduced ATC file> <drug> <days>  >  <fractional abundundance file>
    python3 ../src/transit.py CGI run_model <fractional abundundance file>  >  <CRISPRi DR results file>
    python3 ../src/transit.py CGI visualize <fractional abundance> <gene> <output figure location>
    note: redirect output from stdout to output files as shown above"""


    @classmethod
    def fromargs(self, rawargs): 
        if not hasR:
            print("Error: R and rpy2 (~= 3.0) required to run heatmap.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)
        if len(args)<1: print(self.usage_string())
        self.cmd = args[0]
        self.args = args[1:]
        self.kwargs = kwargs

        return self()

    def Run(self):
        print("Note: CRISPRi-DR (CGI) has been migrated to Transit2.  Please use Transit2 for this.")
        sys.exit(-1)
        cmd,args,kwargs = self.cmd,self.args,self.kwargs

        if cmd=="extract_counts":
            if len(args)<2: 
                print("You have provided incorrect number of args")
                print(self.usage_string())
            fastq_file = args[0]
            ids_file = args[1]
            self.extract_counts(fastq_file, ids_file)
        
        elif cmd=="create_combined_counts":
            if len(args)<2: 
                print("You have provided incorrect number of args")
                print(self.usage_string())
            headers = args[0].split(",")
            counts_file_list = args[1:]
            self.create_combined_counts(headers,counts_file_list)

        elif cmd=="extract_abund":
            if len(args)<7: 
                print("You have provided incorrect number of args")
                print(self.usage_string())
            combined_counts_file = args[0]
            metadata_file = args[1]
            control_condition=args[2]
            sgRNA_efficiency_file = args[3]
            no_dep_abund = args[4]
            drug = args[5]
            days = args[6]
            self.extract_abund(combined_counts_file,metadata_file,control_condition,sgRNA_efficiency_file,no_dep_abund,drug,days)
        elif cmd == "run_model":
            if len(args)<1: 
                print("You have provided incorrect number of args")
                print(self.usage_string())
            ifile_path = args[0] #example frac_abund_RIF_D5.txt
            self.run_model(ifile_path)
        elif cmd == "visualize":
            if len(args)<3: 
                print("You have provided incorrect number of args")
                print(self.usage_string())
            frac_abund_file= args[0]
            gene = args[1]
            fig_location = args[2]
            self.visualize(frac_abund_file, gene, fig_location)
            
        else: 
            print("You have not entered a valid command, here are the options")
            print(self.usage_string())

    def reverse_complement(self, seq):
        complement = {'A':'T','T':'A','C':'G','G':'C'}
        s = list(seq)
        s.reverse()
        for i in range(len(s)):
            s[i] = complement.get(s[i],s[i]) # if unknown, leave as it, e.g > or !
        s = ''.join(s)
        return s

    def extract_counts(self, fastq_file, ids_file):
        IDs = []
        barcodemap = {} # hash from barcode to full ids
        for line in open(ids_file):
            w = line.rstrip().split('\t')
            id = w[0]
            v = id.split('_')
            if len(v)<3: continue
            barcode = v[2]
            IDs.append(id)
            # reverse-complement of barcodes appears in reads, so hash them that way
            barcodemap[self.reverse_complement(barcode)] = id

        counts = {}

        #A,B = "AGCTTCTTTCGAGTACAAAAAC","TCCCAGATTATATCTATCACTGA"
        A,B = "GTACAAAAAC","TCCCAGATTA"
        lenA = len(A)
        cnt,nreads,recognized = 0,0,0
        for line in open(fastq_file):
            cnt += 1
            if cnt%4==2:
                nreads += 1
                if (nreads%1000000==0): sys.stderr.write("reads=%s, recognized barcodes=%s (%0.1f%%)\n" % (nreads,recognized,100.*recognized/float(nreads)))
                seq = line.rstrip()
                a = seq.find(A)
                if a==-1: continue
                b = seq.find(B)
                if b==-1: continue
                sz = b-(a+lenA)
                if sz<10 or sz>30: continue
                barcode = seq[a+lenA:b] # these are reverse-complements, but rc(barcodes) stored in hash too
                if barcode not in barcodemap: continue
                id = barcodemap[barcode]
                if id not in counts: counts[id] = 0
                counts[id] += 1
                recognized += 1

        for id in IDs:
            vals = [id,counts.get(id,0)]
            print('\t'.join([str(x) for x in vals]))
            print("\n")
       

    def create_combined_counts(self,headers, counts_list):
        import pandas as pd
        df_list =[]
        for f in counts_list:
            sys.stderr.write("Adding in file # %s \n"%f)
            counts_df = pd.read_csv(f, sep="\t")
            counts_df["sgRNA"]=counts_df[counts_df.columns[0]].str.split("_v", expand=True)[0]
            counts_df = counts_df.drop(columns=[counts_df.columns[0]])
            counts_df.set_index("sgRNA",inplace=True)
            df_list.append(counts_df)
        combined_df = pd.concat(df_list, axis=1)
        combined_df.columns = headers
        combined_df_text = combined_df.to_csv(sep="\t")
        sys.stderr.write("Number of sgRNAs in combined counts file (present in all counts files): %d \n"%len(combined_df))
        print(combined_df_text)


    def extract_abund(self,combined_counts_file,metadata_file,control_condition,sgRNA_efficiency_file,no_dep_abund,drug,days,PC=1e-8):  
        import pandas as pd
        
        metadata = pd.read_csv(metadata_file, sep="\t")
        metadata = metadata[((metadata["drug"]==drug) | (metadata["drug"]==control_condition)) & (metadata["days_predepletion"]==int(days))]
        if(len(metadata)==0):
            sys.stderr.write("This combination of conditions does not exist in your metadata file. Please select one that does")
            sys.exit(0)
        elif (drug not in metadata["drug"].values.tolist()):
            sys.stderr.write("%s is not found in your metadata. Add the drug's information in the metadata file or select a different drug"%drug)
            sys.exit(0)
        elif (int(days) not in metadata["days_predepletion"].values.tolist()):
            sys.stderr.write("%d is not found in your metadata days of predepletion column. Add the day's information in the metadata file or select a different day"%days)
            sys.exit(0)
        elif (control_condition not in metadata["drug"].values.tolist()):
            sys.stderr.write("%s is not found in your metadata. Add the corresponding information in the metadata file or select a different control"%control_condition)
            sys.exit(0)
        metadata = metadata.sort_values(by=["conc_xMIC"])
        column_names = metadata["column_name"].values.tolist()
        concs_list = metadata["conc_xMIC"].values.tolist()
        
        print("# Condition Tested : "+drug+" D"+days)
        headers = []
        combined_counts_df = pd.read_csv(combined_counts_file,sep="\t", index_col=0)
        combined_counts_df = combined_counts_df[column_names]

        if(len(combined_counts_df.columns)==0):
            sys.stderr.write("The samples assocaited with the selected drugs do not exist in your combined counts file. Please select one that does and check your metadata file has corresponding column names")
            sys.exit(0)
        elif(len(combined_counts_df.columns)<len(metadata)):
            sys.stderr.write("WARNING: Not all of the samples from the metadata based on this criteron have a column in the combined counts file")
      
        sgRNA_efficiency = pd.read_csv(sgRNA_efficiency_file,sep="\t", index_col=0)
        sgRNA_efficiency = sgRNA_efficiency.iloc[:,-1:]
        sgRNA_efficiency.columns = ["sgRNA efficiency"]
        sgRNA_efficiency["sgRNA"] = sgRNA_efficiency.index
        sgRNA_efficiency["sgRNA"]=sgRNA_efficiency["sgRNA"].str.split("_v", expand=True)[0]
        sgRNA_efficiency.set_index("sgRNA",inplace=True)

        no_dep_df = pd.read_csv(no_dep_abund, sep="\t", index_col=0, header=None)
        no_dep_df = no_dep_df.iloc[:,-1:]
        no_dep_df.columns = ["uninduced ATC values"] 
        no_dep_df["uninduced ATC values"] = no_dep_df["uninduced ATC values"]/ no_dep_df["uninduced ATC values"].sum()
        no_dep_df["sgRNA"] = no_dep_df.index
        no_dep_df["sgRNA"]=no_dep_df["sgRNA"].str.split("_v", expand=True)[0]
        no_dep_df.set_index("sgRNA",inplace=True)

        abund_df = pd.concat([sgRNA_efficiency, no_dep_df,combined_counts_df], axis=1)
        abund_df= abund_df[~(abund_df.index.str.contains("Negative") | abund_df.index.str.contains("Empty"))]
        sys.stderr.write("Disregarding Empty or Negative sgRNAs\n")
        sys.stderr.write("%d sgRNAs are all of the following files : sgRNA efficiency metadata, uninduced ATC counts file, combined counts file\n"%len(abund_df))

        headers = ["sgRNA efficiency","uninduced ATC values"]
        for i,col in enumerate(column_names):
            abund_df[col] = abund_df[col]/abund_df[col].sum()
            abund_df[col] = (abund_df[col]+PC)/(abund_df["uninduced ATC values"]+PC)
            headers.append(str(concs_list[i])+"_"+str(i))
            print("# "+str(concs_list[i])+" conc_xMIC"+" - "+col)

        abund_df.columns = headers
        abund_df["sgRNA"] = abund_df.index.values.tolist()
        abund_df[["orf-gene","remaining"]] = abund_df["sgRNA"].str.split('_',n=1,expand=True)
        abund_df[["orf","gene"]]= abund_df["orf-gene"].str.split(':',expand=True)
        abund_df = abund_df.drop(columns=["orf-gene","remaining"])
        abund_df = abund_df.dropna()
        
        abund_df.insert(0, "sgRNA efficiency", abund_df.pop("sgRNA efficiency"))
        abund_df.insert(0, "uninduced ATC values", abund_df.pop("uninduced ATC values"))
        abund_df.insert(0, 'gene', abund_df.pop('gene'))
        abund_df.insert(0, 'orf', abund_df.pop('orf'))
        abund_df.insert(0, 'sgRNA', abund_df.pop('sgRNA'))

        abund_df_text = abund_df.to_csv(sep="\t", index=False)
        print(abund_df_text)


  #####################################################

  # derived from logsigmoidfit.R
  # see heatmap.py for example of how to put data in a pandas.DataFrame and call an R function like make_heatmapFunc()

    def run_model(self, frac_abund_file):
        import pandas as pd
        import numpy as np
        from mne.stats import fdr_correction
        import statsmodels.api as sm
        
        frac_abund_df = pd.read_csv(frac_abund_file, sep="\t",comment='#')

        drug_output = []
        for i,gene in enumerate(set(frac_abund_df["gene"])):
            #print(i,gene)
            sys.stderr.write("Analyzing Gene # %d \n"%i)
            gene_df = frac_abund_df[frac_abund_df["gene"]==gene]
            orf = gene_df["orf"].iloc[0]
            gene_df = gene_df.drop(columns=["orf","gene","uninduced ATC values"])

            melted_df = gene_df.melt(id_vars=["sgRNA","sgRNA efficiency"],var_name="conc",value_name="abund")
            melted_df["conc"] = melted_df["conc"].str.split("_", expand=True)[0].astype(float)
            min_conc = min(melted_df[melted_df["conc"]>0]["conc"])
            melted_df.loc[melted_df["conc"]==0,"conc"] = min_conc/2
            melted_df["abund"] = [0.01+(1-0.01)*(1-np.exp(-2*float(i)))/(1+np.exp(-2*float(i))) for i in melted_df["abund"]]
            melted_df["logsig abund"] = [np.nan if (1-x)== 0 else np.log10(float(x)/(1-float(x))) for x in melted_df["abund"]]
            melted_df["log conc"] = [np.log2(float(x)) for x in melted_df["conc"]]
            

            melted_df = melted_df.dropna()
            if len(melted_df.index)<2:
                drug_output.append([orf,gene,len(gene_df)]+[np.nan]*6)
                continue
            
            Y = melted_df["logsig abund"]
            X = melted_df.drop(columns=["abund", "logsig abund", "sgRNA", "conc"])
            X = sm.add_constant(X)
            model = sm.OLS(Y,X)
            results = model.fit()
            coeffs = results.params
            pvals = results.pvalues
            drug_output.append([orf,gene,len(gene_df)]+coeffs.values.tolist()+pvals.values.tolist())
            sys.stderr.flush()

        drug_out_df = pd.DataFrame(drug_output, columns=["Orf","Gene","Nobs", "intercept","coefficient sgRNA_efficiency","coefficient concentration dependence","pval intercept","pval sgRNA_efficiency","pval concentration dependence"])
        drug_out_df["intercept"] = round(drug_out_df["intercept"],6)
        drug_out_df["coefficient sgRNA_efficiency"] = round(drug_out_df["coefficient sgRNA_efficiency"],6)
        drug_out_df["coefficient concentration dependence"] = round(drug_out_df["coefficient concentration dependence"],6)
        drug_out_df["pval intercept"] = round(drug_out_df["pval intercept"],6)
        drug_out_df["pval sgRNA_efficiency"] = round(drug_out_df["pval sgRNA_efficiency"],6)
        drug_out_df["pval concentration dependence"] = round(drug_out_df["pval concentration dependence"],6)


        mask = np.isfinite(drug_out_df["pval concentration dependence"])
        pval_corrected = np.full(drug_out_df["pval concentration dependence"].shape, np.nan)
        pval_corrected[mask] = fdr_correction(drug_out_df["pval concentration dependence"][mask])[1]
        drug_out_df["qval concentration dependence"] = pval_corrected
        drug_out_df["qval concentration dependence"] = round(drug_out_df["qval concentration dependence"] ,6)
        drug_out_df = drug_out_df.replace(np.nan,1)

        drug_out_df["Z score of concentration dependence"] = (drug_out_df["coefficient concentration dependence"] - drug_out_df["coefficient concentration dependence"].mean())/drug_out_df["coefficient concentration dependence"].std()
        drug_out_df["Z score of concentration dependence"] = round(drug_out_df["Z score of concentration dependence"], 6)
        drug_out_df["Significant Interactions"] = [0] * len(drug_out_df)
        drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) & (drug_out_df["Z score of concentration dependence"]<-2),"Significant Interactions"]=-1
        drug_out_df.loc[(drug_out_df["qval concentration dependence"]<0.05) & (drug_out_df["Z score of concentration dependence"]>2),"Significant Interactions"]=1
        drug_out_df.insert(0, "Significant Interactions", drug_out_df.pop("Significant Interactions"))

        n = len(drug_out_df[drug_out_df["Significant Interactions"]!=0])
        depl_n = len(drug_out_df[drug_out_df["Significant Interactions"]== -1])
        enrich_n = len(drug_out_df[drug_out_df["Significant Interactions"]==1])
        sys.stderr.write("%d Total Significant Gene Interactions\n"%n)
        sys.stderr.write("%d Significant Gene Depletions\n"%depl_n)
        sys.stderr.write("%d Significant Gene Enrichments\n"%enrich_n)
    
        drug_out_df  = drug_out_df.replace(r'\s+',np.nan,regex=True).replace('',np.nan)
        drug_out_txt = drug_out_df.to_csv(sep="\t", index=False)
        print("# Total Significant Gene Interactions : ", n)
        print("# Significant Gene Depletions : ", depl_n)
        print("# Significant Gene Enrichments : ", enrich_n)
        print(drug_out_txt)

    def visualize(self,fractional_abundances_file, gene, fig_location):
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import statsmodels.api as sm

        abund_df = pd.read_csv(fractional_abundances_file,sep="\t", comment="#")
        with open(fractional_abundances_file) as f:
            first_line = f.readline()
            condition = first_line.split(" : ")[1]

        abund_df = abund_df[(abund_df["gene"]==gene)| (abund_df["orf"]==gene)]
        if len(abund_df)==0:
            sys.stderr.write("Gene not found : %d \n"%idx)
            sys.exit(0)
        abund_df = abund_df.reset_index(drop=True)
        all_slopes = []

        df_list = []
        for idx,row in abund_df.iterrows():
            sys.stderr.write("Fitting sgRNA # : %d \n"%idx)
            raw_Y= row[5:].values
            Y = [max(0.01,x) for x in raw_Y]
            Y = [np.log10(x) for x in Y]

            X = abund_df.columns[5:]
            X = [float(i.split("_")[0]) for i in X]
            min_conc = min([i for i in X if i>0])
            X = [min_conc/2 if i==0 else i for i in X ]
            X = [np.log2(float(x)) for x in X]

            data = pd.DataFrame({"Log (Concentration)":X, "Log (Relative Abundance)":Y})
            X = pd.DataFrame({"log concentration":X})
            X_in = sm.add_constant(X, has_constant='add')
            results = sm.OLS(Y,X_in).fit()
            all_slopes.append(results.params[1])
            data["sgRNA efficiency"] = [row["sgRNA efficiency"]] * len(data)
            data["slope"] = [results.params[1]] * len(data)
            df_list.append(data)


        plot_df = pd.concat(df_list)
        plt.figure()
        cmap =  mpl.colors.LinearSegmentedColormap.from_list("", ["#8ecae6","#219ebc","#023047","#ffb703","#fb8500"], N=len(abund_df))
        palette = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
        #print("-----------", bo_palette.as_hex())
        g = sns.lmplot(data=plot_df, x='Log (Concentration)', y='Log (Relative Abundance)', hue="sgRNA efficiency", palette=palette, legend=False,ci=None, scatter=False, line_kws={"lw":0.75})

        sm1 = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=plot_df['sgRNA efficiency'].min(), vmax=0, clip=False), cmap=cmap)
        g.figure.colorbar(sm1, shrink=0.8, aspect=50, label="sgRNA efficiency")
        g.set(ylim=(-2.5, 1.0))
        plt.gca().set_title(gene+"\n"+condition, wrap=True)
        plt.tight_layout()
        plt.savefig(fig_location)
################################

if __name__ == "__main__":

    G = CGI.fromargs(sys.argv[1:])
    G.Run()


