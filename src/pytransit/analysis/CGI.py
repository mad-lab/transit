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
        return """usage (3 sub-commands):
    python3 ../src/transit.py CGI extract_counts <fastq_file> <ids_file> > <counts_file>
    python3 ../src/transit.py CGI create_combined_counts <comma seperated headers> <counts_file_1> <counts_file_2> ... <counts_file_n> > <combined_counts_file>
    python3 ../src/transit.py CGI extract_abund <combined_counts_file> <metadata_file> <reference_condition> <sgRNA_strength_file> <no_depletion_abundances_file> <drug> <days>  >  <frac_abund_file>
    python3 ../src/transit.py CGI run_model <abund_file>  >  <logsigmodfit_file>
    python3 ../src/transit.py CGI post_process <logsigmoidfit_file>  >  <results_file>
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
        cmd,args,kwargs = self.cmd,self.args,self.kwargs

        if cmd=="extract_counts":
           fastq_file = args[0]
           ids_file = args[1]
           self.extract_counts(fastq_file, ids_file)
        
        elif cmd=="create_combined_counts":
           headers = args[0].split(",")
           counts_file_list = args[1:]
           self.create_combined_counts(headers,counts_file_list)

        elif cmd=="extract_abund":
          if len(args)<6: print(self.usage_string())
          combined_counts_file = args[0]
          metadata_file = args[1]
          reference_condition=args[2]
          extrapolated_LFCs_file = args[3]
          no_dep_abund = args[4]
          drug = args[5]
          days = args[6]
          self.extract_abund(combined_counts_file,metadata_file,reference_condition,extrapolated_LFCs_file,no_dep_abund,drug,days)
        elif cmd == "run_model":
          ifile_path = args[0] #example frac_abund_RIF_D5.txt
          self.run_model(ifile_path)
        elif cmd == "post_process":
            logsig_file_path = args[0]
            self.post_process(logsig_file_path)
        else: print(self.usage_string())

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
       

    def create_combined_counts(self,headers, counts_list):
        import pandas as pd
        df_list =[]
        for f in counts_list:
            counts_df = pd.read_csv(f, sep="\t")
            counts_df["sgRNA"]=counts_df[counts_df.columns[0]].str.split("_v", expand=True)[0]
            counts_df = counts_df.drop(columns=[counts_df.columns[0]])
            counts_df.set_index("sgRNA",inplace=True)
            df_list.append(counts_df)
        combined_df = pd.concat(df_list, axis=1)
        combined_df.columns = headers
        combined_df_text = combined_df.to_csv(sep="\t")
        print(combined_df_text)


    def extract_abund(self,combined_counts_file,metadata_file,reference_condition,extrapolated_LFCs_file,no_dep_abund,drug,days,PC=1e-8):  
      import pandas as pd

      metadata = pd.read_csv(metadata_file, sep="\t")
      metadata = metadata[((metadata["drug"]==drug) | (metadata["drug"]==reference_condition)) & (metadata["days_predepletion"]==int(days))]
      metadata = metadata.sort_values(by=["conc_xMIC"])
      column_names = metadata["column_name"].values.tolist()
      concs_list = metadata["conc_xMIC"].values.tolist()
      concs={}
      i=1
      for c in metadata["conc_xMIC"].values.tolist():
          if c not in concs: 
              concs[c] = i
              i=2*i
      headers = []
      combined_counts_df = pd.read_csv(combined_counts_file,sep="\t", index_col=0)
      combined_counts_df = combined_counts_df[column_names]
      
      extrapolated_LFCs = pd.read_csv(extrapolated_LFCs_file,sep="\t", index_col=0)
      extrapolated_LFCs = extrapolated_LFCs.iloc[:,-1:]
      extrapolated_LFCs.columns = ["extrapolated LFCs"]
      extrapolated_LFCs["sgRNA"] = extrapolated_LFCs.index
      extrapolated_LFCs["sgRNA"]=extrapolated_LFCs["sgRNA"].str.split("_v", expand=True)[0]
      extrapolated_LFCs.set_index("sgRNA",inplace=True)

      no_dep_df = pd.read_csv(no_dep_abund, sep="\t", index_col=0, header=None)
      no_dep_df = no_dep_df.iloc[:,-1:]
      no_dep_df.columns = ["-ATC values"]
      no_dep_df["sgRNA"] = no_dep_df.index
      no_dep_df["sgRNA"]=no_dep_df["sgRNA"].str.split("_v", expand=True)[0]
      no_dep_df.set_index("sgRNA",inplace=True)

      abund_df = pd.concat([extrapolated_LFCs, no_dep_df,combined_counts_df], axis=1)
      abund_df= abund_df[~(abund_df.index.str.contains("Negative") | abund_df.index.str.contains("Empty"))]
      headers = ["extrapolated LFCs","-ATC values"]
      for i,col in enumerate(column_names):
         abund_df[col] = abund_df[col]/abund_df[col].sum()
         abund_df[col] = (abund_df[col]+PC)/(abund_df["-ATC values"]+PC)
         headers.append(concs[concs_list[i]])
         print("# "+str(concs[concs_list[i]])+" conc_xMIC"+" - "+col)

      abund_df.columns = headers
      abund_df["sgRNA"] = abund_df.index.values.tolist()
      abund_df[["orf-gene","reminaing"]] = abund_df["sgRNA"].str.split('_',n=1,expand=True)
      abund_df[["orf","gene"]]= abund_df["orf-gene"].str.split(':',expand=True)
      abund_df = abund_df.drop(columns=["orf-gene","reminaing","sgRNA"])
      abund_df = abund_df.dropna()
      
      abund_df.insert(0, "extrapolated LFCs", abund_df.pop("extrapolated LFCs"))
      abund_df.insert(0, "-ATC values", abund_df.pop("-ATC values"))
      abund_df.insert(0, 'gene', abund_df.pop('gene'))
      abund_df.insert(0, 'orf', abund_df.pop('orf'))

      abund_df_text = abund_df.to_csv(sep="\t")
      print(abund_df_text)
      #sys.stderr.write("warning: extrapolated_LFCs values not found for %s gRNAs\n" % a)
      #sys.stderr.write("warning: no_dep values not found for %s gRNAs\n" % b)

  #####################################################

  # derived from logsigmoidfit.R
  # see heatmap.py for example of how to put data in a pandas.DataFrame and call an R function like make_heatmapFunc()

    def run_model(self, frac_abund_file):
        import pandas as pd
        import numpy as np
        from mne.stats import fdr_correction
        import statsmodels.api as sm
        pd.set_option('display.max_columns', 500)
        frac_abund_df = pd.read_csv(frac_abund_file, sep="\t",comment='#')

        drug_output = []
        for i,gene in enumerate(set(frac_abund_df["gene"])):
            #print(i,gene)
            sys.stderr.write("Analyzing Gene # %d \n"%i)
            gene_df = frac_abund_df[frac_abund_df["gene"]==gene]
            orf = gene_df["orf"].iloc[0]
            gene_df = gene_df.drop(columns=["orf","gene","-ATC values"])

            melted_df = gene_df.melt(id_vars=["sgRNA","extrapolated LFCs"],var_name="conc",value_name="abund")
            melted_df["conc"] = melted_df["conc"].str.split(".",expand=True)[0]
            melted_df["abund"] = [0.01+(1-0.01)*(1-np.exp(-2*float(i)))/(1+np.exp(-2*float(i))) for i in melted_df["abund"]]
            melted_df["logsig abund"] = [np.nan if (1-x)== 0 else np.log10(float(x)/(1-float(x))) for x in melted_df["abund"]]
            melted_df["log conc"] = [np.log2(int(x)) for x in melted_df["conc"]]

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

        drug_out_df = pd.DataFrame(drug_output, columns=["Orf","Gene","Nobs", "coef intercept","coef extrapolated_LFCs","coef conc","pval intercept","pval pred_logFC","pval conc"])
    
        mask = np.isfinite(drug_out_df["pval conc"])
        pval_corrected = np.full(drug_out_df["pval conc"].shape, np.nan)
        pval_corrected[mask] = fdr_correction(drug_out_df["pval conc"][mask])[1]
        drug_out_df["qval conc"] = pval_corrected
        drug_out_df = drug_out_df.replace(np.nan,1)

        drug_out_df["Z"] = (drug_out_df["coef conc"] - drug_out_df["coef conc"].mean())/drug_out_df["coef conc"].std()

        n = len(drug_out_df[(drug_out_df["qval conc"]<0.05) & ((drug_out_df["Z"]<-2) | (drug_out_df["Z"]>2))])
        sys.stderr.write("%d Signifincant Genes"%n)
    
        drug_out_df  = drug_out_df.replace(r'\s+',np.nan,regex=True).replace('',np.nan)
        drug_out_txt = drug_out_df.to_csv(sep="\t")
        print(drug_out_txt)



################################

if __name__ == "__main__":

    G = CGI.fromargs(sys.argv[1:])
    G.Run()


