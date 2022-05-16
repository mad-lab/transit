import scipy
import numpy
import heapq
import math
import statsmodels.stats.multitest

import time
import sys
import collections

from pytransit.analysis import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools

############# GUI ELEMENTS ##################

short_name = "anova"
long_name = "anova"
short_desc = "Perform Anova analysis"
long_desc = """Perform Anova analysis"""
EOL = "\n"

transposons = ["", ""]
columns = []

class AnovaAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, AnovaMethod)
def main():
    print("ANOVA example")

class AnovaMethod(base.MultiConditionMethod):
    """
    anova
    """
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, excluded_conditions=[], included_conditions=[], nterm=0.0, cterm=0.0, PC=1, winz=False, refs=[]):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file,
                normalization=normalization, excluded_conditions=excluded_conditions, included_conditions=included_conditions, nterm=nterm, cterm=cterm)
        self.PC = PC
        self.refs = refs
        self.winz = winz

    @classmethod
    def transit_error(self,msg): print("error: %s" % msg) # for some reason, transit_error() in base class or transit_tools doesn't work right; needs @classmethod

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(AnovaMethod.usage_string())
            sys.exit(0)

        combined_wig = args[0]
        annotation = args[2]
        metadata = args[1]
        output_file = args[3]
        normalization = kwargs.get("n", "TTR")
        NTerminus = float(kwargs.get("iN", 0.0))
        CTerminus = float(kwargs.get("iC", 0.0))
        winz = True if "winz" in kwargs else False
        PC = int(kwargs.get("PC", 5))
        refs = kwargs.get("-ref",[]) # list of condition names to use a reference for calculating LFCs
        if refs!=[]: refs = refs.split(',')
        excluded_conditions = list(filter(None, kwargs.get("-exclude-conditions", "").split(",")))
        included_conditions = list(filter(None, kwargs.get("-include-conditions", "").split(",")))

        # check for unrecognized flags
        flags = "-n --exclude-conditions --include-conditions -iN -iC -PC --ref -winz".split()
        for arg in rawargs:
          if arg[0]=='-' and arg not in flags:
            self.transit_error("flag unrecognized: %s" % arg)
            print(AnovaMethod.usage_string())
            sys.exit(0)

        return self(combined_wig, metadata, annotation, normalization, output_file, excluded_conditions, included_conditions, NTerminus, CTerminus, PC, winz, refs)

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, self.unknown_cond_flag) for f in filenamesInCombWig]

    def means_by_condition_for_gene(self, sites, conditions, data):
        """
            Returns a dictionary of {Condition: Mean} for each condition.
            ([Site], [Condition]) -> {Condition: Number}
            Site :: Number
            Condition :: String
        """
        nTASites = len(sites)
        wigsByConditions = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
            wigsByConditions[c].append(i)

        return { c: numpy.mean(self.winsorize(data[wigIndex][:, sites])) if nTASites > 0 else 0 for (c, wigIndex) in wigsByConditions.items() }

    def means_by_rv(self, data, RvSiteindexesMap, genes, conditions):
        """
            Returns Dictionary of mean values by condition
            ([[Wigdata]], {Rv: SiteIndex}, [Gene], [Condition]) -> {Rv: {Condition: Number}}
            Wigdata :: [Number]
            SiteIndex :: Number
            Gene :: {start, end, rv, gene, strand}
            Condition :: String
        """
        MeansByRv = {}
        for gene in genes:
            Rv = gene["rv"]
            MeansByRv[Rv] = self.means_by_condition_for_gene(RvSiteindexesMap[Rv], conditions, data)
        return MeansByRv

    def group_by_condition(self, wigList, conditions):
        """
            Returns array of datasets, where each dataset corresponds to one condition.
            ([[Wigdata]], [Condition]) -> [[DataForCondition]]
            Wigdata :: [Number]
            Condition :: String
            DataForCondition :: [Number]
        """
        countsByCondition = collections.defaultdict(lambda: [])
        countSum = 0
        for i, c in enumerate(conditions):
          countSum += numpy.sum(wigList[i])
          countsByCondition[c].append(wigList[i])

        return (countSum, [numpy.array(v).flatten() for v in countsByCondition.values()])

    # since this is in both ZINB and ANOVA, should move to stat_tools.py?

    def winsorize(self, counts):
      # input is insertion counts for gene: list of lists: n_replicates (rows) X n_TA sites (cols) in gene
      unique_counts = numpy.unique(numpy.concatenate(counts))
      if (len(unique_counts) < 2): return counts
      else:
        n, n_minus_1 = unique_counts[heapq.nlargest(2, range(len(unique_counts)), unique_counts.take)]
        result = [[ n_minus_1 if count == n else count for count in wig] for wig in counts]
        return numpy.array(result)

    def run_anova(self, data, genes, MeansByRv, RvSiteindexesMap, conditions):
        """
            Runs Anova (grouping data by condition) and returns p and q values
            ([[Wigdata]], [Gene], {Rv: {Condition: Mean}}, {Rv: [SiteIndex]}, [Condition]) -> Tuple([Number], [Number])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            Mean :: Number
            SiteIndex: Integer
            Condition :: String
        """
        count = 0
        self.progress_range(len(genes))

        pvals,Rvs,status = [],[],[]
        for gene in genes:
            count += 1
            Rv = gene["rv"]
            if (len(RvSiteindexesMap[Rv]) <= 1):
                status.append("TA sites <= 1")
                pvals.append(1)
            else:
                countSum, countsVec = self.group_by_condition(list(map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data)), conditions)
                if self.winz: countsVec = self.winsorize(countsVec)

                if (countSum == 0):
                    pval = 1
                    status.append("No counts in all conditions")
                    pvals.append(pval)
                else:
                    stat,pval = scipy.stats.f_oneway(*countsVec)
                    status.append("-")
                    pvals.append(pval)
            Rvs.append(Rv)

            # Update progress
            text = "Running Anova Method... %5.1f%%" % (100.0*count/len(genes))
            self.progress_update(text, count)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[1] # BH, alpha=0.05

        p,q,statusMap = {},{},{}
        for i,rv in enumerate(Rvs):
            p[rv],q[rv],statusMap[rv] = pvals[i],qvals[i],status[i]
        return (p, q, statusMap)

    def calcLFCs(self,means,refs=[],PC=1):
      if len(refs)==0: refs = means # if ref condition(s) not explicitly defined, use mean of all
      grandmean = numpy.mean(refs)
      lfcs = [math.log((x+PC)/float(grandmean+PC),2) for x in means]
      return lfcs

    def Run(self):
        self.transit_message("Starting Anova analysis")
        start_time = time.time()

        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        self.transit_message("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)
        if self.winz: self.transit_message("Winsorizing insertion counts")

        conditionsByFile, _, _, orderingMetadata = tnseq_tools.read_samples_metadata(self.metadata)
        conditions = self.wigs_to_conditions(
            conditionsByFile,
            filenamesInCombWig)

        conditionsList = self.select_conditions(conditions,self.included_conditions,self.excluded_conditions,orderingMetadata)

        conditionNames = [conditionsByFile[f] for f in filenamesInCombWig]
        fileNames = filenamesInCombWig
        
        data, fileNames, conditionNames, conditions, _, _ = self.filter_wigs_by_conditions3( # in base.py
                data,
                fileNames, # it looks like fileNames and conditionNames have to be parallel to data (vector of wigs)
                conditionNames, # original Condition column in samples metadata file
                self.included_conditions,
                self.excluded_conditions,
                conditions = conditionNames) # this is kind of redundant for ANOVA, but it is here because condition, covars, and interactions could have been manipulated for ZINB

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap, nterm=self.NTerminus, cterm=self.CTerminus)
        MeansByRv = self.means_by_rv(data, RvSiteindexesMap, genes, conditions)

        self.transit_message("Running Anova")
        pvals,qvals,run_status = self.run_anova(data, genes, MeansByRv, RvSiteindexesMap, conditions)

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")

        heads = ("Rv Gene TAs".split() +
                ["Mean_%s" % x for x in conditionsList] +
                ["LFC_%s" % x for x in conditionsList] +
                "pval padj".split() + ["status"])
        file.write("#Console: python3 %s\n" % " ".join(sys.argv))
        file.write("#parameters: normalization=%s, trimming=%s/%s%% (N/C), pseudocounts=%s\n" % (self.normalization,self.NTerminus,self.CTerminus,self.PC))
        file.write('#'+'\t'.join(heads)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            if Rv in MeansByRv:
              means = [MeansByRv[Rv][c] for c in conditionsList]
              refs = [MeansByRv[Rv][c] for c in self.refs]
              LFCs = self.calcLFCs(means,refs,self.PC)
              vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                      ["%0.2f" % x for x in means] + 
                      ["%0.3f" % x for x in LFCs] + 
                      ["%f" % x for x in [pvals[Rv], qvals[Rv]]] + [run_status[Rv]])
              file.write('\t'.join(vals)+EOL)
        file.close()
        self.transit_message("Finished Anova analysis")
        self.transit_message("Time: %0.1fs\n" % (time.time() - start_time))

    @classmethod
    def usage_string(self):
        usage =  """Usage: python3 transit.py anova <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
 Optional Arguments:
  -n <string>         :=  Normalization method. Default: -n TTR
  --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
  --exclude-conditions <cond1,...> := Comma-separated list of conditions to exclude (Default: none)
  --ref <cond> := which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)
  -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
  -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
  -PC <N> := pseudocounts to use for calculating LFC. Default: -PC 5
  -winz   := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)"""
        return usage

if __name__ == "__main__":
    main()

