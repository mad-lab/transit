import scipy
import numpy
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
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions=[], included_conditions=[], nterm=0.0, cterm=0.0, PC=1):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file,
                normalization=normalization, ignored_conditions=ignored_conditions, included_conditions=included_conditions, nterm=nterm, cterm=cterm)
        self.PC = PC

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
        PC = int(kwargs.get("PC", 5))
        ignored_conditions = list(filter(None, kwargs.get("-ignore-conditions", "").split(",")))
        included_conditions = list(filter(None, kwargs.get("-include-conditions", "").split(",")))

        # check for unrecognized flags
        flags = "-n --ignore-conditions --include-conditions -iN -iC -PC".split()
        for arg in rawargs:
          if arg[0]=='-' and arg not in flags:
            self.transit_error("flag unrecognized: %s" % arg)
            print(AnovaMethod.usage_string())
            sys.exit(0)

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions, included_conditions, NTerminus, CTerminus, PC)

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

        return { c: numpy.mean(data[wigIndex][:, sites]) if nTASites > 0 else 0 for (c, wigIndex) in wigsByConditions.items() }

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

    def calcLFCs(self,means,PC=1):
      grandmean = numpy.mean(means)
      lfcs = [math.log((x+PC)/float(grandmean+PC),2) for x in means]
      return lfcs

    def Run(self):
        self.transit_message("Starting Anova analysis")
        start_time = time.time()

        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        self.transit_message("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        conditionsByFile, _, _, orderingMetadata = tnseq_tools.read_samples_metadata(self.metadata)
        conditions = self.wigs_to_conditions(
            conditionsByFile,
            filenamesInCombWig)

        conditionsList = self.select_conditions(conditions,self.included_conditions,self.ignored_conditions,orderingMetadata)
        data, conditions, _, _ = self.filter_wigs_by_conditions2(data, conditions, conditionsList)

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
        file.write('\t'.join(heads)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            if Rv in MeansByRv:
              means = [MeansByRv[Rv][c] for c in conditionsList]
              LFCs = self.calcLFCs(means,self.PC)
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
        usage =  """Usage: python33 transit.py anova <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
 Optional Arguments:
  -n <string>         :=  Normalization method. Default: -n TTR
  --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
  --ignore-conditions <cond1,...> := Comma-separated list of conditions to ignore (Default: none)
  -iN <int> :=  Ignore TAs within given percentage of N terminus. Default: -iN 0
  -iC <int> :=  Ignore TAs within given percentage of C terminus. Default: -iC 0
  -PC <int> := pseudocounts to use for calculating LFC. Default: -PC 5"""
        return usage

#        return """python3 %s anova <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
#
#        Optional Arguments:
#        -n <string>         :=  Normalization method. Default: -n TTR
#        --ignore-conditions <cond1,cond2> :=  Comma separated list of conditions to ignore, for the analysis. Default --ignore-conditions Unknown
#        --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.
#        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
#        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0
#        """ % (sys.argv[0])

if __name__ == "__main__":
    main()

