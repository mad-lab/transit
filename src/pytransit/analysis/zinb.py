import scipy
import numpy
import statsmodels.stats.multitest

import time
import sys
import collections

from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects import r, globalenv

import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools

############# GUI ELEMENTS ##################

short_name = "ZINB"
long_name = "ZINB"
short_desc = "Perform ZINB analysis"
long_desc = """Perform ZINB analysis"""
EOL = "\n"

transposons = ["", ""]
columns = []

class ZinbAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, ZinbMethod)
def main():
    print("ZINB example")

class ZinbMethod(base.MultiConditionMethod):
    """
    Zinb
    """
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions=set()):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file, normalization=normalization)
        self.ignored_conditions = ignored_conditions

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(ZinbMethod.usage_string())
            sys.exit(0)

        combined_wig = args[0]
        annotation = args[1]
        metadata = args[2]
        output_file = args[3]
        normalization = kwargs.get("n", "TTR")
        ignored_conditions = set(kwargs.get("-ignore-conditions", "Unknown").split(","))

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions)

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, "Unknown") for f in filenamesInCombWig]

    def stats_by_condition_for_gene(self, siteIndexes, conditions, data):
        """
            Returns a dictionary of [{Condition: Stat}] for Mean, NzMean, NzPercentage
            ([SiteIndex], [Condition], [WigData]) -> [{Condition: Number}]
            SiteIndex :: Number
            WigData :: [Number]
            Condition :: String
        """
        wigsByConditions = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
            wigsByConditions[c].append(i)

        nonzero = lambda xs: xs[numpy.nonzero(xs)]
        nzperc = lambda xs: numpy.count_nonzero(xs)/float(xs.size)
        zperc = lambda xs: 1 - numpy.count_nonzero(xs)/float(xs.size)

        means = {}
        nz_means = {}
        nz_percs = {}
        for (c, wigIndex) in wigsByConditions.items():
            arr = data[wigIndex][:, siteIndexes]
            means[c] = numpy.mean(arr) if len(arr) > 0 else 0
            nonzero_arr = nonzero(arr)
            nz_means[c] = numpy.mean(nonzero_arr) if len(nonzero_arr) > 0 else 0

            nz_percs[c] = nzperc(arr)
        return [means, nz_means, nz_percs]

    def filter_by_conditions_blacklist(self, data, conditions, ignored_conditions):
        """
            Filters out wigfiles, with ignored conditions.
            ([[Wigdata]], [Condition]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered = [], [];
        for i, c in enumerate(conditions):
          if c not in ignored_conditions:
            d_filtered.append(data[i])
            cond_filtered.append(conditions[i])

        return (numpy.array(d_filtered), numpy.array(cond_filtered))

    def read_samples_metadata(self, metadata_file):
        """
          Filename -> ConditionMap
          ConditionMap :: {Filename: Condition}
        """
        wigFiles = []
        conditionsByFile = {}
        headersToRead = ["condition", "filename"]
        with open(metadata_file) as mfile:
            lines = mfile.readlines()
            headIndexes = [i
                    for h in headersToRead
                    for i, c in enumerate(lines[0].split())
                    if c.lower() == h]
            for line in lines:
                if line[0]=='#': continue
                vals = line.split()
                [condition, wfile] = vals[headIndexes[0]], vals[headIndexes[1]]
                conditionsByFile[wfile] = condition
        return conditionsByFile

    def stats_by_rv(self, data, RvSiteindexesMap, genes, conditions):
        """
            Returns Dictionary of Stats by condition for each Rv
            ([[Wigdata]], {Rv: SiteIndex}, [Gene], [Condition]) -> {Rv: {Condition: Number}}
            Wigdata :: [Number]
            SiteIndex :: Number
            Gene :: {start, end, rv, gene, strand}
            Condition :: String
        """
        MeansByRv = {}
        NzMeansByRv = {}
        NzPercByRv = {}
        for gene in genes:
            Rv = gene["rv"]

            if len(RvSiteindexesMap[gene["rv"]]) > 1: # skip genes with no or 1 TA site
                (MeansByRv[Rv],
                 NzMeansByRv[Rv],
                 NzPercByRv[Rv]) = self.stats_by_condition_for_gene(RvSiteindexesMap[Rv], conditions, data)
        return [MeansByRv, NzMeansByRv, NzPercByRv]

    def global_stats_for_rep(self, data):
        """
        Returns the logit zero percentage and nz_mean for each replicate.
            [[WigData]] -> [[Number] ,[Number]]
        """

        logit_zero_perc = []
        nz_mean = []
        for wig in data:
            zero_perc = (wig.size - numpy.count_nonzero(wig))/float(wig.size)
            logit_zero_perc.append(numpy.log(zero_perc/(1 - zero_perc)))
            nz_mean.append(numpy.mean(wig[numpy.nonzero(wig)]))
        return [numpy.array(logit_zero_perc), numpy.array(nz_mean)]

    def group_by_condition(self, wigList, conditions):
        """
            Returns array of datasets, where each dataset corresponds to one condition.
            ([[Wigdata]], [Condition]) -> [[DataForCondition]]
            Wigdata :: [Number]
            Condition :: String
            DataForCondition :: [Number]
        """
        countsByCondition = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
          countsByCondition[c].append(wigList[i])

        return [numpy.array(v).flatten() for v in countsByCondition.values()]

    def melt_data(self, readCountsForRv, conditions, NZMeanByRep, LogZPercByRep):
        rvSitesLength = len(readCountsForRv[0])
        repeatAndFlatten = lambda xs: numpy.repeat(xs, rvSitesLength)
        return [
                numpy.concatenate(readCountsForRv).astype(int),
                repeatAndFlatten(conditions),
                repeatAndFlatten(NZMeanByRep),
                repeatAndFlatten(LogZPercByRep)
               ]

    def def_r_zinb_signif(self):
        r('''
            zinb_signif = function(count, condition, NZmean, logitZPerc, sat_adjust=TRUE)
                {
                  suppressMessages(require(pscl))
                  mod1 = tryCatch(
                     {
                         if (sat_adjust) {
                             zeroinfl(count~0+condition+offset(log(NZmean))|0+condition+offset(logitZPerc), dist="negbin")
                         } else {
                             zeroinfl(count~0+condition, dist="negbin")
                         }
                     }, error=function(err) { return(NULL) } )
                  mod0 = tryCatch( # null model, independent of conditions
                     {
                         if (sat_adjust) {
                             zeroinfl(count~1+offset(log(NZmean))|1+offset(logitZPerc), dist="negbin")
                         } else {
                             zeroinfl(count~1, dist="negbin")
                         }
                     }, error=function(err) { return(NULL) } )
                  if (is.null(mod1) | is.null(mod0)) { return(1) }
                  if (sum(is.na(coef(summary(mod1))$count[,4]))>0) { return(1) } # rare failure mode - has coefs, but pvals are NA
                  df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
                  pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()

                  return(pval)
                }
        ''')

        return globalenv['zinb_signif']

    def run_zinb(self, data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions):
        """
            Runs Zinb for each gene across conditions and returns p and q values
            ([[Wigdata]], [Gene], [[Number]], {Rv: [SiteIndex]}, [Condition]) -> Tuple([Number], [Number])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            SiteIndex: Integer
            Condition :: String
        """
        count = 0
        self.progress_range(len(genes))
        pvals,Rvs = [],[]
        r_zinb_signif = self.def_r_zinb_signif()

        for gene in genes:
            count += 1
            Rv = gene["rv"]
            ## TODO :: Option for sat adjustment?
            ([ readCounts,
               condition,
               NZmean,
               logitZPerc]) = self.melt_data(
                       map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data),
                       conditions, NZMeanByRep, LogZPercByRep)
            r_args = map(numpy2ri, [readCounts, condition, NZmean, logitZPerc]) + [True]
            pval = r_zinb_signif(*r_args)[0]
            pvals.append(pval)
            Rvs.append(Rv)
            # Update progress
            text = "Running ZINB Method... %5.1f%%" % (100.0*count/len(genes))
            self.progress_update(text, count)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals)[1] # BH, alpha=0.05

        p,q = {},{}
        for i,rv in enumerate(Rvs):
            p[rv],q[rv] = pvals[i],qvals[i]
        return (p, q)

    def Run(self):
        self.transit_message("Starting ZINB analysis")
        start_time = time.time()

        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        self.transit_message("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        conditions = self.wigs_to_conditions(
            self.read_samples_metadata(self.metadata),
            filenamesInCombWig)
        data, conditions = self.filter_by_conditions_blacklist(data, conditions, self.ignored_conditions)

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap)
        [MeansByRv, NzMeansByRv, NzPercByRv] = self.stats_by_rv(data, RvSiteindexesMap, genes, conditions)
        LogZPercByRep, NZMeanByRep = self.global_stats_for_rep(data)

        self.transit_message("Running ZINB")
        pvals,qvals = self.run_zinb(data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions)

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")
        conditionsList = list(set(conditions))
        vals = "Rv Gene TAs".split() + map(lambda v: "mean_" + v, conditionsList) + map(lambda v: "NZmean_" + v, conditionsList) + map(lambda v: "NZperc_" + v, conditionsList) + "pval padj".split()

        file.write('\t'.join(vals)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            if Rv in MeansByRv:
              vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                      ["%0.1f" % MeansByRv[Rv][c] for c in conditionsList] +
                      ["%0.1f" % NzMeansByRv[Rv][c] for c in conditionsList] +
                      ["%0.1f" % NzPercByRv[Rv][c] for c in conditionsList] +
                      ["%f" % x for x in [pvals[Rv], qvals[Rv]]])
              file.write('\t'.join(vals)+EOL)
        file.close()
        self.transit_message("Finished Zinb analysis")

    @classmethod
    def usage_string(self):
        return """python %s zinb <combined wig file> <annotation .prot_table> <samples_metadata file> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        --ignore-conditions <cond1,cond2> :=  Comma seperated list of conditions to ignore, for the analysis. Default --ignore-conditions Unknown

        """ % (sys.argv[0])

if __name__ == "__main__":
    main()

