import scipy
import numpy
import heapq
import statsmodels.stats.multitest

import time
import sys
import collections

from rpy2.robjects import r, globalenv, IntVector, FloatVector, StrVector, packages as rpackages

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
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions=[], included_conditions=[], winz=False):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file, normalization=normalization)
        self.ignored_conditions = ignored_conditions
        self.included_conditions = included_conditions
        self.unknown_cond_flag = "FLAG-UNMAPPED-CONDITION-IN-WIG"
        self.winz = winz

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
        winz = True if kwargs.has_key("w") else False
        ignored_conditions = filter(None, kwargs.get("-ignore-conditions", "").split(","))
        included_conditions = filter(None, kwargs.get("-include-conditions", "").split(","))
        if len(included_conditions) > 0 and len(ignored_conditions) > 0:
            print("Cannot use both include-conditions and ignore-conditions flags")
            print(ZinbMethod.usage_string())
            sys.exit(0)

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions, included_conditions, winz)

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, self.unknown_cond_flag) for f in filenamesInCombWig]

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
            if (len(siteIndexes) == 0): # If no TA sites, write 0
                means[c] = 0
                nz_means[c] = 0
                nz_percs[c] = 0
            else:
                arr = data[wigIndex][:, siteIndexes]
                means[c] = numpy.mean(arr) if len(arr) > 0 else 0
                nonzero_arr = nonzero(arr)
                nz_means[c] = numpy.mean(nonzero_arr) if len(nonzero_arr) > 0 else 0
                nz_percs[c] = nzperc(arr)
        return [means, nz_means, nz_percs]

    def filter_wigs_by_conditions(self, data, conditions, ignored_conditions, included_conditions):
        """
            Filters conditions that are ignored/included.
            ([[Wigdata]], [Condition], [Condition], [Condition]) -> Tuple([[Wigdata]], [Condition])
        """
        ignored_conditions, included_conditions = (set(ignored_conditions), set(included_conditions))
        d_filtered, cond_filtered = [], [];
        if len(ignored_conditions) > 0 and len(included_conditions) > 0:
            self.transit_error("Both ignored and included conditions have len > 0", ignored_conditions, included_conditions)
            sys.exit(0)
        elif (len(ignored_conditions) > 0):
            self.transit_message("conditions ignored: {0}".format(ignored_conditions))
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag) and (c not in ignored_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
        elif (len(included_conditions) > 0):
            self.transit_message("conditions included: {0}".format(included_conditions))
            d_filtered, cond_filtered = [], [];
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag) and (c in included_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
        else:
            for i, c in enumerate(conditions):
              if (c != self.unknown_cond_flag):
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
            zinb_signif = function(count, condition, NZmean, logitZPerc, sat_adjust=TRUE) {
              suppressMessages(require(pscl))
              suppressMessages(require(MASS))
              melted = data.frame(cnt=count, cond=condition, NZmean = NZmean, logitZperc = logitZPerc)
              sums = aggregate(melted$cnt,by=list(melted$cond),FUN=sum)
              # to avoid model failing due to singular condition, add fake counts of 1 to all conds if any cond is all 0s
              if (0 %in% sums[,2]) {
                # print("adding pseudocounts")
                for (i in 1:length(sums[,1])) {
                  subset = melted[melted$cond==sums[i,1],]
                  newvec = subset[1,]
                  newvec$cnt = 1
                  melted = rbind(melted,newvec) }
              }
              status = "-"
              minCount = min(melted$cnt)
              mod1 = tryCatch(
                { if (minCount == 0) {
                  if (sat_adjust) {
                    zeroinfl(cnt~0+cond+offset(log(NZmean))|0+cond+offset(logitZperc),data=melted,dist="negbin")
                  } else {
                    zeroinfl(cnt~0+cond+offset(log(NZmean)),data=melted,dist="negbin")
                  }
                } else {
                  glm.nb(cnt~0+cond,data=melted)
                }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              mod0 = tryCatch( # null model, independent of conditions
                { if (minCount == 0) {
                  if (sat_adjust) { zeroinfl(cnt~1+offset(log(NZmean))|1+offset(logitZperc),data=melted,dist="negbin") }
                  else { zeroinfl(cnt~1,data=melted,dist="negbin") }
                } else {
                  glm.nb(cnt~1,data=melted)
                }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              if (is.null(mod1) | is.null(mod0)) { return (c(1, paste0("Model Error. ", status))) }
              if ((minCount == 0) && (sum(is.na(coef(summary(mod1))$count[,4]))>0)) { return(c(1, "Has Coefs, pvals are NAs")) } # rare failure mode - has coefs, but pvals are NA
              df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
              pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
              # this gives same answer, but I would need to extract the Pvalue...
              #require(lmtest)
              #print(lrtest(mod1,mod0))
              return (c(pval, "-"))
            }
        ''')

        return globalenv['zinb_signif']

    def winsorize(self, data):
        unique_counts = numpy.unique(numpy.concatenate(data))
        if (len(unique_counts) < 2):
            return data
        else:
            n, n_minus_1 = unique_counts[heapq.nlargest(2, xrange(len(unique_counts)), unique_counts.take)]
            result = [[ n_minus_1 if count == n else count
                        for count in wig] for wig in data]
        return numpy.array(result)

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
        pvals,Rvs, status = [],[], []
        r_zinb_signif = self.def_r_zinb_signif()
        if (self.winz):
            self.transit_message("Winsorizing and running analysis...")

        for gene in genes:
            count += 1
            # Rv = gene["rv"]
            Rv = "Rv3915"
            if (len(RvSiteindexesMap[Rv]) <= 1):
                status.append("TA sites <= 1")
                pvals.append(1)
            else:
                ## TODO :: Option for sat adjustment?
                norm_data = self.winsorize(map(
                    lambda wigData: wigData[RvSiteindexesMap[Rv]], data)) if self.winz else map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data)
                ([ readCounts,
                   condition,
                   NZmean,
                   logitZPerc]) = self.melt_data(
                           norm_data,
                           conditions, NZMeanByRep, LogZPercByRep)
                if (numpy.sum(readCounts) == 0):
                    status.append("No counts in all conditions")
                    pvals.append(1)
                else:
                    r_args = [IntVector(readCounts), StrVector(condition), FloatVector(NZmean), FloatVector(logitZPerc)] + [True]
                    pval, msg = r_zinb_signif(*r_args)
                    status.append(msg)
                    pvals.append(float(pval))
            Rvs.append(Rv)
            # Update progress
            text = "Running ZINB Method... %5.1f%%" % (100.0*count/len(genes))
            self.progress_update(text, count)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals)[1] # BH, alpha=0.05

        p,q,statusMap = {},{},{}
        for i,rv in enumerate(Rvs):
            p[rv],q[rv],statusMap[rv] = pvals[i],qvals[i],status[i]
        return (p, q, statusMap)

    def Run(self):
        self.transit_message("Starting ZINB analysis")
        start_time = time.time()
        packnames = ("MASS", "pscl")
        r_packages_needed = [x for x in packnames if not rpackages.isinstalled(x)]
        if (len(r_packages_needed) > 0):
            self.transit_error(
                    "Error: Following R packages are required: %(0)s. From R console, You can install them using install.packages(c(%(0)s))"
                    % ({'0': '"{0}"'.format('", "'.join(r_packages_needed))}))
            sys.exit(1)


        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        self.transit_message("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        conditions = self.wigs_to_conditions(
            self.read_samples_metadata(self.metadata),
            filenamesInCombWig)
        data, conditions = self.filter_wigs_by_conditions(data, conditions, self.ignored_conditions, self.included_conditions)

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap)
        [MeansByRv, NzMeansByRv, NzPercByRv] = self.stats_by_rv(data, RvSiteindexesMap, genes, conditions)
        LogZPercByRep, NZMeanByRep = self.global_stats_for_rep(data)

        self.transit_message("Running ZINB")
        pvals,qvals,run_status = self.run_zinb(data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions)

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")
        conditionsList = self.included_conditions if len(self.included_conditions) > 0 else list(set(conditions))
        head = ("Rv Gene TAs".split() +
                map(lambda v: "mean_" + v, conditionsList) +
                map(lambda v: "NZmean_" + v, conditionsList) +
                map(lambda v: "NZperc_" + v, conditionsList) +
                "pval padj".split() + ["status"])

        file.write('\t'.join(head)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                    ["%0.1f" % MeansByRv[Rv][c] for c in conditionsList] +
                    ["%0.1f" % NzMeansByRv[Rv][c] for c in conditionsList] +
                    ["%0.1f" % NzPercByRv[Rv][c] for c in conditionsList] +
                    ["%f" % x for x in [pvals[Rv], qvals[Rv]]]) + [run_status[Rv]]
            file.write('\t'.join(vals)+EOL)
        file.close()
        self.transit_message("Finished Zinb analysis")
        self.transit_message("Time: %0.1fs\n" % (time.time() - start_time))

    @classmethod
    def usage_string(self):
        return """python %s zinb <combined wig file> <annotation .prot_table> <samples_metadata file> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        -w                  :=  If set, Winsorize data.
        --ignore-conditions <cond1,cond2> :=  Comma seperated list of conditions to ignore, for the analysis.
        --include-conditions <cond1,cond2> :=  Comma seperated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.

        """ % (sys.argv[0])

if __name__ == "__main__":
    main()

