import scipy
import numpy
import heapq
import statsmodels.stats.multitest

import time
import sys
import collections

hasR = False
try:
    import rpy2.robjects
    hasR = True
except Exception as e:
    hasR = False

if hasR:
    from rpy2.robjects import r, DataFrame, globalenv, IntVector, FloatVector, StrVector, packages as rpackages

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
DEBUG = False
GENE = None

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
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions=[], included_conditions=[], winz=False, nterm=5.0, cterm=5.0, condition="Condition", covars=[], interactions = []):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file,
                normalization=normalization, ignored_conditions=ignored_conditions, included_conditions=included_conditions, nterm=nterm, cterm=cterm)
        self.winz = winz
        self.covars = covars
        self.interactions = interactions
        self.condition = condition

    @classmethod
    def fromargs(self, rawargs):
        if not hasR:
            print("Error: R and rpy2 (< 2.9.0) required to run ZINB analysis.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2<2.9.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(ZinbMethod.usage_string())
            sys.exit(0)

        if (kwargs.get('v', False)):
            global DEBUG
            DEBUG = True

        if (kwargs.get('-gene', False)):
            global GENE
            GENE = kwargs.get('-gene', None)

        combined_wig = args[0]
        metadata = args[1]
        annotation = args[2]
        output_file = args[3]
        normalization = kwargs.get("n", "TTR")
        NTerminus = float(kwargs.get("iN", 5.0))
        CTerminus = float(kwargs.get("iC", 5.0))
        condition = kwargs.get("-condition", "Condition")
        covars = filter(None, kwargs.get("-covars", "").split(","))
        interactions = filter(None, kwargs.get("-interactions", "").split(","))
        winz = True if kwargs.has_key("w") else False
        ignored_conditions = filter(None, kwargs.get("-ignore-conditions", "").split(","))
        included_conditions = filter(None, kwargs.get("-include-conditions", "").split(","))
        if len(included_conditions) > 0 and len(ignored_conditions) > 0:
            self.transit_error("Cannot use both include-conditions and ignore-conditions flags")
            print(ZinbMethod.usage_string())
            sys.exit(0)

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions, included_conditions, winz, NTerminus, CTerminus, condition, covars, interactions)

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, self.unknown_cond_flag) for f in filenamesInCombWig]

    def wigs_to_covariates(self, covariatesMap, filenamesInCombWig):
        """
            Returns list of covariate lists. Each covariate list consists of covariates corresponding to given wigfiles.
            ([{FileName: Covar}], [FileName]) -> [[Covar]]
            Covar :: String
        """
        try:
            return [[covarsByFile[f]
                        for f in filenamesInCombWig]
                        for covarsByFile in covariatesMap]
        except KeyError:
            self.transit_error("Error: Covariates not found for file {0}".format(f))
            sys.exit(0)

    def wigs_to_interactions(self, interactionsMap, filenamesInCombWig):
        """
            Returns list of interaction lists. Each interaction list consists of covariates corresponding to given wigfiles.
            ([{FileName: Interaction}], [FileName]) -> [[Interaction]]
            Interaction :: String
        """
        try:
            return [[covarsByFile[f]
                        for f in filenamesInCombWig]
                        for covarsByFile in interactionsMap]
        except KeyError:
            self.transit_error("Error: Interaction var not found for file {0}".format(f))
            sys.exit(0)

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

    def melt_data(self, readCountsForRv, conditions, covariates, interactions, NZMeanByRep, LogZPercByRep):
        rvSitesLength = len(readCountsForRv[0])
        repeatAndFlatten = lambda xs: numpy.repeat(xs, rvSitesLength)

        return [
                numpy.concatenate(readCountsForRv).astype(int),
                repeatAndFlatten(conditions),
                map(repeatAndFlatten, covariates),
                map(repeatAndFlatten, interactions),
                repeatAndFlatten(NZMeanByRep),
                repeatAndFlatten(LogZPercByRep)
               ]

    def def_r_zinb_signif(self):
        r('''
            zinb_signif = function(df,
                zinbMod1,
                zinbMod0,
                nbMod1,
                nbMod0, DEBUG = F) {
              suppressMessages(require(pscl))
              suppressMessages(require(MASS))
              melted = df

              # filter out genes that have low saturation across all conditions, since pscl sometimes does not fit params well (resulting in large negative intercepts and high std errors)
              NZpercs = aggregate(melted$cnt,by=list(melted$cond),FUN=function(x) { sum(x>0)/length(x) })
              if (max(NZpercs$x)<=0.15) { return(c(pval=1,status="pan-essential, not analyzed")) }

              sums = aggregate(melted$cnt,by=list(melted$cond),FUN=sum)
              # to avoid model failing due to singular condition, add fake counts of 1 to all conds if any cond is all 0s
              if (0 %in% sums[,2]) {
                # print("adding pseudocounts")
                for (i in 1:length(sums[,1])) {
                  subset = melted[melted$cond==sums[i,1],]
                  newvec = subset[1,]
                  newvec$cnt = 1 # note: NZmean and NZperc are copied from last dataset in condition
                  #newvec$cnt = as.integer(mean(subset$cnt))+1 # add the mean for each condition as a pseudocount
                  melted = rbind(melted,newvec) }
              }
              status = "-"
              minCount = min(melted$cnt)
              f1 = ""
              mod1 = tryCatch(
                {
                  if (minCount == 0) {
                    f1 = zinbMod1
                    mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin")
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod1' }
                    mod
                  } else {
                    f1 = nbMod1
                    glm.nb(as.formula(nbMod1),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              f0 = ""
              mod0 = tryCatch( # null model, independent of conditions
                {
                  if (minCount == 0) {
                    f0 = zinbMod0
                    mod = zeroinfl(as.formula(zinbMod0),data=melted,dist="negbin")
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod0' }
                    mod
                  } else {
                    f0 = nbMod0
                    glm.nb(as.formula(nbMod0),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              if (DEBUG) {
                  print("Model 1:")
                  print(f1)
                  print(summary(mod1))
                  print("Model 0:")
                  print(f0)
                  print(summary(mod0))
              }

              if (is.null(mod1) | is.null(mod0)) { return (c(1, paste0("Model Error. ", status))) }
              if ((minCount == 0) && (sum(is.na(coef(summary(mod1))$count[,4]))>0)) { return(c(1, "Has Coefs, pvals are NAs")) } # rare failure mode - has coefs, but pvals are NA
              df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
              pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
              # this gives same answer, but I would need to extract the Pvalue...
              #require(lmtest)
              #print(lrtest(mod1,mod0))
              return (c(pval, status))
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

    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def run_zinb(self, data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions, covariates, interactions):
        """
            Runs Zinb for each gene across conditions and returns p and q values
            ([[Wigdata]], [Gene], [Number], [Number], {Rv: [SiteIndex]}, [Condition], [Covar], [Interaction]) -> Tuple([Number], [Number], [Status])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            SiteIndex: Integer
            Condition :: String
            Covar :: String
            Interaction :: String
            Status :: String
        """

        count = 0
        self.progress_range(len(genes))
        pvals,Rvs, status = [],[], []
        r_zinb_signif = self.def_r_zinb_signif()
        if (self.winz):
            self.transit_message("Winsorizing and running analysis...")

        self.transit_message("Condition: %s" % self.condition)
        covarsFormula = ""
        if (len(self.covars) > 0):
            covarsFormula = "+".join([""] + self.covars)
            self.transit_message("Covariates: {0}".format(covarsFormula))

        condFormula = "cond"
        if (len(self.interactions) > 0):
            condFormula = "+".join(["cond*{0}".format(interaction) for interaction in self.interactions])
            self.transit_message("Interactions: {0}".format(condFormula))

        sat_adjust = True # just for testing False, but we always want to leave it as True
        zinbMod1 = ("cnt~1+{cond}+offset(log(NZmean)){covars}|1+{cond}+offset(logitZperc){covars}" if sat_adjust else "cnt~1+{cond}{covars}").format(cond=condFormula, covars=covarsFormula)
        zinbMod0 = ("cnt~1+offset(log(NZmean)){covars}|1+offset(logitZperc){covars}" if sat_adjust else "cnt~1{covars}").format(covars=covarsFormula)
        nbMod1 = "cnt~1+{cond}{covars}".format(cond=condFormula, covars=covarsFormula)
        nbMod0 = "cnt~1{covars}".format(covars=covarsFormula)


        # proposed new way to do it
        comp1a = "1+cond"
        comp1b = "1+cond"
        # include cond in mod0 only if testing interactions
        comp0a = "1" if len(self.interactions)==0 else "1+cond"
        comp0b = "1" if len(self.interactions)==0 else "1+cond"
        for I in self.interactions: comp1a += "*"+I; comp1b += "*"+I; comp0a += "+"+I; comp0b += "+"+I
        for C in self.covars: comp1a += "+"+C; comp1b += "+"+C; comp0a += "+"+C; comp0b += "+"+C
        zinbMod1 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp1a,comp1b)
        zinbMod0 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp0a,comp0b)

        nbMod1 = "cnt~1+{cond}{covars}".format(cond=condFormula, covars=covarsFormula)
        nbMod0 = "cnt~1{covars}".format(covars=covarsFormula)
        toRFloatOrStrVec = lambda xs: FloatVector(xs) if self.is_number(xs[0]) else StrVector(xs)

        for gene in genes:
            count += 1
            Rv = gene["rv"]
            ## Single gene case for debugging
            if (GENE):
                Rv = None
                if RvSiteindexesMap.has_key(GENE):
                    Rv = GENE
                else:
                    for g in genes:
                        if (g['gene'] == GENE):
                            Rv = g["rv"]
                            break
                if not Rv:
                    self.transit_error("Cannot find gene: {0}".format(GENE))
                    sys.exit(0)

            if (len(RvSiteindexesMap[Rv]) <= 1):
                status.append("TA sites <= 1")
                pvals.append(1)
            else:
                norm_data = self.winsorize(map(
                    lambda wigData: wigData[RvSiteindexesMap[Rv]], data)) if self.winz else map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data)
                ([ readCounts,
                   condition,
                   covarsData,
                   interactionsData,
                   NZmean,
                   logitZPerc]) = self.melt_data(
                           norm_data,
                           conditions, covariates, interactions, NZMeanByRep, LogZPercByRep)
                if (numpy.sum(readCounts) == 0):
                    status.append("No counts in all conditions")
                    pvals.append(1)
                else:
                    df_args = {
                        'cnt': IntVector(readCounts),
                        'cond': toRFloatOrStrVec(condition),
                        'NZmean': FloatVector(NZmean),
                        'logitZperc': FloatVector(logitZPerc)
                        }
                    ## Add columns for covariates and interactions if they exist.
                    df_args.update(map(lambda (i, c): (c, toRFloatOrStrVec(covarsData[i])), enumerate(self.covars)))
                    df_args.update(map(lambda (i, c): (c, toRFloatOrStrVec(interactionsData[i])), enumerate(self.interactions)))

                    melted = DataFrame(df_args)
                    # r_args = [IntVector(readCounts), StrVector(condition), melted, map(lambda x: StrVector(x), covars), FloatVector(NZmean), FloatVector(logitZPerc)] + [True]
                    debugFlag = True if DEBUG or GENE else False
                    pval, msg = r_zinb_signif(melted, zinbMod1, zinbMod0, nbMod1, nbMod0, debugFlag)
                    status.append(msg)
                    pvals.append(float(pval))
                if (DEBUG or GENE):
                    self.transit_message("Pval for Gene {0}: {1}, status: {2}".format(Rv, pvals[-1], status[-1]))
                if (GENE):
                    self.transit_message("Ran for single gene. Exiting...")
                    sys.exit(0)
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

        condition_name = self.condition
        conditionsByFile, covariatesByFileList, interactionsByFileList = tnseq_tools.read_samples_metadata(self.metadata, self.covars, self.interactions, condition_name=condition_name)

        ## [Condition] in the order of files in combined wig
        conditions = self.wigs_to_conditions(
            conditionsByFile,
            filenamesInCombWig)
        ## [Covariate] in the order of files in combined wig
        covariates = self.wigs_to_covariates(
            covariatesByFileList,
            filenamesInCombWig)
        ## [Interaction] in the order of files in combined wig
        interactions = self.wigs_to_interactions(
            interactionsByFileList,
            filenamesInCombWig)
        data, conditions, covariates, interactions = self.filter_wigs_by_conditions(
                data,
                conditions,
                covariates = covariates,
                interactions = interactions,
                ignored_conditions = self.ignored_conditions,
                included_conditions = self.included_conditions)

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap, nterm=self.NTerminus, cterm=self.CTerminus)
        [MeansByRv, NzMeansByRv, NzPercByRv] = self.stats_by_rv(data, RvSiteindexesMap, genes, conditions)
        LogZPercByRep, NZMeanByRep = self.global_stats_for_rep(data)

        self.transit_message("Running ZINB")
        pvals, qvals, run_status = self.run_zinb(data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions, covariates, interactions)

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")
        conditionsList = self.included_conditions if len(self.included_conditions) > 0 else list(set(conditions))
        head = ("Rv Gene TAs".split() +
                map(lambda v: "mean_" + v, conditionsList) +
                map(lambda v: "NZmean_" + v, conditionsList) +
                map(lambda v: "NZperc_" + v, conditionsList) +
                "pval padj".split() + ["status"])

        file.write("#Console: python %s\n" % " ".join(sys.argv))
        file.write('\t'.join(head)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                    ["%0.2f" % MeansByRv[Rv][c] for c in conditionsList] +
                    ["%0.2f" % NzMeansByRv[Rv][c] for c in conditionsList] +
                    ["%0.2f" % NzPercByRv[Rv][c] for c in conditionsList] +
                    ["%f" % x for x in [pvals[Rv], qvals[Rv]]]) + [run_status[Rv]]
            file.write('\t'.join(vals)+EOL)
        file.close()
        self.transit_message("Finished Zinb analysis")
        self.transit_message("Time: %0.1fs\n" % (time.time() - start_time))

    @classmethod
    def usage_string(self):
        return """python %s zinb <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        -w                  :=  If set, Winsorize data.
        --ignore-conditions <cond1,cond2> :=  Comma separated list of conditions to ignore, for the analysis.
        --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.
        -iN <float>     :=  Ignore TAs occuring within given percentage of the N terminus. Default: -iN 5.0
        -iC <float>     :=  Ignore TAs occuring within given percentage of the C terminus. Default: -iC 5.0
        --covars <covar1,...>     :=  Comma separated list of covariates (in metadata file) to include, for the analysis.
        --interactions <covar1,...>     :=  Comma separated list of covariates to include, that interact with the condition for the analysis. Must be factors
        """ % (sys.argv[0])

if __name__ == "__main__":
    main()

