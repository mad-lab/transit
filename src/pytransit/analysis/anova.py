import scipy
import numpy
import statsmodels.stats.multitest

import time
import sys
import collections

import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.norm_tools as norm_tools
import pytransit.utils as utils

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
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file, normalization=normalization, ignored_conditions=ignored_conditions)

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)
        combined_wig = kwargs['-wig']
        metadata = kwargs['-meta']
        annotation = kwargs['-prot']
        normalization = kwargs.get("n", "TTR")
        ignored_conditions = set(kwargs.get("-ignore-conditions", " ").split(","))
        output_file = kwargs['o']

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions)

    def Run(self):
        self.transit_message("Starting Anova analysis")
        print(self.ignored_conditions)
        start_time = time.time()

        def read_samples_metadata(metadata_file):
            conditions = set()
            wigFiles = []
            conditionsByFile = {}
            with open(metadata_file) as mfile:
                next(mfile) # Skip headers
                for line in mfile:
                    if line[0]=='#': continue
                    [id, condition, wfile] = line.split()
                    conditionsByFile[wfile] = condition
            return conditionsByFile

        def read_combined_wig(fname):
          sites,countsByWig,files = [],[],[]
          with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("#File:"):
                  files.append(line.split()[1])
            countsByWig = [[] for _ in files]
            for line in lines:
                if line[0]=='#': continue
                cols = map(int, line.split("\t")[0:-1])
                position, wigCounts = cols[0], cols[1:]
                sites.append(position)
                for i, c in enumerate(wigCounts):
                    countsByWig[i].append(c)

          return (numpy.array(sites), numpy.array(countsByWig), files)

        def wigsToConditions(conditionsByFile, filenamesInCombWig):
            return [conditionsByFile.get(f, "Unknown") for f in filenamesInCombWig]

        def meansByConditionForGene(sites, conditions, data):
          wigsByConditions = collections.defaultdict(lambda: [])
          for i, c in enumerate(conditions):
            wigsByConditions[c].append(i)

          return { c: numpy.mean(data[wigIndex][:, sites]) for (c, wigIndex) in wigsByConditions.items() }

        def filterByConditionBlackList(data, conditions, ignored_conditions):
          d_filtered, cond_filtered = [], [];
          for i, c in enumerate(conditions):
            if c not in ignored_conditions:
              d_filtered.append(data[i])
              cond_filtered.append(conditions[i])

          return (numpy.array(d_filtered), numpy.array(cond_filtered))

        def groupByCondition(wigList, conditions):
          countsByCondition = collections.defaultdict(lambda: [])
          for i, c in enumerate(conditions):
            countsByCondition[c].append(wigList[i])

          return [numpy.array(v).flatten() for v in countsByCondition.values()]

        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = read_combined_wig(self.combined_wig)
        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            (data, factors) = norm_tools.normalize_data(data, self.normalization)

        conditionsByFile = read_samples_metadata(self.metadata)
        conditions = wigsToConditions(conditionsByFile, filenamesInCombWig)
        data, conditions = filterByConditionBlackList(data, conditions, self.ignored_conditions)

        genes = utils.read_genes(self.annotation_path)

        siteshash = {}
        for i,TA in enumerate(sites): siteshash[TA] = i

        TAsites = {}
        for g,(start,end,Rv,gene,strand) in enumerate(genes):
            siteindexes = []
            ## Coords off by one?
            for i in range(start,end): # end+1?
              co = i+1
              if co in siteshash: siteindexes.append(siteshash[co])
            TAsites[Rv] = siteindexes

        Means = {}
        for (start,end,Rv,gene,strand) in genes:
          if len(TAsites[Rv]) > 0: # skip genes with no TA sites
            Means[Rv] = meansByConditionForGene(TAsites[Rv], conditions, data)

        pvals,Rvs, = [],[] # lists
        for (start,end,Rv,gene,strand) in genes:
          if Rv in Means:
            countsvec = groupByCondition(map(lambda wigData: wigData[TAsites[Rv]], data), conditions)
            stat,pval = scipy.stats.f_oneway(*countsvec)
            pvals.append(pval)
            Rvs.append(Rv)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape,numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[1] # BH, alpha=0.05

        p,q = {},{}
        for i,rv in enumerate(Rvs):
           p[rv],q[rv] = pvals[i],qvals[i]
        pvals,qvals = p,q, # hashes on Rv

        file = open(self.output,"w")
        conditionsList = list(set(conditions))
        vals = "Rv Gene TAs".split() + conditionsList + "pval padj".split()
        file.write('\t'.join(vals)+EOL)
        for (start,end,Rv,gene,strand) in genes:
          if Rv in Means:
            vals = [Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % Means[Rv][c] for c in conditionsList]+["%f" % x for x in [pvals[Rv],qvals[Rv]]]
            file.write('\t'.join(vals)+EOL)
        file.close()

if __name__ == "__main__":
    main()

