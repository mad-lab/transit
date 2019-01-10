import time
import sys
import heapq
import numpy as np

import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools

short_name = "winsorize"
long_name = "winsorize"
short_desc = "Perform Winsorization"
long_desc = """Perform Winsorization"""
EOL = "\n"

transposons = ["", ""]
columns = []

class Winsorize(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, WinsorizeMethod)

class WinsorizeMethod(base.MultiConditionMethod):
    """
    winsorize
    """

    def __init__(self, combined_wig, annotation, output_file):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, "", annotation, output_file)

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(WinsorizeMethod.usage_string())
            sys.exit(0)

        combined_wig = args[0]
        annotation = args[1]
        output_file = args[2]

        return self(combined_wig, annotation, output_file)

    def winsorize_data(self, data):
        unique_counts = np.unique(np.concatenate(data))
        if (len(unique_counts) < 2):
            return data
        else:
            n, n_minus_1 = unique_counts[heapq.nlargest(2, xrange(len(unique_counts)), unique_counts.take)]
            result = [[ n_minus_1 if count == n else count
                        for count in wig] for wig in data]
        return np.array(result)

    def Run(self):
        self.transit_message("Starting Winsorization")
        start_time = time.time()

        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        genes = tnseq_tools.read_genes(self.annotation_path)
        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap)

        self.transit_message("Running Winsorize")

        count = 0
        result = np.array([[] for _ in range(len(filenamesInCombWig))])
        self.progress_range(len(genes))
        for gene in genes:
            count += 1
            Rv = gene["rv"]
            # Update progress
            winz_data = self.winsorize_data(map(lambda xs: xs[RvSiteindexesMap[Rv]], data))
            result = np.hstack((result, winz_data))
            text = "Running Winsorize Method... %5.1f%%" % (100.0*count/len(genes))
            self.progress_update(text, count)
        print("******")
        print(result)
        norm_data = np.array(result)

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")
        file.write("#Converted to CombinedWig with TRANSIT.\n")
        file.write("#Console: python %s\n" % " ".join(sys.argv))
        for f in filenamesInCombWig:
            file.write("#File: %s\n" % f)
        for i in range(len(sites)):
            file.write('\t'.join(
                [str(sites[i])] +
                ["%0.1f" % x for x in list(norm_data[...,i])])+"\n")
        file.close()
        self.transit_message("Finished Winsorization")
        self.transit_message("Time: %0.1fs\n" % (time.time() - start_time))

    @classmethod
    def usage_string(self):
        return """python %s winsorize <combined wig file> <annotation .prot_table> <output file>
        """ % (sys.argv[0])


def main():
    print("Winsorize example")
