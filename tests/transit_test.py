import unittest
import os
import sys

basedir = os.path.dirname(__file__)
ctrl_rep1 = basedir + "/../src/pytransit/data/glycerol_H37Rv_rep1.wig"
ctrl_rep2 = basedir + "/../src/pytransit/data/glycerol_H37Rv_rep2.wig"
ctrl_data_txt = ",".join([ctrl_rep1, ctrl_rep2])
mini_wig = basedir + "/data/test.wig"

combined_wig = basedir + "/../src/pytransit/data/cholesterol_glycerol_combined.dat"
samples_metadata = basedir + "/../src/pytransit/data/samples_metadata_cg.txt"
samples_metadata_covariates = basedir + "/../src/pytransit/data/samples_metadata_cg_covar.txt"
samples_metadata_interactions = basedir + "/../src/pytransit/data/samples_metadata_cg_interactions.txt"

exp_rep1 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep1.wig"
exp_rep2 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep2.wig"
exp_rep3 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep3.wig"
exp_data_txt = ",".join([exp_rep1, exp_rep2, exp_rep3])

all_data_list = [ctrl_rep1, ctrl_rep2, exp_rep1, exp_rep2, exp_rep3]

annotation = basedir + "/../src/pytransit/genomes/H37Rv.prot_table"
small_annotation = basedir + "/data/test.prot_table"
output = basedir + "/testoutput.txt"
hist_path = output.rsplit(".", 1)[0] + "_histograms"
tpp_output_base = basedir + "/test_tpp_temp"
tpp_output_paths = [tpp_output_base + i for i in [".counts", ".reads1", ".sam", ".tn_stats", ".trimmed1", ".trimmed1_failed_trim", ".wig", "_a.counts", "_b.counts", "_c.counts", "_1.counts", "_2.counts", "_3.counts"]]

# For tpp
reads1 = basedir + "/data/test.fastq"
test_multicontig = basedir + "/data/test-multicontig.fna"
test_multicontig_reads1 = basedir + "/data/test-multicontig-1.fastq"
test_multicontig_reads2 = basedir + "/data/test-multicontig-2.fastq"
h37fna = basedir + "/../src/pytransit/genomes/H37Rv.fna"


class TransitTestCase(unittest.TestCase):

    def setUp(self):
        # Print header
        self.header()

    def tearDown(self):
        for f in tpp_output_paths:
            if os.path.exists(f):
                print("Removing tpp test file")
                os.remove(f)

        if os.path.exists(hist_path):
            print("Removing histogram files")
            for f in os.listdir(hist_path):
                os.remove(os.path.join(hist_path, f))
            os.rmdir(hist_path)

        # Check if there were output files and remove them
        if os.path.exists(output):
            print("Removing output file...")
            os.remove(output)

        genes_path = output.rsplit(".", 1)[0] + "_genes" + output.rsplit(".", 1)[1]

        if os.path.exists(genes_path):
            print("Removing genes file...")
            os.remove(genes_path)

    def header(self):
        print("\n")
        print("#"*20)
        print(self.id())
        print("#"*20)

def count_hits(path):
    hits = 0
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        if float(tmp[-1]) < 0.05:
            hits+=1
    return hits

# for ANOVA output; assume last 3 columns are pval, qval, and status

def significant_pvals_qvals(fname, pcol=-2, qcol=-1):
    pvals, qvals = [], []
    #with open(fname) as f:
    #    lines = f.readlines()
    #for line in lines[2:]:
    for line in open(fname):
        if line[0]=='#': continue
        if "pval" in line and "padj" in line: continue
        cols = line.split("\t")
        # Read in position as int, and readcounts as float
        pvals.append(float(cols[pcol]))
        qvals.append(float(cols[qcol]))

    return (list(filter(lambda p: p < 0.05, pvals)), list(filter(lambda q: q < 0.05, qvals)))

