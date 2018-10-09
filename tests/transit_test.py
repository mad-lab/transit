import unittest
import os

basedir = os.path.dirname(__file__)
ctrl_rep1 = basedir + "/../src/pytransit/data/glycerol_H37Rv_rep1.wig"
ctrl_rep2 = basedir + "/../src/pytransit/data/glycerol_H37Rv_rep2.wig"
ctrl_data_txt = ",".join([ctrl_rep1, ctrl_rep2])
print(ctrl_data_txt)

combined_wig = basedir + "/../src/pytransit/data/cholesterol_glycerol_combined.dat"
samples_metadata = basedir + "/../src/pytransit/data/samples_metadata_cg.txt"

exp_rep1 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep1.wig"
exp_rep2 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep2.wig"
exp_rep3 = basedir + "/../src/pytransit/data/cholesterol_H37Rv_rep3.wig"
exp_data_txt = ",".join([exp_rep1, exp_rep2, exp_rep3])

all_data_list = [ctrl_rep1, ctrl_rep2, exp_rep1, exp_rep2, exp_rep3]

annotation = basedir + "/../src/pytransit/genomes/H37Rv.prot_table"
small_annotation = basedir + "/test.prot_table"
output = basedir + "/testoutput.txt"


class TransitTestCase(unittest.TestCase):

    def setUp(self):

        # Print header
        self.header()

        # Check if there were output files and remove them
        if os.path.exists(output):
            print "Removing output file..."
            os.remove(output)

        genes_path = output.rsplit(".", 1)[0] + "_genes" + output.rsplit(".", 1)[1]
        if os.path.exists(genes_path):
            print "Removing genes file..."
            os.remove(genes_path)




    def header(self):
        print "\n"
        print "#"*20
        print self.id()
        print "#"*20



def count_hits(path):
    hits = 0
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        if float(tmp[-1]) < 0.05:
            hits+=1
    return hits

