import sys

sys.path.insert(0, '../src/')

import shutil
import unittest
import os

import pytransit.norm_tools as norm_tools

from pytransit.analysis.gumbel import GumbelMethod
from pytransit.analysis.binomial import BinomialMethod
from pytransit.analysis.griffin import GriffinMethod
from pytransit.analysis.hmm import HMMMethod

from pytransit.analysis.resampling import ResamplingMethod
from pytransit.analysis.rankproduct import RankProductMethod


ctrl_rep1 = "../src/pytransit/data/glycerol_H37Rv_rep1.wig"
ctrl_rep2 = "../src/pytransit/data/glycerol_H37Rv_rep2.wig"
ctrl_data_txt = ",".join([ctrl_rep1, ctrl_rep2])

exp_rep1 = "../src/pytransit/data/cholesterol_H37Rv_rep1.wig"
exp_rep2 = "../src/pytransit/data/cholesterol_H37Rv_rep2.wig"
exp_rep3 = "../src/pytransit/data/cholesterol_H37Rv_rep3.wig"
exp_data_txt = ",".join([exp_rep1, exp_rep2, exp_rep3])

annotation = "../src/pytransit/genomes/H37Rv.prot_table"

output = "testoutput.txt"

class TestMethods(unittest.TestCase):
 
    def setUp(self):

        print ""
        print "#"*20
        if os.path.exists(output):
            print "Removing output file..."
            os.remove(output) 

        hist_path = output.rsplit(".", 1)[0] + "_histograms"
        if os.path.exists(hist_path):
            print "Removing histogram dir..."
            shutil.rmtree(hist_path)
    
        genes_path = output.rsplit(".", 1)[0] + "_genes" + output.rsplit(".", 1)[1]
        if os.path.exists(genes_path):
            print "Removing genes file..."
            os.remove(genes_path)

 
    def test_Gumbel(self):
        args = [ctrl_data_txt, annotation, output, "-s", "1000", "-b", "100"]
        G = GumbelMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_Binomial(self):
        args = [ctrl_data_txt, annotation, output, "-s", "1000", "-b", "100"]
        G = BinomialMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_Griffin(self):
        args = [ctrl_data_txt, annotation, output, "-s", "1000", "-b", "100"]
        G = GriffinMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))


    def test_HMM(self):
        args = [ctrl_data_txt, annotation, output, "-s", "1000", "-b", "100"]
        G = HMMMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        genes_path = output.rsplit(".", 1)[0] + "_genes." + output.rsplit(".", 1)[1]
        self.assertTrue(os.path.exists(genes_path))        


    def test_resampling(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))


    def test_resampling_adaptive(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000", "-a"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
    

    def test_resampling_histogram(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000", "-h"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hist_path = output.rsplit(".", 1)[0] + "_histograms"
        self.assertTrue(os.path.isdir(hist_path))


 
if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestMethods)
    #unittest.TextTestRunner(verbosity=2).run(suite)
