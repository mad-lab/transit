import sys

sys.path.insert(0, '../src/')
#sys.path.insert(0, '/home/travis/build/mad-lab/transit/src/')

import shutil
import unittest
import os

from transit_test import *

import pytransit
import pytransit.norm_tools as norm_tools

from pytransit.analysis.gumbel import GumbelMethod
from pytransit.analysis.binomial import BinomialMethod
from pytransit.analysis.griffin import GriffinMethod
from pytransit.analysis.hmm import HMMMethod

from pytransit.analysis.resampling import ResamplingMethod
from pytransit.analysis.rankproduct import RankProductMethod


class TestMethods(TransitTestCase):
 
 
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
    

    #def test_resampling_histogram(self):
    #    args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-s", "1000", "-h"]
    #    G = ResamplingMethod.fromargs(args)
    #    G.Run()
    #    self.assertTrue(os.path.exists(output))
    #    hist_path = output.rsplit(".", 1)[0] + "_histograms"
    #    self.assertTrue(os.path.isdir(hist_path))


 
if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestMethods)
    #unittest.TextTestRunner(verbosity=2).run(suite)
