import sys
import os

basedir = os.path.dirname(__file__)
sys.path.insert(0, basedir + '/../src/')
#sys.path.insert(0, '/home/travis/build/mad-lab/transit/src/')

import shutil
import unittest

from transit_test import *

import pytransit
import pytransit.norm_tools as norm_tools

# Single condition methods
from pytransit.analysis.gumbel import GumbelMethod
from pytransit.analysis.binomial import BinomialMethod
from pytransit.analysis.griffin import GriffinMethod
from pytransit.analysis.hmm import HMMMethod
from pytransit.analysis.anova import AnovaMethod
from pytransit.analysis.zinb import ZinbMethod

# Comparative methods
from pytransit.analysis.resampling import ResamplingMethod
from pytransit.analysis.rankproduct import RankProductMethod
from pytransit.analysis.utest import UTestMethod

# Genetic Interactions
from pytransit.analysis.gi import GIMethod

def significant_pvals_qvals(fname, pcol=-2, qcol=-1):
    print(fname)
    pvals, qvals = [], []
    with open(fname) as f:
        lines = f.readlines()
    for line in lines[1:]:
        if line[0]=='#': continue
        cols = line.split("\t")
        # Read in position as int, and readcounts as float
        pvals.append(float(cols[pcol]))
        qvals.append(float(cols[qcol]))

    return (filter(lambda p: p < 0.05, pvals), filter(lambda q: q < 0.05, qvals))


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

    def test_anova(self):
        args = [combined_wig, annotation, samples_metadata, output, "--ignore-conditions", "Unknown"]
        G = AnovaMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        sig_qvals.sort()
        self.assertEqual(
            len(sig_pvals),
            196,
            "sig_pvals expected: %d, actual: %d" % (196, len(sig_pvals)))
        self.assertEqual(
            len(sig_qvals),
            37,
            "sig_qvals expected: %d, actual: %d" % (37, len(sig_qvals)))

    def test_zinb(self):
        args = [combined_wig, annotation, samples_metadata, output, "--ignore-conditions", "Unknown"]
        G = ZinbMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
        sig_qvals.sort()
        self.assertEqual(
            len(sig_pvals),
            397,
            "sig_pvals expected: %d, actual: %d" % (397, len(sig_pvals)))
        self.assertEqual(
            len(sig_qvals),
            106,
            "sig_qvals expected: %d, actual: %d" % (106, len(sig_qvals)))

    #def test_resampling_histogram(self):
    #    args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-s", "1000", "-h"]
    #    G = ResamplingMethod.fromargs(args)
    #    G.Run()
    #    self.assertTrue(os.path.exists(output))
    #    hist_path = output.rsplit(".", 1)[0] + "_histograms"
    #    self.assertTrue(os.path.isdir(hist_path))


    def test_GI(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000"]
        G = GIMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_utest(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output]
        G = UTestMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))


    def test_GI(self):
        args = [ctrl_data_txt, exp_data_txt, ctrl_data_txt, exp_data_txt, annotation, output,
                    "-s", "1000"]
        G = GIMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))



if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestMethods)
    #unittest.TextTestRunner(verbosity=2).run(suite)
