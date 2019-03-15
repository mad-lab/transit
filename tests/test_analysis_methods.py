import matplotlib
matplotlib.use('Agg')

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

# Comparative methods
from pytransit.analysis.resampling import ResamplingMethod
from pytransit.analysis.rankproduct import RankProductMethod
from pytransit.analysis.utest import UTestMethod

# Genetic Interactions
from pytransit.analysis.gi import GIMethod

class TestMethods(TransitTestCase):
    def test_Gumbel(self):
        args = [ctrl_data_txt, small_annotation, output, "-s", "1000", "-b", "100"]
        G = GumbelMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_Binomial(self):
        args = [ctrl_data_txt, small_annotation, output, "-s", "1000", "-b", "100"]
        G = BinomialMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_Griffin(self):
        args = [ctrl_data_txt, small_annotation, output, "-s", "1000", "-b", "100"]
        G = GriffinMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))


    def test_HMM(self):
        args = [mini_wig, small_annotation, output]
        G = HMMMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        genes_path = output.rsplit(".", 1)[0] + "_genes." + output.rsplit(".", 1)[1]
        self.assertTrue(os.path.exists(genes_path))


    def test_resampling(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        self.assertLessEqual(
                abs(len(sig_pvals) - 37),
                2,
                "sig_pvals expected in range: %s, actual: %d" % ("[35, 39]", len(sig_qvals)))
        self.assertLessEqual(
                abs(len(sig_qvals) - 35),
                1,
                "sig_qvals expected in range: %s, actual: %d" % ("[34, 36]", len(sig_qvals)))

    def test_resampling_adaptive(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-a"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        print(len(sig_pvals))
        print(len(sig_qvals))
        self.assertLessEqual(
                abs(len(sig_pvals) - 37),
                2,
                "sig_pvals expected in range: %s, actual: %d" % ("[35, 39]", len(sig_qvals)))
        self.assertLessEqual(
                abs(len(sig_qvals) - 35),
                1,
                "sig_qvals expected in range: %s, actual: %d" % ("[34, 36]", len(sig_qvals)))

    def test_resampling_histogram(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-s", "1000", "-h"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        self.assertTrue(
                os.path.isdir(hist_path),
                "histpath expected: %s" % (hist_path))

    def test_resampling_multistrain(self):
        args = [ctrl_data_txt, exp_data_txt, ','.join([small_annotation, small_annotation]), output, "-h"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        self.assertTrue(
                os.path.isdir(hist_path),
                "histpath expected: %s" % (hist_path))

    def test_anova(self):
        args = [combined_wig, annotation, samples_metadata, output, "--ignore-conditions", "Unknown"]
        G = AnovaMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
        sig_qvals.sort()
        self.assertEqual(
            len(sig_pvals),
            196,
            "sig_pvals expected: %d, actual: %d" % (196, len(sig_pvals)))
        self.assertEqual(
            len(sig_qvals),
            36,
            "sig_qvals expected: %d, actual: %d" % (36, len(sig_qvals)))

    def test_GI(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-s", "1000"]
        G = GIMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

    def test_utest(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output]
        G = UTestMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))


    def test_GI(self):
        args = [ctrl_data_txt, exp_data_txt, ctrl_data_txt, exp_data_txt, small_annotation, output,
                    "-s", "1000"]
        G = GIMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestMethods)
    #unittest.TextTestRunner(verbosity=2).run(suite)
