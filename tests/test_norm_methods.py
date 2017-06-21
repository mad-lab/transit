import sys

sys.path.insert(0, '/home/travis/build/mad-lab/transit/src/')

import shutil
import unittest
import os
import numpy

import pytransit.norm_tools as norm_tools
import pytransit.tnseq_tools as tnseq_tools

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

all_data_list = [ctrl_rep1, ctrl_rep2, exp_rep1, exp_rep2, exp_rep3]

annotation = "../src/pytransit/genomes/H37Rv.prot_table"

output = "testoutput.txt"



# RAW STATISTICS:

raw_ctrl_rep1 = (0.41855103545338784, 53.913866362844317, 128.81073464420675, 70.0, 3855.0, 4022244.0, 4.022213863266365, 33.02399374787363)

raw_ctrl_rep2 = (0.51571610481871188, 86.141009315729505, 167.03183885640027, 89.0, 5944.0, 6426550.0, 3.975522680364774, 33.45593924344892)


raw_exp_rep1 = (0.43946116212050129, 52.944722203605657, 120.47645336424084, 56.0, 47546.0, 3949941.0, 54.83729352080048, 4237.703504066973)

raw_exp_rep2 = (0.43891160109912203, 53.013082233094295, 120.78305084745763, 46.0, 217960.0, 3955041.0, 105.78284470780153, 14216.199240166581)

raw_exp_rep3 = (0.35898398230681589, 60.667260907445879, 168.99712493465759, 60.0, 102013.0, 4526081.0, 42.17620960361282, 2327.971405759911)


raw_means = [53.913866362844317, 86.141009315729505, 52.944722203605657, 53.013082233094295, 60.667260907445879]


def count_hits(path):
    hits = 0
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        if float(tmp[-1]) < 0.05:
            hits+=1
    return hits

class TestMethods(unittest.TestCase):
 
    def setUp(self):

        print ""
        print "#"*20
        if os.path.exists(output):
            print "Removing output file..."
            os.remove(output) 

        #hist_path = output.rsplit(".", 1)[0] + "_histograms"
        #if os.path.exists(hist_path):
        #    print "Removing histogram dir..."
        #    shutil.rmtree(hist_path)
    
        genes_path = output.rsplit(".", 1)[0] + "_genes" + output.rsplit(".", 1)[1]
        if os.path.exists(genes_path):
            print "Removing genes file..."
            os.remove(genes_path)


    def test_nonorm(self):
        data,position = tnseq_tools.get_data(all_data_list)
        norm_data,factors = norm_tools.normalize_data(data, "nonorm")
        self.assertTrue((factors == numpy.array([ 1.])).all())
        N = len(all_data_list)
        for k in range(N):
           self.assertEqual(numpy.mean(norm_data[k]), raw_means[k])


    def test_TTR(self):
        N = len(all_data_list)
        data,position = tnseq_tools.get_data(all_data_list)
        norm_data,factors = norm_tools.normalize_data(data, "TTR")
        self.assertFalse((factors == numpy.ones(N)).all())
        for k in range(N):
           self.assertNotEqual(numpy.mean(norm_data[k]), raw_means[k])
#    """

    def test_resampling_nonorm(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "nonorm"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 200)
#

    def test_resampling_TTR(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "TTR"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 10)

#

    def test_resampling_NZMean(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "nzmean"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 20)

#

    def test_resampling_TotReads(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "totreads"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 15)

#

    def test_resampling_Quantile(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "quantile"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 50)

#
    
    def test_resampling_ZINFNB(self):
        args = [ctrl_rep1, ctrl_rep2, annotation, output, "-s", "1000", "-n", "zinfnb"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))

#
    """
    def test_resampling_BGC(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000", "-n", "betageom"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 20)
    """


#    """
#
    """
    def test_resampling_aBGC(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "-s", "1000", "-n", "aBGC"]
        G = ResamplingMethod.fromargs(args)
        G.Run()
        self.assertTrue(os.path.exists(output))
    """

    



    


 
if __name__ == '__main__':
    unittest.main()
