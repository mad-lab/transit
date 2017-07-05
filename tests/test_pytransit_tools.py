import sys

sys.path.insert(0, '../src/')

import shutil
import unittest
import os
import numpy

import pytransit.norm_tools as norm_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.stat_tools as stat_tools


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


class TestTnSeqTools(unittest.TestCase):
 
    def setUp(self):

        print "#"*20
        print self.id()
        print "#"*20
        print "\n"
        if os.path.exists(output):
            print "Removing output file..."
            os.remove(output) 
    
        genes_path = output.rsplit(".", 1)[0] + "_genes" + output.rsplit(".", 1)[1]
        if os.path.exists(genes_path):
            print "Removing genes file..."
            os.remove(genes_path)


    def test_read_data(self):
        data,position = tnseq_tools.get_data(all_data_list)
        K,N = data.shape
        
        self.assertEqual(K, 5)
        self.assertGreater(N, 70000)
        
    def test_genes_creation_fromwig(self):
        G = tnseq_tools.Genes(all_data_list, annotation)
        N = len(G)
        test_orf = "Rv0001"
        test_name = "dnaA"
        #Test list creation
        self.assertGreater(N, 3000)

        #Test dictionary lookup + data
        self.assertEqual(G[test_orf].orf, test_orf)
        self.assertEqual(G[test_orf].name, test_name)

        #Test list lookup + data
        self.assertEqual(G[0].orf, test_orf)
        self.assertEqual(G[0].name, test_name)
    

    def test_genes_creation_fromdata(self):
        data,position = tnseq_tools.get_data(all_data_list)
        Kreps,Nsites = data.shape
        G = tnseq_tools.Genes([], annotation, data=data, position=position)
        N = len(G)
        test_orf = "Rv0001"
        test_name = "dnaA"
        #Test list creation
        self.assertGreater(N, 3000)

        #Test dictionary lookup + data
        self.assertEqual(G[test_orf].orf, test_orf)
        self.assertEqual(G[test_orf].name, test_name)

        #Test list lookup + data
        self.assertEqual(G[0].orf, test_orf)
        self.assertEqual(G[0].name, test_name)


    def test_file_types(self):
        types = tnseq_tools.get_file_types(all_data_list)
        types = set(types)
        self.assertEqual(len(types), 1)
        self.assertTrue("himar1" in types)



    def test_normalization(self):
        N = len(all_data_list)
        data,position = tnseq_tools.get_data(all_data_list)
        norm_data,factors = norm_tools.normalize_data(data, "TTR")
        self.assertFalse((factors == numpy.ones(N)).all())
   
    
    



 
if __name__ == '__main__':
    unittest.main()
