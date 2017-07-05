import sys
sys.path.insert(0, '../src/')

import os
import shutil
import unittest
import os
import numpy

from transit_test import *

import pytransit.norm_tools as norm_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.stat_tools as stat_tools



class TestTnSeqTools(TransitTestCase):
 

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
