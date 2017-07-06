import sys

sys.path.insert(0, '../src/')
#sys.path.insert(0, '/home/travis/build/mad-lab/transit/src/')

import shutil
import unittest
import os

from transit_test import *

import pytpp
import pytpp.__main__

tppMain = pytpp.__main__.main


NOFLAG_PRIMER = [float(x) for x in "1000    983 937 0   937 932 932 1.0 74605   907 2   0.0121573621071 4348593 1.02756339581   0.0190915393988 1.0 8   9   0   9".split()]

FLAG_PRIMER = [float(x) for x in "1000  983 935 0   935 931 931 1.0 74605   906 2   0.0121439581797 4348593 1.02759381898   0.0191247691862 1.0 8   9   0   9".split()]

NOFLAG_NOPRIMER = [float(x) for x in "".split()]

FLAG_NOPRIMER = [float(x) for x in "".split()]


def get_stats(path):
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        break
    return [float(x) for x in tmp[2:]]


class TestTPP(TransitTestCase):
 
    def test_tpp_noflag_primer(self):

        arguments = ["-bwa", "/usr/bin/bwa", "-ref", "H37Rv.fna", "-reads1", "test.fastq", "-output",
                     "tpp_output_noflag_primer", "-himar1"]
        tppMain(arguments)
        stats = get_stats("tpp_output_noflag_primer.tn_stats")
        self.assertTrue(NOFLAG_PRIMER == stats) 


    def test_tpp_flag_primer(self):

        arguments = ["-bwa", "/usr/bin/bwa", "-ref", "H37Rv.fna", "-reads1", "test.fastq", "-output",
                     "tpp_output_flag_primer", "-himar1", "-flags", "-k 1"]
        tppMain(arguments)
        stats = get_stats("tpp_output_flag_primer.tn_stats")
        self.assertTrue(FLAG_PRIMER == stats)


    """
    def test_tpp_noflag_noprimer(self):

        arguments = ["-bwa", "/usr/bin/bwa", "-ref", "H37Rv.fna", "-reads1", "test.fastq", "-output",
                     "tpp_output_noflag_noprimer", "-primer", " "]
        tppMain(arguments)
        stats = get_stats("tpp_output_noflag_noprimer.tn_stats")
        self.assertTrue(NOFLAG_NOPRIMER == stats)


    
    def test_tpp_flag_noprimer(self):
        
        arguments = ["-bwa", "/usr/bin/bwa", "-ref", "H37Rv.fna", "-reads1", "test.fastq", "-output",
                     "tpp_output_flag_noprimer", "-flags", "-k 1", "-primer", " "]
        tppMain(arguments)
        stats = get_stats("tpp_output_flag_noprimer.tn_stats")
        self.assertTrue(FLAG_PRIMER == stats)
    """


 
if __name__ == '__main__':
    unittest.main()


