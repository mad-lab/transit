import sys
import os
import inspect

basedir = os.path.dirname(__file__)
sys.path.insert(0, basedir + '/../src/')
#sys.path.insert(0, '/home/travis/build/mad-lab/transit/src/')

import shutil
import unittest

from transit_test import *
from pytpp.tpp_tools import cleanargs

import pytpp.__main__

tppMain = pytpp.__main__.main

def get_bwa():
    if (os.path.exists("/usr/bin/bwa")):
        return "/usr/bin/bwa"
    elif (os.path.exists("/usr/local/bin/bwa")):
        return "/usr/local/bin/bwa"
    return ""

bwa_path = get_bwa()


NOFLAG_NOPRIMER = [
        "# density: 0.000",
        "# NZ_mean (among templates): 1.0",
        "# FR_corr (Fwd templates vs. Rev templates): -0.000"
        ]

FLAG_NOPRIMER = [
        "# density: 0.000",
        "# NZ_mean (among templates): 1.0",
        "# FR_corr (Fwd templates vs. Rev templates): -0.000",
        "# bwa flags: -k 1"
        ]

NOFLAG_PRIMER = [
        "# TA_sites: 74605",
        "# TAs_hit: 914",
        "# mapped_reads (both R1 and R2 map into genome, and R2 has a proper barcode): 967",
        "# density: 0.012",
        "# NZ_mean (among templates): 1.0",
        "# FR_corr (Fwd templates vs. Rev templates): 0.019"
        ]

FLAG_PRIMER = [
        "# TA_sites: 74605",
        "# TAs_hit: 914",
        "# bwa flags: -k 1",
        "# mapped_reads (both R1 and R2 map into genome, and R2 has a proper barcode): 967",
        "# density: 0.012",
        "# NZ_mean (among templates): 1.0",
        "# FR_corr (Fwd templates vs. Rev templates): 0.019"
        ]

def get_stats(path):
    for line in open(path):
        if line.startswith("#"):
            print(line[1:].split(":"))
            continue
        break
    return [float(x) if (type(x) != list) else x for x in tmp[2:]]

def verify_stats(stats_file, expected):
    with open(stats_file) as f:
        lines = set([line.strip() for line in f])
        print(lines)
        print(set(expected) - lines)
        return len(set(expected) - lines) == 0
    return False

class TestTPP(TransitTestCase):

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_noflag_primer(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-himar1"])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), NOFLAG_PRIMER))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_flag_primer(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-himar1", "-flags", "-k 1"])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), FLAG_PRIMER))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_noflag_noprimer(self):
        # with self.assertRaises(SystemExit):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-primer", " "])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), NOFLAG_NOPRIMER))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_flag_noprimer(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-flags", "-k 1", "-primer", " "])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), FLAG_NOPRIMER))

if __name__ == '__main__':
    unittest.main()


