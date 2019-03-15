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


NOFLAG_PRIMER = [
        "# TA_sites: 74605",
        "# TAs_hit: 914",
        "# mapped_reads (both R1 and R2 map into genome, and R2 has a proper barcode): 967",
        "# density: 0.012",
        "# NZ_mean (among templates): 1.0",
        "# FR_corr (Fwd templates vs. Rev templates): 0.019",
        "# transposon type: Himar1",
        "# protocol type: Sassetti",
        "# primer_matches: 8 reads (0.8%) contain CTAGAGGGCCCAATTCGCCCTATAGTGAGT (Himar1)"
        ]

MME1_PROTOCOL = [
        "# TA_sites: 74605",
        "# TAs_hit: 34",
        "# mapped_reads (both R1 and R2 map into genome, and R2 has a proper barcode): 967",
        "# density: 0.000",
        "# NZ_mean (among templates): 1.0",
        "# transposon type: Himar1",
        "# protocol type: Mme1",
        "# primer_matches: 8 reads (0.8%) contain CTAGAGGGCCCAATTCGCCCTATAGTGAGT (Himar1)"
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

MULTICONTIG = [
        "# TA_sites:",
        "#   a: 89994",
        "#   b: 646",
        "#   c: 664",
        "# TAs_hit:",
        "#   a: 63",
        "#   b: 0",
        "#   c: 0",
        "# density:",
        "#   a: 0.001",
        "#   b: 0.000",
        "#   c: 0.000",
        "# max_count (among templates):",
        "#   a: 1",
        "#   b: 0",
        "#   c: 0",
        "# max_site (coordinate):",
        "#   a: 4977050",
        "#   b: 57441",
        "#   c: 38111" ]

MULTICONTIG_AUTO_IDS = [
        "# TA_sites:",
        "#   1: 89994",
        "#   2: 646",
        "#   3: 664",
        "# TAs_hit:",
        "#   1: 63",
        "#   2: 0",
        "#   3: 0",
        "# density:",
        "#   1: 0.001",
        "#   2: 0.000",
        "#   3: 0.000",
        "# max_count (among templates):",
        "#   1: 1",
        "#   2: 0",
        "#   3: 0",
        "# max_site (coordinate):",
        "#   1: 4977050",
        "#   2: 57441",
        "#   3: 38111" ]

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
        diff = set(expected) - lines
        if (len(diff) == 0):
            return True
        print("Diff: ", diff)
        return False

class TestTPP(TransitTestCase):

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_noflag_primer(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-protocol", "sassetti"])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), NOFLAG_PRIMER))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_flag_primer(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-himar1", "-flags", "-k 1"])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), FLAG_PRIMER))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_protocol_mme1(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", h37fna, "-reads1", reads1, "-output", tpp_output_base, "-protocol", "Mme1"])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), MME1_PROTOCOL))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_multicontig_empty_prefix(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", test_multicontig, "-reads1", test_multicontig_reads1, "reads2", test_multicontig_reads2, "-output", tpp_output_base, "-replicon-ids", "a,b,c", "-maxreads", "10000", "-primer", ""])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), MULTICONTIG))

    @unittest.skipUnless(len(bwa_path) > 0, "requires BWA")
    def test_tpp_multicontig_auto_replicon_ids(self):
        (args, kwargs) = cleanargs(["-bwa", bwa_path, "-ref", test_multicontig, "-reads1", test_multicontig_reads1, "reads2", test_multicontig_reads2, "-output", tpp_output_base, "-replicon-ids", "auto", "-maxreads", "10000", "-primer", ""])
        tppMain(*args, **kwargs)
        self.assertTrue(verify_stats("{0}.tn_stats".format(tpp_output_base), MULTICONTIG_AUTO_IDS))

if __name__ == '__main__':
    unittest.main()


