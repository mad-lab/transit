import scipy
import pandas as pd

import time

import base
import pytransit
import pytransit.transit_tools as transit_tools

############# GUI ELEMENTS ##################

short_name = ""
long_name = ""
short_desc = ""
long_desc = """"""

transposons = ["", ""]
columns = []

class AnovaAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, AnovaMethod)
def main():
    print("ANOVA example")

class AnovaMethod(base.MultiConditionMethod):
    """
    anova
    """
    def __init__(self, combined_wig, metadata, annotation, output_file):
        print("callig anova method")
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file)

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        return self("combined_wig", "metadata", "annotation", "output_file")

    def Run(self):
        self.transit_message("Starting Anova analysis")
        start_time = time.time()

if __name__ == "__main__":
    main()

