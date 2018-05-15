#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]


import combined_wig
import igv

# EXPORT METHODS
methods = {}
methods["combined_wig"] = combined_wig.CombinedWigExport()
methods["igv"] = igv.IGVExport()



