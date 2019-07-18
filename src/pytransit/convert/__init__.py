#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]


from pytransit.convert import gff_to_prot_table

# EXPORT METHODS
methods = {}
methods["gff_to_prot_table"] = gff_to_prot_table.GffProtConverter()

