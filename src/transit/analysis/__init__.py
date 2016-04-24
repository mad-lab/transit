#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]

import gumbel
import example
import example2
import globalruns
import binomial
import griffin
import resampling
import hmm
import rankproduct


methods = {}
methods["gumbel"] = {"module":gumbel, "gui":gumbel.GumbelGUI, "method":gumbel.Gumbel}
methods["example"] = {"module":example, "gui":example.ExampleGUI, "method":example.ExampleMethod}
methods["binomial"] = {"module":binomial, "gui":binomial.BinomialGUI, "method":binomial.Binomial}
methods["griffin"] = {"module":griffin, "gui":griffin.GriffinGUI, "method":griffin.Griffin}
methods["hmm"] = {"module":hmm, "gui":hmm.hmmGUI, "method":hmm.HMM}
methods["resampling"] = {"module":resampling, "gui":resampling.resamplingGUI, "method":resampling.Resampling}
methods["rankproduct"] = {"module":rankproduct, "gui":rankproduct.rankproductGUI, "method":rankproduct.Rankproduct}
methods["globalruns"] = {"module":globalruns, "gui":globalruns.globalrunsGUI, "method":globalruns.GlobalGumbel}

def defineGUI(wxobj):
    gui = {}
    for m in methods:
        gui[m] = methods[m]["gui"](wxobj)
    return gui



