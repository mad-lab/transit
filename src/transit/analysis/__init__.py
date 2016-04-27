#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]

import gumbel
import example
import globalruns
import binomial
import griffin
import resampling
import hmm
import rankproduct



methods = {}
methods["gumbel"] = {"module":gumbel, "method":gumbel.Gumbel}
methods["example"] = {"module":example, "method":example.Example}
methods["binomial"] = {"module":binomial, "method":binomial.Binomial}
methods["griffin"] = {"module":griffin, "method":griffin.Griffin}
methods["hmm"] = {"module":hmm, "method":hmm.HMM}
methods["resampling"] = {"module":resampling, "method":resampling.Resampling}
methods["rankproduct"] = {"module":rankproduct, "method":rankproduct.Rankproduct}
methods["tn5gaps"] = {"module":tn5gaps, "method":tn5gaps.Tn5Gaps}



