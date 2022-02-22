#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]

from pytransit.analysis import base

from pytransit.analysis import gumbel
from pytransit.analysis import example
from pytransit.analysis import tn5gaps
from pytransit.analysis import binomial
from pytransit.analysis import griffin
from pytransit.analysis import resampling
from pytransit.analysis import hmm
from pytransit.analysis import rankproduct
from pytransit.analysis import gi
from pytransit.analysis import utest
from pytransit.analysis import normalize
from pytransit.analysis import pathway_enrichment
from pytransit.analysis import anova
from pytransit.analysis import zinb
from pytransit.analysis import tnseq_stats
from pytransit.analysis import corrplot
from pytransit.analysis import heatmap
from pytransit.analysis import ttnfitness

methods = {}
methods["example"] = example.ExampleAnalysis()
methods["gumbel"] = gumbel.GumbelAnalysis()
methods["binomial"] = binomial.BinomialAnalysis()
methods["griffin"] = griffin.GriffinAnalysis()
methods["hmm"] = hmm.HMMAnalysis()
methods["resampling"] = resampling.ResamplingAnalysis()
methods["tn5gaps"] = tn5gaps.Tn5GapsAnalysis()
methods["rankproduct"] = rankproduct.RankProductAnalysis()
methods["utest"] = utest.UTestAnalysis()
methods["GI"] = gi.GIAnalysis()
methods["anova"] = anova.AnovaAnalysis()
methods["zinb"] = zinb.ZinbAnalysis()
#methods["mcce"] = mcce.MCCEAnalysis()
#methods["mcce2"] = mcce2.MCCE2Analysis()
#methods["motifhmm"] = motifhmm.MotifHMMAnalysis()
methods["normalize"] = normalize.Normalize()
methods["pathway_enrichment"]=pathway_enrichment.PathwayAnalysis()
methods["tnseq_stats"]=tnseq_stats.TnseqStats()
methods["corrplot"]=corrplot.Corrplot()
methods["heatmap"]=heatmap.Heatmap()
methods["ttnfitness"]=ttnfitness.TTNFitnessAnalysis()

# EXPORT METHODS
from pytransit.analysis import norm

export_methods = {}
export_methods["norm"] = norm.NormAnalysis()

