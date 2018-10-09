#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]

import base

import gumbel
import example
import tn5gaps
import binomial
import griffin
import resampling
import hmm
import rankproduct
import gi
import utest
import normalize
import pathway_enrichment #08/22/2018 by Ivan
import anova
import tnseq_stats

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
#methods["mcce"] = mcce.MCCEAnalysis()
#methods["mcce2"] = mcce2.MCCE2Analysis()
#methods["motifhmm"] = motifhmm.MotifHMMAnalysis()
methods["normalize"] = normalize.Normalize()
methods["pathway_enrichment"]=pathway_enrichment.GSEAAnalysis()
methods["tnseq_stats"]=tnseq_stats.TnseqStats()

# EXPORT METHODS
import norm
export_methods = {}
export_methods["norm"] = norm.NormAnalysis()

