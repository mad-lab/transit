#__all__ = []

from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]

import pytransit.analysis.base as base

import pytransit.analysis.gumbel as gumbel
import pytransit.analysis.example as example
import pytransit.analysis.tn5gaps as tn5gaps
import pytransit.analysis.binomial as binomial
import pytransit.analysis.griffin as griffin
import pytransit.analysis.resampling as resampling
import pytransit.analysis.hmm as hmm
import pytransit.analysis.rankproduct as rankproduct
import pytransit.analysis.gi as gi
import pytransit.analysis.utest as utest
import pytransit.analysis.normalize as normalize
import pytransit.analysis.pathway_enrichment as pathway_enrichment
import pytransit.analysis.anova as anova
import pytransit.analysis.zinb as zinb
import pytransit.analysis.tnseq_stats as tnseq_stats

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
methods["pathway_enrichment"]=pathway_enrichment.GSEAAnalysis()
methods["tnseq_stats"]=tnseq_stats.TnseqStats()

# EXPORT METHODS
import pytransit.analysis.norm as norm

export_methods = {}
export_methods["norm"] = norm.NormAnalysis()

