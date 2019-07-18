import sys
import os
import io

hasR = False
try:
    import rpy2.robjects
    hasR = True
except Exception as e:
    hasR = False

import base
import numpy
import scipy
import pdfkit
import statsmodels.stats.multitest
import statsmodels.stats.multicomp
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
from multitnseq_helpers import def_r_samples_corrplot, def_r_clustering, def_r_make_heatmap, def_r_conditions_corrplot

###### GUI ELEMENTS ######
short_name = "MultiTnSeq"
long_name = "MultiTnSeq"
short_desc = "Run the MultiTnSeq pipeline"
long_desc = """Run the MultiTnSeq pipeline"""
transposons = ["", ""]
EOL = "\n"

class MultiTnSeqAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, MultiTnSeqMethod)

class MultiTnSeqMethod(base.AnalysisMethod):
    """
    MultiTnSeq pipeline, uses R scripts
    """
    def __init__(self, fna, annotation, combined_wig, metadata, cond_metadata, output_file, K, D, debatch, stdnorm):
        base.AnalysisMethod.__init__(self, short_name, long_name, short_desc, long_desc, output_file, annotation)
        self.ref = fna
        self.annotation = annotation
        self.combined_wig = combined_wig
        self.samples_metadata_fname = metadata
        self.conditions_metadata_fname = cond_metadata
        self.output_file = output_file
        self.K = K
        self.D = D
        self.debatch = False
        self.stdnorm = False
        self.rundir = ""

    @classmethod
    def fromargs(self, rawargs):
        if not hasR:
            print("Error: R and rpy2 (< 2.9.0) required to run ZINB analysis.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2<2.9.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(MultiTnSeqMethod.usage_string())
            sys.exit(0)

        fna = args[0]
        annotation = args[1]
        combined_wig = args[2]
        metadata = args[3]
        cond_metadata = args[4]
        output_file = args[5]
        K = kwargs.get("K", 10)
        D = kwargs.get("D", 5)
        debatch = kwargs.get("-debatch", False)
        stdnorm = kwargs.get("-stdnorm", False)

        return self(fna, annotation, combined_wig, metadata, cond_metadata, output_file, K, D, debatch, stdnorm)

    def read_genes(self, fname, descriptions=False):
      genes = []
      for line in open(fname):
        w = line.split('\t')
        data = [int(w[1])-1,int(w[2])-1,w[8],w[7],w[3]]
        if descriptions==True: data.append(w[0])
        genes.append(data)
      return genes

    def read_genome(self, filename):
      s = ""
      n = 0
      for line in open(filename):
        if n==0: n = 1 # skip first
        else: s += line[:-1]
      return s

    def hash_genes(self, genes,genome):
      hash = {}
      for gene in genes:
        a,b = gene[0],gene[1]
        for i in range(a,b+1):
          hash[i] = gene
      return hash

    def makefname(self, name,dir):
      if dir=="": return name
      else: return "%s/%s" % (dir,name)

    def Run(self):
      print("Running")
      # self.generate_pdf()
      # self.transit_message("Done.")

      Samples = tnseq_tools.Spreadsheet(self.makefname(self.samples_metadata_fname, self.rundir))
      Conditions = tnseq_tools.Spreadsheet(self.makefname(self.conditions_metadata_fname, self.rundir))

      genome = self.read_genome(self.makefname(self.ref,self.rundir))
      genes = self.read_genes(self.makefname(self.annotation,self.rundir))
      ghash = self.hash_genes(genes,genome)

      datasets = Samples.getcol("Filename") # filenames

      (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.makefname(self.combined_wig, self.rundir)) # data: rows are TA sites, cols are datasets
      # "datasets" are columns (indexes) in combined wig (data)
      # SampleIndexes are list of rows in Samples to be analyzed (contains Filenames and Conditions), i.e. which datasets (row indexes in data) represent each condition?
      # wigindexes are corresponding columns in data (from original order, matched by Filename)

      SampleIndexes,wigindexes = [],[]
      for cond in Conditions.keys:
        for row in range(Samples.nrows):
          if Samples.get(row,"Condition") == cond:
            fname = Samples.get(row,"Filename")
            if fname not in filenamesInCombWig:
              print "error: filename '%s' listed in samples metadata not found in combined wig file" % fname; sys.exit(0)
            SampleIndexes.append(row)
            wigindexes.append(filenamesInCombWig.index(fname))

      # reordering data matrix (select columns, order by cond) for normalization; now should be parallel to SampleIndexes
      data = data[wigindexes, :]
      data = data[wigindexes, :]
      neworder = [filenamesInCombWig[i] for i in wigindexes]
      filenamesInData = neworder

      # foreach cond, list of indexes into data
      cond2datasets = {}
      for cond in Conditions.keys:
        cond2datasets[cond] = []
        for row in range(Samples.nrows):
          if Samples.get(row,"Condition")==cond:
            fname = Samples.get(row,"Filename")
            cond2datasets[cond].append(filenamesInData.index(fname))
        if len(cond2datasets[cond]) == 0:
          print "error: no samples found in metadata for condition %s" % cond; sys.exit(0)

      (normed, factors) = norm_tools.normalize_data(data, method='TTR')
      Nds,Nsites = normed.shape # transposed: rows are datasets, cols are TA sites

      file = open(self.makefname("temp_counts_TTR.txt", self.rundir),"w")
      vals = "coord ORF gene".split()+[Samples.get(r,"Id") for r in SampleIndexes] # in order of Conditions
      file.write('\t'.join(vals)+EOL)
      for i,co in enumerate(sites):
        rv,gene = "igr","igr"
        if co in ghash:
          annot = ghash[co]
          rv,gene = annot[2],annot[3]
        vals = [str(co),rv,gene]+["%0.1f" % (int(x)) for x in list(normed[:,i])]
        file.write('\t'.join([str(x) for x in vals])+EOL)
      file.close()

      r_correlation = def_r_samples_corrplot()
      r_correlation("temp_counts_TTR.txt")

      # write table of stats (saturation,NZmean)
      file = open(self.makefname("temp_stats.txt", self.rundir),"w")
      file.write("dataset\tdensity\tmean_ct\tNZmean\tNZmedian\tmax_ct\ttotal_cts\tskewness\tkurtosis\n")
      for i in range(data.shape[0]):
        density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis = tnseq_tools.get_data_stats(data[i,:])
        vals = [filenamesInData[i], "%0.2f" % density, "%0.1f" % meanrd, "%0.1f" % nzmeanrd, "%d" % nzmedianrd, maxrd, totalrd, "%0.1f" % skew, "%0.1f" % kurtosis]
        file.write('\t'.join([str(x) for x in vals])+EOL)
      file.close()

      ################################

      # new version, hashes on Rv

      siteshash = {}
      for i,TA in enumerate(sites): siteshash[TA] = i

      TAsites = {}
      for g,(start,end,Rv,gene,strand) in enumerate(genes):
        siteindexes = []
        for i in range(start,end): # end+1?
          co = i+1
          if co in siteshash: siteindexes.append(siteshash[co])
        TAsites[Rv] = siteindexes



      if False: # print means for each gene for each dataset
        print '\t'.join(filenamesInData)
        print '\t'.join([Samples.get(x,"Id") for x in SampleIndexes])
        for (start,end,Rv,gene,strand) in genes:
          if len(TAsites[Rv])>0: # skip genes with no TA sites
            means = []
            for j in range(Nds):
              obs = normed[j,TAsites[Rv]]
              means.append(numpy.mean(obs))
          print '\t'.join([Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in means])
        #sys.exit(0)


      Counts = {} # sub-arrays of (datasets X sites) for genes X conds, where #TAs>0
      for (start,end,Rv,gene,strand) in genes:
        siteindexes = TAsites[Rv]
        local = normed[:,siteindexes]
        counts = []
        for j in range(Conditions.nrows):
          cond = Conditions.get(j,"Condition")
          obs = local[cond2datasets[cond],:]
          counts.append(obs) # sub-matrices
        Counts[Rv] = counts

      Means = {}
      for (start,end,Rv,gene,strand) in genes:
       if len(TAsites[Rv])>0: # skip genes with no TA sites
        means = []
        for j in range(Conditions.nrows):
          cond = Conditions.get(j,"Condition")
          obs = Counts[Rv][j]
          means.append(numpy.mean(obs))
        Means[Rv] = means

      file = open(self.makefname("temp_gene_means.txt", self.rundir),"w")
      file.write('\t'.join("ORF Gene TAs".split()+Conditions.getcol('Condition'))+EOL)
      for (start,end,Rv,gene,strand) in genes:
        if Rv not in Means: continue # skip genes with no TA sites
        means,sites = Means[Rv],TAsites[Rv]
        file.write('\t'.join([Rv,gene,str(len(sites))]+["%0.1f" % x for x in means])+EOL)
      file.close()

      ################################
      # do ANOVA to identify genes with significant variability
      # it is much faster to do this in python than R (and don't have to create the melted file)

      pvals,Rvs, = [],[] # lists
      for (start,end,Rv,gene,strand) in genes:
       if Rv in Means:
        countsvec = Counts[Rv]
        countsvec = [x.flatten() for x in countsvec]
        stat,pval = scipy.stats.f_oneway(*countsvec)
        pvals.append(pval)
        Rvs.append(Rv)

      pvals = numpy.array(pvals)
      mask = numpy.isfinite(pvals)
      qvals = numpy.full(pvals.shape,numpy.nan)
      qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[1] # BH, alpha=0.05

      p,q = {},{}
      for i,rv in enumerate(Rvs):
         p[rv],q[rv] = pvals[i],qvals[i]
      pvals,qvals = p,q, # hashes on Rv

      # post-hoc analysis to identify "high" and "low" subgroups
      Sortedmeans,Divisions = {},{} # sortedmeans: list of (mean,cond); divisions: list of (lowsubset,highsubset) or None
      for (start,end,Rv,gene,strand) in genes:
        if Rv in qvals and qvals[Rv]<0.05:
          #print Rv,gene
          countsvec = Counts[Rv]
          countsvec = [x.flatten() for x in countsvec]
          allcounts,allconds = numpy.array([]),numpy.array([])
          for j in range(len(Conditions.keys)):
            allcounts = numpy.append(allcounts,countsvec[j])
            for k in range(len(countsvec[j])):
              allconds = numpy.append(allconds,Conditions.keys[j])
          mc = statsmodels.stats.multicomp.MultiComparison(allcounts,allconds)
          tuk = mc.tukeyhsd() # alpha=0.1
          #print tuk
          reject,n = {},0
          groups = tuk.groupsunique.tolist()
          for j,group1 in enumerate(groups):
            for k,group2 in enumerate(groups):
              if j<k: reject[(group1,group2)] = reject[(group2,group1)] = tuk.reject[n]; n += 1
          sortedmeans = [(numpy.mean(countsvec[k]),Conditions.keys[k]) for k in range(len(countsvec))] # list of (mean,condname)
          sortedmeans = sorted(sortedmeans)
          Sortedmeans[Rv] = sortedmeans
          sortedgroups = [x[1] for x in sortedmeans]
          #for (m,c) in sortedmeans: print c,"%0.1f" % m
          candidate_divisions = [] # lists of condition names in 2 subgroups (lower and higher, but not necessarily in that order)
          for j in range(1,len(sortedmeans)):
            lowsubset,highsubset = sortedgroups[:j],sortedgroups[j:]
            alldistinct = True
            for group1 in lowsubset:
              for group2 in highsubset: alldistinct = alldistinct and reject[(group1,group2)]; # print group1,group2,reject[(group1,group2)]
            diff = sortedmeans[j][0]-sortedmeans[j-1][0] # could be negative if conds not sorted
            vals = lowsubset+['<<<--->>>']+highsubset+["%s" % alldistinct,"%0.1f" % diff]
            #print '\t'.join(vals)
            if diff<0: diff,lowsubset,highsubset = -diff,highsubset,lowsubset
            if alldistinct==True: candidate_divisions.append((diff,lowsubset,highsubset))
          if len(candidate_divisions)>0:
            candidate_divisions.sort(reverse=True)
            (diff,lowsubset,highsubset) = candidate_divisions[0]
            Divisions[Rv] = (lowsubset,highsubset)

    #The tukeyhsd of statsmodels doesn't return P value.
    #So, if you want to know P value, calculate from these outputted value or use R.
    #After saving the output into a variable res, you can get p-values by applying
    # psturng(np.abs(res.meandiffs / res.std_pairs), len(res.groupsunique), res.df_total), where psturng comes from from statsmodels.stats.libqsturng import psturng

      # save pvals as temp_anova.txt
      file = open(self.makefname("temp_anova.txt", self.rundir),"w")
      vals = "Rv Gene TAs".split()+Conditions.keys+"pval padj".split()
      file.write('\t'.join(vals)+EOL)
      for (start,end,Rv,gene,strand) in genes:
       if Rv in Means:
        m,s = numpy.mean(Means[Rv]),numpy.std(Means[Rv])
        scv = s/m
        vals = [Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in Means[Rv]]+["%f" % x for x in [pvals[Rv],qvals[Rv]]] # could append SCV, but beware of how multi_TnSeq.py counts hits in temp_anova.txt
        file.write('\t'.join(vals)+EOL)
      file.close()

      # write sorted means and high/low subgroups into temp_subgroups.txt
      file = open(self.makefname("temp_divisions.txt", self.rundir),"w")
      vals = "Rv Gene TAs".split()+Conditions.keys
      file.write('\t'.join(vals)+EOL)
      for (start,end,Rv,gene,strand) in genes:
       if Rv in qvals and qvals[Rv]<0.05:
        vals = [Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in Means[Rv]]
        vals += [x[1] for x in Sortedmeans[Rv]]+["%0.1f" % x[0] for x in Sortedmeans[Rv]]
        if Rv in Divisions:
          (lowsubset,highsubset) = Divisions[Rv]
          vals += [','.join(lowsubset),','.join(highsubset)]
        file.write('\t'.join(vals)+EOL)
      file.close()


      if False: # print gene means for each condition, and batch means
        for (start,end,Rv,gene,strand) in genes:
         vals = []
         if Rv in Means:
          batches = {}
          for r in range(Conditions.nrows):
            batch = Conditions.get(r,"Batch")
            if batch not in batches: batches[batch] = []
            batches[batch].append(Means[Rv][r])
            vals.append(Means[Rv][r]) ###
          for b in "CC1 CC2 CC3 KO".split(): vals.append(numpy.mean(batches[b]))
          print '\t'.join([Rv,gene]+["%0.1f" % x for x in vals]) ###
      ###sys.exit(0)



      # remove batch effects (subtract batch means from mean gene counts; adjust each to global mean)
      if self.debatch and "Batch" in Conditions.headers:
        print "<BR>correcting means for batch effects..."
        for (start,end,Rv,gene,strand) in genes:
         if Rv in Means:
          batches = {}
          for r in range(Conditions.nrows):
            batch = Conditions.get(r,"Batch")
            if batch not in batches: batches[batch] = []
            batches[batch].append(Means[Rv][r])

          globalmean = numpy.mean(Means[Rv])
          keys = sorted(batches.keys())
          batchmeans = {}
          for b in keys: batchmeans[b] = numpy.mean(batches[b])

          for r in range(Conditions.nrows):
            batch = Conditions.get(r,"Batch")
            delta = globalmean-batchmeans[batch] # correction for each batch
            Means[Rv][r] = max(0,Means[Rv][r]+delta)

      if False: # print gene means for each condition (with batch corrections)
        vals = ['ORF','gene']+Conditions.getcol("Condition")
        print '\t'.join(vals)
        for (start,end,Rv,gene,strand) in genes:
           if Rv not in Means: continue
           vals = [Means[Rv][r] for r in range(Conditions.nrows)]
           print '\t'.join([Rv,gene]+["%0.1f" % x for x in vals])
        sys.exit(0)

      ################################
      # calc LFCs

      # consider normalizing by reference conditions?
      # consider only calculating for non-reference conditions?

      lfcs = {}
      for (start,end,Rv,gene,strand) in genes:
       if Rv in qvals and qvals[Rv]<0.05:
        lfcvec = []
        for j in range(Conditions.nrows): # could remove means for reference conditions
          a = Means[Rv][j]
          b = numpy.median(Means[Rv])
          PC = 5
          lfc = numpy.log2((a+5)/float(b+5))
          lfcvec.append(lfc)
        lfcs[Rv] = lfcvec

      if self.stdnorm:
        print "<BR>applying standard normalization to LFCs"
        for i,cond in enumerate(Conditions.keys):
          col = [row[i] for row in lfcs.values()]
          m,s = numpy.mean(col),numpy.std(col) # consider using quantiles(col,[0,0.25,0.5,0.75,1.0])
          for Rv in lfcs.keys():
            lfcs[Rv][i] = (lfcs[Rv][i]-m)/s

      file = open(self.makefname("temp_LFCs.txt", self.rundir),"w")
      vals = "Rv Gene".split()+Conditions.keys # or non-ref-conds
      file.write('\t'.join(vals)+EOL)
      cnt = 0
      for (start,end,Rv,gene,strand) in genes:
       if Rv in qvals and qvals[Rv]<0.05:
        vals = [Rv,gene]+["%0.3f" % x for x in lfcs[Rv]]
        file.write('\t'.join([str(x) for x in vals])+EOL)
        cnt += 1
      file.close()
      if cnt==0: print "error: no significantly varying genes found by ANOVA"; return

      r_conditions_corrplot = def_r_conditions_corrplot()
      r_make_heatmap = def_r_make_heatmap()
      r_clustering = def_r_clustering()

      r_conditions_corrplot("temp_LFCs.txt")
      r_make_heatmap("temp_LFCs.txt")
      r_clustering("temp_LFCs.txt", self.K, self.D)

      self.generate_pdf()
      self.transit_message("Done.")

    def to_html(self):
        stats_table = '<table border="1" cellpadding="5" cellspacing="0">'
        for line in open("temp_stats.txt"):
          w = line.rstrip().split('\t')
          stats_table += "<tr>"
          for val in w:
              stats_table += "<td>" + val + "</td>"
          stats_table += "</tr>"
        stats_table += '</table>'
        lfcs_path = os.getcwd() + "/lfcs_boxplot.png"

        corr_plot = '<p><br>lfcs_boxplot.png:<br><img src="{0}"></p>'.format(lfcs_path)

        anova_varying = 0
        skip = 1
        for line in open("temp_anova.txt"):
            # skip header lines
            if skip > 0: skip -= 1; continue
            w = line.rstrip().split("\t")
            if w[-1] != 'nan' and float(w[-1]) < 0.05: anova_varying += 1

        anova_results = '<p>ANOVA finds {0} significantly varying genes'.format(anova_varying)
        conditions_corrplot_path = os.getcwd() + "/conditions_corrplot.png"
        conditions_corrplot = '<p><br>conditions_corrplot.png:<br><img src="{0}"></p>'.format(conditions_corrplot_path)
        lfcs_heatmap = "<p><br>heatmap.png:<BR><img src='{0}'></p>".format(os.getcwd() + "/heatmap.png")
        cluster_analysis = ''' '''

        return '''
            <p>
            Links
                <ul>
                    <li> <a href="temp_stats.txt">spreadsheet of statistics on TnSeq datasets</a>
                    <li> <a href="temp_counts_TTR.txt">spreadsheet of normalized counts</a>
                    <li> <a href="temp_anova.txt">spreadsheet with gene means and Anova results</a>
                    <li> <a href="temp_clust.txt">spreadsheet with LFCs and clusters</a>
                    <li> <a href="temp_divisions.txt">spreadsheet with divisions of conditions into low and high subgroups for significantly varying genes </a>
                    <li> <a href="temp_scores.txt"> spreadsheet with loadings of genes onto Varimax dimensions, along with normalized associations</a>
                    <li> <a href="PFilteredValues.txt"> text file with enriched functional categories in each clusters</a>
                </ul>
            </p>
            <hr>
            <H3>Quality Control</H3>
            <P>Statistics of TnSeq datasets:<BR>
            {stats_table}
            {corr_plot}
            <hr>
            {anova_results}
            <br>
            Conditions corrplot is based on significantly varying genes
            {conditions_corrplot}
            <hr>
            <h3>Cluster Analysis</h3>
            <p>LFCs are relative to the median across all conditions</p>
            {lfcs_heatmap}
        '''.format(
                stats_table=stats_table,
                corr_plot=corr_plot,
                anova_results=anova_results,
                conditions_corrplot=conditions_corrplot,
                lfcs_heatmap=lfcs_heatmap
            )
# print "<HR>"
# print "<H3>Cluster Analysis</H3>"

# print "<P>LFCs are relative to the median across all conditions"
# newfname = "%s/heatmap.png" % dir
# print "<BR>heatmap.png:<BR><img src='%s'>" % newfname

# newfname = "%s/cluster_opt.png" % dir
# print "<P>cluster_opt.png:<BR><img src='%s'>" % newfname

# newfname = "%s/pca_genes.png" % dir
# print "<P>LD-PCA, colored by Kmeans cluster (K=%s)" % K
# print "<BR>pca_genes.png:<BR><img src='%s'>" % newfname

# newfname = "%s/hclust_genes.png" % dir
# print "<P>hclust_genes.png:<BR><img src='%s'>" % newfname

# newfname = "%s/hclust_conditions.png" % dir
# print "<P>hclust_conditions.png:<BR><img src='%s'>" % newfname


# print "<HR>"

# print "<H3>Pathway Enrichment</H3>"

# print """This analysis shows which functional categories are enriched in each cluster.
# It uses the Sanger annotation (based on Cole's original descriptions of ORF functions,
# with updates.  Only the categories with <I>adjusted p-values</I> < 0.05 are shown below.
# P-values were calculated using the Hypergeometric distribution, and
# the multiple-test correction by the Benjamini-Hochberg procedure
# was applied to control the FDR at 5%."""

# from functionPathWay import sangerClassification
# sangerRoles="H37Rv.sanger_roles2"
# protTable="%s/ref.prot_table" % dir
# clusters="%s/temp_clust.txt" %dir
# resultFileName="%s/PFilteredValues.txt" % dir
# sangerClassification(sangerRoles,protTable,clusters,resultFileName)

# print "<PRE>"
# for line in open(resultFileName):
#   print line,
# print "</PRE>"

# print "<HR>"

# print "<H3>Principle Component Analysis</H3>"


# print "<P>Mapping of activators onto Principle Components"
# newfname = "%s/pca_conditions.png" % dir
# print "<BR>pca_conditions.png:<BR><img src='%s'>" % newfname

# newfname = "%s/condition_PCs.png" % dir
# print "<P>condition_PCs.png:<BR><img src='%s'>" % newfname

# newfname = "%s/varimax.png" % dir
# print "<P>varimax.png:<BR><img src='%s'>" % newfname
# print "</BODY></HTML>"

    def generate_pdf(self):
        html_content = self.to_html()
        file_name = os.path.splitext(os.path.basename(self.output_file))[0]
        with open("{0}.html".format(file_name), 'w') as f:
            f.write(html_content)
        pdfkit.from_string(html_content, "{0}.pdf".format(file_name))

    @classmethod
    def usage_string(self):
        return """
            python %s multitnseq <ref.fna> <ref.prot_table> <combined_wig> <samples_metadata> <conditions_metadata> <output_filename> [-K n] [-D n] [-debatch] [-stdnorm]
        """ % (sys.argv[0])

def main():
    print("MultiTnseq example.")

    return self()

if __name__ == "__main__":
    main()

