import sys

try:
    import wx
    hasWx = True
    #Check if wx is the newest 3.0+ version:
    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False
except Exception as e:
    hasWx = False
    newWx = False

import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

import base
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "gumbel"
long_name = "Bayesian analysis of essentiality based on long gaps."
description = """Bayesian methods of analyzing longest runs of non-insertions in a row. Estimates the parameters using the MCMC sampling, and estimates posterior probabilities of essentiality. 

Reference: DeJesus et al. (2013; Bioinformatics)"""
transposons = ["himar1"]
columns = ["Orf","Name","Desc","k","n","r","s","zbar", "Call"]


############# Analysis Method ##############

class GumbelAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, description, transposons, GumbelMethod, GumbelGUI, [GumbelFile])


################## FILE ###################

class GumbelFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#Gumbel", columns)

    def getHeader(self, path):
        ess=0; unc=0; non=0; short=0
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[-1] == "E": ess+=1
            if tmp[-1] == "U": unc+=1
            if tmp[-1] == "NE": non+=1
            if tmp[-1] == "S": short+=1

        text = """Results:
    Essentials: %s
    Uncertain: %s
    Non-Essential: %s
    Short: %s
        """ % (ess, unc, non, short)
        return text




################# GUI ##################

class GumbelGUI(base.AnalysisGUI):


    def definePanel(self, wxobj):
        self.wxobj = wxobj
        gumbelPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        #wxobj.gumbelPanel.SetMinSize( wx.Size( 50,1 ) )
        #wxobj.gumbelPanel.SetMaxSize( wx.Size( 250,-1 ) )

        gumbelSection = wx.BoxSizer( wx.VERTICAL )

        gumbelLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelLabel.Wrap( -1 )
        gumbelSection.Add( gumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        mainSizer1 = wx.BoxSizer( wx.HORIZONTAL )

        mainSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        labelSizer = wx.BoxSizer( wx.VERTICAL )

        gumbelSampleLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelSampleLabel.Wrap( -1 )
        labelSizer.Add( gumbelSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
 
        gumbelBurninLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Burn-In", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelBurninLabel.Wrap( -1 )
        labelSizer.Add( gumbelBurninLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        gumbelTrimLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Trim", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelTrimLabel.Wrap( -1 )
        labelSizer.Add( gumbelTrimLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        gumbelReadLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelReadLabel.Wrap( -1 )
        labelSizer.Add( gumbelReadLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        gumbelRepLabel = wx.StaticText( gumbelPanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelRepLabel.Wrap( -1 )
        labelSizer.Add( gumbelRepLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

        mainSizer2.Add( labelSizer, 1, wx.EXPAND, 5 )

        widgetSizer = wx.BoxSizer( wx.VERTICAL )

        self.wxobj.gumbelSampleText = wx.TextCtrl( gumbelPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
        widgetSizer.Add( self.wxobj.gumbelSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        self.wxobj.gumbelBurninText = wx.TextCtrl( gumbelPanel, wx.ID_ANY, u"500", wx.DefaultPosition, wx.DefaultSize, 0 )
        widgetSizer.Add( self.wxobj.gumbelBurninText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        self.wxobj.gumbelTrimText = wx.TextCtrl( gumbelPanel, wx.ID_ANY, u"1", wx.DefaultPosition, wx.DefaultSize, 0 )
        widgetSizer.Add( self.wxobj.gumbelTrimText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        gumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
        self.wxobj.gumbelReadChoice = wx.Choice( gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelReadChoiceChoices, 0 )
        self.wxobj.gumbelReadChoice.SetSelection( 0 )
        widgetSizer.Add( self.wxobj.gumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

        gumbelRepChoiceChoices = [ u"Sum", u"Mean" ]
        self.wxobj.gumbelRepChoice = wx.Choice( gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelRepChoiceChoices, 0 )
        self.wxobj.gumbelRepChoice.SetSelection( 0 )
        widgetSizer.Add( self.wxobj.gumbelRepChoice, 0, wx.ALL|wx.EXPAND, 5 )

        mainSizer2.Add( widgetSizer, 1, wx.EXPAND, 5 )

        mainSizer1.Add( mainSizer2, 1, wx.EXPAND, 5 )

        gumbelSection.Add( mainSizer1, 1, wx.EXPAND, 5 )

        gumbelButton = wx.Button( gumbelPanel, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
        gumbelSection.Add( gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        gumbelPanel.SetSizer( gumbelSection )
        gumbelPanel.Layout()
        gumbelSection.Fit( gumbelPanel )


        #Connect events
        gumbelButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )

        self.panel = gumbelPanel






########## METHOD #######################

ALPHA = 1
BETA = 1
class GumbelMethod(base.SingleConditionMethod):
    """   
    Gumbel
 
    """
    def __init__(self,
                ctrldata,
                annotation_path,
                output_file,
                samples=10000,
                burnin=500,
                trim=1,
                minread=1,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, short_name, long_name, description, ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.samples = samples
        self.burnin = burnin
        self.trim = trim
        self.minread = minread
      

        self.cache_nn = {}
 


    @classmethod
    def fromGUI(self, wxobj):
        """ """

        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get selected files
        ctrldata = wxobj.ctrlSelected()
        if not transit_tools.validate_control_datasets(ctrldata):
            return None

        #Validate transposon types
        if not transit_tools.validate_filetypes(ctrldata, transposons):
            return None


        #Read the parameters from the wxPython widgets
        minread = int(wxobj.gumbelReadChoice.GetString(wxobj.gumbelReadChoice.GetCurrentSelection()))
        samples = int(wxobj.gumbelSampleText.GetValue())
        burnin = int(wxobj.gumbelBurninText.GetValue())
        trim = int(wxobj.gumbelTrimText.GetValue())
        replicates = wxobj.gumbelRepChoice.GetString(wxobj.gumbelRepChoice.GetCurrentSelection())
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(ctrldata[0])
        defaultFileName = "gumbel_%s_s%d_b%d_t%d.dat" % (".".join(name.split(".")[:-1]), samples, burnin, trim)
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(ctrldata,
                annotationPath,
                output_file,
                samples,
                burnin,
                trim,
                minread,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self,rawargs):

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")
    
        samples = int(kwargs.get("s", 10000))
        burnin = int(kwargs.get("b", 500))
        trim = int(kwargs.get("t", 1))
        minread = int(kwargs.get("m", 1))
        replicates = kwargs.get("r", "Sum")
        normalization = None
        LOESS = False
        ignoreCodon = True
        NTerminus = float(kwargs.get("iN", 0.0))
        CTerminus = float(kwargs.get("iC", 0.0))

        return self(ctrldata,
                annotationPath,
                output_file,
                samples,
                burnin,
                trim,
                minread,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.status_message("Starting Gumbel Method")

        #Set Default parameter values
        w1 = 0.15
        w0 = 1.0 - w1
        ALPHA = 1
        BETA = 1
        ALPHA_w = 600
        BETA_w = 3400
        mu_c = 0
        acctot = 0.0
        phi_start = 0.3
        sigma_c = 0.01 
        
        start_time = time.time()
       
        self.progress_range(self.samples+self.burnin)
        
        #Get orf data
        self.transit_message("Reading Annotation")
        self.transit_message("Getting Data")

        

        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, minread=self.minread, reps=self.replicates, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)

        ii_good = numpy.array([self.good_orf(g) for g in G]) # Gets index of the genes that can be analyzed

        K = G.local_insertions()[ii_good]
        N = G.local_sites()[ii_good]
        R = G.local_runs()[ii_good]
        S = G.local_gap_span()[ii_good]
        T = G.local_gene_span()[ii_good]

        self.transit_message("Doing Regression")
        mu_s, temp, sigma_s = stat_tools.regress(R, S) # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = stat_tools.regress(S, R) # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(G)
        N_GOOD = sum(ii_good)

        self.transit_message("Setting Initial Class")
        Z_sample = numpy.zeros((N_GOOD, self.samples))
        Z = [self.classify(g.n, g.r, 0.5)   for g in G if self.good_orf(g)]
        Z_sample[:,0] = Z
        N_ESS = numpy.sum(Z_sample[:,0] == 1)
        
        phi_sample = numpy.zeros(self.samples) #[]
        phi_sample[0] = phi_start
        phi_old = phi_start
        phi_new = 0.00
        
        SIG = numpy.array([self.sigmoid(g.s, g.t) * scipy.stats.norm.pdf(g.r, mu_r*g.s, sigma_r) for g in G if self.good_orf(g)])


        i = 1; count = 0;
        while i < self.samples:

            # PHI
            acc = 1.0
            phi_new  = phi_old + random.gauss(mu_c, sigma_c)
            i0 = Z_sample[:,i-1] == 0
            if phi_new > 1 or phi_new <= 0 or (self.F_non(phi_new, N[i0], R[i0]) - self.F_non(phi_old, N[i0], R[i0])) < math.log(random.uniform(0,1)):
                phi_new = phi_old
                acc = 0.0
                flag = 0
            
            # Z
            Z = self.sample_Z(phi_new, w1, N, R, S, T, mu_s, sigma_s, SIG)
            
            # w1
            N_ESS = sum(Z == 1)
            w1 = scipy.stats.beta.rvs(N_ESS + ALPHA_w, N_GOOD - N_ESS + BETA_w)
            
            count +=1
            acctot+=acc
            
            if (count > self.burnin) and (count % self.trim == 0):
                phi_sample[i] = phi_new
                Z_sample[:,i] = Z
                i+=1
            
            phi_old = phi_new
            #Update progress
            text = "Running Gumbel Method... %2.0f%%" % (100.0*(count+1)/(self.samples+self.burnin))
            self.progress_update(text, count)
            self.transit_message_inplace(text)


        ZBAR = numpy.apply_along_axis(numpy.mean, 1, Z_sample)
        (ess_t, non_t) = stat_tools.bayesian_ess_thresholds(ZBAR)

        #Orf    k   n   r   s   zbar
        self.output.write("#Gumbel\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s, samples=%s, minread=%s, trim=%s\n" % (",".join(self.ctrldata), self.annotation_path, self.output, self.samples, self.minread, self.trim))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))

        self.output.write("#Data: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#Annotation path: %s\n" % (",".join(self.ctrldata))) 
        self.output.write("#FDR Corrected thresholds: %f, %f\n" % (ess_t, non_t))
        self.output.write("#MH Acceptance-Rate:\t%2.2f%%\n" % (100.0*acctot/count))
        self.output.write("#Total Iterations Performed:\t%d\n" % count)
        self.output.write("#Sample Size:\t%d\n" % i)
        self.output.write("#phi estimate:\t%f\n" % numpy.average(phi_sample))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        self.output.write("#%s\n" % "\t".join(columns))
        i = 0
        data = []
        for g in G:
            if not self.good_orf(g):
                zbar = -1.0
            else:
                zbar = ZBAR[i]
                i+=1
            if zbar > ess_t:
                call = "E"
            elif non_t <= zbar <= ess_t:
                call = "U"
            elif 0 <= zbar < non_t:
                call = "NE"
            else:
                call = "S"
            data.append("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s\n" % (g.orf, g.name, g.desc, g.k, g.n, g.r, g.s, zbar, call))
        data.sort()
        for line in data:
            self.output.write(line)
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="Gumbel")
        self.finish()
        self.transit_message("Finished Gumbel Method") 

    @classmethod
    def usage_string(self):
        return """python %s gumbel <comma-separated .wig files> <annotation .prot_table> <output file> [Optional Arguments]
    
        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -b <integer>    :=  Number of Burn-in samples. Default -b 500
        -m <integer>    :=  Smallest read-count to consider. Default: -m 1
        -t <integer>    :=  Trims all but every t-th value. Default: -t 1
        -r <string>     :=  How to handle replicates. Sum or Mean. Default: -r Sum
        -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
        -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
        """ % (sys.argv[0])

    def good_orf(self, gene):
        return (gene.n >= 3 and gene.t >= 150)

    def classify(self, n,r,p):
        if n == 0: return 0
        q = 1-p; B = 1/math.log(1/p); u = math.log(n*q,1/p);
        pval = 1 - scipy.stats.gumbel_r.cdf(r,u,B)
        if pval < 0.05: return(1)
        else: return(0)

    def F_non(self, p, N, R):
        q = 1.0 - p
        total = numpy.log(scipy.stats.beta.pdf(p,ALPHA,BETA))
        mu = numpy.log(N*q) / numpy.log(1/p)
        sigma = 1/math.log(1/p);
        total+= numpy.sum(numpy.log(scipy.stats.gumbel_r.pdf(R, mu, sigma)))
        return(total)
    
    def sample_Z(self, p, w1, N, R, S, T, mu_s, sigma_s, SIG):
        G = len(N)
        q = 1.0-p
        mu = numpy.log(N*q) / numpy.log(1.0/p)
        sigma = 1.0/math.log(1.0/p);
        h0 = ((scipy.stats.gumbel_r.pdf(R,mu,sigma)) * scipy.stats.norm.pdf(S, mu_s*R, sigma_s)  * (1-w1))
        h1 = SIG * w1
        p_z1 = h1/(h0+h1)
        return scipy.stats.binom.rvs(1, p_z1, size=G)

    def sigmoid(self,d,n):
        Kn = 0.1
        MEAN_DOMAIN_SPAN = 300

        if d == 0: return(0.00)
        f = 1./(1.+math.exp(Kn*(MEAN_DOMAIN_SPAN-d)))
        #if n in self.cache_nn: return f/self.cache_nn[n]
        tot = 0
        N = int(n+1)
        for i in range(1,N): tot += 1.0/(1.0+math.exp(Kn*(MEAN_DOMAIN_SPAN-i)))
        self.cache_nn[n] = tot
        return f/tot





if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)


    G = GumbelMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


