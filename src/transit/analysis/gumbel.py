import sys
import wx
import os
import time
import math
import random
import numpy
import scipy.stats
import datetime

import base
import transit_tools

import tnseq_tools
import norm_tools
import stat_tools

method_name = "Gumbel"


############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.gumbelPanel.Hide()

def Show(wxobj):
    wxobj.gumbelPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run Gumbel" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.gumbelPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    wxobj.gumbelPanel.SetMinSize( wx.Size( 50,1 ) )
    wxobj.gumbelPanel.SetMaxSize( wx.Size( 250,-1 ) )

    gumbelSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Gumbel Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelLabel.Wrap( -1 )
    gumbelSection.Add( wxobj.gumbelLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    bSizer14 = wx.BoxSizer( wx.HORIZONTAL )

    bSizer15 = wx.BoxSizer( wx.HORIZONTAL )

    bSizer16 = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelSampleLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelSampleLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelSampleLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelBurninLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Burn-In", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelBurninLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelBurninLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelTrimLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Trim", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelTrimLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelTrimLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelReadLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Minimum Read", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelReadLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelReadLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

    wxobj.gumbelRepLabel = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Replicates", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.gumbelRepLabel.Wrap( -1 )
    bSizer16.Add( wxobj.gumbelRepLabel, 1, wx.ALL, 5 )

    bSizer15.Add( bSizer16, 1, wx.EXPAND, 5 )

    bSizer17 = wx.BoxSizer( wx.VERTICAL )

    wxobj.gumbelSampleText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelSampleText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    wxobj.gumbelBurninText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"500", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelBurninText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    wxobj.gumbelTrimText = wx.TextCtrl( wxobj.gumbelPanel, wx.ID_ANY, u"1", wx.DefaultPosition, wx.DefaultSize, 0 )
    bSizer17.Add( wxobj.gumbelTrimText, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    gumbelReadChoiceChoices = [ u"1", u"2", u"3", u"4", u"5" ]
    wxobj.gumbelReadChoice = wx.Choice( wxobj.gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelReadChoiceChoices, 0 )
    wxobj.gumbelReadChoice.SetSelection( 0 )
    bSizer17.Add( wxobj.gumbelReadChoice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )

    gumbelRepChoiceChoices = [ u"Sum", u"Mean" ]
    wxobj.gumbelRepChoice = wx.Choice( wxobj.gumbelPanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gumbelRepChoiceChoices, 0 )
    wxobj.gumbelRepChoice.SetSelection( 0 )
    bSizer17.Add( wxobj.gumbelRepChoice, 0, wx.ALL|wx.EXPAND, 5 )


    bSizer15.Add( bSizer17, 1, wx.EXPAND, 5 )


    bSizer14.Add( bSizer15, 1, wx.EXPAND, 5 )
        
    gumbelSection.Add( bSizer14, 1, wx.EXPAND, 5 )

    wxobj.gumbelButton = wx.Button( wxobj.gumbelPanel, wx.ID_ANY, u"Run Gumbel", wx.DefaultPosition, wx.DefaultSize, 0 )
    gumbelSection.Add( wxobj.gumbelButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.progressLabel1 = wx.StaticText( wxobj.gumbelPanel, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.progressLabel1.Wrap( -1 )
    gumbelSection.Add( wxobj.progressLabel1, 0, wx.ALL, 5 )

    wxobj.gumbelProgress = wx.Gauge( wxobj.gumbelPanel, wx.ID_ANY, 20, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL|wx.GA_SMOOTH )
    wxobj.gumbelProgress.SetValue( 0 )
    gumbelSection.Add( wxobj.gumbelProgress, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


    wxobj.gumbelPanel.SetSizer( gumbelSection )
    wxobj.gumbelPanel.Layout()
    gumbelSection.Fit( wxobj.gumbelPanel )


    #Connect events
    wxobj.gumbelButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.gumbelPanel


def updateProgressBar(wxobj, count):
    wxobj.gumbelProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.gumbelProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.gumbelButton.Enable()





########## CLASS #######################

ALPHA = 1
BETA = 1
class Gumbel(base.SingleConditionMethod):
    """   
    Gumbel
 
    """
    def __init__(self,
                output_file,
                annotation_path,
                ctrldata,
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

        base.SingleConditionMethod.__init__(self, "Gumbel", "Gumbel Method", "Gumbel method from DeJesus et al. (Bioinformatics, 2013)", output_file, annotation_path, ctrldata, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)
        self.samples = samples
        self.burnin = burnin
        self.trim = trim
        self.minread = minread
      

        self.cache_nn = {}
 


    @classmethod
    def fromGUI(self, wxobj):
        """ """

        #Get selected files
        all_selected = wxobj.ctrlSelected()
        if len(all_selected) ==0:
            wxobj.ShowError("Error: No dataset selected.")
            return None

        #Get Annotation file
        annotationPath = wxobj.annotationFilePicker.GetPath()
        if not annotationPath:
            wxobj.ShowError("Error: No annotation file selected.")
            return None


        #Read the parameters from the wxPython widgets
        ctrldata = all_selected
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
        name = transit_tools.basename(all_selected[0])
        defaultFileName = "gumbel_%s_s%d_b%d_t%d.dat" % (".".join(name.split(".")[:-1]), samples, burnin, trim)
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(output_file,
                annotationPath,
                ctrldata,
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
    def fromconsole(self,
                output_path,
                annotation_path,
                ctrldata,
                samples=10000,
                burnin=500,
                trim=1,
                minread=1,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0):

        return self(ctrldata,
                outpath,
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
        
        #Get orf data
        self.status_message("Reading Annotation")
        self.status_message("Getting Data")

        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, minread=self.minread, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)

        ii_good = numpy.array([self.good_orf(g) for g in G]) # Gets index of the genes that can be analyzed

        K = G.local_insertions()[ii_good]
        N = G.local_sites()[ii_good]
        R = G.local_runs()[ii_good]
        S = G.local_gap_span()[ii_good]
        T = G.local_gene_span()[ii_good]

        

        self.status_message("Doing Regression")
        mu_s, temp, sigma_s = stat_tools.regress(R, S) # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = stat_tools.regress(S, R) # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(G)
        N_GOOD = sum(ii_good)

        self.status_message("Setting Initial Class")
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
            #i0 = numpy.logical_and(Z_sample[:,i-1] == 0, ii_good)
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
        


        ZBAR = numpy.apply_along_axis(numpy.mean, 1, Z_sample)
        (ess_t, non_t) = stat_tools.bayesian_ess_thresholds(ZBAR)

        VERBOSE = False
        #Orf    k   n   r   s   zbar
        self.output.write("#Gumbel\n")
        #self.output.write("#Command: python transit.py %s\n" % " ".join(["%s=%s" %(key,val) for (key,val) in kwargs.items()]))
        self.output.write("#FDR Corrected thresholds: %f, %f\n" % (ess_t, non_t))
        self.output.write("#MH Acceptance-Rate:\t%2.2f%%\n" % (100.0*acctot/count))
        self.output.write("#Total Iterations Performed:\t%d\n" % count)
        self.output.write("#Sample Size:\t%d\n" % i)
        self.output.write("#phi estimate:\t%f\n" % numpy.average(phi_sample))
        self.output.write("#Time: %s\n" % (time.time() - start_time))
        if VERBOSE: self.output.write("#Orf\tName\tDesc\tk\tn\tr\ts\tzbar\tCall\tSample\n")
        else: self.output.write("#Orf\tName\tDesc\tk\tn\tr\ts\tzbar\tCall\n")
        i = 0
        data = []
        for g in G:

            if not self.good_orf(g):
                zbar = -1.0
            else:
                zbar = ZBAR[i]
                i+=1
            sample_str = ""
            if VERBOSE: sample_str = "\t"+ ",".join(["%d" % x for x in Z_sample[i,:]])
            if zbar > ess_t:
                call = "E"
            elif non_t <= zbar <= ess_t:
                call = "U"
            elif 0 <= zbar < non_t:
                call = "NE"
            else:
                call = "S"
            data.append("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s%s\n" % (g.orf, g.name, g.desc, g.k, g.n, g.r, g.s, zbar, call, sample_str))
        data.sort()
        for line in data:
            self.output.write(line)
        self.output.close()

        self.status_message("") # Printing empty line to flush stdout 
        self.status_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.status_message("Finished Gumbel Method") 


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
        if n in self.cache_nn: return f/self.cache_nn[n]
        tot = 0
        N = int(n+1)
        for i in range(1,N): tot += 1.0/(1.0+math.exp(Kn*(MEAN_DOMAIN_SPAN-i)))
        self.cache_nn[n] = tot
        return f/tot






if __name__ == "__main__":


    
    G = Gumbel("results_gumbel_test.dat",
                "H37Rv.prot_table",
                ["glycerol_H37Rv_merged.wig"],
                samples=10000,
                burnin=500,
                trim=1,
                minread=1,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0)

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


