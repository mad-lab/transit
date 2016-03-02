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
        
        #Get orf data
        self.status_message("Reading Annotation") 
        orf2info = transit_tools.get_gene_info(self.annotation_path)
        hash = transit_tools.get_pos_hash(self.annotation_path)

        self.status_message("Getting Data")
        (data, position) = transit_tools.get_data(self.ctrldata)
        orf2reads, orf2pos = transit_tools.get_gene_reads(hash, data, position, orf2info, ignoreCodon=self.ignoreCodon, ignoreNTerm=self.NTerminus, ignoreCTerm=self.CTerminus, orf_list=orf2info.keys())

        start_time = time.time()
        (ORF_all, K_all, N_all, R_all, S_all, T_all) = self.get_orf_data(orf2reads, orf2pos, orf2info, self.minread)
        bad_orf_set = set([ORF_all[g] for g in xrange(len(N_all)) if not self.good_orf(N_all[g], T_all[g])]);
        bad_orf_set.add("Rvnr01");
        (ORF, K, N, R, S, T) = self.get_orf_data(orf2reads, orf2pos, orf2info, self.minread, bad=bad_orf_set)
        
        orf2name = {}; orf2desc = {}
        for line in open(self.annotation_path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            orf = tmp[8]; name = tmp[7]; desc = tmp[0];
            orf2name[orf] = name; orf2desc[orf] = desc

        self.status_message("Doing Regression")
        mu_s, temp, sigma_s = self.regress(R,S) # Linear regression to estimate mu_s, sigma_s for span data
        mu_r, temp, sigma_r = self.regress(S, R) # Linear regression to estimate mu_r, sigma_r for run data

        N_GENES = len(N)

        self.status_message("Setting Initial Class")
        Z_sample = numpy.zeros((N_GENES, self.samples))
        Z = [self.classify(N[g], R[g], 0.5)   for g in xrange(N_GENES)]
        Z_sample[:,0] = Z
        N_ESS = numpy.sum(Z_sample[:,0] == 1)
        
        phi_sample = numpy.zeros(self.samples) #[]
        phi_sample[0] = phi_start
        phi_old = phi_start
        phi_new = 0.00
        
        SIG = numpy.zeros(len(S))
        for g in range(len(S)):
            SIG[g] = self.sigmoid(S[g], T[g]) * scipy.stats.norm.pdf(R[g], mu_r*S[g], sigma_r)


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
            w1 = scipy.stats.beta.rvs(N_ESS + ALPHA_w, N_GENES - N_ESS + BETA_w)
            
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
        (ess_t, non_t) = transit_tools.fdr_post_prob(ZBAR)

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
        i = -1
        for g in xrange(len(ORF_all)):
            k = K_all[g]; n = N_all[g];  r = R_all[g]; s = S_all[g]; orf=ORF_all[g];
            if orf not in bad_orf_set:
                i+=1; zbar = ZBAR[i]; sample_str = "\t"+ ",".join(["%d" % x for x in Z_sample[i,:]]);
            else:
                zbar = -1.0; sample_str = "\t" + "-1";
            if not VERBOSE: sample_str = ""
        
            if zbar > ess_t: call = "E"
            elif non_t <= zbar <= ess_t: call = "U"
            elif 0 <= zbar < non_t: call = "NE"
            else: call = "S"

            self.output.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s%s\n" % (orf, orf2name.get(orf,"-"), orf2desc.get(orf,"-"), k, n, r, s, zbar, call, sample_str))

        self.output.close()


        self.status_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.status_message("Finished Gumbel Method") 








    def get_orf_data(self, orf2reads, orf2pos, orf2info, min_read, repchoice="Sum", bad=set()):
        G = len([orf for orf in orf2reads if orf not in bad]); g = 0;
        K = numpy.zeros(G); N = numpy.zeros(G); R = numpy.zeros(G);
        S = numpy.zeros(G); T = numpy.zeros(G); ORF = [];
        
        for orf in sorted(orf2reads):
            if orf in bad: continue
            run = [0,0,0]; maxrun = [0,0,0]; k = 0; n = 0; s = 0;
            start = orf2info.get(orf,[0,0,0,0])[2]; end = orf2info.get(orf,[0,0,0,0])[3]; length = end-start
            for i,rdrow in enumerate(orf2reads[orf]):
                pos = orf2pos[orf][i]
                if repchoice == "Sum":
                    rd = numpy.sum(rdrow)
                elif repchoice == "Mean":
                    rd = numpy.mean(rdrow)
                else:
                    rd = rdrow[0]
                
                n += 1
                if rd < min_read:
                    run[0] +=1
                    if run[0] == 1: run[1] = pos
                    run[2] = pos
                else:
                    k += 1
                    maxrun = max(run,maxrun)
                    run = [0,0,0]
            maxrun = max(run, maxrun)
            r = maxrun[0]
            if r > 0: s = maxrun[2] + 2  - maxrun[1]
            else: s = 0
            t = 0
            if orf2reads[orf]:
                t = max(orf2pos[orf]) + 2  - min(orf2pos[orf])
            
            K[g]=k; N[g]=n; R[g]=r; S[g]=s; T[g]=t;
            ORF.append(orf)
            g+=1
        return((ORF,K,N,R,S,T))


    def good_orf(self, n, s):
        return (n >= 3 and s >= 150)


    def regress(self, X,Y):
        N = len(X)
        xbar = numpy.average(X)
        ybar = numpy.average(Y)
        xybar = numpy.average([X[i]*Y[i] for i in range(N)])
        x2bar = numpy.average([X[i]*X[i] for i in range(N)])
        B = (xybar - xbar*ybar)/(x2bar - xbar*xbar)
        A0 = ybar - B*xbar
        
        yfit = [ A0 + B *X[i] for i in range(N)]
        yres = [Y[i] - (A0 + B *X[i]) for i in range(N)]
        var = sum([math.pow(yres[i],2) for i in range(N) ])/(N-2)
        std = math.sqrt(var)

        return(B, A0, std)


    def classify(self, n,r,p):
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
        G = len(N); h1 = numpy.zeros(G)
        q = 1-p
        mu = numpy.log(N*q) / numpy.log(1/p)
        sigma = 1/math.log(1/p);
        h0 = ((scipy.stats.gumbel_r.pdf(R,mu,sigma)) * scipy.stats.norm.pdf(S, mu_s*R, sigma_s)  * (1-w1))
        h1 = SIG * w1
        p_z1 = h1/(h0+h1)
        new_Z = scipy.stats.binom.rvs(1, p_z1, size=G)
        return(new_Z)



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


    
    G = Gumbel("gumbel_test.dat",
                "H37Rv.prot_table",
                ["glycerol_H37Rv_merged.wig"],
                samples=100,
                burnin=50,
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


