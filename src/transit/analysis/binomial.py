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
import transit.transit_tools as transit_tools
import transit.tnseq_tools as tnseq_tools
import transit.norm_tools as norm_tools
import transit.stat_tools as stat_tools

#method_name = "binomial"


############# GUI ELEMENTS ##################
def Hide(wxobj):
    wxobj.binomialPanel.Hide()

def Show(wxobj):
    wxobj.binomialPanel.Show()

def getInstructions():
        return """Instructions:

1. Make sure you have one control sample selected.
2. Modify the options as desired.
3. Click on the "Run binomial" button.
4. Choose a name for the output file.
5. Wait until the execution finishes and the output is added to the file list at the bottom of the screen.
                """



def getPanel(wxobj):
    wxobj.binomialPanel = wx.Panel( wxobj.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
    #wxobj.binomialPanel.SetMinSize( wx.Size( 50,1 ) )
    #wxobj.binomialPanel.SetMaxSize( wx.Size( 250,-1 ) )

    binomialSection = wx.BoxSizer( wx.VERTICAL )

    wxobj.binomialLabel = wx.StaticText( wxobj.binomialPanel, wx.ID_ANY, u"Binomial Options", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.binomialLabel.Wrap( -1 )
    binomialSection.Add( wxobj.binomialLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    binomialSizer1 = wx.BoxSizer( wx.HORIZONTAL )
    binomialSizer2 = wx.BoxSizer( wx.HORIZONTAL )
    binomialLabelSizer = wx.BoxSizer( wx.VERTICAL )
    binomialControlSizer = wx.BoxSizer( wx.VERTICAL )
    
    #Samples 
    wxobj.binomialSampleLabel = wx.StaticText( wxobj.binomialPanel, wx.ID_ANY, u"Samples", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.binomialSampleLabel.Wrap(-1)
    binomialLabelSizer.Add(wxobj.binomialSampleLabel, 1, wx.ALL, 5)
    wxobj.binomialSampleText = wx.TextCtrl( wxobj.binomialPanel, wx.ID_ANY, u"10000", wx.DefaultPosition, wx.DefaultSize, 0 )
    binomialControlSizer.Add(wxobj.binomialSampleText, 0, wx.ALL|wx.EXPAND, 5)

    #Burnin 
    wxobj.binomialBurnLabel = wx.StaticText( wxobj.binomialPanel, wx.ID_ANY, u"Burn-in", wx.DefaultPosition, wx.DefaultSize, 0 )
    wxobj.binomialBurnLabel.Wrap(-1)
    binomialLabelSizer.Add(wxobj.binomialBurnLabel, 1, wx.ALL, 5)
    wxobj.binomialBurnText = wx.TextCtrl( wxobj.binomialPanel, wx.ID_ANY, u"500", wx.DefaultPosition, wx.DefaultSize, 0 )
    binomialControlSizer.Add(wxobj.binomialBurnText, 0, wx.ALL|wx.EXPAND, 5)


    binomialSizer2.Add(binomialLabelSizer, 1, wx.EXPAND, 5)
    binomialSizer2.Add(binomialControlSizer, 1, wx.EXPAND, 5)
    binomialSizer1.Add(binomialSizer2, 1, wx.EXPAND, 5 )


    binomialSection.Add( binomialSizer1, 1, wx.EXPAND, 5 )

    wxobj.binomialButton = wx.Button( wxobj.binomialPanel, wx.ID_ANY, u"Run binomial", wx.DefaultPosition, wx.DefaultSize, 0 )
    binomialSection.Add( wxobj.binomialButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

    wxobj.binomialPanel.SetSizer( binomialSection )
    wxobj.binomialPanel.Layout()
    binomialSection.Fit( wxobj.binomialPanel )

    #Connect events
    wxobj.binomialButton.Bind( wx.EVT_BUTTON, wxobj.RunMethod )

    return wxobj.binomialPanel


def updateProgressBar(wxobj, count):
    wxobj.binomialProgress.SetValue(count)

def SetProgressRange(wxobj, X):
    wxobj.binomialProgress.SetRange(X)

def enableButton(wxobj):
    wxobj.binomialButton.Enable()



def getColumnNames():
    return ["Orf","Name","Description","Mean Insertion","Sites per Replicate","Total Insertions","Total Sites","thetabar", "zbar", "Call"]


def getFileHeaderText(path):
    ess=0; unc=0; non=0; short=0
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        if tmp[-1] == "Essential": ess+=1
        if tmp[-1] == "Uncertain": unc+=1
        if tmp[-1] == "Non-Essential": non+=1

    text = """Results:
    Essentials: %s
    Uncertain: %s
    Non-Essential: %s
        """ % (ess, unc, non)
    return text


FileTypes = {}
FileTypes["#Binomial"] = (transit_tools.getTabTableData, getColumnNames, [getFileHeaderText])


########## CLASS #######################

class Binomial(base.SingleConditionMethod):
    """   
    binomial
 
    """
    def __init__(self,
                ctrldata,
                annotation_path,
                output_file,
                samples=10000,
                burnin=500,
                replicates="Sum",
                normalization=None,
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0, wxobj=None):

        base.SingleConditionMethod.__init__(self, "binomial", "Binomial Method", "A basic binomial method to show how you could add your own method to TRANSIT.", ctrldata, annotation_path, output_file, replicates=replicates, normalization=normalization, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.samples = samples
        self.burnin = burnin


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
        samples = int(wxobj.binomialSampleText.GetValue())
        burnin = int(wxobj.binomialBurnText.GetValue())
        ignoreCodon = True
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        replicates = "Sum"
        normalization = None
        LOESS = False

        #Get output path
        name = transit_tools.basename(all_selected[0])
        defaultFileName = "binomial_output.dat"
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")



        return self(ctrldata,
                annotationPath,
                output_file,
                samples,
                burnin,
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus, wxobj)

    @classmethod
    def fromargs(self, rawargs): 
        (args, kwargs) = transit_tools.cleanargs(rawargs)

        print "ARGS:", args
        print "KWARGS:", kwargs

        ctrldata = args[0].split(",")
        annotationPath = args[1]
        outpath = args[2]
        output_file = open(outpath, "w")

        samples = int(kwargs.get("s", 10000))
        burnin = int(kwargs.get("b", 500))
        replicates = "Sum"
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
                replicates,
                normalization,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus)

    def Run(self):

        self.transit_message("Starting Binomial Method")
        start_time = time.time()
        
        self.progress_range(self.samples+self.burnin)

        #Get orf data
        self.transit_message("Getting Data")
        G = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus)


        #Parameters
        self.transit_message("Setting Parameters")
        w1 = 0.15
        w0 = 1.0 - w1
        ALPHA = 1
        BETA = 1
        ALPHA_w = 600
        BETA_w = 3400
        mu_c = 0

        Ngenes = len(G)
        sample_size = self.samples+self.burnin
        numReps = len(self.ctrldata)

        theta = numpy.zeros((Ngenes, sample_size))
        theta[:,0] = 0.10

        rho0 = numpy.zeros(sample_size); rho0[0] = 0.5;  Kp0 = numpy.zeros(sample_size); Kp0[0] = 10;
        rho1 = numpy.zeros(sample_size); rho1[0] = 0.10; Kp1 = numpy.zeros(sample_size); Kp1[0] = 3;

        pi0 = 0.5; M0 = 1;
        pi1 = 0.5; M1 = 1;

        a0 = 10; b0 = 1;
        a1 = 10; b1 = 1;

        Z = numpy.zeros((Ngenes, sample_size))
        pz1 = numpy.zeros(sample_size);
        n1 = 0
        alpha_w = 0.5; beta_w = 0.5;


        w1 = scipy.stats.beta.rvs(alpha_w, beta_w)
        W1 = numpy.zeros(sample_size); W1[0] = w1



        #
        self.transit_message("Setting Initial Values")
        K = numpy.array([sum([1 for x in gene.reads.flatten() if x> 0]) for gene in G])
        N = numpy.array([len(gene.reads.flatten()) for gene in G])

        for g,gene in enumerate(G):
            if N[g] == 0: theta[g][0] = 0.5
            elif K[g]/float(N[g]) == 0: theta[g][0] = 0.001
            elif K[g]/float(N[g]) == 1: theta[g][0] = 0.001
            else: theta[g][0] = K[g]/float(N[g])

            #print g, ORF[g], K[g], N[g], theta[g][0]
            Z[g][0] = scipy.stats.bernoulli.rvs(1-theta[g][0])


        acc_p0 = 0; acc_k0 = 0;
        acc_p1 = 0; acc_k1 = 0;


        rho0c_std = 0.010
        kp0c_std = 1.40
        rho1c_std = 0.009
        kp1c_std = 1.1


        for i in range(1, sample_size):

            i0 = Z[:,i-1] == 0; n0 = numpy.sum(i0);
            i1 = Z[:,i-1] == 1; n1 = numpy.sum(i1);

            theta[i0,i] = scipy.stats.beta.rvs(Kp0[i-1]*rho0[i-1] + K[i0],  Kp0[i-1]*(1-rho0[i-1]) + N[i0] - K[i0])
            theta[i1,i] = scipy.stats.beta.rvs(Kp1[i-1]*rho1[i-1] + K[i1],  Kp1[i-1]*(1-rho1[i-1]) + N[i1] - K[i1])
            
            rho0_c = rho0[i-1] + scipy.stats.norm.rvs(0, rho0c_std)
            Kp0_c = Kp0[i-1] + scipy.stats.norm.rvs(0, kp0c_std)


            if rho0_c <= 0: rho0[i] = rho0[i-1]
            else:
                fc = numpy.log(scipy.stats.beta.pdf(rho0_c, M0*pi0, M0*(1-pi0)))
                f0 = numpy.log(scipy.stats.beta.pdf(rho0[i-1], M0*pi0, M0*(1-pi0)))
                fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0_c, Kp0[i-1]*(1-rho0_c))))
                f0 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0[i-1], Kp0[i-1]*(1-rho0[i-1]))))
    
                if numpy.log(scipy.stats.uniform.rvs()) < fc - f0:
                    rho0[i] = rho0_c
                    acc_p0+=1
                else: rho0[i] = rho0[i-1]


            if Kp0_c <= 0: Kp0[i] = Kp0[i-1]
            else:
                fc = numpy.log(scipy.stats.gamma.pdf(Kp0_c, a0, b0));
                f0 = numpy.log(scipy.stats.gamma.pdf(Kp0[i-1], a0, b0));
                fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0_c*rho0[i], Kp0_c*(1-rho0[i]))))
                f0 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0[i], Kp0[i-1]*(1-rho0[i]))))
    
                if numpy.log(scipy.stats.uniform.rvs()) < fc - f0:
                    Kp0[i] = Kp0_c
                    acc_k0+=1
                else: Kp0[i] = Kp0[i-1]

            rho1_c = rho1[i-1] + scipy.stats.norm.rvs(0, rho1c_std)
            Kp1_c = Kp1[i-1] + scipy.stats.norm.rvs(0, kp1c_std)


            if rho1_c <= 0:
                rho1[i] = rho1[i-1]
            else:
                fc = numpy.log(scipy.stats.beta.pdf(rho1_c, M1*pi1, M1*(1-pi1)))
                f1 = numpy.log(scipy.stats.beta.pdf(rho1[i-1], M1*pi1, M1*(1-pi1)))
                fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1_c, Kp1[i-1]*(1-rho1_c))))
                f1 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1[i-1], Kp1[i-1]*(1-rho1[i-1]))))
    
                if numpy.log(scipy.stats.uniform.rvs()) < fc - f1:
                    rho1[i] = rho1_c
                    acc_p1+=1
                else: rho1[i] = rho1[i-1]

            if Kp1_c <= 0: Kp1[i] = Kp1[i-1]
            else:
                fc = numpy.log(scipy.stats.gamma.pdf(Kp1_c, a1, b1));
                f1 = numpy.log(scipy.stats.gamma.pdf(Kp1[i-1], a1, b1));
                fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1_c*rho1[i], Kp1_c*(1-rho1[i]))))
                f1 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1[i], Kp1[i-1]*(1-rho1[i]))))

                if numpy.log(scipy.stats.uniform.rvs()) < fc - f1:
                    Kp1[i] = Kp1_c
                    acc_k1+=1
                else: Kp1[i] = Kp1[i-1]


            g0 = scipy.stats.beta.pdf(theta[:,i], Kp0[i]*rho0[i], Kp0[i]*(1-rho0[i])) * (1-w1)
            g1 = scipy.stats.beta.pdf(theta[:,i], Kp1[i]*rho1[i], Kp1[i]*(1-rho1[i])) * (w1)
            p1 = g1/(g0+g1)
            p1 = numpy.nan_to_num(p1)

            
            try:
                Z[:,i] = scipy.stats.bernoulli.rvs(p1)
            except:
                inan = numpy.isnan(p1)
                print >> sys.stderr, "K=\t", K[inan]
                print >> sys.stderr, "N=\t", N[inan]
                print >> sys.stderr, "theta=", theta[inan,i]
                sys.exit()
            pz1[i] = p1[0]


            i1 = Z[:,i] == 1; n1 = numpy.sum(i1);
            #w1 = 0.15
            w1 = scipy.stats.beta.rvs(alpha_w + n1, beta_w + Ngenes - n1)
            W1[i] = w1


            #Update progress
            text = "Running Gumbel Method... %2.0f%%" % (100.0*(i+1)/(sample_size))
            self.progress_update(text, i)
            self.transit_message_inplace(text)


        z_bar = numpy.apply_along_axis(numpy.mean, 1, Z[:, self.burnin:])
        theta_bar = numpy.apply_along_axis(numpy.mean, 1, theta[:, self.burnin:])
        (ess_threshold, noness_threshold) = transit_tools.fdr_post_prob(z_bar)

        self.output.write("#Binomial\n")
        #output.write("#Command: %s\n" % " ".join(["%s=%s" %(key,val) for (key,val) in kwargs.items()]))
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: ctrldata=%s, annotation=%s, output=%s, samples=%s, burnin=%s\n" % (",".join(self.ctrldata), self.annotation_path, self.output, self.samples, self.burnin))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))

        self.output.write("#Thresholds: (%1.5f, %1.5f)\n" % (ess_threshold,noness_threshold))
        self.output.write("#rho0 Acceptance Rate:\t%f%%\n" % ((100.0*acc_p0)/sample_size))
        self.output.write("#Kp0  Acceptance Rate:\t%f%%\n" % ((100.0*acc_k0)/sample_size))
        self.output.write("#rho1 Acceptance Rate:\t%f%%\n" % ((100.0*acc_p1)/sample_size))
        self.output.write("#Kp1  Acceptance Rate:\t%f%%\n" % ((100.0*acc_k1)/sample_size))

        self.output.write("#%s\n" % "\t".join(getColumnNames()))

        data = []
        for g,gene in enumerate(G):
            c = "Uncertain"
            if z_bar[g] > ess_threshold:
                c = "Essential"
            if z_bar[g] < noness_threshold:
                c = "Non-Essential"    
            data.append("%s\t%s\t%s\t%1.1f\t%d\t%d\t%d\t%f\t%f\t%s" % (gene.orf, gene.name, gene.desc, K[g]/float(numReps), N[g]/numReps, K[g], N[g], theta_bar[g], z_bar[g], c))

        data.sort()
        for row in data:
            self.output.write("%s\n" % row)
        self.output.close()

        self.transit_message("") # Printing empty line to flush stdout 
        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file()
        self.finish()
        self.transit_message("Finished Binomial Method")


 

    @classmethod
    def usage_string(self):
        return """python %s binomial <comma-separated .wig files> <annotation .prot_table> <output file> [Optional Arguments]

        Optional Arguments:
            -s <int>        :=  Number of samples to take. Default: -s 10000
            -b <int>        :=  Number of burn-in samples to take. Default: -b 500
            -iN <float>     :=  Ignore TAs occuring at given fraction of the N terminus. Default: -iN 0.0
            -iC <float>     :=  Ignore TAs occuring at given fraction of the C terminus. Default: -iC 0.0
            """ % (sys.argv[0])


if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    print "ARGS:", args
    print "KWARGS:", kwargs

    G = Binomial.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")   
    G.print_members()

    print ""
    print "Running:"

    G.Run()


