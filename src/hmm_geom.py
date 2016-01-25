# Copyright 2015.
#   Michael A. DeJesus, Chaitra Ambadipudi, and  Thomas R. Ioerger.
#
#
#    This file is part of TRANSIT.
#
#    TRANSIT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#
#    TRANSIT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.


import sys
import math
import datetime
import numpy
import scipy.stats
import hmm_tools
import transit_tools

# Ignores Divide by Zero warnings caused by calculations in log-space.
numpy.seterr(divide='ignore')

hmm_prefix = "[HMM]"

#def runHMM(wigPathList, protPath, repchoice, output, wx, pubmsg):
def runHMM(wx, pubmsg, **kwargs):

    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False


    
    print hmm_prefix, "Running HMM Method"


    readPathList = kwargs.get("readPathList")
    annotationPath = kwargs.get("annotationPath")
    repchoice = kwargs.get("repchoice", "Sum")
    ignoreCodon = kwargs.get("ignoreCodon", True)
    ignoreNTerm = kwargs.get("ignoreNTerm", 0)
    ignoreCTerm = kwargs.get("ignoreCTerm", 0)
    doLOESS = kwargs.get("doLOESS", False)
    output = kwargs.get("output", sys.stdout)
    SEVEN_STATE = False # TRI

    defaultStates = True
    defaultParameters = True
    defaultTransition = True
    N = 4
    if SEVEN_STATE==True: N = 7 # TRI

    # Get Gene Annotation from Prot Table File
    hash = {}; rv2name = {}
    if annotationPath:
        hash = transit_tools.get_pos_hash(annotationPath)
        rv2name = transit_tools.get_gene_name(annotationPath)

    # Read in read count data
    data = []
    for wigPath in readPathList:
        pos_list = []
        C = []
        for line in open(wigPath):
            if line.startswith("#"): continue
            if line.startswith("variable"): continue
            if line.startswith("location"): continue
            tmp = line.split()
            pos = int(tmp[0])
            reads = int(float(tmp[1]))
            pos_list.append(int(pos))
            C.append(reads)
        data.append(numpy.array(C))

    data = numpy.array(data)
    position = numpy.array(pos_list)

    if doLOESS:
        print hmm_prefix, "Performing LOESS Correction"
        for j in range(len(data)):
            data[j] = transit_tools.loess_correction(position, data[j])


    if repchoice == "Sum":
        O = numpy.round(numpy.sum(data,0))
    elif repchoice == "Mean":
        O = numpy.round(numpy.mean(data,0)) 
    else:
        O = data[0,:]

    #Add 1 to work with numpy's geometric function x : [1, inf]
    O = O + 1

    # If default states, set known labels
    if defaultStates: 
        label= {0:"ES", 1:"GD", 2:"NE",3:"GA"}
    else:
        label = dict([(s, "State-%d" % s) for s in range(0,N)])


    reads = O-1
    reads_nz = sorted(reads[reads !=0 ])
    size = len(reads_nz)
    mean_r = numpy.average(reads_nz[:int(0.95 * size)])


    # If parameters haven't changed, set them as specified in paper
    if defaultParameters:
        mu = numpy.array([1/0.99, 0.01 * mean_r + 2,  mean_r, mean_r*5.0])
        L = 1.0/mu
        label= {0:"ES", 1:"GD", 2:"NE",3:"GA"}
    else:
        L = PARAM

    # TRI
    if SEVEN_STATE==True:
      #mu = numpy.array([1/0.99, mean_r/64.0,mean_r/16.0,mean_r/4.0,  mean_r, mean_r*4.0,mean_r*16.0])
      #mu = numpy.array([1/0.99, mean_r/27.0,mean_r/9.0,mean_r/3.0,  mean_r, mean_r*3.0,mean_r*9.0])
      mu = numpy.array([math.exp((i/4.0)*math.log(mean_r)) for i in range(7)]); mu[0] += 0.01
      L = 1.0/mu
      label = dict([(s, "S%d" % (s+1)) for s in range(0,N)])

    B = [] # Emission Probability Distributions
    for i in range(N):
        B.append(scipy.stats.geom(L[i]).pmf)


    pins = hmm_tools.calculate_pins(O-1)
    pins_obs = sum([1 for rd in O if rd >=2])/float(len(O))
    # IF Default Transition, calculate as specified in paper
    if defaultTransition:
        pnon = 1.0 - pins
        pnon_obs = 1.0 - pins_obs

        for r in range(100):
            if pnon ** r < 0.01: break

        A = numpy.zeros((N,N))
        a = math.log1p(-B[int(N/2)](1)**r)
        b = r*math.log(B[int(N/2)](1)) + math.log(1.0/3) # change to N-1?
        for i in range(N):
            A[i] = [b]*N
            A[i][i] = a

        # TRI
        if SEVEN_STATE==True:
          A = numpy.zeros((N,N))
          x = 0.0001
          for i in range(N):
            for j in range(N):
              if i==j: A[i][j] = 1.0-x
              else: A[i][j] = x/float(N-1)
          # note - A is the transition matrix, but in log form
          A = numpy.log(A)

    else:
        A = numpy.log(PROB)
        r = -1

    PI = numpy.zeros(N) # Initial state distribution
    PI[0] = 0.7; PI[1:] = 0.3/(N-1);

    ###############
    ### VITERBI ###
    (Q_opt, delta, Q) = hmm_tools.viterbi(A, B, PI, O, wx, pubmsg, newWx, scaling=True, discrete=False)
    ###############


    ##################
    ### ALPHA PASS ###
    (log_Prob_Obs, alpha, C) = hmm_tools.forward_procedure(numpy.exp(A), B, PI, O, wx, pubmsg, newWx)

    #################
    ### BETA PASS ###
    beta = hmm_tools.backward_procedure(numpy.exp(A), B, PI, O, wx, pubmsg, newWx, C)


    if wx and newWx: wx.CallAfter(pubmsg, "hmm", msg="Creating HMM Sites output...")
    if wx and not newWx: wx.CallAfter(pubmsg, "hmm", "Creating HMM Sites output...")
    
    T = len(O); total=0; state2count = dict.fromkeys(range(N),0)
    for t in xrange(T):
        state = Q_opt[t]
        state2count[state] +=1
        total+=1


    output.write("#HMM - Sites\n")
    output.write("# Tn-HMM\n")
    #output.write("#Command: python transit.py %s\n" % " ".join(["%s=%s" %(key,val) for (key,val) in kwargs.items()]))
    if wx:
        variables = dict([(key,val) for (key,val) in kwargs.items() if not callable(val)])
        variables["output"] = output.name
        variables["readPathList"]= ",".join(readPathList)
        output.write("#GUI with: %s\n" % (" ".join(["%s=%s" % (key,val) for (key,val) in variables.items()])))
    else:
        output.write("#Console: python %s\n" % " ".join(sys.argv))

    output.write("# \n") 
    output.write("# Mean:\t%2.2f\n" % (numpy.average(reads_nz)))
    output.write("# Median:\t%2.2f\n" % numpy.median(reads_nz))
    output.write("# pins (obs):\t%f\n" % pins_obs)
    output.write("# pins (est):\t%f\n" % pins)
    output.write("# Run length (r):\t%d\n" % r)
    output.write("# State means:\n")
    output.write("#    %s\n" % "   ".join(["%s: %8.4f" % (label[i], mu[i]) for i in range(N)]))
    output.write("# Self-Transition Prob:\n")
    output.write("#    %s\n" % "   ".join(["%s: %2.4e" % (label[i], A[i][i]) for i in range(N)]))
    output.write("# State Emission Parameters (theta):\n")
    output.write("#    %s\n" % "   ".join(["%s: %1.4f" % (label[i], L[i]) for i in range(N)]))
    output.write("# State Distributions:")
    output.write("#    %s\n" % "   ".join(["%s: %2.2f%%" % (label[i], state2count[i]*100.0/total) for i in range(N)]))


    T = len(O)
    last_orf = ""
    for t in xrange(T):
        s_lab = label.get(int(Q_opt[t]), "Unknown State")
        gamma_t = (alpha[:,t] * beta[:,t])/numpy.sum(alpha[:,t] * beta[:,t])
        genes_at_site = hash.get(pos_list[t], [""])
        genestr = ""
        if not (len(genes_at_site) == 1 and not genes_at_site[0]):
            genestr = ",".join(["%s_(%s)" % (g,rv2name.get(g, "-")) for g in genes_at_site])

        output.write("%s\t%s\t%s\t%s\t%s\n" % (int(pos_list[t]), int(O[t])-1, "\t".join(["%-9.2e" % g for g in gamma_t]), s_lab, genestr))

    output.close()

    print hmm_prefix, "Finished HMM - Sites Method"
    if not output.name.startswith("<"):
        data = {"path":output.name, "type":"HMM - Sites", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}

        print hmm_prefix, "Adding File:", output.name

        if wx and newWx: wx.CallAfter(pubmsg, "file", data=data)
        if wx and not newWx: wx.CallAfter(pubmsg, "file", data)        

        if wx and newWx: wx.CallAfter(pubmsg, "hmm", msg="Creating HMM Genes output...")
        if wx and not newWx: wx.CallAfter(pubmsg, "hmm", "Creating HMM Genes output...")

        genes_path = ".".join(output.name.split(".")[:-1]) + "_genes." + output.name.split(".")[-1]
        hmm_tools.post_process_genes(output.name, annotationPath, ignoreCodon, ignoreNTerm, ignoreCTerm, output=open(genes_path,"w"))

        data2 = {"path":genes_path, "type":"HMM - Genes", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        #data["path"] =  genes_path
        #data["type"] = "HMM - Genes"
        print hmm_prefix, "Adding File:", genes_path

        if wx:
            if newWx:
                wx.CallAfter(pubmsg, "file", data=data2)
                wx.CallAfter(pubmsg, "hmm", msg="Finished!")
                wx.CallAfter(pubmsg,"finish", msg="hmm")
            else:
                wx.CallAfter(pubmsg, "file", data2)
                wx.CallAfter(pubmsg, "hmm", "Finished!")
                wx.CallAfter(pubmsg,"finish", "hmm")




#                      Copyright 2013; Michael A. DeJesus, Thomas R. Ioerger.


