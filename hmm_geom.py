import sys
import numpy
import scipy.stats
import math
import hmm_tools

# Ignores Divide by Zero warnings caused by calculations in log-space.
numpy.seterr(divide='ignore')


def runHMM(wigPath, protPath, output, wx, pubmsg):


    defaultStates = True
    defaultParameters = True
    defaultTransition = True
    N = 4

    # Get Gene Annotation from Prot Table File
    hash = {}; rv2name = {}
    if protPath:
        hash = hmm_tools.hash_prot_genes(protPath)
        rv2name = hmm_tools.get_prot_names(protPath)



    # Read in read count data
    pos_list = []
    C = []
    for line in open(wigPath):
        if line.startswith("#"): continue
        if line.startswith("variable"): continue
        if line.startswith("location"): continue
        tmp = line.split()
        pos = int(tmp[0]); reads = int(tmp[1]);
        pos_list.append(int(pos))
        C.append(reads+1)

    O = numpy.array(C) # Reads



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
        b = r*math.log(B[int(N/2)](1)) + math.log(1.0/3)
        for i in range(N):
            A[i] = [b]*N
            A[i][i] = a

    else:
        A = numpy.log(PROB)
        r = -1

    PI = numpy.zeros(N) # Initial state distribution
    PI[0] = 0.7; PI[1:] = 0.3/(N-1);


    ###############
    ### VITERBI ###
    (Q_opt, delta, Q) = hmm_tools.viterbi(A, B, PI, O, wx, pubmsg, scaling=True, discrete=False)
    ###############

    

    wx.CallAfter(pubmsg, "hmm", msg="Creating HMM Sites output...")
    T = len(O); total=0; state2count = dict.fromkeys(range(N),0)
    for t in xrange(T):
        state = Q_opt[t]
        state2count[state] +=1
        total+=1


    output.write("# Tn-HMM\n")
    output.write("# Command Used: python %s\n" % " ".join(sys.argv))
    output.write("# \n") 
    output.write("# Mean:\t%2.2f\n" % (numpy.average(reads_nz)))
    output.write("# Median:\t%2.2f\n" % numpy.median(reads_nz))
    output.write("# pins (obs):\t%f\n" % pins_obs)
    output.write("# pins (est):\t%f\n" % pins)
    output.write("# Run length (r):\t%d\n" % r)
    output.write("# Self-Transition Prob:\n")
    output.write("#    %s\n" % "   ".join(["%s: %2.4e" % (label[i], A[i][i]) for i in range(N)]))
    output.write("# State Emission Parameters (theta):\n")
    output.write("#    %s\n" % "   ".join(["%s: %1.4f" % (label[i], L[i]) for i in range(N)]))
    output.write("# State Distributions:")
    output.write("#    %s\n" % "   ".join(["%s: %2.2f%%" % (label[i], state2count[i]*100.0/total) for i in range(N)]))


    T = len(O)
    last_orf = ""
    if not protPath: output.write("\n")
    for t in xrange(T):
        current_orf = hash.get(pos_list[t], "non-coding")
        if (not last_orf or last_orf != current_orf) and protPath:
            orf = hash.get(pos_list[t], "non-coding")
            output.write( "\n%s %s\n" % (orf, rv2name.get(orf,"-")))
            last_orf = current_orf
    
        output.write( "%-7d %-7d %-4s  %-4s\n" % (pos_list[t], O[t]-1, "   ".join( "%-9.2e" % B[i](O[t]) for i in range(N)), label.get(int(Q_opt[t]), "Unknown State")))



    output.close()



#                      Copyright 2013; Michael A. DeJesus, Thomas R. Ioerger.


