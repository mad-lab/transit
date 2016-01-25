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

#Ignores Divide by Zero warnings caused by calculations in log-space.
numpy.seterr(divide='ignore')

dehmm_prefix = "[DE-HMM]"

def cleanlog(X):
    if type(X) == type(numpy.array([1])):
        temp = numpy.zeros(len(X)) -23
        ii = X != float("-inf")
        temp[ii] = X[ii]
        return temp
    else:
        if X == float("-inf"):
            return -23
        else:
            return X




def runDEHMM(wx, pubmsg, **kwargs):

    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False

    print dehmm_prefix, "Running DE-HMM Method"
   
    expList = kwargs.get("expList")
    ctrlList = kwargs.get("ctrlList")
    annotationPath = kwargs.get("annotationPath")
    repchoice = kwargs.get("repchoice", "Sum")
    ignoreCodon = kwargs.get("ignoreCodon", True)
    ignoreNTerm = kwargs.get("ignoreNTerm", 0)
    ignoreCTerm = kwargs.get("ignoreCTerm", 0)
    doLOESS = kwargs.get("doLOESS", False)
    clusterPenalty = kwargs.get("clusterPenalty", 15.0)
    sitesPenalty = kwargs.get("sitesPenalty", 1.0)
    output = kwargs.get("output", sys.stdout)


    N1 = len(ctrlList)
    N2 = len(expList)

    orf2info = transit_tools.get_gene_info(annotationPath)

    hash = transit_tools.get_pos_hash(annotationPath)
    rv2name = transit_tools.get_gene_name(annotationPath)

    (data, position) = transit_tools.get_data(ctrlList + expList)
    factors = []

    if doLOESS:
        print dehmm_prefix, "Performing LOESS Correction"
        for j in range(len(data)):
            data[j] = transit_tools.loess_correction(position, data[j])


    print dehmm_prefix, "Normalizing with", "TTR"
    factors = transit_tools.TTR_factors(data)
    data = factors * data

    #O = numpy.round(numpy.sum(data,0))
    combined_data = numpy.zeros((3, len(data[0])))

    combined_data[0,:] = numpy.round(numpy.mean(data[:N1,:],0)) + 1
    combined_data[1,:] = numpy.round(numpy.mean(data[N1:,:],0)) + 1
    #combined_data[2,:] = numpy.round((combined_data[0,:]+combined_data[1,:])/2.0) + 1
    combined_data[2,:] = numpy.round((numpy.mean(data[:N1,:],0) +  numpy.mean(data[N1:,:],0))/2.0) + 1

    # TRI
    #combined_data[0,:] = numpy.round(numpy.sum(data[:N1,:],0)) + 1
    #combined_data[1,:] = numpy.round(numpy.sum(data[N1:,:],0)) + 1
    #combined_data[2,:] = numpy.round((numpy.sum(data[:N1,:],0) +  numpy.sum(data[N1:,:],0))/2.0) + 1

    ctrldata = combined_data[0,:]
    expdata = combined_data[1,:]
    jointdata = combined_data[2,:]
    
    #Control HMM
    print dehmm_prefix, "Running Control HMM"
    #(ctrl_likelihood_list, ctrl_gamma_list, ctrl_state_list, ctrl_params) = runHMM(ctrldata, wx, pubmsg, newWx)
    (ctrl_likelihood_list, ctrl_gamma_list, ctrl_state_list, ctrl_params,ctrl_mu,ctrl_labels) = runHMM(ctrldata, wx, pubmsg, newWx)

    #Exp HMM
    print dehmm_prefix, "Running Experimental HMM"
    #(exp_likelihood_list, exp_gamma_list, exp_state_list, exp_params) = runHMM(expdata, wx, pubmsg, newWx)
    (exp_likelihood_list, exp_gamma_list, exp_state_list, exp_params,exp_mu,exp_labels) = runHMM(expdata, wx, pubmsg, newWx)

    #Joint HMM
    print dehmm_prefix, "Running Joint HMM"
    #(joint_likelihood_list, joint_gamma_list, joint_state_list, joint_params,joint_mu) = runHMM(jointdata, wx, pubmsg, newWx)
    (joint_likelihood_list, joint_gamma_list, joint_state_list, joint_params,joint_mu,joint_labels) = runHMM(jointdata, wx, pubmsg, newWx)

    print dehmm_prefix, "Finished the Control, Experimental and Joint HMMs"

    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Calculating LLR...")
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Calculating LLR...")

    print dehmm_prefix, "Calculating LLR"

    # Get LLR
    T = len(ctrldata)
    L1 = numpy.zeros(T)
    L2 = numpy.zeros(T)
    LLR = numpy.zeros(T)
    for i in range(T):
        #print exp_gamma_list[i]
        #print ctrl_gamma_list[i]
        #print joint_gamma_list[i]
        #print ""
        indepA = cleanlog(numpy.log(numpy.sum(scipy.stats.geom.pmf(ctrldata[i], ctrl_params) *  ctrl_gamma_list[i])))
        indepB = cleanlog(numpy.log(numpy.sum(scipy.stats.geom.pmf(expdata[i], exp_params) *  exp_gamma_list[i])))
        combA = cleanlog(numpy.log(numpy.sum(scipy.stats.geom.pmf(ctrldata[i], joint_params) * joint_gamma_list[i])))
        combB = cleanlog(numpy.log(numpy.sum(scipy.stats.geom.pmf(expdata[i], joint_params) * joint_gamma_list[i])))

        comb = combA + combB
        indep = indepA + indepB
        L1[i] = indep
        L2[i] = comb
        #LLR[i] = comb - indep
        LLR[i] = indep - comb

    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Segmenting genome...")
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Segmenting genome...")
    segment_bits = segmentLLR(LLR, clusterPenalty, sitesPenalty)



    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Creating individual sites file...")
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Creating individual sites file...")
    output.write("#DE-HMM Sites\n")
    output.write("#Cluster Penalty: %1.2f\n" %  clusterPenalty)
    output.write("#Sites Penalty: %1.2f\n" % sitesPenalty)
    output.write("# Ctrl HMM Paramters: %s\n" % "\t".join(["%1.4f" % p for p in ctrl_params]))
    output.write("#   State means: "+' '.join(["%8.4f" % (p) for p,s in zip(ctrl_mu,ctrl_labels)])+"\n")
    output.write("# Exp HMM Paramters: %s\n" % "\t".join(["%1.4f" % p for p in exp_params]))
    output.write("#   State means: "+' '.join(["%8.4f" % (p) for p,s in zip(exp_mu,exp_labels)])+"\n")
    output.write("# Joint HMM Paramters: %s\n" % "\t".join(["%1.4f" % p for p in joint_params]))
    output.write("#   State means: "+' '.join(["%8.4f" % (p) for p,s in zip(joint_mu,joint_labels)])+"\n")
    for i in range(T):
        genes_at_site = hash.get(position[i], [""])
        genestr = ""
        if not (len(genes_at_site) == 1 and not genes_at_site[0]):
            genestr = ",".join(["%s_(%s)" % (g,rv2name.get(g, "-")) for g in genes_at_site])

        #output.write("# %6.4f %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f %6.4f | %6.4f %6.4f %6.4f %6.4f\n" % (ctrl_gamma_list[i][0],ctrl_gamma_list[i][1],ctrl_gamma_list[i][2],ctrl_gamma_list[i][3],exp_gamma_list[i][0],exp_gamma_list[i][1],exp_gamma_list[i][2],exp_gamma_list[i][3],joint_gamma_list[i][0],joint_gamma_list[i][1],joint_gamma_list[i][2],joint_gamma_list[i][3])) # TRI

        output.write("%7d\t%5d\t%5d\t%8.4f\t%8.4f\t%8.4f\t%s\t%s\t%s\t%5d\t%s\n" % (position[i], ctrldata[i]-1, expdata[i]-1, L1[i], L2[i], LLR[i], ctrl_state_list[i], exp_state_list[i], joint_state_list[i], segment_bits[i], genestr))


    output.close()
    
    data = {"path":output.name, "type":"DE-HMM - Sites", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
    if wx and newWx: wx.CallAfter(pubmsg, "file", data=data)
    if wx and not newWx: wx.CallAfter(pubmsg, "file", data)
    print dehmm_prefix, "Adding File:", output.name


    print dehmm_prefix, "Analyzing Segments"
    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Analyzing Segments...")
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Analyzing Segments...")


    segments_path = ".".join(output.name.split(".")[:-1]) + "_segments." + output.name.split(".")[-1]


    (gene2sites, segment_data) = parseSegments(output.name)
    outputdata = analyzeSegments(segment_data, gene2sites, orf2info)


    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Creating Segments output...")
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Creating Segments output...")

    output2 = open(segments_path, "w")
    output2.write("#DE-HMM - Segments\n")
    output2.write("#ID\tNum. sites\tstart\tend\tLength\tcount genes\tGenes\tsum1\tsum2\tlog2FC\tpval_resampling\tqval_resampling\n")
    for row in outputdata:
        output2.write("%d\t%d\t%d\t%d\t%d\t%d\t%s\t%1.1f\t%1.1f\t%1.2f\t%1.5f\t%1.5f\n" % tuple(row))
    output2.close()


    #[i, n, start, end, span, count_genes, ",".join(["%s_(%s)" % (g, orf2info.get(g, "-")[0]) for g in sorted(set(clean_genes)) if g]), sum1, sum2, log2FC, pval_resampling]

    data2 = {"path":segments_path, "type":"DE-HMM - Segments", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
    print dehmm_prefix, "Adding File:", segments_path
    if wx:
        if newWx:
            wx.CallAfter(pubmsg, "file", data=data2)
            wx.CallAfter(pubmsg, "dehmm", msg="Finished!")
            wx.CallAfter(pubmsg,"finish", msg="hmm")
        else:
            wx.CallAfter(pubmsg, "file", data2)
            wx.CallAfter(pubmsg, "dehmm", "Finished!")
            wx.CallAfter(pubmsg,"finish", "hmm")


    print dehmm_prefix, "Finished DE-HMM"


def runHMM(O, wx, pubmsg, newWx):

    SEVEN_STATE = True # TRI
    ####################################
    #Prepare parameters
    N = 4
    reads = O-1
    reads_nz = sorted(reads[reads !=0 ])
    size = len(reads_nz)
    mean_r = numpy.average(reads_nz[:int(0.95 * size)])
    
    
    mu = numpy.array([1/0.99, 0.01 * mean_r + 2,  mean_r, mean_r*5.0])
    L = 1.0/mu
    label= {0:"ES", 1:"GD", 2:"NE",3:"GA"}

    if SEVEN_STATE==True: # TRI
      N = 7
      mu = numpy.array([math.exp((i/4.0)*math.log(mean_r)) for i in range(7)]); mu[0] += 0.01
      #mu = numpy.array([1.01, 0.25*mean_r , 0.5*mean_r , 0.75*mean_r , mean_r , 4.0*mean_r , 10.0*mean_r ])
      L = 1.0/mu
      label = dict([(s, "S%d" % (s+1)) for s in range(0,N)])


    B = [] # Emission Probability Distributions
    for i in range(N):
        B.append(scipy.stats.geom(L[i]).pmf)


    pins = hmm_tools.calculate_pins(O-1)
    pins_obs = sum([1 for rd in O if rd >=2])/float(len(O))
    # IF Default Transition, calculate as specified in paper
    pnon = 1.0 - pins
    pnon_obs = 1.0 - pins_obs

    for r in range(100):
        #if pnon ** r < 0.01: break
        if pnon ** r < 0.001: break

    A = numpy.zeros((N,N))
    a = math.log1p(-B[int(N/2)](1)**r)
    b = r*math.log(B[int(N/2)](1)) + math.log(1.0/3)
    for i in range(N):
        A[i] = [b]*N
        A[i][i] = a

    # TRI
    if SEVEN_STATE==True:
      A = numpy.zeros((N,N))
      x = 0.001
      for i in range(N):
        for j in range(N):
          if i==j: A[i][j] = 1.0-x
          else: A[i][j] = x/float(N-1)
      # note - A is the transition matrix, but in log form
      A = numpy.log(A)


    PI = numpy.zeros(N) # Initial state distribution
    PI[0] = 0.7; PI[1:] = 0.3/(N-1);


    #############################################


    ###############
    ### VITERBI ###
    (Q_opt, delta, Q) = viterbi(A, B, PI, O, wx, pubmsg, newWx, scaling=True, discrete=False)
    ###############


    ##################
    ### ALPHA PASS ###
    (log_Prob_Obs, alpha, C) = forward_procedure(numpy.exp(A), B, PI, O, wx, pubmsg, newWx)
    ##################

    #################
    ### BETA PASS ###
    beta = backward_procedure(numpy.exp(A), B, PI, O, wx, pubmsg, newWx, C)
    #################

    T = len(O); total=0; state2count = dict.fromkeys(range(N),0)
    for t in xrange(T):
        state = Q_opt[t]
        state2count[state] +=1
        total+=1

    obs_list = []
    gamma_list = []
    likelihood_list = []
    state_list = []
    T = len(O)
    for t in xrange(T):
        gamma_t = (alpha[:,t] * beta[:,t])/numpy.sum(alpha[:,t] * beta[:,t])
        gamma_list.append(gamma_t)
        likelihood_t = [B[j](O[t]) for j in range(len(B))]
        likelihood_list.append(likelihood_t)
        state_t = label.get(int(Q_opt[t]), "Unknown State")
        state_list.append(state_t)


    return (likelihood_list, gamma_list, state_list, L, mu, label)




def forward_procedure(A, B, PI, O, wx, pubmsg, newWx):
    T = len(O)
    N = len(B)
    alpha = numpy.zeros((N,  T))
    C = numpy.zeros(T)

    alpha[:,0] = PI * [B[i](O[0]) for i in range(N)]

    C[0] = 1.0/numpy.sum(alpha[:,0])
    alpha[:,0] = C[0] * alpha[:,0]


    ITERATIONS = T*4
    count = 2*T
    for t in xrange(1, T):
        #B[i](O[:,t])  =>  numpy.prod(B[i](O[:,t]))
        #b_o = numpy.array([numpy.prod(B[i](O[:,t])) for i in range(N)])
        b_o = [B[i](O[t]) for i in range(N)]

        alpha[:,t] = numpy.dot(alpha[:,t-1], A) * b_o

        C[t] = numpy.nan_to_num(1.0/numpy.sum(alpha[:,t]))
        alpha[:,t] = numpy.nan_to_num(alpha[:,t] * C[t])

        if numpy.sum(alpha[:,t]) == 0:
            alpha[:,t] = 0.0000000000001

        if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1
        #print t, O[:,t], alpha[:,t]

    log_Prob_Obs = - (numpy.sum(numpy.log(C)))
    return(( log_Prob_Obs, alpha, C ))




def backward_procedure(A, B, PI, O, wx, pubmsg, newWx, C=None):

    N = len(B)
    T = len(O)
    beta = numpy.zeros((N,T))

    ITERATIONS = T*4
    count = 3*T

    beta[:,T-1] = 1.0
    if C!=None: beta[:,T-1] = beta[:,T-1] * C[T-1]


    for t in xrange(T-2, -1, -1):
        #B[i](O[:,t])  =>  numpy.prod(B[i](O[:,t]))
        #b_o = numpy.array([numpy.prod(B[i](O[:,t])) for i in range(N)])
        b_o = [B[i](O[t]) for i in range(N)]


        beta[:,t] = numpy.nan_to_num(numpy.dot(A, (b_o * beta[:,t+1] ) ))


        if sum(beta[:,t]) == 0:
            beta[:,t] = 0.0000000000001

        if C!=None:
            beta[:,t] = beta[:,t] * C[t]

        if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

        #print t, beta[:,t]

    return(beta)


def viterbi(A, B, PI, O, wx, pubmsg, newWx, scaling=True, discrete=False):
    N=len(B)
    T = len(O)
    delta = numpy.zeros((N, T))

    b_o = [B[i](O[0]) for i in range(N)]
    if scaling:
        delta[:,0] = numpy.log(PI) + numpy.log(b_o)

    Q = numpy.zeros((N, T))


    count = 0
    ITERATIONS = T*4
    for t in xrange(1, T):
        b_o = [B[i](O[t]) for i in range(N)]
        #nus = delta[:, t-1] + numpy.log(A)
        nus = delta[:, t-1] + A
        delta[:,t] = nus.max(1) + numpy.log(b_o)
        Q[:,t] = nus.argmax(1)
        if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

    Q_opt = [numpy.argmax(delta[:,T-1])]
    for t in xrange(T-2, -1, -1):
        Q_opt.insert(0, Q[Q_opt[0],t+1])
        if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Running DE-HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))

        count+=1

    if wx and newWx: wx.CallAfter(pubmsg, "dehmm", msg="Running DE-HMM Method... %2.0f%%" % (100.0*(ITERATIONS)/(ITERATIONS)))
    if wx and not newWx: wx.CallAfter(pubmsg, "dehmm", "Running DE-HMM Method... %2.0f%%" % (100.0*(ITERATIONS)/(ITERATIONS)))

    return((Q_opt, delta, Q))







class Results:
    def __init__(self,S):
        self.S = S # the (sub-)sequence of values
        self.B = [[0,0],[0,0]] # bit-strings
        self.V = [[0,0],[0,0]] # scores

def select(S,mu):
    vals = []
    for i in range(len(mu)):
        if mu[i]==1: vals.append(S[i])
    return vals

def cost(mu,p,q):
    r,n = 0,0
    for i in range(len(mu)):
        if mu[i]==1:
            n += 1
            if i==0 or mu[i-1]==0: r += 1
    return p*r+q*n

def all_bit_vectors(n):
    if n==1: return [[0],[1]]
    sub = all_bit_vectors(n-1)
    return [x+[0] for x in sub]+[x+[1] for x in sub]



def calc_score(S,mu,p,q):
    return sum(select(S,mu))-cost(mu,p,q)

def argmax(choices,vals):
    return choices[vals.index(max(vals))]

def int_list(string):
    return [int(x) for x in list(string)]


def opt_seg(S,p,q):
    if len(S)==2:
        results = Results(S)
        for i in range(2):
            for j in range(2):
                results.B[i][j] = [i,j]
                results.V[i][j] = calc_score(S,results.B[i][j],p,q)
        return results

    if len(S)==3:
        results = Results(S)
        for i in range(2):
            for j in range(2):
                results.B[i][j] = [i,0,j]
                results.V[i][j] = calc_score(S,results.B[i][j],p,q)
                b2 = [i,1,j]; s2 = calc_score(S,b2,p,q)
                if s2>results.V[i][j]: results.B[i][j],results.V[i][j] = b2,s2
        return results

    i = int(len(S)/2)
    L,R = S[:i],S[i:]
    Lres = opt_seg(L,p,q)
    Rres = opt_seg(R,p,q)
    results = Results(S)

    # S = [a...b]+[c...d]
    for a in range(2): # left end
        for d in range(2): # right end
            bestchoice,bestscore = None,None
            for b in range(2): # left side of connection in middle
                for c in range(2): # right side of connection in middle
                    bits = Lres.B[a][b]+Rres.B[c][d]
                    score = Lres.V[a][b]+Rres.V[c][d]
                    if b==1 and c==1: score += p # discount 1 cluster-init penalty by connecting 2 in middle
                    if bestscore==None or score>bestscore: bestscore,bestchoice = score,bits
            results.B[a][d],results.V[a][d] = bestchoice,bestscore
    return results



def segmentLLR(S, p, q):

    results = opt_seg(S,p,q)
    choices,scores = [],[]
    for i in range(2):
        for j in range(2):
            choices.append(results.B[i][j]); scores.append(results.V[i][j])
    mu = argmax(choices,scores)
    #vals = select(S,mu); s,c = sum(vals),cost(mu)
    return mu
    



def resampling(reads, S=10000, N1=1, N2=1):
    count_utail = 0
    count_ltail = 0
    count_2tail = 0

    delta_sum_list = []
    s_performed = 0

    if len(reads) >0:
        sum1,sum2 = numpy.sum(reads[:,:N1]), numpy.sum(reads[:,N1:])
    else:
        sum1 = 0; sum2 = 0;

    try:
        log2FC = math.log(float(sum2/float(N2) if sum2 > 0 else 1)/float(sum1/float(N1) if sum1 > 0 else 1),2)
    except:
        log2FC = 0

    for s in range(S):
        #Reads
        if len(reads) >0:
            perm = numpy.random.permutation(reads.flatten()).reshape(reads.shape)
            sumA,sumB = numpy.sum(perm[:,:N1]), numpy.sum(perm[:,N1:])
        else:
            sumA = 0; sumB = 0;
        #all_perm=np.array((list(itertools.permutations([0,1,2,3]))))
        #(reads.flatten()[(all_perm[numpy.random.randint(0,len(all_perm),size=100)]+3*np.arange(100)[...,numpy.newaxis]).flatten()]).reshape(reads.shape)

        delta_sum_list.append(sumB-sumA)
        if sumB-sumA >= sum2-sum1: count_utail+=1
        if sumB-sumA <= sum2-sum1: count_ltail+=1
        if abs(sumB-sumA) >= abs(sum2-sum1): count_2tail+=1
        s_performed+=1

    pval_utail = count_utail/float(s_performed)
    pval_ltail = count_ltail/float(s_performed)
    pval_2tail = count_2tail/float(s_performed)

    return (sum1, sum2, log2FC, pval_2tail)


def parseSegments(path):
    temp = []
    segment_data = []
    segment_num = 0
    insegment = False
    gene2sites = {}
    genes_in_segment = set()
    for line in open(path):
        if line.startswith("#"): continue
        #60     0       0    -0.1    -0.2     0.0       0   Rv0001_(dnaA)
        tmp = line.split()
        if len(tmp) == 11:
            ta,rd1,rd2,ll1,ll2,llr,ctrlS,expS,jointS,flag,genestr = line.split()
        else:
            ta,rd1,rd2,ll1,ll2,llr,ctrlS,expS,jointS,flag = line.split()
            genestr = ""

        #print ta,rd1,rd2,ll1,ll2,llr,flag,gene
        genes_in_site = [g.split("_")[0] for g in genestr.split(",")]
        for gene in genes_in_site:
            if gene not in gene2sites: gene2sites[gene] = []
            gene2sites[gene].append([ta,flag])

    
    
    
        #new segment
        if flag == "1" and not insegment:
            temp = []
            genes_in_segment = set()
            insegment = True

    
    
        #within segment        
        if flag == "1" and insegment:
            genes_in_segment.update(genes_in_site)
            temp.append([ta, rd1,rd2,ll1,ll2,llr,flag, gene])
    
        #no segment
        if flag == "0" and not insegment:
            pass

        #ending segment
        if flag == "0" and insegment:

            genes_in_segment.discard("")
            #print genes_in_segment
            insegment = False
            temp.append([ta, rd1,rd2,ll1,ll2,llr,flag, gene])
            segment_data.append([temp, genes_in_segment])
            temp = []
            genes_in_segment = set()

    if temp != []:
        insegment = False
        temp.append([ta, rd1,rd2,ll1,ll2,llr,flag, gene])
        segment_data.append([temp, genes_in_segment])

    return (gene2sites, segment_data)



def fdr_corrected_pval(X):
    n = len(X)
    qvalues = numpy.zeros(n)
    pvalues = numpy.array(X)
    pvalues.sort()
    pvalues = pvalues[::-1]

    for i in xrange(n):
        rank = n - i
        qvalues[i] = n/float(rank) * pvalues[i]

    for i in xrange(0, n-1):
        if qvalues[i] < qvalues[i+1]:
            qvalues[i+1] = qvalues[i]

    p2qval = dict([(p,q) for (p,q) in zip(pvalues,qvalues)])
    return numpy.array([p2qval[p] for p in X])



def analyzeSegments(segment_data, gene2sites, orf2info):
    #print segment_data
    #print len(segment_data)
   
    data = []
    pval = [] 
    N = len(segment_data)
    for i in range(N):

        sites_data = segment_data[i][0]
        genes_in_segment = segment_data[i][1]
        n = len(sites_data)
        start = sites_data[0][0]
        end = sites_data[-1][0]
        span = int(end) - int(start)
        reads = numpy.zeros((n, 2))
        LLR = numpy.zeros(n)
        count_resampling = 0
        count_genes = 0
        for j, (ta, rd1,rd2,ll1,ll2,llr,flag,gene) in enumerate(sites_data):
            #gene = gene.strip()
            LLR[j] = float(llr)
            reads[j,:] = [float(rd1), float(rd2)]
            #if gene:
            #    genes_in_segment.add(gene)



        count_genes = 0
        clean_genes = set()
        #print i, genes_in_segment
        for gene in genes_in_segment:
            n_sites = len(gene2sites[gene])
            n_de = sum([1 for (pos,f) in gene2sites[gene] if f == "1"])
            if n_de < (n_sites/2.0): continue
            count_genes += 1
            clean_genes.add(gene)

        #Resampling    
        (sum1, sum2, log2FC, pval_resampling) = resampling(reads, S=10000, N1=1, N2=1)


        #Gamma
        #LLR_N = numpy.mean(LLR)
        #k = 2
        #pval_gamma = 1.0 - scipy.stats.gamma.cdf(LLR_N, n*k/2.0, scale=2.0/n)

        data.append([i, n, int(start), int(end), span, count_genes, ",".join(["%s_(%s)" % (g, orf2info.get(g, "-")[0]) for g in sorted(set(clean_genes)) if g]), sum1, sum2, log2FC, pval_resampling])
        pval.append(pval_resampling)

    qval = fdr_corrected_pval(pval)

    #print data
    #print pval

    return ([data[i] + [qval[i]] for i in range(len(data))])





