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
import numpy
import math
import scipy.stats
import transit_tools



def forward_procedure(A, B, PI, O, wx, pubmsg):
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

        if wx: wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1
        #print t, O[:,t], alpha[:,t]

    log_Prob_Obs = - (numpy.sum(numpy.log(C)))
    return(( log_Prob_Obs, alpha, C ))



def backward_procedure(A, B, PI, O, wx, pubmsg, C=None):

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


        if wx: wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

        #print t, beta[:,t]

    return(beta)




def viterbi(A, B, PI, O, wx, pubmsg, scaling=True, discrete=False):
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
        if wx: wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

    Q_opt = [numpy.argmax(delta[:,T-1])]
    for t in xrange(T-2, -1, -1):
        Q_opt.insert(0, Q[Q_opt[0],t+1])
        if wx: wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

    if wx: wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(ITERATIONS)/(ITERATIONS)))

    return((Q_opt, delta, Q))


def hash_gff_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        if tmp[2] != "gene": continue
        features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
        start, end = int(tmp[3]), int(tmp[4])
        for i in range(start,end+1):
            if i not in hash:
                hash[i] = features.get("ID", "missing")
    return hash

def get_gff_names(path):
    rv2name = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        if tmp[2] != "gene": continue
        features = dict([tuple(f.split("=")) for f in tmp[8].strip().split(";")])
        rv2name[features.get("ID", "missing")] = features.get("Name", "-")
    return(rv2name)


def hash_prot_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        start, end = int(tmp[1]), int(tmp[2])
        for i in range(start,end+1):
            #if i not in hash:
            #    hash[i] = tmp[8]
            hash[i] = tmp[8]
    return hash


def get_prot_names(path):
    orf2name = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf2name[tmp[8]] = tmp[7]
    return(orf2name)




def calculate_pins(reads):

    non_ess_reads = []
    temp = []
    for rd in reads:

        if rd >=1:
            if len(temp) < 10: non_ess_reads.extend(temp)
            non_ess_reads.append(rd)
            temp = []
        else:
            temp.append(rd)

    return(sum([1 for rd in non_ess_reads if rd >= 1])/float(len(non_ess_reads)) )



### Gumbel/EVD Variance and Expected run funnctions
def VarR(n,p):
    """ VarR_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n) (Schilling, 1990)"""

    r2 = getR2(n)
    E2 = getE2(n)

    A = math.pow(math.pi,2.0)/(6* math.pow(math.log(1.0/p),2.0))
    V = A + 1/12.0 + r2 + E2

    return ( V )



def ExpectedRuns(n,p):
    """ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n) (Schilling, 1990)"""

    q = 1-p
    gamma = getGamma()
    r1 = getR1(n)
    E1 = getE1(n)


    A = math.log(n*q,1.0/p)
    B = gamma/math.log(1.0/p)
    ER = A + B -0.5 + r1 + E1

    return( ER )


def getGamma():
    """Euler-Mascheroni constant ~ 0.577215664901 """
    return(0.5772156649015328606)


def getR1(n):
    """Small Correction term. Defaults to 0.000016 for now"""
    return(0.000016)

def getR2(n):
    """Small Correction term. Defaults to 0.00006 for now"""
    return(0.00006)

def getE1(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return(0.01)

def getE2(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return(0.01)




def post_process_genes(path, annotationPath, ignoreCodon=True, ignoreNTerm=5, ignoreCTerm=5, output=sys.stdout):
    
    
    
    S = 0
    NC = False
    C = 2


    gene2stats = []
    count = 0
    reads_list = []
    states_list = []
    run = 0
    gene = ""
    data = []


    orf2info = transit_tools.get_gene_info(annotationPath)
    hash = transit_tools.get_pos_hash(annotationPath)


    data = [] 
    position = []
    reads = []
    pos2state = {}
    
    K = 1
    T = 0
    for line in open(path):
        if line.startswith("#"): continue
        if line.startswith("location"): continue
        if line.startswith("variable"): continue
        T+=1

    i = 0
    data = numpy.zeros((K,T))
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split()
        pos = int(tmp[0])
        rd = float(tmp[1])
        S = tmp[-2]
        data[0,i] = rd
        position.append(pos)
        pos2state[pos] = S
        i+=1

    #position = numpy.array(position)

    (orf2reads, orf2pos) = transit_tools.get_gene_reads(hash, data, position, orf2info, ignoreCodon=ignoreCodon, ignoreNTerm=ignoreNTerm, ignoreCTerm=ignoreCTerm, orf_list=orf2info.keys())


    theta = numpy.mean(data[0,:] > 0)

    output.write("#HMM - Genes\n")
    output.write("#Orf\tName\tDesc\tN\tn0\tn1\tn2\tn3\tAvg. Insertions\tAvg. Reads\tState Call\n")

    for gene in sorted(orf2info):
        if gene.startswith("non-coding") and not NC: continue

        name = orf2info[gene][0]
        desc = orf2info[gene][1]

        # Reads and Insertions
        reads = numpy.array(orf2reads.get(gene, []))
        n = len(reads)
        freq = numpy.average([c > 0 for c in reads])
        reads_nz = [c for c in reads if c > 0]

        avg_read_nz = 0
        if len(reads_nz) > 0:
            avg_read_nz = numpy.average(reads_nz)

        # State
        states = [pos2state[p] for p in orf2pos.get(gene,[])]
        n = len(states)
        statedist = {}
        for st in states:
            if st not in statedist: statedist[st] = 0
            statedist[st] +=1

        # State counts
        n0 = statedist.get("ES", 0); n1 = statedist.get("GD", 0);
        n2 = statedist.get("NE", 0); n3 = statedist.get("GA", 0);

       
         
        #print gene, name, reads, states, n0, n1, n2, n3, theta

        if n > 0:
            E = ExpectedRuns(n,   1.0 - theta)
            V = VarR(n,   1.0 - theta)
            if n0 == n: S = "ES"
            elif n0 >= int(E+(3*math.sqrt(V))): S = "ES"
            else: S = max([(statedist.get(s, 0), s) for s in ["ES", "GD", "NE", "GA"]])[1]
        else:
            E = 0.0
            V = 0.0
            S = "N/A"
        output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.4f\t%1.2f\t%s\n" % (gene, name, desc, n, n0, n1, n2, n3, freq, avg_read_nz, S))

    output.close()

