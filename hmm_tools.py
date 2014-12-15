import sys
import numpy
import math
import scipy.stats


def forward_procedure(A, B, PI, O, scaling = True, discrete=False):
    T = len(O)
    N = len(B)
    alpha = numpy.zeros((N,  T))
    C = numpy.zeros(T)

    if scaling:
        C = numpy.zeros(T)

    alpha[:,0] = PI * [B[i](O[0]) for i in range(N)]

    C[0] = 1.0/numpy.sum(alpha[:,0])
    alpha[:,0] = C[0] * alpha[:,0]

    for t in xrange(1, T):
        b_o = [B[i](O[t]) for i in range(N)]
        alpha[:,t] = numpy.dot(alpha[:,t-1], A) * b_o

        C[t] = 1.0/numpy.sum(alpha[:,t])
        alpha[:,t] = alpha[:,t] * C[t]

    log_Prob_Obs = - (numpy.sum(numpy.log(C)))
    return(( log_Prob_Obs, alpha, C ))



def backward_procedure(A, B, PI, O, C=None, discrete=False):

    N = len(B)
    T = len(O)
    beta = numpy.zeros((N,T))


    beta[:,T-1] = 1.0
    if C!=None: beta[:,T-1] = beta[:,T-1] * C[T-1]


    for t in xrange(T-2, -1, -1):
        b_o = [B[i](O[t]) for i in range(N)]
        beta[:,t] = numpy.dot(A, (b_o * beta[:,t+1] ) )

        if C!=None:
            beta[:,t] = beta[:,t] * C[t]

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
    ITERATIONS = T + T -1
    for t in xrange(1, T):
        b_o = [B[i](O[t]) for i in range(N)]
        #nus = delta[:, t-1] + numpy.log(A)
        nus = delta[:, t-1] + A
        delta[:,t] = nus.max(1) + numpy.log(b_o)
        Q[:,t] = nus.argmax(1)
        wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

    Q_opt = [numpy.argmax(delta[:,T-1])]
    for t in xrange(T-2, -1, -1):
        Q_opt.insert(0, Q[Q_opt[0],t+1])
        wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(count-1)/(ITERATIONS)))
        count+=1

    wx.CallAfter(pubmsg, "hmm", msg="Running HMM Method... %2.0f%%" % (100.0*(ITERATIONS)/(ITERATIONS)))

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




def post_process_genes(path, annotationPath, output=sys.stdout):
    
    
    
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


    orf2desc = {}
    for line in open(annotationPath):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        orf2desc[tmp[8]] = tmp[0]


    for line in open(path):
        if line.startswith("#"): continue
        if line.startswith("\n"): continue

        if line.startswith("non-coding") or line[0] not in "0123456789":
            if gene != "": gene2stats.append((gene, name, reads_list, states_list))
            reads_list = []
            states_list = []
            run = 0
    
            if line.startswith("non-coding"):
                count+=1
                gene = "non-coding_%d" % count
                name = "-"
                continue
            else:
                gene = line.split()[0]
                name = line.split()[1]
                continue

        tmp = line.split()
        pos = int(tmp[0])
        read = int(tmp[1])
        state = tmp[-1]

        reads_list.append(read)
        states_list.append(state)
        data.append(read)


    output.write("#Orf\tName\tDesc\tN\tn0\tn1\tn2\tn3\tAvg. Insertions\tAvg. Reads\tState Call\n")

    for gene, name, reads, states in gene2stats:
        if gene.startswith("non-coding") and not NC: continue
        # Reads and Insertions
        n = len(reads)
        freq = numpy.average([c > 0 for c in reads])
        reads_nz = [c for c in reads if c > 0]
        avg_read_nz = 0
        if len(reads_nz) > 0:
            avg_read_nz = numpy.average(reads_nz)

        # States
        n = len(states)
        statedist = {}
        for st in states:
            if st not in statedist: statedist[st] = 0
            statedist[st] +=1

        # State counts
        n0 = statedist.get("ES", 0); n1 = statedist.get("GD", 0);
        n2 = statedist.get("NE", 0); n3 = statedist.get("GA", 0);

        theta = sum([1 for rd in data if rd > 0])/float(len(data))
        E = ExpectedRuns(n,   1.0 - theta)
        V = VarR(n,   1.0 - theta)
        if n0 == n: S = "ES"
        elif n0 >= int(E+(3*math.sqrt(V))): S = "ES"
        else: S = max([(statedist.get(s, 0), s) for s in ["ES", "GD", "NE", "GA"]])[1]

        output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.4f\t%1.1f\t%s\n" % (gene, name, orf2desc.get(gene,"-"), n, n0, n1, n2, n3, freq, avg_read_nz, S))

