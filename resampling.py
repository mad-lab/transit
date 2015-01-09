import sys
import random


def sum_reads(reads, N1, N2):
    sum1 = 0
    sum2 = 0
    n = len(reads)
    if n > 0:
        for row in reads:
            sum1+=sum(row[:N1])
            sum2+=sum(row[N1:])
        return (sum1/float(N1*n),  sum2/float(N2*n))
    else:
        return(0,0)


def maxrun(lst,item=0):
    best = 0
    i,n = 0,len(lst)
    while i<n:
        if lst[i]==item:
            j = i+1
            while j<n and lst[j]==item: j += 1
            r = j-i
            if r>best: best = r
            i = j
        else: i += 1
    return best



def comb_runs(reads, N1, N2):
    if not reads: return (0,0)
    combined_col1 = [sum(row[:N1]) for row in reads]
    combined_col2 = [sum(row[N1:]) for row in reads]
    return (maxrun(combined_col1),  maxrun(combined_col2))
      
def fdr_corrected_pval(X):
    import numpy
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



def runResampling(ctrlString, expString, annotationPath, output, wx, pubmsg, doNormalize=True):

    FILES1 = ctrlString.split(",")
    FILES2 = expString.split(",")
    ANNOTATION = annotationPath

    orf2coords = {}
    orf2reads = {}
    #hash = {}
    #for line in open(ANNOTATION):
    #    if line.startswith("#"): continue
    #    if line.startswith("geneID"): continue
    #    tmp = line.strip().split("\t")
    #    orf = tmp[0]
    #    if orf == "intergenic": continue
    #    start = int(tmp[1])
    #    end = int(tmp[2])
    #    strand = tmp[3]
    #    TA_list = [int(ta) for ta in tmp[4:] if ta != "no TAs"]
    #    for pos in TA_list:
    #    #for pos in range(start+delta-2, end+1-delta):
    #        if pos not in hash: hash[pos] = []
    #        hash[pos].append(orf)
    #    orf2reads[orf] = []
    #    orf2coords[orf] = (start, end)


    hash = {}
    orf2reads = {}
    orf2data = {}
    for line in open(ANNOTATION):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        name = tmp[7]
        desc = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        for pos in range(start, end+1):
            if pos not in hash: hash[pos] = []
            hash[pos].append(orf)
        orf2reads[orf] = []
        orf2data[orf] = (name, desc)

    data = []
    position = []
    for i,path in enumerate(FILES1+FILES2):
        reads = []

        for line in open(path):
            if line.startswith("#"): continue
            if line.startswith("location"): continue
            if line.startswith("variable"): continue
            tmp = line.split()
            pos = int(tmp[0])
            rd = float(tmp[1])
            reads.append(rd)
            if i == 0: position.append(pos)
        data.append(reads)

    N1 = len(FILES1)
    N2 = len(FILES2)
    N = len(data)
    total_hits = [sum(x) for x in data]
    TAs_hit = [sum([1 for y in x if y > 0])  for x in data]
    mean_hits = [total_hits[i]/float(TAs_hit[i]) for i in range(N)]
    grand_total = sum(mean_hits)
    grand_mean = grand_total/float(N)
    factors = [grand_mean/float(mean_hits[i]) for i in range(N)]


    orf2TAs = {}
    for i,pos in enumerate(position):
        if pos not in hash: continue
        for gene in hash[pos]:
            if gene not in orf2TAs: orf2TAs[gene] = 0
            orf2TAs[gene]+=1
            if doNormalize:
                row = [data[j][i]*factors[j] for j in range(N)] #NORMALIZED
            else:
                row = [data[j][i] for j in range(N)] #NOT NORMALIZED
            orf2reads[gene].append(row)
    
    S = 10000

    output.write("#Resampling\n")
    output.write("#Command: python %s\n" % " ".join(sys.argv))
    output.write("#Control Samples:  %s\n" % ", ".join(FILES1))
    output.write("#Experimental Samples:  %s\n" % ", ".join(FILES2))
    if doNormalize:
        output.write("#Normalization factors: %s\n" % "\t".join(["%1.4f" % f for f in factors]))
    else:
        output.write("#Not normalized\n")

    output.write("#Orf\t%Name\tDescription\tN\tTAs Hit\tAvg Rd 1\tAvg Rd 2\tDelta Rd\tp-value\tp-adj\n")
    count = 0
    G = len(orf2reads)
    orf2out = {}
    pval=[]
    for orf in sorted(orf2reads):
    
        fullreads = orf2reads[orf]
        reads = [row for row in fullreads if sum(row) > 0 ]
        sum1,sum2 = sum_reads(reads, N1, N2)
        run1,run2 = comb_runs(fullreads, N1, N2)

        count_utail = 0
        count_ltail = 0
        count_2tail = 0
   
        for s in range(S):
            #Reads
            flat_reads = [item for sublist in reads for item in sublist]
            flat_perm = sorted(flat_reads, key=lambda k: random.random())
            perm = [flat_perm[ii*(N1+N2):(ii+1)*(N1+N2)]    for ii in range(len(flat_perm)/(N1+N2))]
            sumA,sumB = sum_reads(perm, N1, N2)

            if sumB-sumA >= sum2-sum1: count_utail+=1
            if sumB-sumA <= sum2-sum1: count_ltail+=1
            if abs(sumB-sumA) >= abs(sum2-sum1): count_2tail+=1

        #Reads - pvalues
        pval_utail = count_utail/float(S)
        pval_ltail = count_ltail/float(S)
        pval_2tail = count_2tail/float(S)

        orf2data[orf] = (orf, orf2data[orf][0], orf2data[orf][1],orf2TAs.get(orf,0), len(reads), sum1, sum2, sum2-sum1, pval_2tail)
        pval.append(pval_2tail)
    
        count += 1
        wx.CallAfter(pubmsg, "resampling", msg="Running Resampling Method... %2.0f%%" % (100.0*(count+1)/(G)))

    qval = fdr_corrected_pval(pval)
    count = 0
    for orf in sorted(orf2reads):
        output.write("%s\t%s\t%s\t%d\t%d\t%1.1f\t%1.1f\t%1.1f\t%1.5f" % orf2data[orf])
        output.write("\t%1.5f\n" % qval[count])
        count += 1

    output.close()


