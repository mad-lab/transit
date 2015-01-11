import sys
import random
import numpy
import transit_tools

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



def runResampling(ctrlString, expString, annotationPath, sampleSize, output, wx, pubmsg, doNormalize=True):

    ctrlList = ctrlString.split(",")
    expList = expString.split(",")

    N1 = len(ctrlList)
    N2 = len(expList)

    orf2info = transit_tools.get_gene_info(annotationPath)

    hash = transit_tools.get_pos_hash(annotationPath)
    (data, position) = transit_tools.get_data(ctrlList + expList)
    factors = transit_tools.get_norm_factors(data)
    norm_data = factors * data
    orf2reads,orf2pos = transit_tools.get_gene_reads(hash, norm_data, position, orf2info, orf_list=orf2info.keys())
    
    S = sampleSize

    output.write("#Resampling\n")
    output.write("#Command: python %s\n" % " ".join(sys.argv))
    output.write("#Control Samples:  %s\n" % ", ".join(ctrlList))
    output.write("#Experimental Samples:  %s\n" % ", ".join(expList))
    if doNormalize:
        output.write("#Normalization factors: %s\n" % "\t".join(["%1.4f" % f for f in factors]))
    else:
        output.write("#Not normalized\n")

    output.write("#Orf\t%Name\tDescription\tN\tTAs Hit\tAvg Rd 1\tAvg Rd 2\tDelta Rd\tp-value\tp-adj\n")
    count = 0
    G = len(orf2reads)
    orf2out= {}
    pval=[]
    for orf in sorted(orf2reads):
        if orf == "non-coding": continue

        fullreads = orf2reads[orf]
        reads = numpy.array([row for row in fullreads if numpy.sum(row) > 0 ])
        if len(reads) >0:
            sum1,sum2 = numpy.sum(reads[:,:N1]), numpy.sum(reads[:,N1:])
        else:
            sum1 = 0; sum2 = 0

        count_utail = 0
        count_ltail = 0
        count_2tail = 0
   
        for s in range(S):
            #Reads
            if len(reads) == 0: break
            perm = numpy.random.permutation(reads.flatten()).reshape(reads.shape)
            sumA,sumB = numpy.sum(perm[:,:N1]), numpy.sum(perm[:,N1:])

            #all_perm=np.array((list(itertools.permutations([0,1,2,3]))))
            #(reads.flatten()[(all_perm[numpy.random.randint(0,len(all_perm),size=100)]+3*np.arange(100)[...,numpy.newaxis]).flatten()]).reshape(reads.shape)

            if sumB-sumA >= sum2-sum1: count_utail+=1
            if sumB-sumA <= sum2-sum1: count_ltail+=1
            if abs(sumB-sumA) >= abs(sum2-sum1): count_2tail+=1

        #Reads - pvalues
        pval_utail = count_utail/float(S)
        pval_ltail = count_ltail/float(S)
        pval_2tail = count_2tail/float(S)

        orf2out[orf] = (orf, orf2info[orf][0], orf2info[orf][1], len(fullreads), len(reads), sum1, sum2, sum2-sum1, pval_2tail)
        pval.append(pval_2tail)
    
        count += 1
        wx.CallAfter(pubmsg, "resampling", msg="Running Resampling Method... %2.0f%%" % (100.0*(count+1)/(G)))

    qval = fdr_corrected_pval(pval)
    count = 0
    for orf in sorted(orf2out):
        output.write("%s\t%s\t%s\t%d\t%d\t%1.1f\t%1.1f\t%1.1f\t%1.5f" % orf2out[orf])
        output.write("\t%1.5f\n" % qval[count])
        count += 1

    output.close()


