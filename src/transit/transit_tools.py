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
import os
import math
import ntpath
import numpy
import scipy.optimize
import scipy.stats


def aton(aa):
    return(((aa-1)*3)+1)

def parseCoords(strand, aa_start, aa_end, start, end):
    if strand == "+":
        return((aton(aa_start) + start,  aton(aa_end) + start))
    # Coordinates are Reversed... to match with Trash FILE TA coordinates
    if strand == "-":
        return((end - aton(aa_end), end - aton(aa_start)))



def get_wig_stats(path):
    reads = []
    for line in open(path):
        if line[0] not in "0123456789": continue
        tmp = line.split()
        pos = int(tmp[0])
        rd = float(tmp[1])
        reads.append(rd)
    reads = numpy.array(reads)
    
    density = numpy.mean(reads>0)
    meanrd = numpy.mean(reads)
    nzmeanrd = numpy.mean(reads[reads>0])
    nzmedianrd = numpy.median(reads[reads>0])
    maxrd = numpy.max(reads)
    totalrd = numpy.sum(reads)

    skew = scipy.stats.skew(reads[reads>0])
    kurtosis = scipy.stats.kurtosis(reads[reads>0])
    return (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis)


def get_reads(path):
    data = []
    for line in open(path):
        if line[0] not in "0123456789": continue
        tmp = line.split()
        rd = int(tmp[1])
        data.append(rd)
    return data


def fetch_name(filepath):
    return os.path.splitext(ntpath.basename(filepath))[0]


def basename(filepath):
    return ntpath.basename(filepath)


def cleanargs(rawargs):
    args = []
    kwargs = {}
    count = 0
    while count < len(rawargs):
        if rawargs[count].startswith("-"):
            try:
                kwargs[rawargs[count][1:]] = rawargs[count+1]
                count += 1
            except IndexError as IE:
                kwargs[rawargs[count][1:]] = True
        else:
            args.append(rawargs[count])
        count += 1

    return (args, kwargs)


def get_pos_hash(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        start = int(tmp[1])
        end = int(tmp[2])
        for pos in range(start, end+1):
            if pos not in hash: hash[pos] = []
            hash[pos].append(orf)
    return hash


def get_data(wig_list):
    """ Returns a tuple of (data, position) containing a matrix of raw read counts, and list of coordinates. """
    K = len(wig_list)
    T = 0
    for line in open(wig_list[0]):
        if line.startswith("#"): continue
        if line.startswith("location"): continue
        if line.startswith("variable"): continue
        T+=1
    
    data = numpy.zeros((K,T))
    position = numpy.zeros(T)
    for j,path in enumerate(wig_list):
        reads = []
        i = 0
        for line in open(path):
            if line.startswith("#"): continue
            if line.startswith("location"): continue
            if line.startswith("variable"): continue
            tmp = line.split()
            pos = int(tmp[0])
            rd = float(tmp[1])
            data[j,i] = rd
            position[i] = pos
            i+=1
    return (data, position)


def combine_replicates(data, method="Sum"):

    if method == "Sum":
        combined = numpy.round(numpy.sum(data,0))
    elif method == "Mean":
        combined = numpy.round(numpy.mean(data,0))
    elif method == "TTRMean":
        factors = transit_tools.TTR_factors(data)
        data = factors * data
        target_factors = transit_tools.norm_to_target(data, 100)
        data = target_factors * data
        combined = numpy.round(numpy.mean(data,0))
    else:
        combined = data[0,:]

    return combined



def normalize_data(data, method="nonorm", wigList=[], annotationPath=""):

    factors = []
    if method == "nzmean":
        factors = nzmean_factors(data)
        data = factors * data
    elif method == "totreads":
        factors = totreads_factors(data)
        data = factors * data
    elif method == "TTR":
        factors = TTR_factors(data)
        data = factors * data
    elif method == "zinfnb":
        factors = zinfnb_factors(data)
        data = factors * data
    elif method == "quantile":
        data = quantile_norm(data)
    elif method == "betageom":
        data = betageom_norm(data)
    elif method == "aBGC":
        data = aBGC_norm(data)
    elif method == "emphist":
        assert ctrlList != None, "Control list cannot be empty!"
        assert expList != None, "Experimental list cannot be empty!"
        assert annotationPath != "", "Annotation path cannot be empty!"
        factors = emphist_factors(wigList, annotationPath)
        data = factors * data
    else:
        method = "nonorm"
        pass

    return (data, factors)



def nzmean_factors(data):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs_hit = numpy.sum(data > 0, 1)
    mean_hits = total_hits/TAs_hit
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    factors = numpy.zeros((K,1))
    factors[:,0] = grand_mean/mean_hits
    return factors



def totreads_factors(data):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs = float(N)
    mean_hits = total_hits/TAs
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    factors = numpy.zeros((K,1))
    factors[:,0] = grand_mean/mean_hits
    return factors


def emphist_factors(wig_list, prot_path):
    orf2info = get_gene_info(prot_path)
    hash = get_pos_hash(prot_path)
    (data, position) = get_data(wig_list)
    orf2reads, orf2pos = get_gene_reads(hash, data, position, orf2info)
    K = len(data)
    N = len(data[0])
    temp = []
    for j in range(K):
        reads_per_gene = []
        for orf in sorted(orf2reads.keys()):
            tempdata = numpy.array(orf2reads[orf])
            if len(tempdata) > 0:
                reads_per_gene.append(numpy.sum(tempdata[:,j]))
        temp.append(reads_per_gene)

    temp = numpy.array(temp)

    factors = numpy.ones((K,1))
    for j in range(1, K):
        ii_good  = numpy.logical_and(temp[0,:] > 0,  temp[j,:] > 0)
        logFC = numpy.log(temp[j,ii_good]/temp[0,ii_good])
        mean = numpy.mean(logFC)
        std = numpy.sqrt(numpy.var(logFC))
        X = numpy.linspace(mean - (5*std),  mean + (std*5), 50000)
        R = scipy.stats.gaussian_kde(logFC)
        Y = R(X)
        peakLogFC = X[Y.argmax()]
        if peakLogFC < 0:
            factors[j,0] = numpy.exp(abs(peakLogFC))
        else:
            factors[j,0] = 1.0/numpy.exp(abs(peakLogFC))

    return factors



def aBGC_norm(data, doTotReads = True, bgsamples = 200000):
    K,N = data.shape
    norm_data = numpy.zeros(data.shape)
    S = bgsamples
    F = [i/100.0 for i in range(0,31) if i % 2 == 0]
    BGC = []
    param_list = []
    for j in range(K):

        nzdata = data[j][data[j] > 0]
        nzdata.sort()
        Nall = len(data[j])
        Nnz = len(nzdata)
        GOF_list = []
        for frac in F:
            tQ = numpy.arange(0,Nnz)/float(Nnz)
            rho = 1.0/(scipy.stats.trim_mean(nzdata, frac))
            rho_to_fit = rho

            try:
                A = (numpy.sum(numpy.power(numpy.log(1.0-tQ),2)))/(numpy.sum(nzdata*numpy.log(1.0-tQ)))
                Kp = (2.0 * numpy.exp(A) - 1)   /(numpy.exp(A) + rho - 1)
                temp = scipy.stats.geom.rvs(scipy.stats.beta.rvs(Kp*rho, Kp*(1-rho), size=S), size=S)
            except Except as e:
                print "aBGC Error:", str(e)
                print "%rho=s\tKp=%s\tA=%s" % (rho, Kp, A)
                temp = scipy.stats.geom.rvs(0.01, size=S)

            corrected_nzdata = [cleaninfgeom(scipy.stats.geom.ppf(ecdf(temp, x), rho_to_fit), rho_to_fit) for x in nzdata]
            corrected_nzmean = numpy.mean(corrected_nzdata)

            Fp = scipy.stats.geom.ppf(numpy.arange(1,Nnz+1)/float(Nnz), 1.0/corrected_nzmean)
            ii_inf = Fp == float("inf")
            Fp[ii_inf] = max(Fp[~ii_inf]) + 100
            ch2_indiv = numpy.power(corrected_nzdata- Fp, 2)/ Fp
            GOF = max(ch2_indiv)
            GOF_list.append((GOF, frac, rho_to_fit, Kp))

        gof, frac, best_rho, best_Kp = sorted(GOF_list)[0]
        BGsample = scipy.stats.geom.rvs(scipy.stats.beta.rvs(best_Kp*best_rho, best_Kp*(1-best_rho), size=S), size=S)
        #BGC.append(dict([(x, removeinf(scipy.stats.geom.ppf(ecdf(temp, x), best_rho), best_rho)) for x in data[j]]))
        for i in range(N):
            norm_data[j,i] = cleaninfgeom(scipy.stats.geom.ppf(ecdf(BGsample, data[j,i]), best_rho), best_rho)

    if doTotReads:
        return totreads_factors(norm_data) * norm_data
    return norm_data



def empirical_theta(X):
    return numpy.mean(X > 0)

def trimmed_empirical_mu(X, t=0.05):
    return scipy.stats.trim_mean(X[X > 0], t)

def TTR_factors(data, thetaEst=empirical_theta, muEst=trimmed_empirical_mu):
    K = len(data)
    N = len(data[0])

    factors = numpy.zeros((K,1))
    for j in range(K):
        factors[j] = (thetaEst(data[0]) * muEst(data[0]))/(thetaEst(data[j]) * muEst(data[j]))

    return factors



def Fzinfnb(params, args):
    pi, mu, r = params
    Fdata = args
    temp0 = numpy.nan_to_num(numpy.log(pi + scipy.stats.nbinom.pmf(Fdata[Fdata==0], mu, r)))
    tempnz = numpy.nan_to_num(numpy.log(1.0-pi)+scipy.stats.nbinom.logpmf(Fdata[Fdata>0], mu, r))
    negLL = -(numpy.sum(temp0) + numpy.sum(tempnz))
    return negLL


def zinfnb_factors(data):
   
    N = len(data)
    G = len(data[0])

    factors = numpy.zeros((N, 1))
    for j in range(N):
        initParams = [0.3, 10, 0.5]
        M = "L-BFGS-B"
        Fdata = numpy.array(data[j])
        results = scipy.optimize.minimize(Fzinfnb, initParams, args=(Fdata,), method=M, bounds=[(0.0001, 0.9999),(0.0001, None),(0.0001, 0.9999)])

        
        pi, n, p = results.x
        mu = n*(1-p)/p
        factors[j,0] = 1.0/mu
    

    return factors



def thetanorm_factors(data):

    K = len(data)
    N = len(data[0])
    
    factors = numpy.zeros((K,1))
    for j in range(K):
        factors[j] = (numpy.mean(data[0] > 0) * scipy.stats.trim_mean(data[0][data[0] > 0], 0.05))/(numpy.mean(data[j] > 0) * scipy.stats.trim_mean(data[j][data[j] > 0], 0.05))

    return factors





def quantile_norm(data):
    """Performs Quantile Normalization as described by Bolstad et al. 2003"""
    
    N = len(data)
    G = len(data[0])

    #Sort columns
    s_data = numpy.array([sorted(col) for col in data])

    #Get ranks of original data
    ranks = numpy.zeros(data.shape, dtype=int)
    for j in range(N):
        ranks[j,:] = scipy.stats.rankdata(data[j], method='dense')

    #Get empirical distribution
    ranked_means = numpy.mean(s_data,0)
    
    #Create dictionary of rank to new empirical values
    rank2count = dict([(r,c) for (r,c) in zip(scipy.stats.rankdata(ranked_means, method='dense'), ranked_means)])

    #Assign values
    norm_data = numpy.zeros(data.shape)
    for i in range(G):
        norm_data[:,i] = [rank2count[ranks[j,i]] for j in range(N)]
    return norm_data




def ecdf(S, x):
    return numpy.sum(S<=x)/float(len(S))


def cleaninfgeom(x, rho):
    if x == float('inf'):
        return scipy.stats.geom.ppf(0.9999999999999999, rho)
    else:
        return x



def betageom_norm(data, doNZMean = True, bgsamples=200000):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs_hit = numpy.sum(data > 0,1)
    mean_hits = total_hits/TAs_hit
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    norm_data = numpy.zeros(data.shape)
    for j in range(K):

        tQ = numpy.arange(0,N)/float(N)
        eX = numpy.array([rd for rd in data[j]])
        eX.sort()

        rho = max(1.0/scipy.stats.trim_mean(eX+1, 0.001), 0.0001)
        A = (numpy.sum(numpy.power(numpy.log(1.0-tQ),2)))/(numpy.sum(eX*numpy.log(1.0-tQ)))
        Kp = max((2.0 * numpy.exp(A) - 1)   /(numpy.exp(A) + rho - 1), 10)
    
        try:
            BGsample = scipy.stats.geom.rvs(scipy.stats.beta.rvs(Kp*rho, Kp*(1-rho), size=bgsamples), size=bgsamples)
        except Exception as e:
            print "BGC ERROR with rho=%f, Kp=%f, A=%s" % (rho, Kp, A)
            print str(e)
            BGsample = scipy.stats.geom.rvs(rho, size=bgsamples)

        for i in range(N):
            norm_data[j,i] = cleaninfgeom(scipy.stats.geom.ppf(ecdf(BGsample, data[j,i]), 1.0/grand_mean), 1.0/grand_mean)

        #mapping = dict([(x, cleaninfgeom(scipy.stats.geom.ppf(ecdf(BGsample, x), 1.0/grand_mean), 1.0/grand_mean)) for x in data[j]])
        #for i in range(N):
        #    try:
        #        norm_data[j,i] = mapping[data[j,i]]
        #    except KeyError:
        #        print "Error: %s  | key = %s not found. Using original Data. Notify authors!" % (KeyError, data[j,i])
        #        norm_data[j,i] = data[j,i]
    
    if doNZMean:
        return nzmean_norm(norm_data)
    return norm_data







def nzmean_norm(data):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs_hit = numpy.sum(data > 0,1)
    mean_hits = total_hits/TAs_hit
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    factors = numpy.zeros((K,1))
    factors[:,0] = grand_mean/mean_hits
    return factors * data
    

def norm_to_target(data, target):
    (K,N) = data.shape
    factors = numpy.zeros((K,1))
    factors[:,0] = float(target)/numpy.mean(data,1)
    return factors


def get_gene_reads(hash, data, position, orf2info, ignoreCodon=True, ignoreNTerm=0, ignoreCTerm=0, orf_list=set()):
    (K,N) = data.shape

    orf2reads = dict([(orf,[]) for orf in orf_list])
    orf2pos = dict([(orf,[]) for orf in orf_list])

    for i in range(N):
        coord = position[i]
        genes_with_coord = hash.get(coord, [])
        for gene in genes_with_coord:
            if gene not in orf2reads: orf2reads[gene] = []
            if gene not in orf2pos: orf2pos[gene] = []

            start,end,strand = orf2info.get(gene, [0,0,0,0])[2:5]
            #print gene, coord, start, end, strand, data[:,i].shape

            if strand == "+":
                #Ignore TAs at stop codon
                if ignoreCodon and coord > end-3:
                    continue

            else:
                #ignore TAs at stop codon
                if ignoreCodon and coord < start + 3:
                    continue

            #print "passed first IF"

            #Ignore TAs at beginning n%
            if (coord-start)/float(end-start) < (ignoreNTerm/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and NTerm", ignoreNTerm/100.0
                continue
            
            #print "passed second IF"

            #Ignore TAs at end c%
            if (coord-start)/float(end-start) > ((100-ignoreCTerm)/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and Cterm", ignoreCTerm/100.0
                continue



            #print "adding data"
            #print gene, data[:,i]
            orf2reads[gene].append(data[:,i])
            orf2pos[gene].append(position[i])
    return (orf2reads, orf2pos)
   


def tricube(X):
    result = numpy.zeros(len(X))
    ii = numpy.logical_and(X >= -1, X <= 1)
    result[ii] = numpy.power(1 - numpy.power(numpy.abs(X[ii]), 3), 3)
    return result


def loess(X, Y, h=10000):
    smoothed = numpy.zeros(len(Y))
    for i,x in enumerate(X):
        W = tricube((X-x)/float(h))
        sW = numpy.sum(W)
        wsX = numpy.sum(W*X)
        wsY = numpy.sum(W*Y)
        wsXY = numpy.sum(W*X*Y)
        sXX = numpy.sum(X*X)
        B = (sW * wsXY - wsX * wsY)/(sW * sXX - numpy.power(wsX,2))
        A = (wsY - B*wsX) / sW
        smoothed[i] = B*x + A
    return smoothed




def loess_correction(X, Y, h=10000, window=100):
    Y = numpy.array(Y)
    size = len(X)/window + 1
    x_w = numpy.zeros(size)
    y_w = numpy.zeros(size)
    for i in range(len(X)/window + 1):
        x_w[i] = window*i
        y_w[i] = sum(Y[window*i:window*(i+1)])

    ysmooth = loess(x_w, y_w, h)
    mline = numpy.mean(y_w)
    y_w * (ysmooth/mline)

    normalized_Y = numpy.zeros(len(Y))
    for i in range(size):
        normalized_Y[window*i:window*(i+1)] = Y[window*i:window*(i+1)] * (ysmooth[i]/mline)
    
    return normalized_Y


def fdr_post_prob(Z_raw, ALPHA=0.05):
    Z = numpy.sort(Z_raw)[::-1]
    W = 1 - Z
    N = len(Z)

    ess_threshold = 1.00
    INDEX = range(3, N+1)
    count = 0
    for i in INDEX:
        count +=1
        wi = 1 - Z[i-1]
        ai_n = (ALPHA*i)/N
        mean_wi = numpy.average(W[0:i-2])
        delta_w = wi - mean_wi
        if delta_w > ai_n:
            ess_threshold = Z[i-1]
            break

    noness_threshold = 0.00
    count = 0
    INDEX = range(0, N+1)
    INDEX.sort(reverse=True)
    for i in INDEX:
        wi = Z[N-i+1]
        ai_n = (ALPHA*i)/N
        mean_wi = numpy.average(Z[N-i+1:])
        delta_w = Z[N-i+1] - mean_wi
        count +=1
        if ai_n > delta_w:
            break
        noness_threshold = Z[N-i]

    return(ess_threshold, noness_threshold)



def get_gene_info(path):
    orf2info = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        name = tmp[7]
        desc = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        strand = tmp[3]
        orf2info[orf] = (name, desc, start, end, strand)
    return orf2info

def get_gene_name(path):
    orf2info = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        name = tmp[7]
        desc = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        strand = tmp[3]
        orf2info[orf] = name
    return orf2info
    






