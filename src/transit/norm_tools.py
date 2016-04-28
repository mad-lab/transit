import sys
import numpy
import scipy.stats
import scipy.optimize
import warnings




def normalize_data(data, method="nonorm", wigList=[], annotationPath=""):
    """Normalizes the numpy array by the given normalization method.

    Arguments:
        data (numpy array): (K,N) numpy array defining read-counts at N sites
            for K datasets.
        method (str): Name of the desired normalization method.
        wigList (list): List of paths for the desired wig-formatted datasets.
        annotationPath (str): Path to the prot_table annotation file.
    
    Returns:
        numpy array: Array with the normalized data.
        list: List containing the normalization factors. Empty if not used.
    """
    factors = []
    if method == "nonorm":
        pass
    elif method == "nzmean":
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
        warnstr = "Normalization method '%s' is unknown. Read-counts were not normalized." % (method)
        warnings.warn(warnstr)
    return (data, factors)

def nzmean_factors(data):
    """Returns the normalization factors for the data, using the NZMean method."""
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
    """Returns the normalization factors for the data, using the total reads
    method."""
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
    """Returns the normalized data, using the empirical hist method."""
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
    """Returns the normalized data using the aBGC method."""
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
    """Calculates the observed density of the data."""
    return numpy.mean(X > 0)

def trimmed_empirical_mu(X, t=0.05):
    """Estimates the trimmed mean of the data."""
    return scipy.stats.trim_mean(X[X > 0], t)

def TTR_factors(data, thetaEst=empirical_theta, muEst=trimmed_empirical_mu):
    """Returns the normalization factors for the data, using the TTR method."""
    K = len(data)
    N = len(data[0])

    factors = numpy.zeros((K,1))
    for j in range(K):
        factors[j] = (thetaEst(data[0]) * muEst(data[0]))/(thetaEst(data[j]) * muEst(data[j]))

    return factors

def Fzinfnb(params, args):
    """Objective function for the zero-inflated NB method."""
    pi, mu, r = params
    Fdata = args
    temp0 = numpy.nan_to_num(numpy.log(pi + scipy.stats.nbinom.pmf(Fdata[Fdata==0], mu, r)))
    tempnz = numpy.nan_to_num(numpy.log(1.0-pi)+scipy.stats.nbinom.logpmf(Fdata[Fdata>0], mu, r))
    negLL = -(numpy.sum(temp0) + numpy.sum(tempnz))
    return negLL


def zinfnb_factors(data):
    """Returns the normalization factors for the data using the zero-inflated
    negatibe binomial method."""
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
    """Calculates an empirical CDF of the given data."""
    return numpy.sum(S<=x)/float(len(S))


def cleaninfgeom(x, rho):
    """Returns a 'clean' output from the geometric distribution."""
    if x == float('inf'):
        return scipy.stats.geom.ppf(0.9999999999999999, rho)
    else:
        return x



def betageom_norm(data, doNZMean = True, bgsamples=200000):
    """Returns normalized data according to the BGC method."""
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

    if doNZMean:
        return nzmean_norm(norm_data)
    return norm_data



def norm_to_target(data, target):
    """Returns factors to normalize the data to the given target value."""
    (K,N) = data.shape
    factors = numpy.zeros((K,1))
    factors[:,0] = float(target)/numpy.mean(data,1)
    return factors





