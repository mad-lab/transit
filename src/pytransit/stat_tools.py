import sys,math,random
import numpy
import scipy.stats


def sample_trunc_norm_post(data, S, mu0, s20, k0, nu0):
    n = len(data)
    s2 = numpy.var(data,ddof=1)
    ybar = numpy.mean(data)
    kn = k0+n
    nun = nu0+n
    mun = (k0*mu0 + n*ybar)/float(kn)
    s2n = (1.0/nun) * (nu0*s20 + (n-1)*s2 + (k0*n/float(kn))*numpy.power(ybar-mu0,2))

    s2_post = 1.0/scipy.stats.gamma.rvs(nun/2.0, scale=2.0/(s2n*nun), size=S)

    # Truncated Normal since counts can't be negative
    min_mu = 0
    max_mu = 1000000
    trunc_a = (min_mu-mun)/numpy.sqrt(s2_post/float(kn))
    trunc_b = (max_mu-mun)/numpy.sqrt(s2_post/float(kn))

    mu_post = scipy.stats.truncnorm.rvs(a=trunc_a, b=trunc_b, loc=mun, scale=numpy.sqrt(s2_post/float(kn)), size=S)

    return (mu_post, s2_post)

#

def FWER_Bayes(X):
    ii = numpy.argsort(numpy.argsort(X))
    P_NULL = numpy.sort(X)
    W = 1 - P_NULL
    N = len(P_NULL)
    P_ALT = numpy.zeros(N)
    for i in range(N):
        P_ALT[i] = 1.0 - numpy.prod(W[:i+1])
    return P_ALT[ii]

#

def bFDR(X):
    N = len(X)
    ii = numpy.argsort(numpy.argsort(X))
    P_NULL = numpy.sort(X)
    P_ALT = numpy.zeros(N)
    for i in range(N):
        P_ALT[i] = numpy.mean(P_NULL[:i+1])
    return P_ALT[ii]

#

def HDI_from_MCMC(posterior_samples, credible_mass=0.95):
    # Credit to 'user72564'
    # https://stackoverflow.com/questions/22284502/highest-posterior-density-region-and-central-credible-region
    # Computes highest density interval from a sample of representative values,
    # estimated as the shortest credible interval
    # Takes Arguments posterior_samples (samples from posterior) and credible mass (normally .95)
    sorted_points = sorted(posterior_samples)
    ciIdxInc = numpy.ceil(credible_mass * len(sorted_points)).astype('int')
    nCIs = len(sorted_points) - ciIdxInc
    ciWidth = [0]*nCIs
    for i in range(0, nCIs):
        ciWidth[i] = sorted_points[i + ciIdxInc] - sorted_points[i]
    HDImin = sorted_points[ciWidth.index(min(ciWidth))]
    HDImax = sorted_points[ciWidth.index(min(ciWidth))+ciIdxInc]
    return(HDImin, HDImax)

#

def transformToRange(X, new_min, new_max, old_min=None, old_max=None):

    if old_min == None:
        old_min = min(X)
    if old_max == None:
        old_max = max(X)
    
    old_range = old_max - old_min

    new_range = new_max - new_min
    return [float(x - old_min) / old_range * new_range + new_min for x in X]

#

def fact(n):
    if n == 0: return (1)
    else: return reduce(lambda x,y: x*y, range(1,n+1))

#

def comb1(n,k):
    prod = 1
    for i in range(1,k+1):
        prod = prod * (n - (k-i))/float(i)
    return(prod)

#

def comb(n, k):
    if k < 0 or k > n:
        return 0
    if k > n - k: # take advantage of symmetry
        k = n - k
    c = 1
    for i in range(k):
        c = c * (n - (k - (i+1)))
        c = c // (i+1)
    return c

#

def norm(x, mu,sigma):
    """Normal distribution"""
    sigma = float(sigma)
    return(1/(sigma*(math.sqrt(2*math.pi))) * math.exp( -0.5 * math.pow( (x-mu)/sigma,2)))

#

def binom(k,n,p):
    """Binomial distribution. Uses Normal approximation for large 'n' """
    if n >= 100:
        return(norm(k, n*p, math.sqrt(n*p*(1-p)) ) )
    else:
        return(comb(n,k) * math.pow(p, k) * math.pow(1-p, n-k))

#

def binom_cdf(k,n,p):
    """CDF of the binomial distribution"""
    return(sum([binom(i,n,p) for i in range(0,k+1)]))

#

def binom_test(k,n,p, type="two-sided"):
    """Does a binomial test given success, trials and probability."""
    if type == "less": return(binom_cdf(k,n,p))
    elif type == "greater": return(1-binom_cdf(k-1,n,p))
    else:
        if p == 0: return(1) #return(k == 0)
        elif p == 1: return(1) #return(k == n)
        else:
            relErr = 1 + 1e-7
            d = binom(k,n,p)
            m = n * p
            if k == m: return(1)
            elif (k < m):
                ri = range(int(math.ceil(m)), n+1)
                y = sum([1 for j in ri if binom(j,n,p) <= d*relErr])
                return(binom_cdf(k,n,p) + (1-binom_cdf(int(n-y),n,p)))
            else:
                ri = range(0, int(math.floor(m)))
                y = sum([1 for j in ri if binom(j,n,p) <= d*relErr])
                return(binom_cdf(y-1,n,p) + (1-binom_cdf(k-1,n,p)))


##############################
# Bernoulli Diff Distribution
def dberndiff(d, peq, p01, p10):
    N = numpy.size(d)
    if N == 0:
        return 0.0
    if N == 1:
        if type(d) == type(()):
            d = d[0]
        if d == 0:
            return peq
        else:
            if d == -1:
                return p01
            if d == 1:
                return p10
            return 0.0
#
    else:
        d = numpy.array(d)
        result = numpy.zeros(N)
        result[d == -1] = p01
        result[d == 0] = peq
        result[d == 1] = p10
        return result

#

def qberndiff(d, peq, p01, p10):
    return numpy.sum([ dberndiff(x, peq, p01, p10) for x in range(-1, d + 1) ])

#############################
# Binomial Diff Distribution

def dbinomdiff(d, n, P):
    S = numpy.array(my_perm(d, n))
    return numpy.sum(multinomial(S, P))

#

def qbinomdiff(d, n, peq, p01, p10):
    return numpy.sum([ dbinomdiff(x, n, peq, p01, p10) for x in range(-n, d + 1) ])

#

def my_perm(d, n):
    S = []
    if d == 0:
        for i in range(n + 1):
            r = n - i
            if isEven(r):
                S.append((int(r / 2.0), i, int(r / 2.0)))
#
    if d > 0:
        for i in range(d, n + 1):
            r = n - i
            if i == d:
                S.append((0, n - d, d))
            elif i > d:
                r = n - (i + (i - d))
                if 0 <= r <= n:
                    S.append((i - d, r, i))
#
    if d < 0:
        for i in range(abs(d), n + 1):
            r = n - i
            if i == abs(d):
                S.append((-d, n + d, 0))
            elif i > d:
                r = n - (i + (i + d))
                if 0 <= r <= n:
                    S.append((i, r, i + d))
#
    return S

#

def multinomial(K, P):
    N = numpy.sum(K, 1)
    if K.shape == P.shape:
        return tricoeff(N, K) * numpy.prod([ numpy.power(P[i], K[i]) for i in range(len(K)) ], 1)
    else:
        return tricoeff(N, K) * numpy.prod([ numpy.power(P, K[i]) for i in range(len(K)) ], 1)

#

def log_fac(n):
    return numpy.sum(numpy.log(numpy.arange(2, n + 1)))

#

def tricoeff(N, S):
    try:
        LOG_FAC
    except NameError:
        LOG_FAC = []
        for i in range(numpy.max(N) + 1):
            LOG_FAC.append(log_fac(i))
#
        LOG_FAC = numpy.array(LOG_FAC)
#
    return numpy.exp(LOG_FAC[N] - (LOG_FAC[S[:, 0]] + LOG_FAC[S[:, 1]] + LOG_FAC[S[:, 2]]))

#

def isEven(x):
    return x % 2 == 0

#

def regress(X,Y):
    """Performs linear regression given two vectors, X, Y."""
    N = len(X)
    xbar = numpy.average(X)
    ybar = numpy.average(Y)
    xybar = numpy.average([X[i]*Y[i] for i in range(N)])
    x2bar = numpy.average([X[i]*X[i] for i in range(N)])
    B = (xybar - xbar*ybar)/(x2bar - xbar*xbar)
    A0 = ybar - B*xbar

    yfit = [ A0 + B *X[i] for i in range(N)]
    yres = [Y[i] - (A0 + B *X[i]) for i in range(N)]
    var = sum([math.pow(yres[i],2) for i in range(N) ])/(N-2)
    std = math.sqrt(var)

    return(B, A0, std)

#

def boxcoxtransform(x, lambdax):
    """
    Performs a box-cox transformation to data vector X.
    WARNING: elements of X should be all positive! 
    Fixed: '>' has changed to '<'
    """
    if x <= 0:
        raise ValueError("Nonpositive value(s) in X vector")
    if abs(lambdax) < 1.0e-5:
        return(math.log(x))
    else:
        return((x**lambdax - 1.0)/lambdax)
    #return math.log(x) if abs(lambdax) < 1.0e-5 else (x**lambdax - 1.0)/lambdax

#

def loglik(X, lambdax):
    """
    Computes the log-likelihood function for a transformed vector Xtransform.     
    """
    n = len(X)
    Xtrans = [boxcoxtransform(x, lambdax) for x in X]
    meanX = sum(Xtrans) / float(n)
    S2 = (lambdax - 1.0) * sum([math.log(x) for x in X])
    S = sum([(x-meanX) **2 for x in Xtrans])
    S1= (-n/2.0)*math.log(S/n) 
    return  S2+S1

#

def boxcoxTable(X, minlambda, maxlambda, dellambda):
    """
    Returns a table of (loglik function, lambda) pairs
    for the data.
    """
    # Create a table (lambda, loglik)
    out = []
    vallambda = minlambda
    while vallambda <= maxlambda+1.0e-5:
        llik = loglik(X, vallambda)
        out.append((llik, vallambda))
        vallambda += dellambda
    return out  

#

def phi_coefficient(X,Y):
    """Calculates the phi-coefficient for two bool arrays"""
    N = len(X)
    assert len(X) == len(Y), "Length of arrays must be equal"
    x1y1 = sum([int(X[j]) == int(Y[j])  ==  1 for j in range(N)])
    x1y0 = sum([int(X[j]) == 1 and int(Y[j])  ==  0 for j in range(N)])
    x0y1 = sum([int(X[j]) == 0 and int(Y[j])  ==  1 for j in range(N)])
    x0y0 = sum([int(X[j]) == int(Y[j])  ==  0 for j in range(N)])
    x1 = x1y1 + x1y0
    x0 = x0y1 + x0y0
    y1 = x1y1 + x0y1
    y0 = x1y0 + x0y0
    phi_coeff = (x1y1*x0y0 - x1y0*x0y1)/math.sqrt(x1*x0*y1*y0)
    return phi_coeff

#

def BH_fdr_correction(X):
    """Adjusts p-values using the Benjamini Hochberg procedure"""
    n = len(X)
    qvalues = numpy.zeros(n)
    pvalues = numpy.array(X)
    pvalues.sort()
    pvalues = pvalues[::-1]
    
    for i in range(n):
        rank = n - i
        qvalues[i] = n/float(rank) * pvalues[i]
 
    for i in range(0, n-1):
        if qvalues[i] < qvalues[i+1]:
            qvalues[i+1] = qvalues[i]
        
    p2qval = dict([(p,q) for (p,q) in zip(pvalues,qvalues)])
    return numpy.array([p2qval[p] for p in X])

#

def bayesian_ess_thresholds(Z_raw, ALPHA=0.05):
    """Returns Essentiality Thresholds using a BH-like procedure"""
    Z = numpy.sort(Z_raw)[::-1]
    W = 1 - Z
    N = len(Z)

    ess_threshold = 1.00
    INDEX = list(range(3, N+1))
    count = 0
    for i in INDEX:
        count +=1
        wi = 1 - Z[i-1]
        ai_n = (ALPHA*i)/N
        mean_wi = numpy.average(W[0:i-2])
        delta_w = wi - mean_wi
        #if count < 30: print(i, wi, ai_n, delta_w)
        if delta_w > ai_n:
            ess_threshold = Z[i-1]
            #print("i", i)
            break

    noness_threshold = 0.00
    count = 0
    INDEX = list(range(0, N+1))
    INDEX.sort(reverse=True)
    for i in INDEX:
        wi = Z[N-i+1]
        ai_n = (ALPHA*i)/N
        mean_wi = numpy.average(Z[N-i+1:])
        delta_w = Z[N-i+1] - mean_wi
        count +=1
        #print(count)
        #if count < 20:
        #   print(i, wi, ai_n, mean_wi, delta_w, N-i+1, N-1, W[N-i-1], W[i-1])

        if ai_n > delta_w:
        #   print(i, wi, ai_n, mean_wi, delta_w, N-i+1, N-1, W[N-i-1], W[i-1])
            break
        noness_threshold = Z[N-i]

    return(ess_threshold, noness_threshold)

#

def tricube(X):
    #TODO: Write docstring
    result = numpy.zeros(len(X))
    ii = numpy.logical_and(X >= -1, X <= 1)
    result[ii] = numpy.power(1 - numpy.power(numpy.abs(X[ii]), 3), 3)
    return result

#

def loess(X, Y, h=10000):
    #TODO: Write docstring
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

# X is coords, Y is counts

def loess_correction(X, Y, h=10000, window=100):
    #TODO: Write docstring
    Y = numpy.array(Y)
    size = int(len(X)/window) + 1
    x_w = numpy.zeros(size)
    y_w = numpy.zeros(size)
    for i in range(size):
        x_w[i] = window*i
        y_w[i] = sum(Y[window*i:window*(i+1)])

    ysmooth = loess(x_w, y_w, h)
    mline = numpy.mean(y_w)

    normalized_Y = numpy.zeros(len(Y))
    for i in range(size):
        normalized_Y[window*i:window*(i+1)] = Y[window*i:window*(i+1)] / (ysmooth[i]/mline)

    return normalized_Y

#

def F_mean_diff_flat(*args, **kwargs):
    A = args[0]
    B = args[1]
    return numpy.mean(B) - numpy.mean(A)

# 
   
def F_sum_diff_flat(*args, **kwargs):
    A = args[0]
    B = args[1]
    return numpy.sum(B) - numpy.sum(A)

#

def F_mean_diff_dict(*args, **kwargs):
    D = args[0] 
    data1_total = 0; data2_total = 0;
    data1_size = 0; data2_size = 0;
    for L in D:
        data1_total+= numpy.sum(D[L][0])
        data1_size+= len(D[L][0])
        data2_total+= numpy.sum(D[L][1])
        data2_size+= len(D[L][1])
    return (data2_total/float(data2_size)) - (data1_total/float(data1_size))

# 

def F_sum_diff_dict(*args, **kwargs):
    D = args[0]
    data1_total = 0; data2_total = 0;
    for L in D:
        data1_total+= numpy.sum(D[L][0])
        data2_total+= numpy.sum(D[L][1])
    return data2_total - data1_total

#

# input: an a array with counts pooled across both conditions

def F_shuffle_flat(*args, **kwargs):
    X = args[0]
    return numpy.random.permutation(X)

# data is 2D array of counts (samples X TA sites) including all reps in all libs in both conditions (as rows)
# DESTRUCTIVE; caller must make copy of input arg (data) if they want to preserve it

def site_restricted_permutation(data):
  if len(data)==0: return data
  nsamples,nTAs = data.shape[0],data.shape[1]
  for i in range(nTAs): numpy.random.shuffle(data[:,i])
  return data

# apparently, this takes 1 arg, a dictionary D which maps library strings to pairs of count vectors

def F_shuffle_dict_libraries(*args, **kwargs):
    D = args[0]
    E = {}
    for L in D:
        #print("L", L)
        n1 = len(D[L][0])
        combined = numpy.append(D[L][0], D[L][1])
        #print("combined", combined)
        perm = numpy.random.permutation(combined)
        #print("perm", perm)
        #print("perm[:n1]", perm[:n1])
        E[L] = numpy.array([perm[:n1], perm[n1:]],dtype=object)
        #print("D[L]", D[L])
    return E

#

def resampling(data1, data2, S=10000, testFunc=F_mean_diff_flat,
            permFunc=F_shuffle_flat, adaptive=False, lib_str1="", lib_str2="",PC=1,site_restricted=False):
    """Does a permutation test on two sets of data.

    Performs the resampling / permutation test given two sets of data using a
    function defining the test statistic and a function defining how to permute
    the data.

    Args:
        data1: List or numpy array with the first set of observations.
        data2: List or numpy array with the second set of observations.
        S: Number of permutation tests (or samples) to obtain.
        testFunc: Function defining the desired test statistic. Should accept
                two lists as arguments. Default is difference in means between
                the observations.
        permFunc: Function defining the way to permute the data. Should accept
                one argument, the combined set of data. Default is random
                shuffle.
        adaptive: Cuts-off resampling early depending on significance.

    Data arrays: (data1 and data2)
      Regular resampling used to take 1D arrays of counts pooled (flattened) over replicates.
        Now 2D arrays are passed in and flatten them.
        Uses F_shuffle_flat() and F_sum_diff_flat().
      If using library strings, then inputs are 2D arrays of counts for each sample. 
        Character in lib_str indicates which lib it is in.  Make a dict out of these to pass to permFunc.
        Uses F_shuffle_dict_libraries() and F_sum_diff_dict_libraries().
      If site_restricted, keep input arrays as 2D and pass to site_restricted_permutation() and F_sum_diff_flat().

    Returns:
        Tuple with described values
            - test_obs -- Test statistic of observation.
            - mean1 -- Arithmetic mean of first set of data.
            - mean2 -- Arithmetic mean of second set of data.
            - log2FC -- Normalized log2FC the means.
            - pval_ltail -- Lower tail p-value.
            - pval_utail -- Upper tail p-value.
            - pval_2tail -- Two-tailed p-value.
            - test_sample -- List of samples of the test statistic.
    
    :Example:
        >>> import pytransit.stat_tools as stat_tools
        >>> import numpy
        >>> X = numpy.random.random(100)
        >>> Y = numpy.random.random(100)
        >>> (test_obs, mean1, mean2, log2fc, pval_ltail, pval_utail, pval_2tail, test_sample) = stat_tools.resampling(X,Y)
        >>> pval_2tail
        0.2167
        >>> test_sample[:3]
        [0.076213992904990535, -0.0052513291091412784, -0.0038425140184765172]
    
    """

    # Do basic sanity checks:
    # - Check library strings match in some way
    lib_diff = set(lib_str1) ^ set(lib_str2)
    if lib_diff:
        #raise ValueError("At least one library string has a letter not used by the other:\ %s" % ", ".join(lib_diff))
        raise ValueError("At least one library string has a letter not used by the other: " + ", ".join(lib_diff))

    if lib_str1 and site_restricted:
      raise Exception("Cannot do site_restricted resampling with library strings at same time")

    # - Check input has some data
    assert len(data1) > 0, "Data1 cannot be empty"
    assert len(data2) > 0, "Data2 cannot be empty"

    if isinstance(data1,list): data1 = numpy.array(data1)
    if isinstance(data2,list): data2 = numpy.array(data2)

    #TRI note - now I am switching resampling() so caller passes in NON-flattened arrays of counts
    if not site_restricted and not lib_str1: 
      data1 = data1.flatten()
      data2 = data2.flatten()

    count_ltail = 0
    count_utail = 0
    count_2tail = 0

    test_list = []

    # Calculate basic statistics for the input data:
    # flattened (pooled) if not lib_str, else multiple samples (rows) for each lib
    n1 = len(data1) # number of samples (i.e. rows) for site_restricted or lib_str, or pooled counts if flattened
    n2 = len(data2)

    mean1 = 0
    if n1 > 0:
        mean1 = numpy.mean(data1) # over all counts pooled across reps and libs for cond1
    mean2 = 0
    if n2 > 0:
        mean2 = numpy.mean(data2)

    if PC>0: log2FC = math.log((mean2+PC)/(mean1+PC),2) # as of 3/5/20
    else:
      # Only adjust log2FC if one of the means is zero
      if mean1 > 0 and mean2 > 0: log2FC = math.log((mean2)/(mean1),2)
      else: log2FC = math.log((mean2+1.0)/(mean1+1.0),2)

    # Get stats and info based on whether working with libraries or not:
    nTAs = 0
    if lib_str1:
        # note: returns a generator, not a list
        # Get number of TA sites implied
        nTAs = len(data1.flatten())//len(lib_str1)
        assert len(data2.flatten())//len(lib_str2) == nTAs, "Datasets do not have matching sites; check input data and library strings."

        # Get data
        # for lib_str, perm is a dict mapping letters to pairs of numpy arrays (1 for each cond)
        perm = get_lib_data_dict(data1.flatten(), lib_str1, data2.flatten(), lib_str2, nTAs) 
        test_obs = testFunc(perm) 
    else:
        try:
            # site_retricted use F_sum_diff_flat() as testFunc too
            test_obs = testFunc(data1, data2) # first call, actual value from observed counts
        except Exception as e:
            print("")
            print("!"*100)
            print("Error: Could not apply test function to input data!")
            print("data1", data1)
            print("data2", data2)
            print("")
            print("\t%s" % e)
            print("!"*100)
            print("")
            return None

        if site_restricted:
          data = numpy.concatenate((data1,data2),axis=0) # keep it as a 2D array
          perm = data.copy() # this array will get modified with each permutation
        else: # pool all counts (across conditions) into 1 big array
          perm = numpy.zeros(n1+n2)
          perm[:n1] = data1
          perm[n1:] = data2


    count_ltail = 0
    count_utail = 0
    count_2tail = 0
    test_list = []
    s_performed = 0
    for s in range(S):
        if mean1+mean2 > 0:
#            perm = permFunc(perm) 
            if site_restricted: perm = site_restricted_permutation(perm) #TRI - I could have passed this in as permFunc, but I don't want to require the caller to know this
            else: perm = permFunc(perm) #TRI
            if not lib_str1:
                test_sample = testFunc(perm[:n1], perm[n1:])
            else: # case for lib strings
                test_sample = testFunc(perm) # how do I know how many counts are in cond1 or cond2? perm is a dict over lib strings (and conds?)
        else:
            test_sample = 0

        test_list.append(test_sample)
        if test_sample <= test_obs: count_ltail+=1
        if test_sample >= test_obs: count_utail+=1
        if abs(test_sample) >= abs(test_obs): count_2tail+=1

        s_performed+=1
        if adaptive:
            if s_performed == round(S*0.01) or s_performed == round(S*0.1) or s_performed == round(S*1):
                    if count_2tail >= round(S*0.01*0.10):
                        break



    pval_ltail = count_ltail/float(s_performed)
    pval_utail = count_utail/float(s_performed)
    pval_2tail = count_2tail/float(s_performed)

    return (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, test_list)








#

def cumulative_average(new_x, n, prev_avg):
    return ((new_x + (n*prev_avg))/(n+1.0), n+1)

#

def text_histogram(X, nBins = 20, resolution=200, obs = None):
    MIN = numpy.min(X)
    MAX = numpy.max(X)
    bin_list = numpy.linspace(MIN, MAX, nBins)
    hit_flag = "->"; empty_flag = "  "
    for b_l, b_u in zip(bin_list[:-2], bin_list[1:]):
        Z = numpy.logical_and(b_l <= X,  X < b_u)
        density = numpy.mean(Z)
        if obs != None and (b_l <= obs < b_u):
            flag = hit_flag
        else:
            flag  = empty_flag
        print("%-12f\t%s|%s" % (b_l, flag, "#"*int(resolution*density)))
    Z = numpy.logical_and(bin_list[-1] <= X,  X < float("inf"))
    density = numpy.mean(Z)
    if obs != None and (bin_list[-1] <= obs < float("inf")):
        flag = hit_flag
    else:
        flag = empty_flag
    print("%-12f\t%s|%s" % (bin_list[-1], flag,  "#"*int(resolution*density)))



def parse_lib_index(nData, libstr, nTAs):
    full_index = numpy.arange(nData)
    lib_to_index = {}
    for k,L in enumerate(libstr):
        if L not in lib_to_index: lib_to_index[L] = []
        lib_to_index[L] += list(full_index[k*nTAs:((k+1)*nTAs)])
    for L,index in lib_to_index.items():
        lib_to_index[L] = numpy.array(index)
    return lib_to_index

#

def combine_lib_dicts(L1, L2):
    KEYS = L1.keys()
    DATA = {}
    for K in KEYS:
        DATA[K] = numpy.array([L1[K], L2[K]],dtype=object)
    return DATA

# it looks like data1 is supposed to be pre-flattened (see parse_lib_index())

def get_lib_data_dict(data1, ctrl_lib_str, data2, exp_lib_str, nTAs):
    lib1_index_dict = parse_lib_index(len(data1), ctrl_lib_str, nTAs)
    lib2_index_dict = parse_lib_index(len(data2), exp_lib_str, nTAs)

    lib1_data_dict = dict([(L, data1[lib1_index_dict[L]]) for L in sorted(lib1_index_dict)])
    lib2_data_dict = dict([(L, data2[lib2_index_dict[L]]) for L in sorted(lib2_index_dict)])

    data_dict = combine_lib_dicts(lib1_data_dict, lib2_data_dict)
    return data_dict


#TEST-CASES

if __name__ == "__main__":

    """
    n = 20
    p = 0.5
    k = 14
    print("")
    print("#########################################")
    print("############ BINOM TEST #################")
    print("#########################################")
    print("Coin Tosses: %d" % n)
    print("Success Prob: %3.2f" % p)
    print("Observed: %d" % k)


    print("")
    print("Left-Tail Test:")
    print("%d tosses, p-value = %f" % (k, binom_test(k,n,p,"less")))

    print("")
    print("Right-Tail Test:")
    print("%d tosses, p-value = %f" % (k, binom_test(k,n,p,"greater")))


    print("")
    print("Two-Sided Test:")
    print("%d tosses, p-value = %f" % (k, binom_test(k,n,p,"two-sided")))



    print("")
    print("")
    print("#########################################")
    print("############ RESAMPLING #################")
    print("#########################################")

    data1 = scipy.stats.norm.rvs(100,10, size=1000)
    data2 = scipy.stats.norm.rvs(105,10, size=1000)

    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, test_list) = resampling(data1, data2, S=10000)
    print("Data1:")
    text_histogram(data1, nBins = 20)
    print("")
    print("Data2:")
    text_histogram(data2, nBins = 20)
    print("")
    print("Results:", (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail))
    print("")
    print("Resampling Histogram:"   )
    text_histogram(test_list, nBins = 20, obs=test_obs)

    """

    ## TEST
    import pytransit.transit_tools as transit_tools
    import pytransit.tnseq_tools as tnseq_tools
    import pytransit.norm_tools as norm_tools
    import sys

    ctrldata = ["/pacific/home/mdejesus/transit/tests/GI/H37Rv_day0_rep1.wig",
                "/pacific/home/mdejesus/transit/tests/GI/Rv2680_day0_rep1.wig", 
                "/pacific/home/mdejesus/transit/tests/GI/H37Rv_day0_rep2.wig",
                "/pacific/home/mdejesus/transit/tests/GI/Rv2680_day0_rep2.wig"]


    expdata = ["/pacific/home/mdejesus/transit/tests/GI/H37Rv_day32_rep1.wig", 
               "/pacific/home/mdejesus/transit/tests/GI/H37Rv_day32_rep2.wig",
               "/pacific/home/mdejesus/transit/tests/GI/H37Rv_day32_rep3.wig",
               "/pacific/home/mdejesus/transit/tests/GI/Rv2680_day32_rep1.wig",
               "/pacific/home/mdejesus/transit/tests/GI/Rv2680_day32_rep2.wig",
               "/pacific/home/mdejesus/transit/tests/GI/Rv2680_day32_rep3.wig"]

    annotation = "/pacific/home/mdejesus/transit/tests/GI/H37Rv.prot_table"

    i = 202
    if len(sys.argv) > 1:
        i = int(sys.argv[1])
    DO_LIB = True
    if len(sys.argv) > 2:
        DO_LIB = bool(int(sys.argv[2]))
    
    
    if DO_LIB:
        ctrl_lib_str = "ABAB"
        exp_lib_str = "AAABBB"
    else:
        ctrl_lib_str = ""
        exp_lib_str = ""
    
    Kctrl = len(ctrldata)
    Kexp  = len(expdata)

    (data, position) = transit_tools.get_validated_data(ctrldata+expdata)
    (K,N) = data.shape

    (data, factors) = norm_tools.normalize_data(data, "TTR", ctrldata+expdata, annotation)

    G = tnseq_tools.Genes(ctrldata + expdata, annotation, data=data, position=position)


    gene = G[i]

    print("\n\n")
    print("#"*100)
    print("#  (%s)  NEW TEST:   %s"  % (DO_LIB, gene))
    print("#"*100)
    print("")
       

 
    ii = numpy.ones(gene.n) == 1
        
    data1 = gene.reads[:Kctrl,ii].flatten() #TRI should we not flatten if doing site_restricted?
    data2 = gene.reads[Kctrl:,ii].flatten()
    
    data_dict = get_lib_data_dict(data1, ctrl_lib_str, data2, exp_lib_str, gene.n) # not used?

    if DO_LIB:
        (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  resampling(data1, data2, S=10000, testFunc=F_mean_diff_dict, permFunc=F_shuffle_dict_libraries, adaptive=False, lib_str1=ctrl_lib_str, lib_str2=exp_lib_str)
    else:
        (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  resampling(data1, data2, S=10000, testFunc=F_mean_diff_flat, permFunc=F_shuffle_flat, adaptive=False, lib_str1=ctrl_lib_str, lib_str2=exp_lib_str)
        
    print("Resampling Histogram:")
    text_histogram(testlist, nBins = 20, obs=test_obs)
    
     
        
