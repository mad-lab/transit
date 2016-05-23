import math
import numpy
import scipy.stats

def transformToRange(X, new_min, new_max, old_min=None, old_max=None):

    if old_min == None:
        old_min = min(X)
    if old_max == None:
        old_max = max(X)
    
    old_range = old_max - old_min

    new_range = new_max - new_min
    return [float(x - old_min) / old_range * new_range + new_min for x in X]




def fact(n):
    if n == 0: return (1)
    else: return reduce(lambda x,y: x*y, range(1,n+1))


def comb1(n,k):
    prod = 1
    for i in range(1,k+1):
        prod = prod * (n - (k-i))/float(i)
    return(prod)


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


def norm(x, mu,sigma):
    """Normal distribution"""
    sigma = float(sigma)
    return(1/(sigma*(math.sqrt(2*math.pi))) * math.exp( -0.5 * math.pow( (x-mu)/sigma,2)))


def binom(k,n,p):
    """Binomial distribution. Uses Normal approximation for large 'n' """
    if n >= 100:
        return(norm(k, n*p, math.sqrt(n*p*(1-p)) ) )
    else:
        return(comb(n,k) * math.pow(p, k) * math.pow(1-p, n-k))


def binom_cdf(k,n,p):
    """CDF of the binomial distribution"""
    return(sum([binom(i,n,p) for i in range(0,k+1)]))


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


def boxcoxtransform(x, lambdax):
    """
    Performs a box-cox transformation to data vector X.
    WARNING: elements of X should be all positive! 
    Fixed: '>' has changed to '<'
    """
    if x <= 0:
        raise ArgumentError, "Nonpositive value(s) in X vector"
    if abs(lambdax) < 1.0e-5:
        return(math.log(x))
    else:
        return((x**lambdax - 1.0)/lambdax)
    #return math.log(x) if abs(lambdax) < 1.0e-5 else (x**lambdax - 1.0)/lambdax

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

def BH_fdr_correction(X):
    """Adjusts p-values using the Benjamini Hochberg procedure"""
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

def bayesian_ess_thresholds(Z_raw, ALPHA=0.05):
    """Returns Essentiality Thresholds using a BH-like procedure"""
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
        #if count < 30: print i, wi, ai_n, delta_w
        if delta_w > ai_n:
            ess_threshold = Z[i-1]
            #print "i", i
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
        #print count
        #if count < 20:
        #   print i, wi, ai_n, mean_wi, delta_w, N-i+1, N-1, W[N-i-1], W[i-1]

        if ai_n > delta_w:
        #   print i, wi, ai_n, mean_wi, delta_w, N-i+1, N-1, W[N-i-1], W[i-1]
            break
        noness_threshold = Z[N-i]

    return(ess_threshold, noness_threshold)


def tricube(X):
    #TODO: Write docstring
    result = numpy.zeros(len(X))
    ii = numpy.logical_and(X >= -1, X <= 1)
    result[ii] = numpy.power(1 - numpy.power(numpy.abs(X[ii]), 3), 3)
    return result


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


def loess_correction(X, Y, h=10000, window=100):
    #TODO: Write docstring
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



def F_mean_diff_flat(A, B):
    return numpy.mean(B) - numpy.mean(A)
    
def F_sum_diff_flat(A, B):
    return numpy.sum(B) - numpy.sum(A)


def F_shuffle_flat(X):
    return numpy.random.permutation(X)

def resampling(data1, data2, S=10000, testFunc=F_mean_diff_flat,
            permFunc=F_shuffle_flat, adaptive=False):
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
    count_ltail = 0
    count_utail = 0
    count_2tail = 0

    test_list = []

    n1 = len(data1)
    n2 = len(data2)

    mean1 = 0
    if n1 > 0:
        mean1 = numpy.mean(data1)
    mean2 = 0
    if n2 > 0:
        mean2 = numpy.mean(data2)

    test_obs = testFunc(data1, data2)

    try:
        norm_mean1 = mean1/float(n1) if mean1 > 0.0 else 1.0
        norm_mean2 = mean2/float(n2) if mean2 > 0.0 else 1.0
        log2FC = math.log(norm_mean2/norm_mean1,2)
    except:
        log2FC = 0

    perm = numpy.zeros(n1+n2)
    perm[:n1] = data1
    perm[n1:] = data2

    s_performed = 0
    for s in range(S):
        if len(perm) >0:
            perm = permFunc(perm)
            test_sample = testFunc(perm[:n1], perm[n1:])
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






#TEST-CASES

if __name__ == "__main__":
    sdata = """
.15 .09 .18 .10 .05 .12 .08
.05 .08 .10 .07 .02 .01 .10
.10 .10 .02 .10 .01 .40 .10
.05 .03 .05 .15 .10 .15 .09
.08 .18 .10 .20 .11 .30 .02
.20 .20 .30 .30 .40 .30 .05 
"""
    X = [float(x) for x in sdata.split()]
    out = boxcoxTable(X, -2, 2, 0.10)
    out.sort()
    print "Log likelihood function:"
    for (llik, lambdax) in out:
       print llik, lambdax
    print "best lambda = ", out[-1][1]
    print "with loglik function value = ",out[-1][0]




    print ""
    print ""
    print "#########################################"
    print "############ BINOM TEST #################"
    print "#########################################"
    print "DEFAULT"
    n = 15
    p = 0.5
    for i in range(0,16):
        print i, binom(i,n,p), binom_cdf(i,n,p)


    print ""
    print "LESS"    
    for i in range(0,16):
        print i, binom_test(i,n,p,"less")

    print ""
    print "GREATER"
    for i in range(0,16):
        print i, binom_test(i,n,p,"greater")


    print ""
    print "TWO-SIDED"
    for i in range(0,16):
        print i, binom_test(i,n,p,"two-sided")



    print ""
    print ""
    print "#########################################"
    print "############ RESAMPLING #################"
    print "#########################################"

    data1 = scipy.stats.norm.rvs(100,10, size=100)
    data2 = scipy.stats.norm.rvs(105,10, size=100)


    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail) = resampling(data1, data2, S=10000)
    print "Data1:", data1
    print "Data2:", data2
    print "Results:", (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail)

