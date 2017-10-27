import sys
import numpy
import scipy.stats
import scipy.optimize
import warnings

import tnseq_tools


class NormMethod:
    name = "undefined"
    @staticmethod
    def normalize():
        raise NotImplemented

class NZMeanNorm(NormMethod):
    name = "nzmean"

    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        """Returns the normalization factors for the data, using the NZMean method.

        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.

        Returns:
            numpy array: Array with the normalization factors for the nzmean method.

        :Example:
            >>> import pytransit._tools.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> factors = norm_tools.nzmean_factors(data)
            >>> print factors
            array([[ 1.14836149],
                   [ 0.88558737]])

        .. seealso:: :class:`normalize_data`

        """
        (K,N) = data.shape
        total_hits = numpy.sum(data,1)
        TAs_hit = numpy.sum(data > 0, 1)
        mean_hits = total_hits/TAs_hit
        grand_total = numpy.sum(mean_hits)
        grand_mean = grand_total/float(K)
        factors = numpy.zeros((K,1))
        factors[:,0] = grand_mean/mean_hits
        data = factors * data
        return (data, factors)



class TotReadsNorm(NormMethod):
    name = "totreads"

    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        """Returns the normalization factors for the data, using the total reads
        method.

        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.

        Returns:
            numpy array: Array with the normalization factors for the totreads method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> factors = norm_tools.totreads_factors(data)
            >>> print factors
            array([[ 1.2988762],
                   [ 0.8129396]])

        .. seealso:: :class:`normalize_data`

        """
        (K,N) = data.shape
        total_hits = numpy.sum(data,1)
        TAs = float(N)
        mean_hits = total_hits/TAs
        grand_total = numpy.sum(mean_hits)
        grand_mean = grand_total/float(K)
        factors = numpy.zeros((K,1))
        factors[:,0] = grand_mean/mean_hits
        data = factors * data
        return (data, factors)


class TTRNorm(NormMethod):
    name = "emphist"

    def empirical_theta(X):
        """Calculates the observed density of the data.

        This is used as an estimate insertion density by some normalization methods.
        May be improved by more sophisticated ways later on.

        Arguments:
            data (numpy array): (N) numpy array defining read-counts at N sites.

        Returns:
            float: Density of the given dataset.

        :Example:
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> import pytransit.norm_tools as norm_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> theta = norm_tools.empirical_theta(data)
            >>> print theta
            0.467133570136


        .. seealso:: :class:`TTR_factors`
        """
        return numpy.mean(X > 0)

    def trimmed_empirical_mu(X, t=0.05):
        """Estimates the trimmed mean of the data.

        This is used as an estimate of mean count by some normalization methods.
        May be improved by more sophisticated ways later on.

        Arguments:
            data (numpy array): (N) numpy array defining read-counts at N sites.
            t (float): Float specifying fraction of start and end to trim.

        Returns:
            float: (Trimmed) Mean of the given dataset.

        :Example:
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> import pytransit.norm_tools as norm_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> mu = norm_tools.trimmed_empirical_mu(data)
            >>> print mu
            120.73077107

        .. seealso:: :class:`TTR_factors`
        """
        return scipy.stats.trim_mean(X[X > 0], t)


    @staticmethod
    def normalize(data, wigList=[], annotationPath="", thetaEst=empirical_theta, muEst=trimmed_empirical_mu, target=100.0):
        """Returns the normalization factors for the data, using the TTR method.


        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.
            thetaEst (function): Function used to estimate density. Should take a list
                of counts as input.
            muEst (function): Function used to estimate mean count. Should take a list
                of counts as input.

        Returns:
            numpy array: Array with the normalization factors for the TTR method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> factors = norm_tools.TTR_factors(data)
            >>> print factors
            array([[ 1.        ],
                   [ 0.62862886]])

        .. seealso:: :class:`normalize_data`
        """
        K = len(data)
        N = len(data[0])

        factors = numpy.zeros((K,1))
        for j in range(K):
            factors[j] = float(target)/(thetaEst(data[j]) * muEst(data[j]))
        data = factors * data
        return (data, factors)


class EmpHistNorm(NormMethod):
    name = "emphist"

    @staticmethod
    def Fzinfnb(params, args):
        """Objective function for the zero-inflated NB method."""
        pi, mu, r = params
        Fdata = args
        temp0 = numpy.nan_to_num(numpy.log(pi + scipy.stats.nbinom.pmf(Fdata[Fdata==0], mu, r)))
        tempnz = numpy.nan_to_num(numpy.log(1.0-pi)+scipy.stats.nbinom.logpmf(Fdata[Fdata>0], mu, r))
        negLL = -(numpy.sum(temp0) + numpy.sum(tempnz))
        return negLL

    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        """Returns the normalized data, using the empirical hist method.

        Arguments:
            wigList (list): List of paths to wig formatted datasets.
            annotationPath (str): Path to annotation in .prot_table or GFF3 format.

        Returns:
            numpy array: Array with the normalization factors for the emphist method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> factors = norm_tools.emphist_factors(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"], "transit/genomes/H37Rv.prot_table")
            >>> print factors
            array([[ 1.        ],
                   [ 0.63464722]])

        .. seealso:: :class:`normalize_data`
        """

        G = tnseq_tools.Genes(wigList, annotationPath)
        K = len(wigList)
        temp = []
        for j in range(K):
            reads_per_gene = []
            for gene in G:
                tempdata = numpy.array(gene.reads)
                if len(tempdata[0]) > 0:
                    reads_per_gene.append(numpy.sum(tempdata[j,:]))
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

        data = factors * data
        return (data, factors)


class AdaptiveBGCNorm(NormMethod):
    name = "aBGC"

    def ecdf(S, x):
        """Calculates an empirical CDF of the given data."""
        return numpy.sum(S<=x)/float(len(S))

    def cleaninfgeom(x, rho):
        """Returns a 'clean' output from the geometric distribution."""
        if x == float('inf'):
            return scipy.stats.geom.ppf(0.9999999999999999, rho)
        else:
            return x

    @staticmethod
    def normalize(data, wigList=[], annotationPath="", doTotReads = True, bgsamples = 200000):
        """Returns the normalized data using the aBGC method.


        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.
            doTotReads (bool):  Boolean specifying whether to do TTR normalization as well.
            bgsamples (int): Integeer specifying how many samples to take.

        Returns:
            numpy array: Array with the normalized data.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> normdata = norm_tools.aBGC_norm(data)
            >>> print normdata
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

        .. seealso:: :class:`normalize_data`
        """

        K,N = data.shape
        norm_data = numpy.zeros(data.shape)
        S = bgsamples
        F = [i/100.0 for i in range(0,31) if i % 2 == 0]
        BGC = []
        param_list = []
        bgc_factors = []
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
                    bgc_factors.append((rho, Kp))
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
            (norm_data, factors) = TTRNorm.normalize(norm_data)
        return (norm_data, bgc_factors)



class ZeroInflatedNBNorm(NormMethod):
    name = "zinfb"

    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        """Returns the normalization factors for the data using the zero-inflated
        negative binomial method.


        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.

        Returns:
            numpy array: Array with the normalization factors for the zinfnb method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> factors = norm_tools.zinfnb_factors(data)
            >>> print factors
            [[ 0.0121883 ]
             [ 0.00747111]]

        .. seealso:: :class:`normalize_data`
        """
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
        data = factors * data
        return (data, factors)


class QuantileNorm(NormMethod):
    name = "quantile"

    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        """Performs Quantile Normalization as described by Bolstad et al. 2003

        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.

        Returns:
            numpy array: Array with the data normalized by the quantile normalization method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> normdata = norm_tools.quantile_norm(data)
            >>> print normdata

        .. seealso:: :class:`normalize_data`

        """
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
        return (norm_data, numpy.ones(1))





class BetaGeomNorm(NormMethod):
    name = "betageom"

    def ecdf(S, x):
        """Calculates an empirical CDF of the given data."""
        return numpy.sum(S<=x)/float(len(S))

    def cleaninfgeom(x, rho):
        """Returns a 'clean' output from the geometric distribution."""
        if x == float('inf'):
            return scipy.stats.geom.ppf(0.9999999999999999, rho)
        else:
            return x

    @staticmethod
    def normalize(data, wigList=[], annotationPath="", doTTR = True, bgsamples=200000):
        """Returns normalized data according to the BGC method.

        Arguments:
            data (numpy array): (K,N) numpy array defining read-counts at N sites
                for K datasets.
            doTTR (bool): Boolean specifying whether to do TTR norm as well.
            bgsamples (int): Integer specifying how many samples to take.

        Returns:
            numpy array: Array with the data normalized using the betageom method.

        :Example:
            >>> import pytransit.norm_tools as norm_tools
            >>> import pytransit.tnseq_tools as tnseq_tools
            >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
            >>> print data
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                   [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
            >>> normdata = norm_tools.betageom_norm(data)
            >>> print normdata
            [[ 0.  0.  0. ...,  0.  0.  0.]
             [ 0.  0.  0. ...,  0.  0.  0.]]

        .. seealso:: :class:`normalize_data`
        """

        (K,N) = data.shape
        total_hits = numpy.sum(data,1)
        TAs_hit = numpy.sum(data > 0,1)
        mean_hits = total_hits/TAs_hit
        grand_total = numpy.sum(mean_hits)
        grand_mean = grand_total/float(K)
        norm_data = numpy.zeros(data.shape)
        bgc_factors = []
        for j in range(K):

            tQ = numpy.arange(0,N)/float(N)
            eX = numpy.array([rd for rd in data[j]])
            eX.sort()

            rho = max(1.0/scipy.stats.trim_mean(eX+1, 0.001), 0.0001)
            A = (numpy.sum(numpy.power(numpy.log(1.0-tQ),2)))/(numpy.sum(eX*numpy.log(1.0-tQ)))
            Kp = max((2.0 * numpy.exp(A) - 1)   /(numpy.exp(A) + rho - 1), 10)

            bgc_factors.append((rho,Kp))
            try:
                BGsample = scipy.stats.geom.rvs(scipy.stats.beta.rvs(Kp*rho, Kp*(1-rho), size=bgsamples), size=bgsamples)
            except Exception as e:
                print "BGC ERROR with rho=%f, Kp=%f, A=%s" % (rho, Kp, A)
                print str(e)
                BGsample = scipy.stats.geom.rvs(rho, size=bgsamples)

            for i in range(N):
                norm_data[j,i] = cleaninfgeom(scipy.stats.geom.ppf(ecdf(BGsample, data[j,i]), 1.0/grand_mean), 1.0/grand_mean)

        if doTTR:
            (norm_data, factors) = TTRNorm.normalize(norm_data)
        return (norm_data, bgc_factors)


class NoNorm(NormMethod):
    name = "nonorm"
    @staticmethod
    def normalize(data, wigList=[], annotationPath=""):
        return (data, numpy.ones(1))


methods = {}
methods["nonorm"] = NoNorm
methods["TTR"] = TTRNorm
methods["nzmean"] = NZMeanNorm
methods["totreads"] = TotReadsNorm
methods["betageom"] = BetaGeomNorm
methods["zinfnb"] = ZeroInflatedNBNorm
methods["quantile"] = QuantileNorm
methods["aBGC"] = AdaptiveBGCNorm
methods["emphist"] = EmpHistNorm


#########################
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

    :Example:
        >>> import pytransit.norm_tools as norm_tools
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        (normdata, normfactors) = norm_tools.normalize_data(data, "TTR")   # Some methods require annotation and path to wig files.
        >>> print normfactors
        array([[ 1.        ],
               [ 0.62862886]])
        >> print normdata
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

    .. note:: Some normalization methods require the wigList and annotationPath arguments.

    """
    factors = []
    if method in methods:
        return methods[method].normalize(data, wigList, annotationPath)
    else:
        warnstr = "Normalization method '%s' is unknown. Read-counts were not normalized." % (method)
        warnings.warn(warnstr)
    return methods["nonorm"].normalize(data, wigList, annotationPath)


def empirical_theta(X):
    """Calculates the observed density of the data.

    This is used as an estimate insertion density by some normalization methods.
    May be improved by more sophisticated ways later on.

    Arguments:
        data (numpy array): (N) numpy array defining read-counts at N sites.

    Returns:
        float: Density of the given dataset.

    :Example:
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> import pytransit.norm_tools as norm_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> theta = norm_tools.empirical_theta(data)
        >>> print theta
        0.467133570136


    .. seealso:: :class:`TTR_factors`
    """
    return numpy.mean(X > 0)

def trimmed_empirical_mu(X, t=0.05):
    """Estimates the trimmed mean of the data.

    This is used as an estimate of mean count by some normalization methods.
    May be improved by more sophisticated ways later on.

    Arguments:
        data (numpy array): (N) numpy array defining read-counts at N sites.
        t (float): Float specifying fraction of start and end to trim.

    Returns:
        float: (Trimmed) Mean of the given dataset.

    :Example:
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> import pytransit.norm_tools as norm_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> mu = norm_tools.trimmed_empirical_mu(data)
        >>> print mu
        120.73077107

    .. seealso:: :class:`TTR_factors`
    """

    return scipy.stats.trim_mean(X[X > 0], t)


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
    negative binomial method.


    Arguments:
        data (numpy array): (K,N) numpy array defining read-counts at N sites
            for K datasets.

    Returns:
        numpy array: Array with the normalization factors for the zinfnb method.

    :Example:
        >>> import pytransit.norm_tools as norm_tools
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> factors = norm_tools.zinfnb_factors(data)
        >>> print factors
        [[ 0.0121883 ]
         [ 0.00747111]]

    .. seealso:: :class:`normalize_data`
    """
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
    return numpy.array(factors)

#

def ecdf(S, x):
    """Calculates an empirical CDF of the given data."""
    return numpy.sum(S<=x)/float(len(S))


def cleaninfgeom(x, rho):
    """Returns a 'clean' output from the geometric distribution."""
    if x == float('inf'):
        return scipy.stats.geom.ppf(0.9999999999999999, rho)
    else:
        return x

#

def norm_to_target(data, target):
    """Returns factors to normalize the data to the given target value.

    Arguments:
        data (numpy array): (K,N) numpy array defining read-counts at N sites
            for K datasets.
        target (float): Floating point specifying the target for the mean of the data/

    Returns:
        numpy array: Array with the factors necessary to normalize mean to target.

    :Example:
        >>> import pytransit.norm_tools as norm_tools
        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_data(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])
        >>> factors = norm_tools.norm_to_target(data, 100)
        >>> print factors
        [[ 1.8548104 ]
         [ 1.16088726]]


    .. seealso:: :class:`normalize_data`
    """
    (K,N) = data.shape
    factors = numpy.zeros((K,1))
    factors[:,0] = float(target)/numpy.mean(data,1)
    return factors
