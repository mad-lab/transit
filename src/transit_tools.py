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



def aton(aa):
    return(((aa-1)*3)+1)

def parseCoords(strand, aa_start, aa_end, start, end):
    if strand == "+":
        return((aton(aa_start) + start,  aton(aa_end) + start))
    # Coordinates are Reversed... to match with Trash FILE TA coordinates
    if strand == "-":
        return((end - aton(aa_end), end - aton(aa_start)))



def get_wig_stats(path):
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        if line.startswith("variable"): continue
        if line.startswith("location"): continue
        tmp = line.split()
        pos = int(tmp[0])
        rd = float(tmp[1])
        data.append(rd)
    return (sum(data), sum([1 for rd in data if rd > 0])/float(len(data)),  sum(data)/float(len(data)), max(data))


def get_reads(path):
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        if line.startswith("variable"): continue
        if line.startswith("location"): continue
        tmp = line.split()
        rd = int(tmp[1])
        data.append(rd)
    return data


def fetch_name(filepath):
    return os.path.splitext(ntpath.basename(filepath))[0]




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


def normalize_data(data):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs_hit = numpy.sum(data > 0,1)
    mean_hits = total_hits/TAs_hit
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    factors = numpy.zeros((K,1))
    factors[:,0] = grand_mean/mean_hits
    return factors * data
    

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

            if strand == "+":
                #Ignore TAs at stop codon
                if ignoreCodon and coord > end-3:
                    continue

            else:
                #ignore TAs at stop codon
                if ignoreCodon and coord < start + 3:
                    continue


            #Ignore TAs at beginning n%
            if (coord-start)/float(end-start) <= (ignoreNTerm/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and NTerm", ignoreNTerm/100.0
                continue
            
            #Ignore TAs at end c%
            if (coord-start)/float(end-start) >= ((100-ignoreCTerm)/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and Cterm", ignoreCTerm/100.0
                continue


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
    






