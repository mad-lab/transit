import sys
import os
import math
import ntpath
import numpy


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
        rd = int(tmp[1])
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



def get_norm_factors(data):
    (K,N) = data.shape
    total_hits = numpy.sum(data,1)
    TAs_hit = numpy.sum(data > 0, 1)
    mean_hits = total_hits/TAs_hit
    grand_total = numpy.sum(mean_hits)
    grand_mean = grand_total/float(K)
    factors = numpy.zeros((K,1))
    factors[:,0] = grand_mean/mean_hits

    return factors

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


            #Ignore TAs at beginning 5%
            if (coord-start)/float(end-start) <= (ignoreNTerm/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and NTerm", ignoreNTerm/100.0
                continue
            
            #Ignore TAs at end 5%
            if (coord-start)/float(end-start) >= ((100-ignoreCTerm)/100.0):
                #print "Ignoring", coord, "from gene", gene, "with", (coord-start)/float(end-start), "perc and Cterm", ignoreCTerm/100.0
                continue


            orf2reads[gene].append(data[:,i])
            orf2pos[gene].append(position[i])
    return (orf2reads, orf2pos)
    

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
    






