import sys
import math
import warnings
import numpy
import scipy.stats
from functools import total_ordering

try:
    import norm_tools
    noNorm = False
except ImportError:
    noNorm = True
    warnings.warn("Problem importing the norm_tools.py module. Read-counts will not be normalized.")

@total_ordering
class Gene:
    """This is a longer explanation, which may include math with latex syntax
    Then, you need to provide optional subsection in this order (just to be
    consistent and have a uniform documentation. Nothing prevent you to
    switch the order):

        - parameters using ``:param <name>: <description>``
        - type of the parameters ``:type <name>: <description>``
        - returns using ``:returns: <description>``
        - examples (doctest)
        - seealso using ``.. seealso:: text``
        - notes using ``.. note:: text``
        - warning using ``.. warning:: text``
        - todo ``.. todo:: text``


    Attributes:
        orf (str): ORF ID. Must be unique.
        name (str): Human readable name of the ORF.
        reads (list): Read-counts data for the ORF.
        position (list): Position of TA sites for the ORF.
        start (int): Start coordinate for the ORF.
        end (int): End coordinate for the ORF.
        strand (str): Strand for the ORF.

    :Example:

        >>> import tnseq_tools
        >>> G = tnseq_tools.Gene("Rv0001", "dnaA", [[0,0,0,0,1,0,32]], start=1, end=1500, strand="+")
        >>> print G
        Rv0001 (dnaA): 2,7,4

        .. warning:: orf must be unique.
        .. seealso:: :class:`Genes`
    """

    def __init__(self, orf, name, desc, reads, position, start=0, end=0, strand=""):
        """Initializes the Gene object."""

        self.orf = orf
        self.name = name
        self.desc = desc
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = numpy.array(reads)
        self.position = numpy.array(position, dtype=int)
        self.tosses = tossify(self.reads)
        try:
            self.runs = runs(self.tosses)
        except Exception as e:
            print orf, name, self.tosses
            raise e

        self.k = int(numpy.sum(self.tosses))
        self.n = len(self.tosses)
        try:
            self.r = numpy.max(self.runs)
        except Exception as e:
            print orf, name, self.tosses
            print self.runs
            raise e

        self.s = self.get_gap_span()
        self.t = self.get_gene_span()

    def __getitem__(self, i):
        """Return read-counts at position i."""
        return self.reads[:, i]

    def __str__(self):
        """Return a string representation of the object."""
        return "%s\t(%s)\tk=%d\tn=%d\tr=%d\ttheta=%1.5f" % (self.orf, self.name, self.k, self.n, self.r, self.theta())

    def __eq__(self, other):
        return self.orf == other.orf

    def __lt__(self, other):
        return self.orf < other.orf 

    def get_gap_span(self):
        """Returns the span of the maxrun of the gene (i.e. number of nucleotides)."""
        if len(self.position) > 0:
            if self.r == 0:
                return 0
            index = runindex(self.runs)
            #maxii = numpy.argmax(self.runs)
            maxii = numpy.argwhere(self.runs == numpy.max(self.runs)).flatten()[-1]
            runstart = index[maxii]
            runend = runstart + max(self.runs) - 1
            return self.position[runend] - self.position[runstart] + 2
        else:
            return 0
        
    
    def get_gene_span(self):
        """Returns the number of nucleotides spanned by the gene."""
        if len(self.position) > 0:
            return self.position[-1] - self.position[0] + 2
        return 0


    def theta(self):
        """Return the insertion density ("theta") for the gene."""
        if self.n:
            return float(self.k)/self.n
        else:
            return 0.0

    def phi(self):
        """ Return the non-insertion density ("phi") for the gene."""
        return 1.0 - self.theta()

    def total_reads(self):
        """ Return the total reads for the gene."""
        return numpy.sum(self.reads, 1)


    def calculate_span(self):
        """Caclulates the span based on the coordinates"""
        # TODO: Check if it works.
        runs = self.runs()
        if len(self.raw_data) > 0:
            runs = self.runs(include_pos=True) or [(0,0)]
            maxrun = max(runs)
            return self.get_TA_coord(maxrun[1] + maxrun[0]-1) + 2 - self.get_TA_coord(maxrun[1]) 
        else:
            return -1

    def calculate_length(self, raw_data):
        """Caclulates the length based on the coordinates"""
        # TODO: Check if it works.
        if len(self.raw_data) > 0:
            return self.raw_data[-1][0] + 2 - self.raw_data[0][0]
        else:
            return -1



class Genes:

    def __getitem__(self, i):
        """Defines __getitem__ method so that it works as dictionary and list."""
        if isinstance(i, int):
            return(self.genes[i])

        if isinstance(i, basestring):
            return self.genes[self.orf2index[i]]


    def __contains__(self, item):
        """Defines __contains__ to check if gene exists in the list."""
        return item in self.orf2index


    def __len__(self):
        """Defines __len__ returning number of genes."""
        return len(self.genes)


    def __str__(self):
        """Defines __str__ to print a generic str with the size of the list."""
        return "Genes Object (N=%d)" % len(self.genes)
    
    def __init__(self, wigList, protTable, norm="nonorm", reps="All", minread=1, ignoreCodon = True, nterm=0.0, cterm=0.0, include_nc = False, data=[], position=[]):
        """Initializes the gene list based on the list of wig files and a prot_table."""
        self.wigList = wigList
        self.protTable = protTable
        self.norm = norm
        self.reps = reps
        self.minread = minread
        self.ignoreCodon = ignoreCodon
        self.nterm = nterm
        self.cterm = cterm
        self.include_nc = include_nc


        self.orf2index = {}
        self.genes = []
        
        orf2info = get_gene_info(self.protTable)
        if not numpy.any(data):
            (data, position) = get_data(self.wigList)
        hash = get_pos_hash(self.protTable)
       

        if not noNorm:
            (data, factors) = norm_tools.normalize_data(data, norm, self.wigList, self.protTable)
        else:
            factors = []
       
        if reps.lower() != "all":
            data = numpy.array([combine_replicates(data, method=reps)])

        K,N = data.shape
        
        orf2posindex = {}
        visited_list = []
        for i in range(N):
            genes_with_coord = hash.get(position[i], [])
            for gene in genes_with_coord:
                if gene not in orf2posindex: visited_list.append(gene)
                if gene not in orf2posindex: orf2posindex[gene] = []

                name,desc,start,end,strand = orf2info.get(gene, ["", "", 0, 0, "+"])
                
                if strand == "+":
                    if self.ignoreCodon and position[i] > end - 3:
                        continue
                else:
                    if self.ignoreCodon and position[i] < start + 3:
                        continue

                if (position[i]-start)/float(end-start) < (self.nterm/100.0):
                    continue
                
                if (position[i]-start)/float(end-start) > ((100-self.cterm)/100.0):
                    continue

                orf2posindex[gene].append(i)

        count = 0
        for line in open(self.protTable):
            tmp = line.split("\t")
            gene = tmp[8]
            name,desc,start,end,strand = orf2info.get(gene, ["", "", 0, 0, "+"])
            posindex = orf2posindex.get(gene, [])
            if posindex:
                pos_start = orf2posindex[gene][0]
                pos_end = orf2posindex[gene][-1]
                self.genes.append(Gene(gene, name, desc, data[:, pos_start:pos_end+1], position[pos_start:pos_end+1], start, end, strand))
            else:
                self.genes.append(Gene(gene, name, desc, numpy.array([[]]), numpy.array([]), start, end, strand))
            self.orf2index[gene] = count
            count += 1


    def local_insertions(self):
        """Returns numpy array with the number of insertions, 'k', for each gene."""
        G = len(self.genes)
        K = numpy.zeros(G)
        for i in xrange(G):
            K[i] = self.genes[i].k
        return K


    def local_sites(self):
        """Returns numpy array with total number of TA sites, 'n', for each gene."""
        G = len(self.genes)
        N = numpy.zeros(G)
        for i in range(G):
            N[i] = self.genes[i].n
        return N


    def local_runs(self):
        """Returns numpy array with maximum run of non-insertions, 'r', for each gene."""
        G = len(self.genes)
        R = numpy.zeros(G)
        for i in xrange(G):
            R[i] = self.genes[i].r
        return R 

    def local_gap_span(self):
        """Returns numpy array with the span of nucleotides of the largest gap,
        's', for each gene."""
        G = len(self.genes)
        S = numpy.zeros(G)
        for i in xrange(G):
            S[i] = self.genes[i].s
        return S
  
    def local_gene_span(self):
        """Returns numpy array with the span of nucleotides of the gene,
        't', for each gene."""
        G = len(self.genes)
        T = numpy.zeros(G)
        for i in xrange(G):
            T[i] = self.genes[i].t
        return T

    def local_reads(self):
        """Returns numpy array of lists containing the read counts for each gene."""
        all_reads = []
        G = len(self.genes)
        for i in xrange(G):
            all_reads.extend(self.genes[i].reads)
        return numpy.array(all_reads)


    def local_thetas(self):
        """Returns numpy array of insertion frequencies, 'theta', for each gene."""
        G = len(self.genes)
        theta = numpy.zeros(G)
        for i in xrange(G):
            theta[i] = self.genes[i].theta()
        return theta


    def local_phis(self):
        """Returns numpy array of non-insertion frequency, 'phi', for each gene."""
        return 1.0 - self.theta()


    ######

    def global_insertion(self):
        """Returns total number of insertions, i.e. sum of 'k' over all genes."""
        G = len(self.genes)
        total = 0
        for i in xrange(G):
            total += self.genes[i].k
        return total

    def global_sites(self):
        """Returns total number of sites, i.e. sum of 'n' over all genes."""
        G = len(self.genes)
        total = 0
        for i in xrange(G):
            total += self.genes[i].n
        return total

    def global_run(self):
        """Returns the run assuming all genes were concatenated together."""
        return maxrun(self.tosses())

    def global_reads(self):
        """Returns the reads among the library."""
        return self.data

    def global_theta(self):
        """Returns global insertion frequency, of the library."""
        return float(self.global_insertion())/self.global_sites()

    def global_phi(self):
        """Returns global non-insertion frequency, of the library."""
        return 1.0 - self.global_theta()

    def total_reads(self):
        """Returns total reads among the library."""
        reads_total = 0
        for g in self.genes:
            reads_total += g.total_reads()
        return reads_total

    def tosses(self):
        """Returns list of bernoulli trials, 'tosses', representing insertions in the gene."""
        all_tosses = []
        for g in self.genes:
            all_tosses.extend(g.tosses)
        return all_tosses




def tossify(data):
    """Reduces the data into Bernoulli trials (or 'tosses') based on whether counts were observed or not."""
    K,N = data.shape
    reduced = numpy.sum(data,0)
    return numpy.zeros(N) + (numpy.sum(data, 0) > 0)

def runs(data):
    """Return list of all the runs of consecutive non-insertions."""
    runs = []
    current_r = 0
    for read in data:
        if read > 0: # If ending a run of zeros
            if current_r > 0: # If we were in a run, add to list
                runs.append(current_r)
            current_r = 0
            runs.append(current_r)
        else:
            current_r += 1
    # If we ended in a run, add it
    if current_r > 0:
        runs.append(current_r)

    if not runs:
        return [0]
    return runs

def runindex(runs):
    """Returns a list of the indexes of the start of the runs; complements runs()."""
    index = 0
    index_list = []
    runindex = 0
    for r in runs:
        for i in range(r):
            if i == 0:
                runindex = index
            index+=1
        if r == 0:
            runindex = index
            index+=1
        index_list.append(runindex)
    return index_list


def get_data(wig_list):
    """ Returns a tuple of (data, position) containing a matrix of raw read counts, and list of coordinates. """
    K = len(wig_list)
    T = 0

    if not wig_list:
        return (numpy.zeros((1,0)), numpy.zeros(0))

    for line in open(wig_list[0]):
        if line[0] not in "0123456789": continue
        T+=1

    data = numpy.zeros((K,T))
    position = numpy.zeros(T)
    for j,path in enumerate(wig_list):
        reads = []
        i = 0
        for line in open(path):
            if line[0] not in "0123456789": continue
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



def get_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate."""
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

def get_coordinate_map(galign_path, reverse=False):
    c1Toc2 = {}
    for line in open(galign_path):
        if line.startswith("#"): continue
        tmp = line.split()
        star = line.strip().endswith("*")
        leftempty = tmp[0].startswith("-")
        rightempty = tmp[1].endswith("-")
        if leftempty:
            left = -1
        else:
            left = int(tmp[0])
        if rightempty:
            right = -1
        elif leftempty:
            right = int(tmp[1])
        else:
             right = int(tmp[2])
            
        if not reverse:
            if not leftempty:
                c1Toc2[left] = right
        else:
            if not rightempty:
                 c1Toc2[right] = left    
    return c1Toc2

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


def getGamma():
    """Euler-Mascheroni constant ~ 0.577215664901 """
    return(0.5772156649015328606)

    
def ExpectedRuns(n,p):
    """Expected value of the run of non=insertions (Schilling, 1990):
    
        ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n)

    """   
    q = 1-p
    gamma = getGamma()
    r1 = getR1(n)
    E1 = getE1(n)
    A = math.log(n*q,1.0/p)
    B = gamma/math.log(1.0/p)
    ER = A + B -0.5 + r1 + E1
    return ER 
    

def VarR(n,p):
    """Variance of the expected run of non-insertons (Schilling, 1990):
 
        VarR_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n)

    """
    r2 = getR2(n)
    E2 = getE2(n)    
    A = math.pow(math.pi,2.0)/(6* math.pow(math.log(1.0/p),2.0))
    V = A + 1/12.0 + r2 + E2
    return V

    
def GumbelCDF(x,u,B):
    """CDF of the Gumbel distribution:

        e^(-e^( (u-x)/B))

    """
    return (math.exp( -1 * math.exp((u-x)/B )))
    

def griffin_analysis(genes_obj, pins):
    """Implements the basic Gumbel analysis of runs of non-insertion, described
     in Griffin et al. 2011.

    This analysis method calculates a p-value of observing the maximun run of
    TA sites without insertions in a row (i.e. a "run", r). Unusually long
    runs are indicative of an essential gene or protein domain. Assumes that
    there is a constant, global probability of observing an insertion
    (tantamount to a Bernoulli probability of success).

    Args:
        genes_obj (Genes): An object of the Genes class defining the genes.
        pins (float): The probability of insertion.

    Returns:
        list. List of lists with results and information for the genes. The
            elements of the list are as follows:
            - ORF ID.
            - Gene Name.
            - Gene Description.
            - Number of TA sites with insertions.
            - Number of TA sites.
            - Length of largest run of non-insertion.
            - Expected run for a gene this size.
            - p-value of the observed run.
    """

    pnon = 1.0 - pins
    results = []
    for gene in genes_obj:
        if gene.n == 0:
            results.append([gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, 0.0, 1.000])
        else:
            B = 1.0/math.log(1.0/pnon)
            u = math.log(gene.n*pins, 1.0/pnon)
            exprun = ExpectedRuns(gene.n, pnon)
            pval = 1.0 - GumbelCDF(gene.r, u, B)
            results.append([gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, exprun, pval])
    return(results)





if __name__ == "__main__":

    G = Genes(sys.argv[1].split(","), sys.argv[2], norm="TTR")
    print "#Insertion: %s" % G.global_insertion()
    print "#Sites: %s" % G.global_sites()
    print "#Run: %s" % G.global_run()
    print "#Theta: %1.4f" % G.global_theta()
    print "#Phi: %1.4f" % G.global_phi()
    print "#"

    #g = G["Rv1968"]
    #print g
    #print g.runs
    #print runindex(g.runs)
    
    #sys.exit()

    griffin_results = griffin_analysis(G, G.global_theta())
    for i,gene in enumerate(sorted(G)):
        pos = gene.position
        exprun, pval = griffin_results[i][-2:]
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.1f\t%1.5f" % (gene.orf, gene.name, gene.k, gene.n, gene.r, gene.s, gene.t, exprun, pval)



