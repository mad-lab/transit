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
import datetime
import random
import numpy
import matplotlib.pyplot as plt
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



resampling_prefix = "[Resampling]"

#def runResampling(ctrlString, expString, annotationPath, sampleSize, histPath, doAdaptive, ignoreCodon, ignoreNTerm, ignoreCTerm, output, wx, pubmsg, doNormalize=True):
def runResampling(wx, pubmsg, **kwargs):
    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False


    print resampling_prefix, "Running Resampling Method"

    arguments = locals().items()


    ctrlList = kwargs.get("ctrlList")
    expList = kwargs.get("expList")
    annotationPath = kwargs.get("annotationPath")
    sampleSize = kwargs.get("sampleSize", 10000)
    histPath = kwargs.get("histPath")
    doAdaptive = kwargs.get("doAdaptive", False)
    normalize = kwargs.get("normalize", "nzmean")
    doLOESS = kwargs.get("doLOESS", False)
    ignoreCodon = kwargs.get("ignoreCodon", True)
    ignoreNTerm = kwargs.get("ignoreNTerm", 0)
    ignoreCTerm = kwargs.get("ignoreCTerm", 0)
    output = kwargs.get("output")


    N1 = len(ctrlList)
    N2 = len(expList)

    orf2info = transit_tools.get_gene_info(annotationPath)

    hash = transit_tools.get_pos_hash(annotationPath)
    (data, position) = transit_tools.get_data(ctrlList + expList)

    print resampling_prefix, "Normalizing with", normalize
    if normalize == "nzmean":
        factors = transit_tools.nzmean_factors(data)
        data = factors * data
    if normalize == "totreads":
        factors = transit_tools.totreads_factors(data)
        data = factors * data
    elif normalize == "zinfnb":
        factors = transit_tools.zinfnb_factors(data)
        data = factors * data
    elif normalize == "quantile":
        data = transit_tools.quantile_norm(data)
    else:
        pass

    if doLOESS:
        for j in range(len(data)):
            data[j] = transit_tools.loess_correction(position, data[j])
        

    orf2reads,orf2pos = transit_tools.get_gene_reads(hash, data, position, orf2info, ignoreCodon=ignoreCodon, ignoreNTerm=ignoreNTerm, ignoreCTerm=ignoreCTerm, orf_list=orf2info.keys())
    
    S = sampleSize

    output.write("#Resampling\n")
    output.write("#Command: python transit.py %s\n" % " ".join(["%s=%s" %(key,val) for (key,val) in kwargs.items()]))
    output.write("#Control Samples:  %s\n" % ", ".join(ctrlList))
    output.write("#Experimental Samples:  %s\n" % ", ".join(expList))
    if  normalize in ["nzmean", "totreads", "zinfnb"]:
        output.write("#%s factors: %s\n" % (normalize ,"\t".join(["%1.4f" % f for f in factors])))
    elif normalize == "quantile":
        output.write("#quantile factors: no factors used in quantile normalization\n")
    else:
        output.write("#Not normalized\n")

    output.write("#Orf\t%Name\tDescription\tN\tTAs Hit\tSum Rd 1\tSum Rd 2\tDelta Rd\tlog2 FC\tp-value\tp-adj\n")
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
            sum1 = 0; sum2 = 0;

        try:
            log2FC = math.log(float(sum2)/float(sum1),2)
        except:
            log2FC = 0


        count_utail = 0
        count_ltail = 0
        count_2tail = 0
   
        #delta_sum_list = numpy.zeros(S)
        delta_sum_list = []
        s_performed = 0
        for s in range(S):
            #Reads
            #if len(reads) == 0: breaki
            if len(reads) >0:
                perm = numpy.random.permutation(reads.flatten()).reshape(reads.shape)
                sumA,sumB = numpy.sum(perm[:,:N1]), numpy.sum(perm[:,N1:])
            else:
                 sumA = 0; sumB = 0;
            #all_perm=np.array((list(itertools.permutations([0,1,2,3]))))
            #(reads.flatten()[(all_perm[numpy.random.randint(0,len(all_perm),size=100)]+3*np.arange(100)[...,numpy.newaxis]).flatten()]).reshape(reads.shape)

            delta_sum_list.append(sumB-sumA)
            if sumB-sumA >= sum2-sum1: count_utail+=1
            if sumB-sumA <= sum2-sum1: count_ltail+=1
            if abs(sumB-sumA) >= abs(sum2-sum1): count_2tail+=1


            s_performed+=1
            if doAdaptive:
                if s_performed == round(S*0.01) or s_performed == round(S*0.1) or s_performed == round(S*1):
                    if count_2tail >= round(S*0.01*0.10):
                        break
            


        #Reads - pvalues
        pval_utail = count_utail/float(s_performed)
        pval_ltail = count_ltail/float(s_performed)
        pval_2tail = count_2tail/float(s_performed)




        orf2out[orf] = (orf, orf2info[orf][0], orf2info[orf][1], len(fullreads), len(reads), sum1, sum2, sum2-sum1, log2FC, pval_2tail)
        pval.append(pval_2tail)


        count += 1
        # the histogram of the data
        if wx and newWx: wx.CallAfter(pubmsg, "histogram", msg=(delta_sum_list, orf, histPath, sum2-sum1))
        if wx and not newWx: wx.CallAfter(pubmsg, "histogram", (delta_sum_list, orf, histPath, sum2-sum1))

        # Update Progress
        if wx and newWx: wx.CallAfter(pubmsg, "resampling", msg="Running Resampling Method... %2.0f%%" % (100.0*(count+1)/(G)))
        if wx and not newWx: wx.CallAfter(pubmsg, "resampling", "Running Resampling Method... %2.0f%%" % (100.0*(count+1)/(G)))

    qval = fdr_corrected_pval(pval)
    count = 0
    for orf in sorted(orf2out):
        output.write("%s\t%s\t%s\t%d\t%d\t%1.1f\t%1.1f\t%1.1f\t%1.2f\t%1.5f" % orf2out[orf])
        output.write("\t%1.5f\n" % qval[count])
        count += 1

    output.close()



    print resampling_prefix, "Finished Resampling Method"
    if not output.name.startswith("<"):
        data = {"path":output.name, "type":"Resampling", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print resampling_prefix, "Adding File:", output.name
        if wx and newWx: wx.CallAfter(pubmsg, "file", data=data)
        if wx and not newWx: wx.CallAfter(pubmsg, "file", data)

    if wx:
        if newWx:
            wx.CallAfter(pubmsg, "resampling", msg="Finished!")
            wx.CallAfter(pubmsg,"finish", msg="resampling")
        else:
            wx.CallAfter(pubmsg, "resampling", "Finished!")
            wx.CallAfter(pubmsg,"finish", "resampling")








