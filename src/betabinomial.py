import sys
import time
import datetime
#import trash_tools
import numpy
import scipy.stats
import transit_tools





binomial_prefix = "[Binomial]"


def get_orf_data_transit(orf2reads, orf2pos, orf2info, REALREPS):
    G = len(orf2reads)
    K = numpy.zeros(G,dtype=int); N = numpy.zeros(G,dtype=int); ORF = [];

    g = 0
    for orf in sorted(orf2reads):
        if REALREPS:
            fullreads = numpy.array(orf2reads[orf]).flatten()
            K[g]=numpy.sum(fullreads > 0)
            N[g]=numpy.sum(fullreads >= 0)
            ORF.append(orf)
        else:
            fullreads = numpy.array(orf2reads[orf])
            #print fullreads
            #print numpy.sum(numpy.sum(fullreads,1) > 0)
            #print len(numpy.sum(fullreads,1))
            if fullreads != []:
                K[g] = numpy.sum(numpy.sum(fullreads,1) > 0)
                N[g] = len(numpy.sum(fullreads,1))
            ORF.append(orf)
        g+=1
    return((ORF,K,N))



def runBinomial(wx, pubmsg, **kwargs):

    try:
        from wx.lib.pubsub import pub
        pub.subscribe
        newWx = True
    except AttributeError as e:
        from wx.lib.pubsub import Publisher as pub
        newWx = False


    print binomial_prefix, "Running Binomial Method"


    wigList = kwargs.get("readPathList")
    annotationPath = kwargs.get("annotationPath")
    desiredSamples = kwargs.get("samples", 10000)
    BURNIN = kwargs.get("burnin", 500)
    IGNORECODON = kwargs.get("ignoreCodon", True)
    ignoreNTerm = kwargs.get("ignoreNTerm", 0)
    ignoreCTerm = kwargs.get("ignoreCTerm", 0)
    output = kwargs.get("output")
    
    
    w1 = 0.15
    w0 = 1.0 - w1
    ALPHA = 1
    BETA = 1
    ALPHA_w = 600
    BETA_w = 3400
    mu_c = 0

    acctot = 0.0

    LANES = [1]
    #MINIMUM_READ = min_read
    
    SAMPLE_SIZE = BURNIN+desiredSamples
    
    MID = False
    REALREPS = True
    numReps = len(wigList)
    VERBOSE = False

    orf2info = transit_tools.get_gene_info(annotationPath)
    hash = transit_tools.get_pos_hash(annotationPath)
    (data, position) = transit_tools.get_data(wigList)
    orf2reads, orf2pos = transit_tools.get_gene_reads(hash, data, position, orf2info, ignoreCodon=IGNORECODON, ignoreNTerm=ignoreNTerm, ignoreCTerm=ignoreCTerm, orf_list=orf2info.keys())

    (ORF, K, N) = get_orf_data_transit(orf2reads, orf2pos, orf2info, REALREPS)


    G = len(K)
    theta = numpy.zeros((G, SAMPLE_SIZE))
    theta[:,0] = 0.10

    tc = numpy.zeros(G); tp = numpy.zeros(G)

    rho0 = numpy.zeros(SAMPLE_SIZE); rho0[0] = 0.5;  Kp0 = numpy.zeros(SAMPLE_SIZE); Kp0[0] = 10;
    rho1 = numpy.zeros(SAMPLE_SIZE); rho1[0] = 0.10; Kp1 = numpy.zeros(SAMPLE_SIZE); Kp1[0] = 3;

    pi0 = 0.5; M0 = 1;
    pi1 = 0.5; M1 = 1;

    a0 = 10; b0 = 1;
    a1 = 10; b1 = 1;

    Z = numpy.zeros((G, SAMPLE_SIZE))
    pz1 = numpy.zeros(SAMPLE_SIZE);
    n1 = 0
    alpha_w = 0.5; beta_w = 0.5;

    #w1 = 0.15
    w1 = scipy.stats.beta.rvs(alpha_w, beta_w)
    W1 = numpy.zeros(SAMPLE_SIZE); W1[0] = w1



    for g in range(G):
        if N[g] == 0: theta[g][0] = 0.5
        elif K[g]/float(N[g]) == 0: theta[g][0] = 0.001
        elif K[g]/float(N[g]) == 1: theta[g][0] = 0.001
        else: theta[g][0] = K[g]/float(N[g])

        #print g, ORF[g], K[g], N[g], theta[g][0]
        Z[g][0] = scipy.stats.bernoulli.rvs(1-theta[g][0])



    acc_p0 = 0; acc_k0 = 0;
    acc_p1 = 0; acc_k1 = 0;


    rho0c_std = 0.010
    kp0c_std = 1.40
    rho1c_std = 0.009
    kp1c_std = 1.1


    for i in xrange(1, SAMPLE_SIZE):
        
        i0 = Z[:,i-1] == 0; n0 = numpy.sum(i0);
        i1 = Z[:,i-1] == 1; n1 = numpy.sum(i1);
        
        theta[i0,i] = scipy.stats.beta.rvs(Kp0[i-1]*rho0[i-1] + K[i0],  Kp0[i-1]*(1-rho0[i-1]) + N[i0] - K[i0])
        theta[i1,i] = scipy.stats.beta.rvs(Kp1[i-1]*rho1[i-1] + K[i1],  Kp1[i-1]*(1-rho1[i-1]) + N[i1] - K[i1])
        
        if VERBOSE: print >> sys.stderr, "i=%d\ttheta[0,i]=%f\tn0=%d\tn1=%d\tw1=%f" %(i, theta[0,i], n0, n1, w1)
        
        #rho0_c = rho0[i-1] + scipy.stats.norm.rvs(0,0.010)
        #Kp0_c = Kp0[i-1] + scipy.stats.norm.rvs(0,1.68)
        rho0_c = rho0[i-1] + scipy.stats.norm.rvs(0, rho0c_std)
        Kp0_c = Kp0[i-1] + scipy.stats.norm.rvs(0, kp0c_std)
        
        
        if rho0_c <= 0: rho0[i] = rho0[i-1]
        else:
            fc = numpy.log(scipy.stats.beta.pdf(rho0_c, M0*pi0, M0*(1-pi0)))
            f0 = numpy.log(scipy.stats.beta.pdf(rho0[i-1], M0*pi0, M0*(1-pi0)))
            fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0_c, Kp0[i-1]*(1-rho0_c))))
            f0 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0[i-1], Kp0[i-1]*(1-rho0[i-1]))))

            if numpy.log(scipy.stats.uniform.rvs()) < fc - f0:
                rho0[i] = rho0_c
                acc_p0+=1
            else: rho0[i] = rho0[i-1]



        if Kp0_c <= 0: Kp0[i] = Kp0[i-1]
        else:
            fc = numpy.log(scipy.stats.gamma.pdf(Kp0_c, a0, b0));
            f0 = numpy.log(scipy.stats.gamma.pdf(Kp0[i-1], a0, b0));
            fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0_c*rho0[i], Kp0_c*(1-rho0[i]))))
            f0 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i0,i], Kp0[i-1]*rho0[i], Kp0[i-1]*(1-rho0[i]))))

            if numpy.log(scipy.stats.uniform.rvs()) < fc - f0:
                Kp0[i] = Kp0_c
                acc_k0+=1
            else: Kp0[i] = Kp0[i-1]

        if VERBOSE: print >> sys.stderr, "rho0_c=%f\trho0=%f\tKp0_c=%f\tKp0=%f" %(rho0_c, rho0[i], Kp0_c, Kp0[i])


        rho1_c = rho1[i-1] + scipy.stats.norm.rvs(0, rho1c_std)
        Kp1_c = Kp1[i-1] + scipy.stats.norm.rvs(0, kp1c_std)


        if rho1_c <= 0:
            rho1[i] = rho1[i-1]
        else:
            fc = numpy.log(scipy.stats.beta.pdf(rho1_c, M1*pi1, M1*(1-pi1)))
            f1 = numpy.log(scipy.stats.beta.pdf(rho1[i-1], M1*pi1, M1*(1-pi1)))
            fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1_c, Kp1[i-1]*(1-rho1_c))))
            f1 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1[i-1], Kp1[i-1]*(1-rho1[i-1]))))

            if numpy.log(scipy.stats.uniform.rvs()) < fc - f1:
                rho1[i] = rho1_c
                acc_p1+=1
            else: rho1[i] = rho1[i-1]


        if Kp1_c <= 0: Kp1[i] = Kp1[i-1]
        else:
            fc = numpy.log(scipy.stats.gamma.pdf(Kp1_c, a1, b1));
            f1 = numpy.log(scipy.stats.gamma.pdf(Kp1[i-1], a1, b1));
            fc += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1_c*rho1[i], Kp1_c*(1-rho1[i]))))
            f1 += numpy.sum(numpy.log(scipy.stats.beta.pdf(theta[i1,i], Kp1[i-1]*rho1[i], Kp1[i-1]*(1-rho1[i]))))

            if numpy.log(scipy.stats.uniform.rvs()) < fc - f1:
                Kp1[i] = Kp1_c
                acc_k1+=1
            else: Kp1[i] = Kp1[i-1]



        if VERBOSE: print >> sys.stderr, "rho1_c=%f\trho1=%f\tKp1_c=%f\tKp1=%f" %(rho1_c, rho1[i], Kp1_c, Kp1[i])

        g0 = scipy.stats.beta.pdf(theta[:,i], Kp0[i]*rho0[i], Kp0[i]*(1-rho0[i])) * (1-w1)
        g1 = scipy.stats.beta.pdf(theta[:,i], Kp1[i]*rho1[i], Kp1[i]*(1-rho1[i])) * (w1)
        p1 = g1/(g0+g1)
        p1 = numpy.nan_to_num(p1)
        if VERBOSE: print >> sys.stderr, "nan=%d\tinf=%d" % ( numpy.sum(numpy.isnan(p1)), numpy.sum(numpy.isinf(p1)) )
        if VERBOSE: print >> sys.stderr, "g0[0]=%f\tg1[0]=%f\tp1[0]=%f" % (g0[0],g1[0],  p1[0])
        try:
            Z[:,i] = scipy.stats.bernoulli.rvs(p1)
        except:
            inan = numpy.isnan(p1)
            print >> sys.stderr, "K=\t", K[inan]
            print >> sys.stderr, "N=\t", N[inan]
            print >> sys.stderr, "theta=", theta[inan,i]
            sys.exit()
        pz1[i] = p1[0]

        i1 = Z[:,i] == 1; n1 = numpy.sum(i1);
        #w1 = 0.15
        w1 = scipy.stats.beta.rvs(alpha_w + n1, beta_w + G - n1)
        W1[i] = w1

        
        #Update 
        if wx:
            if newWx:
                wx.CallAfter(pubmsg, "binomial", msg="Running Binomial Method... %2.0f%%" % (100.0*(i+1)/(SAMPLE_SIZE)))
            else:
                wx.CallAfter(pubmsg, "binomial", "Running Binomial Method... %2.0f%%" % (100.0*(i+1)/(SAMPLE_SIZE)))



    z_bar = numpy.apply_along_axis(numpy.mean, 1, Z[:, BURNIN:])
    theta_bar = numpy.apply_along_axis(numpy.mean, 1, theta[:, BURNIN:])
    (ess_threshold, noness_threshold) = transit_tools.fdr_post_prob(z_bar)

    output.write("#Binomial\n")
    output.write("#Command: %s\n" % " ".join(["%s=%s" %(key,val) for (key,val) in kwargs.items()]))
    output.write("#Thresholds: (%1.5f, %1.5f)\n" % (ess_threshold,noness_threshold))
    output.write("#rho0 Acceptance Rate:\t%f%%\n" % ((100.0*acc_p0)/SAMPLE_SIZE))
    output.write("#Kp0  Acceptance Rate:\t%f%%\n" % ((100.0*acc_k0)/SAMPLE_SIZE))
    output.write("#rho1 Acceptance Rate:\t%f%%\n" % ((100.0*acc_p1)/SAMPLE_SIZE))
    output.write("#Kp1  Acceptance Rate:\t%f%%\n" % ((100.0*acc_k1)/SAMPLE_SIZE))


    output.write("#Orf\tName\tDescription\tMean Insertion\tSites per Replicate\tTotal Insertions\tTotal Sites\tthetabar\tzbar\tCall\n")

    for g in xrange(G):

        c = "Uncertain"
        if z_bar[g] > ess_threshold:
            c = "Essential"
        if z_bar[g] < noness_threshold:
            c = "Non-Essential"

        output.write("%s\t%s\t%s\t%1.1f\t%d\t%d\t%d\t%f\t%f\t%s\n" % (ORF[g], orf2info[ORF[g]][0], orf2info[ORF[g]][1], K[g]/float(numReps), N[g]/numReps, K[g], N[g], theta_bar[g], z_bar[g], c))


    output.close()

    if not output.name.startswith("<"):
        data = {"path":output.name, "type":"Binomial", "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p")}
        print binomial_prefix, "Adding File:", output.name
        if wx:
            if newWx:
                wx.CallAfter(pubmsg, "file", data=data)
            else:
                wx.CallAfter(pubmsg, "file", data)

    if wx:
        if newWx:
            wx.CallAfter(pubmsg, "binomial", msg="Finished!")
            wx.CallAfter(pubmsg,"finish", msg="binomial")
        else:
            wx.CallAfter(pubmsg, "binomial", "Finished!")
            wx.CallAfter(pubmsg,"finish", "binomial")

    print binomial_prefix, "Finished Binomial Method"




