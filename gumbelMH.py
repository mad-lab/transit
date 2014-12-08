import sys
import time
#import trash_tools
from MH_tools import *

def runGumbel(readpath, annotationPath, min_read, samples):


    BURNIN = 1
    TRIM = 1
    w1 = 0.15
    w0 = 1.0 - w1
    ALPHA = 1
    BETA = 1
    ALPHA_w = 600
    BETA_w = 3400
    mu_c = 0
    
    acctot = 0.0
    
    LANES = [1]
    MINIMUM_READ = min_read
    SAMPLE_SIZE = samples
    MID = False
    PATH = readpath
    PROT_PATH = annotationPath
    IGV = False; WIG = False; TRASH = False;
    
    if readpath.endswith(".igv"):
        IGV = True
    elif readpath.endswith(".wig"):
        WIG = True
    elif readpath.endswith(".txt"):
        WIG = True

    if IGV:
        orf_to_reads = read_IGV_file(PATH)
    elif WIG:
        orf_to_reads = read_WIG_file(PATH, PROT_PATH)
    else:
        orf_to_reads = read_TRASH_file(PATH, LANES, MINIMUM_READ)

    start_time = time.time()
    (ORF_all, K_all, N_all, R_all, S_all, T_all) = get_orf_data(orf_to_reads, MINIMUM_READ, mid=MID, prot=PROT_PATH)
    bad_orf_set = set([ORF_all[g] for g in xrange(len(N_all)) if not good_orf(N_all[g], T_all[g])]);
    bad_orf_set.add("Rvnr01");
    (ORF, K, N, R, S, T) = get_orf_data(orf_to_reads, MINIMUM_READ, mid=MID, prot=PROT_PATH, bad=bad_orf_set)

    orf2name = {}; orf2desc = {}
    for line in open(PROT_PATH):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]; name = tmp[7]; desc = tmp[0];
        orf2name[orf] = name; orf2desc[orf] = desc

    #rv = "Rv0001"
    #print orf_to_reads[rv]
    #print K
    #print K[0]
    #print N[0]
    #sys.exit()


    mu_s, temp, sigma_s = regress(R,S) # Linear regression to estimate mu_s, sigma_s for span data
    mu_r, temp, sigma_r = regress(S, R) # Linear regression to estimate mu_r, sigma_r for run data


    N_GENES = len(N)
    Z_sample = numpy.zeros((N_GENES, SAMPLE_SIZE))
    Z = [classify(N[g], R[g], 0.5)   for g in xrange(N_GENES)]
    Z_sample[:,0] = Z
    N_ESS = numpy.sum(Z_sample[:,0] == 1)

    phi_sample = numpy.zeros(SAMPLE_SIZE) #[]
    phi_sample[0] = PHI_START
    phi_old = PHI_START
    phi_new = 0.00

    SIG = numpy.zeros(len(S))
    for g in range(len(S)):
        SIG[g] = sigmoid(S[g], T[g]) * scipy.stats.norm.pdf(R[g], mu_r*S[g], sigma_r) 



    i = 1; count = 0;
    while i < SAMPLE_SIZE:

        # PHI
        acc = 1.0
        phi_new  = phi_old + random.gauss(mu_c, sigma_c)
        i0 = Z_sample[:,i-1] == 0
        if phi_new > 1 or phi_new <= 0 or (F_non(phi_new, N[i0], R[i0]) - F_non(phi_old, N[i0], R[i0])) < math.log(random.uniform(0,1)):
            phi_new = phi_old
            acc = 0.0
            flag = 0

        # Z
        Z = sample_Z(phi_new, w1, N, R, S, T, mu_s, sigma_s, SIG)

        # w1
        N_ESS = sum(Z == 1)
        w1 = scipy.stats.beta.rvs(N_ESS + ALPHA_w, N_GENES - N_ESS + BETA_w)
    
        count +=1
        acctot+=acc

        if (count > BURNIN) and (count % TRIM == 0):
            phi_sample[i] = phi_new
            Z_sample[:,i] = Z
            i+=1

    
        phi_old = phi_new


    ZBAR = numpy.apply_along_axis(numpy.mean, 1, Z_sample)
    
    (ess_t, non_t) = fdr_post_prob(ZBAR)

    #Orf	k	n	r	s	zbar
    print "#FDR Corrected thresholds:", ess_t, non_t
    print "#Command used:\tpython %s" % (" ".join(sys.argv))
    print "#MH Acceptance-Rate:\t%2.2f%%" % (100.0*acctot/count)
    print "#Total Iterations Performed:\t%d" % count
    print "#Sample Size:\t%d" % i
    print "#phi estimate:\t%f" % numpy.average(phi_sample)
    print "#Time:", time.time() - start_time
    if VERBOSE: print "#Orf\tName\tDesc\tk\tn\tr\ts\tzbar\tCall\tSample"
    else: print "#Orf\tName\tDesc\tk\tn\tr\ts\tzbar\tCall"
    i = -1
    for g in xrange(len(ORF_all)):
        k = K_all[g]; n = N_all[g];  r = R_all[g]; s = S_all[g]; orf=ORF_all[g];
        if orf not in bad_orf_set: 
            i+=1; zbar = ZBAR[i]; sample_str = "\t"+ ",".join(["%d" % x for x in Z_sample[i,:]]);
        else:
            zbar = -1.0; sample_str = "\t" + "-1";
        if not VERBOSE: sample_str = ""

        if zbar > ess_t: call = "E"
        elif non_t <= zbar <= ess_t: call = "U"
        elif 0 <= zbar < non_t: call = "NE"
        else: call = "S"

        print "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%s%s" % (orf, orf2name.get(orf,"-"), orf2desc.get(orf,"-"), k, n, r, s, zbar, call, sample_str )
            
