import sys,numpy
import scipy.stats

STATES = "ES GD NE GA".split()

# first pass...

headers = []
data = []
all = []
Sats = {}
NZmeans = {}
Means = {}
TAs = {}
Calls = {}
for line in open(sys.argv[1]):
  line = line.strip()
  if line[0]=='#': headers.append(line); continue
  w = line.split('\t')
  nTA = int(w[3])
  if nTA==0: continue
  votes = [int(x) for x in w[4:8]]
  m = max(votes)
  consistency = m/float(nTA)
  id = w[0]
  Sats[id] = float(w[-3])
  NZmeans[id] = float(w[-2])
  Means[id] = Sats[id]*NZmeans[id]
  TAs[id] = nTA
  Calls[id] = w[-1]
  all.append(consistency)
  data.append(w)

Info = ["# HMM confidence info:"]
Info.append("# avg gene-level consistency of HMM states: %s" % (round(numpy.mean(all),4)))

cons = numpy.mean(all)

IDs = [x[0] for x in data]

meanSats,meanNZmeans = {},{}
stdSats,stdNZmeans = {},{}
SatParams = {}
NZmeanParams = {}
MeanParams = {}

Info.append("# state posterior probability distributions:")
for st in ["ES","GD","NE","GA"]:
  sub = list(filter(lambda id: Calls[id]==st,IDs))
  nzmeans = [NZmeans[id] for id in sub]
  sats = [Sats[id] for id in sub]
  means = [Means[id] for id in sub]

  meanSat = numpy.mean(sats)
  meanNZmean = numpy.mean(nzmeans)
  stdSat = numpy.std(sats)
  stdNZmean = numpy.std(nzmeans)
  medNZmeans = numpy.median(nzmeans)
  iqrNZmeans = scipy.stats.iqr(nzmeans)
  meanMeans = -999 if len(means)==0 else numpy.median(means) # -999 if there are no GA genes, for example
  stdMeans = max(1.0,0.7314*scipy.stats.iqr(means)) # don't let stdev collapse to 0 for ES

  # model NZmean with robust Normal distribution: use median and IQR
  sigma = 0.7413*iqrNZmeans
  sigma = max(0.01,sigma) # don't let it collapse to 0
  NZmeanParams[st] = (medNZmeans,sigma)

  # model saturation as Beta distribution, fit by method of moments:
  # https://real-statistics.com/distribution-fitting/method-of-moments/method-of-moments-beta-distribution/
  alpha = meanSat*( (meanSat*(1.0-meanSat)/(stdSat*stdSat) - 1.0) )
  beta = alpha*(1.0-meanSat)/meanSat
  SatParams[st] = (alpha,beta)
  MeanParams[st] = (meanMeans,stdMeans)

  #Info.append("#   Sat[%s]:    Beta(alpha=%s,beta=%s), E[sat]=%s" % (st,round(alpha,2),round(beta,2),round(alpha/(alpha+beta),3)))
  #Info.append("#   NZmean[%s]: Norm(mean=%s,stdev=%s)" % (st,round(NZmeanParams[st][0],2),round(NZmeanParams[st][1],2)))
  Info.append("#   Mean[%s]:   Norm(mean=%s,stdev=%s)" % (st,round(MeanParams[st][0],2),round(MeanParams[st][1],2)))


def normalize(L):
  tot = sum(L)
  return [x/float(tot) for x in L]

# prob of each state is joint prob over saturation (Beta) and NZmean (Normal)
# for all 4 states; uses globals

def calc_probs3(sats,nzmean):
  A,B = [],[]
  for st in STATES:
    alpha,beta = SatParams[st]
    meanNZmean,stdNZmean = NZmeanParams[st]
    A.append(scipy.stats.beta.pdf(sat,a=alpha,b=beta))
    B.append(scipy.stats.norm.pdf(nzmean,loc=meanNZmean,scale=stdNZmean))
  A = normalize(A)
  B = normalize(B)
  probs = [a*b for a,b in zip(A,B)]
  return normalize(probs)

# prob of each state is based Gaussian density for Mean count for each gene (combines Sat and NZmean)

def calc_probs4(sats,nzmean):
  probs = []
  for st in STATES:
    meanMeans,stdMeans = MeanParams[st]
    probs.append(0 if meanMeans<0 else scipy.stats.norm.pdf(sat*nzmean,loc=meanMeans,scale=stdMeans))
  return normalize(probs)

##################################

# second pass...

num_low_conf,num_ambig = 0,0
results = []
for line in open(sys.argv[1]):
  line = line.strip()
  if line[0]=='#': continue
  w = line.split('\t')
  id,Call = w[0],w[-1]
  nTA = int(w[3])
  if nTA==0: continue; # print (line); continue
  votes = [int(x) for x in w[4:8]]
  m = max(votes)
  consistency = m/float(nTA)
  sat,NZmean = float(w[-3]),float(w[-2])
  PC = 0.01 # shrink range of saturation form [0,1] to [0.01,0.99] to prevent nan's from Beta
  sat = max(PC,min(1.0-PC,sat))
  probs = calc_probs4(sat,NZmean) # normalized
  conf = probs[STATES.index(Call)]

  flag=""

  #if conf<0.5: flag="low-confidence"
  #if max(probs)<0.7:
  #  if (Call=="ES" or Call=="GD") and (probs[0]>0.25 and probs[1]>0.25): flag="ambiguous"
  #  if (Call=="GD" or Call=="NE") and  (probs[1]>0.25 and probs[2]>0.25): flag="ambiguous"
  #  if (Call=="NE" or Call=="GA") and (probs[2]>0.25 and probs[3]>0.25): flag="ambiguous"

  if conf<0.2: flag = "low-confidence"
  if conf>=0.2 and conf!=max(probs): flag = "ambiguous"

  if flag=="ambiguous": num_ambig += 1
  if flag=="low-confidence": num_low_conf += 1
  
  vals = w+[round(sat*NZmean,1),round(consistency,3)]+[round(x,6) for x in probs]
  vals += [round(conf,4),flag]
  results.append(vals)

Info.append("# num low-confidence genes=%s, num ambiguous genes=%s" % (num_low_conf,num_ambig))

###############
# print results

for line in headers[:-1]: print(line)
for line in Info: print(line)
newheaders = headers[-1].split('\t') # assume last comment line is column headers
newheaders += "Mean consis probES probGD probNE probGA conf flag".split()
print('\t'.join(newheaders))

for vals in results:
  print ('\t'.join([str(x) for x in vals]))

