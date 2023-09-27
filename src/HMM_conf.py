import sys,numpy
import scipy.stats

data = []
all = []
Sats = {}
NZmeans = {}
TAs = {}
Calls = {}
for line in open(sys.argv[1]):
  line = line.strip()
  if line[0]=='#': continue
  w = line.split('\t')
  nTA = int(w[3])
  if nTA==0: continue
  votes = [int(x) for x in w[4:8]]
  m = max(votes)
  consistency = m/float(nTA)
  id = w[0]
  Sats[id] = float(w[-3])
  NZmeans[id] = float(w[-2])
  TAs[id] = nTA
  Calls[id] = w[-1]
  all.append(consistency)
  data.append(w)

print("# avg gene-level consistency of HMM states: %s" % (round(numpy.mean(all),4)))

cons = numpy.mean(all)

IDs = [x[0] for x in data]

meanSats,meanNZmeans = {},{}
stdSats,stdNZmeans = {},{}

print("# state posterior probability distributions:")
for st in ["ES","GD","NE","GA"]:
  sub = list(filter(lambda id: Calls[id]==st,IDs))
  meanSat = numpy.mean([Sats[id] for id in sub])
  meanNZmean = numpy.mean([NZmeans[id] for id in sub])
  stdSat = numpy.std([Sats[id] for id in sub])
  stdNZmean = numpy.std([NZmeans[id] for id in sub])
  print("#  %s: genes=%s, meanSat=%s, stdSat=%s, meanNZmean=%s, stdNZmean=%s" % (st,len(sub),round(meanSat,3),round(stdSat,3),round(meanNZmean,1),round(stdNZmean,1)))
  meanSats[st] = meanSat
  meanNZmeans[st] = meanNZmean
  stdSats[st] = meanSat
  stdNZmeans[st] = meanNZmean

def calc_prob(sat,NZmean,meanSat,stdSat,meanNZmean,stdNZmean):
  a = scipy.stats.norm.pdf(sat,loc=meanSat,scale=stdSat)
  b = scipy.stats.norm.pdf(NZmean,loc=meanNZmean,scale=stdNZmean)
  return a*b

# second pass...

for line in open(sys.argv[1]):
  line = line.strip()
  if line[0]=='#': print (line); continue
  w = line.split('\t')
  id,Call = w[0],w[-1]
  nTA = int(w[3])
  if nTA==0: print (line,"<skip>"); continue
  votes = [int(x) for x in w[4:8]]
  m = max(votes)
  consistency = m/float(nTA)
  #prob = scipy.stats.binom.cdf(m,nTA,cons)
  probs = []
  STATES = "ES GD NE GA".split()
  for st in STATES:
    probs.append(calc_prob(float(w[-3]),float(w[-2]),meanSats[st],stdSats[st],meanNZmeans[st],stdNZmeans[st]))
  totprob = sum(probs)
  relprobs = [x/float(totprob) for x in probs]
  flag=0
  if max(relprobs)<0.7:
    if (Call=="ES" or Call=="GD") and (relprobs[0]>0.25 and relprobs[1]>0.25): flag=1
    if (Call=="GD" or Call=="NE") and  (relprobs[1]>0.25 and relprobs[2]>0.25) :flag=1
    if (Call=="NE" or Call=="GA") and (relprobs[2]>0.25 and relprobs[3]>0.25):flag=1
  conf = relprobs[STATES.index(Call)]
  #vals = w+[nTA,m,round(consistency,3),round(prob,4)]
  vals = w+[round(consistency,3)]+[round(x,6) for x in relprobs]
  vals += [round(conf,4),flag]
  print ('\t'.join([str(x) for x in vals]))

