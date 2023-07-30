import sys

FASTQ_FILE = sys.argv[1]
IDS_FILE = sys.argv[2]

###################################

complement = {'A':'T','T':'A','C':'G','G':'C'}

def reverse_complement(seq):
  s = list(seq)
  s.reverse()
  for i in range(len(s)):
    s[i] = complement.get(s[i],s[i]) # if unknown, leave as it, e.g > or !
  s = ''.join(s)
  return s

###################################

# example: RVBD0015c:pknA_Essential_GACGCGCGCGAATGCGGTGTCG_22mer_v2PAMscore22 

IDs = []
barcodemap = {} # hash from barcode to full ids

for line in open(IDS_FILE):
  w = line.rstrip().split('\t')
  id = w[0]
  v = id.split('_')
  if len(v)<3: continue
  barcode = v[2]
  IDs.append(id)
  # reverse-complement of barcodes appears in reads, so hash them that way
  barcodemap[reverse_complement(barcode)] = id

###################################

# example: AGCTTCTTTCGAGTACAAAAAC xxxx TCCCAGATTATATCTATCACTGA , where xxxx is rev. compl. of a barcode

counts = {}

A,B = "AGCTTCTTTCGAGTACAAAAAC","TCCCAGATTATATCTATCACTGA"
#A,B = "GTACAAAAAC","TCCCAGATTA"
lenA = len(A)
cnt,nreads,recognized = 0,0,0
for line in open(FASTQ_FILE):
  cnt += 1
  if cnt%4==2:
    nreads += 1
    if (nreads%1000000==0): sys.stderr.write("reads=%s, recognized barcodes=%s (%0.1f%%)\n" % (nreads,recognized,100.*recognized/float(nreads)))
    seq = line.rstrip()
    a = seq.find(A)
    if a==-1: continue
    b = seq.find(B)
    if b==-1: continue
    sz = b-(a+lenA)
    if sz<10 or sz>30: continue
    barcode = seq[a+lenA:b] # these are reverse-complements, but rc(barcodes) stored in hash too
    if barcode not in barcodemap: continue
    id = barcodemap[barcode]
    if id not in counts: counts[id] = 0
    counts[id] += 1
    recognized += 1

for id in IDs:
  vals = [id,counts.get(id,0)]
  print('\t'.join([str(x) for x in vals]))
