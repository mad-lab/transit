#!/usr/bin/env python

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


import glob,os,sys,time,math
import sys, re, shutil
import platform
import gzip
import subprocess
from collections import defaultdict

def cleanargs(rawargs):
    #TODO: Write docstring
    args = []
    kwargs = {}
    count = 0
    while count < len(rawargs):
        if rawargs[count].startswith("-"): #and len(rawargs[count].split(" ")) == 1:
            if count + 1 < len(rawargs) and (not rawargs[count+1].startswith("-") or len(rawargs[count+1].split(" ")) > 1):
                kwargs[rawargs[count][1:]] = rawargs[count+1]
                count += 1
            else:
                kwargs[rawargs[count][1:]] = True
        else:
            args.append(rawargs[count])
        count += 1
    return (args, kwargs)



def analyze_dataset(wigfile):
  data = []
  TAs,ins,reads = 0,0,0
  for line in open(wigfile):
    if line[0]=='#': continue
    if line[:3]=='var': continue # variableStep
    w = line.rstrip().split()
    TAs += 1
    cnt = int(w[1])
    if cnt>1: ins += 1
    reads += cnt
    data.append((cnt,w[0]))

  output = open(wigfile+".stats","w")
  output.write("total TAs: %d, insertions: %d (%0.1f%%), total reads: %d\n" % (TAs,ins,100*(ins/float(TAs)),reads))
  output.write("mean read count per non-zero site: %0.1f\n" % (reads/float(ins)))
  output.write("5 highest counts:\n")
  data.sort(reverse=True)
  for cnt,coord in data[:5]:
    output.write("coord=%s, count=%s\n" % (coord,cnt))
  output.close()

#############################################################################

def fastq2reads(infile,outfile,maxreads):
  output = open(outfile,"w")
  cnt,tot = 0,0
  for line in open(infile):
    if cnt==0 and line[0]=='@':
        tot += 1
        if tot%1000000==0: message("%s reads processed" % tot)
    if maxreads > -1:
        if tot > maxreads:
            break
    if cnt==0:
      h = line[1:] # strip off '@'
      #h = h.replace(' ','_')
      output.write(">%s" % h)
    if cnt==1: output.write(line)
    cnt = (cnt+1)%4
  output.close()

# the headers for each pair must be identical up to /1 and /2 at the ends
# if the variable character with the read number occurs in the middle, move it to the end

def fix_paired_headers_for_bwa(reads1,reads2):
  a = open(reads1)
  b = open(reads2)
  temp1 = reads1+".temp"
  temp2 = reads2+".temp"
  c = open(temp1,"w")
  d = open(temp2,"w")
  tot = 0
  try:
   while True:
    e = a.readline().rstrip()
    f = b.readline().rstrip()
    if len(e)<=2 or len(f)<=2: break
    if e[0]=='>':
      tot += 1
      if tot%1000000==0: message("%s reads processed" % tot)
      # find first position where there is a difference
      i,n = 0,len(e)
      if len(f)!=n: raise Exception('Error: unexpected format of headers in .fastq files')
      while i<n and e[i]==f[i]: i += 1
      if i==n: raise Exception('Error: unexpected format of headers in .fastq files')
      e = e.replace(' ','_')
      f = f.replace(' ','_')
      e = e.replace('/','_') # this was neceesary for bwa 0.7.10 but not 0.7.12
      f = f.replace('/','_')
      #if i<n-1:
      if e[i+1:]!=f[i+1:]: raise Exception('Error: unexpected format of headers in .fastq files')
      e,f = e[:-1],f[:-1] # strip EOL
      # needed for bwa 0.7.12? which apparently trims off last 2 chars to make ids identical
      # "/[1|2]" are automatically trimmed, but not /3
      e = e[:i]+e[i+1:]
      f = f[:i]+f[i+1:]
    c.write(e+"\n")
    d.write(f+"\n")
  except Exception as ex:
    a.close(); b.close(); c.close(); d.close()
    error(ex.args[0])
  a.close(); b.close(); c.close(); d.close()

  if(os.path.exists(reads1)): os.remove(reads1)
  if(os.path.exists(reads2)): os.remove(reads2)
  os.rename(temp1, reads1)
  os.rename(temp2, reads2)

  '''
  if platform.platform == 'win32':
	  os.system("move %s %s" % (temp1,reads1))
	  os.system("move %s %s" % (temp2,reads2))
  else:
	  os.system("mv %s %s"%(temp1, reads1))
	  os.system("mv %s %s" % (temp2, reads2))
  '''


# checks for a match allowing 1 or 2 mismatches
# if not a match returns -1,-1. If match occurs, returns 1 and the start index of match
def bit_parallel_with_max_2_error(text, pattern, m):
    S_table = defaultdict(int)
    for i, c in enumerate(pattern):
        S_table[c] |= 1 << i
    R0 = 0
    R1 = 0
    R2 = 0
    mask = 1 << (m - 1)
    for j, c in enumerate(text):
        S = S_table[c]
        shR0 = (R0 << 1) | 1
        shR1 = (R1 << 1) | 1
        R0 = shR0 & S
        R1 = ((R1 << 1) | 1) & S | shR0
        R2 = ((R2 << 1) | 1) & S | shR1
        # first if-statement commented because we have already checked for exact match using find() in mmfind()
        #if R0 & mask:  #exact match
            #return 1, j - m + 1
        if R1 & mask:   # 1 mismatch
            return 1, j - m + 1
        if R2 & mask:   # 2 mismatches
            return 1, j - m + 1
    return -1,-1


# checks for a match allowing 1 mismatch
# if not a match returns -1,-1. If match occurs, returns 1 and the start index of match
def bit_parallel_with_max_1_error(text, pattern, m):
    S_table = defaultdict(int)
    for i, c in enumerate(pattern):
        S_table[c] |= 1 << i
    R0 = 0
    R1 = 0
    mask = 1 << (m - 1)
    for j, c in enumerate(text):
        S = S_table[c]
        shR0 = (R0 << 1) | 1
        R0 = shR0 & S
        R1 = ((R1 << 1) | 1) & S | shR0
        # first if-statement commented because we have already checked for exact match using find() in mmfind()
        #if R0 & mask:  #exact match
            #return 1, j - m + 1
        if R1 & mask:   # 1 mismatch
            return 1, j - m + 1
    return -1,-1


# this function is a replacement for below mmfind() for speedup
# it assumes the length of H is <32
# find index of H[1..m] in G[1..n] with up to max (1 or 2) mismatches
def mmfind(G,n,H,m,max): # lengths; assume n>m
    a = G.find(H)
    if a!=-1: return a # shortcut for perfect matches
    a,b = -1,-1
    if max==1 :
        a,b = bit_parallel_with_max_1_error(G, H, m)
    elif max==2 :
        a,b = bit_parallel_with_max_2_error(G, H, m)
    if a==1: return b
    return -1


''' replaced with above mmfind()
# find index of H[1..m] in G[1..n] with up to max mismatches

def mmfind(G,n,H,m,max): # lengths; assume n>m
  a = G[:n].find(H[:m])
  if a!=-1: return a # shortcut for perfect matches
  for i in range(0,n-m):
    cnt = 0
    for k in range(m):
      if G[i+k]!=H[k]: cnt += 1
      if cnt>max: break
    if cnt<=max: return i
  return -1
'''


def extract_staggered(infile,outfile,vars):
  Tn = vars.prefix
  message("prefix sequence: %s" % vars.prefix)
  lenTn = len(Tn)
  ADAPTER2 = "TACCACGACCA"
  lenADAP = len(ADAPTER2)

  #P,Q = 5,10 # 1-based inclusive positions to look for start of Tn prefix
  P,Q = 0,15
  if vars.barseq_catalog_out!=None: Q = 100 # relax for barseq

  vars.tot_tgtta = 0
  vars.truncated_reads = 0
  output = open(outfile,"w")
  tot = 0
  #print infile
  if vars.barseq_catalog_out!=None:
    barcodes_file = vars.base+".barseq" # I could define this in vars
    catalog = open(barcodes_file,"w")
    barseq1 = "TGCAGGGATGTCCACGAGGTCTCT" # const regions surrounding barcode
    barseq2 = "CGTACGCTGCAGGTCGACGGCCGG"
    barseq1len,barseq2len = len(barseq1),len(barseq2)
  for line in open(infile):
    #print line
    line = line.rstrip()
    if not line: continue
    if line[0]=='>': header = line; continue
    tot += 1
    if tot%1000000==0: message("%s reads processed" % tot)
    readlen = len(line)
    a = mmfind(line,readlen,Tn,lenTn,vars.mm1) # allow some mismatches
    b = mmfind(line,readlen,ADAPTER2,lenADAP, 1) # look for end of short frags
    if a>=P and a<=Q:
      gstart,gend = a+lenTn,readlen
      if b!=-1: gend = b; vars.truncated_reads += 1
      if gend-gstart<20: continue # too short # I should make this a param
      output.write(header+"\n")
      output.write(line[gstart:gend]+"\n")
      vars.tot_tgtta += 1
    if vars.barseq_catalog_out!=None:
      n = max(a,readlen)
      c = mmfind(line,n,barseq1,barseq1len,vars.mm1) # only have to search as far as Tn prefix
      d = mmfind(line,n,barseq2,barseq2len,vars.mm1)
      seq = "XXXXXXXXXXXXXXXXXXXX"
      if c!=-1 and d!=-1:
        size = d-c-barseq1len
        if size>=15 and size<=25: seq = line[c+barseq1len:d]
      catalog.write(header+"\n")
      catalog.write(seq+"\n")
  if vars.barseq_catalog_out!=None: catalog.close()
  output.close()
  if vars.tot_tgtta == 0:
    raise ValueError("Error: Input files did not contain any reads matching prefix sequence with %d mismatches" % vars.mm1)


def message(s):
  #print "[tn_preprocess]",s
  #sys.stdout.flush()
  sys.stderr.write("[tn_preprocess] "+s+"\n")

def get_id(line):
  a,b = line.find(":")+1,line.rfind("#")
  if b==-1: b = line.rfind("_")
  return line[a:b]

# select the reads from infile that have headers occuring in goodreads

def select_reads(goodreads,infile,outfile):
  hash = {}
  for line in open(goodreads):
    if line[0]=='>':
      #id = line[line.find(":")+1:line.rfind("#")]
      id = get_id(line)
      hash[id] = 1

  output = open(outfile,"w")
  for line in open(infile):
    if line[0]=='>':
      header = line
      id = get_id(line)
    else:
      if hash.has_key(id):
        output.write(header)
        output.write(line)
  output.close()

def replace_ids(infile1,infile2,outfile):
  f = open(infile1)
  g = open(infile2)
  h = open(outfile,"w")

  while True:
    a = f.readline()
    b = g.readline()
    if len(a)<2: break
    if a[0]=='>': header = a
    else:
      h.write(header)
      h.write(b)
  f.close()
  g.close()
  h.close()

# indexes i and j are 1-based and inclusive (could be -1)

def select_cycles(infile,i,j,outfile):
  output = open(outfile,"w")
  for line in open(infile):
    if line[0]=='>': header = line
    else:
      output.write(header)
      output.write(line[i-1:j]+"\n")
  output.close()

def read_genome(filename):
  s = ""
  for line in open(filename):
    if line[0]=='>': continue # skip fasta header
    else: s += line[:-1]
  return s

# convert to bistring (8 bits; bit 0 is low-order bit)
#
# Bit Description
# 0 0x1 template having multiple segments in sequencing
# 1 0x2 each segment properly aligned according to the aligner
# 2 0x4 segment unmapped
# 3 0x8 next segment in the template unmapped
# 4 0x10 SEQ being reverse complemented
# 5 0x20 SEQ of the next segment in the template being reversed
# 6 0x40 the first segment in the template
# 7 0x80 the last segment in the template
#
# code[6]=1 means read1
# code[4]=1 means reverse strand

def samcode(num): return bin(int(num))[2:].zfill(8)[::-1]

def template_counts(ref,sam,bcfile,vars):
  genome = read_genome(ref)
  barcodes = {}


  fil1 = open(bcfile)
  fil2 = open(sam)

  idx=1
  for line in fil1:
    if idx==1: break
    idx+=1

  idx=1
  for line in fil2:
    if idx==2: break
    idx+=1

  '''
  for line in open(bcfile):
    line = line.rstrip()
    if line[0]=='>': id = line[1:]
    else: barcodes[id] = line
  '''
  hits = {}
  vars.tot_tgtta,vars.mapped = 0,0
  vars.r1 = vars.r2 = 0

  #for line in open(sam):
  bcline=''
  for line in fil2:
    try:
      bcline = fil1.next().rstrip()
      if bcline[0] !='>': bc = bcline
    except StopIteration:
      pass
    if line[0]=='@': continue
    else:
      w = line.split('\t')
      code = samcode(w[1])
      if 'S' in w[5]: continue #elimate softclipped reads
      if code[6]=="1": # previously checked for for reads1's via w[1]<128
        vars.tot_tgtta += 1
        if code[2]=="0": vars.r1 += 1
      if code[7]=="1" and code[2]=="0": vars.r2 += 1
      # include "improperly mapped reads, which might just be short frags
      #if w[1]=="99" or w[1]=="83" or w[1]=="97" or w[1]=="81":
      if code[6]=="1" and code[2]=="0" and code[3]=="0": # both reads mapped (proper or not)
        vars.mapped += 1
        readlen = len(w[9])
        pos,size = int(w[3]),int(w[8]) # note: size could be negative
        strand,delta = 'F',-2
        if code[4]=="1": strand,delta = 'R',readlen

        pos += delta
        #bc = barcodes[w[0]]
        if pos not in hits: hits[pos] = []
        hits[pos].append((strand,size,bc))

  sites = []
  for i in range(len(genome)-1):
    if genome[i:i+2].upper()=="TA":
      pos = i+1
      h = hits.get(pos,[])
      f = filter(lambda x: x[0]=='F',h)
      r = filter(lambda x: x[0]=='R',h)
      h.sort()
      unique = {}
      for (strand,size,bc) in h:
        #print strand,bc,size
        s = "%s-%s-%s" % (strand,bc,size)
        unique[s] = 1
      u = unique.keys()
      uf = filter(lambda x: x[0]=='F',u)
      ur = filter(lambda x: x[0]=='R',u)
      data = [pos,len(f),len(uf),len(r),len(ur),len(f)+len(r),len(uf)+len(ur)]
      sites.append(data)

  return sites # (coord, Fwd_Rd_Ct, Fwd_Templ_Ct, Rev_Rd_Ct, Rev_Templ_Ct, Tot_Rd_Ct, Tot_Templ_Ct)

# pretend that all reads count as unique templates

def increase_counts(pos,sites, strand):
        if strand == "F":
            sites[pos][1] += 1  #if read has been found before, tally 1 more in R reads
            sites[pos][2] += 1  #if read has been found before, tally 1 more in R reads
        if strand == "R":
            sites[pos][3] += 1  #if read has been found before, tally 1 more in R reads
            sites[pos][4] += 1  #if read has been found before, tally 1 more in R reads
        sites[pos][5] += 1  #if read has been found before, tally 1 more in R reads
        sites[pos][6] += 1  #if read has been found before, tally 1 more in R reads


def read_counts(ref,sam,vars):
    genome = read_genome(ref)
    sites = {}
    for i in range(len(genome)-1):
        if genome[i:i+2]=="TA" or vars.transposon=='Tn5':
          pos = i+1
          sites[pos] = [pos,0,0,0,0,0,0]

    hits = {}
    vars.tot_tgtta,vars.mapped = 0,0
    vars.r1 = vars.r2 = 0
    for line in open(sam):
        if line[0]=='@': continue
        else:
            w = line.split('\t')
            code,icode = samcode(w[1]),int(w[1])
            vars.tot_tgtta += 1
            if icode==0 or icode==16:
                vars.r1 += 1
                vars.mapped += 1
                readlen = len(w[9])
                pos = int(w[3])
                if vars.protocol.lower() == "mme1":
                    strand,delta = 'F',readlen
                    if code[4]=="1": strand,delta = 'R',1
                    site1 = pos + delta - 2 #if on + strand, take column 3 position and add 1bp,
                    site2 = pos + delta - 1 #check one off just in case it enzyme chewed too much
                    if site1 in sites:
                        increase_counts(site1, sites, strand)
                    if site2 in sites:
                        increase_counts(site2, sites, strand)
                else:
                    strand,delta = 'F',-2
                    if code[4]=="1": strand,delta = 'R',readlen
                    site1 = pos + delta #if on + strand, take column 3 position and add 1bp)
                    if site1 in sites:
                        increase_counts(site1, sites, strand)

    results = []
    for key in sorted(sites.keys()):
        results.append(sites[key])
    return results # (coord, Fwd_Rd_Ct, Fwd_Templ_Ct, Rev_Rd_Ct, Rev_Templ_Ct, Tot_Rd_Ct, Tot_Templ_Ct)



def driver(vars):
  vars.reads1 = vars.base+".reads1"
  vars.reads2 = vars.base+".reads2"
  vars.trimmed1 = vars.base+".trimmed1"
  vars.trimmed2 = vars.base+".trimmed2"
  vars.barcodes1 = vars.base+".barcodes1"
  vars.barcodes2 = vars.base+".barcodes2"
  vars.genomic2 = vars.base+".genomic2"
  vars.sai1 = vars.base+".sai1"
  vars.sai2 = vars.base+".sai2"
  vars.sam = vars.base+".sam"
  vars.tc = vars.base+".counts"
  vars.wig = vars.base+".wig"
  vars.stats = vars.base+".tn_stats"

  if not vars.prefix:
    if vars.transposon=="Tn5": vars.prefix = "TAAGAGACAG"
    elif vars.transposon=="Himar1": vars.prefix = "ACTTATCAGCCAACCTGTTA"
    else: vars.prefix = ""

  try:
     extract_reads(vars)

     run_bwa(vars)

     generate_output(vars)

  except ValueError as err:
    message("")
    message("%s" % " ".join(err.args))
    message("Exiting.")
    sys.exit()

  except IOError as err:
    message("")
    message("%s" % " ".join(err.args))
    message("Make sure you have read/write access in the directories containing the necessary files.")
    message("Note: If TPP cannot find index files for the FASTA sequence (i.e. *.fna.bwt, *.fna.pac, *.fna.ann, *.fna.sa), it will attempt to create them.")
    message("Exiting.")
    sys.exit()


  message("Done.")

def uncompress(filename):
   outfil = open(filename[0:-3], "w+")
   for line in gzip.open(filename):
      outfil.write(line)
   return filename[0:-3]

def extract_reads(vars):
    message("extracting reads...")

    flag = ['','']
    for idx, name in enumerate([vars.fq1, vars.fq2]):
        if idx==1 and vars.single_end==True: continue
        fil = open(name)
        for line in fil:
            if line[0] == '>':
                flag[idx] = 'FASTA'
                break
            flag[idx] = 'FASTQ'
            break
        fil.close()

    if vars.fq1.endswith('.gz'):
       vars.fq1 = uncompress(vars.fq1)

    if vars.fq2.endswith('.gz'):
       vars.fq2 = uncompress(vars.fq2)

    if(flag[0] == 'FASTQ'):
        message("fastq2reads: %s -> %s" % (vars.fq1,vars.reads1))
        fastq2reads(vars.fq1,vars.reads1,vars.maxreads)
    else:
        shutil.copyfile(vars.fq1, vars.reads1)

    if vars.single_end==True:
      message("assuming single-ended reads")
      message("creating %s" % vars.trimmed1)
      extract_staggered(vars.reads1,vars.trimmed1,vars)

      return

    if(flag[1] == 'FASTQ'):
        message("fastq2reads: %s -> %s" % (vars.fq2,vars.reads2))
        fastq2reads(vars.fq2,vars.reads2,vars.maxreads)
    else:
        shutil.copyfile(vars.fq2, vars.reads2)

    message("fixing headers of paired reads for bwa...")
    fix_paired_headers_for_bwa(vars.reads1,vars.reads2)

    message("extracting barcodes and genomic parts of reads...")

    message("creating %s" % vars.trimmed1)
    extract_staggered(vars.reads1,vars.trimmed1,vars)

    message("creating %s" % vars.trimmed2)
    select_reads(vars.trimmed1,vars.reads2,vars.trimmed2)
    #message("creating %s" % vars.barcodes2)
    #select_cycles(vars.trimmed2,22,30,vars.barcodes2)
    #message("creating %s" % vars.genomic2)
    #select_cycles(vars.trimmed2,43,-1,vars.genomic2)

    # instead of using select_cycles, do these both in one shot by looking for constant seqs
    message("creating %s" % vars.barcodes2)
    message("creating %s" % vars.genomic2)
    extract_barcodes(vars.trimmed2,vars.barcodes2,vars.genomic2, vars.mm1)

    message("creating %s" % vars.barcodes1)
    replace_ids(vars.trimmed1,vars.barcodes2,vars.barcodes1)

#  pattern for read 2...
#    TAGTGGATGATGGCCGGTGGATTTGTG GTAATTACCA TGGTCGTGGTAT CCCAGCGCGACTTCTTCGGCGCACACACC TAACAGGTTGGCTGATAAGTCCCCG?AGAT AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGT
#    -----const1---------------- --barcode- ---const2--- ------genomic---------------- ------const3--------------------------------------------------------------
#    const suffix might appear if fragment is shorter than read length; if so, truncate
#    if genomic part is too short, just output at least 20bp of const so as not to mess up BWA
#    could the start of these be shifted slightly?

def extract_barcodes(fn_tgtta2,fn_barcodes2,fn_genomic2,mm1):
  const1 = "GATGGCCGGTGGATTTGTG"
  const2 = "TGGTCGTGGTAT"
  const3 = "TAACAGGTTGGCTGATAAG"
  nconst1,nconst2,nconst3 = len(const1),len(const2),len(const3)
  fl_barcodes2 = open(fn_barcodes2,"w")
  fl_genomic2 = open(fn_genomic2,"w")
  tot,DEBUG = 0,0
  for line in open(fn_tgtta2):
    line = line.rstrip()
    if line[0]=='>': header = line
    else:
      tot += 1
      if tot%1000000==0: message("%s reads processed" % tot)
      #a  = line.find(const1)
      #b  = line.find(const2)
      #c  = line.find(const3)
      a  = mmfind(line,len(line),const1,nconst1, mm1)
      b  = mmfind(line,len(line),const2,nconst2, mm1)
      c  = mmfind(line,len(line),const3,nconst3, mm1)
      bstart,bend = a+nconst1,b
      gstart,gend = b+nconst2,len(line)
      if c!=-1 and c-gstart>20: gend = c
      if a==-1 or bend<bstart+5 or bend>bstart+15:
        # you can't just reject these, beacuse they are paired with R1
        # but setting the genomic part to the first 20 cycles should prevent it from mapping
        bstart,bend = 0,10
        gstart,gend = 0,20
        barcode,genomic = "XXXXXXXXXX","XXXXXXXXXX"
      else: barcode,genomic = line[bstart:bend],line[gstart:gend]
      if DEBUG==1:
        fl_barcodes2.write(header+"\n")
        fl_barcodes2.write(line+"\n")
        #fl_barcodes2.write((" "*bstart)+line[bstart:bend]+"\n")
        fl_barcodes2.write((" "*bstart)+barcode+"\n")
        fl_genomic2.write(header+"\n")
        fl_genomic2.write(line+"\n")
        #fl_genomic2.write((" "*gstart)+line[gstart:gend]+"\n")
        fl_genomic2.write((" "*gstart)+genomic+"\n")
      else:
        fl_barcodes2.write(header+"\n")
        #fl_barcodes2.write(line[bstart:bend]+"\n")
        fl_barcodes2.write(barcode+"\n")
        fl_genomic2.write(header+"\n")
        #fl_genomic2.write(line[gstart:gend]+"\n")
        fl_genomic2.write(genomic+"\n")
  fl_barcodes2.close()
  fl_genomic2.close()
  if DEBUG==1: sys.exit(0)



def bwa_subprocess(command, outfile):
    commandstr = " ".join(command)
    if outfile.name != "<stdout>":
        commandstr += " > %s" % outfile.name
    message(commandstr)
    process = subprocess.Popen(command, stdout=outfile, stderr=subprocess.PIPE)
    for line in iter(process.stderr.readline, ''):
        if "Permission denied" in line:
            raise IOError("Error: BWA encountered a permissions error: \n\n%s" % line)
        if "invalid option" in line:
            raise ValueError("Error: Unrecognized flag for BWA: %s" % (line.split()[-1]))
        sys.stderr.write("%s\n" % line.strip())




def run_bwa(vars):
    message("mapping reads using BWA...(this takes a couple of minutes)")

    if not os.path.exists(vars.ref+".amb"):
      cmd = [vars.bwa, "index", vars.ref]
      bwa_subprocess(cmd, sys.stdout)


    cmd = [vars.bwa, "aln"]
    if vars.flags.strip():
        cmd.extend( vars.flags.split(" "))
    cmd.extend([vars.ref, vars.trimmed1])
    outfile = open(vars.sai1, "w")
    bwa_subprocess(cmd, outfile)


    if vars.single_end==True:
        cmd = [vars.bwa, "samse", vars.ref, vars.sai1, vars.trimmed1]
        outfile = open(vars.sam, "w")
        bwa_subprocess(cmd, outfile)

    else:

        cmd = [vars.bwa, "aln"]
        if vars.flags.strip():
            cmd.extend(vars.flags.split(" "))
        cmd.extend([vars.ref, vars.genomic2])
        outfile = open(vars.sai2, "w")
        bwa_subprocess(cmd, outfile)

        cmd = [vars.bwa, "sampe", vars.ref, vars.sai1, vars.sai2, vars.trimmed1, vars.genomic2]
        outfile = open(vars.sam, "w")
        bwa_subprocess(cmd, outfile)



def stats(vals):
  sum,ss = 0,0
  for x in vals: sum += x; ss += x*x
  N = float(len(vals))
  mean = sum/N
  var = ss/N-mean*mean
  stdev = math.sqrt(var)
  return mean,stdev

def corr(X,Y):
  muX,sdX = stats(X)
  muY,sdY = stats(Y)

  if sdX == 0 or sdY == 0:
    raise ValueError("Warning: Standard deviations of counts is zero.")

  cX = [x-muX for x in X]
  cY = [y-muY for y in Y]
  s = sum([x*y for (x,y) in zip(cX,cY)])
  return s/(float(len(X))*sdX*sdY)

def get_read_length(filename):
   fil = open(filename)
   i = 0
   for line in fil:
      if i == 1:
         #print "reads1 line: " + line
         return len(line.strip())
      i+=1

def get_genomic_portion(filename):
   fil = open(filename)
   i = 0
   tot_len = 0.0
   n = 1
   for line in fil:
      if i%2 == 1:
         tot_len += len(line.strip())
         n += 1
      i+=1
   return tot_len/n


# return list of (item,cnt) sorted by counts

def popularity(lst):
  hash = {}
  for x in lst:
    if x not in hash: hash[x] = 0
    hash[x] += 1
  data = [(hash[x],x) for x in hash.keys()]
  data.sort(reverse=True)
  data = [(y,x) for (x,y) in data]
  return data


def create_barseq_catalog(vars):
  barcodes = {} # headers->barcodes
  for line in open(vars.base+".barseq"):
    if line[0]=='>': 
      header = line.rstrip()[1:]
      header = header.split()[0] # in case it has a space, which is dropped by bwa in sam file
    else: barcodes[header] = line.split('\n')[0]

  sites,nreads = {},0
  for line in open(vars.sam):
    if line[0]=='@': continue
    w = line.split('\t')
    samcode = int(w[1])
    nreads += 1
    if samcode in [0,16]: # what about PE reads?
      header,coord = w[0],int(w[3])
      # see how I adjust coords in read_counts()
      strand,delta = 'F',-2
      if samcode==16: strand,delta = 'R',len(w[9]); 
      coord += delta
      if strand=='R': coord *= -1
      if coord not in sites: sites[coord] = [] # use -co for minus strand
      sites[coord].append(barcodes[header])

  mapsto = {} # barcodes->sites
  totbc,maptoTAs = 0,0
  genome = read_genome(vars.ref)
  for site,bclist in sites.items():
    for bc in bclist:
      if bc not in mapsto: mapsto[bc] = []
      if 'X' not in bc: 
        mapsto[bc].append(site)
        totbc += 1
        if site<0: site *= -1
        if genome[site-1:site+1]=="TA": maptoTAs += 1

  # good barcodes are those that are associated with only 1 site (must be TA?)
  # what about redundant sites like IS elements?
  goodbc = {}
  for bc,sites2 in mapsto.items(): 
    pop = popularity(sites2) # make a table of locations at which the barcode appears
    #print bc,pop
    if len(pop)==1: goodbc[bc] = 1
  #for x in sorted(sites.items()): print x[0],genome[x[0]-1:x[0]+1],popularity(x[1])

  file = open(vars.barseq_catalog_out,"w")  
  file.write("# Barseq (stats are at the bottom): reads = %s, ref = %s\n" % (vars.fq1,vars.ref))
  n = len(genome)
  a,b = 0,0
  for i in range(n-1):
    if genome[i:i+2]=="TA":
      a += 1
      co = i+1
      barcodes_pos = [(x,'+') for x in sites.get(co,[])]
      barcodes_neg = [(x,'-') for x in sites.get(-co,[])]
      pop = popularity(barcodes_pos+barcodes_neg)
      if len(pop)>0: b += 1
      else: file.write("# %s %s\n" %(co,genome[i:i+2]))
      for (bc,strand),cnt in pop:
        if bc in goodbc: 
          file.write("%s %s %s %s %s\n" % (co,genome[i:i+2],strand,bc,cnt))
  file.write("# total_reads: %s, total_barcodes: %s, map_to_TAs: %s, distinct_bc: %s, unimapped: %s, TA_sites_hit: %s/%s\n" % (nreads,totbc,maptoTAs,len(mapsto.keys()),len(goodbc.keys()),b,a))
  file.close()


def generate_output(vars):
  if vars.barseq_catalog_out!=None: 
    message("creating Barseq catalog file: "+vars.barseq_catalog_out)
    create_barseq_catalog(vars)

  message("tabulating template counts and statistics...")
  if vars.single_end==True: counts = read_counts(vars.ref,vars.sam,vars) # return read counts copied as template counts
  else: counts = template_counts(vars.ref,vars.sam,vars.barcodes1,vars)
  tcfile = open(vars.tc,"w")
  tcfile.write('\t'.join("coord Fwd_Rd_Ct Fwd_Templ_Ct Rev_Rd_Ct Rev_Templ_Ct Tot_Rd_Ct Tot_Templ_Ct".split())+"\n")
  for data in counts: tcfile.write('\t'.join([str(x) for x in data])+"\n")
  tcfile.close()

  if vars.mapped == 0:
    raise ValueError('Error: BWA was unable to map any reads to the genome.')

  message("writing %s" % vars.wig)
  output = open(vars.wig,"w")

  read1 = os.path.basename(vars.fq1)
  read2 = os.path.basename(vars.fq2)
  fi = re.split(r'\.', os.path.basename(vars.ref))[0]
  output.write("# Generated by tpp from " + read1 + " and " + read2 + "\n")
  output.write("variableStep chrom="+ fi + "\n")
  for data in counts: output.write("%s %s\n" % (data[0],data[-1]))
  output.close()

  primer = "CTAGAGGGCCCAATTCGCCCTATAGTGAGT"
  vector = "CTAGACCGTCCAGTCTGGCAGGCCGGAAAC"
  adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  Himar1 = "ACTTATCAGCCAACCTGTTA"

  tot_reads,nprimer,nvector,nadapter,misprimed = 0,0,0,0,0
  for line in open(vars.reads1):
    if line[0]=='>': tot_reads += 1; continue
    if primer in line: nprimer += 1
    if vector in line: nvector += 1
    if adapter in line: nadapter += 1
    if Himar1[:-5] in line and Himar1 not in line: misprimed += 1




  rcounts = [x[5] for x in counts]
  tcounts = [x[6] for x in counts]
  rc,tc = sum(rcounts),sum(tcounts)
  ratio = rc/float(tc) if (rc != 0 and tc !=0) else 0
  ta_sites = len(rcounts)
  tas_hit = len(filter(lambda x: x>0,rcounts))
  density = tas_hit/float(ta_sites)
  counts.sort(key=lambda x: x[-1])
  max_tc = counts[-1][6]
  max_coord = counts[-1][0]
  NZmean = tc/float(tas_hit)

  try:
    FR_corr = corr([x[1] for x in counts],[x[3] for x in counts])
  except ValueError:
    FR_corr = float("nan")
  try:
    BC_corr = corr([x for x in rcounts if x!=0],[x for x in tcounts if x!=0])
  except ValueError:
    BC_corr = float("nan")


  read_length = get_read_length(vars.base + ".reads1")
  mean_r1_genomic = get_genomic_portion(vars.base + ".trimmed1")
  if vars.single_end==False: mean_r2_genomic = get_genomic_portion(vars.base + ".genomic2")

  output = open(vars.stats,"w")
  version = "1.0"
  #output.write("# title: Tn-Seq Pre-Processor, version %s\n" % vars.version)
  output.write("# title: Tn-Seq Pre-Processor\n")
  output.write("# date: %s\n" % time.strftime("%m/%d/%Y %H:%M:%S"))
  output.write("# command: python ")
  output.write(' '.join(sys.argv)+"\n")
  output.write('# transposon type: %s\n' % vars.transposon)
  output.write('# protocol type: %s\n' % vars.protocol)
  output.write('# bwa flags: %s\n' % vars.flags)
  output.write('# read1: %s\n' % vars.fq1)
  output.write('# read2: %s\n' % vars.fq2)
  output.write('# ref_genome: %s\n' % vars.ref)
  output.write("# total_reads %s (or read pairs)\n" % tot_reads)
  #output.write("# truncated_reads %s (fragments shorter than the read length; ADAP2 appears in read1)\n" % vars.truncated_reads)
  output.write("# trimmed_reads %s (reads with valid Tn prefix, and insert size>20bp)\n" % vars.tot_tgtta)
  output.write("# reads1_mapped %s\n" % vars.r1)
  output.write("# reads2_mapped %s\n" % vars.r2)
  output.write("# mapped_reads %s (both R1 and R2 map into genome)\n" % vars.mapped)
  output.write("# read_count %s (TA sites only, for Himar1)\n" % rc)
  output.write("# template_count %s\n" % tc)
  output.write("# template_ratio %0.2f (reads per template)\n" % ratio)
  output.write("# TA_sites %s\n" % ta_sites)
  output.write("# TAs_hit %s\n" % tas_hit)
  output.write("# density %0.3f\n" % density)
  output.write("# max_count %s (among templates)\n" % max_tc)
  output.write("# max_site %s (coordinate)\n" % max_coord)
  output.write("# NZ_mean %0.1f (among templates)\n" % NZmean)
  output.write("# FR_corr %0.3f (Fwd templates vs. Rev templates)\n" % FR_corr)
  output.write("# BC_corr %0.3f (reads vs. templates, summed over both strands)\n" % BC_corr)
  output.write("# primer_matches: %s reads (%0.1f%%) contain %s (Himar1)\n" % (nprimer,nprimer*100/float(tot_reads),primer))
  output.write("# vector_matches: %s reads (%0.1f%%) contain %s (phiMycoMarT7)\n" % (nvector,nvector*100/float(tot_reads),vector))
  output.write("# adapter_matches: %s reads (%0.1f%%) contain %s (Illumina/TruSeq index)\n" % (nadapter,nadapter*100/float(tot_reads),adapter))
  output.write("# misprimed_reads: %s reads (%0.1f%%) contain Himar1 prefix but don't end in TGTTA\n" % (misprimed,misprimed*100/float(tot_reads)))
  output.write("# read_length: %s bp\n" % read_length)
  output.write("# mean_R1_genomic_length: %0.1f bp\n" % mean_r1_genomic)
  if vars.single_end==False: output.write("# mean_R2_genomic_length: %0.1f bp\n" % mean_r2_genomic)

  #output.write("# most_abundant_prefix: %s reads start with %s\n" % (temp[0][1],temp[0][0]))
  # since these are reads (within Tn prefix stripped off), I expect ~1/4 to match Tn prefix
  vals = [vars.fq1,vars.fq2,tot_reads,vars.tot_tgtta,vars.r1,vars.r2,vars.mapped,rc,tc,ratio,ta_sites,tas_hit,max_tc,density,max_coord,NZmean,FR_corr,BC_corr,nprimer,nvector,nadapter,misprimed]
  output.write('\t'.join([str(x) for x in vals])+"\n")
  output.close()

  message("writing %s" % vars.stats)
  #os.system("grep '#' %s" % vars.stats)
  infile = open(vars.stats)
  for line in infile:
      if '#' in line:
          print line.rstrip()
  infile.close()

#############################################################################

def error(s):
  print "error:",s
  sys.exit(0)

def warning(s):
  print "warning:",s





def set_protocol_defaults(vars, protocol):
    #protocol = kwargs.get("protocol", "sassetti")
    if protocol == "sassetti":
        set_sassetti_defaults(vars)
    elif protocol == "mme1":
        set_mme1_defaults(vars)
    elif protocol == "tn5":
        set_tn5_defaults(vars)
    else:
        set_sassetti_defaults(vars)

def set_attributes(vars, attributes_list, override=False):
    for (attr, value) in attributes_list:
        if override:
            setattr(vars, attr, value)
        else:
            if not hasattr(vars, attr):
                setattr(vars, attr, value)


def set_sassetti_defaults(vars):
    attributes_list = []
    attributes_list.append(("transposon", "Himar1"))
    attributes_list.append(("protocol", "Sassetti"))
    attributes_list.append(("prefix", "ACTTATCAGCCAACCTGTTA"))
    attributes_list.append(("maxreads", -1))
    attributes_list.append(("mm1", 100))
    set_attributes(vars, attributes_list)


def set_mme1_defaults(vars):
    attributes_list = []
    attributes_list.append(("transposon", "Himar1"))
    attributes_list.append(("protocol", "Mme1"))
    attributes_list.append(("prefix", ""))
    attributes_list.append(("maxreads", -1))
    attributes_list.append(("mm1", 2))
    set_attributes(vars, attributes_list)


def set_tn5_defaults(vars):
    attributes_list = []
    attributes_list.append(("transposon", "Tn5"))
    attributes_list.append(("protocol", "Tn5"))
    attributes_list.append(("prefix", ""))
    attributes_list.append(("maxreads", -1))
    attributes_list.append(("mm1", 2))
    set_attributes(vars, attributes_list)




def verify_inputs(vars):

    if not os.path.exists(vars.fq1): error("reads1 file not found: "+vars.fq1)
    vars.single_end = False
    if vars.fq2=="": vars.single_end = True
    elif not os.path.exists(vars.fq2): error("reads2 file not found: "+vars.fq2)
    if not os.path.exists(vars.ref): error("reference file not found: "+vars.ref)
    if vars.base == '': error("prefix cannot be empty")
    if vars.fq1 == vars.fq2: error('fastq files cannot be identical')
    if vars.barseq_catalog_in!=None and vars.barseq_catalog_out!=None: error('barseq catalog input and output files cannot both be defined at the same time')

    # If Mme1 protocol, warn that we don't use read2 file
    if vars.protocol.lower() == "mme1" and not vars.single_end:
        warning("Ignoring Read 2 file. TPP assumes Mme1 protocol runs in single-end mode.")
        vars.single_end = True
        vars.fq2 = ""

    if os.path.isdir(vars.bwa):
        bwaexec_unix = os.path.join(vars.bwa, "bwa")
        bwaexec_win = os.path.join(vars.bwa, "bwa.exe")
        if os.path.exists(bwaexec_unix) and not os.path.isdir(bwaexec_unix):
            warning("did not include BWA executable name. Assuming BWA executable is named 'bwa'")
            vars.bwa = bwaexec_unix
        elif os.path.exists(bwaexec_win) and not os.path.isdir(bwaexec_win):
            warning("did not include BWA executable name. Assuming BWA executable is named 'bwa.exe'")
            vars.bwa = bwaexec_win
        else:
            error('cannot find BWA executable. Please include the full executable name as well as its directory.')
    elif not os.path.exists(vars.bwa):
        error('cannot find BWA executable. Please include the full executable name as well as its directory.')

def initialize_globals(vars, args=[], kwargs={}):
    vars.fq1,vars.fq2,vars.ref,vars.bwa,vars.base,vars.maxreads = "","","","","temp",-1
    vars.mm1 = 1 # mismatches allowed in Tn prefix
    vars.transposon = 'Himar1'
    vars.protocol = "Sassetti"
    vars.prefix = "ACTTATCAGCCAACCTGTTA"
    vars.flags = ""
    vars.barseq_catalog_in = vars.barseq_catalog_out = None
    
    # Update defaults
    protocol = kwargs.get("protocol", "").lower()
    if protocol:
        set_protocol_defaults(vars, protocol)
    elif not kwargs:
        read_config(vars)

    # If running in console mode with flags
    if "protocol" in kwargs:
        vars.protocol = kwargs["protocol"]
    if "himar1" in kwargs:
        vars.transposon = "Himar1"
    if "tn5" in kwargs:
        vars.transposon = "Tn5"
    if "protocol" in kwargs:
        vars.protocol = kwargs["protocol"]
    if "primer" in kwargs:
        vars.prefix = kwargs["primer"]
    if "reads1" in kwargs:
        vars.fq1 = kwargs["reads1"]
    if "reads2" in kwargs:
        vars.fq2 = kwargs["reads2"]
    if "bwa" in kwargs:
        vars.bwa = kwargs["bwa"]
    if "ref" in kwargs:
        vars.ref = kwargs["ref"]
    if "maxreads" in kwargs:
        vars.maxreads = int(kwargs["maxreads"])
    if "output" in kwargs:
        vars.base = kwargs["output"]
    if "mismatches" in kwargs:
        vars.mm1 = int(kwargs["mismatches"])
    if "barseq_catalog_in" in kwargs:
        vars.barseq_catalog_in = kwargs["barseq_catalog_in"]
    if "barseq_catalog_out" in kwargs:
        vars.barseq_catalog_out = kwargs["barseq_catalog_out"]
    if "flags" in kwargs:
        vars.flags = kwargs["flags"]
    # note: if last flag expected an arg but was end of list, it gets value True ; check for this and report as missing # TRU, 10/28/17


def read_config(vars):
  if not os.path.exists("tpp.cfg"): return
  for line in open("tpp.cfg"):
    w = line.split()
    if len(w)>=2 and w[0]=='reads1': vars.fq1 = w[1]
    if len(w)>=2 and w[0]=='reads2': vars.fq2 = w[1]
    if len(w)>=2 and w[0]=='ref': vars.ref = w[1]
    if len(w)>=2 and w[0]=='bwa': vars.bwa = w[1]
    if len(w)>=2 and w[0]=='prefix': vars.base = w[1]
    if len(w)>=2 and w[0]=='mismatches1': vars.mm1 = int(w[1])
    if len(w)>=2 and w[0]=='transposon': vars.transposon = w[1]
    if len(w)>=2 and w[0]=='protocol': vars.protocol = " ".join(w[1:])
    if len(w)>=2 and w[0]=='primer': vars.prefix = w[1]
    if len(w)>=2 and w[0]=='flags': vars.flags = " ".join(w[1:])
    if len(w)>=2 and w[0]=='barseq_catalog_in': vars.barseq_catalog_in = w[1]
    if len(w)>=2 and w[0]=='barseq_catalog_out': vars.barseq_catalog_out = w[1]


def save_config(vars):
  f = open("tpp.cfg","w")
  f.write("reads1 %s\n" % vars.fq1)
  f.write("reads2 %s\n" % vars.fq2)
  f.write("ref %s\n" % vars.ref)
  f.write("bwa %s\n" % vars.bwa)
  f.write("prefix %s\n" % vars.base)
  f.write("mismatches1 %s\n" % vars.mm1)
  f.write("transposon %s\n" % vars.transposon)
  f.write("protocol %s\n" % vars.protocol)
  f.write("primer %s\n" % vars.prefix)
  f.write("flags %s\n" % vars.flags)
  if vars.barseq_catalog_in!=None: f.write("barseq_catalog_in %s\n" % vars.barseq_catalog_in)
  if vars.barseq_catalog_out!=None: f.write("barseq_catalog_out %s\n" % vars.barseq_catalog_out)
  f.close()

def show_help():
  print 'usage: python PATH/src/tpp.py -bwa <EXECUTABLE_WITH_PATH> -ref <REF_SEQ> -reads1 <FASTQ_OR_FASTA_FILE> [-reads2 <FASTQ_OR_FASTA_FILE>] -output <BASE_FILENAME> [-maxreads <N>] [-mismatches <N>] [-flags "<STRING>"] [-tn5|-himar1] [-primer <seq>] [-barseq_catalog_in|_out <file>]'

class Globals:
  pass
