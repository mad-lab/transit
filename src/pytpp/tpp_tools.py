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
        # Special case for handling multiple entries for "-ref" and "-replicon_ids"
        if rawargs[count] == "-ref":
          if count+1>=len(rawargs): error("must give comma-separated list as arg for -ref")
          kwargs['ref'] = rawargs[count+1].split(',')
          count += 1
        elif rawargs[count] == "-replicon-ids":
          if count+1>=len(rawargs): error("must give comma-separated list as arg for -replicon-ids")
          kwargs['replicon-ids'] = rawargs[count+1].split(',')
          count += 1
        elif rawargs[count].startswith("-"): #and len(rawargs[count].split(" ")) == 1:
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
   # Check to make sure that headers differ only by 1 character. (Read number. 1 can be single read or Read 2 of paired-end)
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
      # if i==n: raise Exception('Error: unexpected format of headers in .fastq files')
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

##############################

# original implementation
# find index of H[1..m] in G[1..n] with up to max mismatches
# note: this find first match, not necessarily the best (with min mismatches)

def mmfind1(G,n,H,m,max): # lengths; assume n>m
  a = G[:n].find(H[:m])
  if a!=-1: return a # shortcut for perfect matches
  for i in range(0,n-m):
    cnt = 0
    for k in range(m):
      if G[i+k]!=H[k]: cnt += 1
      if cnt>max: break
    if cnt<=max: return i
  return -1


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
def mmfind2(G,n,H,m,max): # lengths; assume n>m
    a = G.find(H)
    if a!=-1: return a # shortcut for perfect matches
    a,b = -1,-1
    if max==1 :
        a,b = bit_parallel_with_max_1_error(G, H, m)
    elif max==2 :
        a,b = bit_parallel_with_max_2_error(G, H, m)
    if a==1: return b
    return -1


# TRI (10/20/2021): I switched back to mmfind1(), since 
#   mmfind2() wasn't working right; increasing -mismatches caused fewer reads to be recognized with prefix and trimmed


def mmfind(G,n,H,m,max): return mmfind1(G,n,H,m,max)

##############################

def windowize(origin, window_size):                         # Generate P,Q values based on a window size (tolerance)
    lower_bound = origin - ((window_size + 0)/2)            # Example, if we assume origin = 28:
    upper_bound = origin + ((window_size + 1)/2)            # for window_size = 0: 28, 28
                                                            #                   1: 28, 29
    return lower_bound, upper_bound                         #                   2: 27, 29 etc.

def extract_staggered(infile,outfile,vars):
  Tn = vars.prefix
  message("prefix sequence: %s" % vars.prefix)
  lenTn = len(Tn)
  ADAPTER2 = "TACCACGACCA"
  lenADAP = len(ADAPTER2)

  #P,Q = 0,15
  #P,Q = 0,50 # relax this, because it has caused problems for various users; shouldn't matter, if prefix is long enough to make random occurences unlikely
  P,Q = vars.primer_start_window

  if vars.window!=None: P,Q = vars.window[0],vars.window[1]
 
  # if -window-size specified, automatically set P,Q by tolerance [RJ]
  if vars.window_size!=-1:
    origin = 28 - len(Tn)                                     
    P,Q = windowize(origin, vars.window_size)

  if vars.barseq_catalog_out!=None: P,Q = 0,100 # relax for barseq
  vars.tot_tgtta = 0
  vars.truncated_reads = 0
  output = open(outfile,"w")
  output_failed = open(outfile+'_failed_trim',"w")                          # [WM] [add]
  if vars.window_size!=-1: message("Looking for start of Tn prefix with P,Q = %d,%d (origin = %d, window size = %d)" % (P,Q,origin,vars.window_size)) # [RJ] Outputting P,Q values and origin/window size
  else: message("Looking for start of Tn prefix within positions [%d,%d]" % (P,Q))
  tot = 0
  if vars.barseq_catalog_out!=None:
    barcodes_file = vars.base+".barseq" # I could define this in vars
    catalog = open(barcodes_file,"w")
    barseq1 = "TGCAGGGATGTCCACGAGGTCTCT" # const regions surrounding barcode
    barseq2 = "CGTACGCTGCAGGTCGACGGCCGG"
    barseq1len,barseq2len = len(barseq1),len(barseq2)
  for line in open(infile):
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
      minReadLen = 15 if vars.protocol.lower() == "mme1" else 20
      if gend-gstart < minReadLen: continue # too short # I should make this a param
      output.write(header+"\n")
      output.write(line[gstart:gend]+"\n")
      vars.tot_tgtta += 1
    else:                                               # [WM] [add]
      #Output reads that failed to be trimmed.          # [WM] [add]
      output_failed.write(header+"\n")                  # [WM] [add]
      output_failed.write(line+"\n")                    # [WM] [add]
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
  output_failed.close()                                 # [WM] [add]
  if vars.tot_tgtta == 0:
    raise ValueError("Error: Input files did not contain any reads matching prefix sequence with %d mismatches" % vars.mm1)
  vars.tot_reads = tot

def message(s):
  #print("[tn_preprocess]",s)
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
      if id in hash:
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

def read_genome(filename, replicon_index):
  s = ""
  cur_index = 0
  first_iteration = True
  for line in open(filename):
    if line[0]=='>': 
      if first_iteration:
        continue # skip fasta header
      else:
        cur_index += 1
    else:
      if cur_index == replicon_index:
        s += line.strip()
    first_iteration = False
  return s.upper()


def parse_sam_header(sam_filename):
    parsed_header = {}
    
    f = open(sam_filename, "r")
    for line in f:

        # Stop parsing once end of SAM header is reached
        if not line[0] == '@':
            break
        
        header_line = line.split()
        at_sign_sam_tag = header_line[0]
        if at_sign_sam_tag not in parsed_header:
            parsed_header[at_sign_sam_tag] = []

        tag_line = {}
        for tag_entry in header_line[1:]:
            tag_pair = tag_entry.split(':')
            tag_name = tag_pair[0]
            tag_value = ':'.join(tag_pair[1:])  # in case the value of the tag contains ':', which we used as delimiter for split
            tag_line[tag_name] = tag_value
        parsed_header[at_sign_sam_tag].append(tag_line)
    
    f.close()

    return parsed_header


def get_replicon_names_from_sam_header(sam_header):
    
    if "@SQ" not in sam_header:
        raise ValueError("No @SQ tags found in SAM header")

    replicon_names = []
    
    for header_line in sam_header["@SQ"]:
        if "SN" not in header_line:
            raise ValueError("No SN tag found in an entry of the SAM header")
        else:
            replicon_names.append(header_line["SN"])
    
    return replicon_names


def count_sam_header_lines(sam_header):
    
    num_sam_header_lines = 0

    for at_sign_sam_tag in sam_header:
        for line in sam_header[at_sign_sam_tag]:
            num_sam_header_lines += 1

    return num_sam_header_lines


# convert to bistring (16 bits; bit 0 is low-order bit)
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

def samcode(num): 
    num = int(num)
    binary_num = bin(num)
    truncated_0b_binary_num = binary_num[2:]
    zero_padded_binary_num = truncated_0b_binary_num.zfill(16)
    rev = zero_padded_binary_num[::-1]
    return rev

def template_counts(ref,sam,bcfile,vars):
  vars.mapped = vars.r1 = vars.r2 = 0

  barcodes = {}
  for line in open(bcfile):
    line = line.rstrip()
    if line[0]=='>': id = line[1:]
    else: barcodes[id] = line
  
  hits = {}
  sam_header = parse_sam_header(sam)
  replicon_names = get_replicon_names_from_sam_header(sam_header)
  for replicon in replicon_names:
    hits[replicon] = {}
  
  skip = count_sam_header_lines(sam_header)
  for line in open(sam):
    #if line[0]=='@': continue
    if skip>0: skip -= 1; continue
    else:
      w = line.split('\t')
      bc = barcodes[w[0]]
      code = samcode(w[1])
      if 'S' in w[5]: continue # eliminate softclipped reads
      if code[6]=="1" and code[2]=="0": vars.r1 += 1
      if code[6]=="1" and code[3]=="0": vars.r2 += 1
      if bc=="XXXXXXXXXX": continue
      if code[6]=="1" and code[1]=="1": # both reads map properly (83 or 99) and has legit barcode
        vars.mapped += 1
        readlen = len(w[9])
        pos,size = int(w[3]),int(w[8]) # note: size could be negative
        strand,delta = 'F',-2
        if code[4]=="1": strand,delta = 'R',readlen

        pos += delta
        replicon_name = w[2]
        if pos not in hits[replicon_name]: hits[replicon_name][pos] = []
        hits[replicon_name][pos].append((strand,size,bc))

  sites_list = []
  for replicon_index in range(vars.num_replicons):
    sites = []
    genome = read_genome(ref, replicon_index)
    for i in range(len(genome)-1):
      #if genome[i:i+2].upper()=="TA":
      if vars.transposon=="Himar1" and genome[i:i+2].upper()!="TA": continue 
      else:
        pos = i+1
        h = hits[replicon_names[replicon_index]].get(pos,[])
        f = list(filter(lambda x: x[0]=='F',h))
        r = list(filter(lambda x: x[0]=='R',h))
        u = list(set(h))
        uf = list(filter(lambda x: x[0]=='F',u))
        ur = list(filter(lambda x: x[0]=='R',u))
        data = [pos,len(f),len(uf),len(r),len(ur),len(f)+len(r),len(uf)+len(ur)]
        sites.append(data)    
    sites_list.append(sites)

  return sites_list     # list of (coord, Fwd_Rd_Ct, Fwd_Templ_Ct, Rev_Rd_Ct, Rev_Templ_Ct, Tot_Rd_Ct, Tot_Templ_Ct)

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
    sites_dict = {}
    sam_header = parse_sam_header(sam)
    replicon_names = get_replicon_names_from_sam_header(sam_header)
    
    for replicon_names_index in range(vars.num_replicons):
        sites = {}
        genome = read_genome(ref, replicon_names_index)
        for i in range(len(genome)-1):
            #if genome[i:i+2]=="TA" or vars.transposon=='Tn5':
            if vars.transposon=='Himar1' and genome[i:i+2]!="TA": continue
            pos = i+1
            sites[pos] = [pos,0,0,0,0,0,0]
            sites_dict[replicon_names[replicon_names_index]] = sites

    vars.tot_tgtta = 0
    vars.mapped = 0
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
                replicon_name = w[2]

                if vars.protocol.lower() == "mme1":
                    strand,delta = 'F',readlen
                    if code[4]=="1": strand,delta = 'R',1
                    site1 = pos + delta - 2 #if on + strand, take column 3 position and add 1bp,
                    site2 = pos + delta - 1 #check one off just in case it enzyme chewed too much
                    if site1 in sites_dict[replicon_name]:
                        increase_counts(site1, sites_dict[replicon_name], strand)
                    if site2 in sites_dict[replicon_name]:
                        increase_counts(site2, sites_dict[replicon_name], strand)
                else:
                    strand,delta = 'F',-2
                    if code[4]=="1": strand,delta = 'R',readlen
                    site1 = pos + delta #if on + strand, take column 3 position and add 1bp)
                    if site1 in sites_dict[replicon_name]:
                        increase_counts(site1, sites_dict[replicon_name], strand)

    results_list = []
    for replicon_index in range(vars.num_replicons):
        results = []
        for key in sorted(sites_dict[replicon_names[replicon_index]].keys()):
            results.append(sites_dict[replicon_names[replicon_index]][key])
        results_list.append(results)
    return results_list # list of (coord, Fwd_Rd_Ct, Fwd_Templ_Ct, Rev_Rd_Ct, Rev_Templ_Ct, Tot_Rd_Ct, Tot_Templ_Ct)


def driver(vars):
  # [RJ] These variables are for the extract_reads() step (no reference genome involved)
  vars.reads1 = vars.base+".reads1"
  vars.reads2 = vars.base+".reads2"
  vars.trimmed1 = vars.base+".trimmed1"         # Final fastq for read 1 is stored in this file
  vars.trimmed2 = vars.base+".trimmed2"
  vars.barcodes1 = vars.base+".barcodes1"
  vars.barcodes2 = vars.base+".barcodes2"
  vars.genomic2 = vars.base+".genomic2"         # Final fastq for read 2 is stored in this file

  # [RJ] These variables are for the run_bwa() step (reference genome involved)
  vars.sai1 = vars.base+".sai1"                 # [RJ] SAI only used when using bwa_alg 'aln'
  vars.sai2 = vars.base+".sai2"
  vars.sam = vars.base+".sam"

  # [RJ] These variables are for the generate_output() step
  vars.tc = []
  vars.wig = []
  vars.stats = vars.base+".tn_stats"

  # [RJ] Handle multi-line fastas by making a multiline fasta out of all specified references
  if len(vars.ref) > 1:
    total_num_records = 0
    fasta_names_combined = ""
    contents_of_multiline_fasta = ""
    # Create multi-line FASTA if more than one reference file is specified
    for reference_genome in vars.ref:
      fasta_names_combined += ("%s%s" % ("_", os.path.splitext(os.path.basename(reference_genome))[0])) 
      with open(reference_genome, "r") as ref:
        contents = ref.read()
      if not contents.endswith('\n'):
        contents += '\n'
      contents_of_multiline_fasta += contents
      num_records = contents.count('>')
      total_num_records += num_records
    fasta_names_combined = os.path.join(os.path.dirname(vars.ref[0]), fasta_names_combined[1:] + ".fa")    # Remove first _
    with open(fasta_names_combined, "w") as multiline_fasta:
      multiline_fasta.write(contents_of_multiline_fasta)
    message("%d records from %s combined into single .fa file: %s" % (total_num_records, ", ".join(vars.ref), fasta_names_combined))
    vars.ref = fasta_names_combined
  else:
    with open(vars.ref[0], "r") as ref:
      contents = ref.read()
    num_records = contents.count('>')
    total_num_records = num_records
    vars.ref = vars.ref[0]
    message("One reference genome specified: %s, containing %d records." % (vars.ref, num_records))
  vars.num_replicons = total_num_records

  if vars.num_replicons != len(vars.replicon_ids):
    if vars.num_replicons is 1:
      vars.replicon_ids = ['']
    # Autogenerate ids if 'auto' flag present
    elif len(vars.replicon_ids) is 1 and vars.replicon_ids[0].strip() == 'auto':
      message("Autogenerating replicon_ids...")
      vars.replicon_ids = [str(i) for i in range(1, vars.num_replicons + 1)]
    else:
      raise error("%d replicons detected in reference genome, but only %d replicon_ids specified" % (vars.num_replicons, len(vars.replicon_ids)))

  if len(vars.replicon_ids) > 1:
   for name in vars.replicon_ids:
    vars.tc.append(vars.base+"_"+name+".counts")
    vars.wig.append(vars.base+"_"+name+".wig")
  else:
    vars.tc.append(vars.base+".counts")
    vars.wig.append(vars.base+".wig")
  
  # Handled independently for CLI and GUI.
  # if not vars.prefix:
  #   if vars.transposon=="Tn5": vars.prefix = "TAAGAGACAG"
  #   elif vars.transposon=="Himar1": vars.prefix = "ACTTATCAGCCAACCTGTTA" # [ORIGINAL]
  #   else: vars.prefix = ""

  try:
    extract_reads(vars)
    
    run_bwa(vars)
    
    generate_output(vars)

  except ValueError as err:
    message("")
    message(err.args)
    message("Exiting.")
    sys.exit()

  except IOError as err:
    message("")
    message("%s" % " ".join(str(v) for v in err.args))     # [RJ] Fixed to prevent erroring out on numeric arguments
    message("Make sure you have read/write access in the directories containing the necessary files.")
    message("Note: If TPP cannot find index files for the FASTA sequence (i.e. *.fna.bwt, *.fna.pac, *.fna.ann, *.fna.sa), it will attempt to create them.")
    message("Exiting.")
    sys.exit()

  message("Done.")

# gunzip <fastq>.gz file, written to <fastq>

def uncompress(filename):
   newfname = filename[0:-3]
   print("uncompressing %s" % filename)
   outfil = open(newfname, "wb+")
   for line in gzip.open(filename):
      outfil.write(line)
   return newfname

def copy_fasta(infile,outfile,maxreads=-1):
  a = open(infile)
  b = open(outfile,"w")
  cnt = 0
  while True:
    hdr = a.readline()
    seq = a.readline()
    if len(hdr)==0: break
    b.write(hdr)
    b.write(seq)
    cnt += 1
    if maxreads>0 and cnt>=maxreads:break
  a.close()
  b.close()

def extract_reads(vars):
    message("extracting reads...")

    if vars.fq1.endswith('.gz'):
       vars.fq1 = uncompress(vars.fq1)

    if vars.fq2.endswith('.gz'):
       vars.fq2 = uncompress(vars.fq2)

    flag = ['','']
    for idx, name in enumerate([vars.fq1, vars.fq2]):
        print(name)
        if idx==1 and vars.single_end==True: continue
        fil = open(name)
        for line in fil:
            if line[0] == '>':
                flag[idx] = 'FASTA'
                break
            flag[idx] = 'FASTQ'
            break
        fil.close()

    if(flag[0] == 'FASTQ'):
        message("fastq2reads: %s -> %s" % (vars.fq1,vars.reads1))
        fastq2reads(vars.fq1,vars.reads1,vars.maxreads)
    else:
        #shutil.copyfile(vars.fq1, vars.reads1)
        copy_fasta(vars.fq1, vars.reads1, vars.maxreads)

    if vars.single_end==True:
      message("assuming single-ended reads")
      message("creating %s" % vars.trimmed1)
      extract_staggered(vars.reads1,vars.trimmed1,vars)

      return

    if(flag[1] == 'FASTQ'):
        message("fastq2reads: %s -> %s" % (vars.fq2,vars.reads2))
        fastq2reads(vars.fq2,vars.reads2,vars.maxreads)
    else:
        #shutil.copyfile(vars.fq2, vars.reads2)
        copy_fasta(vars.fq2, vars.reads2, vars.maxreads)

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
      if len(genomic)==0:
          genomic = "XXX" #Necessary to avoid a bizarre error with bwa when there is an empty line.
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
    #process.wait()
    (pout,perr) = process.communicate()
    #for line in iter(process.stderr.readline, ''):
    for line in perr.split(b'\n'): # returned by communicate()
        if b"Permission denied" in line:
            raise IOError("Error: BWA encountered a permissions error: \n\n%s" % line)
        if b"invalid option" in line:
            raise ValueError("Error: Unrecognized flag for BWA: %s" % (line.split()[-1]))
        sys.stderr.write("%s\n" % line.strip())



def run_bwa(vars):
    message("mapping reads using BWA...(this takes a couple of minutes)")

    # Create index if it doesn't already exist (.amb file is created)
    if not os.path.exists(vars.ref+".amb"):
      cmd = [vars.bwa, "index", vars.ref]
      bwa_subprocess(cmd, sys.stdout)
    
    
    cmd = [vars.bwa, vars.bwa_alg]
    if vars.flags.strip():
        cmd.extend(vars.flags.split(" "))
    cmd.extend([vars.ref, vars.trimmed1])
    

    if vars.bwa_alg == "mem":
        if vars.single_end == False:
            cmd.extend([vars.genomic2])
        outfile = open(vars.sam, "w")
        bwa_subprocess(cmd, outfile)
        outfile.close()

    elif vars.bwa_alg == "aln":
        outfile = open(vars.sai1, "w")
        bwa_subprocess(cmd, outfile)
        outfile.close()

        if vars.single_end==True:
            cmd = [vars.bwa, "samse", vars.ref, vars.sai1, vars.trimmed1]
            outfile = open(vars.sam, "w")
            bwa_subprocess(cmd, outfile)
            outfile.close()

        else:
            cmd = [vars.bwa, vars.bwa_alg]
            if vars.flags.strip():
                cmd.extend(vars.flags.split(" "))
            
            cmd.extend([vars.ref, vars.genomic2])
            outfile = open(vars.sai2, "w")
            bwa_subprocess(cmd, outfile)
            outfile.close()
            
            cmd = [vars.bwa, "sampe", vars.ref, vars.sai1, vars.sai2, vars.trimmed1, vars.genomic2]
            outfile = open(vars.sam, "w")
            bwa_subprocess(cmd, outfile)
            outfile.close()
    else:
        raise ValueError("Error: Invalid BWA algorithm '%s' specified, acceptable algorithms are 'aln' and 'mem'" % vars.bwa_alg)


def stats(vals):
  N = float(len(vals))
  tot,ss = 0,0
  if N == 0:
      return 0, 0
  for x in vals: tot += x; ss += x*x
  mean = tot/N
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
         #print("reads1 line: " + line)
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

#TODO: fix this...?
def create_barseq_catalog(vars, replicon_index):
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
  genome = read_genome(vars.ref, replicon_index)
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
    #print(bc,pop)
    if len(pop)==1: goodbc[bc] = 1
  #for x in sorted(sites.items()): print(x[0],genome[x[0]-1:x[0]+1],popularity(x[1]))

  file = open(vars.barseq_catalog_out[replicon_index],"w")  
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
    message("creating Barseq catalog file: "+vars.barseq_catalog_out) # [RJ] TODO: This should be indexed by replicon_index too... maybe
    create_barseq_catalog(vars)
  
  rcounts,tcounts,rc,tc,ratio,ta_sites,tas_hit,density,max_tc,max_coord,NZmean,FR_corr,BC_corr = [],[],[],[],[],[],[],[],[],[],[],[],[]

  message("tabulating template counts and statistics for reference genome %s" % vars.ref)
  # message("tabulating template counts and statistics...")
  if vars.single_end==True: counts = read_counts(vars.ref,vars.sam,vars) # return read counts copied as template counts
  else: counts = template_counts(vars.ref,vars.sam,vars.barcodes1,vars)
  for replicon_index in range(vars.num_replicons):
    tcfile = open(vars.tc[replicon_index],"w")
    tcfile.write('\t'.join("coord Fwd_Rd_Ct Fwd_Templ_Ct Rev_Rd_Ct Rev_Templ_Ct Tot_Rd_Ct Tot_Templ_Ct".split())+"\n")
    for data in counts[replicon_index]: tcfile.write('\t'.join([str(x) for x in data])+"\n")
    tcfile.close()

    if vars.mapped == 0:
      raise ValueError('Error: BWA was unable to map any reads to the genome.')

    message("writing %s" % vars.wig[replicon_index])
    output = open(vars.wig[replicon_index],"w")

    read1 = os.path.basename(vars.fq1)
    read2 = os.path.basename(vars.fq2)
    fi = re.split(r'\.', os.path.basename(vars.ref))[0]
    output.write("# Generated by tpp from " + read1 + " and " + read2 + "\n")
    output.write("variableStep chrom="+ fi)
    if vars.num_replicons > 1:
      output.write(", replicon=%d" % replicon_index)
    output.write("\n")
    for data in counts[replicon_index]: output.write("%s %s\n" % (data[0],data[-1]))    # This is where the wig is actually written
    output.close()

    cur_rcounts = [x[5] for x in counts[replicon_index]]
    cur_tcounts = [x[6] for x in counts[replicon_index]]
    cur_rc = sum(cur_rcounts)
    cur_tc = sum(cur_tcounts)
    cur_ratio = cur_rc/float(cur_tc) if (cur_rc != 0 and cur_tc !=0) else 0
    cur_ta_sites = len(cur_rcounts)
    cur_tas_hit = len(list(filter(lambda x: x>0,cur_rcounts)))
    cur_density = cur_tas_hit/float(cur_ta_sites) if cur_tas_hit != 0 else 0
    counts[replicon_index].sort(key=lambda x: x[-1])
    cur_max_tc = counts[replicon_index][-1][6]
    cur_max_coord = counts[replicon_index][-1][0]
    cur_NZmean = cur_tc/float(cur_tas_hit) if cur_tas_hit != 0 else 0

    try:
      cur_FR_corr = corr([x[1] for x in counts[replicon_index]],[x[3] for x in counts[replicon_index]])
    except ValueError:
      cur_FR_corr = float("nan")
    try:
      cur_BC_corr = corr([x for x in cur_rcounts if x!=0],[x for x in cur_tcounts if x!=0])
    except ValueError:
      cur_BC_corr = float("nan")
    
    rcounts.append(cur_rcounts)
    tcounts.append(cur_tcounts)
    rc.append(cur_rc)
    tc.append(cur_tc)
    ratio.append(cur_ratio)
    ta_sites.append(cur_ta_sites)
    tas_hit.append(cur_tas_hit)
    density.append(cur_density)
    max_tc.append(cur_max_tc)
    max_coord.append(cur_max_coord)
    NZmean.append(cur_NZmean)
    FR_corr.append(cur_FR_corr)
    BC_corr.append(cur_BC_corr)

  tot_reads = vars.tot_reads
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
  output.write('# replicon_ids: %s\n' % ','.join(vars.replicon_ids))
  output.write("# total_reads (or read pairs): %s\n" % tot_reads)
  output.write("# truncated_reads %s (genomic inserts shorter than the read length; ADAP2 appears in read1)\n" % vars.truncated_reads)
  output.write("# trimmed_reads (reads with valid Tn prefix, and insert size>20bp): %s\n" % vars.tot_tgtta)
  output.write("# reads1_mapped: %s\n" % vars.r1)
  output.write("# reads2_mapped: %s\n" % vars.r2)
  output.write("# mapped_reads (both R1 and R2 map into genome, and R2 has a proper barcode): %s \n" % vars.mapped)

  if vars.num_replicons>1:
    output.write("# read_count (TA sites only, for Himar1):\n")
    for replicon_ids,read_count in zip(vars.replicon_ids, rc):
      output.write("#   %s: %s\n" % (replicon_ids, read_count))
    output.write("# template_count:\n")
    for replicon_ids,template_count in zip(vars.replicon_ids, tc):
      output.write("#   %s: %s\n" % (replicon_ids, template_count))
    output.write("# template_ratio (reads per template):\n")
    for replicon_ids,template_ratio in zip(vars.replicon_ids, ratio):
      output.write("#   %s: %0.2f\n" % (replicon_ids, template_ratio))
    output.write("# TA_sites:\n")
    for replicon_ids,num_ta_sites in zip(vars.replicon_ids, ta_sites):
      output.write("#   %s: %s\n" % (replicon_ids, num_ta_sites))
    output.write("# TAs_hit:\n")
    for replicon_ids,num_tas_hit in zip(vars.replicon_ids, tas_hit):
      output.write("#   %s: %s\n" % (replicon_ids, num_tas_hit))
    output.write("# density:\n")
    for replicon_ids,dens in zip(vars.replicon_ids, density):
      output.write("#   %s: %0.3f\n" % (replicon_ids, dens))
    output.write("# max_count (among templates):\n")
    for replicon_ids,max_template_counts in zip(vars.replicon_ids, max_tc):
      output.write("#   %s: %s\n" % (replicon_ids, max_template_counts))
    output.write("# max_site (coordinate):\n")
    for replicon_ids,max_site in zip(vars.replicon_ids, max_coord):
      output.write("#   %s: %s\n" % (replicon_ids, max_site))
    output.write("# NZ_mean (among templates):\n")
    for replicon_ids,nzmean in zip(vars.replicon_ids, NZmean):
      output.write("#   %s: %0.1f\n" % (replicon_ids, nzmean))
    output.write("# FR_corr (Fwd templates vs. Rev templates):\n")
    for replicon_ids,frcorr in zip(vars.replicon_ids, FR_corr):
      output.write("#   %s: %0.3f\n" % (replicon_ids, frcorr))
    output.write("# BC_corr (reads vs. templates, summed over both strands):\n")
    for replicon_ids,bccorr in zip(vars.replicon_ids, BC_corr):
      output.write("#   %s: %0.3f\n" % (replicon_ids, bccorr))

  else: # just one replicon (contig); this format is read to populate table of stats for datasets in TPP GUI ("# <var> <val>")
    output.write("# read_count (TA sites only, for Himar1): %s\n" % rc[0])
    output.write("# template_count: %s\n" % tc[0])
    output.write("# template_ratio (reads per template): %0.2f\n" % ratio[0])
    output.write("# TA_sites: %s\n" % ta_sites[0])
    output.write("# TAs_hit: %s\n" % tas_hit[0])
    output.write("# density: %0.3f\n" % density[0])
    output.write("# max_count (among templates): %s\n" % max_tc[0])
    output.write("# max_site (coordinate): %s\n" % max_coord[0])
    output.write("# NZ_mean (among templates): %0.1f\n" % NZmean[0])
    output.write("# FR_corr (Fwd templates vs. Rev templates): %0.3f\n" % FR_corr[0])
    output.write("# BC_corr (reads vs. templates, summed over both strands): %0.3f\n" % BC_corr[0])

  primer = "CTAGAGGGCCCAATTCGCCCTATAGTGAGT"
  vector = "CTAGACCGTCCAGTCTGGCAGGCCGGAAAC"
  adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  ADAPTER2 = "TACCACGACCA" # rc of const2 region of R2, between barcode and genomic; these reads will be truncated here
  Himar1 = "ACTTATCAGCCAACCTGTTA"
  trimmed_reads,nprimer,nvector,nadapter,misprimed,ntruncated = 0,0,0,0,0,0
  for line in open(vars.trimmed1):
    if line[0]=='>': trimmed_reads += 1; continue
    if primer in line: nprimer += 1
    if vector in line: nvector += 1
    if adapter in line: nadapter += 1
    #if "TGTTA" in line and Himar1 not in line: misprimed += 1
    # basically, these should correspond to insertions at non-TA sites (so the terminal TA of ...TGTTA will be different)
    if Himar1[:-5] in line and Himar1 not in line: misprimed += 1 

  output.write("# Break-down of total reads (%s):\n" % tot_reads)
  output.write("#  %s reads (%0.1f%%) lack the expected Tn prefix\n" % (tot_reads-vars.tot_tgtta,(tot_reads-vars.tot_tgtta)*100/float(tot_reads)))

  output.write("# Break-down of trimmed reads with valid Tn prefix (%s):\n" % trimmed_reads)
  output.write("#  primer_matches: %s reads (%0.1f%%) contain %s (Himar1)\n" % (nprimer,nprimer*100/float(trimmed_reads),primer))
  output.write("#  vector_matches: %s reads (%0.1f%%) contain %s (phiMycoMarT7)\n" % (nvector,nvector*100/float(trimmed_reads),vector))
  output.write("#  adapter_matches: %s reads (%0.1f%%) contain %s (Illumina/TruSeq index)\n" % (nadapter,nadapter*100/float(trimmed_reads),adapter))
  output.write("#  misprimed_reads: %s reads (%0.1f%%) contain Himar1 prefix but don't end in TGTTA\n" % (misprimed,misprimed*100/float(trimmed_reads)))

  output.write("# read_length: %s bp\n" % read_length)
  output.write("# mean_R1_genomic_length: %0.1f bp\n" % mean_r1_genomic)
  if vars.single_end==False: output.write("# mean_R2_genomic_length: %0.1f bp\n" % mean_r2_genomic)

  output.close()

  message("writing %s" % vars.stats)
  #os.system("grep '#' %s" % vars.stats[replicon_index])
  infile = open(vars.stats)
  for line in infile:
      if '#' in line:
          print(line.rstrip())
  infile.close()

#############################################################################

def error(s):
  print("error:",s)
  sys.exit(0)

def warning(s):
  print("warning:",s)





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
    
    if not vars.ref.__class__.__name__ == "list":
        vars.ref = [vars.ref]
    for ref_genome in vars.ref:
        if not os.path.exists(ref_genome): error("reference file not found: "+ref_genome)
    if vars.base == '': error("prefix cannot be empty")
    # Empty prefix is allowed (03/04/2019 - When reads are pretrimmed)
    # if len(vars.prefix)==0: error("primer sequence cannot be empty")
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
    vars.fq1,vars.fq2,vars.ref,vars.bwa,vars.bwa_alg,vars.replicon_ids,vars.base,vars.maxreads = "","","","","","","",-1
    vars.mm1 = 1 # mismatches allowed in Tn prefix AND adapter prefix on read2
    vars.transposon = 'Himar1'
    vars.protocol = "Sassetti"
    vars.prefix = "ACTTATCAGCCAACCTGTTA"
    vars.flags = ""
    vars.barseq_catalog_in = vars.barseq_catalog_out = None
    vars.window_size = -1
    vars.primer_start_window = 0,20
    vars.window = None
    vars.bwa_alg = "aln" # changing from mem back to aln because of /dev/shm error on Windows machines [TRI,9/14/24]
    
    # Update defaults
    protocol = kwargs.get("protocol", "").lower()
    if protocol:
        set_protocol_defaults(vars, protocol)
    elif not kwargs:
        read_config(vars)

    # If running in console mode with flags
    if "protocol" in kwargs:
        protocol = kwargs["protocol"].lower().capitalize()
        if protocol not in "Sassetti Tn5 Mme1".split(): error("protocol must be one of: Sassetti, Tn5, Mme1")
        vars.protocol = protocol
        if vars.protocol in ["Sassetti","Mme1"]: vars.transposon = "Himar1"
        if vars.protocol=="Tn5": vars.transposon = "Tn5"
    if "himar1" in kwargs:
        vars.transposon = "Himar1"
    if "tn5" in kwargs:
        vars.transposon = "Tn5"
    if "primer" in kwargs:
        vars.prefix = kwargs["primer"].strip()
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

    if "window-size" in kwargs:                             # [RJ] Adding support for window-size, which is the tolerance of positions for the Tn prefix
        vars.window_size = int(kwargs["window-size"])
        if vars.window_size < 6:
            raise ValueError("Error: window-size cannot be less than 6")

    if "primer-start-window" in kwargs:
        w = kwargs["primer-start-window"]
        w = w.split(',')
        vars.window = (int(w[0]),int(w[1]))

    if "bwa-alg" in kwargs:
        if kwargs["bwa-alg"] not in [ "mem", "aln" ]:
            raise ValueError("Error: bwa-alg can only be 'mem' or 'aln'")
        else:
            vars.bwa_alg = kwargs["bwa-alg"]
   
    if "replicon-ids" in kwargs:
        vars.replicon_ids = kwargs["replicon-ids"]

    # Handle no primer but Tn5 protocol case
    # Handled in initialize_globals and tpp_gui
    if "primer" not in kwargs and vars.transposon == "Tn5":
        vars.prefix = "TAAGAGACAG"
    if "primer" not in kwargs and vars.transposon == "Himar1":
        vars.prefix = "ACTTATCAGCCAACCTGTTA"

    # note: if last flag expected an arg but was end of list, it gets value True ; check for this and report as missing # TRI, 10/28/17

def read_config(vars):
  if not os.path.exists("tpp.cfg"): return
  for line in open("tpp.cfg"):
    w = line.split()
    if len(w)>=2 and w[0]=='reads1': vars.fq1 = w[1]
    if len(w)>=2 and w[0]=='reads2': vars.fq2 = w[1]
    if len(w)>=2 and w[0]=='ref': vars.ref = ' '.join(w[1:])
    if len(w)>=2 and w[0]=='ids': vars.replicon_ids = w[1] #vars.replicon_ids = ','.join(w[1:])
    if len(w)>=2 and w[0]=='bwa': vars.bwa = w[1]
    if len(w)>=2 and w[0]=='bwa-alg': vars.bwa_alg = w[1]
    if len(w)>=2 and w[0]=='flags': vars.flags = " ".join(w[1:])
    #if len(w)>=2 and w[0]=='prefix': vars.base = w[1]
    if len(w)>=2 and w[0]=='mismatches1': vars.mm1 = int(w[1])
    if len(w)>=2 and w[0]=='maxreads': vars.maxreads = int(w[1])
    if len(w)>=2 and w[0]=='window_size': vars.window_size = int(w[1])
    if len(w)>=2 and w[0]=='primer_start_window': v = w[1].split(','); vars.primer_start_window = (int(v[0]),int(v[1]))
    if len(w)>=2 and w[0]=='transposon': vars.transposon = w[1]
    if len(w)>=2 and w[0]=='protocol': vars.protocol = " ".join(w[1:])
    if len(w)>=2 and w[0]=='primer': vars.prefix = w[1]
    if len(w)>=2 and w[0]=='barseq_catalog_in': vars.barseq_catalog_in = w[1]
    if len(w)>=2 and w[0]=='barseq_catalog_out': vars.barseq_catalog_out = w[1]


def save_config(vars):
  f = open("tpp.cfg","w")
  f.write("reads1 %s\n" % vars.fq1)
  f.write("reads2 %s\n" % vars.fq2)
  f.write("ref %s\n" % ' '.join(vars.ref))
  f.write("ids %s\n" % ','.join(vars.replicon_ids))
  f.write("bwa %s\n" % vars.bwa)
  f.write("bwa_alg %s\n" % vars.bwa_alg)
  f.write("flags %s\n" % vars.flags)
  f.write("mismatches1 %s\n" % vars.mm1)
  f.write("primer_start_window %s,%s\n" % (vars.primer_start_window[0],vars.primer_start_window[1]))
  f.write("window_size %s\n" % vars.window_size)
  if vars.maxreads>-1: f.write("maxreads %s\n" % vars.maxreads)
  f.write("transposon %s\n" % vars.transposon)
  f.write("protocol %s\n" % vars.protocol)
  f.write("primer %s\n" % vars.prefix)
  if vars.barseq_catalog_in!=None: f.write("barseq_catalog_in %s\n" % vars.barseq_catalog_in)
  if vars.barseq_catalog_out!=None: f.write("barseq_catalog_out %s\n" % vars.barseq_catalog_out)
  f.close()

def show_help():
  #print('usage: python PATH/src/tpp.py -bwa <EXECUTABLE_WITH_PATH> -ref <fasta-file|comma_separated_list> -reads1 <FASTQ_OR_FASTA_FILE> [-reads2 <FASTQ_OR_FASTA_FILE>] -output <BASE_FILENAME> [-maxreads <N>] [-mismatches <N>] [-flags "<STRING>"] [-tn5|-himar1] [-primer <seq>] [-primer-start-window INT,INT] [-window-size INT] [-barseq_catalog_in|_out <file>] [-replicon-ids <comma_separated_list_of_names>]')

  print('usage: python PATH/src/tpp.py -bwa <EXECUTABLE_WITH_PATH> -ref <fasta-file|comma_separated_list> -reads1 <FASTQ_OR_FASTA_FILE> [-reads2 <FASTQ_OR_FASTA_FILE>] -output <BASE_FILENAME> [OPTIONAL ARGS]')
  print('  OPTIONAL ARGS:')
  print('    -protocol [Sassetti|Tn5|Mme1] # which sample prep protocol was used?; sassetti protocol is the default; this sets the default transposon and primer sequence')
  print('    -primer <seq>      # prefix of reads corresponding to end of transposon at junction with genomic sequence; can override default seq' )
  print('    -maxreads <INT>')
  print('    -mismatches <INT>  # when searching for constant regions in reads 1 and 2; default is 1')
  print('    -flags "<STRING>"  # args to pass to BWA')
  print('    -bwa-alg [aln|mem]  # Algorithm to use for mapping reads with bwa; default is \'aln\'' )
  print('    -primer-start-window INT,INT # position in read to search for start of primer; default is: [0,20]')
  print('    -window-size INT   # automatic method to set window')
  #print('    -barseq_catalog_in|-barseq_catalog_out <file>')
  print('    -replicon-ids <comma_separated_list_of_names> # if multiple replicons/genomes/contigs/sequences were provided in -ref, give them names.')
  print('                                                  # Enter \'auto\' for autogenerated ids.')

class Globals:
  pass
