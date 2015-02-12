# Copyright 2015 by Thomas R. Ioerger
# Texas A&M University
# ioerger@cs.tamu.edu

import wx
import glob,os,sys,time,math
import sys, re, shutil
import platform

class MyForm(wx.Frame):
 
    def __init__(self,vars):
        self.vars = vars
        initialize_globals(self.vars)

        wx.Frame.__init__(self, None, wx.ID_ANY, "Tn-Seq PreProcessor") # v%s" % vars.version)
        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, wx.ID_ANY)
        #panel.SetFont(wx.Font(14,wx.DECORATIVE,wx.NORMAL,wx.NORMAL,False,u'times'))
        # panel.SetFont(wx.Font(14,wx.DECORATIVE,wx.NORMAL,wx.NORMAL,False,u'fixed'))
        #panel.SetFont(wx.Font(14,wx.DECORATIVE,wx.NORMAL,wx.BOLD,False,u'courier'))

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.list_ctrl = None    
        self.InitMenu()
        self.InitFiles(panel,sizer)

        buttonrow = wx.BoxSizer(wx.HORIZONTAL)

        btn = wx.Button(panel, label="Start")
        btn.Bind(wx.EVT_BUTTON, self.map_reads)
        buttonrow.Add(btn,0,0,0,10)

        btn = wx.Button(panel, label="Quit")
        btn.Bind(wx.EVT_BUTTON, self.OnQuit)
        buttonrow.Add(btn,0,0,0,10)
        sizer.Add(buttonrow,0,0,0)

        self.InitList(panel,sizer)

        panel.SetSizer(sizer)
        # self.SetSize((1305, 700))
        self.SetSize((900, 700))
        #self.SetTitle('Simple menu')
        self.Centre()
        #self.Show(True)

        self.pid = None

    '''
    def initialize_globals(self):
      vars = self.vars
      vars.fq1,vars.fq2,vars.ref,vars.bwa,vars.base,vars.maxreads = "","","","","temp",-1
      read_config(vars)
    '''

    def InitFiles(self,panel,sizer):
        vars = self.vars
        sizer0 = wx.BoxSizer(wx.HORIZONTAL)
        label0 = wx.StaticText(panel, label='BWA executable:',size=(350,-1))
        sizer0.Add(label0,0,0,0)
        print vars.bwa
        self.picker0 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="path to BWA",size=(400,30),path="/pacific/home/cambadipudi")#os.path.abspath(vars.bwa))
        sizer0.Add(self.picker0, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer0,0,wx.EXPAND,0)

        sizer3 = wx.BoxSizer(wx.HORIZONTAL)
        label3 = wx.StaticText(panel, label='Choose a reference genome (FASTA):',size=(350,-1))
        sizer3.Add(label3,0,0,0)
        self.picker3 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the reference genome", wildcard='*.fna;*.fasta;*.fa', size=(400,30),path=vars.ref)
        sizer3.Add(self.picker3, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer3,0,wx.EXPAND,0)
       
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        label1 = wx.StaticText(panel, label='Choose the Fastq file for read 1:',size=(350,-1))
        sizer1.Add(label1,0,0,0)
        self.picker1 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the .fastq file for read 1", wildcard='*.fastq;*.fq;*.reads;*.fasta;*.fa', size=(400,30),path=vars.fq1)
        sizer1.Add(self.picker1, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer1,0,wx.EXPAND,0)
       
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        label2 = wx.StaticText(panel, label='Choose the Fastq file for read 2:',size=(350,-1))
        sizer2.Add(label2,0,0,0)
        self.picker2 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the .fastq file for read 2", wildcard='*.fastq;*.fq;*.reads;*.fasta;*.fa', size=(400,30),path=vars.fq2)
        sizer2.Add(self.picker2, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer2,0,wx.EXPAND,0)

        sizer4 = wx.BoxSizer(wx.HORIZONTAL)
        label4 = wx.StaticText(panel, label='Prefix to use for filenames:',size=(350,-1))
        sizer4.Add(label4,0,0,0)
        self.base = wx.TextCtrl(panel,value=vars.base,size=(400,30))
        sizer4.Add(self.base, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer4,0,wx.ALL,0)

        sizer5 = wx.BoxSizer(wx.HORIZONTAL)
        label5 = wx.StaticText(panel, label='Max reads:',size=(350,-1))
        sizer5.Add(label5,0,0,0)
        self.maxreads = wx.TextCtrl(panel,size=(400,30))
        sizer5.Add(self.maxreads, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
        sizer.Add(sizer5,0,wx.ALL,0)


    def InitList(self,panel,sizer):
        self.list_ctrl = wx.ListCtrl(panel, size=(-1,-1), style=wx.LC_HRULES|wx.LC_VRULES|wx.LC_REPORT|wx.BORDER_SUNKEN)
        self.list_ctrl.InsertColumn(0, 'dataset',width=300)
        self.list_ctrl.InsertColumn(1, 'total reads',wx.LIST_FORMAT_RIGHT,width=125)
        self.list_ctrl.InsertColumn(2, 'TGTTA prefix', wx.LIST_FORMAT_RIGHT,width=125)
        self.list_ctrl.InsertColumn(3, 'mapped\nreads', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(4, 'template\ncount', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(5, 'TAs hit', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(6, 'insertion\ndensity',wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(7, 'NZmean', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(8, 'maxcount', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(9, 'primer', wx.LIST_FORMAT_RIGHT,width=90)
        self.list_ctrl.InsertColumn(10, 'vector',wx.LIST_FORMAT_RIGHT,width=90)
        #self.list_ctrl.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL))
        #btn = wx.Button(panel, label="Add Line")
        #btn.Bind(wx.EVT_BUTTON, self.add_line)
 
        sizer.Add(self.list_ctrl, 0, wx.ALL|wx.EXPAND, 10)
        #sizer.Add(btn, 0, wx.ALL|wx.CENTER, 5)

    def InitMenu(self):    
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()

        #dataset_menuitem = fileMenu.Append(wx.ID_ANY, 'Add New Dataset', 'Analyze New Dataset')
        #self.Bind(wx.EVT_MENU, self.addNewDataset, dataset_menuitem)

        quit_menuitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
        self.Bind(wx.EVT_MENU, self.OnQuit, quit_menuitem)

        menubar.Append(fileMenu, '&File')
        self.SetMenuBar(menubar)

    def addNewDataset(self, event):
      dlg = wx.FileDialog(
          self, message="Choose a file",
          defaultDir=".",
          defaultFile="",
          wildcard="*.wig",
          style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
          )
      if dlg.ShowModal() == wx.ID_OK:
          paths = dlg.GetPaths()
          for path in paths:
            print "analyzing dataset:",path
            analyze_dataset(path)
      dlg.Destroy()
      self.update_dataset_list()

    def update_dataset_list(self):
      if self.list_ctrl==None: return
      self.list_ctrl.DeleteAllItems()
      self.index = 0
      for fname in glob.glob("*.tn_stats"):
        stats = self.read_stats_file(fname)
        
        #self.add_data(fname,

        vals = [stats.get("total_reads","?"),stats.get("TGTTA_reads","?"),stats.get("mapped_reads","?"),stats.get("template_count","?"), stats.get("TAs_hit","?"), stats.get("density", "?"), stats.get("NZ_mean", "?"), stats.get("max_count", "?"), stats.get("primer_matches:","?"),stats.get("vector_matches:","?")]

        self.add_data(fname, vals)

    def read_stats_file(self,fname):
      stats = {}
      for line in open(fname):
        w = line.rstrip().split()
        stats[w[1]] = w[2]
      return stats

    def add_data(self, dataset,vals):
        self.list_ctrl.InsertStringItem(self.index, dataset)
        for i in range(1, len(vals)+1):
            self.list_ctrl.SetStringItem(self.index, i, vals[i-1])
        self.index += 1

    def OnQuit(self, e):
        print "Quitting TPP.  Good bye."
        self.vars.action = "quit"
        self.Close()
        return 0

    def map_reads(self,event):
      # add bwa path, prefix
      bwapath = self.picker0.GetPath()
      fq1,fq2,ref,base,maxreads = self.picker1.GetPath(),self.picker2.GetPath(),self.picker3.GetPath(),self.base.GetValue(),self.maxreads.GetValue()

      self.vars.bwa = bwapath
      self.vars.fq1 = fq1
      self.vars.fq2 = fq2
      self.vars.ref = ref
      self.vars.base = base  
      if maxreads == '': self.vars.maxreads = -1
      else: self.vars.maxreads = int(maxreads)

      self.vars.action = "start"
      self.Close()
      return 0

# http://www.blog.pythonlibrary.org/2013/09/04/wxpython-how-to-update-a-progress-bar-from-a-thread/\

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
    if maxreads > -1:
        if tot > maxreads:
            break
    if cnt==0: 
      h = line[1:] # strip off '@'
      #h = h.replace(' ','_')
      output.write(">%s" % h)
    if cnt==1: output.write(line)
    cnt = (cnt+1)%4
  print cnt,tot
  output.close()

# the headers for each pair must be identical up to /1 and /2 at the ends
# if the variable character with the read number occurs in the middle, move it to the end

def fix_paired_headers_for_bwa(reads1,reads2):
  a = open(reads1)
  b = open(reads2)
  temp1 = "temp."+reads1
  temp2 = "temp."+reads2
  c = open(temp1,"w")
  d = open(temp2,"w")
  try:
   while True:
    e = a.readline().rstrip()
    f = b.readline().rstrip()
    if len(e)<=2 or len(f)<=2: break
    if e[0]=='>':
      # find first position where there is a difference
      i,n = 0,len(e)
      if len(f)!=n: raise Exception('unexpected format of headers in .fastq files')
      while i<n and e[i]==f[i]: i += 1
      if i==n: raise Exception('unexpected format of headers in .fastq files')
      e = e.replace(' ','_')
      f = f.replace(' ','_')
      e = e.replace('/','_') # this was neceesary for bwa 0.7.10 but not 0.7.12
      f = f.replace('/','_')
      #if i<n-1:
      if e[i+1:]!=f[i+1:]: raise Exception('unexpected format of headers in .fastq files')
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

def mmfind(G,n,H,m): # lengths; assume n>m
  for i in range(0,n-m):
    cnt,max = 0,1
    for k in range(m):
      if G[i+k]!=H[k]: cnt += 1
      if cnt>max: break
    if cnt<=max: return i
  return -1

def extract_staggered(infile,outfile):
  Tn = "ACTTATCAGCCAACCTGTTA"
  lenTn = len(Tn)

  LEN = -1
  P,Q = 5,10 # 1-based inclusive positions to look for start of Tn prefix

  n,m = -1,len(Tn)
  output = open(outfile,"w")
  for line in open(infile):
    if line[0]=='>': h = line; continue
    elif n==-1: n = len(line.rstrip()) # readlen
    if LEN==-1: LEN = n-lenTn-(Q+1) # genomic suffix len
    #a = line.find(Tn)
    a = mmfind(line,n,Tn,m)
    if a>=P and a<=Q:
      w = line[a+lenTn:a+lenTn+LEN]
      output.write(h)
      output.write(w+"\n")
  output.close()

def message(s):
  print "[tn_preprocess]",s
  sys.stdout.flush()

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

def template_counts(ref,sam,bcfile):
  genome = read_genome(ref)

  barcodes = {}
  for line in open(bcfile):
    line = line.rstrip()
    if line[0]=='>': id = line[1:]
    else: barcodes[id] = line

  hits = {}
  readlen = -1
  tot,mapped = 0,0
  for line in open(sam):
    if line[0]=='@': continue
    else:
      w = line.split('\t')
      if int(w[1])>=128: tot += 1 # just count read1's
      if w[1]=="99" or w[1]=="83":
        mapped += 1
        if readlen==-1: readlen = len(w[9])
        pos,size = int(w[3]),int(w[8])
        strand,delta = 'F',-2
        if w[1]=="83": strand,delta = 'R',readlen
        pos += delta
        bc = barcodes[w[0]]
        if pos not in hits: hits[pos] = []
        hits[pos].append((strand,size,bc))

  sites = []
  for i in range(len(genome)-1):
    if genome[i:i+2]=="TA":
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

def driver(vars):
  vars.reads1 = vars.base+".reads1"
  vars.reads2 = vars.base+".reads2"
  vars.tgtta1 = vars.base+".tgtta1"
  vars.tgtta2 = vars.base+".tgtta2"
  vars.barcodes1 = vars.base+".barcodes1"
  vars.barcodes2 = vars.base+".barcodes2"
  vars.genomic2 = vars.base+".genomic2"
  vars.sai1 = vars.base+".sai1"
  vars.sai2 = vars.base+".sai2"
  vars.sam = vars.base+".sam"
  vars.tc = vars.base+".counts"
  vars.wig = vars.base+".wig"
  vars.stats = vars.base+".tn_stats"

  extract_reads(vars)

  run_bwa(vars)

  generate_output(vars)

  message("Done.")


def extract_reads(vars):
    message("extracting reads...")
    
    flag = ['','']
    for idx, name in enumerate([vars.fq1, vars.fq2]):
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
        shutil.copyfile(vars.fq1, vars.reads1)
    if(flag[1] == 'FASTQ'): 
        message("fastq2reads: %s -> %s" % (vars.fq2,vars.reads2))
        fastq2reads(vars.fq2,vars.reads2,vars.maxreads)
    else:
        shutil.copyfile(vars.fq2, vars.reads2)

    message("fixing headers of paired reads for bwa...")
    fix_paired_headers_for_bwa(vars.reads1,vars.reads2)

    message("extracting barcodes and genomic parts of reads...")

    message("creating %s" % vars.tgtta1)
    extract_staggered(vars.reads1,vars.tgtta1)
    message("creating %s" % vars.tgtta2)
    select_reads(vars.tgtta1,vars.reads2,vars.tgtta2)
    message("creating %s" % vars.barcodes2)
    select_cycles(vars.tgtta2,22,30,vars.barcodes2)
    message("creating %s" % vars.barcodes1)
    replace_ids(vars.tgtta1,vars.barcodes2,vars.barcodes1)
    message("creating %s" % vars.genomic2)
    select_cycles(vars.tgtta2,43,-1,vars.genomic2)


def run_bwa(vars):
    message("mapping reads using BWA...(this takes a couple of minutes)")

    if not os.path.exists(vars.ref+".amb"):
      cmd = vars.bwa+" index "+vars.ref
      message(cmd)
      os.system(cmd)

    cmd = "%s aln %s %s > %s" % (vars.bwa,vars.ref,vars.tgtta1,vars.sai1)
    message(cmd)
    os.system(cmd)

    cmd = "%s aln %s %s > %s" % (vars.bwa,vars.ref,vars.genomic2,vars.sai2)
    message(cmd)
    os.system(cmd)

    cmd = "%s sampe %s %s %s %s %s > %s" % (vars.bwa,vars.ref,vars.sai1,vars.sai2,vars.tgtta1,vars.genomic2,vars.sam)
    message(cmd)
    os.system(cmd)


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
  cX = [x-muX for x in X]
  cY = [y-muY for y in Y]
  s = sum([x*y for (x,y) in zip(cX,cY)])
  return s/(float(len(X))*sdX*sdY)


def generate_output(vars):
  message("tabulating template counts and statistics...")
  counts = template_counts(vars.ref,vars.sam,vars.barcodes1)

  tcfile = open(vars.tc,"w")
  tcfile.write('\t'.join("coord Fwd_Rd_Ct Fwd_Templ_Ct Rev_Rd_Ct Rev_Templ_Ct Tot_Rd_Ct Tot_Templ_Ct".split())+"\n")
  for data in counts: tcfile.write('\t'.join([str(x) for x in data])+"\n")
  tcfile.close()

  message("writing %s" % vars.wig)
  output = open(vars.wig,"w")

  read1 = os.path.basename(vars.fq1)
  read2 = os.path.basename(vars.fq2)
  fi = re.split(r'\.', os.path.basename(vars.ref))[0]
  output.write("# Generated by tpp from " + read1 + " and " + read2 + "\n")
  output.write("variableStep chrom="+ fi + "\n")
  for data in counts: output.write("%s %s\n" % (data[0],data[-1]))
  output.close()

  tot_tgtta,mapped = 0,0
  r1,r2=0,0
  for line in open(vars.sam):
    if line[0]=='@': continue
    else:
      w = line.split('\t')
      if int(w[1])>=128: tot_tgtta += 1 # just count read1's
      if w[1]=="99" or w[1]=="83":
        mapped += 1
      bin_str = bin(int(w[1]))[2:].zfill(8)[::-1]
      if bin_str[2]=='0' and bin_str[6]=='1': r1 +=1
      if bin_str[2]=='0' and bin_str[7]=='1': r2 +=1
      

  primer = "CTAGAGGGCCCAATTCGCCCTATAGTGAGT"
  vector = "CTAGACCGTCCAGTCTGGCAGGCCGGAAAC"

  tot_reads,nprimer,nvector = 0,0,0
  prefixes = {}
  for line in open(vars.reads1):
    if line[0]=='>': tot_reads += 1; continue
    if primer in line: nprimer += 1
    if vector in line: nvector += 1
    prefix = line[:30]
    if prefix not in prefixes: prefixes[prefix] = 0
    prefixes[prefix] += 1
  temp = prefixes.items()
  temp.sort(key=lambda x: x[1],reverse=True)

  rcounts = [x[5] for x in counts]
  tcounts = [x[6] for x in counts]
  rc,tc = sum(rcounts),sum(tcounts)
  ratio = rc/float(tc)
  ta_sites = len(rcounts)
  tas_hit = len(filter(lambda x: x>0,rcounts))
  density = tas_hit/float(ta_sites)
  counts.sort(key=lambda x: x[-1])
  max_tc = counts[-1][5]
  max_coord = counts[-1][0]
  NZmean = tc/float(tas_hit)
  FR_corr = corr([x[1] for x in counts],[x[3] for x in counts])
  BC_corr = corr([x for x in rcounts if x!=0],[x for x in tcounts if x!=0])

  output = open(vars.stats,"w")
  version = "1.0"
  output.write("# title: Tn-Seq Pre-Processor, version %s\n" % vars.version)
  output.write("# date: %s\n" % time.strftime("%m/%d/%Y %H:%M:%S"))
  output.write("# command: python ")
  output.write(' '.join(sys.argv)+"\n")
  output.write('# read1: %s\n' % vars.fq1)
  output.write('# read2: %s\n' % vars.fq2)
  output.write('# ref_genome: %s\n' % vars.ref)
  output.write("# total_reads %s (read pairs)\n" % tot_reads)
  output.write("# TGTTA_reads %s (reads with valid Tn prefix)\n" % tot_tgtta)
  output.write("# mapped_reads %s (both R1 and R2 map into genome)\n" % mapped)
  output.write("# reads1_mapped %s\n" % r1)
  output.write("# reads2_mapped %s\n" %r2)

  output.write("# read_count %s (TA sites only)\n" % rc)
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

  output.write("# primer_matches: %s reads contain %s\n" % (nprimer,primer))
  output.write("# vector_matches: %s reads contain %s\n" % (nvector,vector))
  #output.write("# most_abundant_prefix: %s reads start with %s\n" % (temp[0][1],temp[0][0]))
  # since these are reads (within Tn prefix stripped off), I expect ~1/4 to match Tn prefix
  vals = [vars.fq1,vars.fq2,tot_reads,tot_tgtta,mapped,r1,r2,rc,tc,ratio,ta_sites,tas_hit,max_tc,max_coord,NZmean,FR_corr,BC_corr,nprimer,nvector]
  output.write('\t'.join([str(x) for x in vals]))
  output.close()

  message("writing %s" % vars.stats)
  os.system("grep '#' %s" % vars.stats)

#############################################################################

def error(s):
  print "error:",s
  sys.exit(0)

def verify_inputs(vars):
  if not os.path.exists(vars.fq1): error("file not found: "+vars.fq1)
  if not os.path.exists(vars.fq2): error("file not found: "+vars.fq2)
  if not os.path.exists(vars.ref): error("file not found: "+vars.ref)

def initialize_globals(vars):
      vars.fq1,vars.fq2,vars.ref,vars.bwa,vars.base,vars.maxreads = "","","","","temp",-1
      read_config(vars)

def read_config(vars):
  if not os.path.exists("tpp.cfg"): return
  for line in open("tpp.cfg"):
    w = line.split()
    if len(w)>=2 and w[0]=='reads1': vars.fq1 = w[1]
    if len(w)>=2 and w[0]=='reads2': vars.fq2 = w[1]
    if len(w)>=2 and w[0]=='ref': vars.ref = w[1]
    if len(w)>=2 and w[0]=='bwa': vars.bwa = w[1]
    if len(w)>=2 and w[0]=='prefix': vars.base = w[1]

def save_config(vars):
  f = open("tpp.cfg","w")
  f.write("reads1 %s\n" % vars.fq1)
  f.write("reads2 %s\n" % vars.fq2)
  f.write("ref %s\n" % vars.ref)
  f.write("bwa %s\n" % vars.bwa)
  f.write("prefix %s\n" % vars.base)
  f.close()
    
class Globals:
  pass

if __name__ == "__main__":
    # if -nowin is command-line arg, skip the GUI and set filenames in vars
    
    vars = Globals()
    #vars.version = "$Revision: 1.5 $".split()[1]
    
    if(len(sys.argv) <= 1):
        app = wx.App(False)
        form = MyForm(vars)
        form.update_dataset_list()

        form.Show()
        app.MainLoop()

        # vars.action not defined, quit...

        if vars.action=="start":
            print "running pre-processing on %s and %s" % (vars.fq1,vars.fq2)
            verify_inputs(vars)
            save_config(vars)
            driver(vars)

    else:
        initialize_globals(vars)
        for i in range(0, len(sys.argv)):
            if sys.argv[i] == '-reads1': 
                vars.fq1 = sys.argv[i+1]
            elif sys.argv[i] == '-reads2':
                vars.fq2 = sys.argv[i+1]
            elif sys.argv[i] == '-bwa':
                vars.bwa = sys.argv[i+1]
            elif sys.argv[i] == '-ref':
                vars.ref = sys.argv[i+1]
            elif sys.argv[i] == '-maxreads':
                vars.maxreads = sys.argv[i+1]
            elif sys.argv[i] == '-prefix':
                vars.base = sys.argv[i+1]
        print 'running pre-processing on %s and %s' % (vars.fq1, vars.fq2)
        verify_inputs(vars)
        save_config(vars)
        driver(vars)
