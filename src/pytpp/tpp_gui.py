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

import sys
import glob
import os
import time
import math
import re
import shutil
import platform
import gzip

try:
    import wx
    import wx.lib.filebrowsebutton
    hasWx = True
except Exception as e:
    hasWx = False

from tpp_tools import *

if hasWx:
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
            label0 = wx.StaticText(panel, label='BWA executable:',size=(340,-1))
            sizer0.Add(label0,0,0,0)
            print os.path.dirname(vars.bwa)
            #self.picker0 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="path to BWA",size=(400,30))#,path=os.path.abspath(vars.bwa))
            #self.picker0.SetDirName('/pacific/home/cambadipudi/chaitra/tpp/')
            self.picker0 = wx.lib.filebrowsebutton.FileBrowseButton(panel, id = wx.ID_ANY, size=(400,30), dialogTitle='Path to BWA', fileMode=wx.OPEN, fileMask='bwa*', startDirectory=os.path.dirname(vars.bwa), initialValue=vars.bwa, labelText='')
            sizer0.Add(self.picker0, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer0,0,wx.EXPAND,0)

            sizer3 = wx.BoxSizer(wx.HORIZONTAL)
            label3 = wx.StaticText(panel, label='Choose a reference genome (FASTA):',size=(340,-1))
            sizer3.Add(label3,0,0,0)
            #self.picker3 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the reference genome", wildcard='*.fna;*.fasta;*.fa', size=(400,30),path=vars.ref)
            self.picker3 = wx.lib.filebrowsebutton.FileBrowseButton(panel, id=wx.ID_ANY, dialogTitle='Please select the reference genome', fileMode=wx.OPEN, fileMask='*.fna;*.fasta;*.fa', size=(400,30), startDirectory=os.path.dirname(vars.ref), initialValue=vars.ref, labelText='')
            sizer3.Add(self.picker3, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer3,0,wx.EXPAND,0)
       
        
            sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            label1 = wx.StaticText(panel, label='Choose the Fastq file for read 1:',size=(340,-1))
            sizer1.Add(label1,0,0,0)
            # self.picker1 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the .fastq file for read 1", wildcard='*.fastq;*.fq;*.reads;*.fasta;*.fa', size=(400,30),path=vars.fq1)
            self.picker1 = wx.lib.filebrowsebutton.FileBrowseButton(panel, id=wx.ID_ANY, dialogTitle='Please select the .fastq file for read 1', fileMode=wx.OPEN, fileMask='*.fastq;*.fq;*.reads;*.fasta;*.fa;*.fastq.gz', size=(400,30), startDirectory=os.path.dirname(vars.fq1), initialValue=vars.fq1, labelText='',changeCallback=self.OnChanged2)
            #self.picker1.OnChanged = self.OnChanged(self.picker1.GetValue(), self.base)
            #self.Bind(wx.EVT_TEXT, self.OnChanged, id=self.picker1.GetId())
            #self.picker1.Bind(wx.EVT_TEXT, self.OnChanged)
            sizer1.Add(self.picker1, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer1,0,wx.EXPAND,0)
       
            sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            label2 = wx.StaticText(panel, label='Choose the Fastq file for read 2:',size=(340,-1))
            sizer2.Add(label2,0,0,0)
            #self.picker2 = wx.FilePickerCtrl(panel, wx.ID_ANY,message="Please select the .fastq file for read 2", wildcard='*.fastq;*.fq;*.reads;*.fasta;*.fa', size=(400,30),path=vars.fq2)
            self.picker2 = wx.lib.filebrowsebutton.FileBrowseButton(panel, id=wx.ID_ANY, dialogTitle='Please select the .fastq file for read 2', fileMode=wx.OPEN, fileMask='*.fastq;*.fq;*.reads;*.fasta;*.fa;*.fastq.gz', size=(400,30), startDirectory=os.path.dirname(vars.fq2), initialValue=vars.fq2, labelText='', changeCallback=self.OnChanged2)
            sizer2.Add(self.picker2, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer2,0,wx.EXPAND,0)

            sizer4 = wx.BoxSizer(wx.HORIZONTAL)
            label4 = wx.StaticText(panel, label='Prefix to use for output filenames:',size=(350,-1))
            sizer4.Add(label4,0,0,0)
            self.base = wx.TextCtrl(panel,value=vars.base,size=(400,30))
            sizer4.Add(self.base, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer4,0,wx.ALL,0)

            sizer7 = wx.BoxSizer(wx.HORIZONTAL)
            label7 = wx.StaticText(panel, label='Transposon used:',size=(350,-1))
            sizer7.Add(label7,0,0,0)
            self.transposon = wx.ComboBox(panel,choices=['Himar1','Tn5'],size=(400,30))
            self.transposon.SetStringSelection(vars.transposon)
            sizer7.Add(self.transposon, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer7,0,wx.ALL,0)    

            sizer5 = wx.BoxSizer(wx.HORIZONTAL)
            label5 = wx.StaticText(panel, label='Max reads (leave blank to use all):',size=(350,-1))
            sizer5.Add(label5,0,0,0)
            self.maxreads = wx.TextCtrl(panel,size=(400,30))
            sizer5.Add(self.maxreads, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer5,0,wx.ALL,0)

            sizer6 = wx.BoxSizer(wx.HORIZONTAL)
            label6 = wx.StaticText(panel, label='Mismatches allowed in Tn prefix:',size=(350,-1))
            sizer6.Add(label6,0,0,0)
            self.mismatches = wx.TextCtrl(panel,value=str(vars.mm1),size=(400,30))
            sizer6.Add(self.mismatches, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
            sizer.Add(sizer6,0,wx.ALL,0)    

            #self.picker1.OnChanged = self.OnChanged(self.picker1.GetValue())


        def OnChanged(self, str_path):
            print "changed"
            value = os.path.basename(str_path).split('.')[0]
            if '_R1' in value or '_R2':
                value = value.split('_')[0]
            self.base.SetValue(value)

        def OnChanged2(self, event):
            value2 = os.path.basename(self.picker2.GetValue()).split('.')[0]
            value1 = os.path.basename(self.picker1.GetValue()).split('.')[0]
            value = os.path.commonprefix([value1, value2])
            self.base.SetValue(value)
            self.base.Refresh()

        def InitList(self,panel,sizer):
            self.list_ctrl = wx.ListCtrl(panel, size=(500,500), style=wx.LC_HRULES|wx.LC_VRULES|wx.LC_REPORT|wx.BORDER_SUNKEN)
            self.list_ctrl.InsertColumn(0, 'Dataset (*.tn_stats)',width=300)
            self.list_ctrl.InsertColumn(1, 'total reads',wx.LIST_FORMAT_RIGHT,width=125)
            self.list_ctrl.InsertColumn(2, 'TGTTA prefix', wx.LIST_FORMAT_RIGHT,width=125)
            self.list_ctrl.InsertColumn(3, 'R1_mapped', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(4, 'R2_mapped', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(5, 'mapped\nreads', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(6, 'template\ncount', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(7, 'TAs hit', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(8, 'insertion\ndensity',wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(9, 'NZmean', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(10, 'maxcount', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(11, 'primer', wx.LIST_FORMAT_RIGHT,width=90)
            self.list_ctrl.InsertColumn(12, 'vector',wx.LIST_FORMAT_RIGHT,width=90)
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
          datasets = []
          for fname in glob.glob("*.tn_stats"):
            filedate = os.path.getmtime(fname)
            datasets.append((filedate,fname))
          datasets.sort(reverse=True)
          for (filedate,fname) in datasets:
            stats = self.read_stats_file(fname)
        
            vals = [stats.get("total_reads","?"),stats.get("TGTTA_reads","?"),stats.get("reads1_mapped", "?"),stats.get("reads2_mapped","?"),stats.get("mapped_reads","?"),stats.get("template_count","?"), stats.get("TAs_hit","?"), stats.get("density", "?"), stats.get("NZ_mean", "?"), stats.get("max_count", "?"), stats.get("primer_matches:","?"),stats.get("vector_matches:","?")]

            #dataset = "[%s] %s" % (time.strftime("%m/%d/%y",time.localtime(os.path.getmtime(fname))),fname[:fname.rfind('.')])
            dsname = "[%s] %s" % (time.strftime("%m/%d/%y",time.localtime(filedate)),fname[:fname.rfind('.')])
            self.add_data(dsname, vals)

        def read_stats_file(self,fname):
          stats = {}
          for line in open(fname):
            w = line.rstrip().split()
            val = ""
            if len(w)>2: val = w[2]
            stats[w[1]] = val
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
          # bwapath = self.picker0.GetPath()
          bwapath = self.picker0.GetValue()
          #fq1,fq2,ref,base,maxreads = self.picker1.GetPath(),self.picker2.GetPath(),self.picker3.GetPath(),self.base.GetValue(),self.maxreads.GetValue()
          fq1, fq2, ref, base, maxreads = self.picker1.GetValue(), self.picker2.GetValue(), self.picker3.GetValue(), self.base.GetValue(), self.maxreads.GetValue()
    
          mm1 = self.mismatches.GetValue()
          try: mm1 = int(mm1)
          except Exception: mm1 = 1

          self.vars.transposon = self.transposon.GetStringSelection()

          self.vars.bwa = bwapath
          self.vars.fq1 = fq1
          self.vars.fq2 = fq2
          self.vars.ref = ref
          self.vars.base = base  
          self.vars.mm1 = mm1
          if maxreads == '': self.vars.maxreads = -1
          else: self.vars.maxreads = int(maxreads)

          self.vars.action = "start"
          self.Close()
          return 0

