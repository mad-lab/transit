import wx
import wx.xrc
import wx.lib.mixins.listctrl as listmix
import ntpath

########################################################################
class SortableListCtrl(wx.ListCtrl):

    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)




class FileFrame(wx.Frame, listmix.ColumnSorterMixin):
    #constructor
    def __init__(self, parent, filePath, method="Gumbel"):

        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "%s - %s" % (method, ntpath.basename(filePath)), pos = wx.DefaultPosition, size = wx.Size( 1019,740 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        bSizer1 = wx.BoxSizer( wx.VERTICAL )

        sbSizer1 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Stats" ), wx.HORIZONTAL )

        self.headerText1 = wx.StaticText( self, wx.ID_ANY, u" ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.headerText1.Wrap( -1 )
        sbSizer1.Add( self.headerText1, 0, wx.ALL, 5 )


        self.headerText2 = wx.StaticText( self, wx.ID_ANY, u" ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.headerText2.Wrap( -1 )
        sbSizer1.Add( self.headerText2, 0, wx.ALL, 5 )


        self.headerText3 = wx.StaticText( self, wx.ID_ANY, u" ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.headerText3.Wrap( -1 )
        sbSizer1.Add( self.headerText3, 0, wx.ALL, 5 )



        bSizer1.Add( sbSizer1, 0, wx.EXPAND, 5 )

        bSizer3 = wx.BoxSizer( wx.VERTICAL )


        self.index_data = 0
        if method == "Gumbel":
            self.list_data = SortableListCtrl(self, size=(-1,100),
                         style=wx.LC_REPORT
                         |wx.BORDER_SUNKEN
                         |wx.LC_SORT_ASCENDING
                         )
            self.initializeGumbel()
            self.populateGumbel(filePath)
            listmix.ColumnSorterMixin.__init__(self, len(self.itemDataMap[0]))
            self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list_data)

        elif method == "HMM - Sites":
            self.list_data = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
            self.initializeHMMSites()
            self.populateHMMSites(filePath)

        elif method == "HMM - Genes":
            self.list_data = SortableListCtrl(self, size=(-1,100),
                         style=wx.LC_REPORT
                         |wx.BORDER_SUNKEN
                         |wx.LC_SORT_ASCENDING
                         )
            self.initializeHMMGenes()
            self.populateHMMGenes(filePath)
            listmix.ColumnSorterMixin.__init__(self, len(self.itemDataMap[0]))
            self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list_data)

        elif method == "Resampling":
            try:
                self.list_data = SortableListCtrl(self, size=(-1,100),
                         style=wx.LC_REPORT
                         |wx.BORDER_SUNKEN
                         |wx.LC_SORT_ASCENDING
                         )
                self.initializeResampling()
                self.populateResampling(filePath)
                listmix.ColumnSorterMixin.__init__(self, len(self.itemDataMap[0]))
                self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list_data)
            except Exception as e:
                print "Error", e



        #self.list_data = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        #self.list_data = SortableListCtrl( self) #, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.SUNKEN_BORDER )
        #self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list_data)

        bSizer3.Add( self.list_data, 1, wx.ALL|wx.EXPAND, 5 )



        bSizer1.Add( bSizer3, 1, wx.EXPAND, 5 )


        self.SetSizer( bSizer1 )
        self.Layout()

        self.Centre( wx.BOTH )
 
        """
        self.index_data = 0
        if method == "Gumbel":
            self.initializeGumbel()
            self.populateGumbel(filePath)

        elif method == "HMM - Sites":
            self.initializeHMMSites()
            self.populateHMMSites(filePath)

        elif method == "HMM - Genes":
            self.initializeHMMGenes()
            self.populateHMMGenes(filePath)
            
        elif method == "Resampling":
            self.initializeResampling()
            self.populateResampling(filePath)
        """

            

    def initializeGumbel(self):
        self.list_data.InsertColumn(0, 'Orf', width=100)
        self.list_data.InsertColumn(1, 'Name', width=85)
        self.list_data.InsertColumn(2, 'Description', width=220)
        self.list_data.InsertColumn(3, 'k', width=75)
        self.list_data.InsertColumn(4, 'n', width=75)
        self.list_data.InsertColumn(5, 'r', width=75)
        self.list_data.InsertColumn(6, 's', width=75)
        self.list_data.InsertColumn(7, 'zbar', width=100)
        self.list_data.InsertColumn(8, 'Call', width=50)
    

    def populateGumbel(self, path):

        self.itemDataMap = {}
        ess=0; unc=0; non=0; short=0
        
        data = []
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            data.append([tmp[0].upper()]+tmp)
            
        data.sort()
        
        for tmp in data:
            if not tmp: continue
            if len(tmp) < 3: continue
            self.list_data.InsertStringItem(self.index_data, tmp[1])
            for i,cell in enumerate(tmp[2:]):
                self.list_data.SetStringItem(self.index_data, i+1, cell)
            self.list_data.SetItemData(self.index_data, self.index_data)
            self.itemDataMap[self.index_data] = tmp[1:]
            self.index_data+=1
            if tmp[8] == "E": ess+=1
            if tmp[8] == "U": unc+=1
            if tmp[8] == "NE": non+=1
            if tmp[8] == "S": short+=1

        text = """Results:
    Essentials: %s
    Uncertain: %s
    Non-Essential: %s
    Short: %s
        """ % (ess, unc, non, short)
        self.headerText1.SetLabel(text)

    def initializeHMMSites(self):
        self.list_data.InsertColumn(0, 'Location', width=100)
        self.list_data.InsertColumn(1, 'Read', width=85)
        self.list_data.InsertColumn(2, 'Likelihood - ES', width=120)
        self.list_data.InsertColumn(3, 'Likelihood - GD', width=120)
        self.list_data.InsertColumn(4, 'Likelihood - NE', width=120)
        self.list_data.InsertColumn(5, 'Likelihood - GA', width=120)
        self.list_data.InsertColumn(6, 'State', width=75)


    def populateHMMSites(self, path):
        #self.itemDataMap = {}
        T = 0; es=0; gd=0; ne=0; ga=0;
        for line in open(path, "r"):
            if line.startswith("#"): continue
            tmp = line.strip().split()
            if not tmp: continue

            n = len(tmp)
            if n < 3:
                tmp + [""] * (7-n)
            self.list_data.InsertStringItem(self.index_data, tmp[0])
            for i,cell in enumerate(tmp[1:]):
                self.list_data.SetStringItem(self.index_data, i+1, cell)
            #self.list_data.SetItemData(self.index_data, self.index_data)
            #self.itemDataMap[self.index_data] = tmp
            self.index_data+=1

            if len(tmp) < 5: continue
            if tmp[-1] == "ES": es+=1
            if tmp[-1] == "GD": gd+=1
            if tmp[-1] == "NE": ne+=1
            if tmp[-1] == "GA": ga+=1
            T+=1
            text = """Results:
    Essential: %1.1f%%
    Growth-Defect: %1.1f%%
    Non-Essential: %1.1f%%
    Growth-Advantage: %1.1f%%
        """ % (100.0*es/T, 100.0*gd/T, 100.0*ne/T, 100.0*ga/T)
        self.headerText1.SetLabel(text)
    


    def initializeHMMGenes(self):
        self.list_data.InsertColumn(0, 'Orf', width=100)
        self.list_data.InsertColumn(1, 'Name', width=85)
        self.list_data.InsertColumn(2, 'Description', width=140)
        self.list_data.InsertColumn(3, 'N', width=75)
        self.list_data.InsertColumn(4, 'n0', width=75)
        self.list_data.InsertColumn(5, 'n1', width=75)
        self.list_data.InsertColumn(6, 'n2', width=75)
        self.list_data.InsertColumn(7, 'n3', width=75)
        self.list_data.InsertColumn(8, 'Avg. Insertions', width=75)
        self.list_data.InsertColumn(9, 'Avg. Reads', width=75)
        self.list_data.InsertColumn(10, 'State Call', width=75)
        

    def populateHMMGenes(self, path):

        self.itemDataMap = {}
        es=0; gd=0; ne=0; ga=0;
        data=[]
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            data.append([tmp[0].upper()] + tmp)
            
        data.sort()
        for tmp in data:
            self.list_data.InsertStringItem(self.index_data, tmp[1])
            for i,cell in enumerate(tmp[2:]):
                self.list_data.SetStringItem(self.index_data, i+1, cell)

            self.list_data.SetItemData(self.index_data, self.index_data)
            self.itemDataMap[self.index_data] = tmp[1:]
            self.index_data+=1
            if len(tmp) < 5: continue
            if tmp[-1] == "ES": es+=1
            if tmp[-1] == "GD": gd+=1
            if tmp[-1] == "NE": ne+=1
            if tmp[-1] == "GA": ga+=1

        text = """Results:
    Essential: %s
    Growth-Defect: %s
    Non-Essential: %s
    Growth-Advantage: %s
        """ % (es, gd, ne, ga)
        self.headerText1.SetLabel(text)


    def initializeResampling(self):
        self.list_data.InsertColumn(0, 'Orf', width=100)
        self.list_data.InsertColumn(1, 'Name', width=85)
        self.list_data.InsertColumn(2, 'Description', width=140)
        self.list_data.InsertColumn(3, 'N', width=85)
        self.list_data.InsertColumn(4, 'TAs Hit', width=100)
        self.list_data.InsertColumn(5, 'Avg. Read 1', width=100)
        self.list_data.InsertColumn(6, 'Avg. Read 2', width=100)
        self.list_data.InsertColumn(7, 'Delta Read', width=100)
        self.list_data.InsertColumn(8, 'p-value', width=75)
        self.list_data.InsertColumn(9, 'q-value', width=75)


    def populateResampling(self, path):
        self.itemDataMap = {}
        de05=0; de01 = 0; count = 0;
        data=[]
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            data.append([tmp[0].upper()] + tmp)
            
        data.sort()
        
        for tmp in data:
            self.list_data.InsertStringItem(self.index_data, tmp[1])
            for i,cell in enumerate(tmp[2:]):
                self.list_data.SetStringItem(self.index_data, i+1, cell)
            self.list_data.SetItemData(self.index_data, self.index_data)
            self.itemDataMap[self.index_data] = tmp[1:]
            self.index_data+=1
            if float(tmp[-1]) < 0.05: de05+=1
            if float(tmp[-1]) < 0.01: de01+=1
            count +=1
        
        text = """Results:
    Significant Hits (q<0.05): %s
    Significant Hits (q<0.01): %s
        """ % (de05, de01)
        self.headerText1.SetLabel(text)
        text2 = """         Notes:
            TAs Hit:   Number of TA sites with insertions combined across conditions.
            Avg Reads: Average Read count, normalized with the Non-Zero Mean normalization method. """
        self.headerText2.SetLabel(text2)


#----------------------------------------------------------------------
    # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
    def GetListCtrl(self):
        return self.list_data
 
    #----------------------------------------------------------------------
    def OnColClick(self, event):
        #print "column clicked self:", self
        #print "column clicked event:", event
        event.Skip()
 
