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

try:
    import wx
    import wx.xrc
    import wx.lib.mixins.listctrl as listmix
    hasWx = True
except Exception as e:
    hasWx = False


import ntpath
import subprocess
import os
import sys

from functools import partial

import trash

import wx.grid

import pytransit
import pytransit.analysis


########################################################################
class SortableListCtrl(wx.ListCtrl):

    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)



#########################################
#menu_titles = [ "Display Histogram",
#                "Display Tracks",]
#
#menu_title_by_id = {}
#for title in menu_titles:
#    menu_title_by_id[ wx.NewId() ] = title
#
##########################################


class ImgFrame(wx.Frame):

    def __init__(self, parent, filePath):
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "%s" % (filePath), pos = wx.DefaultPosition, size = wx.Size( 1150,740 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        bSizer1 = wx.BoxSizer( wx.VERTICAL )
       
        self.m_bitmap1 = wx.StaticBitmap( self, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer1.Add( self.m_bitmap1, 1, wx.ALL|wx.EXPAND, 5 )
        self.SetSizer( bSizer1 )
        self.Layout()
        self.Centre( wx.BOTH )

        img = wx.Image(filePath, wx.BITMAP_TYPE_ANY)
        self.m_bitmap1.SetBitmap(wx.BitmapFromImage(img))
        
        self.Refresh()
        self.Fit()
 




def unknownColNames(path):
    colNames = []
    final_line = ""
    for line in open(path):
        if line.startswith("#"):
            colNames = line.split("\t")
        else:
            final_line = line
    if not final_line:
        print "Error: file appears to be empty"
    tmp = final_line.split("\t")
    if len(colNames) < len(tmp):
        colNames = ["Col%d" % (i) for i in range(len(tmp))]
    return colNames


def unknownTableData(path, colnames):
    row = 0
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        tmp[-1] = tmp[-1].strip()
        rowdict = dict([(colnames[i], tmp[i]) for i in range(len(colnames))])
        data.append((row, rowdict))
        row+=1
    return data


def unknownFileHeaderText(path):
    return "Unknown results file."


def getInfoFromFileType(X):
    for method in pytransit.analysis.methods:
        for filetype in pytransit.analysis.methods[method].filetypes:
            FT = filetype()
            if X == FT.identifier:
                return (method, FT)

    return ("unknown", pytransit.analysis.base.TransitFile())

    

class TransitTable(wx.grid.GridTableBase):
    """
    A custom wx.Grid Table using user supplied data
    """
    def __init__(self, data, colnames):
        """data is a list of the form
        [(rowname, dictionary),
        dictionary.get(colname, None) returns the data for column
        colname
        """
        # The base class must be initialized *first*
        wx.grid.GridTableBase.__init__(self)
        self.data = data
        self.colnames = colnames
        # XXX
        # we need to store the row length and column length to
        # see if the table has changed size
        self._rows = self.GetNumberRows()
        self._cols = self.GetNumberCols()
        self.sorted_col = None
        self.sorted_dir = None



    def GetNumberCols(self):
        return len(self.colnames)

    def GetNumberRows(self):
        return len(self.data)

    def GetColLabelValue(self, col):
        return self.colnames[col]

    def GetRowLabelValue(self, row):
        return "%d" % int(self.data[row][0])

    def GetValue(self, row, col):
        return str(self.data[row][1].get(self.GetColLabelValue(col), ""))

    def GetRawValue(self, row, col):
        return self.data[row][1].get(self.GetColLabelValue(col), "")

    def SetValue(self, row, col, value):
        self.data[row][1][self.GetColLabelValue(col)] = value


    def SortColumn(self, col):
        if self.sorted_col == col:
            self.sorted_dir = not self.sorted_dir
        else:
            self.sorted_col = col
            self.sorted_dir = False

        name = self.colnames[col]
        tempdata = []

        for row in self.data:
            rowname, entry = row
            try:
                tempval = float(entry.get(name, None))
            except:
                tempval = entry.get(name, None)
            tempdata.append((tempval, row))

        tempdata.sort(reverse=self.sorted_dir)
        self.data = []

        for sortvalue, row in tempdata:
            self.data.append(row)





class TransitGridFrame(wx.Frame):

    def __init__(self, parent, path, size=(-1,-1)):

        wx.Frame.__init__(self, parent, size=size)


        self.SetTitle(path)

        bSizer1 = wx.BoxSizer( wx.VERTICAL )

        sbSizer1 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Information" ), wx.HORIZONTAL )


        self.parent = parent
        self.path = path
        self.col = 0
        self.row = 0
        self.ctrldata = self.parent.ctrlSelected()
        self.expdata = self.parent.expSelected()
        self.annotation = self.parent.annotation


        line = open(self.path).readline().strip()
        (method, FT) = getInfoFromFileType(line)

          
        self.filetype = FT 
        
        if self.filetype.identifier == "#Unknown":
            self.columnlabels = unknownColNames(self.path)
        else:
            self.columnlabels = self.filetype.colnames

        data = self.filetype.getData(self.path, self.columnlabels)

        wxheader_list = []
        text = self.filetype.getHeader(self.path)
        wxheader_list.append(wx.StaticText( self, wx.ID_ANY, text, wx.DefaultPosition, wx.DefaultSize, 0 ))
        wxheader_list[-1].Wrap( -1 )
        sbSizer1.Add( wxheader_list[-1], 0, wx.ALL, 5 )


        self.grid = wx.grid.Grid(self, -1)

        bSizer1.Add( sbSizer1, 0, wx.EXPAND, 5 )
        bSizer1.Add( self.grid, 1, wx.EXPAND, 5 )
        self.SetSizer( bSizer1 )
        self.Centre( wx.BOTH )

        self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelDoubleClicked)
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClicked)

        mytable = TransitTable(data, self.columnlabels)
        self.grid.SetTable(mytable, True)

        self.grid.EnableEditing(False)
        self.grid.AdjustScrollbars()
        self.grid.SetColLabelSize(wx.grid.GRID_AUTOSIZE)

        self.grid.AutoSizeColumns()
        self.AutoResizeCols()
        self.grid.ForceRefresh()

        (width, height) = bSizer1.GetMinSize()
        max_width = 1500
        max_height = 800
        width = min(width+50, max_width)
        height = min(height, max_height)
        
        self.SetMinSize((width, height))
        

        self.Layout()
        #self.Show()


    def AutoResizeCols(self):
        max_column_size = 200
        min_column_size = 100
        self.grid.AutoSizeColumns(False)
        for i,label in enumerate(self.columnlabels):
            size = self.grid.GetColSize(i)            
            if size > max_column_size:
                self.grid.SetColSize(i, max_column_size)
            elif size < min_column_size:
                self.grid.SetColSize(i, min_column_size)

    def OnLabelDoubleClicked(self, evt):
        col = evt.GetCol()
        if col != -1:
            self.grid.GetTable().SortColumn(col)
            self.grid.ForceRefresh()



    def OnCellRightClicked(self, evt):

        menu = wx.Menu()
        id1 = wx.NewId()
        sortID = wx.NewId()

        xo, yo = evt.GetPosition()
        self.row = self.grid.YToRow(yo) - 1
        self.col = 0
        val = self.grid.GetCellValue(self.row, 0)

        self.Refresh()
        for (menuname, menufunc) in self.filetype.getMenus():
            newid = wx.NewId()
            menu.Append(newid, menuname)
            newmenufunc = partial(menufunc,  self)
            self.Bind(wx.EVT_MENU, newmenufunc, id=newid)

        self.PopupMenu(menu)
        menu.Destroy()



