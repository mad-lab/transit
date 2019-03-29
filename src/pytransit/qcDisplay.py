# -*- coding: utf-8 -*-

###########################################################################
## Python code generated with wxFormBuilder (version Jun  6 2014)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

try:
    import wx
    import wx.xrc
    hasWx = True
except Exception as e:
    hasWx = False


import sys
import os
import io
from PIL import Image
import numpy
import scipy.stats
import matplotlib.pyplot as plt
#matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools




def WxBitmapToPilImage( myBitmap ) :
    return WxImageToPilImage( WxBitmapToWxImage( myBitmap ) )

def WxBitmapToWxImage( myBitmap ) :
    return wx.ImageFromBitmap( myBitmap )

#-----

def PilImageToWxBitmap( myPilImage ) :
    return WxImageToWxBitmap( PilImageToWxImage( myPilImage ) )

def PilImageToWxImage( myPilImage ):
    myWxImage = wx.EmptyImage( myPilImage.size[0], myPilImage.size[1] )
    myWxImage.SetData( myPilImage.convert( 'RGB' ).tobytes() )
    return myWxImage

def WxImageToWxBitmap( myWxImage ) :
    return myWxImage.ConvertToBitmap()







###########################################################################
## Class qcFrame
###########################################################################

class qcFrame ( wx.Frame ):

    def __init__( self, parent, datasets):

        try:
            self.qc_prefix = "[QualityControl]"
            self.index_stats = 0
            self.plots_list = []

            self.wigList = datasets

            wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Quality Control", pos = wx.DefaultPosition, size = wx.Size( 1560, 900 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

            #self.SetSizeHints( wx.DefaultSize, wx.DefaultSize )

            bSizer9 = wx.BoxSizer( wx.VERTICAL )

            self.m_scrolledWindow1 = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.m_scrolledWindow1.SetScrollRate( 5, 5 )
            bSizer10 = wx.BoxSizer( wx.VERTICAL )

            self.plotsScrolledWindow = wx.ScrolledWindow( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.plotsScrolledWindow.SetScrollRate( 5, 5 )
            self.plotsScrolledWindow.SetMinSize( wx.Size( -1, 515 ) )

            #plotsSizer = wx.BoxSizer( wx.VERTICAL )
            plotsSizer = wx.BoxSizer( wx.HORIZONTAL )

            self.plotsBitmap1 = wx.StaticBitmap( self.plotsScrolledWindow, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
            self.plotsBitmap2 = wx.StaticBitmap( self.plotsScrolledWindow, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
            self.plotsBitmap3 = wx.StaticBitmap( self.plotsScrolledWindow, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )

            plotsSizer.Add( self.plotsBitmap1, 0, wx.ALL, 5 )
            plotsSizer.Add( self.plotsBitmap2, 0, wx.ALL, 5 )
            plotsSizer.Add( self.plotsBitmap3, 0, wx.ALL, 5 )

            #self.plotsBitmap.SetMaxSize( wx.Size( 400,400 ) )
            #self.plotsFigure = Figure()
            #self.plotsAxes = self.plotsFigure.add_subplot(111)
            #self.plotsCanvas = FigureCanvas(self, -1, self.plotsFigure)
            #plotsSizer.Add( self.plotsCanvas, 0, wx.ALL, 5 )

            self.plotsScrolledWindow.SetSizer( plotsSizer )
            self.plotsScrolledWindow.Layout()
            plotsSizer.Fit( self.plotsScrolledWindow )
            bSizer10.Add( self.plotsScrolledWindow, 0, wx.EXPAND |wx.ALL, 5 )

            self.statsScrolledWindow = wx.ScrolledWindow( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.statsScrolledWindow.SetScrollRate( 5, 5 )
            self.statsScrolledWindow.SetMaxSize( wx.Size( -1, -1 ) )
            self.statsScrolledWindow.SetMinSize( wx.Size( -1, 300 ) )


            NoteText = """*Note: Plot 1 and 2 truncate the top 1% of reads for readability.
 Selecting a normalization method from the drop down will normalize the data and refresh the figures and table.
 This may take a long time depending on the normalization method chosen."""

            #self.noticeLabel = wx.StaticText( self.statsScrolledWindow, wx.ID_ANY, u"*Note: Plot 1 and 2 truncate the top 1% of reads for readability.", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE)
            self.noticeLabel = wx.StaticText( self.statsScrolledWindow, wx.ID_ANY, NoteText, wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_LEFT)

            normChoices = sorted(norm_tools.methods.keys()) #[ u"nonorm", "TTR", "betageom"]
            self.normChoice = wx.Choice( self.statsScrolledWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, normChoices, 0 )
            #self.normChoice.SetSelection( 0 )
            self.normChoice.SetStringSelection("nonorm")
            #noteSizer = wx.BoxSizer( wx.VERTICAL )
            #noteSizer.Add(self.noticeLabel, wx.ALL|wx.EXPAND, 5 )


            self.normLabel = wx.StaticText( self.statsScrolledWindow, wx.ID_ANY, u"Normalization:", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE)
            statsSizer = wx.BoxSizer( wx.VERTICAL )
            normSizer = wx.BoxSizer( wx.HORIZONTAL )

            self.statsListCtrl = wx.ListCtrl( self.statsScrolledWindow, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1, 140 ), wx.LC_REPORT |wx.LC_SINGLE_SEL )


            normSizer.Add(self.normLabel, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
            normSizer.Add(self.normChoice, 0, wx.ALL, 5)

            statsSizer.Add( self.noticeLabel, 0, wx.EXPAND, 5 )
            #statsSizer.Add( self.normChoice, 0, wx.ALL, 5 )
            statsSizer.Add(normSizer, 0, wx.ALL, 5)
            statsSizer.Add( self.statsListCtrl, 1, wx.ALL|wx.EXPAND, 5 )

            

            self.statsScrolledWindow.SetSizer( statsSizer )
            self.statsScrolledWindow.Layout()
            statsSizer.Fit( self.statsScrolledWindow )
            bSizer10.Add( self.statsScrolledWindow, 0, wx.EXPAND |wx.ALL, 5 )


            self.m_scrolledWindow1.SetSizer( bSizer10 )
            self.m_scrolledWindow1.Layout()
            bSizer10.Fit( self.m_scrolledWindow1 )
            bSizer9.Add( self.m_scrolledWindow1, 1, wx.EXPAND |wx.ALL, 5 )


            self.SetSizer( bSizer9 )
            self.Layout()

            self.Centre( wx.BOTH )

            ########################
            # Connect Events
            self.statsListCtrl.Bind( wx.EVT_LIST_ITEM_SELECTED, self.onStatsItemSelect)
            self.normChoice.Bind( wx.EVT_CHOICE, self.onNormSelect )
            self.Bind(wx.EVT_CLOSE, self.OnClose)


            #######
            self.index_stats = 0
            self.statsListCtrl.InsertColumn(0, 'File', width=250)
            self.statsListCtrl.InsertColumn(1, 'Density', wx.LIST_FORMAT_CENTRE, width=85)
            self.statsListCtrl.InsertColumn(2, 'Mean Read', wx.LIST_FORMAT_CENTRE, width=85)
            self.statsListCtrl.InsertColumn(3, 'NZMean Read', wx.LIST_FORMAT_CENTRE, width=115)
            self.statsListCtrl.InsertColumn(4, 'NZMedian Read', wx.LIST_FORMAT_CENTRE, width=125)
            self.statsListCtrl.InsertColumn(5, 'Max Read', wx.LIST_FORMAT_CENTRE, width=85)
            self.statsListCtrl.InsertColumn(6, 'Total Reads', wx.LIST_FORMAT_CENTRE, width=85)
            self.statsListCtrl.InsertColumn(7, 'Skew', wx.LIST_FORMAT_CENTRE, width=85)
            self.statsListCtrl.InsertColumn(8, 'Kurtosis', wx.LIST_FORMAT_CENTRE, width=85)




            ############################
            self.norm = "nonorm"
            (self.data, self.position) = tnseq_tools.get_data(self.wigList)
        

            self.refresh()
            #self.updateFiles()
            #self.addPlots()
            #self.statsListCtrl.Select(0)
            #self.onStatsItemSelect(None)
            ###########################
            #self.bSizer9.Fit()

        except Exception as e:
            print(self.qc_prefix, "Error:", e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)




    def __del__( self ):
        pass


    def refresh(self):
        try:
            #(self.data, self.position) = tnseq_tools.get_data(self.wigList)
            self.plots_list = []
            self.statsListCtrl.DeleteAllItems()
            (self.normdata, factors) = norm_tools.normalize_data(self.data, self.norm)
            self.updateFiles()
            self.addPlots()
            self.statsListCtrl.Select(0)
            self.refreshPlots()
        except Exception as e:
            print(self.qc_prefix, "Error:", e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def updateFiles(self):
        self.index_stats = 0

        try:
            for j,row in enumerate(self.normdata):
                name = transit_tools.basename(self.wigList[j])
                (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = tnseq_tools.get_data_stats(row)
                self.statsListCtrl.InsertStringItem(self.index_stats, name)
                self.statsListCtrl.SetStringItem(self.index_stats, 1, "%1.1f" % (density*100.0))
                self.statsListCtrl.SetStringItem(self.index_stats, 2, "%1.1f" % (meanrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 3, "%1.1f" % (nzmeanrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 4, "%1.1f" % (nzmedianrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 5, "%d" % (maxrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 6, "%d" % (totalrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 7, "%1.1f" % (skew))
                self.statsListCtrl.SetStringItem(self.index_stats, 8, "%1.1f" % (kurtosis))
                print(self.qc_prefix, "Adding dataset (%d): %s" % (self.index_stats, name))
                self.index_stats+=1
        except Exception as e:
            print(self.qc_prefix, "Error:", e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)


    def addPlots(self):
        try:
            for i,reads in enumerate(self.normdata):
                #Data
                name = transit_tools.basename(self.wigList[i])
                #reads,position = tnseq_tools.get_data([])
                nzreads = reads[reads>0]
                n_nz = len(nzreads)
                truncnzreads = sorted(nzreads)[:int(n_nz*0.99)]

                fig = plt.figure(facecolor='white', figsize=(5, 5), dpi=100)
                ax = fig.add_subplot(111, frame_on=False)

                #Plot 1
                n, bins, patches = ax.hist(truncnzreads, density=1, facecolor='c', alpha=0.75, bins=100)
                plt.xlabel('Reads')
                plt.ylabel('Probability')

                plt.title('Histogram of Non-Zero Reads\nDataset: %s' % name)
                plt.grid(True)
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                hist_pil_img = Image.open(buf)
                hist_wx_img = PilImageToWxImage( hist_pil_img )
                hist_wx_bitmap = wx.BitmapFromImage(hist_wx_img)

                #Plot 2
                plt.clf()
                ax = fig.add_subplot(111, frame_on=False)
                ((qtheoretical, qdata), (slope, intercept, r)) = scipy.stats.probplot(truncnzreads, dist="geom", sparams=(1.0/numpy.mean(truncnzreads),))

                n_qdata = len(qdata)
                #qtheoretical = qtheoretical[:int(n_qdata*0.99)]
                #qdata = qdata[:int(n_qdata*0.99)]
                ax.plot(qdata, qtheoretical, "ob")

                maxval = max(numpy.max(qtheoretical), numpy.max(qdata))
                ax.plot((0, maxval), (0, maxval), ls="-", c="r")

                plt.xlabel("Data Quantiles")
                plt.ylabel("Theoretical Quantiles")
                plt.title('QQ-Plot with Theoretical Geom\nDataset: %s' % name)
                plt.grid(True)
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                qq_pil_img = Image.open(buf)
                qq_wx_img = PilImageToWxImage(qq_pil_img)
                qq_wx_bitmap = wx.BitmapFromImage(qq_wx_img)

                #Plot 3
                plt.clf()
                ax = fig.add_subplot(111, frame_on=False)
                ax.plot(sorted(reads), "ob", linewidth=3)
                plt.title('Rank Reads\nDataset: %s' % name)
                plt.xlabel("Reads")
                plt.ylabel("Rank")
                plt.grid(True)
                plt.xticks(rotation=45)
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                sorted_pil_img = Image.open(buf)
                sorted_wx_img = PilImageToWxImage(sorted_pil_img)
                sorted_wx_bitmap = wx.BitmapFromImage(sorted_wx_img)

                self.plots_list.append([fig, hist_wx_bitmap, qq_wx_bitmap, sorted_wx_bitmap])


        except Exception as e:
            print(self.qc_prefix, "Error:", e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)


    def onNormSelect(self, event):
        self.norm = self.normChoice.GetString(self.normChoice.GetCurrentSelection())
        print(self.qc_prefix, "Normalizing data using '%s' method" % (self.norm))
        self.refresh() 


    def refreshPlots(self):
        ii = self.statsListCtrl.GetFirstSelected()
        fig, hist_wx_bitmap, qq_wx_bitmap, sorted_wx_bitmap = self.plots_list[ii]
        self.plotsBitmap1.SetBitmap(hist_wx_bitmap)
        self.plotsBitmap2.SetBitmap(qq_wx_bitmap)
        self.plotsBitmap3.SetBitmap(sorted_wx_bitmap)


    def onStatsItemSelect(self, event):
        ii = self.statsListCtrl.GetFirstSelected()
        print(self.qc_prefix, "Showing plots for ", self.statsListCtrl.GetItem(ii, 0).GetText())
        self.refreshPlots()


    def statsSelected(self, col=0):
        selected_stats = []
        current = -1
        while True:
            next = self.statsListCtrl.GetNextSelected(current)
            if next == -1:
                break
            path = self.statsListCtrl.GetItem(next, col).GetText()
            selected_stats.append(path)
            current = next
        return selected_stats

    def OnClose(self, event):
        for (fig, a,b,c) in self.plots_list:
            plt.close(fig)
        self.Destroy()


