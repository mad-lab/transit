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




def WxBitmapToPilImage( myBitmap ) :
    return WxImageToPilImage( WxBitmapToWxImage( myBitmap ) )

def WxBitmapToWxImage( myBitmap ) :
    return wx.ImageFromBitmap( myBitmap )

#-----

def PilImageToWxBitmap( myPilImage ) :
    return WxImageToWxBitmap( PilImageToWxImage( myPilImage ) )

def PilImageToWxImage( myPilImage ):
    myWxImage = wx.EmptyImage( myPilImage.size[0], myPilImage.size[1] )
    myWxImage.SetData( myPilImage.convert( 'RGB' ).tostring() )
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

            wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Quality Control", pos = wx.DefaultPosition, size = wx.Size( 1560,720 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
        
            self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        
            bSizer9 = wx.BoxSizer( wx.VERTICAL )
        
            self.m_scrolledWindow1 = wx.ScrolledWindow( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.m_scrolledWindow1.SetScrollRate( 5, 5 )
            bSizer10 = wx.BoxSizer( wx.VERTICAL )
        
            self.plotsScrolledWindow = wx.ScrolledWindow( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.plotsScrolledWindow.SetScrollRate( 5, 5 )
            self.plotsScrolledWindow.SetMinSize( wx.Size( -1, 150 ) )
            
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
            bSizer10.Add( self.plotsScrolledWindow, 1, wx.EXPAND |wx.ALL, 5 )
        
            self.statsScrolledWindow = wx.ScrolledWindow( self.m_scrolledWindow1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
            self.statsScrolledWindow.SetScrollRate( 5, 5 )
            self.statsScrolledWindow.SetMaxSize( wx.Size( -1, 400 ) )
            self.statsScrolledWindow.SetMinSize( wx.Size( -1, 150 ) )
            statsSizer = wx.BoxSizer( wx.VERTICAL )
        
            #self.statsListCtrl = wx.ListCtrl( self.statsScrolledWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT )
            self.statsListCtrl = wx.ListCtrl( self.statsScrolledWindow, wx.ID_ANY, wx.DefaultPosition, wx.Size( -1, 140 ), wx.LC_REPORT |wx.LC_SINGLE_SEL )
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



            ############################
            self.addFiles(datasets)
            self.addPlots(datasets)
            self.statsListCtrl.Select(0)
            #self.onStatsItemSelect(None)
            ###########################
            #self.bSizer9.Fit() 

        except Exception as e:
            print self.qc_prefix, "Error:", e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



 
    def __del__( self ):
        pass
    

    def addFiles(self, datasets):
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


        try:
            for path in datasets:
                name = transit_tools.basename(path)
                (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis) = tnseq_tools.get_wig_stats(path)
                self.statsListCtrl.InsertStringItem(self.index_stats, name)
                self.statsListCtrl.SetStringItem(self.index_stats, 1, "%1.1f" % (density*100.0))
                self.statsListCtrl.SetStringItem(self.index_stats, 2, "%1.1f" % (meanrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 3, "%1.1f" % (nzmeanrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 4, "%1.1f" % (nzmedianrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 5, "%d" % (maxrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 6, "%d" % (totalrd))
                self.statsListCtrl.SetStringItem(self.index_stats, 7, "%d" % (skew))
                self.statsListCtrl.SetStringItem(self.index_stats, 8, "%d" % (kurtosis))
                print self.qc_prefix, "Adding dataset (%d): %s" % (self.index_stats, name)
                self.index_stats+=1
        except Exception as e:
            print self.qc_prefix, "Error:", e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)


    def addPlots(self, datasets):
        try:
            for i,path in enumerate(datasets):
                #Data
                name = transit_tools.basename(path)
                reads,position = tnseq_tools.get_data([path])
                reads = reads[0]
                reads = numpy.array(reads)
                nzreads = reads[reads>0]

                fig = plt.figure(facecolor='white', figsize=(5, 5), dpi=100)
                ax = fig.add_subplot(111, frame_on=False)
                #Plot 1
                n, bins, patches = ax.hist(reads, normed=1, facecolor='c', alpha=0.75, bins=100)
                plt.xlabel('Reads')
                plt.ylabel('Probability')
                plt.title('Histogram of Reads\nDataset: %s' % name)
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
                ((qtheoretical, qdata), (slope, intercept, r)) = scipy.stats.probplot(nzreads, dist="geom", sparams=(1.0/numpy.mean(nzreads),))
                ax.plot(qdata, qtheoretical, "ob")
                ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="r")

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
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                sorted_pil_img = Image.open(buf)
                sorted_wx_img = PilImageToWxImage(sorted_pil_img)
                sorted_wx_bitmap = wx.BitmapFromImage(sorted_wx_img)                

                self.plots_list.append([hist_wx_bitmap, qq_wx_bitmap, sorted_wx_bitmap])


        except Exception as e:
            print self.qc_prefix, "Error:", e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)


    def onStatsItemSelect(self, event):
        ii = self.statsListCtrl.GetFirstSelected()
        print self.qc_prefix, "Showing plots for ", self.statsListCtrl.GetItem(ii, 0).GetText()
        hist_wx_bitmap, qq_wx_bitmap, sorted_wx_bitmap = self.plots_list[ii]
        self.plotsBitmap1.SetBitmap(hist_wx_bitmap)
        self.plotsBitmap2.SetBitmap(qq_wx_bitmap)
        self.plotsBitmap3.SetBitmap(sorted_wx_bitmap)


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

