import wx
import view_trash
import draw_trash
from math import *
import os
import ntpath
import Image
import ImageDraw
import ImageFont




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



def fetch_name(filepath):
    return os.path.splitext(ntpath.basename(filepath))[0]



#inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TrashFrame(view_trash.MainFrame):
    #constructor
    def __init__(self,parent, dataset_list=["H37Rv_Sassetti_glycerol.wig"], annotation="H37Rv.prot_table"):



        self.size = wx.Size( 1500,800 )
        self.start = 1
        self.end = 10000

        self.orf2data = draw_trash.read_prot_table(annotation)
        self.hash = draw_trash.hash_prot_genes(annotation)


        self.labels = [fetch_name(d) for d in dataset_list]


        self.fulldata = []
        for dataset in dataset_list:
            temp = []
            for line in open(dataset):
                if line.startswith("#"): continue
                if line.startswith("variable"): continue
                if line.startswith("location"): continue
                tmp = line.split()
                pos = int(tmp[0])
                read = int(tmp[1])
                temp.append((pos,read))
            self.fulldata.append(temp)
            


        #initialize parent class
        view_trash.MainFrame.__init__(self,parent)

        self.updateFunc(parent)   
        self.Fit()
 
    def updateFunc(self,event):
        try:
            start = int(self.startText.GetValue())
            end = int(self.endText.GetValue())

            min_read = int(self.minText.GetValue())
            max_read = int(self.maxText.GetValue())

    
            image_pil = draw_trash.draw_canvas(self.fulldata, self.hash, self.orf2data, labels=self.labels, min_read=min_read, max_read=max_read, start=start, end=end)

            #image_pil = draw_trash.draw_canvas(start, end, min_read, max_read, self.data, self.hash, self.orf2data)
            #image_wxBit = PilImageToWxBitmap( image_pil )
            image_wxImg = PilImageToWxImage( image_pil )
            #self.m_bitmap1 = wx.StaticBitmap(self, -1, image_wxBit, (0, 0))
            self.m_bitmap1.SetBitmap(wx.BitmapFromImage(image_wxImg))
            #self.m_bitmap1 = image_wxBit
            self.Refresh()
            image_pil = ""

        except Exception, e:
            print 'error', e
    #put a blank string in text when 'Clear' is clicked

    def leftFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end-start
        new_start = int(start - delta*0.10)
        new_end = int(end - delta*0.10)

        #new_start = start - 500
        #new_end = end - 500
        
        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.updateFunc(event) 


    def rightFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end-start
        new_start = int(start + delta*0.10)
        new_end = int(end + delta*0.10)

        #new_start = start + 500
        #new_end = end + 500

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.updateFunc(event)

    
    def zoomInFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end-start
        new_start = int(start + delta*0.10)
        new_end = int(end - delta*0.10)

        #new_start = start + 500
        #new_end = end - 500

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.updateFunc(event)


    def zoomOutFunc(self, event):
        start = int(self.startText.GetValue())
        end = int(self.endText.GetValue())
        delta = end-start
        new_start = int(start - delta*0.10)
        new_end = int(end + delta*0.10)

        #new_start = start - 500
        #new_end = end + 500

        self.startText.SetValue(str(new_start))
        self.endText.SetValue(str(new_end))
        self.updateFunc(event)    


    def resetFunc(self, event):
        self.startText.SetValue("1")
        self.endText.SetValue("10000")
        self.minText.SetValue("0")
        self.maxText.SetValue("150")
        self.updateFunc(event)


    def saveImageFunc(self, event):
        finished_image = self.m_bitmap1.GetBitmap()
        output_path = self.SaveFile(DIR=".", FILE="reads_canvas.png")
        if output_path:
            finished_image.SaveFile(output_path, wx.BITMAP_TYPE_PNG)
        


    def SaveFile(self, DIR=".", FILE="", WC=""):
        """
        Create and show the Save FileDialog
        """
        path = ""
        dlg = wx.FileDialog(
            self, message="Save file as ...",
            defaultDir=DIR,
            defaultFile=FILE, wildcard=WC, style=wx.SAVE|wx.OVERWRITE_PROMPT
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print "You chose the following output filename: %s" % path
        dlg.Destroy()
        return path



#mandatory in wx, create an app, False stands for not deteriction stdin/stdout
#refer manual for details
if __name__ == "__main__":
    app = wx.App(False)
     
    #create an object of CalcFrame
    frame = TrashFrame(None)
    #show the frame
    frame.Show(True)
    #start the applications
    app.MainLoop()


