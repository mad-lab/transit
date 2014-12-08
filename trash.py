import wx
import view_trash
from math import *
import os
import Image
import ImageDraw
import ImageFont



def normalize(X, old_min, old_max, new_min, new_max):
    old_range = (old_max - old_min)
    new_range = (new_max - new_min)
    if old_range == 0:
        return new_min
    else:
        return (((X - old_min) * new_range) / old_range) + new_min
 

def read_prot_table(path):
    orf2data = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf2data[tmp[8]] = [int(tmp[1]), int(tmp[2]), tmp[3], tmp[7]]
    return(orf2data)


def hash_prot_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        start, end = int(tmp[1]), int(tmp[2])
        for i in range(start,end+1):
            #if i not in hash:
            #    hash[i] = tmp[8]
            hash[i] = tmp[8]
    return hash




def draw_canvas(start, end, min_read, max_read, data, hash, orf2data):

    GENES = []
    TA_SITES = []
    READS = []
    nc_count = 1
    for pos,read in data:
        if start <= pos <= end:
            gene = hash.get(pos,"non-coding")
            if gene == "non-coding" and len(GENES) > 0 and not GENES[-1].startswith("non-coding"):
                gene+="_%d" % nc_count
                nc_count +=1
            if gene not in GENES: GENES.append(gene)
            TA_SITES.append(pos)
            READS.append(read)



    #MAXREAD = 500
    MAXREAD = max_read
    canvas_h = 300
    canvas_w = 700
    drawBorder = False
    VERBOSE = False
    
    ######### Variables ##########
    fontpath="/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans-Bold.ttf"

    #Genes
    genebar_h = 30
    genebar_w = 200


    # Headings menu to left (i.e. read counts, genes, etc/)
    heading_w = 100
    heading_h = canvas_h
    read_counts_w = 5
    read_counts_h = 60
    ta_sites_w = 5
    ta_sites_h = 10
    genes_w = 5
    genes_h = 50
    lwd =2


    #Display
    padding = 5
    padding_r = 40
    fontsize = 16



    ######### Drawing #########
    image = Image.new("RGB",(canvas_w, canvas_h),"white")
    draw = ImageDraw.Draw(image)
    font1 = ImageFont.truetype(fontpath,fontsize)

    if drawBorder:
        draw.line(((0,0),(0, canvas_h)), fill="black", width=4)
        draw.line(((0,0),(canvas_w, 0)), fill="black", width=4)
        draw.line(((canvas_w,canvas_h),(0, canvas_h)), fill="black", width=4)
        draw.line(((canvas_w,canvas_h),(canvas_w, 0)), fill="black", width=4)

    # READ COUNTS - Header
    draw.text((read_counts_w, read_counts_h),'Read Counts', font=font1, fill="black")
    read_text_w, read_text_h = draw.textsize('Read Counts', font=font1)
    draw.line([(2*read_counts_w+read_text_w, 2*read_counts_h+read_text_h), (canvas_w, 2*read_counts_h+read_text_h)], width=lwd, fill="black")
    draw.line([(2*read_counts_w+read_text_w, 0), (2*read_counts_w+read_text_w, heading_h)], width=lwd, fill="black")


    #TA SITES - Header
    read_counts_finish_h = 2*read_counts_h+read_text_h+lwd
    ta_text_w, ta_text_h = draw.textsize('TA Sites', font=font1)
    draw.text(( (2*read_counts_w+read_text_w - ta_text_w)/2.0, read_counts_finish_h +ta_sites_h),'TA Sites', font=font1, fill="black")
    draw.line([(2*read_counts_w+read_text_w, read_counts_finish_h + 2*ta_sites_h + ta_text_h), (canvas_w, read_counts_finish_h + 2*ta_sites_h + ta_text_h)], width=lwd, fill="black")



    #GENES - Header
    ta_sites_finish_h = read_counts_finish_h + 2*ta_sites_h + ta_text_h + 2
    gene_text_w, gene_text_h = draw.textsize('Genes', font=font1)
    draw.text(((2*read_counts_w+read_text_w - gene_text_w)/2.0, ta_sites_finish_h +genes_h),'Genes', font=font1, fill="black")


    #Reads
    TRUNC_READS = [min(rd, MAXREAD) for rd in READS]
    NORM_READS = [normalize(rd, 0, MAXREAD, 0, MAXREAD) for rd in TRUNC_READS]


    new_min = 2*read_counts_w + read_text_w + lwd + padding
    new_max = canvas_w - padding_r

    new_min2 =  padding
    new_max2 = 2*read_counts_h+read_text_h - padding

    for i,TA in enumerate(TA_SITES):

        TApos = normalize(TA, start, end, new_min, new_max)

        #Draw TA line
        draw.line([(TApos, 2*read_counts_h+read_text_h+5), (TApos, read_counts_finish_h + 2*ta_sites_h + ta_text_h-5)], width=lwd, fill="black")


        if NORM_READS[i] == 0: continue
        #Draw Read Line
        read_h = normalize(NORM_READS[i], 0, MAXREAD, new_min2, new_max2)
        draw.line([(TApos, 2*read_counts_h+read_text_h - padding), (TApos, read_counts_finish_h-padding-read_h)], width=lwd, fill=(255,0,0))



    #SCALE:

    MIDREAD = int(MAXREAD/2.0)
    max_text_w, max_text_h = draw.textsize(str(MAXREAD), font=font1)
    draw.text((canvas_w-max_text_w-5, 0+ padding), str(MAXREAD), font=font1, fill="black")

    mid_text_w, mid_text_h = draw.textsize(str(MIDREAD), font=font1)
    draw.text((canvas_w-mid_text_w-5, (0+ padding + 2*read_counts_h+read_text_h - padding-15)/2.0   ), str(MIDREAD), font=font1, fill="black")

    draw.text((canvas_w-20, 2*read_counts_h+read_text_h - padding-15), "0", font=font1, fill="black")



    #GENES
    triangle_size = 10
    for gene in GENES:


        if gene not in orf2data: continue
        gene_start = orf2data[gene][0]
        gene_end = orf2data[gene][1]
        strand = orf2data[gene][2]
        name = orf2data[gene][3]

        norm_start = normalize(max(orf2data[gene][0], start), start, end, new_min, new_max)
        norm_end = normalize(min(orf2data[gene][1], end), start, end, new_min, new_max)


        if VERBOSE:
            print gene, name, start, end, strand

        if strand == "-":

            if gene_start >= start:
                draw.rectangle(((norm_start+triangle_size,ta_sites_finish_h +genes_h),(norm_end,ta_sites_finish_h +genes_h+20)), fill="blue")
                draw.polygon([(norm_start+triangle_size,ta_sites_finish_h +genes_h-5),(norm_start+triangle_size,ta_sites_finish_h +genes_h+20+5), (norm_start,ta_sites_finish_h +genes_h+10)], fill="blue" )

            else:
                draw.rectangle(((norm_start,ta_sites_finish_h +genes_h),(norm_end,ta_sites_finish_h +genes_h+20)), fill="blue")

        else:
            if gene_end <= end:
                draw.rectangle(((norm_start,ta_sites_finish_h +genes_h),(norm_end-triangle_size,ta_sites_finish_h +genes_h+20)), fill="blue")
                draw.polygon([(norm_end-triangle_size,ta_sites_finish_h +genes_h-5),(norm_end-triangle_size,ta_sites_finish_h +genes_h+20+5), (norm_end,ta_sites_finish_h +genes_h+10)], fill="blue" )
            else:
                draw.rectangle(((norm_start,ta_sites_finish_h +genes_h),(norm_end,ta_sites_finish_h +genes_h+20)), fill="blue")


        if name == "-": name = gene
        if not name.startswith("non-coding"):
            name_text_w, name_text_h = draw.textsize(name, font=font1)
            if abs(norm_start-norm_end) >= name_text_w:
                draw.text(( norm_start + (abs(norm_start-norm_end) - name_text_w)/2.0 , ta_sites_finish_h +genes_h+20+10), name, font=font1, fill="black")



    return(image)



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


#start = 1
#end = 2000
#wig_path = "H37Rv_Sassetti_glycerol.wig"
#pt_path = "H37Rv.prot_table"

#orf2data = read_prot_table(pt_path)
#hash = hash_prot_genes(pt_path)
#data = []
#for line in open(wig_path):
#    if line.startswith("#"): continue
#    if line.startswith("variable"): continue
#    tmp = line.split()
#    pos = int(tmp[0])
#    read = int(tmp[1])
#    data.append((pos,read))




#inherit from the MainFrame created in wxFowmBuilder and create CalcFrame
class TrashFrame(view_trash.MainFrame):
    #constructor
    def __init__(self,parent, dataset="H37Rv_Sassetti_glycerol.wig", annotation="H37Rv.prot_table"):



        self.start = 1
        self.end = 2000
        #self.wig_path = "H37Rv_Sassetti_glycerol.wig"
        #pt_path = "H37Rv.prot_table"

        self.orf2data = read_prot_table(annotation)
        self.hash = hash_prot_genes(annotation)
        self.data = []
        for line in open(dataset):
            if line.startswith("#"): continue
            if line.startswith("variable"): continue
            if line.startswith("location"): continue
            tmp = line.split()
            pos = int(tmp[0])
            read = int(tmp[1])
            self.data.append((pos,read))


        #initialize parent class
        view_trash.MainFrame.__init__(self,parent)
        #self.startText.SetValue("1")
        #self.endText.SetValue("2000")

        self.updateFunc(parent)    
        # pick an image file you have in the working folder
        # you can load .jpg  .png  .bmp  or .gif files
        #image_file = 'temp.png'
        #bmp1 = wx.Image(image_file, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        # image's upper left corner anchors at panel coordinates (0, 0)
        #self.m_bitmap1 = wx.StaticBitmap(self, -1, bmp1, (0, 0))

 
    #what to when 'Solve' is clicked
    #wx calls this function with and 'event' object
    def updateFunc(self,event):
        try:
            #command = "python /pacific/HomeFrozen/michael.dejesus/MISC/draw_trash/draw_trash.py -f %s -pt %s -s %s -e %s -o temp.png" % (wig_path, pt_path, self.startText.GetValue(), self.endText.GetValue())

            #os.system(command)

            
            #image_file = 'temp.png'
            #bmp1 = wx.Image(image_file, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
            ## image's upper left corner anchors at panel coordinates (0, 0)
            #self.m_bitmap1 = wx.StaticBitmap(self, -1, bmp1, (0, 0))

            start = int(self.startText.GetValue())
            end = int(self.endText.GetValue())

            min_read = int(self.minText.GetValue())
            max_read = int(self.maxText.GetValue())

            image_pil = draw_canvas(start, end, min_read, max_read, self.data, self.hash, self.orf2data)
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
        self.startText.SetValue(str(int(self.startText.GetValue())-500))
        self.endText.SetValue(str(int(self.endText.GetValue())-500))
        self.updateFunc(event) 

    def rightFunc(self, event):
        self.startText.SetValue(str(int(self.startText.GetValue())+500))
        self.endText.SetValue(str(int(self.endText.GetValue())+500))
        self.updateFunc(event)
    
    def zoomInFunc(self, event):
        self.startText.SetValue(str(int(self.startText.GetValue())+500))
        self.endText.SetValue(str(int(self.endText.GetValue())-500))
        self.updateFunc(event)

    def zoomOutFunc(self, event):
        self.startText.SetValue(str(int(self.startText.GetValue())-500))
        self.endText.SetValue(str(int(self.endText.GetValue())+500))
        self.updateFunc(event)    

    def resetFunc(self, event):
        self.startText.SetValue("1")
        self.endText.SetValue("2000")
        self.minText.SetValue("0")
        self.maxText.SetValue("150")
        self.updateFunc(event)

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


