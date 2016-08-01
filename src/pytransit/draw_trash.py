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


import view_trash
from math import *
import os
import platform

try:
    import Image
    import ImageDraw
    import ImageFont
except ImportError:
    import PIL.Image as Image
    import PIL.ImageDraw as ImageDraw
    import PIL.ImageFont as ImageFont


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


linuxFonts = []
linuxFonts.append("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans-Bold.ttf")
linuxFonts.append("/usr/share/fonts/dejavu-lgc/DejaVuLGCSerifCondensed-Bold.ttf")
linuxFonts.append("/usr/share/fonts/dejavu-lgc/DejaVuLGCSansCondensed-Bold.ttf")

winFonts = []
winFonts.append("consolab.ttf")
winFonts.append("courb.ttf")
winFonts.append("arial.ttf")
fontsize = 16

font = ImageFont.load_default()
if platform.system() == "Linux":
    for fontpath in linuxFonts:
        if os.path.isfile(fontpath):
            font = ImageFont.truetype(fontpath, fontsize)
            break

elif platform.system() == "Windows":
    for fontpath in winFonts:
        try:
            font = ImageFont.truetype(fontpath, fontsize)
            break    
        except:
            pass


def draw_reads(draw, reads, ta_sites, start_x=0, start_y=0, width=400, height=100, start=0, end=500, min_read=0, max_read=500, lwd=2):


    TRUNC_READS = [min(rd, max_read) for rd in reads]
    NORM_READS = [normalize(rd, 0, max_read, 0, max_read) for rd in TRUNC_READS]

    new_min_w = start_x
    new_max_w = start_x + width #- self.padding_r

    new_min_h = start_y
    new_max_h = start_y + height

    for i,TA in enumerate(ta_sites):

        TApos = normalize(TA, start, end, new_min_w, new_max_w)


        if NORM_READS[i] == 0: continue
           
        read_h = normalize(NORM_READS[i], 0, max_read, new_min_h, new_max_h) # height of read line

        Y1 = start_y
        Y2 = start_y - (read_h-start_y)
        draw.line([(TApos, Y1), (TApos, Y2)], width=lwd, fill=(255,0,0))



def draw_ta_sites(draw, ta_sites, start_x=0, start_y=0, width=200, height=0, start=0, end=500, lwd=2):

    new_min_w = start_x
    new_max_w = start_x + width #- self.padding_r
    for i,TA in enumerate(ta_sites):
        TApos = normalize(TA, start, end, new_min_w, new_max_w)
        draw.line([(TApos, start_y+0), (TApos, start_y + height)], width=lwd, fill="black")




def draw_scale(draw, start_x, start_y, height, max_read):

    #print "scale", start_x, start_y, height
    MIDREAD = int(max_read/2.0)
    top_text_w, top_text_h = draw.textsize(str(max_read), font=font)
    draw.text((start_x, start_y), str(max_read), font=font, fill="black")
  

    draw.text((start_x, start_y + height/2.0), str(MIDREAD), font=font, fill="black")

    bottom_text_w, bottom_text_h = draw.textsize(str(MIDREAD), font=font)
    draw.text((start_x+bottom_text_w-(top_text_w/2.0), start_y+height), "0", font=font, fill="black")



def draw_features(draw, features, start, end, start_x, start_y, width, height):
    pass


def draw_genes(draw, GENES, orf2data, start, end, start_x, start_y, width, height):

    padding_h = 3
    text_w, text_h = draw.textsize("RV0001", font=font)        
    gene_h = height - text_h

    #print "GENES height", height
    #print "GENES text_h", text_h
    #print "GENES gene_h", gene_h

    triangle_size = 10
    for gene in GENES:

        if gene not in orf2data: continue
        gene_start = orf2data[gene][0]
        gene_end = orf2data[gene][1]
        strand = orf2data[gene][2]
        name = orf2data[gene][3]

        new_min = start_x
        new_max = start_x + width

        norm_start = normalize(max(orf2data[gene][0], start), start, end, new_min, new_max)
        norm_end = normalize(min(orf2data[gene][1], end), start, end, new_min, new_max)


        #if True:
        #    print gene, name, gene_start, gene_end, strand #, norm_start, norm_end

        if strand == "-":
    
            if gene_start >= start:
                draw.rectangle(((norm_start+triangle_size, start_y+5),(norm_end,start_y+gene_h-5)), fill="blue")
                #draw.polygon([(norm_start,start_y+gene_h/2.0),(norm_start, gene_h+20+5), (norm_start,ta_sites_finish_h +gene_h+10)], fill="blue" )
                draw.polygon([(norm_start+triangle_size, start_y),(norm_start+triangle_size,start_y+gene_h), (norm_start,start_y+gene_h/2.0)], fill="blue" )
    
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end,start_y+gene_h-5)), fill="blue")
                #draw.rectangle(((norm_start, start_y),(norm_end,start_y+gene_h)), fill="blue")

        else:
            if gene_end <= end:
                draw.rectangle(((norm_start, start_y+5),(norm_end-triangle_size, start_y+gene_h-5)), fill="blue")
                draw.polygon([(norm_end-triangle_size, start_y),(norm_end-triangle_size,start_y+gene_h), (norm_end,start_y+gene_h/2.0)], fill="blue" )
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end, start_y+gene_h-5)), fill="blue")
                #draw.rectangle(((norm_start,start_y ),(norm_end, start_y+gene_h)), fill="blue")


        if name == "-": name = gene
        if not name.startswith("non-coding"):
            name_text_w, name_text_h = draw.textsize(name, font=font)
            if abs(norm_start-norm_end) >= name_text_w:
                draw.text(( norm_start + (abs(norm_start-norm_end) - name_text_w)/2.0 , start_y+gene_h+text_h), name, font=font, fill="black")





def get_dynamic_height(N):

    #Set rest of heights and widths
    read_h = 100
    gene_h = 50
    ta_h = 20
    padding_h = 3
    canvas_h = read_h*N + ta_h + gene_h + padding_h + padding_h + 80
    return (canvas_h)


def draw_canvas(fulldata, hash, orf2data, labels=[], min_read=0, max_read=2000, start=1, end=500, canvas_h=-1, canvas_w=1000):
    

    temp_image = Image.new("RGB",(200, 200),"white")
    temp_draw = ImageDraw.Draw(temp_image)
    #Set main draw object

    N = len(fulldata)
    #Set Labels
    if not labels:
        labels= ["Read Counts"]*N
    


    #Get dynamic text widths
    #print "Labels:"
    max_label_w = 0
    for L in labels:
        label_text_w, label_text_h = temp_draw.textsize(L, font=font)
        max_label_w = max(label_text_w, max_label_w)
        #print L

    scale_text_w, scale_text_h = temp_draw.textsize(str(max_read), font=font)
    


    #Set rest of heights and widths
    read_h = 100
    gene_h = 50
    ta_h = 20
    padding_w = 3
    padding_h = 3
    read_w = canvas_w - (max_label_w + scale_text_w + padding_w + padding_w + 30)

    if canvas_h == -1:
        canvas_h = read_h*N + ta_h + gene_h + padding_h + padding_h + 80
    


    image = Image.new("RGB",(canvas_w, canvas_h),"white")
    draw = ImageDraw.Draw(image)

    lwd = 2

    GENES = []
    TA_SITES = []
    READS = []
    nc_count = 1
    for j,data in enumerate(fulldata):
        #print j
        temp = []
        for pos,read in data:
            if start <= pos <= end:
                gene = hash.get(pos,"non-coding")
                if gene == "non-coding" and len(GENES) > 0 and not GENES[-1].startswith("non-coding"):
                    gene+="_%d" % nc_count
                    nc_count +=1
                if j ==0:
                    if gene not in GENES: GENES.append(gene)
                    TA_SITES.append(pos)
                temp.append(read)
        READS.append(temp)



    #print READS
    #print "start", start
    #print "end", end
    #print len(READS), len(TA_SITES)
    #print ""
    #for rd in READS:
    #    print rd
    #print ""

    start_x = max_label_w + padding_w + 21
    draw.line([(start_x, 0), (start_x, canvas_h)], width=lwd, fill="black")
    start_y = 0
    half = 100*0.5
    start_x += 5
    for j in range(len(fulldata)):
        start_y+=read_h+padding_h
        draw.text((10, start_y - half), labels[j], font=font, fill="black")
        draw_reads(draw, READS[j], TA_SITES, start_x, start_y, read_w, read_h, start, end, min_read, max_read)
        draw_scale(draw, start_x+read_w+padding_w+2, start_y-100+10, 70, max_read)
            



    start_y+=10
    #start_x+=5

    #TA sites
    draw.text((30, start_y),'TA Sites', font=font, fill="black")
    draw_ta_sites(draw, TA_SITES, start_x, start_y, read_w, ta_h, start, end)

    #Genes
    start_y += 50
    draw.text((30, start_y+10),'Genes', font=font, fill="black")
    width = read_w
    draw_genes(draw, GENES, orf2data, start, end, start_x, start_y, width, gene_h)

    return(image)







