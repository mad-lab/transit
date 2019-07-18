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


import pytransit.view_trash as view_trash
from math import *
import os
import platform
import numpy
from PIL import Image, ImageDraw, ImageFont

def normalize(X, old_min, old_max, new_min, new_max):
    old_range = (old_max - old_min)
    new_range = (new_max - new_min)
    if old_range == 0:
        return new_min
    else:
        return (((X - old_min) * new_range) / old_range) + new_min


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

    #print("scale", start_x, start_y, height)
    MIDREAD = int(max_read/2.0)
    top_text_w, top_text_h = draw.textsize(str(max_read), font=font)
    draw.text((start_x, start_y), str(max_read), font=font, fill="black")
  

    draw.text((start_x, start_y + height/2.0), str(MIDREAD), font=font, fill="black")

    bottom_text_w, bottom_text_h = draw.textsize(str(MIDREAD), font=font)
    draw.text((start_x+bottom_text_w-(top_text_w/2.0), start_y+height), "0", font=font, fill="black")



 
def draw_features(draw, GENES, orf2data, start, end, start_x, start_y, width, height): 

    padding_h = 3
    text_w, text_h = draw.textsize("RV0001", font=font)
    gene_h = height - text_h

    triangle_size = 10
    for gene in GENES:

        if gene not in orf2data: continue
        gene_start = orf2data[gene][2]
        gene_end = orf2data[gene][3]
        strand = orf2data[gene][4]
        name = orf2data[gene][0]

        new_min = start_x
        new_max = start_x + width

        norm_start = normalize(max(gene_start, start), start, end, new_min, new_max)
        norm_end = normalize(min(gene_end, end), start, end, new_min, new_max)

        color = "gray"
        if gene.startswith("ES-"):
            color = "red"
        elif gene.startswith("GD-"):
            color = "yellow"
        elif gene.startswith("NE-"):
            color = "blue"
        elif gene.startswith("GA-"):
            color = "green"

        if strand == "-":    
            if gene_start >= start:
                draw.rectangle(((norm_start, start_y+5),(norm_end,start_y+gene_h-5)), fill=color)
        
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end,start_y+gene_h-5)), fill=color)
                
        else:
            if gene_end <= end:
                draw.rectangle(((norm_start, start_y+5),(norm_end, start_y+gene_h-5)), fill=color)
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end, start_y+gene_h-5)), fill=color)


        if name == "-": name = gene
        if not name.startswith("non-coding"):
            name_text_w, name_text_h = draw.textsize(name, font=font)
            if abs(norm_start-norm_end) >= name_text_w:
                draw.text(( norm_start + (abs(norm_start-norm_end) - name_text_w)/2.0 , start_y+gene_h+text_h), name, font=font, fill="black")












def draw_genes(draw, GENES, orf2data, start, end, start_x, start_y, width, height, doTriangle=True):

    padding_h = 3
    text_w, text_h = draw.textsize("RV0001", font=font)        
    gene_h = height - text_h


    triangle_size = 10
    if not doTriangle:
        triangle_size = 0
    for gene in GENES:

        if gene not in orf2data: continue
        gene_start = orf2data[gene][2]
        gene_end = orf2data[gene][3]
        strand = orf2data[gene][4]
        name = orf2data[gene][0]

        new_min = start_x
        new_max = start_x + width

        norm_start = normalize(max(gene_start, start), start, end, new_min, new_max)
        norm_end = normalize(min(gene_end, end), start, end, new_min, new_max)


        if strand == "-":
    
            if gene_start >= start:
                draw.rectangle(((norm_start+triangle_size, start_y+5),(norm_end,start_y+gene_h-5)), fill="blue")
                if doTriangle:
                    draw.polygon([(norm_start+triangle_size, start_y),(norm_start+triangle_size,start_y+gene_h), (norm_start,start_y+gene_h/2.0)], fill="blue" )
    
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end,start_y+gene_h-5)), fill="blue")

        else:
            if gene_end <= end:
                draw.rectangle(((norm_start, start_y+5),(norm_end-triangle_size, start_y+gene_h-5)), fill="blue")
                if doTriangle:
                    draw.polygon([(norm_end-triangle_size, start_y),(norm_end-triangle_size,start_y+gene_h), (norm_end,start_y+gene_h/2.0)], fill="blue" )
            else:
                draw.rectangle(((norm_start, start_y+5),(norm_end, start_y+gene_h-5)), fill="blue")


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


def draw_canvas(fulldata, position, hash, orf2data, feature_hashes, feature_data, labels=[], min_read=0, scale=[500], globalScale = False, start=1, end=500, canvas_h=-1, canvas_w=1000):
    

    temp_image = Image.new("RGB",(200, 200),"white")
    temp_draw = ImageDraw.Draw(temp_image)
    #Set main draw object

    N = len(fulldata)
    Nfeat = len(feature_hashes)
    #Set Labels
    if not labels:
        labels= ["Read Counts"]*N
   

    GENES = []
    FEATURES = [[] for j in range(len(feature_hashes))]
    TA_SITES = []
    READS = []
    nc_count = 1
    for j,data in enumerate(fulldata):
        #print(j)
        temp = []
        for i,read in enumerate(data):
            pos = position[i]
            if start <= pos <= end:
                gene = hash.get(pos,["non-coding"])[0]
                if gene == "non-coding" and len(GENES) > 0 and not GENES[-1].startswith("non-coding"):
                    gene+="_%d" % nc_count
                    nc_count +=1
                if j ==0:
                    if gene not in GENES: GENES.append(gene)
                    TA_SITES.append(pos)
                    for f,f_hash in enumerate(feature_hashes):
                        feat = f_hash.get(pos,["non-coding"])[0]
                        if feat not in FEATURES[f]: FEATURES[f].append(feat)
                temp.append(read)
        READS.append(temp)

    max_reads = []
    if globalScale:
        max_reads = [int(numpy.max(READS))] * len(READS)

    else:
        for j,s in enumerate(scale):
            #print(j,s)
            if s < 0:
                max_reads.append(int(numpy.max(READS[j])))
            else:
                max_reads.append(s)

    #Get dynamic text widths
    #print("Labels:")
    max_label_w = 0
    for L in labels:
        label_text_w, label_text_h = temp_draw.textsize(L, font=font)
        max_label_w = max(label_text_w, max_label_w)
        #print(L)

    scale_text_w, scale_text_h = temp_draw.textsize(str(max(max_reads)), font=font)
    


    #Set rest of heights and widths
    read_h = 100
    gene_h = 50
    ta_h = 20
    padding_w = 3
    padding_h = 3
    read_w = canvas_w - (max_label_w + scale_text_w + padding_w + padding_w + 30)

    if canvas_h == -1:
        canvas_h = read_h*N + ta_h + gene_h + padding_h + padding_h + 80 + (gene_h+padding_h+50)*(Nfeat)
    


    image = Image.new("RGB",(canvas_w, canvas_h),"white")
    draw = ImageDraw.Draw(image)

    lwd = 2


    #print(READS)
    #print("start", start)
    #print("end", end)
    #print(len(READS), len(TA_SITES))
    #print("")
    #for rd in READS:
    #    print(rd)
    #print("")

    start_x = max_label_w + padding_w + 21
    draw.line([(start_x, 0), (start_x, canvas_h)], width=lwd, fill="black")
    start_y = 0
    half = 100*0.5
    start_x += 5
    for j in range(len(fulldata)):
        temp_label_text_w, temp_label_text_h = temp_draw.textsize(labels[j], font=font)
        label_text_x = (start_x/2.0) - (temp_label_text_w/2.0)
        start_y+=read_h+padding_h
        #draw.text((10, start_y - half), labels[j], font=font, fill="black")
        draw.text((label_text_x, start_y - half), labels[j], font=font, fill="black")
        draw_reads(draw, READS[j], TA_SITES, start_x, start_y, read_w, read_h, start, end, min_read, max_reads[j])
        draw_scale(draw, start_x+read_w+padding_w+2, start_y-100+10, 70, max_reads[j])
            



    start_y+=10
    #start_x+=5

    #TA sites
    temp_label_text_w, temp_label_text_h = temp_draw.textsize('TA Sites', font=font)
    label_text_x = (start_x/2.0) - (temp_label_text_w/2.0)
    #draw.text((30, start_y),'TA Sites', font=font, fill="black")
    draw.text((label_text_x, start_y),'TA Sites', font=font, fill="black")
    draw_ta_sites(draw, TA_SITES, start_x, start_y, read_w, ta_h, start, end)

    #Genes
    temp_label_text_w, temp_label_text_h = temp_draw.textsize('Genes', font=font)
    label_text_x = (start_x/2.0) - (temp_label_text_w/2.0)
    start_y += 50
    #draw.text((30, start_y+10),'Genes', font=font, fill="black")
    draw.text((label_text_x, start_y+10),'Genes', font=font, fill="black")
    width = read_w
    draw_genes(draw, GENES, orf2data, start, end, start_x, start_y, width, gene_h)

    start_y += gene_h -20#+ padding_h 
    #Features:
    for f in range(len(FEATURES)):
        start_y += gene_h + padding_h + 25
        temp_label_text_w, temp_label_text_h = temp_draw.textsize('Feature-%d' % (f+1), font=font)
        label_text_x = (start_x/2.0) - (temp_label_text_w/2.0)
        draw.text((label_text_x, start_y+10),'Feature-%d' % (f+1), font=font, fill="black")
        width = read_w
        #print(FEATURES[f])
        #draw_genes(draw, FEATURES[f], feature_data[f], start, end, start_x, start_y, width, gene_h))
        draw_features(draw, FEATURES[f], feature_data[f], start, end, start_x, start_y, width, gene_h)
        start_y +=10

    return(image)







