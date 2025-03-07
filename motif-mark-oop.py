#!/usr/bin/env python
import re
import math
import cairo
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="test")
    parser.add_argument("-m", help="motif list file") #type: str
    parser.add_argument("-f", help="fasta file") #type: str
    return parser.parse_args()

args = get_args()
motif = args.m 
fasta = args.f 

####################
# GLOBAL VARIABLES #
####################

#define line metrics
line_width = 1
#define feature variables 
feature_height = 25
#color list for drawing motifs: 
# list of 5 tupols containg RGB values divided 255
#pink, purple, light blue, green, orange
COLOR_LIST = [(255/255, 96/255, 208/255),
              (160/255, 32/255, 255/255),
              (80/255, 208/255, 255/255),
              (96/255, 255/255, 128/255),
              (255/255, 106/255, 16/255)]
#motif variations dict key: given motif, value: regex motif
MOTIF_DICT = {}
#key: given letter value: regex possible letters for motifs
IUPAC_DICT ={
    "U":"[TU]",
    "W":"[ATU]",
    "S":"[CG]",
    "M":"[AC]",
    "K":"[GTU]",
    "R":"[AG]",
    "Y":"[CTU]",
    "B":"[CGTU]",
    "D":"[AGTU]",
    "H":"[ACTU]",
    "V":"[ACG]",
    "N":"[ACGTU]",
}

##################
# DEFINE CLASSES #
##################
 
class Motif: 
    def __init__(self, start, length, gene_number, color):
        self.start = start
        self.length = length 
        self.gene_number = gene_number
        self.color = color
    #draw motif
    def draw(self, context):
        context.set_source_rgb(self.color[0], self.color[1], self.color[2])       #set color to color list
        y = (gene_num*100) + 137.5
        x = self.start + 50
        length = self.length
        height = feature_height
        context.rectangle(x,y,length,height)
        context.fill()
    
    # def legend(self, context):
    #     context.set_source_rgb = (Motif.draw)



class Gene: 
    def __init__(self, length, gene_name, gene_number):
        self.length = length 
        self.gene_name = gene_name
        self.gene_number = gene_number


    def draw(self, context):
        #draw the line
        y = (self.gene_number *100) + 150
        context.set_source_rgb(0,0,0)        #set color to black 
        end_point = self.length + 50
        context.set_line_width(line_width)      #Line density
        context.move_to(50,y)        #(x,y) START POINT
        context.line_to(end_point,y)       #END POINT
        context.stroke()                #DRAW

        #draw header 
        title_height = (self.gene_number * 100) + 100
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(10)
        context.move_to(50, title_height)
        context.show_text(self.gene_name)


class Exon: 
    def __init__(self, start, length, gene_number):
        self.start = start
        self.length = length
        self.gene_number = gene_number

    def draw(self, context):
        context.set_source_rgb(0,0,0)        #set color to black
        y = (gene_num*100) + 137.5
        x = self.start + 50
        length = self.length
        height = feature_height
        context.rectangle(x,y,length,height)
        context.fill()

#############
# FUNCTIONS #
#############

#def find and replace the regex expressions of the
def reg_ex_replace(motif) -> str: 
    '''Takes in motif and outputs seq containing possible variations in regex'''
    motif = motif.upper()
    for char in motif: 
        if char in IUPAC_DICT:
            motif = motif.replace(char, IUPAC_DICT[char])
    return(motif)

#function to find the exon start and length
def find_exon(sequence: str) -> tuple:
    '''Takes a DNA seq and returns the start position and total length'''
    for i in range(len(sequence)):
        if sequence[i].isupper():
            start = int(i)
            break
    for i in range(start, len(sequence)):
        if sequence[i].islower():
            end = int(i)
            break
    length = int(end) - int(start)
    return (start, length)

# open the fasta file and read in the fasta records: 
#fasta dict: header = key seq = value 
fasta_dict = {}
with (open(fasta, "r") as fasta_file,
      open (motif, "r") as motif_file):
    #OPEN MOTIF FILE AND POP MOTIF DICT
    for motif in motif_file: 
        motif = motif.strip()
        new_motif = reg_ex_replace(motif)
        MOTIF_DICT[motif] = new_motif

    for line in fasta_file: 
        #set the header and seq variables equal to the header and seq of that record 
        if line.startswith(">"):
            header = line.strip()
            fasta_dict[header] = ""
        else: 
            fasta_dict[header] += line.strip()
    #iterate through fasta dict to find and draw genes

##########################
# CREATE CANVAS + LEGEND #
##########################

height = (len(fasta_dict) * 100) + 100
length = 0 
for header in fasta_dict:
    if len(fasta_dict[header]) > length:
        length = len(fasta_dict[header]) 
length += 250

#set the coordinates of the canvas
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, length, height)
#make the backgroud
context = cairo.Context(surface)
context.set_source_rgb(1,1,1)
context.paint()


################################
# CREATE GENE + EXONS + MOTIFS #
################################

for i, header in enumerate(fasta_dict):
    #FIND OVERALL GENE INFO: 
    gene_name = header
    gene_num = i 
    gene_length =  len(fasta_dict[header])
    gene = Gene(gene_length,gene_name,gene_num)
    gene.draw(context)

    #FIND EXON INFO: 
    exon_touple = (find_exon(fasta_dict[header]))
    exon_start = exon_touple[0]
    exon_len = exon_touple[1]
    gene_num = i
    exon = Exon(exon_start, exon_len, gene_num)
    exon.draw(context)


    #ITERATE THROUGH FASTA SEQ TO FIND MOTIFS: 
    #MAKE SEQ ALL UPPER CASE:
    seq = fasta_dict[header].upper()
    for j, motif in enumerate(MOTIF_DICT):
        matches = re.finditer(MOTIF_DICT[motif],seq)
        for match in matches:
            start = match.start()
            length = len(match.group())
            gene_number = i
            motif = Motif(start, length, gene_number, COLOR_LIST[j])
            motif.draw(context)

#########################
# CREATE TITLE + LEGEND #
#########################

# TITLE: 
title_height = 25
context.set_source_rgb(0,0,0)
context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
context.set_font_size(20)
context.move_to(350, title_height)
context.show_text("Motif Locations Within 4 Genes")

#LEGEND:
# for motif in motif_dict:
for j, motif in enumerate(MOTIF_DICT):
    context.set_source_rgb(COLOR_LIST[j][0], COLOR_LIST[j][1], COLOR_LIST[j][2])
    y = j*30 + 75
    x = 950
    length = 25
    height = feature_height
    context.rectangle(x,y,length,height)
    context.fill()
    title_height = y+10
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(15)
    context.move_to(820, title_height)
    context.show_text(motif)

#save to png 
name = fasta.split(".")
surface.write_to_png(f"{name[0]}.png")
