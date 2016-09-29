# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:13:18 2016

@author: lukasendler
"""

from Bio import SeqIO
import sys,re
import os 
import collections
import argparse
#########################################################   HELP   #########################################################################
parser = argparse.ArgumentParser(description="""
H E L P:
____________

Needs Biopython to be installed. Filters out sequences of a length below a certain threshold 
""") 
parser.add_argument("--input","-i", dest="input", help="fasta file input, default: STDIN",default=False)
parser.add_argument("--len", "-l",  dest="len", help="minimal length to be kept, integer, default: 350",default="350")
parser.add_argument("--out", dest="out", help="basename for output file (default SDTOUT)", default=False)

options = parser.parse_args()

if options.input ==  False:
    input = sys.stdin
else:
    input = options.input

length = int(options.len)
out_seqs = []

for seq_record in SeqIO.parse(input, "fasta"):
    if len(seq_record.seq) < length:
        continue
    out_seqs.append(seq_record)

if options.out == False:
    outf = sys.stdout
else:
    outf = open(options.out, "w")
SeqIO.write(out_seqs,outf, "fasta")
   
    