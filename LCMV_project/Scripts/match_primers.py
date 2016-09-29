# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 13:11:33 2015

@author: lukasendler
"""

def find_substr(v,ps,pid,strand):
    vs = str(v.seq)
    pos=vs.find(ps)
    if pos != -1:
        print "{}\t{}\t{}\t{}\t{}".format(v.id,pos,pos+len(ps),pid,strand)

    
import sys
from Bio import SeqIO
#fasta1="/Volumes/vetgrid01/LCMV_project/References/viruses.fasta"
#fasta2="/Volumes/vetgrid01/LCMV_project/References/primers.fna"
fasta1=sys.argv[1]
fasta2=sys.argv[2]


viruses = [ x for x in SeqIO.parse(fasta1,"fasta") ]
primer = [ x for x in SeqIO.parse(fasta2, "fasta")]

for p in primer:
    pf = str(p.seq)
    pr = str(p.reverse_complement().seq)
    for v in viruses:
        find_substr(v,pf,p.id,"+")
        find_substr(v,pr,p.id,"-")


