# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 10:09:25 2016

@author: lukasendler
"""
# back to the future
from __future__ import print_function

import sys, re
import os 
import argparse
import gzip

def read_fai(fai_file):
    """
    reads a fasta index file into a dictionary with chromosome naem -> length
    """
    chroms=dict()
    with open(fai_file, 'r') as f:
        for line in f:
            if (re.match("^\s*\#+",line)): # entry is comment/header
                continue
            entries=line.split("\t")
            chroms[entries[0]] = int(entries[1])
    return chroms


parser = argparse.ArgumentParser(description="""Go through a a fasta index file and a gtf file and remove all entries that are not within chromosome bounds""")

parser.add_argument("--gtf","-g", dest="gtffile", help="gtf-file, tries stdin if not set use \"STDIN\" or nothing for piping; default \"False\"", default=False)
parser.add_argument("--fai","-f", dest="faifile", help="fasta index file as created by eg. samtools", required=True)

args = parser.parse_args()
gtf_file = vars(args)['gtffile']
fai_file = vars(args)['faifile']

# read fai file
chroms = read_fai(fai_file)

# open vcf file
if not gtf_file or gtf_file == "STDIN":
    if not sys.stdin.isatty():
        inf = sys.stdin
    else:
        sys.exit("No gtf file or stdinput given")
else:
    if re.search("\.b?gz",gtf_file):
        inf = gzip.open(gtf_file,'rb')
    else:
        inf = open(gtf_file,"r")
        
for line in inf:
   if (re.match("^\s*\#+",line)): # entry is comment/header
       print(line,end="")
       continue
   entries = line.split("\t")
   if ( entries[0] in chroms ) and ( int(entries[4]) <= chroms[entries[0]] and int(entries[3]) >= 0 ):
       print(line,end="")

       