
def read_snp_file(infile,both=False):
    """
    reads a cmh or sync like file and creates a dict with chrom->set(bps) or chrom->bps->line
    """
    if (both):
        snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    else:
        snp_dict=defaultdict(set)  #dictionary of chroms with positions and pVals of snps
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        if both:
            snp_dict[fields[0]][fields[1]]="\t".join(fields[2:])
        else:
            snp_dict[fields[0]].add(fields[1])        
    inf.close()
    return snp_dict


import re
#import numpy as np
import gzip

#from scipy import stats
from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description="""
reads snps from file A and only gets the lines with the same positions from file B. writes to stdout.
if both print line from A and attach B minus the first three fields
if -m, add a 1 or zero at the end of line of file B if the coordinates exist/do not exist in A.
if -d only give the lines in B that do not exist in A
""") 

parser.add_argument("--inA","-a", dest="inA", help="file with snps, can be a (b)gzipped vcf or sync file", required=True)
parser.add_argument("--inB","-b", dest="inB", help="sync or cmh file to get positions from", required=True)
parser.add_argument("--both", dest="both", action="store_true", help="print CHR,POS,rest from A and rest from B ", default=False)
parser.add_argument("-d", dest="diff", action="store_true", help="only get lines in B that are not in A", default=False)
parser.add_argument("-m", dest="mark", action="store_true", help="if line in A add 1 to B else add 0", default=False)
parser.add_argument("--header", dest="header", action="store_true", help="print comment lines in file B", default=False)

args = parser.parse_args()
inA = vars(args)['inA']
inB = vars(args)['inB']
both = vars(args)['both']
diff = vars(args)['diff']
mark = vars(args)['mark']
header = vars(args)['header']

snp_dict=read_snp_file(inA,both)
if re.search("\.b?gz",inB):
    inf = gzip.open(inB,'rb')
else:
    inf = open(inB,"r")
for line in inf:
    if re.match("\#",line):
        if header:
            line=line.rstrip()
            print line
            continue
        else:
            continue
    line=line.rstrip()
    fields=line.split()
    if fields[0] == "":
        continue
    if diff:
       if not ( fields[0] in snp_dict and fields[1] in snp_dict[ fields[0] ]):
           print line
    elif ( fields[0] in snp_dict and fields[1] in snp_dict[ fields[0] ]):
        if both:
            print "\t".join([ fields[0],fields[1],snp_dict[fields[0]][fields[1]] ]),"\t", "\t".join(fields[3:])
        elif mark:
            print line+"\t1"
        else:
            print line
    elif mark:
        print line+"\t0"
inf.close()

