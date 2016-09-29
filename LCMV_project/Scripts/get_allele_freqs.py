import sys
sys.path.append("/Volumes/Temp/Lukas/LCMV_project/Scripts/python")
import Sync_parser
import  os, re
import gzip
import numpy as np
from scipy.stats.stats import nanmean
import argparse
#from scipy import stats
#Author: Lukas Endler
parser = argparse.ArgumentParser(description='read a sync or cmhout file and create a file with the allele frequencies. the major/minor alleles are calculated  over all populations. the resulting file will contain the major allele freq and the coverage for each population.') 
parser.add_argument("-i","--infile", dest="infile", help="sync/cmhout file", required=True)
parser.add_argument("--all", dest="two",  action="store_false", help="output all allelefreqs and not just the major allele freq (calculated from just the major and minor counts) (default: True)", default=True)
parser.add_argument("-c","--cov", dest="cov",  action="store_true", help="output coverages too (default: False)", default=False)
parser.add_argument("-o","--out", dest="outfile", help="output sync file, if \"stdout\" redircet to STDOUT (default: infile.af)", default=None)
parser.add_argument("-p","--pV", dest="pV",  action="store_true", help="output P values (default: False)", default=False)
parser.add_argument("-a","--avg", dest="avg", help="create averages of afs in populations in each files eg. 1-3,4-6 (default: None)", default=None)
parser.add_argument("--pol", dest="pol", help="polarise afs. according rising allele according to population pairs, only effect if not averaging or outputting all allelefrwqs,  eg \"1:4,2:5,3:6\" (default: None)", default=None)
parser.add_argument("--ref", dest="ref", action="store_true", help="also print reference allele (default: False)", default=False)
parser.add_argument("--mc","--min_count", dest="min_count", help="integer indicating the minimal count over all populations above which to consider allels (default: 0)", default=0)
parser.add_argument("--mf", "--min_frac", dest="min_fract", help="float indicating the minimal allele frequency above which an allele is considered in each individual population (default: 0)", default=0)

args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
two = vars(args)['two']
pV = vars(args)['pV']
avg= vars(args)['avg']
ref= vars(args)['ref']
pol= vars(args)['pol']
cov= vars(args)['cov']
min_count= int(vars(args)['min_count'])
min_fract= float(vars(args)['min_fract'])


if pol:
    pol=[ x.split(":") for x in pol.split(",")]
    pol=[ [int(x[0])-1,int(x[1])-1]  for x in pol]

if avg:
    avg=[ x.split("-") for x in avg.split(",")]
    avg=[ [int(x[0])-1,int(x[1])-1]  for x in avg]

if outfile == None:
    outfile = infile + ".af"
if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
if outfile == "stdout" or outfile == "STDOUT":
    out=sys.stdout
else:
    out = open(outfile,"w")
nucs = [x for x in "ATCG"]
for line in inf:
    entry = Sync_parser.SyncLineParser(line,min_count=min_count,min_fract=min_fract)
    afs = entry.get_pop_allele_freqs(two)
    maj_min = "".join([ nucs[x] for x in  entry.get_two_major_alleles() ])
    all_pos=0 # 0: major, 1: minor
    if cov:
        coverages = entry.get_pop_coverages(two) 
    if pol:
        # see which direction of change is most frequent over replicates in the major allele (1: increasing,0,-1: decreasing)
        polarisation = np.median(np.sign(np.array([ afs[x[0]][0]-afs[x[1]][0]  for x in pol ])))
        #freqs = "\t".join([ ":".join([ str(y) for y in x ]) for x in afs ])
        if polarisation < 0:
            maj_min=maj_min[::-1]
            all_pos=1                
    if avg:
        # get the averages by first creating a matrix, averaging, and then recreating vectors
        if cov:
            coverages2 =  [ np.matrix(coverages[pos[0]:(pos[1]+1)]) for pos in avg]
            coverages =  np.array([ x.mean() for x in coverages2])
        #print afs
        try:
            afs2 = [ np.matrix(afs[pos[0]:(pos[1]+1)]) for pos in avg]
        except:
            print >> sys.stderr, "Problem with generating freq matrix for\n"+line+"\n"+str(afs)+"\nSkipping line"
            continue
        #print afs2
        afs = [ np.array(nanmean(x,0)).flatten() for x in afs2]
        #print afs
        #print avg
    if two:
        if cov: 
            freqs = "\t".join([ str(x[all_pos]) + "\t" + str(y) for x,y in zip(afs,coverages) ])
        else:
            freqs = "\t".join([ str(x[all_pos]) for x in afs ])
    else:
        freqs = "\t".join([ ":".join([ str(y) for y in x ]) for x in afs ])
    if pV and entry.cmhp:
        freqs += "\t"+str(entry.cmhp)
    if ref:
        print >>out, "{0}\t{1}\t{2}\t{3}\t{4}".format(entry.chr,entry.pos,entry.ref,maj_min,freqs)
    else:
        print >>out, "{0}\t{1}\t{2}\t{3}".format(entry.chr,entry.pos,maj_min,freqs)
inf.close()
out.close()

