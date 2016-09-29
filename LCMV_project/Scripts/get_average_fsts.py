import sys, os, re
import argparse
import numpy as np
parser = argparse.ArgumentParser(description='read an fst output file and create averages (mean and stdv) of the given comparisons ') 

parser.add_argument("--reps", dest="replicates", help="the replicates to calculate the fst averages for, multiple averages separated by \"|\" (def. \'1:4,2:5,3:6,7:10,8:11,9:12\')",default='1:4,2:5,3:6,7:10,8:11,9:12')
parser.add_argument("--in", dest="infile", help="fst file to be read", required=True)
parser.add_argument("--app", dest="appendix", help="outfile appendix (optional)",default="")
parser.add_argument("--bed","-b", dest="bed", action='store_true',help="outfile in bed graph format, only last mean (optional)",default=False)

args = parser.parse_args()
infile = vars(args)['infile']
replicates = vars(args)['replicates'].split('|')
appendix=vars(args)['appendix']
bed=vars(args)['bed']
replicates=[ x.split(",") for x in replicates]
out=sys.stdout
inf = open(infile,'r')
#outfile=infile+"_".join([ ".".join(x) for x in replicates ])+".avg.fsts"+appendix
#out=open(outfile,'w')
if bed:
    out_head="track\ttype=bedGraph"
else:
    out_head="CHR\tBP" + "\tMEAN\tSTDV"*len(replicates)
print >>out,out_head
field_nums=[]
for line in inf:
    line.rstrip()
    entries=line.split()
    if (len(field_nums) == 0):        
        for reps in replicates:
            fields = []
            for rep in reps:
                for i in range(5,len(entries)):
                    if (re.match(rep+'\=',entries[i])):
                        fields.append(i-5)
                        break
            field_nums.append(fields)
    fsts=np.array([ float(x.split('=')[1]) for x in entries[5:]])
    # for i in fields:
    #     fsts.append(float(entries[i].split('=')[1]))
    # fsts=np.array(fsts)
    means=""
    mean=0.0
    for fields in field_nums:
        mean = fsts[fields].mean()
        stdv = fsts[fields].std()
        means +="\t{}\t{}".format(mean,stdv)
    if bed:
        outstring='{0}\t{1}\t{2}\t{3}'.format(entries[0],int(entries[1])-1,entries[1],mean)
    else:
        outstring='{0}\t{1}\t{2}'.format(entries[0],entries[1],means)
    print >>out, outstring
out.close()
inf.close()
