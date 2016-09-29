import sys, os, re
import argparse
parser = argparse.ArgumentParser(description='read an fst output file and create a tabular format only containing certain comparisons') 

parser.add_argument("--reps", dest="replicates", help="the replicates to obtain from the fst output (def. \'1:4,2:5,3:6\')",default='1:4,2:5,3:6')
parser.add_argument("--in", dest="infile", help="fst file to be read", required=True)
parser.add_argument("--app", dest="appendix", help="outfile appendix (optional)",default="")

args = parser.parse_args()
infile = vars(args)['infile']
replicates = vars(args)['replicates'].split(',')
appendix=vars(args)['appendix']
inf = open(infile,'r')
outfile=infile+".tab"+appendix
out=open(outfile,'w')
out_head="CHR\tBP\t"+"\t".join(replicates)
print >>out,out_head
for line in inf:
    line.rstrip()
    entries=line.split()
    outstring='{0}\t{1}'.format(entries[0],entries[1])
    fsts=[]
    for rep in replicates:
        for entry in entries[4:]:
            b=re.subn(rep+'\=','',entry)
            if (b[1] != 0):
                fsts.append(b[0])
                break
    outstring += "\t"+"\t".join(fsts)
    print >>out,outstring
out.close()
inf.close()

        
            
    
