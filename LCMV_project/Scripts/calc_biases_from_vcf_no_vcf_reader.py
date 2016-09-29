#import vcf
import sys, re
import os 
import numpy as np
from scipy import stats
import argparse
import gzip

def welch_t_test(a, b):
    """
    gets two populations given by number, sum, and sum of squares and performs Welch's T-test (assuming independent variances).
    gives the t, degrees of freedom and the pValue for the mean of b being greater than the mean of a (onesided t test)
    """
    pop1=np.array(a,dtype=float)
    pop2=np.array(b,dtype=float)
    (n,u,v)=(0,1,2)
    if pop1[n] <= 1 or pop2[n] <= 1 or pop1[n] + pop2[n] < 3:
        # not enough reads to say much
        return np.array([0,0,1.0])
    pop1[1]=pop1[1]/pop1[0] # get population1 mean
    pop2[1]=pop2[1]/pop2[0] # get population2 mean
    if pop1[u] <= pop2[u]:
        # alternative allele has greater tail distance mean than reference
        return np.array([0,0,1.0])
    pop1[2]=(pop1[2]-pop1[0]*pop1[1]**2)/(pop1[0]-1) # get pop1 variance
    pop2[2]=(pop2[2]-pop2[0]*pop2[1]**2)/(pop2[0]-1) # get pop2 variance
    # calculate t value
    t = (pop1[u] - pop2[u])/ np.math.sqrt(pop1[v]/pop1[n] + pop2[v]/pop2[n])
    # calculate v (degree of freedoms)
    nd = (pop1[v]/pop1[n] + pop2[v]/pop2[n])**2/(pop1[v]**2/(pop1[n]**2*(pop1[n]-1)) + pop2[v]**2/(pop2[n]**2*(pop2[n]-1)) )
    pV = 1 - stats.t.cdf(t,nd) # onesided, Null Hyp: pop2 has higher or equal tail distance
    return np.array([t,nd,pV])
    

parser = argparse.ArgumentParser(description="""Go through a VCF file and add fisher strand bias (FS), ratio of fw to rev (SB), mean tail distance for the ref and alt allele (RTD and ATD),tail distance bias (TDB), and the absolute value of the RPB calculated by samtools (t test that alt allele more likely to lie in the last or first 11 bases than the ref).  
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file, tries stdin if not set use \"STDIN\" or nothing for piping; default \"False\"", default=False)

args = parser.parse_args()
vcf_file = vars(args)['vcffile']
# open vcf file
if not vcf_file or vcf_file == "STDIN":
    if not sys.stdin.isatty():
        inf = sys.stdin
    else:
        sys.exit("No vcf file or stdinput given")
else:
    if re.search("\.b?gz",vcf_file):
        inf = gzip.open(vcf_file,'rb')
    else:
        inf = open(vcf_file,"r")

print "#CHR\tBPS\tALL\tFS\tSB\tRTD\tATD\tTDB\taRPB"

for line in inf:
    if (re.match("^\s*\#+",line)): # entry is comment/header
        continue
    entries = line.split("\t")
    alleles = entries[3]+entries[4][0]
    if not alleles[1] in ["A","C","T","G"]:
        continue
    INFO = {} # need fields I16 and RPB
    for entry in entries[7].split(";"):
        (field_id,values)=entry.split("=")
        if field_id == "I16":
            INFO[field_id] = [ x for x in values.split(",")]
        if field_id == "RPB":
            INFO[field_id] = float(values)
    read_nums=np.array(INFO['I16'][0:4],dtype=float)
    read_nums=read_nums.reshape((2,2))
    try:
        minexp=np.min(read_nums.sum(axis=0))*np.min(read_nums.sum(axis=0))/read_nums.sum()
        if minexp > 50:
            FS = stats.chi2_contingency(read_nums)[1]    
        elif minexp > 5:
            FS = stats.chi2_contingency(read_nums,correction=True)[1]
        elif read_nums.mean > 1e4:
            FS = stats.chi2_contingency(read_nums,correction=True)[1]
        else:
            #print >>sys.stderr, str(matr)
            FS = stats.fisher_exact(read_nums)[1]
        FS = -10*np.math.log10(FS)
    except:
        FS = 0.0
    (n1,n2)=read_nums.sum(axis=1)
    try:
        SB=float(read_nums[1].min())/read_nums[1].max()
    except:
        SB=np.NaN    
    try:
        RPB=INFO['RPB']
    except:
        RPB=np.NaN
    pop1=[n1,float(INFO['I16'][12]),float(INFO['I16'][13]) ]
    if n1 > 0:
        RTD=pop1[1]/pop1[0]
    else:
        RTD=np.NaN
    pop2=[n2,float(INFO['I16'][14]),float(INFO['I16'][15]) ]
    if n2 > 0: 
        ATD=pop2[1]/pop2[0]
    else:
        ATD = np.NaN
    try:
        (t,nd,TDB)=welch_t_test(pop1, pop2)
        if TDB == 0.0:
            TDB = 1e-20
        TDB=(-10*np.math.log10(TDB))
    except:
        TDB=0.0
    print "{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{}\t{:.3f}\t{:.3f}".format(entries[0],entries[1],alleles,FS,SB,RTD,ATD,TDB,abs(RPB))
