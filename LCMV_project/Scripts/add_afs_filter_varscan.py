
import vcf
import sys, re
import os 
import argparse
from scipy import stats
import select
import gzip
import numpy as np
from collections import defaultdict,OrderedDict,namedtuple

def read_afs(filename):
    # read an allele frequency file in sync format
    inf = open(filename,"r")
    snp_dict=defaultdict(lambda: defaultdict( list ) )
    for line in inf:
        entries = line.rstrip().split()
        for i in range(3,len(entries)):
            freqs=np.array(entries[i].split(":"),dtype=float)
            snp_dict[entries[0]][int(entries[1])].append(freqs)
    inf.close()
    return snp_dict

def get_afs(snp_dict, chrom, bps,smp,alleles,nucs):
    # get allele freeqs from snp_dict
    return [ snp_dict[chrom][bps][smp][ nucs[x] ] for x in alleles ]
    
parser = argparse.ArgumentParser(description="""Go through a varscan created vcf file and add the alternative allele frequency to each sample
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file", default="STDIN")
parser.add_argument("--afs","-a", dest="afsfile", help="afs-file", default=False)
parser.add_argument("--minfreq",'-m', dest="minfreq",type=float, help="minimal minor allele frequency to consider (if not reached in at least one sample, allele will be removed)", default=0.01)
parser.add_argument("--mindepth",'-d', dest="mindepth",type=int, help="minimal sample depth to consider sample freq", default=75)
parser.add_argument("-v", dest="verb", action="store_true", help="verbose (default: FALSE)", default=False)
args = parser.parse_args()
vcf_file = vars(args)['vcffile']
afs_file = vars(args)['afsfile']
minfreq=vars(args)['minfreq']
mindepth=vars(args)['mindepth']
verb=vars(args)['verb']

# nucleotide positions as in popoolation2 sync format
nucs=dict(zip(["A","T","C","G"],[0,1,2,3]))
if afs_file != False:
    snp_dict=read_afs(afs_file)
#maxfreq=1-minfreq
# open vcf file
if vcf_file == "STDIN":
    inf = sys.stdin
    vcf_reader = vcf.Reader(inf)
else:
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    # if AF field not in formats, add it
if not ('AF' in vcf_reader.formats.keys()):
    newForm=vcf_reader.formats[vcf_reader.formats.keys()[0]]
    newForm=newForm._replace(id='AF', num=-1, type='Float', desc='Alternative allele frequency')
    vcf_reader.formats['AF'] = newForm

#vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
entries=[]
count=0
for entry in vcf_reader:    
    if verb:
        count += 1
        if count%500 == 0:
            sys.stderr.write(" processed " +str(count)+" variants at "+ str(entry.POS)+"\n" )
    #vcf_writer.write_record(entry)
    #entry.FORMAT=re.sub(":GL","",entry.FORMAT)
    if len(entry.REF) > 1:
        # cut all variants to SNPs
        entry.REF = entry.REF[0]
    entry.ALT=[ str(x)[0] for x in entry.ALT if str(x)[0] != entry.REF ]       
    afs=np.zeros(len(entry.ALT))
    acs=np.zeros(len(entry.ALT))
    b=list(entry.samples[0].data._fields)   
    if not re.search("(?<![^:])AF",entry.FORMAT):
        entry.add_format('AF')  
        b.append('AF')
    newdata=namedtuple('newdata',b)
    for snmr in range(0,len(entry.samples)):
        smp=entry.samples[snmr]
        nd=smp.data._asdict()
        if afs_file != False:
            nd['AF']=get_afs(snp_dict,entry.CHROM,entry.POS,snmr,entry.ALT,nucs)
        if nd['DP'] > mindepth:
            try:
                safs=np.array(nd['AF'],dtype=float,ndmin=1)
                safs = [ 0.0 if np.isnan(x) else x for x in safs]
                afs=np.maximum(afs,safs)
            except:
                print >>sys.stderr, str(afs) + " " + str(entry.POS)
                sys.exit()
        else:
            nd['AF']=[ "." for x in entry.ALT ]
        nd=newdata(**nd)        
        smp.data=nd
    keep = np.nonzero(afs >= minfreq)[0]
    if len(keep) == 0:
        continue
#    elif len(keep) == len(entry.ALT) or entry.ALT != list:
#        entries.append(entry)
#        continue
    if len(entry.ALT) > 1 and len(entry.ALT) > len(keep):
        for smp in entry.samples:
            nd=smp.data._asdict()
            for field in [ "AF" ]:
                #remove alternative alleles with too low a frequency
                nd[field]=[ nd[field][x] for x in keep ]
            try:
                nd['GT']=re.sub("[^0/.]","0",nd['GT'])
            except:
                pass
            nd=newdata(**nd)
            smp.data=nd
        # remove alternative alleles with too low a freq
        entry.ALT=[ entry.ALT[x] for x in keep]
    if len(keep) > 0:
       entries.append(entry)

vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
for entry in entries:
    vcf_writer.write_record(entry)
#sys.exit(0)

    
            
        
#if not entry.is_snp or len(entry.alleles) > 2:
#    # not snp or more than two alleles
#    #continue
##   
#alleles = "".join([ str(x) for x in entry.alleles])
#if maf:
#    # check if SNP in dict and if same major/minor
#    if not( entry.CHROM in maf_dict.keys() and entry.POS in maf_dict[entry.CHROM].keys()):
#        continue
#    if not ( alleles[0] in maf_dict[entry.CHROM][entry.POS][0] and alleles[1] in maf_dict[entry.CHROM][entry.POS][0]):
#        continue
#    # check if major is reference allele:
#    if alleles != maf_dict[entry.CHROM][entry.POS][0]:
#        maf_dict[entry.CHROM][entry.POS][1]=1-maf_dict[entry.CHROM][entry.POS][1]
#homR=0 # homozygous reference
#homA=0 # homozygous alternative
#nogt=0 # not called
#for sam in samples:
#    gt=entry.genotype(sam).data.GT
#    if not gt:
#        nogt+=1
#    elif gt == '0/0':
#        homR+=1
#    elif gt == '1/1':
#        homA+=1
#tot_cal=len(samples)-nogt # total number of called samples
#if tot_cal == 0:
#    continue
#het=tot_cal-(homR+homA)
#fR=(homR+het/2.0)/tot_cal
#fA=(homA+het/2.0)/tot_cal
#assert abs(1 -(fR + fA)) < 1e-3, "fR and fA do not sum to 1 at"+ entry.CHROM+"\t"+str(entry.POS)
## reduction of heterozygosity from hardy weinberg, F coefficient of inbreeding
#if het == 0:
#    Fi = 1.0
#else:
#    Fi=1.0-(float(het)/tot_cal)/(2.0*fR*fA)
#tot_pop=""
#if maf:
#    # calculate population total
#    fRt=maf_dict[entry.CHROM][entry.POS][1]
#    fAt=1-fRt
#    HomExp=tot_cal*(fRt**2+fAt**2)
#    HetExp=tot_cal*(2*fRt*fAt)
#    try:
#        Fit=1.0-(float(het)/tot_cal)/(2.0*fRt*fAt)
#    except:
#        Fit=1.0
#    tot_pop="\t{:.2}\t{:.2}\t{:.3}\t{:.3}\t{:.3}".format(HomExp,HetExp,fRt,fAt,Fit)
#print entry.CHROM+"\t"+str(entry.POS)+"\t"+alleles+"\t"+"{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}".format(tot_cal,homR,het,homA,Fi,fR,fA)+tot_pop




    

