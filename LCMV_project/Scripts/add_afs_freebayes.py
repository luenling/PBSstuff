
import vcf
import sys, re
import os 
import argparse
from scipy import stats
import select
import gzip
import numpy as np
from collections import defaultdict,OrderedDict,namedtuple

parser = argparse.ArgumentParser(description="""Go through a freebayes created vcf file and add the alternative allele frequency array to each sample
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file", required=True)
parser.add_argument("--minfreq","-m", dest="minfreq",type=float, help="minimal minor allele frequency to consider (if not reached in at least one sample, allele will be removed)", default=0.01)
parser.add_argument("--sb", dest="minsb",type=float, help="minmal ratio of supporting strands", default=0.1)
parser.add_argument("--fsb", dest="maxfsb",type=float, help="maximal FSB", default=30)
parser.add_argument("--rpp", dest="maxrpp",type=float, help="maximal RPP", default=30)
parser.add_argument("--mc", dest="mc",type=float, help="min. supporting reads", default=3)
parser.add_argument("-v", dest="verb", action="store_true", help="verbose (default: FALSE)", default=False)
args = parser.parse_args()
vcf_file = vars(args)['vcffile']
minfreq=vars(args)['minfreq']
minsb=vars(args)['minsb']
maxfsb=vars(args)['maxfsb']
maxrpp=vars(args)['maxrpp']
mc=vars(args)['mc']
verb=vars(args)['verb']
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
if not ('FSB' in vcf_reader.infos.keys()):
    newForm=vcf_reader.infos[vcf_reader.infos.keys()[0]]
    newForm=newForm._replace(id='FSB', num=-1, type='Float', desc='PHRED scaled Chi2 test ( or if expected < 50 > 5 with Yates correction or if minimum expected <=5 Fisher test) value for strand bias')
    vcf_reader.infos['FSB'] = newForm
if not ('SB' in vcf_reader.infos.keys()):
    newForm=vcf_reader.infos[vcf_reader.infos.keys()[0]]
    newForm=newForm._replace(id='SB', num=-1, type='Float', desc='ratio of less common strand to more common strand of reads supporting alternative alleles')
    vcf_reader.infos['SB'] = newForm

#vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
entries=[]
count=0
for entry in vcf_reader:
    if verb:
        count += 1
        if count%500 == 0:
            sys.stderr.write(" processed " +str(count)+" variants at "+ str(entry.POS)+"\n" )
    #vcf_writer.write_record(entry)
    entry.FORMAT=re.sub(":GL","",entry.FORMAT)
    entry.add_format('AF')
    afs=np.zeros(len(entry.ALT))
    acs=np.zeros(len(entry.ALT))    
    b=list(entry.samples[0].data._fields)
    b.remove('GL')
    b.append('AF')
    newdata=namedtuple('newdata',b)
    for smp in entry.samples:
        nd=smp.data._asdict()
        if type(nd['DP']) == int:
            nd['AF']=np.array(nd['AO'],dtype=float,ndmin=1)
            nd['AF']= [ x/nd['DP'] for x in nd['AF'] ]
            try:
                afs=np.maximum(afs,nd['AF'])
                acs=np.maximum(acs,np.array(nd['AO'],dtype=float,ndmin=1))
            except:
                print >>sys.stderr, str(afs) + " " + str(entry.POS)
                sys.exit()
        else: # none type, not called
            #print nd
            nd['AF']=nd['AO']
        del(nd['GL'])
        nd=newdata(**nd)        
        smp.data=nd
    keep = np.nonzero(np.logical_and(afs >= minfreq,acs >= mc))[0]
    if len(keep) == 0:
        continue
#    elif len(keep) == len(entry.ALT) or entry.ALT != list:
#        entries.append(entry)
#        continue
    if len(entry.ALT) > 1 and len(entry.ALT) > len(keep):
        for smp in entry.samples:
            nd=smp.data._asdict()
            #print >>sys.stderr, str(nd)
            if type(nd['DP']) != int:
                continue
            for field in ( "AO","QA","AF" ):
                #remove alternative alleles with too low a frequency
                nd[field]=[ nd[field][x] for x in keep ]
            nd=newdata(**nd)
            smp.data=nd
        # remove alternative alleles with too low a freq
        entry.ALT=[ entry.ALT[x] for x in keep ]
        for field in entry.INFO.keys():
            # remove low freq alleles from info fields
            if type(entry.INFO[field]) == list and len(entry.INFO[field]) == len(afs):
                entry.INFO[field] = [entry.INFO[field][x] for x in keep]
        entry.INFO["NUMALT"]=len(keep)
    FSBs=[]
    SBs=[]
    for i in range(0,len(entry.INFO['SAF'])):
        try:
            matr=np.array([[entry.INFO['SRF'],entry.INFO['SRR']],
                                     [entry.INFO['SAF'][i],entry.INFO['SAR'][i]]])
            minexp=np.min(matr.sum(axis=0))*np.min(matr.sum(axis=0))/matr.sum()
            if minexp > 50:
                FS = stats.chi2_contingency(matr)[1]    
            elif minexp > 5:
                FS = stats.chi2_contingency(matr,correction=True)[1]    
            else:
                #print >>sys.stderr, str(matr)
                FS = stats.fisher_exact(matr)[1]
            FS = -10*np.math.log10(FS)
        except:
            FS = 0.0
        FSBs.append(FS)
        try:
            SB=float(np.min([entry.INFO['SAF'][i],entry.INFO['SAR'][i]]))/np.max([entry.INFO['SAF'][i],entry.INFO['SAR'][i]])
        except:
            SB=np.NaN
        SBs.append(SB)    
    entry.INFO['FSB']=FSBs
    entry.INFO['SB']=SBs
    FSBs=np.array(FSBs,dtype=float,ndmin=1)
    SBs=np.array(SBs,dtype=float,ndmin=1)
    keep = np.nonzero(np.logical_or(FSBs <= maxfsb,SBs >= minsb))[0]
    if len(keep) == 0:
        continue
    if len(entry.ALT) > 1 and len(entry.ALT) > len(keep):
        for smp in entry.samples:
            nd=smp.data._asdict()
            if type(nd['DP']) != int:
                continue
            for field in ("AO","QA","AF"):
                #remove alternative alleles with too low a frequency
                nd[field]=[ nd[field][x] for x in keep ]
            nd=newdata(**nd)
            smp.data=nd
        # remove alternative alleles with too low a freq
        entry.ALT=[ entry.ALT[x] for x in keep]
        for field in entry.INFO.keys():
            # remove low freq alleles from info fields
            if type(entry.INFO[field]) == list and len(entry.INFO[field]) == len(afs):
                entry.INFO[field] = [entry.INFO[field][x] for x in keep]
        entry.INFO["NUMALT"]=len(keep)
    #print entry.INFO['FSB']
    #vcf_writer.write_record(entry)
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




    

