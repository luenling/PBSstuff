######################### filtering a bam file by the ID's 

import pysam 
import sys,re
import os 
import collections
import argparse
# edited by Lukas to support greater numbers of reads (> 200 Mio)
#from Bio import trie
#########################################################   HELP   #########################################################################
parser = argparse.ArgumentParser(description="""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this. This version also uses the trie data structure from the Biopython package, to support bigger numbers of reads. This script takes a set of reads from a text file (one read ID per line) of reads to be filtered out and  as the input and splits a BAM file accordingly into a filtered and not filtered BAM which are subsets of the original BAM.

""") 
parser.add_argument("--input", dest="input", help="A BAM file",required=True)
parser.add_argument("--filter", dest="filt", help="reads to filtered out, can be bam or fastq/fasta also gzipped",required=True)
parser.add_argument("--out", dest="out", help="basename for output file (default input wo bam)", default=False)
parser.add_argument("--unmapped", dest="unmp", action="store_true", help="only get unfiltered reads with at least one mate unmapped", default=False)

options = parser.parse_args()

## read ID's from the fastq file
def get_ids(x):
    ids=list()
    bam=False
    if ( re.search(".*\.bam",x) ):
        bam=True
        filt_file=pysam.AlignmentFile(x,"rb")
    else:
        filt_file= pysam.FastxFile(x,persist=False)
    for entry in filt_file:
        if bam:
            id = entry.query_name
        else:
            id = entry.name
            if id.endswith('/1') or id.endswith('/2'):
                id = id[:-2]
        ids.append(id)
    filt_file.close()
    return set(ids)

print "reading "+ options.filt
filt_reads=get_ids(options.filt)
print "created set for filtering"

## index BAM file if necessary
if not os.path.exists(options.input+".bai"):
	print "indexing "+options.input
	os.system("samtools index "+options.input)

print "filtering " + options.input

samfile=pysam.AlignmentFile(options.input,"rb")
inf = options.input
# remove .bam at end
if options.input[-4:] == ".bam":
    inf = inf[:-4]

if not options.out:
    options.out = inf

filtout=pysam.AlignmentFile(options.out+"_filt.bam","wb",template=samfile)
cleanout=pysam.AlignmentFile(options.out+"_clean.bam","wb",template=samfile)
print "Filtered: "+options.out+"_filt.bam"
print "Clean: "+options.out+"_clean.bam"
### split BAM file
for l in samfile.fetch(until_eof=True):

	if l.query_name in filt_reads:
		filtout.write(l)
	else:
         if options.unmp and  not ( ( l.flag & 0x004 ) or ( l.flag & 0x008 ) ):
             continue
         cleanout.write(l)
filtout.close()
cleanout.close()
samfile.close() 
