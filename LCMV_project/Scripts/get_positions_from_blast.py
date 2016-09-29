def get_coords(pos,start,stop,subjct):
    # gets position to retrieve and 
    strand=1 # if reverse: -1
    idx=0
    idcs=[False,False,False]
    if start >= stop:
        strand=-1
        idx=-1 # reverse strand
        a=start # exchange start and stop values
        start=stop
        stop=a
    rpos=start-1
    while (rpos <= stop):
        #print str(rpos)+" "+str(idx)+" "+ subjct[idx]
        if subjct[idx] != "-":
            rpos += 1
        if rpos == pos[0]:
            idcs[0] = idx
        elif rpos == pos[1]:
            idcs[1] = idx
            break
        idx += strand
    idcs[1] += strand # to get whole intervall
    idcs[2]=strand
    return idcs


# this script should get all positions aligned to a certain position of a subject sequence
from Bio.Blast import NCBIXML
import sys
import os 
import collections
import argparse
from string import maketrans

parser = argparse.ArgumentParser(description="""
 this script should get all positions aligned to a certain position of a subject sequence
""") 

parser.add_argument("--blast","-b", dest="bfile", help="xml file with blast results", required=True)
parser.add_argument("--subjpos","-p", dest="subjpos",  help="subject postions eg. \"3R:17064230-17064232\"", required=True)
args = parser.parse_args()
bfile = vars(args)['bfile']
( subjct, pos)=vars(args)['subjpos'].split(":")
pos = [int(x) for x in pos.split("-")]
bfh = open(bfile)
blast_handles=NCBIXML.parse(bfh)
trantab=maketrans("ACGT","TGCA")

for blast_handle in blast_handles:
    for algn in blast_handle.alignments:
        if algn.accession != subjct:
            continue
        for hsp in algn.hsps:
            # check if in hsp intervall
            if not (pos[0] >= min(hsp.sbjct_end,hsp.sbjct_start) and pos[1] <=  max(hsp.sbjct_end,hsp.sbjct_start)):
                continue
            idcs=get_coords(pos,hsp.sbjct_start,hsp.sbjct_end,hsp.sbjct)
            sbj_seq=str(hsp.sbjct[idcs[0]:idcs[1]:idcs[2]])
            q_seq=str(hsp.query[idcs[0]:idcs[1]:idcs[2]])
            match=str(hsp.match[idcs[0]:idcs[1]:idcs[2]])
            if idcs[2] <= 0: # reverse strand
                sbj_seq=sbj_seq.translate(trantab)
                q_seq=q_seq.translate(trantab)
            print "*******"
            print "{}\t{}\t{}\t{}".format(pos[0],sbj_seq,pos[1],algn.accession)
            print "{}\t{}".format(" "*len(str(pos[0])),match)
            print "{}\t{}\t{}\t{}".format(" "*len(str(pos[0])),q_seq," "*len(str(pos[1])),blast_handle.query)
            print"*******"
            
    

