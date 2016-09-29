# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:44:55 2015

@author: lukasendler
"""

def get_reads_and_nucs_new(chrom,bps,snp_indx,bamfh,read_dict,map_qual=20,base_qual=20):
    """
    function for the new pysam synthax (pysam 0.8.0 and up)
    gets a SNP position (chrom & bps), SNP indx, a bam file handle and returns the numbers of the fw and reverse strands for each nucleotide found    
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101) 
    # have to go through loop, I think
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups
            break
    #b=[x.alignment.seq[x.qpos] for x in pile_col ]
    #counts=[ b.count(x) for x in ["A","C","G","T"]]
    #read_dict=defaultdict(list)
    snp_nucs=[]
    
    for pile_read in pile_col:
        # check whether matched, mapping qual and base qual alright
        if  pile_read.is_del or (pile_read.alignment.mapping_quality < map_qual) or pile_read.alignment.query_qualities[pile_read.query_position] < base_qual :
            continue
        # assign read_read
        read_dict[pile_read.alignment.query_name].append([snp_indx,pile_read.alignment.query_sequence[pile_read.query_position]])
        snp_nucs.append(pile_read.alignment.query_sequence[pile_read.query_position])
    return np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)

# get forward/reverse: pile_read.alignment.is_reverse
# get postion in read: pile_col.query_position
# from end: a.is_tail or a.alignment.alen - pile_col.query_position
