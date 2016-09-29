# usage python fst_calc.py -c trial_r1_r2_r3_b1_b2_b3.cmhout -s 1,2,3 -n 4,5,6 > 2.txt

# -s: serial number of selected populations, e.g -s. 1,2,3
# -n: serial number of not-selected populations (or base populations), e.g. 4,5,6
# important suggestion: pls input all (not-)selected populations one by one, e.g, 1,2,3,4,5 or 4,5,6 
# an order of 1,3,5 is a disaster for this cscript because you will get a wrong estimation of fst2............

# output file: 2.txt, order: chrom	positon	main_allele	fst1	fst2
# fst1: fst1= (Ht-Hs)/Ht, where Ht is the total heterozygosity and Hs is heterozygosity of subpopulations.
# here subpopulaitons are all selected populations (as one subpopulation) and all base populations (as another subpopulation)
# fst2: fst2= 1/k * sigma (fst(i)), here fst(i) is pairwise fst calculated from a base-selected population pair.
  
from optparse import OptionParser
import random

parser=OptionParser()

# import name of sync file
parser.add_option('-c', dest='cmhfile')
parser.add_option('-s', dest='selectedpop')
parser.add_option('-n', dest='notselectedpop')
# import name of coverage file
#parser.add_option('-c', dest='coverage')

(options, args) = parser.parse_args()

fcmh=options.cmhfile
ser_num_s=options.selectedpop
ser_num_n=options.notselectedpop
#nc=options.coverage


def cal_fst_het(l_a_f, l_s_p, l_n_p):
    # l_allele_freq, l_selected_pop, l_notselected_pop
    
    if len(l_s_p) != len(l_n_p):
        print 'Error! The number of selected populations is not equal to the number of notselected populations!'
        return 
    '''
    print l_a_f, l_s_p, l_n_p
    ['2L', '5390', 'T', '8:45', '13:58', '12:36', '26:65', '63:186', '75:217', '0.09577668'] ['1', '2', '3'] ['4', '5', '6']
    input()
    '''
    ######################## fst1
    l_af_s=[l_a_f[2+int(x)] for x in l_s_p]
    l_af_n=[l_a_f[2+int(x)] for x in l_n_p]
    l_s=[]
    l_n=[]
    for x in l_af_s:
        t_a=0
        t_x=x.split(':')
        t_all=int(t_x[0])+int(t_x[1])
        if t_all != 0:
            l_s.append((int(t_x[0])*1.0)/(t_all*1.0))
            # l_s.append((int(t_x[1])*1.0)/(t_all*1.0))
    for x in l_af_n:
        t_a=0
        t_x=x.split(':')
        t_all=int(t_x[0])+int(t_x[1])
        if t_all != 0:
            l_n.append((int(t_x[0])*1.0)/(t_all*1.0))
    
    p_s=0.0
    for x in l_s:
        p_s=p_s+x
    p_s=p_s/(len(l_s)*1.0)
    
    p_n=0.0
    for x in l_n:
        p_n=p_n+x
    p_n=p_n/(len(l_n)*1.0)    
    
    p_t=(p_s+p_n)/2.0
    # (1.0-p_n*p_n-(1.0-p_n)*(1.0-p_n))
    fst1=((1.0-p_t*p_t-(1.0-p_t)*(1.0-p_t))-0.5*((1.0-p_s*p_s-(1.0-p_s)*(1.0-p_s)))-0.5*((1.0-p_n*p_n-(1.0-p_n)*(1.0-p_n))))/(1.0-p_t*p_t-(1.0-p_t)*(1.0-p_t))
    
    ################### fst2
    num_pair=len(l_s_p)
    start_s=int(l_s_p[0])
    start_n=int(l_n_p[0])
    fst2 = 0.0

    for i in range(num_pair):
        vafn=l_a_f[2+i+start_n]
        vafs=l_a_f[2+i+start_s]
        l_f_locus_n=vafn.split(':')
        l_f_locus_s=vafs.split(':')
        #print l_f_locus_n, l_f_locus_s
        fp1=0.5*float(l_f_locus_n[0])/(float(l_f_locus_n[0])+float(l_f_locus_n[1]))+0.5*float(l_f_locus_s[0])/(float(l_f_locus_s[0])+float(l_f_locus_s[1]))
        fpp1=float(l_f_locus_n[0])/(float(l_f_locus_n[0])+float(l_f_locus_n[1]))
        fpp2=float(l_f_locus_s[0])/(float(l_f_locus_s[0])+float(l_f_locus_s[1]))
        #print fp1, fpp1, fpp2
        #raw_input()
        ht2=1.0-fp1*fp1-(1.0-fp1)*(1.0-fp1)
        hs2=0.5*(1.0-fpp1*fpp1-(1.0-fpp1)*(1.0-fpp1))+0.5*(1.0-fpp2*fpp2-(1.0-fpp2)*(1.0-fpp2))
        if ht2 == 0:
            fst2=fst2+0.0
        else:
            fst2=fst2+(1.0 - hs2/ht2)
    fst2=fst2/(float(num_pair))
    
    vt=[]
    vt.append(fst1)
    vt.append(fst2)
    
    return vt
    
def filter_af(lc):
    lc_1=lc[:]
    lc_t=lc_1[3:9]
    #print lc_t
    #input()
    l_0=[0]*6
    for x in lc_t:
        lx=x.split(':')
        ly=[int(y) for y in lx]
        for i in range(len(l_0)):
            l_0[i]=l_0[i]+ly[i]
    i=0
    for x in l_0:
        if x > 0:
            i=i+1
    if i <= 2:
        return lc_1
    else:
        lc_af=[]
        # lc_t: ['0:0:50:2:0:0', '0:0:61:2:0:0', '0:0:34:5:0:0', '0:1:95:2:0:0', '0:7:293:15:0:0', '1:1:313:19:0:0']
        for x in lc_t:
            x1=x.split(':')
            j = 0
            for x2 in x1:
                j = j + int(x2)
            if j != 0:
                lt1=[]
                for x in x1:
                    lt1.append(int(x)/float(j))
                lc_af.append(lt1)
            else:
                lc_af.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        # here is a vector to judge main alleles
        l_j=[0.0]*6
        for x in lc_af:
            for i1 in range(6):
                l_j[i1] = l_j[i1] + x[i1]
        vector_main = []
        l_sort=l_j[:]
        l_j.sort()
        
        for i in range(len(l_sort)):
            if l_j[-1] == l_sort[i]:
                vector_main.append(i)
                break
        for i in range(len(l_sort)):
            if l_j[-2] == l_sort[i]:
                vector_main.append(i)
                break        
        '''
        print l_sort, l_j, vector_main
        [0.0029940119760479044, 0.03542031583092319, 5.638259535351066, 0.3233261368419623, 0.0, 0.0]
        [0.0, 0.0, 0.0029940119760479044, 0.03542031583092319, 0.3233261368419623, 5.638259535351066]
        [2, 3]
        raw_input()
        '''
        for i in range(3,len(lc_1)-1):
            #print lc_1
            x01=lc_1[i]
            x02=x01.split(':')
            #print x02, vector_main
            x03=['0']*6
            x03[vector_main[0]]=x02[vector_main[0]]
            x03[vector_main[1]]=x02[vector_main[1]]
            y01=':'.join(x03)
            lc_1[i] = y01

        return lc_1
       
def transfer(lc):
    '''
    print lc
    #['2L', '75770', 'C', '0:0:50:2:0:0', '0:0:61:2:0:0', '0:0:34:5:0:0', '0:1:95:2:0:0', '0:7:293:15:0:0', '1:1:313:19:0:0', '0.461779']
    input()
    '''
    lc=filter_af(lc)
    
    l_0=[0]*6
    #print len(lc) # 24

    for x in lc[3:(len(lc)-4+3)]:
        lx=x.split(':')
        ly=[int(y) for y in lx]
        for i in range(len(l_0)):
            l_0[i]=l_0[i]+ly[i]
    '''
    print l_0
    [16, 0, 5184, 0, 0, 0]
    '''
    l_1=[]
    for i in range(len(l_0)):
        if l_0[i] != 0:
            l_1.append(i)
    #print l_1 [0,2]

    l_af=lc[0:3]
    if len(l_1)  == 0:
        l_af.append('0:0')
    elif len(l_1) == 1:
        for i in range(3, (len(lc)-4+3)):
            for x in lc[i].split(':'):
                if int(x) > 0:
                    ap= '0'+':'+x
                    l_af.append(ap)
    else:
        for i in range(3, (len(lc)-4+3)):
            lft=lc[i].split(':')
            ap=lft[l_1[0]]+':'+lft[l_1[1]]
            l_af.append(ap)
    l_af.append(lc[-1])
    '''
    print l_af
    ['2L', '6580', 'C', '0:260', '0:260', '0:260', '0:260', '2:258', '0:260', '2:258',
     '0:260', '0:260', '0:260', '0:260', '0:260', '2:258', '0:260', '4:256', '0:260',
     '2:258', '0:260', '4:256', '0:260', '0.0001717933']
    raw_input()
    '''
    
    return l_af


if __name__== "__main__":
    readcmh=open(fcmh,'r')
    print 'chrom\tpositon\tmain_allele\tfst1\tfst2'
    l_selected_pop=ser_num_s.split(',')
    l_notselected_pop=ser_num_n.split(',')
    lcmh=readcmh.readline().split()
    while lcmh:
        
        '''
        print lcmh
        ['2L', '6580', 'C', '0:0:260:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '2:0:258:0:0:0', '0:0:260:0:0:0', '2:0:258:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '0:0:260:0:0:0', '2:0:258:0:0:0', '0:0:260:0:0:0', '4:0:256:0:0:0', '0:0:260:0:0:0', '2:0:258:0:0:0', '0:0:260:0:0:0', '4:0:256:0:0:0', '0:0:260:0:0:0', '0.0001717933']
        raw_input()
        '''
        l_allele_freq=transfer(lcmh)
        #diffstat=cal_diffstat(l_allele_freq, l_selected_pop, l_notselected_pop)
        #assocstat=cal_assocstat(l_allele_freq, l_selected_pop, l_notselected_pop)
        fst_het_vector=cal_fst_het(l_allele_freq, l_selected_pop, l_notselected_pop)
        fst1=fst_het_vector[0]
        fst2=fst_het_vector[1]
        #hs_over_d=fst_het_vector[1]
        #hs_over_hc=fst_het_vector[2]
        
        loutput=lcmh[0:3]
        #loutput.append(str(diffstat))
        #loutput.append(str(assocstat))
        loutput.append(str(fst1))
        loutput.append(str(fst2))
        #loutput.append(str(hs_over_d))
        #loutput.append(str(hs_over_hc))
        #loutput.append(lcmh[-1])
        
        lcmh=readcmh.readline().split()
        print '\t'.join(loutput)
        
        # print lcmh
        # raw_input()
        # print diffstat, assocstat, fst, hs_over_d, hs_over_hc, lcmh[-1]
        
        
    readcmh.close
