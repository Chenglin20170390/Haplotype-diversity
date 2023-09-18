#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author   : WN
# @Time     : 20230811
# @File     : Design_less_dSVs_in_A_line_IPH.py
# @Project  : scripts/Designing Ideal potato haplotypes (IPH)

import itertools
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Find best groups for ideal haplotype for A lines', add_help = False, usage = '\npython3 -chr [chr01] -out [output.file]')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-chr', '--chr', metavar = '[chr01]', help = 'input', required = True)
required.add_argument('-out', '--out', metavar = '[output.txt]', help = 'output', required = True)
optional.add_argument('-h', '--help', action = 'help', help = 'help')
args = parser.parse_args()

#DM reference genome
chromosomes={'chr01':88591686,'chr02':46102915,'chr03':60707570,'chr04':69236331,'chr05':55599697,'chr06':59091578,'chr07':57639317,'chr08':59226000,'chr09':67600300,'chr10':61044151,'chr11':46777387,'chr12':59670755}
##Heterotic group A (STN and GON groups)
HAPs_line=["C005.H1","C025.H1","C056.H1","C091.H1","C093.H1","C151.H1","C098.H1","C121.H1","C174.H1","C005.H2","C025.H2","C056.H2","C091.H2","C093.H2","C151.H2","C098.H2","C121.H2","C174.H2","AE.H1"]

def num_of_dSVs_n_recombinations(sample,chr,pos1,pos2):
    df = pd.read_csv(sample+".dSV.bed",header=None,sep='\t')###from dir:/Input_files_from_two_heterotic_groups
    mask = df[(df[0] == chr) & (df[1] >= pos1) & (df[1] <= pos2)]
    return len(mask)

#Split into 50 bins for each chromosome
chr=args.chr
out=open(args.out,'w')
pos_list=[]
for i in range(1,50):
    pos_list.append(chromosomes[chr]*0.02*i)

# Output less dSVs IPH A lines, the output could be used to detected less dSNPs for the next step.
for combo_pos in itertools.combinations(pos_list,4):
    point=[]
    point.append(0)
    for i in combo_pos:
        point.append(i)
    point.append(chromosomes[chr])
    point=sorted(point)
    target=0
    name=''
    for i in range(0,len(point)-1):
        pos1=point[i]
        pos2=point[i+1]
        Part_sample=[]
        Part_name=[]
        for sample in HAPs_line:
            SVs_num=num_of_dSVs_n_recombinations(sample,chr,pos1,pos2)
            Part_sample.append(SVs_num)
        min_SVs=min(Part_sample)
        for i,v in enumerate(Part_sample):
            if v==min_SVs:
                name+=HAPs_line[i]+'_'
        target+=int(min_SVs)
        name=name[:-1]
        name+='\t'
    name=name[:-1]
    ll=''
    for one in combo_pos:
        ll+=str(one)+'\t'
    ll+=name+'\t'
    out.write(str(target)+'\t'+ll[:-1]+'\n')
    out.flush()
