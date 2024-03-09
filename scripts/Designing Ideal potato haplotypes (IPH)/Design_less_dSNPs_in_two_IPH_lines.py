#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author   : WN
# @Time     : 20230821
# @File     : Design_less_dSNPs_in_two_IPH_lines.py
# @Project  : scripts/Designing Ideal potato haplotypes (IPH)

import itertools
import argparse
import gc
import pandas as pd
import numpy as np

##for example:8.Five-recmobinations-SVs-stat.E454-Elines.py -chr chr01 -file chr01.${file}_from_Design_less_dSVs_in_A_line_IPH.py -out chr01.results
parser = argparse.ArgumentParser(description = 'Find ideal haplotype for E lines', add_help = False, usage = '\npython3 -chr chr01 -file [chr01.FILE_from_Design_less_dSVs_in_A_line_IPH.py] -out [output.file]')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-chr', '--chr', metavar = '[chr01', help = 'input', required = True)
required.add_argument('-file', '--file', metavar = '[chr01.FILE_from_Design_less_dSVs_in_A_line_IPH.py', help = 'input', required = True)
required.add_argument('-out', '--out', metavar = '[output.txt]', help = 'output', required = True)
optional.add_argument('-h', '--help', action = 'help', help = 'help')
args = parser.parse_args()

chr=args.chr
out=open(args.out,'w')

#DM chromosomes size
chromosomes={'chr01':88591686,'chr02':46102915,'chr03':60707570,'chr04':69236331,'chr05':55599697,'chr06':59091578,'chr07':57639317,'chr08':59226000,'chr09':67600300,'chr10':61044151,'chr11':46777387,'chr12':59670755}

#function for find less dSNPs in IPH
def num_of_dSNPs_n_recombinations(sample,chr,pos1,pos2):
    df = pd.read_csv(sample+".dSNP.bed",header=None,sep='\t')
    mask = df[(df[0] == chr) & (df[1] >= pos1) & (df[1] <= pos2)]
    return len(mask)

#centremeres in 12 chromosomes,'chr02':'0_1' indicate undetectable centromeres in Chr02
#filter recmobinations in INVs
centremere={}
list_exam={'chr01':'33000000_33300000','chr02':'0_1','chr03':'11900000_12200000','chr04':'26300000_27700000','chr05':'27000000_27500000','chr06':'15500000_16800000','chr07':'16800000_17500000','chr08':'18600000_18700000','chr09':'18300000_20300000','chr10':'24500000_25100000','chr11':'21000000_22600000','chr12':'18400000_19400000'}
def INVs_cetremere(sample,chr,pos1,pos2):
    test_re='NA'
    list_aaa=[]
    list_aaa.append(list_exam[chr])
    with open(sample+".dSV.bed",'r') as fff:
        for line in fff:
            if 'INV' in line:
                line=line.replace('\n','').split('\t')
                if line[0]==chr:
                    target_aaa=line[1]+'_'+line[2]
                    list_aaa.append(target_aaa) 
    for one in list_aaa:
        a=int(one.split('_')[0])
        b=int(one.split('_')[1])
        if pos1 >a and pos1<b:
            test_re='NO'
        if pos2 >a and pos2<b:
            test_re='NO'
    return test_re

SVs_in_file=[]
with open(args.file,'r') as fd:
    for line in fd:
        line=line.replace('\n','').split('\t')
        SVs_in_file.append(int(line[0]))
with open(args.file,'r') as fd:
    for line in fd:
        line=line.replace('\n','').split('\t') 
        if int(line[0])==min(SVs_in_file):
            out_line=line[0]+'\t'
            out_point=''
            out_line_sample=''
            out_SNPs=0
            point=[]
            point.append(0)
            type_target=[]
            for ppp in line[1:5]:
                out_point+=ppp+'\t'
                point.append(float(ppp))
            point.append(chromosomes[chr])
            point=sorted(point)
            for num,part_sample in enumerate(line[5:]):
                part_sample=part_sample.split('_')
                pos1=point[num]
                pos2=point[num+1]
                SNPs_part=[]
                for each in part_sample:
                    each_SNPs=num_of_dSNPs_n_recombinations(each,chr,pos1,pos2)
                    type_target.append(INVs_cetremere(each,chr,pos1,pos2))
                    SNPs_part.append(each_SNPs)
                min_SNPs_sample=[]
                min_SNPs_value=min(SNPs_part)
                for v,loci in enumerate(SNPs_part):
                    if loci==min_SNPs_value:
                        min_SNPs_sample.append(part_sample[v])
                s=''
                for recylcy in min_SNPs_sample:
                    s+=recylcy+'_'
                out_line_sample+=s[:-1]+'\t'
                out_SNPs+=min_SNPs_value
            if not 'NO' in type_target:
                ll_final=args.chr+'\t'+out_line+str(out_SNPs)+'\t'+out_point[:-1]+'\t'+out_line_sample[:-1]
                out.write(ll_final+'\n')
                out.flush()
