#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author   : WN
# @Time     : 20230821
# @File     : Calculate_hybrids_from_ideal_inbred_lines.py
# @Project  : scripts/Designing Ideal potato haplotypes (IPH)

import numpy as np
import pandas as pd
chromosomes={'chr01':88591686,'chr02':46102915,'chr03':60707570,'chr04':69236331,'chr05':55599697,'chr06':59091578,'chr07':57639317,'chr08':59226000,'chr09':67600300,'chr10':61044151,'chr11':46777387,'chr12':59670755}
chromosomes_list=[88591686,46102915,60707570,69236331,55599697,59091578,57639317,59226000,67600300,61044151,46777387,59670755]

def detail_of_dSNPs_n_recombinations(sample,chr,pos1,pos2):
    df = pd.read_csv(sample+".dSNP.bed",header=None,sep='\t')
    mask = df[(df[0] == chr) & (df[1] >= pos1) & (df[1] <= pos2)]
    mask.to_csv('mask.csv')
    list=mask[2].tolist()
    list=sorted(list)
    return list

#Analysis of 12 chromosomes
for select_line in range(0,12):
    #A-lines indicate the IPH E with 12 chromosomes
    list_all_A=[]
    with open("A-lines.txt","r") as fd:
        line = fd.readlines()[select_line]
        line=line.replace("\n","").split('\t')
        chr=line[0]
        sample_list=line[7:12]
        point=[]
        point.append(0)
        for one in line[3:7]:
            one=float(one)
            point.append(one)
        point.append(chromosomes[chr])
        for ie in range(5):
            pos1=point[ie]
            pos2=point[ie+1]
            sample=sample_list[ie]
            out_list=detail_of_dSNPs_n_recombinations(sample,chr,pos1,pos2)
            for qw in out_list:
                list_all_A.append(qw)
        #print(len(list_all_A))
        #screen output


    #A-lines indicate the IPH E with 12 chromosomes
    list_all_E=[]
    with open("E-lines.txt","r") as fd:
        line = fd.readlines()[select_line]
        line=line.replace("\n","").split('\t')
        chr=line[0]
        sample_list=line[7:12]
        point=[]
        point.append(0)
        for one in line[3:7]:
            one=float(one)
            point.append(one)
        point.append(chromosomes[chr])
        for ie in range(5):
            pos1=point[ie]
            pos2=point[ie+1]
            sample=sample_list[ie]
            out_list=detail_of_dSNPs_n_recombinations(sample,chr,pos1,pos2)
            for qw in out_list:
                list_all_E.append(qw)
        set1 = set(list_all_A)
        set2 = set(list_all_E)
        common = set1 & set2
        #print(len(common))
        #print(len(list_all_E))
        #screen output

