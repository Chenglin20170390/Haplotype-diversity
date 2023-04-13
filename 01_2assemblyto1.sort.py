#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   01_2assemblyto1.py
@Time    :   2021/12/22 12:42:26
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   phasing process with 
'''

import datetime, sys
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile_1=sys.argv[1]  ##assembly 1
infile_2=sys.argv[2]   ##assembly 2
outfile=sys.argv[3]  ##merge.assembly

f_in_1=open(infile_1,'r')
f_in_2=open(infile_2,'r')
f_out=open(outfile,'w')

dic_1={}
list_1=[]
num_1=1
for line in f_in_1:
    if line.startswith('>'):
        in_line=line.strip().split()    
        f_out.write(' '.join(in_line)+'\n')
    else :
        in_line=line.strip().split()
        if num_1 <= 12:
              
            dic_1[str(num_1)] = in_line
            num_1 += 1
        else:
            list_1.append(in_line[0])
f1_num=int(in_line[0])
print('Contig number in file 1 :  ' + str(f1_num))

dic_2={}
list_2=[]
num_2=1
for line in f_in_2:
    if line.startswith('>'):
        in_line=line.strip().split()  
        in_line[1]=str(int(in_line[1])+f1_num)  
        f_out.write(' '.join(in_line)+'\n')
    else :
        in_line=line.strip().split()
        tmp_list=[]
        for num in range(len(in_line)):
            ind_num=int(in_line[num])+f1_num
            if ind_num < f1_num:  ##判断为负值 方向为反方向
                ind_num=-f1_num+int(in_line[num])
            tmp_list.append(str(ind_num))
        if num_2 <= 12:
            dic_2[str(num_2)] = tmp_list
            num_2 += 1
        else:
            list_2.append(str(ind_num))
f2_num=int(in_line[0])
print('Contig number in file 2 :  ' + str(f2_num))


for num in range(1,13):
    f_out.write(' '.join(dic_1[str(num)])+'\n')
    f_out.write(' '.join(dic_2[str(num)])+'\n')
f_out.write('\n'.join(list_1)+'\n')
f_out.write('\n'.join(list_2))


f_in_1.close()
f_in_2.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')