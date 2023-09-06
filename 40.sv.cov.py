#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   40.sv.cov.py
@Time    :   2022/12/06 17:11:08
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin_solab@163.com
@License :   (C)Copyright 2022-2023, CAAS ShenZhen
@Desc    :   python  .py   samtools_depth .file  sv_bed  sv_ft.bed
'''

import datetime, sys
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile=sys.argv[1]
infile1=sys.argv[2]
outfile=sys.argv[3]
f_in=open(infile,'r')
f_in1=open(infile1,'r')
f_out=open(outfile,'w')

dic={}

for line in f_in:
    in_line=line.strip().split()
    pos=in_line[0]+"_"+in_line[1]
    if int(in_line[2]) >= 3:
        dic[pos]=""

for line in f_in1:
    in_line=line.strip().split()
    pos1=in_line[0]+"_"+in_line[1]
    pos2=in_line[0]+"_"+in_line[2]
    if pos1 and pos2 in dic:
        f_out.write(line)


f_in.close()
f_in1.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')