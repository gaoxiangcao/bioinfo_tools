#!/home/caogaoxiang/miniconda3/bin/python
# -*- coding: UTF-8 -*-

import sys

fun_file = sys.argv[1] #分类文件输入
anno_list = sys.argv[2] #注释文件输入
Describe_file = open(fun_file,'r')
Class_file = open(anno_list,'r')
Describe_line = Describe_file.readlines()
Class_line = Class_file.readlines()
newgtf_name = 'DrawAnnotationPic.R'
desktop_path = './'
file_path = desktop_path+newgtf_name+'.txt'
file = open(file_path,'w')
for i in Describe_line[1:]:
    anno_ID = i.split('\t')
    describe = anno_ID[1].replace('\n','')
    count = 0
    for j in Class_line:
        COG_ID = j.split('\t')
        #print(COG_ID)
        if len(COG_ID) > 2:
            #print(COG_ID[-2])
            #print(anno_ID[0])
            if anno_ID[0] == COG_ID[-2]:
                count += 1
            else:
                continue
    COG_list = anno_ID[0]+"\t"+describe+"\t"+str(count)+"\n"
    #print(COG_list)
    file.write(COG_list)
file.close()
Class_file.close()
Describe_file.close()
