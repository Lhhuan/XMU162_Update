
# coding: utf-8
out_file = open('./output/test.txt','w') #write

for lines in open ('../kinase.txt','r'): #文件迭代器, 效果下面三行的综合及格过是一样的
#------------------open file
# f = open('../kinase.txt','r')#read
# lines = f.readline()#read file 
# for line in f.readlines():
    line=line.strip('\n') # skip \n 跳过换行符
    line1= line.split('\t') #split the line by '\t'
    gene= line1[0]
    value = line1[2]
    result = gene+'\t'+value+'\n'
    out_file.write(result)

