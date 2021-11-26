#download from  https://drive.google.com/drive/folders/1NnUmxxSBKh6E4szkhTJ0fSJj3v6bpMgE 
# #data: 2019.11.11
#Release 1.1 (2019-09-08)

#unzip all.txt.zip
gzip all.txt

nohup zless all.txt.gz |awk 'NR>1'| cut -f3| sort -u > QTLbase_gene_list.txt &