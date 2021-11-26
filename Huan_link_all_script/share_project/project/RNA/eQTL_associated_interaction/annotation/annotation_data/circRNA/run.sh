wget -c http://circbase.org/download/hsa_hg19_circRNA.txt
# #----------------------------------
# less hsa_hg19_circRNA.txt |awk 'NR>1' |cut -f1,2,3,5,8 |sort -k1,1 -k2,2n |gzip >needed_hsa_hg19_circRNA_sorted.bed.gz

less hsa_hg19_circRNA.txt |awk 'NR>1' |cut -f1,2,3 |sort -u |sort -k1,1 -k2,2n |gzip >needed_hsa_hg19_circRNA_sorted.bed.gz
perl 03_filter_chr1_22.pl #过滤出needed_hsa_hg19_circRNA_sorted.bed.gz中的chr1_22,得03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz
#将03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz左右各扩500bp得03_needed_hsa_hg19_circRNA_sorted_chr1_22_extend_500bp.bed.gz
#将03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz左右各扩100bp得03_needed_hsa_hg19_circRNA_sorted_chr1_22_extend_100bp.bed.gz

# wget -c http://circatlas.biols.ac.cn/: human_bed_v2.0.txt.zip hg38_hg19_v2.0.txt.zip
# unzip human_bed_v2.0.txt.zip 
# gzip human_bed_v2.0.txt
# unzip hg38_hg19_v2.0.txt.zip
# gzip hg38_hg19_v2.0.txt
# rm human_bed_v2.0.txt.zip
# rm hg38_hg19_v2.0.txt.zip