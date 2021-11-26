wget -c http://promoter.bx.psu.edu/hi-c/downloads/loops-hg19.zip
unzip loops-hg19.zip
perl 02_merge.pl #将hg19/下面的.loops文件合在一起，由这些文件判断，该文件是 0-based,判断cis-trans_interaction 并去重得02_all_unique_hic_loops_tran_cis_1.bed.gz(interaction1_interaction2),
#得02_all_unique_hic_loops_tran_cis_2.bed.gz， 得不去重的文件得02_all_hic_loops_source.bed.gz

zless 02_all_unique_hic_loops_tran_cis_1.bed.gz |awk 'NR>1'|sort -u |sort -k1,1 -k2,2n > 02_all_unique_hic_loops_tran_cis_1_sorted.bed
zless 02_all_unique_hic_loops_tran_cis_2.bed.gz |awk 'NR>1'|sort -u |sort -k1,1 -k2,2n > 02_all_unique_hic_loops_tran_cis_2_sorted.bed
cat 02_all_unique_hic_loops_tran_cis_1_sorted.bed 02_all_unique_hic_loops_tran_cis_2_sorted.bed | sort -u |sort -k1,1 -k2,2n >02_all_unique_hic_loops_tran_cis_1_2_sorted.bed

gzip 02_all_unique_hic_loops_tran_cis_1_sorted.bed
gzip 02_all_unique_hic_loops_tran_cis_2_sorted.bed
gzip 02_all_unique_hic_loops_tran_cis_1_2_sorted.bed
perl 03_filter_chr1_22.pl #过滤出02_all_unique_hic_loops_tran_cis_1_2_sorted.bed.gz中的chr1_22,得02_all_unique_hic_loops_tran_cis_1_2_sorted_chr1_22.bed.gz
# zless 02_all_hic_loops_tran_cis_1.bed.gz | head -1 > 02_all_unique_hic_loops_tran_cis_1.bed
# zless 02_all_hic_loops_tran_cis_1.bed.gz |awk 'NR>1'|sort -u >> 02_all_unique_hic_loops_tran_cis_1.bed
# gzip 02_all_unique_hic_loops_tran_cis_1.bed

# zless 02_all_hic_loops_tran_cis_2.bed.gz | head -1 > 02_all_unique_hic_loops_tran_cis_2.bed
# zless 02_all_hic_loops_tran_cis_2.bed.gz |awk 'NR>1'|sort -u >> 02_all_unique_hic_loops_tran_cis_2.bed
# gzip 02_all_unique_hic_loops_tran_cis_2.bed