# download the file.txt from encode 
#把file.txt第一行放进get_meta.txt
xargs -L 1 curl -O -J -L < get_meta.txt #得metadata.tsv 
perl 01_filter_need_download_file.pl ##将metadata.tsv中选一个数据最大的出来，对于同一样本有重复的，筛选出文件大的哪个，对于并生成download link,得download.sh,和need_metadata.txt
mkdir data
cd data
bash download.sh
cd ..
perl 02_merge_all_data.pl #将data/下面的.gz文件合在一起，并去重得02_all_unique_merge_RBP_narrow_peak.bed.gz,得具有标记的文件02_all_merge_RBP_narrow_peak.bed.gz
zless 02_all_unique_merge_RBP_narrow_peak.bed.gz | sort -k1,1 -k2,2n |gzip >02_all_unique_merge_RBP_narrow_peak_sorted.bed.gz
perl 03_filter_chr1_22.pl ###将02_all_unique_merge_RBP_narrow_peak_sorted.bed.gz中的chr1-22提取出来得03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz
#将03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz左右各扩500bp得03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22_extend_500bp.bed.gz
#将03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz左右各扩100bp得03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz