perl 02_filter_normal_cell_line.pl #将"../01_need_cell.txt" 从已有的细胞信息中对应选出来，得02_exists_cell_line_info.txt，没出现的cell得02_NO_info_cell_line.txt
cp 02_NO_info_cell_line.txt 02_NO_info_cell_line_mannual_info.txt #对 02_NO_info_cell_line_mannual_info.txt 手动填写cell line info
cat 02_exists_cell_line_info.txt 02_NO_info_cell_line_mannual_info.txt >02_all_cell_line_info.txt
perl 03_merge_normal_cell_line_profile.pl #将02_all_cell_line_info.txt 中正常细胞的bed合并，得03_normal_cell_line_profile.bed.gz

liftOver 03_normal_cell_line_profile.bed.gz "/home/huanhuan/reference/hg19ToHg38.over.chain.gz"  ./hg38/03_normal_cell_line_profile.bed  03_normal_cell_line_profile_unmap.bed

gzip 03_normal_cell_line_profile.bed 
cd hg38 
less 03_normal_cell_line_profile.bed |cut -f1-3|sort -u |sort -k1,1 -k2,2n |gzip >03_normal_cell_line_profile_sort.bed.gz
less 03_normal_cell_line_profile.bed |sort -u |sort -k1,1 -k2,2n |gzip >03_normal_cell_line_profile_info_sort.bed.gz
