#Unlike many conservation-scoring programs, phastCons does not rely on a sliding window of fixed size; therefore, short highly-conserved regions and long moderately conserved regions can both obtain high scores.
#score high, highly-conserved
# bash download.sh
perl 01_normalized_data.pl #将chr1-22 merge在一起， #将chr${i}.phastCons100way.wigFix.gz normalized成正常的文件，得merge_phastCons100way.bed.gz
zless merge_phastCons100way.bed.gz | sort -k1,1 -k2,2n | gzip > merge_phastCons100way_sorted.bed.gz
rm merge_phastCons100way.bed.gz
perl 01_normalized_data_per_chr.pl #将chr${i}.phastCons100way.wigFix.gz normalized成正常的文件，得./normalized_per_chr/chr${i}.normalized_merge_phastCons100way.bed.gz
bash sort_per_chr.sh #对/normalized_per_chr/chr$i.normalized_merge_phastCons100way.bed.gz 进行排序，得./normalized_per_chr/chr$i.sorted_merge_phastCons100way.bed.gz
bash rm_replicate.sh #rm /normalized_per_chr/chr$i.normalized_merge_phastCons100way.bed.gz


# sortBed -i merge_phastCons100way.bed.gz | gzip > merge_phastCons100way_sorted.bed.gz

# less merge_phastCons100way.bed | sort -k1,1 -k2,2n > merge_phastCons100way_sorted.bed