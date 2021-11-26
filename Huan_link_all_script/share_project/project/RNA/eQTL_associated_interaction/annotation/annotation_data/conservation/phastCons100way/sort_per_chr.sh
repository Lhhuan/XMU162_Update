perl 01_normalized_data_per_chr.pl
for i in $( seq 22);
do
	echo $i;
    zless ./normalized_per_chr/chr$i.normalized_merge_phastCons100way.bed.gz | sort -k1,1 -k2,2n | gzip > ./normalized_per_chr/chr$i.sorted_merge_phastCons100way.bed.gz
done