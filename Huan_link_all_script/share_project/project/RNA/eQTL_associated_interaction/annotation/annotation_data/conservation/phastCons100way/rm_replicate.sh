for i in $( seq 22);
do
	echo $i;
    rm ./normalized_per_chr/chr$i.normalized_merge_phastCons100way.bed.gz
done