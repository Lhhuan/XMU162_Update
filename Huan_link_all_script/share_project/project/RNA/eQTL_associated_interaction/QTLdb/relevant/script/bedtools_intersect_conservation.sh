for i in $( seq 22);
do
	echo $i;
    bedtools intersect -a ../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed.gz -b ../../annotation_data/conservation/phastCons100way/normalized_per_chr/chr$i.sorted_merge_phastCons100way.bed.gz -wa -wb | gzip > ../output/conservation_per_chr/QTLbase_chr$i.all_eQTL_phastCons100way.bed.gz
done