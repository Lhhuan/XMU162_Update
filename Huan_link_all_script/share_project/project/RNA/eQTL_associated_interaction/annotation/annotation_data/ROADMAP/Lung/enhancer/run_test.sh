cd split_complement_chr1_22_files
echo -e "cd split_complement_chr1_22_files\n"
zless ../split_enhancer_complement_chr1_22.bed.gz | split -l 1000000 -a 5
echo -e " finish zless ../split_enhancer_complement_chr1_22.bed.gz | split -l 1000000 -5\n"
cd /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/non_factor_split/
echo -e "cd /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/non_factor_split/\n"
perl 07_annotation_non_factor_interval_18_different_cutoff_non_factor_split_chr1_22_enhancer.pl
echo -e "finish all\n"