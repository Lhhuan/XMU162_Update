# perl 01_Completion_snp_for_cis_trans_eQTL_by_1kg.pl | sort -u
# echo -e "01_Completion_snp_for_cis_trans_eQTL_by_1kg\n"
# Rscript 02_huan_NHPoisson.R
# echo -e "huan_NHPoisson\n"
# perl 03_filter_emplambda_in_eQTL.pl
# echo -e "03_filter_emplambda_in_eQTL\n"

perl 06_overlap_count_factor_in_and_out_hotspot.pl 
echo -e "06_overlap_count_factor_in_and_out_hotspotn\n"
perl 07_overlap_prepare_for_fisher.pl 
echo -e "07_overlap_prepare_for_fisher\n"
Rscript 08_overlap_fisher_exact_test.R 
echo -e "08_overlap_fisher_exact_test\n"
Rscript 09_overlap_boxplot_for_factor_score.R
echo -e "09_overlap_boxplot_for_factor_score\n"

# perl 06_interval_18_count_factor_in_and_out_hotspot.pl 
# echo -e "06_interval_18_count_factor_in_and_out_hotspotn\n"
# perl 07_interval_18_prepare_for_fisher.pl 
# echo -e "07_interval_18_prepare_for_fisher\n"
# Rscript 08_interval_18_fisher_exact_test.R 
# echo -e "08_interval_18_fisher_exact_test\n"
# Rscript 09_interval_18_boxplot_for_factor_score.R
# echo -e "09_interval_18_boxplot_for_factor_score\n"

