perl 01_filter_significant_eQTL_in_tissue.pl
echo -e "01_filter_significant_eQTL_in_tissue\n"
perl 02_annotation_non_factor_split_interval_18_different_cutoff.pl
echo -e "02_annotation_non_factor_split_interval_18_different_cutoff\n"
perl 03_count_eQTL_hit_factor.pl 
echo -e "03_count_eQTL_hit_factor\n"
# Rscript 10_plot_point_factor_non_factor.R
# echo -e "10_plot_point_factor_non_factor\n"
# Rscript 10_fisher_exact_test_factor_and_non_factor.R
# echo -e "10_fisher_exact_test_factor_and_non_factor\n"