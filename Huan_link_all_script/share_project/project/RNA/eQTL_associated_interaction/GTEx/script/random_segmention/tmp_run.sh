# perl 04_split_chr_and_select_random_chr1_22.pl
# echo -e "04_split_chr_and_select_random_chr1_22\n"
cd non_factor_split
perl 07_annotation_non_factor_split_interval_18_different_cutoff.pl
echo -e "07_annotation_non_factor_split_interval_18_different_cutoff\n"
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl
echo -e "08_count_factor_and_non_factor_annotation_interval18_different_cutoff\n"
# Rscript 10_fisher_exact_test_factor_and_non_factor.R
# echo -e "10_fisher_exact_test_factor_and_non_factor.R\n"
# Rscript 10_plot_point_factor_non_factor.R
# echo -e "10_plot_point_factor_non_factor\n"
# perl 081_sum_random_result.pl
# echo -e "081_sum_random_result\n"
# Rscript 10_fisher_exact_test_factor_and_non_factor_summary.R
# echo -e "10_fisher_exact_test_factor_and_non_factor_summary\n"
# Rscript 10_plot_point_factor_non_factor_sum.R
# echo -e "10_plot_point_factor_non_factor_sum\n"

# cd ..

# Rscript 04_histgram_density_interval.R
# echo -e "04_histgram_density_interval\n"
# Rscript 061_histgram_plot_the_length_of_hotspot.R
# echo -e "061_histgram_plot_the_length_of_hotspot\n"