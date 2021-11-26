
perl 02_random_sampling.pl
echo -e "02_random_sampling\n"
Rscript 03_NHP_big_par.R
echo -e "03_NHP_big_par\n"
gzip /share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/Lung/NHP/*
gzip /share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/Whole_Blood/NHP/*
perl 06_filter_hotspot_for_interval18.pl
echo -e "06_filter_hotspot_for_interval18\n"
cd non_factor_split
perl 07_annotation_non_factor_split_interval_18_different_cutoff.pl
echo -e "07_annotation_non_factor_split_interval_18_different_cutoff\n"
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl
echo -e "08_count_factor_and_non_factor_annotation_interval18_different_cutoff\n"
# Rscript 10_fisher_exact_test_factor_and_non_factor.R
# echo -e "10_fisher_exact_test_factor_and_non_factor.R\n"
# Rscript 10_plot_point_factor_non_factor.R
# echo -e "10_plot_point_factor_non_factor\n"
perl 081_sum_random_result.pl
echo -e "081_sum_random_result\n"
Rscript 10_fisher_exact_test_factor_and_non_factor_summary.R
echo -e "10_fisher_exact_test_factor_and_non_factor_summary\n"
Rscript 10_plot_point_factor_non_factor_sum.R
echo -e "10_plot_point_factor_non_factor_sum\n"

cd ..

Rscript 04_histgram_density_interval.R
echo -e "04_histgram_density_interval\n"
Rscript 061_histgram_plot_the_length_of_hotspot.R
echo -e "061_histgram_plot_the_length_of_hotspot\n"

