#----------------------- #ALL caQTL, 全部的caQTL，不分种族
Rscript huan_NHPoisson_all_lncRNAQTL.R #cutoff =7,cutoff = 24.25(5e-8)
Perl count_number_by_emplambda_in_different_interval_cutoff7_all_caQTL.pl ##统计../output/NHPoisson_emplambda_interval_${i}_all_caQTLbase.txt中 中每个阶段caQTL的数目，
#得../output/count_number_by_emplambda_in_different_interval_cutoff_7_all_caQTL.txt
Rscript density_NHPoisson_batch_all_caQTLbase.R
Rscript bar_plot_NHPoisson_all_caQTL.R