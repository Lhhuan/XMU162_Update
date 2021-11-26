#----------------------- #ALL eQTL, 全部的eQTL，不分种族
Rscript huan_NHPoisson_all_eQTLbase1.R #cutoff =7,
Rscript huan_NHPoisson_all_eQTLbase2.R #cutoff =7,cutoff = 24.25(5e-8)
Rscript huan_NHPoisson_all_eQTLbase_p.R  #huan_NHPoisson_all_eQTLbase.R 的多线程
Perl count_number_by_emplambda_in_different_interval_cutoff7_all_eQTL.pl ##统计../output/NHPoisson_emplambda_interval_${i}_all_eQTLbase.txt中 中每个阶段eqtl的数目，
#得../output/count_number_by_emplambda_in_different_interval_cutoff_7_all_eQTL.txt
perl count_number_by_emplambda_in_different_interval_cutoff_all_eQTL.pl ###统计../../output/ALL_${xQTL}/NHPoisson_emplambda_interval_${i}cutoff_${cutoff}_all_${xQTL}.txt中 中每个阶段eqtl的数目，
#得../../output/ALL_${xQTL}/count_number_by_emplambda_in_different_interval_${cutoff}_all_${xQTL}.txt
Rscript density_NHPoisson_batch_all_eQTLbase.R
Rscript bar_plot_NHPoisson_all_eQTL.R

#--------------------------new data
Rscript huan_NHPoisson_all_eQTLbase2.R
perl count_number_by_emplambda_in_different_interval_cutoff_all_eQTL.pl
Rscript bar_plot_NHPoisson_all_eQTL2.R
Rscript bar_plot_NHPoisson_all_eQTL2_log.R

Rscript point_Plot_emplambda_distance.R