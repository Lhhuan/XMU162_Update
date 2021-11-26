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
qsub eQTL.job
#qdel nubmer #是kill这个进程（比如 qdel 113514）
perl 04_normalized_bed.pl #将"../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz
#得"../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL_normalied.bed.gz"

zless ../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.bed.gz |awk 'NR>1' | sort -k1,1 -k2,2n > ../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed
gzip ../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed


qsub eQTL_supp.job

#-------------
Rscript filter_chr1_QTLbase.R
perl filter_chr1_QTLbase.pl
perl filter_all_QTLbase.pl 
perl 02_filter_emplambda_in_and_out_QTLbase.pl

Rscript 031_histgram_density_interval.R
Rscript 032_histgram_density_interval_chr_specific.R


Rscript qqplot.R

Rscript ks_test.R

Rscript qqplot_specific.R

Rscript hotspot_overlap.R
perl hotspot_overlap.pl

#--------------------------------

perl filter_hotspot_eQTL.pl ###QTLs =eQTL在@interval = (6,7,8,9,12,15)时的hotspot(segment),得"../output/ALL_${QTL}/hotspot/interval_${j}_segment_hotspot.txt.gz"; 得point 热点../../output/ALL_${QTL}/hotspot/interval_${j}_point_hotspot.txt.gz
perl filter_hotspot_eQTL_18.pl ###@QTLs =("hQTL","mQTL","sQTL")在interval = 18时的hotspot(segment),得"../output/ALL_${QTL}/hotspot/interval_${j}_segment_hotspot.txt.gz", 得point 热点../../output/ALL_${QTL}/hotspot/interval_${j}_point_hotspot.txt.gz
bash bedtools_intersect_hotspot.sh #寻找interval 为6,7,8,9,12,15,18时共同的hotspot,得../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz

#---------------------
perl filter_NON_hotspot_eQTL.pl ###QTLs =eQTL在@interval = (6,7,8,9,12,15)时的NON hotspot(segment),得"../output/ALL_${QTL}/un_hotspot/interval_${j}_segment_UNhotspot.txt.gz"; 得point 热点../../output/ALL_${QTL}/un_hotspot/interval_${j}_point_UNhotspot.txt.gz
perl filter_NON_hotspot_eQTL_18.pl ###QTLs =eQTL在interval = 18时的NON hotspot(segment),得"../output/ALL_${QTL}/un_hotspot/interval_${j}_segment_UNhotspot.txt.gz", 得point 热点../../output/ALL_${QTL}/un_hotspot/interval_${j}_point_UNhotspot.txt.gz

bash bedtools_intersect_non_hotspot.sh #寻找interval 为6,7,8,9,12,15,18时共同的 non hotspot,得../../output/ALL_eQTL/un_hotspot/segment_overlap/15_6_7_8_9_12_18_segment_un_hotspot.bed.gz