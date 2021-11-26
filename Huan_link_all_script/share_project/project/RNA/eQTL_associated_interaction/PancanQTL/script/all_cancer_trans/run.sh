perl 02_Completion_snp_for_PancanQTL_by_1kg.pl
Rscript 03_NHP_big_par.R
perl 04_filter_hotspot_for_interval18.pl
perl 05_count_cancer_all_hotspot.pl # 统计trans qtl的数目 ../../output/cancer_total/trans/05_number_of_cancer_all_hotspot.txt.gz
Rscript 06_Cleveland_dot_plot_number_of_cancer_trans_hotspot.R #
perl 07_split_hotspot.pl ##将hotspot makewindows 1MB  ../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}_makewin_${win}.bed.gz
perl 08_extract_cancer_qtl.pl #对"../../data/trans_eQTLs_all_re.gz"进行整理，得cancer 对应的文件""../../output/${cancer}/Trans_eQTL/hotspot_trans_eQTL/interval_${j}/qtl_gene.bed.gz";,得hotspot_gene文件${out_dir}/hotspot_gene.bed.gz
Rscript 09_circos_heatmap_and_link.R 