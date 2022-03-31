perl 01_merge_all_tissue_eQTL.pl #合并部分数据集位点并利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz"筛选落在eur maf >0.05的eqtl,同时利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.01.txt.gz" 补全得"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz"
#并得进入筛选的数据集信息"../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv"
perl 01_1_split_chr1_22.pl #将"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz" split 成chr1-22,得"../output/chr_split/01_merge_all_tissue_cis_eQTL_1kg_Completion_chri.txt.gz"
#-------------egene
    perl 01_merge_all_tissue_egene.pl #利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv" 部分文件，提取pos，p,gene得"../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz"
    perl 02_filter_sig_egene.pl #
    zless "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz"
    zgrep -v "SNP_chr" "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz"


#-----------------
Rscript 02_NHP_big_par.R #
Rscript 02_NHP_big_par_1_11.R
Rscript 02_NHP_big_par_12_22.R

perl 02_1_merge_nhp_chr1_22.pl #将"../output/all_tissue_status/NHP/CHR_${i}_NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz" 合并得"../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz"

perl 03_filter_hotspot_for_interval18.pl ###"../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${tissue}.txt.gz" 得hotspot(segment),"../output/all_tissue_status/hotspot/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz

perl 04_extend_hotspot_18snp.pl ##用"../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz"对"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz"扩展，只对<18 个snp的hotspot进行扩展，以hotspot中心开始，扩到18个snp,得"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18snp.bed.gz",对其进行排序得"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted.bed.gz"
bedtools merge -i "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted.bed.gz" |gzip >"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"

Rscript 041_plot_distribution_the_length_of_hotspot_ori.R
Rscript 041_barplot_hotspot_length_distribution_ori.R

Rscript 041_plot_distribution_the_length_of_hotspot_extend_18snp.R
Rscript 041_barplot_hotspot_length_distribution_extend_18snp.R



perl 07_annotation_markers.pl

perl 072_count_anno_histone_mark.pl ##对../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}*的marker进行count,得$out_dir/${group}_histone_marker.txt.gz
perl 073_calculate_jaccard_index_mark.pl ##对 $input_dir/${mark}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18snp.bed.gz 计算 jaccard index得$out_dir/${group}_cutoff_${cutoff}_marker_jaccard_index.txt.gz




perl count_eqtl_and_variant.pl #