perl 01_merge_all_tissue_eQTL_greater0.pl #合并部分数据集位点并利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz"筛选落在eur maf >0.05的eqtl,同时利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_greater0.txt.gz" 补全得"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz"
#并得进入筛选的数据集信息"../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv"
perl 01_1_split_chr1_22.pl #将"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz" split 成chr1-22,得"../output/chr_split/01_merge_all_tissue_cis_eQTL_1kg_Completion_chri.txt.gz"
perl 012_merge_all_tissue_eQTL_tissue_label.pl  #将"../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv" 中的文件merge在一起，并转bed,得../output/01_merge_all_tissue_cis_eQTLtissue_label.bed.gz
#-------------egene
    perl 01_merge_all_tissue_egene.pl #利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv" 部分文件，提取pos，p,gene得"../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz"
    perl 02_filter_sig_egene.pl #
    zless "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz"
    zgrep -v "SNP_chr" "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz"
    zgrep -v "SNP_chr" "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_1e_5.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_1e_5_sorted.bed.gz"
    cd 

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

Rscript 042_circos_density.R

perl 07_annotation_markers.pl
Rscript 071_heatmap_annotation.R
perl 072_count_anno_histone_mark.pl ##对../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}*的marker进行count,得$out_dir/${group}_histone_marker.txt.gz
perl 073_calculate_jaccard_index_mark.pl ##对 $input_dir/${mark}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18snp.bed.gz 计算 jaccard index得$out_dir/${group}_cutoff_${cutoff}_marker_jaccard_index.txt.gz

bedtools intersect -a "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz" -wa -wb |cut -f1-3,7 |sort -u|sort -k1,1 -k2,2n |gzip > "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"

bedtools intersect -a "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz" -wa -wb |cut -f1-3,7 |sort -u|sort -k1,1 -k2,2n |gzip > "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"

bedtools intersect -a "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_1e_5_sorted.bed.gz" -wa -wb |cut -f1-3,7 |sort -u|sort -k1,1 -k2,2n |gzip > "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_1e_5.bed.gz"

bedtools intersect -a "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz" -wa -wb |cut -f1-7 |sort -u|sort -k1,1 -k2,2n |gzip > "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz"

#---带有tissue特征
bedtools intersect -a "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz" -wa -wb |sort -k1,1 -k2,2n |gzip > "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8_all_info.bed.gz"


perl 08_anno_hotspot_marker_signalValue.pl #
perl 09_anno_gene_marker_signalValue.pl # P: 5e-8
#--------------------------------------- 09_egene_pos_5e_8 kmer
bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/09_egene_pos_5e_8.bed.gz" -fo "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/09_egene_pos_5e_8.fa"

cd /home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/

source activate seekr_source
seekr_kmer_counts 09_egene_pos_5e_8.fa  -o  egene_pos_5e_8_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip egene_pos_5e_8_6mers_uc_us_no_log.csv

#--------------------------------------- 09_egene_pos_0.05 kmer

perl merge_egene0.05_and_pos.pl ##"/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 注释ensg位置得"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot_0.05egene_pos.bed.gz",排序得"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos_sorted.bed.gz"

bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos_sorted.bed.gz" -fo "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene_pos_0.05.fa"

cd /home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/

source activate seekr_source
seekr_kmer_counts egene_pos_0.05.fa  -o  egene_pos_0.05_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip egene_pos_0.05_6mers_uc_us_no_log.csv

#
perl anno_gene0.05_marker_signalValue.pl