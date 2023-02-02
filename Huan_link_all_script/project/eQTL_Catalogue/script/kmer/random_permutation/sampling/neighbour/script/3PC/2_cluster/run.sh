Rscript 11_whole_genome_kpca_louvain_leiden_pcn3.R
perl 12_3pcn_whole_genome_homer_200.pl
#=============================eQTL
Rscript 13_count_cluster_specific_eqtl_gene.R
Rscript 13_3_eqtl_count_hotspot_length.R
Rscript 13_4_eqtl_pvalue.R 
#=============================gc content
perl 14_gc_content.pl
perl 15_boxplot_gc_content.R
#=========================
cp "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/whole_genome/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/random_from_whole_genome.bed.gz

gunzip /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/random_from_whole_genome.bed.gz

cp "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/random_from_cold_region.bed.gz
gunzip /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/random_from_cold_region.bed.gz





perl 16_1_plot_marker_peak_point.pl
perl 16_1_plot_marker_peak_point_miss_as_zero.pl
perl 16_1_plot_marker_peak_point_miss_as_zero_4line.pl
cp "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/
gunzip /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/markers/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz


perl 16_1_plot_marker_peak_point_miss_as_zero_3line.pl
perl 17_merge_gwas_and_hotspot.pl
Rscript 18_boxplot_gwas_loci.R
Rscript 18_boxplot_gwas_trait.R

Rscript 19_barplot_gc_island.R

#=================cCREs
bedtools intersect -a ../../../output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz -b "/share/data0/QTLbase/huan/Cis_Regulatory_Elements/encodeCcreCombined_sorted.bed.gz" -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/hotspot_cluster_cCREs.bed.gz

Rscript 22_barplot_cCREs.R 
Rscript 22_2_heatmap_cCREs.R 
perl 22_3_merge_cCREs_enrichment.pl
Rscript 22_4_plot_heatmap_cCREs_enrichment.R 

#=====================phastCons100way

bedtools intersect -a "/share/data0/QTLbase/huan/Conservation/hg38.phastCons100way_transform.bed.gz" -b ../../../output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz  -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_2cluster/Conservation/hotspot_cluster_phastCons100way.bed.gz
Rscript 23_barplot_Conservation_phastCons100way.R 

#========================hic
Rscript 25_boxplot_chr_intra.R
Rscript 26_boxplot_chr_interaction_frequency.R



#===========================chrom state 
# perl 29_split_chrom_state.pl 
# perl 29_1_generate_noGap.pl 
perl 29_2_merge_chrom_state_enrichment.pl
cp "/share/data0/QTLbase/huan/Adaptive_cell_script/20220901_HAQER_Screen_Local/chromHmm/chromHmm_1581.xlsx" /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/ChromHMM/
Rscript 29_3_plot_heatmap.R 




#================HAQER
perl 30_1_merge_HAQER.pl 
Rscript 30_2_HAQER_plot_heatmap.R 