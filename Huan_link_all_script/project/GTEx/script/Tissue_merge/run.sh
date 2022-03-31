perl 01_transform_sQTL_varint_ID_hg38_to_hg19.pl # #利用"/share/data0/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" 将/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL/*.v8.sqtl_signifpairs.txt.gz转换为hg19, 得/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL_hg19/*.v8.sqtl_signifpairs.txt.gz
perl 011_merge_all_tissue_sQTL_sgene.pl 
Rscript 01_fiter_hotspot_for_homer_and_bar_plot_distribution.R 
#./figure/01_hotspot_distribution.txt,../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed
less /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed | sort -k1,1 -k2,2n > /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed
#---------homer
    source activate huan_py3
    findMotifsGenome.pl ../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed hg19 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/homer/ -size given

#----------

perl 02_calculate_marker_annotation_fraction.pl #"../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz"

bedtools nuc -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa" -bed /share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz |gzip >/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833_GC.bed.gz

Rscript 03_plot_gc_content_distrbution.R 
Rscript 04_pca_cluster.R  #04_10_markers_euclidean_dist.Rdata
Rscript 04_kpca_cluster.R
Rscript 05_Mclust.R

/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz


#-----------------
Rscript 04_pca_cluster_new.R #04_new_pca_tsne.Rdata,04_new_pca_umap.Rdata
Rscript 05_leiden_new.R
Rscript 05_kmeans_new.R 
Rscript 05_Mclust_new.R
Rscript 05_kkl_new.R 
#-------------------------plot
#mart_export.txt from ensembl hg19 1-based
perl 06_merge_sig_QTL_and_gene_info.pl ##为"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz" 找到gene TSS并转bed得../output/Tissue_merge/Cis_eQTL/06_sig_eQTL_gene_TSS.bed.gz,排序得../output/Tissue_merge/Cis_eQTL/06_sig_eQTL_gene_TSS_sorted.bed.gz
perl 07_get_hotspot_to_plot.pl

zless ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_hotspot_eqtl_egene.bed.gz |awk -v OFS="\t" '{print $4,$5,$6,$7}' |sort -u |sort -k1,1 -k2,2n | gzip > ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL.bed.gz

zless ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_hotspot_eqtl_egene.bed.gz |awk -v OFS="\t" '{print $8,$9,$10,$7}' |sort -u |sort -k1,1 -k2,2n |gzip > ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL_TSS.bed.gz

bedtools intersect -a ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL.bed.gz -b ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_100_win.bed.gz -wa -wb |cut -f1-4,8|gzip >../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_eqtl_in_win100.bed.gz  #eqtl\tensg\twinnum

bedtools intersect -a ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL_TSS.bed.gz -b ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_100_win.bed.gz -wa -wb |cut -f1-4,8|gzip >../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_TSS_in_win100.bed.gz #tss\tensg\twinnum

#------------------------------------------win 1000
bedtools makewindows -b "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region.bed" -w 1000 -i winnum | gzip >../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_1kb_win.bed.gz

bedtools intersect -a ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL.bed.gz -b ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_1kb_win.bed.gz -wa -wb |cut -f1-4,8|gzip >../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_eqtl_in_win1kb.bed.gz

bedtools intersect -a ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_eQTL_TSS.bed.gz -b ../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_1kb_win.bed.gz -wa -wb |cut -f1-4,8|gzip >../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_TSS_in_win1kb.bed.gz

#--------------------------------------win10000
perl 06_merge_sig_QTL_and_gene_info.pl
perl 06_merge_sig_QTL_and_gene_info_p_unlimit.pl

perl 07_get_hotspot_to_plot.pl
perl 07_get_hotspot_to_plot_broad_p_unlimit.pl

Rscript 08_interaction_map.R 
Rscript 08_interaction_map_100bp.R
Rscript 08_interaction_map_100kb.R
Rscript 08_interaction_map_100kb_p_unlimit.R

#----sQTL 
perl 06_merge_sQTL_and_gene_info.pl 

#--------
perl merge_all_eqtl.pl 