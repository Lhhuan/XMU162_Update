Rscript 01_id_for_hotspot_and_gene.R #利用"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" 为chr1 hotspot和gene 建立index分别得"../output/01_hotspot_egene_idx.txt.gz"
perl 02_add_index_for_interaction.pl ##用"../output/01_hotspot_idx.txt.gz" 和"../output/01_gene_idx.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"添加idx,得"../output/02_hotspot_egene_idx.txt.gz"
Rscript 03_max_marker_signalValue_for_egene_hotspot.R 
Rscript 04_merge_egene_hotspot_annotation_and_idx.R 


python 03_data_processing_huan.py
python 04_random_walk.py ../output/val_graph.dgl 60
python 04_random_walk.py ../output/test_graph.dgl 60




#--------------------------




perl 07_get_train_val_test_for_GATNE.pl #
Rscript 07_get_train_val_test_for_GATNE.R #