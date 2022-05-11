perl 001_adjust_hotspot.pl #调整"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"的格式得"../output/01_hotspot_egene_h.bed.gz"

python 01_graph_gen_huan.py
python 02_random_walk.py ../output/01_train_graph.dgl 50
python 02_random_walk.py ../output/01_test_graph.dgl 50
python 02_random_walk.py ../output/01_val_graph.dgl 50