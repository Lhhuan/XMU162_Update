perl 001_merge_hotspot_and_egene.pl #利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz"为 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed(../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz)寻找egene,得 "../output/01_hotspot_length_small_than_2290_sorted_egene_h.bed.gz"
python 01_graph_gen_huan.py

python 02_random_walk.py ../output/01_train_graph.dgl 50
python 02_random_walk.py ../output/01_test_graph.dgl 50
python 02_random_walk.py ../output/01_val_graph.dgl 50

source activate huan_py3
python 03_skipgram3_gcn_huan.py
Rscript 04_transform_and_filter_predict.R
perl 05_filter_ori_true.pl #利用 "/home/huanhuan/reference/grch37_ensgID_position_from_ensembl104.txt" 和"../output/01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt"为"../output/04_transform_and_filter_predict_top0.1.txt" 提取hotspot 和 egene 的位置，得 "../output/05_all_predict_top0.1.bed.gz"，"../output/05_trans_predict_top0.1.bed.gz"，"../output/05_cis_predict_top0.1.bed.gz"，并得"../output/01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt" egene 位置得"../output/05_groundtruth_link_pos.bed.gz"
Rscript 06_transform_h5_to_segment.R
zless ../output/06_Hi_C_result_top25.bed.gz |sort -k1,1 -k2,2n |gzip> ../output/06_Hi_C_result_top25_sorted.bed.gz
 
perl 07_merge_hic_left_prediction.pl #"../output/07_all_predict_top0.1_need_left_overlap_hi_c.bed.gz" header 为$h_chr\t$h_start\t$h_end\t$egene\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end\thi-c_value
perl 07_merge_hic_right_prediction.pl

