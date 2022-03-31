perl 001_merge_hotspot_and_egene.pl #利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz"为 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed(../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz)寻找egene,得 "../output/01_hotspot_length_small_than_2290_sorted_egene_h.bed.gz"
python 01_graph_gen_huan.py
python 02_random_walk.py ../output/01_train_graph.dgl 50
python 02_random_walk.py ../output/01_test_graph.dgl 50
python 02_random_walk.py ../output/01_val_graph.dgl 50

hicConvertFormat --matrix ENCFF497EDU.h5 \ --inputFormat h5
-o GInteration_example --outputFormat GInteractions


hicConvertFormat --matrix ENCFF497EDU.h5 \ --inputFormat h5 -o GInteration_example --outputFormat GInteractions

hicConvertFormat --matrix ENCFF497EDU.h5 --inputFormat h5 -o GInteration_example --outputFormat GInteractions