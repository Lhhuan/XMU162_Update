perl 01_id_for_hotspot_and_gene.pl #利用"../../../../../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz" 为chr1 hotspot和gene 建立index分别得"../output/01_hotspot_idx.txt.gz","../output/01_gene_idx.txt.gz"
perl 02_add_index_for_interaction.pl #利用"../output/01_hotspot_idx.txt.gz" 和"../output/01_gene_idx.txt.gz" 为"../../../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz"和"../../../output/nodes_annotation/hotspot_annotation.txt.gz" 增加index,既用index 代替原来hotspot或者gene,得 ../output/02_hotspot_egene_idx.txt.gz, ../output/02_egene_interaction_idx.txt.gz, ../output/02_egene_co-expression_idx.txt.gz, ../output/02_hotspot_egene_reactomeFI_co-expression_idx.bed.gz, ../output/02_hotspot_anno_idx.txt.gz
python 03_data_processing_huan.py
python 

