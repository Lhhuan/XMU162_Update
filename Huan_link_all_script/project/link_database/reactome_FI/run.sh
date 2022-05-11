wget -c https://reactome.org/download/tools/ReatomeFIs/FIsInGene_122921_with_annotations.txt.zip #Version 2021
unzip FIsInGene_122921_with_annotations.txt.zip
Rscript 01_trans_symbol_to_ENSG.R #./output/01_symbol_ENSG.txt
perl 02_unique_symbol_ENSG.pl #unique ./output/01_symbol_ENSG.txt,得"./output/02_query_gene_entrezgene_symbol_ensembl.txt" 和"./output/query_gene_ensembl.txt"
Rscript 03_merge_reactomeFI_ensg.R #./output/03_FIsInGene_122921_with_annotations_ENSG.txt, ./output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt
perl 04_split_network_by_direaction.pl #将 ./output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt 按照方向分为三部分，start-end: ./output/04_ReactomeFI_start_end.txt end-start: ./output/04_ReactomeFI_end_start.txt 双向: ./output/04_ReactomeFI_both.txt