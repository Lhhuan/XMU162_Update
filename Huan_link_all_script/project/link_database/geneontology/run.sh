wget -c http://geneontology.org/gene-associations/goa_human.gaf.gz
perl 01_extract_gene_name.pl #提取goa_human.gaf.gz gene 名字得 "01_unique_gene.csv.gz"

Rscript 02_transfrom_symbol_ensg.R #02_symbol_ENSG.txt

perl 03_adjust_edgepredict_format.pl # transform goa_human.gaf.gz 的symbol to ensg得 /home/huanhuan/project/eQTL_Catalogue/script/prediction/Edgeprediction/output/train/03_go_edgepredict.csv