#download hg19 reference geneome data from ensembl BioMart 
# Hg19_reference_geneome_data_from_ensembl_BioMart.txt
Rscript 01_transfrom_ensg_to_entrez.R #将../output/Hg19_reference_geneome_data_from_ensembl_BioMart.txt 中的ensg转为entrezgene,得../output/01_transfrom_ensg_to_entrez.txt
perl 02_merge_entrez_drug_info.pl #将../output/Hg19_reference_geneome_data_from_ensembl_BioMart.txt 和../output/01_transfrom_ensg_to_entrez.txt merge 到一起，
#得../output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt