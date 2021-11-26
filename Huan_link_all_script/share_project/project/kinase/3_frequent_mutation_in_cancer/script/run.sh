perl 01_filter_01_kinase_gene.pl #从/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/mc3.v0.2.8.PUBLIC.maf.gz 中筛选出
#../../2_strong_gene-impact/output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt基因相关信息，得../output/01_TCGA_gene_in_01_kinase.txt.gz
perl 02_pancancer_count_mutation_occur_in_gene.pl #pancancer 统计../output/01_TCGA_gene_in_01_kinase.txt.gz中 每个gene的突变情况，
#得../output/02_pancancer_count_mutation_occur_in_gene.txt
#并同时根据 /share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/TCGA-Clinical-Data-Resource.csv
#patients对应的具体cancer_type信息，得../output/02_TCGA_gene_in_01_kinase_cancer.txt
perl 03_cancer_specific_count_mutation_occur_in_gene.pl #cancer specific 统计../output/02_TCGA_gene_in_01_kinase_cancer.txt中 每个gene的突变情况，
#得../output/03_cancer_specific_count_mutation_occur_in_gene.txt 
Rscript 04_plot_pancancer_mutation_occur_in_gene.R #为../output/02_pancancer_count_mutation_occur_in_gene.txt 画图
#得../figure/04_pancancer_mutation_occur_in_gene.pdf
Rscript 05_plot_cancer_specific_mutation_occur_in_gene.R # cancer specific 为../output/03_cancer_specific_count_mutation_occur_in_gene.txt 
#画图，得../figure/05_*_mutation_occur_in_gene.pdf *代表cancer的名字