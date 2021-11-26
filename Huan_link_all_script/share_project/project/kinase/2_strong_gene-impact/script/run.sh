#On average, six of these filtered genes per cell line scored as essential below a threshold of −0.6 after CERES correction(pmid:29083409)
# so, the essential gene that CERES <-0.6
Rscript 01_count_col.R #统计../data/DepMap/Achilles_gene_effect.csv的列数,共18334列
perl 02_arrange_gene_effect_to_normal_type_and_filter.pl # 把../data/DepMap/Achilles_gene_effect.csv转换成常用的格式，即，gene cell_line probability
#得../output/02_arrange_gene_effect_to_normal_type.txt, 
#过滤得dependency probability>0的文件得../output/02_arrange_gene_effect_to_normal_type_dependency.txt
Rscript 03_transfrom_kinase_gene_to_entrez_ENSG.R #把../../kinase_ppi/data/kinase.txt 转成entrez和ENSG
#得../output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt
perl 04_filter_kinase_from_Achilles.pl #用../output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt中的kinase提取
#../output/02_arrange_gene_effect_to_normal_type.txt中的kinase
#得kinase文件../output/04_kinase_Achilles.txt
#得dependent_probability 的kinase文件 ../output/04_kinase_Achilles_dependent.txt
#得非 dependent_probability 的kinase文件 ../output/04_kinase_Achilles_non_dependent.txt
#得非kinase 文件 ../output/04_non_kinase_Achilles.txt
perl 05_merge_essential_gene_cell_line_cancer.pl #用../data/DepMap/Cell_lines_annotations_20181226.txt 和../data/DepMap/sample_info.csv为../output/04_kinase_Achilles_dependent.txt
#寻找cell line的信息。首先用../data/DepMap/Cell_lines_annotations_20181226.txt寻找，然后用../output/04_kinase_Achilles_dependent.txt寻找
#得../output/05_merge_essential_gene_cell_line_cancer.txt
perl 06_filter_tcga_code_essential_gene.pl #在../output/05_merge_essential_gene_cell_line_cancer.txt中筛选出有TCGA_code的essential gene
#得../output/06_tcga_code_essential_gene_annotation.txt 和 非tcga code 的文件../output/06_NON-tcga_code_essential_gene_annotation.txt
#得每个tcga code 下的essential_gene文件，../output/06_tcga_code_essential_gene.txt
