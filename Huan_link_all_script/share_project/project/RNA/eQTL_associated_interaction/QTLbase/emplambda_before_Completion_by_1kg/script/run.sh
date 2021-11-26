#------------------------------------------------------------
perl 01_merge_all_kinds_QTL.pl #将../data/中的各种QTL merge到一起，得../output/01_all_kinds_QTL.txt
perl 02_annotate_gene_type.pl #为../output/01_all_kinds_QTL.txt添加genetype 得../output/02_annotate_gene_type.txt.gz
#用"/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt"， "/home/huanhuan/ref_data/RNA_position/data/hg19_v19_miRNA_position.txt"，
#"/home/huanhuan/ref_data/RNA_position/data/hg19_mRNA_position.txt"和"/home/huanhuan/ref_data/RNA_position/data/gencode.v32lift37.long_noncoding_RNAs.gff3.gz"为../output/01_all_kinds_QTL.txt 注释gene type
perl 03_count_snp_QTL_number_in_each_gene_type_and_QTL_type.pl #count ../output/02_annotate_gene_type.txt.gz 各类QTL的数量得文件../output/03_QTL_type_number.txt  在特定QTL种类下的gene type及数量。
#及所有的QTL对应的QTL的数量。
perl count_snp_trait_pos_distance.pl # 计算 ../output/01_all_kinds_QTL.txt中snp 与trait的位置差，得../output/count_snp_trait_pos_distance.txt
Rscript snp_trait_pos_distance_density.R #

#-------------------------------------------------
perl clump_by_sourceid.pl ##分source进行clump,将不同类型QTL按照sourceid分割成"../output/used_to_clump/$used_to_clump"，然后用对应的人种文件/share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID
#进行clump, 得../output/clump/${used_to_clump}.clump
#--------------------------------------------------------------------------------------------------
perl 01_clump_by_tissue_pop_0.5.pl ##分source进行clump,将不同类型QTL按照tissue 和 population 分割成"../output/tissue_pop_used_to_clump/$used_to_clump"，然后用对应的人种文件/share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID
#进行clump, 得 ../output/tissue_pop_used_to_clump/$used_to_clump.clump
perl 02_merge_tag_snp_0.5.pl #将../output/tissue_pop_used_to_clump/$used_to_clump 和../output/tissue_pop_used_to_clump/${used_to_clump}.clump merge 起来，并且计算每个LD的长度得../output/clump_region/$used_to_clump
#得总文件../output/all_qtl_clump_locus_r_square0.5.txt.gz
#----------------------------------------------QTL hotspot

perl filter_chr22.pl ##将../output/01_all_kinds_QTL.txt 中chr22过滤得../output/chr22.txt,过滤chr22中eur得../output/chr22_eur.txt
#得所有的具有 pop 和tissue的文件得../output/merge_QTL_all_QTLtype_pop.txt.gz
perl judge_xQTL_cis_trans.pl ##将../output/merge_QTL_all_QTLtype_pop.txt.gz 按照1MB和10MB划分cis和trans， 分别得../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz


Rscript Manhattan_emplambda_QTL_plot_pa.R #画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure下面
Rscript Manhattan_emplambda_QTL_plot_chr_pa.R #分染色体画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure/Manhattan/下面
Rscript Manhattan_qqman_emplambda_QTL_plot_chr_pa.R  #用qqman分染色体画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure/Manhattan/下面
#--------------------------
perl merge_emplambda_and_org_xQTL.pl #将../output/xQTL/count_number_by_emplambda_in_different_interval_7.3_all_xQTL.txt 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
perl merge_emplambda_and_org_all_QTLbase.pl #对整个QTLbase进行操作 #xQTL 是全部的QTLbase
#将../output/xQTL/count_number_by_emplambda_in_different_interval_7.3_all_xQTL.txt 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
#------------------------------------
perl overlap_emplambda_xQTL.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt按照相同染色体位置，两两取交集，
#得../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz
perl overlap_emplambda_xQTL_cis.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz按照相同染色体位置，两两取交集，
#得 cis ../output/xQTL_merge/Cis_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_1MB_${QTL1}_${QTL2}.txt.gz
#和../output/xQTL_merge/Cis_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz
perl overlap_emplambda_xQTL_cis_trans.pl ##("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz按照相同染色体位置，两两取交集，
#得 QTL1为cis, QTL2为trans文件 ../output/xQTL_merge/Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_cis_${QTL1}_trans_${QTL2}.txt.gz
#和../output/xQTL_merge/Trans_10MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_cis_${QTL1}_trans_${QTL2}.txt.gz
perl overlap_emplambda_xQTL_trans_cis.pl ##("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz按照相同染色体位置，两两取交集，
#得 QTL1为trans, QTL2为cis文件 ../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_trans_${QTL1}_cis_${QTL2}.txt.gz
#和../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_trans_${QTL1}_cis_${QTL2}.txt.gz

perl overlap_emplambda_xQTL_trans.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz按照相同染色体位置，两两取交集，
#得 trans ../output/xQTL_merge/Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_trans_1MB_${QTL1}_${QTL2}.txt.gz
#和../output/xQTL_merge/Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_trans_10MB_${QTL1}_${QTL2}.txt.gz
Rscript Point_plot_emplambda_relevance_pdf.R 
Rscript Point_plot_emplambda_relevance_png.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript Point_plot_emplambda_relevance_per_chr_pdf.R #分染色体画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript Point_plot_emplambda_relevance_per_chr_png.R
Rscript Point_plot_emplambda_relevance_cis_pdf.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript Point_plot_emplambda_relevance_cis_png.R
Rscript Point_plot_emplambda_relevance_cis_trans_pdf.R #../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_trans_${QTL1}_cis_${QTL2}.txt.gz
#和../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_trans_${QTL1}_cis_${QTL2}.txt.gz
# ../output/xQTL_merge/Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_cis_${QTL1}_trans_${QTL2}.txt.gz
#和../output/xQTL_merge/Trans_10MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_cis_${QTL1}_trans_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript Point_plot_emplambda_relevance_cis_trans_lm_pdf.R
Rscript Point_plot_emplambda_relevance_cis_trans_lm_png.R
Rscript Point_plot_emplambda_relevance_cis_trans_glm_pdf.R
Rscript Point_plot_emplambda_relevance_cis_trans_glm_png.R
Rscript Point_plot_emplambda_relevance_cis_trans_loess_png.R
Rscript Point_plot_emplambda_relevance_trans_pdf_fitting.R
Rscript Point_plot_emplambda_relevance_cis_pdf_fitting.R
Rscript count_cis_trans_number.R
perl merge_emplambda_xQTL.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt.gz
#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"，"QTL"#中../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz
perl min_max_pos_in_emplambda.pl #将../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 取出每条染色体上的最大值和最小值。得最大值最小值文件../output/emplambda_in_chr_min_max_position.txt
#得最大值文件../output/emplambda_in_chr_max_position.txt

#---------------------------
#---------------------------------------------
Rscript find_max_min_emplambda_in_pos.R #寻找每条染色体的 max pos和min pos 得../output/find_max_min_emplambda.txt ../output/find_max_emplambda.txt
bedtools makewindows -g ../output/find_max_emplambda.txt -w 1000 | gzip > ../output/interval_1kb.txt.gz
perl adjust_all_qtl_form.pl #调整../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 的格式，得
# ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 
bedtools intersect -a ../output/interval_1kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_1000_win.bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 5000 | gzip > ../output/interval_5kb.txt.gz
bedtools intersect -a ../output/interval_5kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_5kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 10000 | gzip > ../output/interval_10kb.txt.gz
bedtools intersect -a ../output/interval_10kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_10kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 20000 | gzip > ../output/interval_20kb.txt.gz
bedtools intersect -a ../output/interval_20kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_20kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 50000 | gzip > ../output/interval_50kb.txt.gz
bedtools intersect -a ../output/interval_50kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_50kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 100000 | gzip > ../output/interval_100kb.txt.gz
bedtools intersect -a ../output/interval_100kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_100kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 500000 | gzip > ../output/interval_500kb.txt.gz
bedtools intersect -a ../output/interval_500kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_500kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 600000 | gzip > ../output/interval_600kb.txt.gz
bedtools intersect -a ../output/interval_600kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_600kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 700000 | gzip > ../output/interval_700kb.txt.gz
bedtools intersect -a ../output/interval_700kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_700kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 800000 | gzip > ../output/interval_800kb.txt.gz
bedtools intersect -a ../output/interval_800kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_800kb_bed.gz

Rscript Point_plot_emplambda_Manhattan_pdf.R #因为cmplot画单个染色体只显示一般坐标轴，qqman又不能将多种QTL画在一条染色体上，所以尝试用ggplot模拟Manhattan
Rscript Point_plot_emplambda_Manhattan_png.R 
#------------------------------------------------
perl merge_emplambda_count.pl #将#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中../output/ALL_${QTL1}/count_number_by_emplambda_in_different_interval_7.3_all_${QTL1}.txt merge到一起，得
#得../output/count_number_by_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt
Rscript multiplot_histogram.R #画各种qtl emplambda的分布,柱状图。
Rscript circlize_plot.R #画各种qtl emplambda的分布，circlize的分布。
Rscript Point_plot_distance_emplambda.R #将 distance为x,emplambda 为y做point, ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz

#-----------------------------------------length of locus distance 
    zless /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/all_qtl_clump_locus.txt.gz |cut -f9,10,11,12,13 | sort -u -r | gzip > ../output/unique_by_qtl_tissue_pop_locus.txt.gz
    zless ../output/unique_by_qtl_tissue_pop_locus.txt.gz | cut -f1 | sort -u | wc -l #2509339 -1 = 2509338
    zless ../output/unique_by_qtl_tissue_pop_locus.txt.gz | cut -f1,2,3,5 | sort -u -r | gzip > ../output/unique_by_qtl_pop_locus.txt.gz
#
cat ../output/01_all_kinds_QTL.txt | cut -f7,9|gzip > QTL_specific_Pvalue.txt.gz







#--------------- 修改 emplambdaB.fun中plotEmp=FALSE
ALL_metaQTL
ALL_pQTL
ALL_miQTL
ALL_cerQTL
ALL_edQTL
ALL_riboQTL
52374.lncRNAQTL
ALL_caQTL
53957.reQTL
