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

perl merge_1kg_sub_pop.pl #将/share/data0/1kg_phase3_v5_hg19/五个人种的all_SNP merge在一起得../output/all_1kg_phase3_v5_hg19_snp.txt.gz

perl 01_Completion_snp_for_xQTL_by_1kg.pl # 用../output/all_1kg_phase3_v5_hg19_snp.txt.gz 补全../output/merge_QTL_all_QTLtype_pop.txt.gz 得../output/merge_QTL_all_QTLtype_pop_1kg_Completion.txt.gz
#各个文件夹进行NHPoisson
perl 011_unique_1kg_phase3_v5_snp.pl ##将"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz" 的 SNP 提出来得 在一起得../output/011_eur_1kg_phase3_v5_hg19_snp.txt.gz
Rscript 012_random_sample_snp_for_ld_block.R #../output/011_eur_1kg_phase3_v5_hg19_snp.txt.gz 每条染色体随机2个snp,重复10次,得../output/012_random_sample_snp_for_ld_block.txt.gz
gunzip ../output/012_random_sample_snp_for_ld_block.txt.gz
bash test_snp_number_in_LD.sh
perl 013_count_averge_snp_in_ld.pl # 计算../output/012_random_sample_snp_for_ld_block_r2_${i}.tags.list $i =0.5,0.7,0.8 中ld中平均包含的snp的个数得../output/013_count_averge_snp_in_ld.txt,所以取init=6
# Rscript 014_NHP.R
# Rscript 014_NHP_small.R
Rscript 014_NHP_small_par.R
Rscript 014_NHP_big_par.R
Rscript 014_NHP_eQTL.R


perl 015_filter_eQTL_chr_in_parameters.pl #将../output/ALL_eQTL/NHPoisson_emplambda_interval_${i}cutoff_7.3_all_eQTL.txt.gz 中的 chr1,chr2中的数据提出来得../output/ALL_eQTL/chr1/NHPoisson_emplambda_interval_${i}cutoff_7.3_all_eQTL.bedgraph
scp -r -P 22109 /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/  huanhuan@121.192.180.20:/public/home/huanhuan/project/bedgraph_to_wg/
Rscript 016_fitdistrbution.R

#--------------------
perl 017_filter_hotspot.pl ##@QTLs =("hQTL","mQTL","sQTL")在@interval = (6,7,8,9,12,15)时的hotspot(segment),得"../output/ALL_${QTL}/hotspot/interval_${j}_segment_hotspot.txt.gz"; 得point 热点../output/ALL_${QTL}/hotspot/interval_${j}_point_hotspot.txt.gz
perl 017_filter_hotspot_18.pl ###@QTLs =("hQTL","mQTL","sQTL")在interval = 18时的hotspot(segment),得"../output/ALL_${QTL}/hotspot/interval_${j}_segment_hotspot.txt.gz";得point 热点../output/ALL_${QTL}/hotspot/interval_${j}_point_hotspot.txt.gz
#--------------------
bash 018_bedtools_makewindows_intersect.sh 
perl 02_filter_emplambda_in_QTLbase.pl #将../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt 中的QTLbase中../output/merge_QTL_all_QTLtype_pop.txt.gz的数据过滤出来，
#得../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz
perl 02_filter_emplambda_in_QTLbase_all_par.pl # #将../output/ALL_${xQTL}/NHPoisson_emplambda_interval_${interval}cutoff_${cutoff}_all_${xQTL}.txt.gz 中的QTLbase中../output/merge_QTL_all_QTLtype_pop.txt.gz的数据过滤出来，
#得../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${interval}_cutoff_${cutoff}_all_${xQTL}.txt.gz
perl 03_merge_emplambda_xQTL.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL","cerQTL","lncRNAQTL","metaQTL","miQTL","riboQTL")#中../output/ALL_${QTL1}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL1}.txt.gz merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz
#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"，"QTL"#中../output/ALL_${QTL1}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL1}.txt merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz
#-------------------------------
perl 031_count_number_by_emplambda_in_different_interval_cutoff_all_xQTL.pl ###统计各种QTL../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${i}_cutoff_${cutoff}_all_${xQTL}.txt.gz中 中每个阶段eqtl的数目，
#得../output/ALL_${xQTL}/count_number_by_emplambda_in_different_interval_${cutoff}_all_${xQTL}.txt
#----------------------------------
Rscript 031_bar_plot_NHPoisson_all_xQTL.R #为../output/ALL_${xQTL}/count_number_by_emplambda_in_different_interval_${cutoff}_all_${xQTL}.txt画图

Rscript 04_find_max_min_emplambda_in_pos.R #根据../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz寻找每条染色体的 max pos和min pos 得../output/find_max_min_emplambda.txt ../output/find_max_emplambda.txt
perl 05_adjust_all_qtl_form.pl #调整../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz 的格式，得
# ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz

bash 06_bedtools_makewindows_intersect.sh #将../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz 按照../output/find_max_emplambda.txt区域划分../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_${X}_bed.gz
Rscript 07_multiplot_histogram.R #画各种qtl emplambda的分布,柱状图。
OmicCircos_circos_refine_revise_color.R #修改OmicCircos 的代码，可以自由定义heatmap 和heatmap2中的颜色。下次使用要进行source("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/OmicCircos_circos_refine_revise_color.R")
Rscript 08_circlize_plot_average_1MB.R #画emplambda的heatmap,用 orginal, times scale 和 max min scale, less (用的是mean emplambda)
Rscript 08_circlize_plot_average_1MB_test.R #画ALL QTL,eQTL,mQTL,sQTL,caQTL,hQTL,pQTL,edQTL(用的是mean emplambda)
Rscript 08_circlize_plot_1MB_max_replace_na.R #画全部的QTL,spread的后的na 替换为0，再进行max(用的是max emplambda)
#-----------------------------
Rscript 08_circlize_plot_1MB_max_not_replace_na.R ##画全部的QTL,spread的后的na 不替换为0，直接进行(用的是max emplambda in group), max min scale,outlier: 将前10个的值归并为11点个的值
Rscript 08_heatmap_average_1MB_per_chr_not_repalce_na.R ##画emplambda的heatmap, max min scale,outlier: 将前10个的值归并为11点个的值


perl 09_merge_emplambda_and_org_xQTL.pl #将../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
perl 09_merge_emplambda_and_org_all_QTLbase.pl  #对整个QTLbase进行操作 #xQTL 是全部的QTLbase
#将../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
#--------------------------------
perl 091_
#----------------------------

perl 10_overlap_emplambda_xQTL.pl # #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${QTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL}.txt.gz按照相同染色体位置，两两取交集，
#得../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz
perl 10_overlap_emplambda_eQTL_cis_and_trans.pl #/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz 在"cis_1MB","cis_10MB","trans_1MB","trans_10MB"
#按照相同染色体位置，与"caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"取交集得，../output/xQTL_merge/${final_type}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${final_type}_${QTL1}_${QTL2}.txt.gz

perl 10_overlap_emplambda_xQTL_cis_and_trans.pl #/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz 在"cis_1MB","cis_10MB","trans_1MB","trans_10MB"
#按照相同染色体位置，与"caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","lncRNAQTL","miQTL","metaQTL","riboQTL","cerQTL"取交集得，../output/xQTL_merge/${final_type}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${final_type}_${QTL1}_${QTL2}.txt.gz


Rscript 11_Point_plot_eQTL_cis_trans_emplambda_relevance.R #
Rscript 11_Point_plot_emplambda_relevance_pdf.R 
Rscript 11_Point_plot_emplambda_relevance_png.R #画../output/xQTL_merge/ALL/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript 11_Point_plot_emplambda_relevance_per_chr_pdf.R #分染色体画../output/xQTL_merge/ALL/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript 11_Point_plot_emplambda_relevance_per_chr_png.R
Rscript 11_Point_plot_emplambda_relevance_cis_pdf.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
Rscript 11_Point_plot_emplambda_relevance_cis_png.R
Rscript 11_Point_plot_emplambda_relevance_trans_pdf.R
Rscript 11_Point_plot_emplambda_relevance_cis_trans_pdf.R
Rscript 11_Point_plot_emplambda_relevance_cis_trans_png.R
#../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_trans_${QTL1}_cis_${QTL2}.txt.gz
#和../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_trans_${QTL1}_cis_${QTL2}.txt.gz
# ../output/xQTL_merge/Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_cis_${QTL1}_trans_${QTL2}.txt.gz
#和../output/xQTL_merge/Trans_10MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_cis_${QTL1}_trans_${QTL2}.txt.gz中QTL1和QTL2的关联


#---------------------------------------------
Rscript Manhattan_emplambda_QTL_plot_pa.R #画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure下面
Rscript Manhattan_emplambda_QTL_plot_chr_pa.R #分染色体画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure/Manhattan/下面
Rscript Manhattan_qqman_emplambda_QTL_plot_chr_pa.R  #用qqman分染色体画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure/Manhattan/下面
#--------------------------
#---------------------------
#---------------------------------------------

Rscript Point_plot_emplambda_Manhattan_pdf.R #因为cmplot画单个染色体只显示一般坐标轴，qqman又不能将多种QTL画在一条染色体上，所以尝试用ggplot模拟Manhattan
Rscript Point_plot_emplambda_Manhattan_png.R 
#--------------------------------------------------------





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



#----------------------