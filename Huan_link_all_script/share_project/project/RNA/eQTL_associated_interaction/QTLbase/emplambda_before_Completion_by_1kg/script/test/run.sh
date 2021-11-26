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
#
perl count_qtl_per_source.pl # #统计../output/all_qtl_clump_locus_r_square0.5.txt.gz中pop_tissue_qtl的数目，得../output/number_of_pop_tissue_qtl.txt, 排序得../output/sorted_number_of_pop_tissue_qtl.txt
#得mQTL    Blood   EUR数目最多为470179
cat ../output/number_of_pop_tissue_qtl.txt | sort -k4nr >../output/sorted_number_of_pop_tissue_qtl.txt #
perl count_per_region_QTL_number.pl #统计../output/all_qtl_clump_locus_r_square0.5.txt.gz 中每特定距离中QTL(tag snp的数量)，得../output/count_per_100kb_QTL_number.txt
perl count_per_region_snp_number.pl #统计../output/01_all_kinds_QTL.txt 中每特定距离中snp的数量，得../output/count_per_500bp_snp_number.txt, ../output/count_per_750bp_snp_number.txt, ../output/count_per_1000bp_snp_number.txt
#../output/count_per_500bp_snp_number_region.txt, ../output/count_per_750bp_snp_number_region.txt, ../output/count_per_1000bp_snp_number_region.txt
zless "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/data/all.txt.gz" | perl -ane "print $f1"
perl filter_chr22.pl ##将../output/01_all_kinds_QTL.txt 中chr22过滤得../output/chr22.txt,过滤chr22中eur得../output/chr22_eur.txt
#得所有的具有 pop 和tissue的文件得../output/merge_QTL_all_QTLtype_pop.txt.gz
perl judge_xQTL_cis_trans.pl ##将../output/merge_QTL_all_QTLtype_pop.txt.gz 按照1MB和10MB划分cis和trans， 分别得../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz
perl region_QTL.pl #

Rscript per_region_QTL_number.R #
Rscript QTL_density.R
Rscript NHPoisson.R 
Rscript huan_NHPoisson.R
Rscript Manhattan_plot.R #
Rscript density_NHPoisson.R
perl count_number_in_different_interval.pl 

#------------------- #EUR eQTL 
Rscript huan_NHPoisson.R
Rscript density_NHPoisson_batch.R
perl count_number_by_emplambda_in_different_interval_cutoff7.pl # #统计../output/NHPoisson_POTevents_inddat_cutoff_${i}.txt中每个阶段QTL的数目，得../output/count_number_by_emplambda_in_different_interval_cutoff_7.txt
perl count_number_by_emplambda_in_different_interval_cutoff9.pl
Rscript bar_plot_NHPoisson.R

#--------------------------#ALL QTL, population specific
Rscript huan_NHPoisson_all_QTLbase_pop_specific.R
perl count_number_by_emplambda_in_different_interval_cutoff7_all_QTL_pop_specific.pl # #统计../output/population_specific/NHPoisson_POTevents_inddat_cutoff_${cutoff}_all_QTLbase_${pop}_specific.txt中每个阶段qtl的数目，
#得../output/population_specific/count_number_by_emplambda_in_different_interval_cutoff_7_all_QTL_${pop}_specific.txt
Rscript density_NHPoisson_batch_all_QTLbase_pop_specific.R
Rscript bar_plot_NHPoisson_all_QTL_pop_specific.R

#----------------------- #ALL QTL, 全部的QTL，不分种族 
Rscript huan_NHPoisson_all_QTLbase.R 
Perl count_number_by_emplambda_in_different_interval_cutoff7_all_QTL.pl ##统计../output/NHPoisson_emplambda_interval_${i}_all_QTLbase.txt中 中每个阶段qtl的数目，
#得../output/count_number_by_emplambda_in_different_interval_cutoff_7_all_QTL.txt
Rscript density_NHPoisson_batch_all_QTLbase.R
Rscript bar_plot_NHPoisson_all_QTL.R

#----------------------- #ALL eQTL, 全部的eQTL，不分种族
Rscript huan_NHPoisson_all_eQTLbase.R 
Rscript huan_NHPoisson_all_eQTLbase_p.R  #huan_NHPoisson_all_eQTLbase.R 的多线程
Perl count_number_by_emplambda_in_different_interval_cutoff7_all_eQTL.pl ##统计../output/NHPoisson_emplambda_interval_${i}_all_eQTLbase.txt中 中每个阶段eqtl的数目，
#得../output/count_number_by_emplambda_in_different_interval_cutoff_7_all_eQTL.txt
Rscript density_NHPoisson_batch_all_eQTLbase.R
Rscript bar_plot_NHPoisson_all_eQTL.R
#---------------------------------------

perl merge_emplambda_xQTL.pl #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中按照位置，两两取交集，得../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz
Rscript Manhattan_emplambda_QTL_plot_pa.R #画xQTL 中 emplambda的Manhattan，分别存在各个xQTL figure下面
perl merge_emplambda_and_org_xQTL.pl #将../output/xQTL/count_number_by_emplambda_in_different_interval_7.3_all_xQTL.txt 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
perl merge_emplambda_and_org_all_QTLbase.pl #对整个QTLbase进行操作 #xQTL 是全部的QTLbase
#将../output/xQTL/count_number_by_emplambda_in_different_interval_7.3_all_xQTL.txt 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
Rscript point_Plot_emplambda_distance.R #

#-----------------------------------------length of locus distance 
    zless /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/all_qtl_clump_locus.txt.gz |cut -f9,10,11,12,13 | sort -u -r | gzip > ../output/unique_by_qtl_tissue_pop_locus.txt.gz
    zless ../output/unique_by_qtl_tissue_pop_locus.txt.gz | cut -f1 | sort -u | wc -l #2509339 -1 = 2509338
    zless ../output/unique_by_qtl_tissue_pop_locus.txt.gz | cut -f1,2,3,5 | sort -u -r | gzip > ../output/unique_by_qtl_pop_locus.txt.gz

    Rscript locus_distance_density_and_point.R  # 画../output/all_qtl_clump_locus.txt.gz 的distnace 的density plot 和point plot,画point plot失败
    Rscript unique_locus_distance_density_and_point.R  # 画../output/unique_by_qtl_tissue_pop_locus.txt.gz 和 ../output/unique_by_qtl_pop_locus.txt.gz 的distance 的density plot 和point plot, 画point plot失败

    #----------------------------------------

    perl count_of_different_distance_of_region.pl #统计../output/unique_by_qtl_tissue_pop_locus.txt.gz中不同区段长度的数目，得../output/count_of_different_distance_of_region.txt

    Rscript bar_plot_distance_of_region.R 






#--------------------------------------------------------------count p value
Rscript multi_unique_locus_distance_density_and_point.R
#
cat ../output/01_all_kinds_QTL.txt | cut -f7,9|gzip > QTL_specific_Pvalue.txt.gz
Rscript pvalue_density.R


#
perl analysis_json.pl






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



#--------------------------running 
55665.ALL_eQTL
59165.mQTL


#--------------------------finish
29022.sQTL
9146.reQTL #数据少
9316.metaQTL # no data
11732.pQTL
57860.cerQTL # no data
61861.ALL_edQTL #data 很少
31551.miQTL # no data
14626.riboQTL # no data
46226.ALL_caQTL
52374.lncRNAQTL # no data
ALL_hQTL