perl 01_merge_all_eQTL_association_interaction.pl ##将../data/* 下面normalizeda好的文件合在一起,得../output/01_all_eQTL_associated_interaction.txt.gz
zless ~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz | wc -l # 670774443
zless ~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz | cut -f3 |sort -u| wc -l # 645507394
zless ~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz | cut -f3|sort |uniq -d #
perl 02_add_pos_info_for_rsid.pl #将用"~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz" 为 ../output/01_all_eQTL_associated_interaction.txt.gz填充position和ref 和alt的信息
#得../output/02_add_pos_info_for_rsid.txt.gz
perl 03_add_rsid.pl #将用"~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz" 为 ../output/02_add_pos_info_for_rsid.txt.gz填充position和ref 和alt的信息 #得../output/03_add_rsid.txt.gz
perl 04_add_ENSG_Entrez.pl #用"~/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt" 为../output/03_add_rsid.txt.gz添加 ENSG和Entrez,得../output/04_add_gene_ENSG_and_Entrez.txt.gz
perl 05_add_gene_position.pl #用"/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt"， "/home/huanhuan/ref_data/RNA_position/data/hg19_v19_miRNA_position.txt"，
#"/home/huanhuan/ref_data/RNA_position/data/hg19_mRNA_position.txt"和"/home/huanhuan/ref_data/RNA_position/data/gencode.v32lift37.long_noncoding_RNAs.gff3.gz"为../output/04_add_gene_ENSG_and_Entrez.txt.gz 寻找gene的
#位置信息，得../output/05_add_gene_position.txt.gz
zless ../output/05_add_gene_position.txt.gz | cut -f2,3,21,22,24,25,26 >../output/05_used_clump.txt
perl 06_plink_clump.pl #将../output/05_used_clump.txt 按照source，tissue和cis or trans进行分割，得../output/used_to_clump/${database}_${tissue}_cis.eSNP， ../output/used_to_clump/${database}_${tissue}_trans.eSNP
#../output/used_to_clump/${database}_${tissue}_unknown_cis_trans.eSNP 
#然后进行plink clunp,得../output/clump/${database}_${tissue}_cis.eSNP, ../output/clump/${database}_${tissue}_trans.eSNP, ../output/clump/${database}_${tissue}_unknown_cis_trans.eSNP
#并得database与tissue对应的文件../output/06_database_and_tissue_info.txt
    #----------------------------------------因为数据库内组织较多，需要多次循环，而又不确定多线程结果是否准确，所以分为三个sub进行
    perl 06_plink_clump_sub1.pl #进行除PancanQTL和ncRNA_eQTL的其他数据库的计算
    perl 06_plink_clump_sub2.pl ##进行PancanQTL的计算
    perl 06_plink_clump_sub3.pl #进行ncRNA_eQTL的计算
#--------------------------------
perl 07_cut_clump.pl #将../output/clump/*.clumped cut 出snp和与其相关的SNP，得../output/final_clump/*.clumped #得全部的final_clump文件 ../output/07_all_clump_result.txt


    plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq --threads 20 --clump-r2 0.2 --clump-kb 500 --clump /home/chaoqun/phase3/CLUMP/trans/lncRNA_trans_THCA_22.eSNP --out ../output/clump/test_chaoqun_lncRNA_trans_THCA_22.eSNP