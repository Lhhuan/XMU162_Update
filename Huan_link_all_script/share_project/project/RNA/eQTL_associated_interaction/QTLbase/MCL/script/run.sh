Rscript 01_filter_TOP_QTL.R # #挑出top的QTL(hotspot),得../output/01_top_QTL_eQTL_mQTL_miQTL.txt
gzip ../output/01_top_QTL_eQTL_mQTL_miQTL.txt 
perl 02_filter_hotspot_in_1kg_phase3_v5_eur.pl #将"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz" 中的../output/01_top_QTL_eQTL_mQTL_miQTL.txt.gz 提取出来，
#得../output/02_hotspot_in_1kg_phase3_v5_eur.vcf.gz
qsub -I -l nodes=1:compute-0-6:ppn=5
cd /state/partition1/huan
plink --vcf /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/02_hotspot_in_1kg_phase3_v5_eur.vcf.gz --recode --out 02_hotspot_in_1kg_phase3_v5_eur
plink --file 02_hotspot_in_1kg_phase3_v5_eur --ld-window 1000 --ld-window-r2 0.2 --r2 --out 02_hotspot_in_1kg_phase3_v5_eur_1MB_100_0.2
plink --file 02_hotspot_in_1kg_phase3_v5_eur --ld-window 10000 --ld-window-r2 0.2 --r2 --out 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2
gzip 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld
# gzip 02_hotspot_in_1kg_phase3_v5_eur_1MB_100_0.2.ld
# zless 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz | awk -F " "  {print'$3,"\t",$6,"\t",$7'} |gzip >02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld_used_info.gz
#---compute-0-6
#-------------


cd /state/partition1/huan/tmp/

zcat ../02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz | split -l 10800000 -d -a 3 &&ls| grep x | xargs -n1 -i mv {} 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld_{}
gzip *
cd /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/script/
#将大文件分割
perl 03_annotation_ld_band_par.pl ##将用 ../data/cytoBand.txt.gz对compute-0-6下的/state/partition1/huan/tmp/下的文件 
#进行annotation,得/state/partition1/huan/tmp_output/${file}_annotation_bind.txt.gz
perl 031_cat_all_split_annotation.pl #将 "SNP_A\tSNP_B\tR2\tbind\n" 输出到 ../output/03_ALL_annotation_ld_band.txt.gz 并将 /state/partition1/huan/tmp_output/* >>../output/03_ALL_annotation_ld_band.txt.gz
# perl 03_annotation_ld_band.pl #将用 ../data/cytoBand.txt.gz对compute-0-6下的/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz 
# #进行annotation,得../output/03_annotation_ld_band.txt.gz

Rscript 04_MCL_per_band_par.R
zcat ../output/MCL_split_result/* | gzip > ../output/04_part_mcl_result.txt.gz
perl 05_merge_MCL_qtltype.pl #将../output/04_part_mcl_result.txt.gz 和../../../output/01_all_kinds_QTL.txt.gz中特定QTL取overlap得 ../output/05_merge_MCL_qtltype_result.txt.gz
--threads 20


plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID  --threads 5 /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/01_top_QTL_eQTL_mQTL_miQTL_no_header.txt.gz   --out hotspot


plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --threads 5 --blocks --out block

#-------------CLUMP
plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --threads 5 --clump-r2 0.2 --clump-kb 1000 --clump  /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/01_top_QTL_eQTL_mQTL_miQTL.txt --out hotspot



#---------------- make .map and .ped file
plink --vcf /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz --recode --out eur_1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID

#---------------make .map and .ped file by vcftools
# vcftools --gzvcf /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz --plink --out 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID

#---file plink R2
plink --file eur_1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --ld-window 100000 --r2 --out LD_result
plink --file eur_1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --ld-window 1000000 --r2 --out LD_result_1MB_0.2_1000000
plink --file eur_1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --ld-window 1000000  --ld-window-r2 0 --r2 --out LD_result_1MB_0_1000000