# #------------from mulin

#     # bcftools concat 1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr2.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr3.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr4.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr5.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr6.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr7.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr8.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr9.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr10.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr11.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr12.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr13.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr14.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr15.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr16.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr17.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr18.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr19.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr20.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr21.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr22.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chrX.vcf.gz -o 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz -O z

#     # bcftools view -v snps 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz
#     # bcftools annotate --set-id '%CHROM\_%POS' 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.IDS.vcf.gz
#     # bcftools norm -d snps 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.IDS.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.IDS.uniq.vcf.gz

#     # /f/Tools/plink2/plink-1.90/plink --bfile 1kg.phase3.v5.afr.hg19.all --exclude dup.remove --make-bed --out 1kg.phase3.v5.afr.hg19.rdup.all

#     # /f/Tools/plink2/plink-1.90/plink --vcf 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.IDS.uniq.vcf.gz --maf 0.05 --make-bed --out 1kg.phase3.v5.eur.hg19.SNPs
#     # /f/Tools/plink2/plink-1.90/plink --bfile 1kg.phase3.v5.eur.hg19.SNPs --show-tags demo/test.proxy --list-all --tag-r2 0.8

#     # /f/Tools/plink2/plink-1.90/plink --bfile 1kg.phase3.v5.eur.hg19.SNPs.0.05 --clump-p1 1e-5 --clump-p2 0.01 --clump-r2 0.8 --clump-kb 250 --clump test.proxy.signal --clump-snp-field ID -clump-field PVAL --out test.proxy.signal
# #-------------------------------------------------------------
perl 01_filter_EUR_sample_name.pl #将/home/chaoqun/phase3/integrated_call_samples_v3.20130502.ALL.panel 中的EUR提出来，得EUR_sample_list.txt
# perl 02_extract_EUR_samples_in_vcf_format.pl #将/home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 中提取出EUR 的sample info，得1kg.phase3.v5.shapeit2.eur.hg19.chr${i}.vcf.gz
# wget -c http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz
# wget -c http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz

# bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr1._test_vcf.gz -O z
#_------------------------------------------------------------------------------------- #perl 02_extract_EUR_samples_in_vcf_format.pl  中多线程和bash test.sh 多线程 提EUR局费时且不准确，so,手动多线程
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr2.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr3.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr4.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr5.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr6.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr7.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr8.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr9.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr10.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr11.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr12.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr13.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr14.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr15.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr16.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr17.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr18.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr19.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr20.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr21.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "/home/chaoqun/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -o 1kg.phase3.v5.shapeit2.eur.hg19.chr22.vcf.gz -O z &

#---------------


bcftools concat --threads 20 1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr2.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr3.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr4.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr5.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr6.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr7.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr8.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr9.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr10.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr11.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr12.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr13.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr14.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr15.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr16.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr17.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr18.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr19.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr20.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr21.vcf.gz 1kg.phase3.v5.shapeit2.eur.hg19.chr22.vcf.gz -o 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz -O z 
#--------------------------------------------------snp
bcftools view -v snps --threads 20 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz
bcftools norm -d snps --threads 20 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq.vcf.gz  #去重

plink --vcf 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq.vcf.gz --make-bed --out 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq
Rscript makeSNPnamesUnique.R 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq.bim
#------------------------------------------------------------------------all
bcftools norm -d snps --threads 20 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all_uniq.vcf.gz

plink --vcf 1kg.phase3.v5.shapeit2.eur.hg19.all_uniq.vcf.gz --make-bed --out 1kg.phase3.v5.shapeit2.eur.hg19.all_uniq

Rscript makeSNPnamesUnique.R 1kg.phase3.v5.shapeit2.eur.hg19.all_uniq.bim


# #------------------------------------------------------- chr pos as ID
bcftools annotate --set-id '%CHROM\_%POS'  1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz
plink --vcf 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz --make-bed --out 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID
Rscript makeSNPnamesUnique.R 1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.bim
# plink --vcf 1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz --make-bed --out 1kg.phase3.v5.shapeit2.eur.hg19.all
# plink --bfile 1kg.phase3.v5.eur.hg19.SNPs 