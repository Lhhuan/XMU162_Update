bcftools concat 1kg.phase3.v5.shapeit2.eas.hg19.chr1.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr2.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr3.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr4.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr5.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr6.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr7.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr8.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr9.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr10.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr11.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr12.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr13.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr14.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr15.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr16.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr17.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr18.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr19.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr20.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr21.vcf.gz 1kg.phase3.v5.shapeit2.eas.hg19.chr22.vcf.gz -o /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.vcf.gz -O z 
echo -e "bcftools concat\n"
bcftools view -v snps  /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.vcf.gz -O z -o /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz
echo -e "bcftools view -v snps  1kg.phase3.v5.shapeit2.eas.hg19.all.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz\n"
bcftools norm -d snps  /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz -O z -o /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.vcf.gz  #去重
echo -e "bcftools norm -d snps  1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.vcf.gz\n"
bcftools annotate --set-id '%CHROM\_%POS' /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.vcf.gz -O z -o /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.vcf.gz
echo -e "bcftools annotate --set-id '%CHROM\_%POS' 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.vcf.gz -O z -o 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.vcf.gz\n"
plink --vcf /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.vcf.gz --make-bed --out /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID
echo -e "plink --vcf 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.vcf.gz --make-bed --out 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID\n"
Rscript makeSNPnamesUnique.R /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.bim
echo -e "Rscript makeSNPnamesUnique.R 1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID.bim\n"