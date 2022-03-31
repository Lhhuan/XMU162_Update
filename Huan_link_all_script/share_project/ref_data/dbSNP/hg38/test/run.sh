wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz
wget -c ftp://ftp.ncbi.nih.gov/snp/redesign/archive/b152/VCF/GCF_000001405.25.bgz

#----------------------------
mv common_all_20180423.vcf.gz dbsnp_151.hg19.vcf.gz

perl 01_extract_b152_vcf.pl #将GCF_000001405.25.bgz 中所需要的 id,chr,pos,ref,alt 提取出来，得01_extract_b152_vcf.txt.gz




#----------------------
