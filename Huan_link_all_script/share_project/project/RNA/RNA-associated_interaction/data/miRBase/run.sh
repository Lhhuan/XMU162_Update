wget -c ftp://mirbase.org/pub/mirbase/22.1/genomes/hsa.gff3 #gff文件是1 based
source activate py2
CrossMap.py gff "/home/huanhuan/tools/liftOver/hg38ToHg19.over.chain.gz" hsa.gff3 >hsa_hg38_to_hg19.gff3 #将hsa.gff3转换成hg19
source deactivate
perl 01_extract_hg19_has.pl # 提取hsa_hg38_to_hg19.gff3中的hg19,得文件01_hsa_hg19.gff3; 得转换hg19失败的hg38文件01_hsa_hg38_to_hg19_fail.gff3, gff3为1 based文件
perl 02_transform_gff_to_bed.pl #将01_hsa_hg19.gff3 转换为01_hsa_hg19.bed.gz，bed为0-based文件
