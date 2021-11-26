perl 01_normal_format.pl ##将../raw_data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt 和../raw_data/2012-12-21-TransEQTLsFDR0.5.txt 进行初步normal, 得gene name 和部分RSid从对应的v8.egenes.txt.gz中获得，得文件../output/01_normal_format.txt.gz
#并得unique的varint POS的vcf文件../output/01_unique_variant-id_position_hg18.vcf.gz
source activate py2
CrossMap.py vcf "/home/huanhuan/tools/liftOver/hg38ToHg19.over.chain.gz" ../output/01_unique_variant-id_position_hg18.vcf.gz /home/huanhuan/tools/liftOver/hg19.fa ../output/01_unique_variant-id_position_hg19.vcf
source deactivate
perl 02_merge_hg19_position_and_normal_data.pl #将../output/01_unique_variant-id_position_hg19.vcf和../output/01_normal_format.txt.gz merge 在一起，得../normalized/02_normalized.txt.gz
zless ../normalized/02_normalized.txt.gz | wc -l #926837 -1