perl 01_normal_format.pl #将../raw_data/Cis-eQTLs.txt 和../raw_data/Trans-eQTLs.txt 进行初步normal, ../output/01_normal_format.txt.gz
#并得unique的gene POS的bed文件../output/01_unique_gene_position_hg38.bed.gz
source activate py2
CrossMap.py bed "/home/huanhuan/tools/liftOver/hg38ToHg19.over.chain.gz"  ../output/01_unique_gene_position_hg38.bed.gz ../output/01_unique_gene_position_hg19.bed
source deactivate
perl 02_merge_hg19_position_and_normal_data.pl #将../output/01_unique_gene_position_hg19.bed和../output/01_normal_format.txt.gz merge 在一起，得../normalized/02_normalized.txt.gz
zless ../normalized/02_normalized.txt.gz | wc -l #6661690 -1