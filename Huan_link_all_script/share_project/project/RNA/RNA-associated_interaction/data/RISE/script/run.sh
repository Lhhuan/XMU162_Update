cat header.txt > ../output/rise_human_all.txt
zless ../raw_data/rise_human_all.txt.gz >> ../output/rise_human_all.txt 
less ../output/rise_human_all.txt | perl -ane 'chomp; @f = split/\t/;unless($f[0]=~/\./){print "$_\n";}' |cut -f1,2,3> ../output/rise_human_not_null.txt
less ../output/rise_human_all.txt | cut -f1,2,3> ../output/RNA1_hg38.bed
less ../output/rise_human_all.txt | cut -f4,5,6> ../output/RNA2_hg38.bed
source activate py2
CrossMap.py bed  "/home/huanhuan/tools/liftOver/hg38ToHg19.over.chain.gz" ../output/RNA1_hg38.bed >../output/RNA1_hg38_to_hg19.bed 
CrossMap.py bed  "/home/huanhuan/tools/liftOver/hg38ToHg19.over.chain.gz" ../output/RNA2_hg38.bed >../output/RNA2_hg38_to_hg19.bed
source deactivate

perl 01_transform_rise_hg38_to_hg19.pl #利用../output/RNA1_hg38_to_hg19.bed和../output/RNA2_hg38_to_hg19.bed  将../output/rise_human_all.txt的gene version 转为hg19,得../output/01_rise_human_hg19.txt;
#得因为位置信息不全而不能转出hg19的文件../output/01_rise_human_position_null.txt
perl 02_normal_RNA-RNA.pl # 将../output/01_rise_human_hg19.txt 和../output/01_rise_human_position_null.txt  normal成 ../normalized/02_rise_human_hg19_normalized.txt


cat ../normalized/02_rise_human_hg19_normalized.txt| cut -f6,15 | sort -u | wc -l #112886 -1
