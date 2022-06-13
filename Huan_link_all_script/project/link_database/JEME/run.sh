wget -c http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_lasso.zip
wget -c http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_elasticnet.zip
wget -c http://yiplab.cse.cuhk.edu.hk/jeme/fantom5_lasso.zip
wget -c http://yiplab.cse.cuhk.edu.hk/jeme/fantom5_elasticnet.zip

unzip encoderoadmap_lasso.zip
unzip encoderoadmap_elasticnet.zip
unzip fantom5_lasso.zip
unzip fantom5_elasticnet.zip

#ENCODE_Roadmap_info.txt famtom5_info.txt from the website
#手动annotation cell line cancer 信息得ENCODE_Roadmap_info_man.txt， famtom5_info_man.txt

perl 01_merge_normal_cell_line_enhancer_target.pl ##将$dir通过cell info 筛选得./output/01_${dir}.bed.gz，汇总得 "./output/01_merge_enhancer_target_sample.bed.gz"，将./output/01_${dir}.bed.gz 转为hg38得./output/hg38/01_${dir}_enhancer_target_sample.bed，汇总得"./output/hg38/01_merge_enhancer_target_sample.bed.gz"
