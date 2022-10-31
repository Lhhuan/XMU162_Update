wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/files.txt
perl 01_filter_need_file.pl #将files.txt中未经其他处理的bed选出来，对于同一样本有重复的，筛选出文件大的哪个，对于并生成download link,得download.sh,和01_need_file.txt 和对应的细胞文件01_need_cell.txt
cd raw_data
bash ../download.sh 
cd ..
perl 02_filter_normal_cell_line.pl 