perl 01_count_number_and_length_of_hotspot_chr_in_lung_and_whole_blood.pl ##统计/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.g中每条染色体的hotspot的长度及对应的数目，得../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_lung_and_whole_blood.txt.gz  ../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz
perl 02_random_sampling.pl


perl 03_count_average_length_in_diff_cutoff.pl ##计算"../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz" 中cutoff下每条染色体平均hotspot的数目和长度得"../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz"
perl 04_split_chr_and_select_random.pl #根据"../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz"和"/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes" split 染色体并random取hotspot的均长，得"/share/data0/QTLbase/huan/GTEx/random_segment/mean_hotspot/${tissue}/hotspot/${i}/interval_18/${tissue}_segment_hotspot_cutoff_${Cutoff}.bed.gz"
#----单线程太慢，并行又会random同样的结果，所以先用041_prepare.pl准备好切割文件，然后在./tmp_script 进行手动多线程
# perl 041_prepare.pl



