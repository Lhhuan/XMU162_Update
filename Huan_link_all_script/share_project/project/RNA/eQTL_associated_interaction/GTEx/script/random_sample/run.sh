perl 01_count_number_of_eQTL_chr_in_lung_and_whole_blood.pl #统计lung和 whole_blood ${dir}/${tissue}${suffix}中每条染色体的eQTL的数目，得../../output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz
Rscript 02_random_sampling.pl #根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz"，对应染色体从"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz"随机抽取相同数目的位置作为新的QTL,并将其他位置定义为非qtl得"$output_dir/${tissue}_random_select_result_${i}.txt.gz"
Rscript 03_NHP_big_par.R

Rscript 04_histgram_density_interval.R
# Rscript 05_histgram_density_interval_chr_specific.R

#-------------------以hotspot 为基准计数
perl 06_filter_hotspot_for_interval18.pl ###"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/random_select_result/${tissue}_random_select_result_${number}.txt.gz" 时得不同cutoff下的hotspot(segment),得"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/hotspot/${number}/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz"
Rscript 061_histgram_plot_the_length_of_hotspot.R 


perl 20_count_random_hit_eQTL_ratio.pl # 统计"$output_dir/${tissue}_random_select_result_${i}.txt.gz" 中，hit住eQTL的数目，得count_random_hit_eQTL_ratio.txt
Rscript 21_count_segment_length.R