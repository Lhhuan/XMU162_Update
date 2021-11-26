# perl 01_filter_emplambda_0_0.176_segment.pl ###过滤出../../../../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 中emplambda$emplambda >0 && $emplambda<0.176，得 "../../../../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/sampling/whole_blood_segment_hotspot_cutoff_more_0_and_less_than_0.176.bed.gz"
# perl 02_random_genomic_resemble_hotspot.pl #产生10000个与../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}",
# #whole_blood_segment_hotspot_cutoff_more_0_and_less_than_0.176.bed.gz

perl 01_filter_emplambda_0_0.176_segment_fact.pl ###过滤出../../../../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 中$emplambda<0.176，得 "../../../../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/sampling/whole_blood_segment_hotspot_cutoff_0_0.176.bed.gz

perl 02_random_genomic_resemble_hotspot_fact.pl  #产生10000个与../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}",
#whole_blood_segment_hotspot_cutoff_0_0.176.bed.gz
perl 03_annotation_marker.pl #对10000个../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2/*进行histone marker的注释，得/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2/*
perl 04_count_annotation_marker.pl #对/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/$type/$cutoff2/*的marker进行count,得"$out_dir/${type}_${cutoff2}_histone_marker.txt.gz";

perl 05_calculate_jaccard_index.pl #对$input_dir/${mark}_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz 进行jaccard index 进行计算,得的marker进行计算,得"$out_dir/${cutoff2}_jaccard_index_histone_marker.txt.gz";
perl 05_calculate_jaccard_index_1000.pl 


perl 06_annotation_chrom_state_emplambda_0_0.176_par.pl #对10000background_original_random/0_0.176进行chromatin_states的注释，得"/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2"/${i}_15_state_resemble_${input_file_base_name}
perl 07_filter_annotation_chrom_state_par.pl #因为每个片段会被多个chrom state 注释，根据对hotspot的覆盖比例，选出覆盖比例最高的state, interval_18 时，对"/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2/${state}_state_resemble_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"进行过滤，得$dir/filter_$name
"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/random_permutation/34_annotation_chrom_state_original_random_and_permutation_hotspot_par.pl"
zless 


