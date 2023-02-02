perl 061_annotation_for_predicted.pl
# bash annotation_marker_interval18_signalValue_mean.sh
Rscript 062_mean_marker_signalValue_for_predicted.R
Rscript 07_exact_feature_for_predicted.R
python 08_predict_warm_region.py

perl 09_eqtl_to_bed.pl
bedtools intersect -a ../output/09_eqtl_sorted.bed.gz -b "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/warm_hotspot_region_win4518_large_than6.bed.gz" -wa -wb |gzip > ../output/eqtl_warm_hotspot_region_win4518_large_than6.bed.gz

Rscript 10_plot_eqtl_distribution.R
perl 11_plot_marker_peak_point.pl