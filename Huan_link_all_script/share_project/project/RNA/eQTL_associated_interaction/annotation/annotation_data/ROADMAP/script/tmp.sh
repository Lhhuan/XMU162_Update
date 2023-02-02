perl 04_merge_eQTL_Catalogue_marker.pl 
echo -e "04_merge_eQTL_Catalogue_marker\n"
perl 05_tranform_hg19_to_hg38.pl 
echo -e "05_tranform_hg19_to_hg38\n"

perl 04_merge_eQTL_Catalogue_marker_sample.pl 
echo -e "04_merge_eQTL_Catalogue_marker_sample\n"
perl 05_tranform_hg19_to_hg38_sample.pl
echo -e "05_tranform_hg19_to_hg38_sample\n"

perl 04_merge_eQTL_Catalogue_marker_signalValue.pl 
echo -e "04_merge_eQTL_Catalogue_marker_signalValue\n"
perl 05_tranform_hg19_to_hg38_signalValue.pl 
echo -e "05_tranform_hg19_to_hg38_signalValue\n"

# source activate huan_py3
# Rscript 06_intersect_merge_bed_and_signalValue.R
# echo -e "05_tranform_hg19_to_hg38_signalValue\n"