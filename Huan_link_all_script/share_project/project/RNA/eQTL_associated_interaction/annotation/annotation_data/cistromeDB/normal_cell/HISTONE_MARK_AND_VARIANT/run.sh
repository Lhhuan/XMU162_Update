#------------------------------------------unsed
#--------------------------------------
#----------------------------------------
bedtools makewindows -b  merge_pos_info_narrow_peak_complement.bed.gz -w 379 | gzip >merge_pos_info_narrow_peak_complement_split.bed.gz


ln merge_pos_info_narrow_peak_complement_split.bed.gz  /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/non_HISTONE_modification_split.bed.gz
ln merge_pos_info_narrow_peak_sort.bed.gz  /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/HISTONE_modification.bed.gz


bedtools makewindows -b  merge_pos_info_narrow_peak_union_complement.bed.gz -w 713 | gzip >merge_pos_info_narrow_peak_union_complement_split.bed.gz
ln merge_pos_info_narrow_peak_union_complement_split.bed.gz  /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/non_HISTONE_modification_split_union.bed.gz
ln merge_pos_info_narrow_peak_sort_union.bed.gz  /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/HISTONE_modification_union.bed.gz