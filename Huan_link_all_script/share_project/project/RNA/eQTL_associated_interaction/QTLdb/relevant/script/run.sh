bash eQTL_overlap_hotspot_bedtools_intersect.sh
bash eQTL_interval_15_non_hotspot_bedtools_intersect.sh
bash eQTL_overlap_non_hotspot_bedtools_intersect.sh
bash eQTL_interval_15_hotspot_bedtools_intersect.sh


Rscript 01_merge_all_annotation.R #将../output/annotation 下的文件merge 到一起