perl 01_adjust_tissue_merge.pl #判断"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" 中来自不同study 的相同tissue id,组织名字，roadmap,cistromeDB等信息是否完全相同,判断结果是相同
perl 02_tissue_level_Marker_source.pl #根据"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" 的 tissue_ontology_id在"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_roadmap_info.tsv" 和"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv"提取tissue 的annotation info,得./output/02_tissue_level_Marker_source.txt
perl 03_tissue_specific_marker_merge.pl #
Rscript 04_intersect_merge_bed_and_signalValue.R
perl 05_computer_marker_matrix_na0.pl
Rscript 06_extract_tissue_feature.R 
python 07_predict_normal.py
perl 08_split_prediction_result.py 
perl 09_plot_marker_matrix_na0_predict.pl 
perl 10_merge_normal_tissue_eqtl.pl
perl 11_bedtools_eqtl.pl 
Rscript 12_barplot_eqtl_distribution.R
Rscript 13_boxplot_gc_content_distribution.R
perl 14_enrichment_ChromHMM.pl 
perl 15_enrichment_cCREs.pl
perl 16_enrichment_HAQER.pl 
#========================overlapenrichment
source activate huan_py3
cd /home/huanhuan/tools/gonomics/src/github.com/vertgenlab/gonomics-main/cmd/overlapEnrichments
go run overlapEnrichments.go exact ./testdata/elements1.bed ./testdata/elements2.bed ./testdata/tinyNoGap.bed ./testdata/123.bed




Rscript 09_eqtl_distribution.R
perl 10_normal_tissue_distribution_in_marker.pl 
perl 