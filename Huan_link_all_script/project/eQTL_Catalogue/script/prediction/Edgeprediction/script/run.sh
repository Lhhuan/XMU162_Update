perl 01_adjust_reactomeFI_format.pl #调整"~/project/link_database/reactome_FI/output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt"格式得 ../output/train/01_reactome.csv.gz
perl 02_adjust_all_hotspot_format.pl ##调整"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"格式得 "../output/01_all_hotspot.csv.gz"
Rscript 021_split_train_and_valid.R  #"./test/021_hotsopt_valid.csv.gz"
perl 022_adjust_train_valid.pl #adjust ../output/test/021_hotsopt_valid.csv.gz to json得"../output/test/02_valid_hotspot.json.gz"，"../output/02_all_hotspot.csv.gz"减去"../output/test/021_hotsopt_valid.csv.gz"得 "../output/train/02_train_hotspot.csv.gz"
less "../output/train/02_train_hotspot.csv" |head -563046 > "../output/train/02_train1_hotspot.csv"
python 03_train_validation_models.py 
Rscript 04_plot_AUC_distribution.R 






