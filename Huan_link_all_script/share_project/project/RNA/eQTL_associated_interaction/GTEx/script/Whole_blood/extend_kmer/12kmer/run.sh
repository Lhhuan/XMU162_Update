
cd  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/12/kmer/extend/

bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/extend_kmer/12kmer/extend_Whole_Blood_segment_hotspot_cutoff_0.176.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/12/kmer/extend/extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa"
source activate /home/huanhuan/miniconda3
seekr_norm_vectors  -k 12 extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa 
# seekr_kmer_counts extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa -o 12mers.csv -k 12 -mv mean.npy -sv std.npy

# seekr_pearson 12mers.csv 12mers.csv -o example_vs_self.csv 
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv

# # seekr_kmer_counts extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa -o 6mers_uc_us_no_log.csv --log2 none -uc -us
# # gzip 6mers_uc_us_no_log.csv
# conda deactivate


# # perl 02_random_genomic_resemble_extend_hotspot.pl
# perl 03_transform_kmer_result_and_anno.pl
# Rscript 05_Pie.R


# # Rscript 06_heatmap.R
# # Rscript 07_count_kmer_count.R

