perl 01_extend_segment_6.pl #将"../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz"; #长度<6的片段为扩展为6或者7，../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6/extend_whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz
bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/extend_Whole_Blood_segment_hotspot_cutoff_0.176.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa"
cd  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/

source activate /home/huanhuan/miniconda3
seekr_norm_vectors extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa
seekr_kmer_counts extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa -o 6mers.csv -mv mean.npy -sv std.npy

seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv

seekr_kmer_counts extend_Whole_Blood_segment_hotspot_cutoff_0.176.fa -o 6mers_uc_us_no_log.csv --log2 none -uc -us
gzip 6mers_uc_us_no_log.csv
conda deactivate


# perl 02_random_genomic_resemble_extend_hotspot.pl
perl 03_transform_kmer_result_and_anno.pl
Rscript 05_Pie.R
Rscript 06_heatmap.R
Rscript 07_count_kmer_count.R

