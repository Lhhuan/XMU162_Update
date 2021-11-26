perl 01_get_hotspot_fa.pl ##用 ~/ref_data/genecode/GRCh37.primary_assembly.genome.fa 为"../../../../output/ALL_eQTL/cis_trans/hotspot/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.bed.gz" 提取序列得 ,用a"${output_dir}/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.fa"

source activate /home/huanhuan/miniconda3
seekr_norm_vectors /home/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa
seekr_kmer_counts ../../../../../../output/ALL_eQTL/cis_trans/kmer/interval_18/cis_1MB/interval_18_cutoff_7.3_cis_1MB_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.fa -o 6mers.csv -mv mean.npy -sv std.npy
seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv

