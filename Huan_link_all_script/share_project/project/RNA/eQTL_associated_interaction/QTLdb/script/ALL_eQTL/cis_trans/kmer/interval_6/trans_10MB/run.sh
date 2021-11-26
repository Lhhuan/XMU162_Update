perl 01_get_hotspot_fa.pl ##用 ~/ref_data/genecode/GRCh37.primary_assembly.genome.fa 为"../../../../output/ALL_eQTL/cis_trans/hotspot/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.bed.gz" 提取序列得 ,用a"${output_dir}/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.fa"

seekr_norm_vectors /home/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa
seekr_kmer_counts ../../../../../output/ALL_eQTL/cis_trans/kmer/interval_6/trans_10MB/interval_6_cutoff_7.3_trans_10MB_eQTL_segment_hotspot_length_more_than6.fa -o 6mers.csv -mv mean.npy -sv std.npy
seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv

