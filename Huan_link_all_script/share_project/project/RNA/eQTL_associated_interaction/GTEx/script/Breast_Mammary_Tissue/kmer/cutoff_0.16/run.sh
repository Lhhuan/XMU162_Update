perl 06_filter_segment_6.pl
bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/Lung_segment_hotspot_cutoff_0.16.bed.gz" -fo "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/Lung_segment_hotspot_cutoff_0.16.fa"


source activate /home/huanhuan/miniconda3

seekr_kmer_counts /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/Lung_segment_hotspot_cutoff_0.16.fa -o 6mers.csv -mv mean.npy -sv std.npy

seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv