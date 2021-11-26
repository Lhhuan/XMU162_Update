perl 01_filter_cluster0.pl  #/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/
bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster0/cluster0.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster0/cluster0.fa"
cd  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster0/

source activate /home/huanhuan/miniconda3
seekr_norm_vectors ../Whole_Blood_segment_hotspot_cutoff_0.176.fa
seekr_kmer_counts cluster0.fa -o 6mers.csv -mv mean.npy -sv std.npy
seekr_pearson 6mers.csv 6mers.csv -o 6example_vs_self.csv
seekr_graph 6example_vs_self.csv 0.13 -g 6example_vs_self.gml -c 6communities.csv

#---------

seekr_norm_vectors ../Whole_Blood_segment_hotspot_cutoff_0.176.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy
seekr_kmer_counts cluster0.fa -o 5mers.csv -k 5 -mv mean_5mers.npy -sv std_5mers.npy
seekr_pearson 5mers.csv 5mers.csv -o 5example_vs_self.csv
seekr_graph 5example_vs_self.csv 0.13 -g 5example_vs_self.gml -c 5communities.csv


#------------------------------------
conda deactivate


perl  02_transform_kmer_result.pl

Rscript  03_heatmap.R