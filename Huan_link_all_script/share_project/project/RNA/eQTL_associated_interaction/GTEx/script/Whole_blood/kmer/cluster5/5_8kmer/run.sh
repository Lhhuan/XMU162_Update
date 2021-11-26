perl 01_filter_segment_8.pl
perl 01_filter_hotspot_segment_8.pl


bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/8kmer/cluster5.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/8kmer/cluster5.fa"
bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/8kmer/Whole_Blood_segment_hotspot_cutoff_0.176.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/8kmer/hotspot_0.176.fa"


cd  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/8kmer
source activate /home/huanhuan/miniconda3
seekr_norm_vectors hotspot_0.176.fa -k 8 -mv mean_8mers.npy -sv std_8mers.npy
# seekr_norm_vectors hotspot_0.176.fa -k 7 -mv mean_7mers.npy -sv std_7mers.npy
seekr_kmer_counts cluster5.fa -o 8mers.csv -k 8 -mv mean_8mers.npy -sv std_8mers.npy
seekr_pearson 8mers.csv 8mers.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv
conda deactivate

#----------6mers
seekr_norm_vectors hotspot_0.176.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy
# seekr_norm_vectors hotspot_0.176.fa -k 7 -mv mean_7mers.npy -sv std_7mers.npy
seekr_kmer_counts cluster5.fa -o 6mers.csv -k 6 -mv mean_6mers.npy -sv std_6mers.npy
seekr_pearson 6mers.csv 6mers.csv -o 6example_vs_self.csv
seekr_graph 6example_vs_self.csv 0.13 -g 6example_vs_self.gml -c 6communities.csv

#----------7mers
seekr_norm_vectors hotspot_0.176.fa -k 7 -mv mean_7mers.npy -sv std_7mers.npy
# seekr_norm_vectors hotspot_0.176.fa -k 7 -mv mean_7mers.npy -sv std_7mers.npy
seekr_kmer_counts cluster5.fa -o 7mers.csv -k 7 -mv mean_7mers.npy -sv std_7mers.npy
seekr_pearson 7mers.csv 7mers.csv -o 7example_vs_self.csv
seekr_graph 7example_vs_self.csv 0.13 -g 7example_vs_self.gml -c 7communities.csv

#----------5mers
cd /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/5mer/ 


seekr_norm_vectors hotspot_0.176.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy
seekr_kmer_counts cluster5.fa -o 5mers.csv -k 5 -mv mean_5mers.npy -sv std_5mers.npy
seekr_pearson 5mers.csv 5mers.csv -o 5example_vs_self.csv
seekr_graph 5example_vs_self.csv 0.13 -g 5example_vs_self.gml -c 5communities.csv
conda deactivate
#---------

mv 5* ../5mer
#------------

cd /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/cluster5/5_8kmer/
perl 02_transform_kmer_result.pl 


Rscript 06_heatmap.R