bedtools getfasta -fi /share/Projects/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed  "/home/huanhuan/project/GTEx/output/Tissue_merge/hotspot_quantile/org_hotspot_quantile25_75.bed" -fo "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/quantile25_75/kmer/6/org_hotspot_quantile25_75.fa"

cd /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/quantile25_75/kmer/6/

source activate seekr_source
seekr_norm_vectors org_hotspot_quantile25_75.fa
seekr_kmer_counts org_hotspot_quantile25_75.fa -o 6mers.csv -mv mean.npy -sv std.npy

seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv

seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv
cp communities.csv communities_5.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 10 -s 0 -c communities_10.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 9 -s 0 -c communities_9.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 8 -s 0 -c communities_8.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 7 -s 0 -c communities_7.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 6 -s 0 -c communities_6.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 4 -s 0 -c communities_4.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 3 -s 0 -c communities_3.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 2 -s 0 -c communities_2.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 5 -s 0 -c communities_5_specify.csv


conda deactivate


cd /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/