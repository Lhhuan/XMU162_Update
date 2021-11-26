# bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed  "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/chr1_11/Chr1_11_Tissue_merge_0.176_extend_sorted_merge.bed.gz" -fo "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.fa"

cd /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/test_kkl/

source activate /home/huanhuan/miniconda3
# seekr_norm_vectors Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.fa
# seekr_kmer_counts Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.fa -o 6mers.csv -mv mean.npy -sv std.npy

seekr_pearson 6kmer_chr1_11_pca.csv 6kmer_chr1_11_pca.csv -o example_vs_self.csv
seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv
# cp communities.csv communities_5.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 10 -s 0 -c communities_10.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 9 -s 0 -c communities_9.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 8 -s 0 -c communities_8.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 7 -s 0 -c communities_7.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 6 -s 0 -c communities_6.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 4 -s 0 -c communities_4.csv
# seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -n 3 -s 0 -c communities_3.csv


# seekr_graph example_vs_self.csv 0.13 -g ./3_community/example_vs_self.gml -n 3 -s 0 -c ./3_community/communities_3.csv
# seekr_graph example_vs_self.csv 0.13 -g ./4_community/example_vs_self.gml -n 4 -s 0 -c ./4_community/communities_4.csv
# seekr_graph example_vs_self.csv 0.13 -g ./6_community/example_vs_self.gml -n 6 -s 0 -c ./6_community/communities_6.csv
# seekr_graph example_vs_self.csv 0.13 -g ./7_community/example_vs_self.gml -n 7 -s 0 -c ./7_community/communities_7.csv
# seekr_graph example_vs_self.csv 0.13 -g ./8_community/example_vs_self.gml -n 8 -s 0 -c ./8_community/communities_8.csv

# conda deactivate


