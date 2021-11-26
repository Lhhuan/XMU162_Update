perl 07_filter_chr1_11.pl #过滤"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz 中的chr1_11,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/chr1_11/Chr1_11_Tissue_merge_0.176_extend_sorted_merge.bed.gz";

bash 6kmer.sh 
perl 08_transform_kmer_result.pl
Rscript 09_pie.R

perl 11_split_cluster_and_homer.pl
perl 11_split_cluster_and_homer_other_n_community.pl

cp -r /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/homer/ /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/5_community/homer
Rscript 12_heatmap_feature.R 
Rscript 12_heatmap_feature_other_n_community.R

source activate  /home/huanhuan/miniconda3
Rscript 13_tf_enrichment.R 
Rscript 14_tf_enrichment_Plot.R 
perl extract_

#test
findMotifsGenome.pl "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/6_community/homer/communities_1.bed" /home/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/6kmer/test/ -size 200 -mset  vertebrates


/share/apps/Homer/bin/findMotifsGenome.pl "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/6_community/homer/communities_1.bed" hg19 /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/6kmer/test2/ -size 200  -preparse -preparsedDir /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/6kmer/test3/



findMotifsGenome.pl "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/6_community/homer/communities_1.bed" /home/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/6kmer/test/ -size 200 -mset  vertebrates