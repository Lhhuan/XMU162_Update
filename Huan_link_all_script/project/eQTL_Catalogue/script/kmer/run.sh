bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -fo "/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/Tissue_merge_hotspot_extend_18_snp.fa"

cd /share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/

source activate seekr_source
seekr_kmer_counts Tissue_merge_hotspot_extend_18_snp.fa  -o  6mers_uc_us_no_log.csv --log2 none -uc -us
gzip 6mers_uc_us_no_log.csv


Rscript 07_count_kmer_count.R
