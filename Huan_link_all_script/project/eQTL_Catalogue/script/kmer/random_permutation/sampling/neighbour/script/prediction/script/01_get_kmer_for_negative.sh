bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -fo "../output/01_all_negative.fa"
cd ../output/

source activate seekr_source
seekr_kmer_counts 01_all_negative.fa  -o  01_all_negative_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip 01_all_negative_6mers_uc_us_no_log.csv
rm 01_all_negative.fa

conda deactivate 
bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" |cut -f1-3,5 |gzip >01_all_negative.gc_content.bed.gz