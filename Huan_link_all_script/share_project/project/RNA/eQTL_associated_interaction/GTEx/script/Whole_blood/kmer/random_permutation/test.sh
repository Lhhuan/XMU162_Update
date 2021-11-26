
source activate /home/huanhuan/miniconda3
cd /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/original_random/
# $cutoff = 0.176
# $group ="hotspot"
# $tissue = "Whole_Blood"
for i in $(seq 1 2)
do 

# sorted_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.fa
    seekr_kmer_counts sorted_${i}_resemble_Whole_Blood_segment_hotspot_cutoff_0.176.fa -o ${i}_6mers_uc_us_no_log.csv --log2 none -uc -us
    echo $i
done