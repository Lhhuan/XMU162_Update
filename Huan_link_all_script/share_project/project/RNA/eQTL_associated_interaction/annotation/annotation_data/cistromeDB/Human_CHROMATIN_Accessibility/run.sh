wget -c http://cistrome.org/db/batchdata/R56Q7GGRZEY7L4PH4RA9.tar.gz
tar -zxvf R56Q7GGRZEY7L4PH4RA9.tar.gz
tar -zxvf human_ca.tar.gz
rm R56Q7GGRZEY7L4PH4RA9.tar.gz
rm human_ca.tar.gz
cd human_ca
gzip *

zless merge_pos_info_narrow_peak.bed.gz | sort -k1,1 -k2,2n |gzip >merge_pos_info_narrow_peak_sorted.bed.gz
rm merge_pos_info_narrow_peak.bed.gz

perl 03_filter_chr1_22.pl