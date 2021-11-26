wget -c  http://cistrome.org/db/batchdata/24KRO157XZ5Y204IEVFN.tar.gz
tar -zxvf 24KRO157XZ5Y204IEVFN.tar.gz
tar -zxvf human_factor.tar.gz
rm 24KRO157XZ5Y204IEVFN.tar.gz
cd human_factor
gzip *
rm ../human_factor.tar.gz

zless merge_pos_info_narrow_peak.bed.gz | sort -k1,1 -k2,2n | gzip >merge_pos_info_narrow_peak_sorted.bed.gz
rm merge_pos_info_narrow_peak.bed.gz