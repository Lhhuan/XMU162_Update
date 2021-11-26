wget -c http://cistrome.org/db/batchdata/GTYPP2KEMBOVQL3DDGS2.tar.gz
tar -zxvf GTYPP2KEMBOVQL3DDGS2.tar.gz
tar -zxvf human_hm.tar.gz
rm GTYPP2KEMBOVQL3DDGS2.tar.gz
rm human_hm.tar.gz
cd human_hm
gzip *

zless merge_pos_info_narrow_peak.bed.gz | sort -k1,1 -k2,2n |gzip >merge_pos_info_narrow_peak_sorted.bed.gz
rm merge_pos_info_narrow_peak.bed.gz
perl 03_filter_chr1_22.pl #