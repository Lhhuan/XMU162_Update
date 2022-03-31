bedtools intersect -a $fo1 -b $hotspot -wo |cut -f4-6  >$hotspot_in_region
bedtools intersect -a $hotspot_in_region -b $eqtl_gene -wo |gzip>$qtl_in_hotspot_region #hotspot snp tss
bedtools makewindows -b $fo1 -w $win -i winnum | gzip >$win_region
zless $qtl_in_hotspot_region |awk -v OFS='\t' '{print $4,$5,$6,$7}' |sort -u |sort -k1,1 -k2,2n | gzip > $qtl_in_r
zless $qtl_in_hotspot_region |awk -v OFS='\t' '{print $8,$9,$10,$7}' |sort -u |sort -k1,1 -k2,2n | gzip > $tss_in_r
bedtools intersect -a $qtl_in_r -b $win_region -wa -wb |cut -f1-4,8|gzip > $qtl_in_win
bedtools intersect -a $tss_in_r -b $win_region -wa -wb |cut -f1-4,8|gzip > $tss_in_win