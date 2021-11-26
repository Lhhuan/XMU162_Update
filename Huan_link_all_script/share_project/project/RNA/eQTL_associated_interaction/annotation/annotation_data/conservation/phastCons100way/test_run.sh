#/share/apps/R_depends/bin/bash
perl 01_normalized_data.pl

zless merge_phastCons100way.bed.gz | sort -k1,1 -k2,2n | gzip > merge_phastCons100way_sorted.bed.gz

rm merge_phastCons100way.bed.gz

