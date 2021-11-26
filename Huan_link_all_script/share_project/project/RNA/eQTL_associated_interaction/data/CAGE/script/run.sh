perl 01_normalized.pl #将../raw_data/cage_association_summary_statistics.txt normalized 成../normalized/01_normalized.txt.gz
zless ../normalized/01_normalized.txt.gz |sort-u| wc -l #3533718-1
