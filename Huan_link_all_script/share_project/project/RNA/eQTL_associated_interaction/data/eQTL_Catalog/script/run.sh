perl 01_normalized.pl #将../raw_data/cis_eQTL_table_conditional_ALL ../raw_data/trans_eQTL_table_conditional_ALL normalized 成../normalized/01_normalized.txt.gz
zless ../normalized/01_normalized.txt.gz |sort -u| wc -l #3642317-1
