#raw_data中的数据，rsid 的position与dbsnp中记录的差1，比如rs10003238	chr4	47548775，而dbSNP中记录的是GRCh37.p13 chr 4	NC_000004.11:g.47548776A>G
perl 01_normalized.pl #将../raw_data/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt 和../raw_data/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt normalized 成../normalized/01_SEEQTL_normalized.txt.gz
zless ../normalized/01_SEEQTL_normalized.txt.gz | wc -l #127149-1
