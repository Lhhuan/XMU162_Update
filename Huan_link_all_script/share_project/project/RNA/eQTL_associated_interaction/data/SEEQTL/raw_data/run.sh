    #Processed genotype and expression data
    wget -c https://www.seeqtl.org/gbrowse2/data/genotype_expression_data.tar.gz
    tar zxvf genotype_expression_data.tar.gz
    #-------------------------------
    #-----------this maybe need
    # eQTLs of meta-analysis (consensus eQTLs) from HapMap human lymphoblastoid cell lines (sample summary)
    wget -c http://www.seeqtl.org/gbrowse2/data/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt
    wget -c https://www.seeqtl.org/gbrowse2/data/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt
    #raw_data中的数据，rsid 的position与dbsnp中记录的差1，比如rs10003238	chr4	47548775，而dbSNP中记录的是GRCh37.p13 chr 4	NC_000004.11:g.47548776A>G