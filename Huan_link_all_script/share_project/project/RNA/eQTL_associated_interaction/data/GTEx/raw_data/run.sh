GTEx V8
#-------------
    #Multi-Tissue eQTL Data
    #eGENE
    #Metasoft results for all tissues in the V8 release. For information about the column headers, please see
    wget -c https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_eqtl_data/GTEx_Analysis_v8.metasoft.txt.gz
    #------------------------------
    #Single-Tissue cis-QTL Data
    #eGene and significant variant-gene associations based on permutations. The archive contains a *.egenes.txt.gz and *.signif_variant_gene_pairs.txt.gz file for each tissue. 
    #Note that the *.egenes.txt.gz files contain data for all genes tested; to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
    wget -c https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
    tar xvf GTEx_Analysis_v8_eQTL.tar
    #---------------------------
    #Lookup table for all variants genotyped in GTEx, with chromosome positions, REF and ALT alleles, RS IDs from dbSNP 151, GTEx variant IDs (constructed as chr_pos_alt_ref_build), 
    #and hg19 liftover variant ID, for all variants in release V8. #hg38转hg37
    wget -c https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
