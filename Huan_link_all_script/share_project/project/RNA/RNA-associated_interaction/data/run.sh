#------data from starbase    
    #RNA-RNA #Get the Data for NcRNA-RNA Interaction Network.
    #Retrieve the interaction network of ncRNA-RNA identified from high-throughput sequencing data of RNA-RNA interactome
    nohup curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=mRNA&RNA=all&interNum=1&expNum=1&cellType=all' > ENCORI_hg19_mRNA-RNA_all.txt &
    nohup curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=sncRNA&RNA=all&interNum=1&expNum=1&cellType=all' > ENCORI_hg19_sncRNA-RNA_all.txt &
    nohup curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=lncRNA&RNA=all&interNum=1&expNum=1&cellType=all' > ENCORI_hg19_lncRNA-RNA_all.txt &
    nohup curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=miRNA&RNA=all&interNum=1&expNum=1&cellType=all' > ENCORI_hg19_miRNA-RNA_all.txt &
    nohup curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=pseudogene&RNA=all&interNum=1&expNum=1&cellType=all' > ENCORI_hg19_pseudogene-RNA_all.txt &
    #-------------------------
    #miRNA-Target #Get MiRNA Target Data.
    #Retrieve data of the miRNA-target interactions by intersecting the predicting target sites of miRNAs with binding sites of Ago protein
    curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=mRNA&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=PITA,RNA22&target=all&cellType=HeLa' > ENCORI_hg19_CLIP-seq_all_miRNA_mRNA.txt
    #---------------------------------------
    #Degradome-RNA #Get Data for MiRNAs Cleavage Events.
    #Retrieve data for cleavage events of miRNAs on genes supported by degradome-seq data
    curl 'http://starbase.sysu.edu.cn/api/degradomeRNA/?assembly=hg19&geneType=mRNA&miRNA=all&degraExpNum=1&target=all&cellType=all' > ENCORI_hg19_degradome-seq_all.txt
    #-------------------------
    #CeRNA
    #Get Data of the CeRNA Networks.
    #Retrieve the ceRNA networks from thousands of interactions of miRNA-targets supported by CLIP-seq data.
    #download failed
    curl 'http://starbase.sysu.edu.cn/api/ceRNA/?assembly=hg19&geneType=mRNA&ceRNA=all&miRNAnum=2&family="all"&pval=0.01&fdr=0.01&pancancerNum=0' >ENCORI_hg19_ceRNA-network_all.txt
    #--------------------------
    #RBP-Target
    #Get RBP Target Data.
    #Retrieve data of RBP-RNA interactions supported by the binding sites of RBPs derived from CLIP-seq data
    curl 'http://starbase.sysu.edu.cn/api/RBPTarget/?assembly=hg19&geneType=mRNA&RBP=all&clipExpNum=5&pancancerNum=0&target=all&cellType=all' > ENCORI_hg19_RBPTarget_all.txt
#-----------------data from RISE
    #----------------human Whole RNA-RNA interactome
    #The RRI data are in the BEDPE file format, defined in BEDTOOLS.  Each entry denote an RRI(RNA-RNA Interaction).
        # The following are the details of each column.   Here RNA1 and RNA2 are the two interacting RNAs :
        # 1. chr1:  The name of the chromosome of the first RNA.
        # 2. start1:  The zero-based starting position of the interacting region on chr1.
        # 3. end1:  The one-based ending position of interacting region on chr1.
        # 4. chr2:  The name of the chromosome of the second RNA.
        # 5. start2:  The zero-based starting position of the interacting region on chr2.
        # 6. end2:  The one-based ending position of the interacting region on chr2.
        # 7. rise_id:  The interaction entry ID in RISE.
        # 8. score:  All " . " here (just adapt to file format).
        # 9. strand1:  Defines the strand for the RNA1. Either '+'' or '-'.   Note: if some entry source provides no strand information, we use " . " instead.
        # 10. strand2:  Defines the strand for the RNA2. Either '+' or '-'.   Note: if some entry source provides no strand information, we use " . " instead.
        # 11. gene_id1:   Ensembl gene id of RNA1.   Note: if there is no available Ensembl gene id, use " . " instead.
        # 12. gene_name1.
        # 13. gene_id2:   Ensembl gene id of RNA2.  Note: if there is no available Ensembl gene id, use " . " instead.
        # 14. gene_name2.
        # 15. gene_type1:   Gene type of RNA1.
        # 16. gene_type2:   Gene type of RNA2.
        # 17. genomic_context1:   Genomic context of RNA1 interaction region.   Note: if there is no interacting region information, use " . " instead.
        # 18. genomic_context2:   Genomic context of RNA2 interaction region.  Note: if there is no interacting region information, use " . " instead.
        # 19. method:   RNA interaction detected method.
        # 10. species.
        # 21. cell line:   Cell line of the detected interaction.   Note: for interaction from other databases, there is no cell line information, use " . " instead.
        # 22. pubmed_id:   Pubmed id of the reference paper where the interaction from.
    wget http://rise.life.tsinghua.edu.cn/download-data/rise_human_all.txt.gz
#-----------------------------------------------------
#data from NPIDB(Nucleic acid–Protein Interaction DataBase ) #http://npidb.belozersky.msu.ru/download.html
#
    #---------- 
    wget -c http://npidb.belozersky.msu.ru/data/pdb_new/archives/all.tar
#--------------------------------------------------

#----------------------------------------
#data from RNAinter：a comprehensive resource for human RNA-associated (RNA–RNA/RNA–protein) interaction
    #header
    # - RNAInter ID: unique identifier for each entry in RNAInter database
    # - Interactor1: interactor1 in current entry
    # - ID1: ID of interactor1
    # - Category1: category of interactor1
    # - Species1: organism name of interactor1
    # - Interactor2: interactor2 in current entry
    # - ID2: ID of interactor2
    # - Category2: category of interactor2
    # - Species2: organism name of interactor2
    # - Score: the integrative confidence score of the current entry
    #-------All RNA-RNA interactions:
    #header: RNAInter ID	Interactor1	ID1	Category1	Species1	Interactor2	ID2	Category2	Species2	Score
    wget -c http://www.rna-society.org/raid/download/RNA-RNA.zip
    unzip RNA-RNA.zip #RNA-RNA.txt
    #------All RNA-Protein interactions
    wget -c http://www.rna-society.org/raid/download/RNA-Protein.zip
    unzip RNA-Protein.zip #RNA-Protein.txt
    #------All RNA-DNA interactions
    wget -c http://www.rna-society.org/raid/download/RNA-DNA.zip
    unzip RNA-DNA.zip #RNAInter_interaction_RD.txt