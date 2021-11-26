#------download data, data from starbase   
    #miRNA-Target #Get MiRNA Target Data.
    #Retrieve data of the miRNA-target interactions by intersecting the predicting target sites of miRNAs with binding sites of Ago protein
    perl 01_download_miRNA-Target.pl 
    #-----------------------------------------
    #Degradome-RNA #Get Data for MiRNAs Cleavage Events.
    #Retrieve data for cleavage events of miRNAs on genes supported by degradome-seq data
    perl 02_download_Degradome-RNA.pl
    #-------------------------------------------------
    #RNA-RNA #Get the Data for NcRNA-RNA Interaction Network.
    #Retrieve the interaction network of ncRNA-RNA identified from high-throughput sequencing data of RNA-RNA interactome
    perl 03_download_RNA-RNA.pl
    #-------------------------
    #CeRNA
    #Get Data of the CeRNA Networks.
    #Retrieve the ceRNA networks from thousands of interactions of miRNA-targets supported by CLIP-seq data.
    perl 04_download_CeRNA.pl
    #--------------------------
    #RBP-Target
    #Get RBP Target Data.
    #Retrieve data of RBP-RNA interactions supported by the binding sites of RBPs derived from CLIP-seq data
    perl 05_download_RBP-Target.pl
#--------------------------------------
#-------------------normalize data interaction_name1 为miRNA name或者gene name;interaction_name2 为RBP或者gene name
    perl 06_normal_miRNA-Target.pl  ##normalize "../raw_data/ENCORI_hg19_CLIP-seq_all_miRNA-${geneType}_${program}_miRNA-Target.txt", 得../normalized/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized.txt,
    #得normalized fail的文件../output/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized_fail.txt
    perl 07_normal_Degradome-RNA.pl #normalize ../raw_data/ENCORI_hg19_degradome-seq_all_Degradome-${geneType}_Degradome-RNA.txt, 得../normalized/07_ENCORI_hg19_degradome-seq_all_normalized.txt,
    #得normalized fail的文件../output/07_ENCORI_hg19_degradome-seq_all_normalized_fail.txt
    perl 08_normal_RNA-RNA.pl ##normalize "../raw_data/ENCORI_hg19_${geneType}-RNA_RNA-RNA_all.txt", 得../normalized/08_ENCORI_hg19_RNA-RNA_normalized.txt;
    perl 09_normal_CeRNA.pl #normalize "../raw_data/ENCORI_hg19_ceRNA-${geneType}-network_all.txt", 得../normalized/09_ENCORI_hg19_CeRNA_normalized.txt;
    #CeRNA 信息量不大
    perl 10_normal_RBP-Target.pl ##normalize "../raw_data/ENCORI_hg19_RBP-${geneType}_RBP-Target.txt", 得../normalized/10_ENCORI_hg19_RBP-target_normalized.txt;

    cat ../normalized/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized.txt | cut -f6,15 | sort -u | wc -l  #327289-1
    cat ../normalized/07_ENCORI_hg19_degradome-seq_all_normalized.txt | cut -f6,15 | sort -u | wc -l  #162929-1
    cat ../normalized/08_ENCORI_hg19_RNA-RNA_normalized.txt | cut -f6,15 | sort -u | wc -l #80001-1
    cat ../normalized/09_ENCORI_hg19_CeRNA_normalized.txt |cut -f6,15 | sort -u | wc -l #1125238 -1
    cat ../normalized/10_ENCORI_hg19_RBP-Target_normalized.txt | cut -f6,15 | sort -u | wc -l #47397 -1

