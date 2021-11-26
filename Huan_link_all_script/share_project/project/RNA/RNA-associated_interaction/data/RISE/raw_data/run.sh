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