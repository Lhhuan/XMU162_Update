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