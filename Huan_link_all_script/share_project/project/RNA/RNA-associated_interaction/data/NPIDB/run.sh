wget -c http://npidb.belozersky.msu.ru/data/pdb_new/lists/ListOfComplexes.txt

#data from NPIDB(Nucleic acid–Protein Interaction DataBase ) #http://npidb.belozersky.msu.ru/download.html
#
    #---------- 
    wget -c http://npidb.belozersky.msu.ru/data/pdb_new/archives/all.tar
    tar -xf all.tar
    mkdir sub_pdb
    mv *.pdb ./sub_pdb
    tar -xzf all.tar.gz
    mkdir sub_pdb2
    mv *.pdb ./sub_pdb2
#--------------------------------------------------