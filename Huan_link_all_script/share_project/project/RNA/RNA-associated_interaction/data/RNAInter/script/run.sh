cat header.txt >../output/RNA-Protein.txt
cat ../raw_data/RNA-Protein.txt >> ../output/RNA-Protein.txt
cat header.txt > ../output/RNA-RNA.txt
cat ../raw_data/RNA-RNA.txt >> ../output/RNA-RNA.txt
cat header.txt > ../output/RNA-DNA.txt
cat ../raw_data/RNAInter_interaction_RD.txt >> ../output/RNA-DNA.txt
perl 01_normalized_interaction.pl ##normalize ../output/${geneType}.txt, 得../normalized/01_${geneType}_normalized.txt;
cat ../normalized/01_RNA-Protein_normalized.txt| cut -f6,15 | sort -u | wc -l  #23577247 -1
cat ../normalized/01_RNA-RNA_normalized.txt| cut -f6,15 | sort -u | wc -l #2361458 -1
cat ../normalized/01_RNA-DNA_normalized.txt |cut -f6,15 | sort -u | wc -l  #93447 -1