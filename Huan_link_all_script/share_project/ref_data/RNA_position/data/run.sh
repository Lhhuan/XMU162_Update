wget -c ftp://mirbase.org/pub/mirbase/19/genomes/hsa.gff3
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.long_noncoding_RNAs.gff3.gz
wget -c https://tcga.xenahubs.net/download/probeMap/hugo_gencode_good_hg19_V24lift37_probemap
mv hsa.gff3 hg19_v19_miRNA_position.txt
mv hugo_gencode_good_hg19_V24lift37_probemap hg19_mRNA_position.txt