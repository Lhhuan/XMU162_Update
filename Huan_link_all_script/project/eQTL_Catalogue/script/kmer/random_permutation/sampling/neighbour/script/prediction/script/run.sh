bash 01_get_kmer_for_negative.sh
"bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed $fo1 |gzip >$gc_file