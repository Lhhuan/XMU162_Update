Rscript 01_random_sapmling.R
perl 02_plot_marker_peak_point.pl

bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed ./output/markers/01_random_c1.bed |gzip >./output/GC_content/01_random_c1_gc_content.bed.gz
bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed ./output/markers/01_random_c2.bed |gzip >./output/GC_content/01_random_c2_gc_content.bed.gz
bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed ./output/markers/01_random_c3.bed |gzip >./output/GC_content/01_random_c3_gc_content.bed.gz

Rscript 03_boxplot_gc_content.R