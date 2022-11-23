#=================hi-c H1-hESC cross chr
    perl 24_adjust_hic_format.pl 
    zless "/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_same_chr_adjust.bed.gz" |sort -k1,1 -k2,2n |gzip >"/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_same_chr_adjust_sorted.bed.gz"
    bedtools intersect -a "/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_same_chr_adjust_sorted.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster4DNFIQYQWPF5_5000.bed.gz
    Rscript 25_boxplot_Hic.R

    #=================hi-c H1-hESC cross chr
    perl 24_adjust_hic_format_cross_chr.pl
    zless "/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_cross_chr_adjust.bed.gz" |sort -k1,1 -k2,2n |gzip >"/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_cross_chr_adjust_sorted.bed.gz"
    bedtools intersect -a "/share/data0/QTLbase/huan/hic/4DNFIQYQWPF5_5000.ginteractions_cross_chr_adjust_sorted.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster4DNFIQYQWPF5_5000_cross_chr.bed.gz

    Rscript 25_boxplot_Hic_cross_chr.R

    #==============hi-c HFFc6
    perl 24_adjust_hic_format_cross_chr_HFFc6.pl
    zless "/share/data0/QTLbase/huan/hic/4DNFIFLJLIS5_5000.ginteractions_cross_chr_adjust.bed.gz" |sort -k1,1 -k2,2n |gzip >"/share/data0/QTLbase/huan/hic/4DNFIFLJLIS5_5000.ginteractions_cross_chr_adjust_sorted.bed.gz"
    bedtools intersect -a "/share/data0/QTLbase/huan/hic/4DNFIFLJLIS5_5000.ginteractions_cross_chr_adjust_sorted.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster_4DNFIFLJLIS5_5000_cross_chr.bed.gz
    Rscript 25_boxplot_Hic_cross_chr_HFFc6.R
    #==============micro-C XL H1-hESC cross chr
    perl 24_adjust_microc_format_cross_chr_H1-hESC.pl

    bedtools intersect -a "/share/data0/QTLbase/huan/hic/4DNFI2TK7L2F_5000.ginteractions_cross_chr_adjust_sorted.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster_4DNFI2TK7L2F_5000_cross_chr.bed.gz
    Rscript 25_boxplot_microc_cross_chr_H1-hESC.R
    #==================micro-C XL HFFc6 cross chr
    perl 24_adjust_microc_format_cross_chr_HFFc6.pl

    bedtools intersect -a "/share/data0/QTLbase/huan/hic/4DNFI18Q799K_5000.ginteractions_cross_chr_adjust_sorted.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster_4DNFI18Q799K_5000_cross_chr.bed.gz
    Rscript 25_boxplot_microc_cross_chr_HFFc6.R