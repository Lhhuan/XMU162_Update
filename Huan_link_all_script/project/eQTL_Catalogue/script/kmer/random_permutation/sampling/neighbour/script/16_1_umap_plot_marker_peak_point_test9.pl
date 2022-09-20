#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;
system "cp /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/homer_200/*.bed /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/markers/";
my @levels=("max");
my $cluster_inputdir = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/markers/";

# my @markers=("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
my @markers=("H3K27ac");
foreach my $level(@levels){
     my $output_dir ="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/markers/point/${level}";
    foreach my $marker(@markers){
        my $marker_bw = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/${marker}_merge_${level}_signalvalue.bw";
        $ENV{'cluster_inputdir'} = $cluster_inputdir ;
        $ENV{'marker'} = $marker ;
        $ENV{'level'} = $level ;
        $ENV{'marker_bw'} = $marker_bw ;
        $ENV{'output_dir'} = $output_dir ;
        system "bash plot_chip_seq_point_umap_cluster9.sh";
        print "$marker\t$level\n";
    }
}