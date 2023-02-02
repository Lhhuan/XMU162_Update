#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;
my $dir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output";
my @levels=("mean","max");
my $cluster_inputdir = "$dir";

# my @markers=("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
my @markers=("H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3");
foreach my $level(@levels){
     my $output_dir ="$dir/figure/markers/${level}";
    foreach my $marker(@markers){
        my $marker_bw = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/${marker}_merge_${level}_signalvalue.bw";
        $ENV{'cluster_inputdir'} = $cluster_inputdir ;
        $ENV{'marker'} = $marker ;
        $ENV{'level'} = $level ;
        $ENV{'marker_bw'} = $marker_bw ;
        $ENV{'output_dir'} = $output_dir ;
        system "bash plot_chip_seq_point_umap.sh";
        print "$marker\t$level\n";
    }
}

my @cistromes = ("Human_FACTOR","Human_CHROMATIN_Accessibility");
foreach my $level(@levels){
    my $output_dir ="$dir/figure/markers/${level}";
    foreach my $marker1(@cistromes){
        my $marker_bw ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/${marker1}/merge_${level}_signalvalue.bw";
        my $marker= $marker1;
        $marker =~ s/Human_FACTOR/TFBS/g;
        $marker =~ s/Human_CHROMATIN_Accessibility/CA/g;
        $ENV{'cluster_inputdir'} = $cluster_inputdir ;
        $ENV{'marker'} = $marker ;
        $ENV{'level'} = $level ;
        $ENV{'marker_bw'} = $marker_bw ;
        $ENV{'output_dir'} = $output_dir ;
        system "bash plot_chip_seq_point_umap.sh";
        print "$marker\t$level\n";
    }
}

foreach my $level(@levels){
    my $output_dir ="$dir/figure/markers/${level}";
    my $marker_bw ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38/normal_cell_line_ctcf_merge_${level}_signalvalue.bw";
    my $marker="CTCF";
    $ENV{'cluster_inputdir'} = $cluster_inputdir ;
    $ENV{'marker'} = $marker ;
    $ENV{'level'} = $level ;
    $ENV{'marker_bw'} = $marker_bw ;
    $ENV{'output_dir'} = $output_dir ;
    system "bash plot_chip_seq_point_umap.sh";
    print "$marker\t$level\n";
}
# # my @pol2=()
#=============pol2
foreach my $level(@levels){
    my $output_dir ="$dir/figure/markers/${level}";
    my $marker_bw ="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/Po2/data/PolII_hg38_chr1_22_merge_${level}_signalvalue.bw";
    my $marker="PolII(CD4+_T_cells)";
    $ENV{'cluster_inputdir'} = $cluster_inputdir ;
    $ENV{'marker'} = $marker ;
    $ENV{'level'} = $level ;
    $ENV{'marker_bw'} = $marker_bw ;
    $ENV{'output_dir'} = $output_dir ;
    system "bash plot_chip_seq_point_umap.sh";
    print "$marker\t$level\n";
}

foreach my $level(@levels){
    my $output_dir ="$dir/figure/markers/${level}";
    my $marker_bw ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/POL2/merge_${level}_signalvalue.bw";
    my $marker="PolII(cistromeDB)";
    $ENV{'cluster_inputdir'} = $cluster_inputdir ;
    $ENV{'marker'} = $marker ;
    $ENV{'level'} = $level ;
    $ENV{'marker_bw'} = $marker_bw ;
    $ENV{'output_dir'} = $output_dir ;
    system "bash plot_chip_seq_point_umap.sh";
    print "$marker\t$level\n";
}