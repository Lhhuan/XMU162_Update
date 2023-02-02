#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use File::Path;
use Parallel::ForkManager;
use List::MoreUtils ':all';
my $pm = Parallel::ForkManager->new(36); ## 设置最大的线程数目

my $f1 = "../output/02_cancer_cellLine_marker.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^TCGA/){
        my @f = split/\t/;
        my $TCGA =$f[0];
        my $marker =$f[1];
        my $cell_line_refine =$f[2];
        my $tmp_tissue =$f[3];
        push @{$hash1{$tmp_tissue}},$marker;
    }
}

my $level="mean";
my @all_tissue = sort keys %hash1;
# foreach my $tissue(sort keys %hash2){
foreach my $k(@all_tissue[180..185]){
    # my $pid = $pm->start and next; #开始多线程
    my $tissue=$k;
    my $cluster_inputdir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/$tissue";
    my @markers =@{$hash1{$k}};
    my @uni_markers=uniq(@markers);
    foreach my $marker(@uni_markers){
        my $output_dir ="/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/${tissue}/markers_plot/na0/${level}";
        if (-e $output_dir){
            print "$output_dir exist\n";
        }
        else{
            mkpath($output_dir);
        }
        my $marker_bw = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/cell_type_level/$k/${marker}_merge_${level}_signalvalue.bw";
        my $marker1 = $marker;
        $marker1 =~s/Human_FACTOR/TFBS/g;
        $marker1 =~s/Human_CHROMATIN_Accessibility/CHROMATIN_Accessibility/g;
        # $ENV{'marker'} = $marker ;
        $ENV{'marker1'} = $marker1 ;
        $ENV{'level'} = $level ;
        $ENV{'marker_bw'} = $marker_bw ;
        $ENV{'output_dir'} = $output_dir ;
        $ENV{'cluster_inputdir'} = $cluster_inputdir ;
        my $out_file = "${output_dir}/${marker}_${level}.bed.gz";
        unless(-e $out_file){
            system "bash plot_chip_seq_point_umap_miss_as_zero.sh";
            print "$marker\t$level\n";
        }
    }
    # $pm->finish;  #多线程结束
}

