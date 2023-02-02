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

my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/082_final_ca_tf_TCGA_sample.txt" ;
my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/092_final_histone_TCGA_mark.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 


my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^TCGA/){
        my @f = split/\t/;
        my $TCGA =$f[0];
        my $factor = $f[1];
        my $cancer_cell = $f[2];
        my $cell_line_refine =$f[3];
        my $k = "$TCGA/$cell_line_refine";
        push @{$hash1{$k}},$factor;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^TCGA/){
        my @f = split/\t/;
        my $TCGA =$f[0];
        my $factor = $f[1];
        my $cancer_cell = $f[2];
        my $cell_line_refine =$f[3];
        my $k = "$TCGA/$cell_line_refine";
        push @{$hash1{$k}},$factor;
    }
}


my $input_file="../../predicted_region/predicted_regions_win5000_large_than6.bed";
my $level="mean";
my @ks = sort keys %hash1;
foreach my $k(@ks[156..161]){
    # my $pid = $pm->start and next; #开始多线程
    my $output_dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/${k}/marker/markers_na0/${level}";
    if (-e $output_dir){
        print "$output_dir exist\n";
    }
    else{
        mkpath($output_dir);
    }
    my @markers =@{$hash1{$k}};
    my @uni_markers=uniq(@markers);
    foreach my $marker(@uni_markers){

        my $marker_bw = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/cell_type_level/$k/${marker}_merge_${level}_signalvalue.bw";
        $ENV{'marker'} = $marker ;
        $ENV{'level'} = $level ;
        $ENV{'marker_bw'} = $marker_bw ;
        $ENV{'output_dir'} = $output_dir ;
        $ENV{'input_file'} = $input_file ;
        my $out_file = "${output_dir}/${marker}_${level}.bed.gz";
        unless(-e $out_file){
            system "bash computeMatrix_sample_miss_as0.sh";
            print "$k\t$marker\t$level\n";
        }
    }
    # $pm->finish;  #多线程结束
}

