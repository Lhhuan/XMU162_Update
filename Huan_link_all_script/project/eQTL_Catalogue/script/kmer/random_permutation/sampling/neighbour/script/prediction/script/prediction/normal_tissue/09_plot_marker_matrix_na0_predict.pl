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

my $f1 = "./output/02_tissue_level_Marker_source.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^tissue_ontology_id/){
        my @f = split/\t/;
        my $tissue_ontology_id =$f[0];
        my $refine_tissue_label2 =$f[1];
        my $marker =$f[2];
        my $ID =$f[3];
        my $k="$refine_tissue_label2\t$marker";
        unless($marker=~/HISTONE_MARK_AND_VARIANT/){
            push @{$hash1{$refine_tissue_label2}},$marker;
        }
    }
}
# my $input_file="../predicted_region/predicted_regions_win5000_large_than6.bed";
my @histone_markers = ("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
foreach my $histone_marker(@histone_markers){
    $hash2{$histone_marker}=1;
}
my $level="mean";
foreach my $k(sort keys %hash1){
    # my $pid = $pm->start and next; #开始多线程
    my $tissue=$k;
    my $cluster_inputdir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/$tissue";
    my @markers =@{$hash1{$k}};
    my @uni_markers=uniq(@markers);
    foreach my $marker(@uni_markers){
        my $marker1=$marker;
        $marker1 =~s/Human_FACTOR/TFBS/g;
        $marker1 =~s/Human_CHROMATIN_Accessibility/CA/g;
        if(exists $hash2{$marker}){
            my $dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/marker/histone_hg38";
            my $output_dir ="/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/markers_plot/na0/${level}";
            if (-e $output_dir){
                print "$output_dir exist\n";
            }
            else{
                mkpath($output_dir);
            }
            my $marker_bw = "$dir/${marker}_merge_${level}_signalvalue.bw";
            $ENV{'marker1'} = $marker1 ;
            $ENV{'level'} = $level ;
            $ENV{'marker_bw'} = $marker_bw ;
            $ENV{'output_dir'} = $output_dir ;
            $ENV{'cluster_inputdir'} = $cluster_inputdir ;
            my $out_file = "${output_dir}/${marker1}_${level}.gz";
            unless(-e $out_file){
                system "bash plot_chip_seq_point_umap_miss_as_zero.sh";
                print "$marker\t$level\n";
            }
        }
        else{
            my $dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/marker";
            my $output_dir ="/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/markers_plot/na0/${level}";
            if (-e $output_dir){
                print "$output_dir exist\n";
            }
            else{
                mkpath($output_dir);
            }
            my $marker_bw = "$dir/${marker}_merge_${level}_signalvalue.bw";
            $ENV{'marker1'} = $marker1 ;
            $ENV{'level'} = $level ;
            $ENV{'marker_bw'} = $marker_bw ;
            $ENV{'output_dir'} = $output_dir ;
            $ENV{'cluster_inputdir'} = $cluster_inputdir ;
            my $out_file = "${output_dir}/${marker1}_${level}.gz";
            unless(-e $out_file){
                system "bash plot_chip_seq_point_umap_miss_as_zero.sh";
                print "$marker\t$level\n";
            } 
        }
    }
    # $pm->finish;  #多线程结束
}

