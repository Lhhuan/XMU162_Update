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
            $hash1{$refine_tissue_label2}=1;
        }
    }
}

foreach my $tissue(sort keys %hash1){
    print "$tissue\tstart\n";
    my $dir ="/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}";
    my $f1 = "$dir/predicted_key_info.txt.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $fo1 = "$dir/predicted_Class0.bed";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    my $fo2 = "$dir/predicted_Class1.bed";
    open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
    my $fo3 = "$dir/predicted_Class2.bed";
    open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
    while(<$I1>)
    {
        chomp;
        my @f=split/\t/;
        unless(/^hotspot/){
            my $hotspot = $f[0];
            my $predict_class =$f[1];
            my @t=split/:/,$hotspot;
            my $chr= $t[0];
            my @cc = split/-/,$t[1];
            my $start =$cc[0];
            my $end=$cc[1];
            my $output="$chr\t$start\t$end";
            if($predict_class=~/0/){
                print $O1 "$output\n";
            }elsif($predict_class=~/1/){
                print $O2 "$output\n";                
            }else{
                print $O3 "$output\n";
            }
        }
    }
    print "$tissue\tfinish\n";
}