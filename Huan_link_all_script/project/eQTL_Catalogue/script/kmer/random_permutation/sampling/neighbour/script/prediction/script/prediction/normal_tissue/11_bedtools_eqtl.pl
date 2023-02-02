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
# my $pm = Parallel::ForkManager->new(36); ## 设置最大的线程数目

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
    my $f1 = "$dir/08_eqtl_sorted.bed.gz";
    my $hotspot= "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/predicted_region/predicted_regions_win5000_large_than6.bed";
    system "bedtools intersect -a $f1 -b $hotspot -wa -wb |gzip > $dir/11_predicted_region_eqtl.bed.gz";
    
    print "$tissue\tfinish\n";
}