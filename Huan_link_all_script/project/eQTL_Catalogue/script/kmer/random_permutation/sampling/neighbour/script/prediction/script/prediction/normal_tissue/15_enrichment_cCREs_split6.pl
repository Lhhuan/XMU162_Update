#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use File::Path;
# use List::MoreUtils ':all';
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(36); ## 设置最大的线程数目

my $f1 = "./output/02_tissue_level_Marker_source.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^tissue_ontology_id/){
        my @f = split/\t/;
        my $refine_tissue_label2 =$f[1];
        $hash1{$refine_tissue_label2}=1;
    }
}
my $no_gap_file = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted_nogap.bed.gz";
my @types=("CTCF","enhD","enhP","K4m3","prom");

# foreach my $tissue(sort keys %hash1){
my @all_tissue = sort keys %hash1;
foreach my $tissue(@all_tissue[30..35]){
    # my $pid = $pm->start and next; #开始多线程
    my $output_dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/$tissue/Enrichment/cCREs";
    if(! -d $output_dir){
        mkpath($output_dir);
    }
    else{
        print "$output_dir exists\n";
    }
    #================
    foreach my $type(@types){
        my $cre  = "/share/data0/QTLbase/huan/Cis_Regulatory_Elements/encodeCcreCombined_sorted_${type}.bed.gz";
        foreach my $i(0..2){
            my $hotspot = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/$tissue/predicted_Class${i}.bed";
            my $output = "$output_dir/${type}_C${i}_enrichment.txt";
            chdir("/home/huanhuan/tools/gonomics/src/github.com/vertgenlab/gonomics-main/cmd/overlapEnrichments") ;
            system "go run overlapEnrichments.go exact $hotspot $cre $no_gap_file $output";
            print "$tissue\t$i\t$type\n";
        }
    }
    # $pm->finish;  #多线程结束
}

