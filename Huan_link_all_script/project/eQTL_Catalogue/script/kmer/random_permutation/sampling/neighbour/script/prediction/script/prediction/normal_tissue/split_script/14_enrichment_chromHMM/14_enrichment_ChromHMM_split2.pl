#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
use File::Path;
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
        my $marker =$f[2];
        my $ID = $f[3];
        unless($marker =~/HISTONE_MARK_AND_VARIANT|Human_CHROMATIN_Accessibility|Human_FACTOR/){
            push @{$hash1{$refine_tissue_label2}},$ID;
            # print "$marker\n";
        }

    }
}
my $no_gap_file = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted_nogap.bed.gz";
my @states=("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies");
my @all_tissue = sort keys %hash1;
foreach my $tissue(@all_tissue[6..11]){
# foreach my $tissue(sort keys %hash1){
    # my $pid = $pm->start and next; #开始多线程
    my $output_dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/$tissue/Enrichment/ChromHMM";
    if(! -d $output_dir){
        mkpath($output_dir);
    }
    else{
        print "$output_dir exists\n";
    }
    my @IDs = @{$hash1{$tissue}};
    my @uniq_IDs =uniq(@IDs);
    foreach my $state(@states){
        my $state1 = $state;
        $state1 =~s/\//_/g;
        foreach my $ID(@uniq_IDs){
            my $chromHMM = "/share/data0/QTLbase/huan/chromatin_states/15_coreMarks/$state1/${ID}_15_${state1}_coreMarks_hg38lift_dense.bed.gz";
            foreach my $i(0..2){
                my $hotspot = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/$tissue/predicted_Class${i}.bed";
                my $output = "$output_dir/${ID}_${state1}_C${i}_enrichment.txt";
                chdir("/home/huanhuan/tools/gonomics/src/github.com/vertgenlab/gonomics-main/cmd/overlapEnrichments") ;
                system "go run overlapEnrichments.go exact $hotspot $chromHMM $no_gap_file $output";
                print "$tissue\t$ID\t$i\t$state1\n";
            }
        }
    }
    # $pm->finish;  #多线程结束
}
