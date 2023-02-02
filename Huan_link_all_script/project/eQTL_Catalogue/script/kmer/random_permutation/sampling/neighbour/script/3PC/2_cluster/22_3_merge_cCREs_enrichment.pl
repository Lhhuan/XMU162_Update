#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use List::MoreUtils ':all';
# use Parallel::ForkManager;
# my $pm = Parallel::ForkManager->new(70); ## 设置最大的线程数目
my @types=("CTCF","enhD","enhP","K4m3","prom");
my $fo1 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/22_3_merge_cCREs_enrichment.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Method\tFilename1\tFilename2\tLenElements1\tLenElements2\tOverlapCount\tDebugCheck\tExpectedOverlap\tEnrichment\tEnrichPValue\tDepletePValue\ttype\tcluster\n";


my $dir = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/enrichment/";

foreach my $type (@types) {
    # print "$file\n";
    foreach my $i(1..2){
        my $f1 = "$dir/${type}_C${i}_enrichment.txt";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        # open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
        while(<$I1>)
        {
            chomp;
            my @f = split/\t/;
            unless(/^#/){
                print $O1 "$_\t$type\tC$i\n";
            }
        }
    }
}




