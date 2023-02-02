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
my @states=("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF_Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies");
my $fo1 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/ChromHMM/merge_15_state_enrichment.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Method\tFilename1\tFilename2\tLenElements1\tLenElements2\tOverlapCount\tDebugCheck\tExpectedOverlap\tEnrichment\tEnrichPValue\tDepletePValue\tstate\tsample\tcluster\n";

foreach my $state(@states){
    # print "$state\n";
    my $dir = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/ChromHMM/2cluster/${state}";
    opendir (DIR, $dir) or die "can't open the directory!";
    my @files = readdir DIR;
    foreach my $file (@files) {
        # print "$file\n";
        if ( $file =~ /txt/) {
            my @t= split/_/,$file;
            my $sample = $t[0];
            my $cluster=$t[-2];
            # print "$cluster\n";
            my $f1 = "$dir/$file";
            open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
            # open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            while(<$I1>)
            {
                chomp;
                my @f = split/\t/;
                unless(/^#/){
                    print $O1 "$_\t$state\t$sample\t$cluster\n";
                }
            }
        }
    }
}




