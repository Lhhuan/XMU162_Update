#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
# use Parallel::ForkManager;

my $f1 = "/share/data0/GTEx/data/GTEx_Analysis_v8_trans_eQTL/GTEx_Analysis_v8_trans_eGenes_fdr05.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo3 = "/home/huanhuan/project/GTEx/output/Tissue_merge/trans_eQTL.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

while(<$I1>)
{
    chomp;
    unless(/tissue_id/){
        my @f = split/\t/;
        my $variant_id =$f[6];
        my @t =split/\_/,$variant_id;
        my $chr =$t[0];
        my $start =$t[1];
        my $end =$start +1;
        print $O3 "$chr\t$start\t$end\n";
    }
}

