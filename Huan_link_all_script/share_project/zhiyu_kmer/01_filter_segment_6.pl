#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $f1 ="nodes_feats.tsv";

open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#----------------------
my $fo1 = "01_filter_nodes.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
while(<$I1>)
{
    chomp;
    unless(/^chr/){
        my $a= $_;
        $a =~ s/\.0//g;
        my @f= split/\t/,$a;

        my $chr =$f[0];
        my $start =$f[1];
        my $end =$f[2];
        print $O1 "chr${chr}\t$start\t$end\n";
    }

}
# print "$group\t$cutoff\n";
      