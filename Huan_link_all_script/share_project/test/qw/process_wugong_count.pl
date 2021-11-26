#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
# use List::Util qw/max min/;
# use Parallel::ForkManager;


my $f1 = "dominant.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "final_dominant.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";



my %hash1;
while(<$I1>)
{
    chomp;
    if (/^D19GP00025/){
        print $O1 "Gene\t$_\n";
    }
    else{
        my @f = split/\t/;
        my $gene= $f[0];
        my $gene2= $gene;
        $gene2 =~ s/\:.*//g;
        print "$gene2\n";
        my $
        # my $sample = 
        # push @{$hash1{$}}
    }
  
}
