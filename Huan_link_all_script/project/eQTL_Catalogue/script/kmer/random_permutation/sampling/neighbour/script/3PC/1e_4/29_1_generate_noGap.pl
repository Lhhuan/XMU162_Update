#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;


my $f1 = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo2 = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted_nogap.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $chr =$f[0];
    my $end=$f[1];
    # print $O2 "$chr\t0\t$end\tname$.\n";
    print $O2 "$chr\t0\t$end\n";
}
close($O2);



