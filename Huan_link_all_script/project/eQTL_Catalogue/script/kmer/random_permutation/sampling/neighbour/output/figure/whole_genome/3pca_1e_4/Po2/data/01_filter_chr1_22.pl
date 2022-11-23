#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
# use List::Util qw/max min/;
# use List::Util qw/sum/;
# use Parallel::ForkManager;


my $f1 = "PolII_hg38.bed";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); 
my $fo1 = "PolII_hg38_chr1_22.bedgraph";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my %hash1;
foreach my $i(1..22){
    $i="chr${i}";
    $hash1{$i}=1;
}
my %hash2;

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr=$f[0];  
    my $start =$f[1];
    my $end =$f[2];
    my $id= $f[3];
    my $signal=$f[4];
    if(exists $hash1{$chr}){
        print $O1 "$chr\t$start\t$end\t$signal\n";
        # unless(exists $hash2{$k1}){
        #     print $O1 "$_\n";
        #     $hash2{$k2}=1;
        # }
    }
}
