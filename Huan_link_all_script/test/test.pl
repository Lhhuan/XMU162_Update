#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "D21GP02757.hg19_multianno.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $fo1 = "filter.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


while(<$I1>)
{
    chomp;
    if(/^Chr/){
        print $O1 "$_\n";
    }
    else{
        my @f=split/\t/;
        my $a= $f[91];
        my $b= $f[96];
        if($b=~/A|T|C|G/ && $a >80){
            print $O1 "$_\n";
        }
    }
}


