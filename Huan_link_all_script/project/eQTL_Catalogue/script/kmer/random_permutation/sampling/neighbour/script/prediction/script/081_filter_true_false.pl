#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;


my $f1 = "../output/08_warm_region_predict_key_info.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "../output/081_warm_region_predict_hotspot_true.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/081_warm_region_predict_hotspot_false.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my %hash1;

while(<$I1>)
{
    chomp;
    unless(/^Chr/){
        my @f= split/\t/;
        my $chr =$f[0];
        my $start=$f[1];
        my $end =$f[2];
        my $class =$f[4];
        if($class >0){
            print $O1 "$chr\t$start\t$end\n";
        }
        else{
            print $O2 "$chr\t$start\t$end\n";
        }
    }

}
# close($O1);



