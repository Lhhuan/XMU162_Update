#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;
my $f1 = "files.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "01_need_file_info.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $file_name= $f[0];
    my $detail = $f[1];
    my @ts=split/;/,$detail;
    if($detail=~/ treatment=None/ && $detail=~/ type=bed/){
        print $O1 "$_\n";
        foreach my $t(@ts){
            if($t =~/size/){
                print "$t\n";
            }
        }
        
    }
    # # if($file_name=~/\*.bed.gz$/){
    #     foreach my $t(@ts){
    #         print "$t\n";
    #         if($t=~/ treatment=None/ && $t=~/ type=bed/){
    #             print $O1 "$_\n";
                
    #         }
    #     }
    # } 
}

