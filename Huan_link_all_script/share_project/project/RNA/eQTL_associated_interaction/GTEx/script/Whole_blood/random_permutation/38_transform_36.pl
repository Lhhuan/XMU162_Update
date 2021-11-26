#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Whole_Blood";

# my @d_types = ("up","down");
my @states = (25);

# my @types = ("random_permutation","original_random");


my @types = ("original_random","hotspot");
foreach my $type(@types){
    foreach my $state(@states){
        my $output_dir = "../../../output/${tissue}/Cis_eQTL/enrichment/interval_18/ALL/" ;
        my $f1 = "$output_dir/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}.txt";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        # my $fo1 = "$output_dir/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}_refine.txt";
        # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        my $fo2 = "$output_dir/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}.csv";
        open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
        # open my $O1, "| gzip >$fo1" or die $!;
        # print $O1 "chromatin_state\tnumber\trandom_number\n";    
        my %hash1;

        while(<$I1>)
        {
            chomp;
            my @f = split/\t/;
            my $ot = join(",",@f);
            print $O2 "$ot\n";
            #-----
        }
    }
}

