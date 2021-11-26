#将../output/09_fisher_test_result.txt和../output/08_start_neighbour_events.txt merge到一起
#得../output/10_start_neighbour_events_fisher.txt 
#得显著的文件../output/10_start_neighbour_events_fisher_significant.txt  
#得不显著的文件 ../output/10_start_neighbour_events_fisher_not_significant.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use File::Copy;

my $f1 ="../output/08_start_neighbour_events.txt";#
my $f2 ="../output/09_fisher_test_result.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 ="../output/10_start_neighbour_events_fisher.txt ";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 ="../output/10_start_neighbour_events_fisher_significant.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 ="../output/10_start_neighbour_events_fisher_not_significant.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my $header = "start_symbol\tp_value\tevents_in_neighbour";
print $O1 "$header\n";
print $O2 "$header\n";
print $O3 "$header\n";

my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^symbol/){
        my $symbol = $f[0];
        my $events = $f[1];
        $hash1{$symbol}=$events;
    }
}


while(<$I2>)
{
    chomp;
    my @f= split /\t/;
    unless(/^symbol/){
        my $symbol = $f[0];
        my $p_value = $f[1];
        if (exists $hash1{$symbol}){
            my $events = $hash1{$symbol};
            my $output = "$symbol\t$p_value\t$events";
            print $O1 "$output\n";
            if ($p_value<0.05){
                print $O2 "$output\n";
            }
            else{
                print $O3 "$output\n";
            }
        }
    }
}

