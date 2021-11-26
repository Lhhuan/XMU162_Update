#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中../output/ALL_${QTL1}/count_number_by_emplambda_in_different_interval_7.3_all_${QTL1}.txt merge到一起，得
#得../output/count_number_by_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my %hash1;


# my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
my @QTL1s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");

my $fo1 = "../output/count_number_by_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
foreach my $QTL1(@QTL1s){
    my $f1 = "../output/ALL_${QTL1}/count_number_by_emplambda_in_different_interval_7.3_all_${QTL1}.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    print $O1 "emplambda\tnumber\txQTL\n";
    my %hash2;
    # print "$QTL1\n";
    while(<$I1>)
    {
        chomp;
        unless(/emplambda/){
            my @f = split/\t/;
            my $emplambda =$f[0];
            my $number =$f[1]; 
            my $interval = $f[2];
            if ($interval =~/1000/){ 
                print $O1 "$emplambda\t$number\t$QTL1\n";
                print "$interval\n";
            }
        }
    }
}