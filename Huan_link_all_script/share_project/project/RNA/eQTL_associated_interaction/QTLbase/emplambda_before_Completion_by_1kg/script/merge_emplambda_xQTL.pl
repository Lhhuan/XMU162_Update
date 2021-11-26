#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt.gz

#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"，"QTL"#中../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt merge到一起，得
#得../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my %hash1;


# my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
my @QTL1s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");

my $fo1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "emplambda\tt\tchr\txQTL\n";
# my $fo1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
foreach my $QTL1(@QTL1s){
    my $f1 = "../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

    my %hash2;
    print "$QTL1\n";
    while(<$I1>)
    {
        chomp;
        unless(/emplambda/){
            my @f = split/\t/;
            my $emplambda =$f[0];
            my $t =$f[1]; 
            my $pos = $t; #snp pos
            my $chr =$f[2]; #snp chr
            unless ($emplambda =~/NA/){ 
                print $O1 "$_\t$QTL1\n";
            }
        }
    }

    close($I1);
}
close($O1);


my @QTL2s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
my $fo2 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "emplambda\tt\tchr\txQTL\n";
# my $fo1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
foreach my $QTL2(@QTL2s){
    my $f1 = "../output/ALL_${QTL2}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL2}.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

    my %hash2;
    print "$QTL2\n";
    while(<$I1>)
    {
        chomp;
        unless(/emplambda/){
            my @f = split/\t/;
            my $emplambda =$f[0];
            my $t =$f[1]; 
            my $pos = $t; #snp pos
            my $chr =$f[2]; #snp chr
            unless ($emplambda =~/NA/){ 
                print $O2 "$_\t$QTL2\n";
            }
        }
    }
}